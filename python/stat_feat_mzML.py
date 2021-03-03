#!/usr/bin/env python

import sys
from math import fabs, log
from pyteomics import mzml
import pandas as pd
#from lshash.lshash import LSHash
import collections
import numpy as np

fn_feats = sys.argv[1] + ".features.tsv"
fn_mzML = sys.argv[1] + ".mzML"

ppm = 1e-6
tol = 10.0 * ppm #>2
rt_gap = 0.5 #min
mz_cutoff = 0.01
iso_shifts = [0, 1, 2, 3] #mono, +1 and +2
delta_C13 = 1.0033548
bin_size = 20.0

class feature:
    def __init__(self):
        self.rtStart = 0.0
        self.rtEnd = 0.0
        self.mz = 0.0
        self.rtApex = 0.0
        self.mzApex = 0.0
        self.intApex = 0.0
        self.intSum = 0.0
        #add
        self.charge = 0
        self.id = -1
        self.ms1_peak = []
        self.ms2_peak = []
        self.ms2_scan = []

    def output(self):
        return "(MZ: %6.4f, RT: %6.4f, Apex: %6.4f, rtStart: %4.3f, rtEnd: %4.3f)" \
        % (self.mz, self.rtApex*60.0, self.intApex, self.rtStart, self.rtEnd)

#load all features (mz and rt)
lines = open(fn_feats, 'r').readlines()
tags = { k:n for n,k in enumerate(lines[0].strip().split('\t')) }

def get_hash(mz):
    return round(mz*bin_size)

feats = []
dict_mz_feat = collections.defaultdict(list)

for n, l in enumerate(lines[1:]):
   es = l.strip().split('\t')
   c = int(es[tags['charge']])
   feat = feature()
   feat.rtStart = float(es[tags['rtStart']])
   feat.rtEnd = float(es[tags['rtEnd']])
   feat.mz = float(es[tags['mz']])
   feat.mzApex = float(es[tags['mostAbundantMz']])
   feat.rtApex = float(es[tags['rtApex']])
   feat.intApex = float(es[tags['intensityApex']])
   feat.intSum = float(es[tags['intensitySum']])
   feat.charge = c
   #save
   feat.id = n
   feats.append(feat)
   assert( n+1 == len(feats) )

   for s in iso_shifts:
       dM = delta_C13 *s /c
       hash_code = get_hash(feat.mz+dM)
       dict_mz_feat[hash_code-1].append(n)
       dict_mz_feat[hash_code].append(n)
       dict_mz_feat[hash_code+1].append(n)

num = len(feats)
print("Number_of_feats:", num)

N_ms2 = 0
N_hit = 0
print("Parsing mzML file")
for spectrum in mzml.read(fn_mzML):
    level = int(spectrum["ms level"])
    if level == 1: continue #skip ms1
    N_ms2 += 1

    base_mz = spectrum.get('base peak m/z')
    scanList = spectrum.get('scanList')
    mz_lst = spectrum.get('m/z array')
    inten = spectrum.get('intensity array')
    index = int(spectrum.get('index'))
    scan_id = index+1

    scan = scanList['scan']
    rt = float(scan[0]['scan start time'])
    precursor = float(scan[0]['[Thermo Trailer Extra]Monoisotopic M/Z:'])

    #features within range
    hash_key = get_hash(precursor)
    cands = dict_mz_feat[hash_key]
    N_t = 0
    if len(cands) > 0:
        hits = []
        for s in iso_shifts:
            hits += [c for c in cands if fabs(feats[c].mz - precursor + delta_C13*s/feats[c].charge) < mz_cutoff]
        for c in hits:
            rt0_s = feats[c].rtStart
            rt0_e = feats[c].rtEnd
            if rt0_s <= rt and rt <= rt0_e:
                feats[c].ms2_scan.append(scan_id)
                N_t += 1
    if N_t>0: N_hit += 1
    #print(N_t, len(cands))
print(N_hit, "ms2 out of", N_ms2)

#counts = []
N_feat = 0
N_hit = 0
for f in feats:
    N_feat += 1
    if len(f.ms2_scan)==0: continue
    N_hit += 1
    print("Feature", f.id+1)
    print(f.output())
    for scan in f.ms2_scan:
        print("MS2:", scan)
    #counts.append(len(f.ms2_scan))
print(N_hit, "features out of", N_feat)
