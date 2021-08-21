#!/usr/bin/env python
import sys
from math import fabs
import numpy as np
from scipy.stats import hypergeom as hg
from scipy.stats import zscore

diffm = []
if len(sys.argv) >= 2:
  fn = sys.argv[1]
else:
  fn = "psm.tsv"
lines = open( fn, 'r' ).readlines()
t = {}
for n, tag in enumerate(lines[0].split('\t')):
  t[tag] = n
for line in lines[1:]:
  es = line.split('\t')
  diffm.append(float(es[t["Delta Mass"]]))

#tags
nAA = 21
tags = "ACDEFGHIKLMNPQRSTUVWY"
pos = {}
for i in range(nAA):
  pos[tags[i]] = i

flag_loc = False
loc_db = []
all_loc_AA = np.zeros([nAA])
all_AA = np.zeros([nAA])
for fn in sys.argv[2:]:
  #fn = sys.argv[2]
  flag_loc = True
  lines = open( fn, 'r' ).readlines()
  #parse tags
  t = {}
  for n,tag in enumerate(lines[0].split('\t')):
    t[tag] = n
  for l in lines[1:]:
    es = l.split('\t')
    best_loc_str = es[t["best_locs"]]
    loc_db.append( (float(es[t["massdiff"]]), best_loc_str) )
    for i in range(len(best_loc_str)-1): #ignore enzymetic bias
      c = best_loc_str[i]
      AA = c.upper()
      if c.islower():
        all_loc_AA[pos[AA]] += 1
      else:
        all_AA[pos[AA]] += 1
    #if 'u' in best_loc_str:
    #  print best_loc_str
  #print all_loc_AA / all_AA

sortm = sorted(diffm)
ndat = len(sortm)
nmid = 4
nwin = 9
dtol = 0.02

def get_loc( query_mass ):
  this_loc_AA = np.zeros([nAA])
  this_AA = np.zeros([nAA])
  for dm, best_loc_str in loc_db:
    if fabs( dm - query_mass ) < dtol:
      for i in range(len(best_loc_str)):
        c = best_loc_str[i]
        AA = c.upper()
        if c.islower():
          this_loc_AA[pos[AA]] += 1
        else:
          this_AA[pos[AA]] += 1
  #print this_loc_AA
  #print this_AA
  #print all_loc_AA
  #print all_AA
  vals = []
  sum_all = np.sum(all_AA)
  sum_loc = np.sum(this_AA)
  for i in range(nAA):
    c = tags[i]
    hit = this_loc_AA[i]
    hit_bg = this_AA[i]
    ref = all_loc_AA[i]
    ref_bg = all_AA[i] + 1
    #val = hg.pmf( hit, ref_bg, ref, hit_bg )
    #val1 = (hit/hit_bg) / ((ref-hit)/(ref_bg-hit_bg))
    val2 = (hit/ref_bg) * (sum_all/sum_loc)

    #use val2
    vals.append((c, val2))
  return sorted( vals, key=lambda x:x[1], reverse=True )

#for i in xrange(ndat-nwin):
#    dwin = sortm[i:i+nwin]
#    if dwin[-1]-dwin[0]<dtol:
#        print dwin[nmid]
#print sortm
ps = 0
pe = nwin
while ps < ndat and pe < ndat:
  if sortm[pe-1] - sortm[ps] > dtol:
    ps += 1
    pe += 1
    if pe == ndat: break
  else:
    dmid = np.median( sortm[ps:pe] )
    #extend window
    while sortm[pe-1]-dmid < dtol and pe < ndat:
      #extend window
      pe += 1
      dmid = np.median( sortm[ps:pe] )
    #adjust
    width = np.std( sortm[ps:pe] )
    while True:
      new_width = np.std( sortm[ps+1:pe+1] )
      if new_width < width:
        ps += 1
        pe += 1
        width = new_width
      else:
        break
    #for gnuplot
    if dmid <= 3.1 and dmid >= -0.1:
      print "%.4f" % dmid, pe - ps
    else:
      vals = get_loc( dmid )
      zsc = zscore([v[1] for v in vals])
      #print zsc
      output_str = "%.4f %d" % (dmid, pe-ps+1)
      for i, val in enumerate(vals):
        if i < 1 or zsc[i] > 1.2:
          output_str += " %s %.3f %.3f (%.5f)" % (val[0], val[1], zsc[i], width)
      print output_str

    ps = pe
    pe = ps + nwin

