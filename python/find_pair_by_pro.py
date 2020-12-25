#!/usr/bin/env python
import sys
from math import fabs
from collections import defaultdict

fn = sys.argv[1]
shift = float(sys.argv[2])

def find_pair(peak_lst, gap):
    tol = 0.01
    peak_lst = sorted(peak_lst)

    bank = []
    bank.append([peak_lst.pop(0),1])
    while len(peak_lst)>0:
        if fabs(peak_lst[0]-bank[-1][0])<tol:
            #average
            p, n = bank[-1]
            bank[-1][0] = (peak_lst.pop(0) + p*n)/(n+1)
            bank[-1][1] += 1
        else:
            bank.append([peak_lst.pop(0),1])

    tmp = []
    tmp.append( bank.pop(0) )
    hits = []
    while len(bank)>0:
        #discard p in tmp that has no pair
        while len(tmp)>0 and fabs(bank[0][0]-tmp[0][0])>(gap+tol*2):
            tmp.pop(0)
        if len(tmp)>0 and fabs(bank[0][0]-tmp[0][0]-gap)<(tol*2):
            #hits
            hits.append([tmp.pop(0), bank[0]])
            tmp.append(bank.pop(0))
        else:
            tmp.append(bank.pop(0))
    return hits

lines = open(fn, 'r').readlines()
peaks = defaultdict(list)
for l in lines:
    elems = l.strip().split()
    pep = elems[0]+":"+elems[1]
    mass = float(elems[2])
    peaks[pep].append(mass)

for pep in peaks.keys():
    if len(peaks[pep])>1:
        hits = find_pair(peaks[pep], shift)
        if len(hits)>0:
            for hit in hits:
                print pep, hit[0][0], hit[0][1], hit[1][0], hit[1][1]
