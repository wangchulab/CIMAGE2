#!/usr/bin/env python
import sys
from collections import defaultdict

flt = sys.argv[1]
res = sys.argv[2]
label = None
if len(sys.argv)>3:
  label = sys.argv[3]

lines = open( flt, 'r' ).readlines()
sel = defaultdict( int )
for l in lines:
  es = l.split('\t')
  if label == None:
    sel[es[1]] = 1
  else:
    if es[12] == label:
      sel[es[1]] = 1

lines = open( res, 'r' ).readlines()
print lines[0].strip()
for l in lines[1:]:
  l = l.strip()
  es = l.split()
  if sel[es[0]] == 1:
    print l
