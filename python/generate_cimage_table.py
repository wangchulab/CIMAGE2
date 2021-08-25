#!/usr/bin/env python
import sys, os
from collections import defaultdict
from math import fabs

tabs = {}
fp_lst = {}
lines = open(sys.argv[1], 'r').readlines()
for l in lines:
  es = l.strip().split()
  fp_lst[es[0]] = open(es[0]+".table.txt", 'w')
  tabs[es[0]] = defaultdict(list)
  for tag in es[1:]:
    ts = tag.split(':')
    for t in ts[1:]:
      t2 = t.split(',')
      tabs[es[0]][ts[0]].append( (t2[0], int(t2[1])) )
    
#print tabs
def clean_header( lines ):
  while lines[0][0] == "!":
    lines = lines[1:]
  return lines

map_tag = {}
map_ndx = {}
scr_path =  os.path.split(os.path.realpath(__file__))[0]
lines = open(scr_path+"/light.default.table.txt", 'r').readlines()
lines = clean_header(lines)
for n, tag in enumerate(lines[0].split()):
  map_tag[n] = tag
  map_ndx[tag] = n

def output( fp, aa_tag, title="" ):
  n = len(aa_tag.keys())
  out_str = title
  for i in xrange(n):
    out_str += '\t' + str(aa_tag[i])
  fp.write(out_str+'\n')

#output col name
for k in fp_lst.keys():
  output(fp_lst[k], map_tag)

for l in lines[1:]:
  es = l.strip().split('\t')
  title = es[0]
  ref_tag = {}
  for n, a in enumerate(es[1:]):
    ref_tag[n] = int(a)
  for k in fp_lst.keys():
    #copy ref tab
    aa_tag = {}
    for i in ref_tag.keys():
      aa_tag[i] = ref_tag[i]
    #modify
    if title in tabs[k].keys():
      tab = tabs[k][title]
      for aa, nn in tab:
        #print aa, nn
        i = map_ndx[aa]
        aa_tag[i] += nn

    output(fp_lst[k], aa_tag, title)

#close
for k in fp_lst.keys():
  fp_lst[k].close()

