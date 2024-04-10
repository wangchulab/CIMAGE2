#!/usr/bin/env python
import sys
import Bio
from Bio import SeqIO

#get mod pos and peptide/motif
#Usage: cmd tsv/txt fasta marks

Lwin = 20

db_seq = {}
for rec in SeqIO.parse(sys.argv[2], "fasta"):
  cur_id = rec.id
  cur_seq = rec.seq
  if cur_id[:8] == "Reverse_": continue
  k = cur_id.split('|')[1]
  db_seq[k] = cur_seq

#marks = [c for c in sys.argv[3]]
marks = []
map_mark = {}
for c in sys.argv[3].split('|'):
  if "," in c:
    mark, label = c.split(',')
    map_mark[label] = mark
    marks.append(mark)
  else:
    marks.append(c)

def load_cimage(fn):
  #print "#Loading from cimage_combine output"
  lines = open(fn, 'r').readlines()
  #tags = lines[0].split()
  for l in lines[1:]:
    es = l.split('\t')
    if es[0] != " ":
      #skip title
      continue
    pro = es[1]
    if pro[:8] == "Reverse_": continue
    pep = es[4]
    R = float(es[6])
    R2 = float(es[7])
    if R<0.0001:
      R=R2
      #for BAU
      #if R2>=0.07: R=1.0/R2
      #elif R2>0.0001: R=15.0
      #else:
      #  continue
    if R<0.0001: continue
    yield (pro, pep, R)
  return

def load_pfind(fn):
  print "Loading from pFind all_result"
  lines = open(fn, 'r').readlines()
  t = {}
  nm = 0
  for n,tag in enumerate(lines[0].strip().split('\t')):
    t[tag] = n
  for l in lines[1:]:
    es = l.split('\t')
    pro = es[t["Protein AC"]].split('/')[0]
    if pro[:4] == "REV_": continue
    pro = pro.split('|')[1]
    pep = es[t["Sq"]]
    #print pep
    mods = es[t["Mod_Sites"]]
    R = 0.0
    #print mods
    for mod in mods.split(';')[:-1]:
      #print mod
      loc, mod = mod.split(',')
      loc = int(loc)
      for l in map_mark.keys():
        marker = map_mark[l]
        if mod.split('#')[0] == l:
          #match
          #print loc, pep[loc+nm]
          pep = pep[:loc+nm]+marker+pep[loc+nm:]
          yield (pro, pep, R)
  return

def load_msfragger(fn):
  print "Loading from MSFragger peptide.tsv(Ionquant TBD)"
  lines = open(fn, 'r').readlines()
  t = {}
  nm = 0
  for n,tag in enumerate(lines[0].strip().split('\t')):
    t[tag] = n
  for l in lines[1:]:
    es = l.split('\t')
    pro = es[t["Protein ID"]]
    pep_raw = es[t["Peptide"]]
    R = 0.0
    mods = es[t["Assigned Modifications"]]
    for mod in mods.split(','):
      mod = mod.strip()
      if len(mod)<3: continue
      tmp = mod.split('(')
      try:
        loc = int(tmp[0][:-1])
      except:
        loc = 0
      mod = tmp[0][-1]+"("+tmp[1]
      for l in map_mark.keys():
        marker = map_mark[l]
        if mod == l:
          pep = pep_raw[:loc+nm]+marker+pep_raw[loc+nm:]
          yield (pro, pep, R)
  return

def extract_motif_pos(pro, pep):
  #print "-", pro, pep
  lst = []
  try:
    seq = db_seq[pro]
  except:
    print "Warning: protein", pro, "not found"
    return lst

  locs = []
  L = len(pep)
  clean_pep = ""
  j = 0
  for i, c in enumerate(pep):
    if c.isupper():
      clean_pep += c
      j += 1
    elif c in marks:
      locs.append(j)

  L_all = len(seq)
  start = seq.find(clean_pep)
  if start<0:
    print "Warning: peptide", pep, "can't be located in", pro
  else:
    for loc in locs:
      pt = loc+start
      motif = ""
      #first half
      if pt <= Lwin:
        motif += "-"*(Lwin-pt+1)
        motif += seq[:pt]
      else:
        motif += seq[pt-Lwin-1:pt]
      #second half
      if pt+Lwin > L_all:
        motif += seq[pt:]
        motif += "-"*(pt+Lwin-L_all)
      else:
        motif += seq[pt:pt+Lwin]
      lst.append( (motif, pt) )

  return lst

fn = sys.argv[1]
if fn[:10] == "all_result":
  #reading pfind results
  for (pro, pep, R) in load_pfind(fn):
    #print pro, pep
    for motif, pos in extract_motif_pos(pro, pep):
      print pro, motif, pos, R
elif fn[-4:] == ".txt":
  #cimage_combine results
  for (pro, pep, R) in load_cimage(fn):
    for motif, pos in extract_motif_pos(pro, pep):
      print pro, motif, pos, R
elif fn[-4:] == ".tsv":
  #msfragger results
  for (pro, pep, R) in load_msfragger(fn):
    for motif, pos in extract_motif_pos(pro, pep):
      print pro, motif, pos, R
else:
  print "Wrong input file!"

