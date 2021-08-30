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

marks = [c for c in sys.argv[3]]

def load_cimage(fn):
  #print "#Loading from cimage_combine output"
  lines = open(fn, 'r').readlines()
  tags = lines[0].split()
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
  lines = open(fn, 'r').readlines()
  for l in lines:
    yield (pro, pep, R)
  return

def load_msfragger(fn):
  assert(0)
  print "Loading from MSFragger peptide.tsv"
  lines = open(fn, 'r').readlines()
  for l in lines:
    es = l.split('\t')
    pro = "pro"
    pep = "XXX"
    R = 1.0
    yield (pro, pep, R)
  return

def extract_motif_pos(pro, pep):
  #return 
  lst = []
  try:
    seq = db_seq[pro]
  except:
    print "Warning: protein", pro, "not found"

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
    for motif, pos in extract_motif_pos(pro, pep):
      print pro, motif, pos, R
elif fn[-4:] == ".txt":
  #cimage_combine results
  for (pro, pep, R) in load_cimage(fn):
    for motif, pos in extract_motif_pos(pro, pep):
      print pro, motif, pos, R
elif fn[-4:] == ".tsv":
  #msfragger results
  #for (pro, pep, R) in load_msfragger(fn):
  #  break
  print "MSF to be done"
else:
  print "Wrong input file!"

