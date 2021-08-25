#!/usr/bin/env python
import atomParams as at
import sys, os
from math import fabs

def clean_header( lines ):
  while lines[0][0] == "!":
    lines = lines[1:]
  return lines

scr_path = os.path.split(os.path.realpath(__file__))[0]

map_tag = {}
lines = open(scr_path+"/light.default.table.txt", 'r').readlines()
lines = clean_header(lines)
for n, tag in enumerate(lines[0].split()):
  map_tag[n] = tag

map_ref_mass = {}
for l in lines[1:]:
  es = l.strip().split()
  name = es[0]
  mass = 0.0
  for n, a in enumerate(es[1:]):
    mass += at.mass[map_tag[n]] * int(a)
  map_ref_mass[name] = mass

lines = open(sys.argv[1], 'r').readlines()
lines = clean_header(lines)
for l in lines[1:]:
  es = l.strip().split()
  name = es[0]
  mass = 0.0
  for n, a in enumerate(es[1:]):
    mass += at.mass[map_tag[n]] * int(a)
  ref_mass = map_ref_mass[name]
  dm = fabs( ref_mass - mass )
  if dm > 0.01:
    print name, dm

