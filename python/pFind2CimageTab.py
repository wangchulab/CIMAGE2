#!/usr/bin/env python
import sys
from collections import defaultdict

# Usage: input: pFind_light.spectra, pFind_heavy.spectra
# pFind2CimageTab.py name L-mods H-mods other-mods

map_mod_mass = {}
map_mod_sign = {}
L_labels = []
H_labels = []

ipi_lst = defaultdict(None)
scan_lst = []
cross_tab = defaultdict(list)

sample_name = sys.argv[1]
for m in sys.argv[2].split('|'):
  es = m.split(':')
  L_labels.append( es[0] )
  map_mod_mass[ es[0] ] = float(es[1])
  map_mod_sign[ es[0] ] = es[2]
for m in sys.argv[3].split('|'):
  es = m.split(':')
  H_labels.append( es[0] )
  map_mod_mass[ es[0] ] = float(es[1])
  map_mod_sign[ es[0] ] = es[2]
for m in sys.argv[4].split('|'):
  es = m.split(':')
  map_mod_mass[ es[0] ] = float(es[1])
  map_mod_sign[ es[0] ] = es[2]

#check
#print map_mod_mass, map_mod_sign, L_labels, H_labels
for tag in [ "light", "heavy" ]:#
  fn = "pFind_" + tag + ".spectra"
  lines = open(fn, 'r').readlines()
  #make tag
  ts = {}
  for n, k in enumerate(lines[0].split('\t')):
    ts[k] = n
  #parse data
  #print ts
  for l in lines[1:]:
    label = None #"light"
    es = l.split('\t')
    if len(es)>10:
      #fraction
      fn = es[ts["File_Name"]]
      frac = fn.split('.')[0]
      frac = frac.replace( sample_name+"_", "" )
      #params
      charge = es[ts["Charge"]]
      scan = es[ts["Scan_No"]]
      mass = es[ts["Exp.MH+"]] ##
      mods = es[ts["Modification"]]
      score = float(es[ts["Final_Score"]])
      #id
      mprots = es[ts["Proteins"]][:-1]
      prots = mprots.split('/')
      p_mark = 0
      while prots[p_mark][:4] == "REV_" or prots[p_mark][:4] == "CON_":
        p_mark += 1
        if p_mark == len(prots):
          p_mark -= 1
          break
      prot = prots[p_mark]
      if prot[:4] == "CON_" or prot[:4] == "REV_":
        uni = prot.split('.')[0]
        syb = "decoy"
        disc = mprots
      else:
        if len(prot.split('|'))<3: continue
        uni = prot.split('|')[1]
        syb = prot.split('|')[2].split('_')[0]
        disc = mprots
      #seq
      seq = es[ts["Sequence"]]
      pos = es[ts["Positions"]]
      pos = pos.split('/')[0]
      pos = pos.split(',')
      #mod
      skip = False
      for i, mod in enumerate(mods.split(';')[:-1]):
        mod = mod.split(',')
        if tag=="light" and mod[1] in L_labels:
          label = "light"
        if tag=="heavy" and mod[1] in H_labels:
          label = "heavy"
        #if mod[1] not in map_mod_sign.keys():
        #  skip = True
        #  break
        sign = map_mod_sign[mod[1]]
        seq = seq[:int(mod[0])+i] + sign + seq[int(mod[0])+i:]
      #if skip: continue
      seq = seq.replace("-", "")
      #check
      if label != tag: continue
      #cap
      seq = pos[1] + '.' + seq + '.' + pos[2]
      #output
      key = uni+":"+seq+":"+charge+":"+frac
      scan_str = key + " " + sample_name + " " + scan + " " + label + "\n"
      scan_lst.append( scan_str )
      ipi_lst[ uni ] = ( syb, disc )
      cross_tab[ key ].append( (mass, scan, score ) )

#save ipi_name.table
ipiout = open("ipi_name.table", 'w')
ipiout.write("name\n")
for uni in ipi_lst.keys():
  syb, disc = ipi_lst[ uni ]
  ipiout.write( uni )
  ipiout.write( "\t" )
  ipiout.write( syb )
  ipiout.write( " " )
  ipiout.write( disc )
  ipiout.write( "\n" )
ipiout.close()

#save all_scan.table
scanout = open("all_scan.table", 'w')
scanout.write("key run scan HL\n")
for scan in scan_lst:
    scanout.write(scan)
scanout.close()

#save cross_scan.table
crossout = open("cross_scan.table", 'w')
crossout.write( "key mass %s\n" % (sample_name) )
for scan_core in cross_tab.keys():
    uni = scan_core.split(":")[0]
    rank = sorted( cross_tab[scan_core], key=lambda x: x[2] )
    neutral_mass = rank[-1][0]
    id_scan = rank[-1][1]
    crossout.write(scan_core+" "+str(neutral_mass)+" "+str(id_scan)+"\n")
crossout.close()

