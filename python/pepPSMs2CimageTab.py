#!/usr/bin/env python

import sys
from math import fabs
import numpy as np
from collections import defaultdict
import sys, glob
from Bio import SeqIO

mass_tol = 0.01

mass_norm_AA = {
"G" : 57.02146,
"A" : 71.03711,
"S" : 87.03203,
"P" : 97.05276,
"V" : 99.06841,
"T" :101.04768,
"C" :103.00918, #+57.021464
"L" :113.08406,
"I" :113.08406,
"N" :114.04293,
"D" :115.02694,
"Q" :128.05858,
"K" :128.09496,
"E" :129.04259,
"M" :131.04048,
"H" :137.05891,
"F" :147.06841,
"R" :156.10111,
"Y" :163.06333,
"W" :186.07931,
"n" :1.0078
}

L_labels = []
H_labels = []
norm_labels = []
norm_nterm_labels = []

ipi_lst = defaultdict(None)
scan_lst = []
cross_tab = defaultdict(list)

sample_name = sys.argv[1]
db_fn = sys.argv[2]

db_seq = {}
for record in SeqIO.parse(db_fn, "fasta"):
    if record.id[:8] == "Reverse_": continue
    pro = record.id.split()[0]
    db_seq[pro] = record.seq

#light
for mod in sys.argv[3].split('|'):
    elems = mod.split(':')
    if len(elems) == 3:
        AA = elems[0]
        mod_mass = float(elems[1])
        tag = elems[2]
        L_labels.append((AA, mod_mass, tag))
#heavy
for mod in sys.argv[4].split('|'):
    elems = mod.split(':')
    if len(elems) == 3:
        AA = elems[0]
        mod_mass = float(elems[1])
        tag = elems[2]
        H_labels.append((AA, mod_mass, tag))
#normal
for mod in sys.argv[5].split('|'):
    elems = mod.split(':')
    if len(elems) == 3:
        AA = elems[0]
        mod_mass = float(elems[1])
        tag = elems[2]
        if AA[0] == "n":
            norm_nterm_labels.append((AA, mod_mass, tag))
        else:
            norm_labels.append((AA, mod_mass, tag))

print "quant labels"
print L_labels
print H_labels
print "normal labels"
print norm_labels
print norm_nterm_labels

pepxml_dir = {}
pepxml_dir["light"] = sys.argv[6]
pepxml_dir["heavy"] = sys.argv[7]

def mark_seq(uAA, seq, dAA, markers):
    ans = ""
    if len(markers)==0:
        ans = uAA + "." + seq + "." + dAA
    else:
        #all_seq = uAA + "." + seq[:ndx+1] + "*" + seq[ndx+1:] + "." + dAA
        ans = uAA + "."
        prev = 0
        for marker, ndx in markers:
            ans = ans + seq[prev:ndx+1] + marker
            prev = ndx+1
        ans = ans + seq[prev:] + "." + dAA
    return ans

print "Parsing pepXML files..."
for tag in [ "light", "heavy" ]:
    #build mod list
    std_AA = [] #mark AA that's modified but not marked (SILAC)
    mod_list = []
    nterm_labels = []
    if tag == "light":
        for m in L_labels:
            if m[1] < mass_tol and m[2]=="-": #??
                std_AA.append(m[0])
            elif m[0][0]=="n":
                nterm_labels.append(m)
            else:
                mod_list.append(m)
    elif tag =="heavy":
        for m in H_labels:
            if m[1] < mass_tol and m[2]=="-": #??
                std_AA.append(m[0])
            elif m[0][0]=="n":
                nterm_labels.append(m)
            else:
                mod_list.append(m)
    else:
        print "Warning tag error!"
        break
    print "## ", tag
    print "mod list"
    print mod_list
    print "nter labels"
    print nterm_labels
    print "others"
    print std_AA

    print "Loading ..."
    fn = pepxml_dir[tag]+"/psm.tsv"
    print fn
    lines = open(fn, 'r').readlines()
    elems = lines[0].strip().split('\t')
    tags = {}
    for n, e in enumerate(elems):
        tags[e] = n
    #print tags
    #assert(0)
    for l in lines[1:]:
        elems = l.strip().split('\t')

        pro = elems[tags["Protein"]]
        disc = elems[tags["Protein Description"]].strip()
        disc = disc.replace("'", "", 10)
        full_seq = db_seq[pro]
        seq = elems[tags["Peptide"]]
        pro = pro + "\t" + disc
        #locate
        loc = full_seq.find(seq)
        uAA = '-'
        dAA = '-'
        if loc > 0: uAA = full_seq[loc-1]
        if loc + len(seq) < len(full_seq): dAA = full_seq[loc+len(seq)]
        #print uAA +'.'+seq+'.'+dAA
        #assert(0)

        spects = elems[tags["Spectrum"]].split('.')
        raw_fn = spects[0]
        scan_id = spects[1]
        charge = int(spects[3])
        score = float(elems[tags["Expectation"]])

        pep_mass = float(elems[tags["Calculated M/Z"]])
        mass = float(elems[tags["Observed M/Z"]]) * charge

        ts = raw_fn.split('_')
        com_fn = ""
        for t in ts[:-1]:
            com_fn = com_fn + t + "_"
        com_fn = com_fn[:-1]
        fn_ndx = ts[-1]

        rt = float(elems[tags["Retention"]])
        #ori_dM = float(elems[tags["Original Delta Mass"]])
        #adj_dM = float(elems[tags["Adjusted Delta Mass"]])

        current_label = None
        nterm_marker = "-"
        markers = []

        #debug
        print pro, seq, loc
        print pep_mass, uAA, dAA, com_fn, fn_ndx, rt, charge

        #check mod
        mods = elems[tags["Assigned Modifications"]]
        if len(mods)>1:
            mods = mods.split(',')
        else:
            mods = []
        for assigned_mod in mods:
            assigned_mod = assigned_mod.strip()
            tmps = assigned_mod[:-1].split('(')
            posA = tmps[0].strip()
            print assigned_mod, posA
            if posA != "N-term":
                mod_mass = float(tmps[1])
                ndx = int(posA[:-1])-1
                ideal_mass = mass_norm_AA[seq[ndx]]
                diff_mass = mod_mass #differ from xml file
                #check label mod
                for m in mod_list:
                    if m[0] == seq[ndx] and fabs(diff_mass-m[1]) <= mass_tol: #AA
                        current_label = tag
                        if m[2]!="-": markers.append((m[2],ndx))
                #check normal mod
                for m in norm_labels:
                    if m[0] == seq[ndx] and fabs(diff_mass-m[1]) <= mass_tol: #AA
                        if m[2]!="-": markers.append((m[2],ndx))
            else:
                #check nterm
                mod_nterm_mass = 0.0

        print current_label
        #after check
        if current_label == None and len(std_AA)>0:
            for aa in std_AA:
                if aa in seq:
                    current_label = tag
                    break
        #finsh mod

        if current_label != None:
            # to be quantified
            if pro not in ipi_lst.keys():
                ipi_lst[pro] = 1
            else:
                ipi_lst[pro] = ipi_lst[pro]+1
            gene = pro.split()[0].split('|')[1]
            if pro[:7] == "Reverse": gene = "Reverse_" + gene
            all_seq = mark_seq(uAA, seq, dAA, markers)
            if nterm_marker!="-":
                all_seq = all_seq[:2] + nterm_marker + all_seq[2:]
            scan_key = gene + ":" + all_seq + ":" + str(charge) + ":" + fn_ndx
            scan_lst.append(scan_key + " " + com_fn + " " + str(scan_id) + " " + current_label + "\n")
            cross_tab[scan_key].append((mass, scan_id, score))

print "STAT:", len(ipi_lst), len(scan_lst)

#output
def get_symbol(disc):
  gene = ""
  elems = disc.split()
  if "GN=" in disc:
    #get from record
    for elem in elems:
      if elem[:3] == "GN=":
        gene = elem[3:]
  elif len(elems)>0:
    gene = elems[0]
  else:
    gene = "NULL"
    print disc
  return gene

#save ipi_name.table
ipiout = open("ipi_name.table", 'w')
ipiout.write("name\n")

out_ipi_lst = [ p for p in ipi_lst.keys() if ipi_lst[p]>=1 ]

out_uni_lst = []
for pro in out_ipi_lst:
    #skip rev decoys ?
    #elems = pro.split('|')
    #ipi = elems[1] #uniprot actually
    #disc = elems[2]
    #tmp = disc.split()
    #gene = tmp[0]
    #disc = disc[len(gene)+1:]
    #symbol = get_symbol(disc)
    elems = pro.split('\t')
    disc = elems[1]
    ipi = elems[0].split('|')[1]
    out_uni_lst.append(ipi)
    if elems[0][:7] == "Reverse":
        ipi = "Reverse_" + ipi
    #disc = disc.split("[REVIEWED]")[0]
    #label = disc.find("[NOT_REVIEWED]")
    label = elems[0].split('|')[2]
    symbol = elems[0].split('|')[2]
    #if label>0:
    #     disc = disc[:label]
    #elems = disc.split()
    #disc = ""
    #for elem in elems:
    #  if "=" in elem: break
    #  disc = disc + " " + elem
    #disc = disc[:min(60,len(disc))]
    #disc = disc.replace("'", "")
    #ipiout.write(pro)
    #ipiout.write("\n")
    ipiout.write(ipi)
    ipiout.write("\t")
    ipiout.write(symbol)
    ipiout.write(" ")
    ipiout.write(disc)
    ipiout.write("\n")
ipiout.close()

#save all_scan.table
scanout = open("all_scan.table", 'w')
scanout.write("key run scan HL\n")
for scan in scan_lst:
    uni = scan.split(":")[0]
    if uni in out_uni_lst: scanout.write(scan)
scanout.close()

#save cross_scan.table
crossout = open("cross_scan.table", 'w')
crossout.write( "key mass %s\n" % (com_fn) )
for scan_core in cross_tab.keys():
    uni = scan_core.split(":")[0]
    if uni not in out_uni_lst: continue
    rank = sorted( cross_tab[scan_core], key=lambda x: x[2] )
    neutral_mass = rank[-1][0]
    id_scan = rank[-1][1]
    crossout.write(scan_core+" "+str(neutral_mass)+" "+str(id_scan)+"\n")
crossout.close()
