#!/usr/bin/env python

import sys, glob
import argparse
from math import fabs
import numpy as np
from collections import defaultdict
from Bio import SeqIO

parse = argparse.ArgumentParser(prog = 'csv2CimageTab', add_help=False)
parse.add_argument('-c', '--config', help='config file')
parse.add_argument('-t', '--title', help='prefix of raw data')
parse.add_argument('-l', '--light', help='modifications of L label')
parse.add_argument('-m', '--medium', default="", help='modifications of M label')
parse.add_argument('-h', '--heavy', help='modifications of H label')
parse.add_argument('-L', '--light_file', help='csv of L')
parse.add_argument('-M', '--medium_file', help='csv of M')
parse.add_argument('-H', '--heavy_file', help='csv of H')
parse.add_argument('-n', '--normal', help="modifications of ALL")
parse.add_argument('-f', '--fasta', help="fasta db")
parse.add_argument('-r', '--rev', default="Reverse_", help="discard rev seq")
parse.add_argument('-j', '--jump', type=int, default=1, help="skip N lines of csv")
parse.add_argument('-s', '--space', default='\t', help="sep col")

args = parse.parse_args()
#parse.print_help()

mass_tol = 0.01

mass_norm_AA = {
"G" : 57.02146,
"A" : 71.03711,
"S" : 87.03203,
"P" : 97.05276,
"V" : 99.06841,
"T" :101.04768,
"C" :160.03064,
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
M_labels = []
H_labels = []
norm_labels = []
norm_nterm_labels = []

ipi_lst = defaultdict(None)
scan_lst = []
cross_tab = defaultdict(list)

sample_name = args.title
db_fn = args.fasta

db_seq = {}
for record in SeqIO.parse(db_fn, "fasta"):
    if record.id[:8] == args.rev: continue
    pro = record.id.split()[0].split('|')[1]
    db_seq[pro] = record.seq

#light
for mod in args.light.split('|'):
    elems = mod.split(':')
    if len(elems) == 3:
        AA = elems[0]
        mod_mass = float(elems[1])
        tag = elems[2]
        L_labels.append((AA, mod_mass, tag))
#medium
for mod in args.medium.split('|'):
    elems = mod.split(':')
    if len(elems) == 3:
        AA = elems[0]
        mod_mass = float(elems[1])
        tag = elems[2]
        M_labels.append((AA, mod_mass, tag))
#heavy
for mod in args.heavy.split('|'):
    elems = mod.split(':')
    if len(elems) == 3:
        AA = elems[0]
        mod_mass = float(elems[1])
        tag = elems[2]
        H_labels.append((AA, mod_mass, tag))
#normal
for mod in args.normal.split('|'):
    elems = mod.split(':')
    if len(elems) == 3:
        AA = elems[0]
        mod_mass = float(elems[1])
        tag = elems[2]
        if AA[0] == "n":
            norm_nterm_labels.append((AA, mod_mass, tag))
        else:
            norm_labels.append((AA, mod_mass, tag))

print("quant labels")
print(L_labels)
print(M_labels)
print(H_labels)
print("normal labels")
print(norm_labels)
print(norm_nterm_labels)

csv_dir = {}
csv_dir["light"] = args.light_file
csv_dir["medium"] = args.medium_file
csv_dir["heavy"] = args.heavy_file

def mark_seq(uAA, seq, dAA, markers):
    ans = ""
    if len(markers)==0:
        ans = uAA + "." + seq + "." + dAA
    else:
        #all_seq = uAA + "." + seq[:ndx+1] + "*" + seq[ndx+1:] + "." + dAA
        markers = sorted(markers, key=lambda x:x[1])
        ans = uAA + "."
        prev = 0
        for marker, ndx in markers:
            ans = ans + seq[prev:ndx+1] + marker
            prev = ndx + 1
            #print(ans, prev)
        ans = ans + seq[prev:] + "." + dAA
    return ans

print("Parsing CSV files...")
for tag in [ "light", "medium", "heavy" ]:
    if csv_dir[tag] is None: continue
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
    elif tag =="medium":
        for m in M_labels:
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
        print("Warning tag error!")
        break
    print("## ", tag)
    print("mod list")
    print(mod_list)
    print("nter labels")
    print(nterm_labels)
    print("others")
    print(std_AA)

    print("Loading ...")
    fn = csv_dir[tag]
    print(fn)

    lines = open(fn, 'r').readlines()
    elems = lines[0].strip().split(args.space)
    tags = {}
    for n, e in enumerate(elems):
        tags[e] = n
    #print(tags)
    for l in lines[args.jump:]:
        elems = l.strip().split(args.space)
        #print(elems)

        pro = elems[tags['"Master Protein Accessions"']]
        pro = pro.strip('"').split(';')[0]
        if len(pro) == 0:
            pro = elems[tags['"Protein Accessions"']]
            pro = pro.strip('"').split(';')[0]
        disc = ""
        mod_seq = elems[tags['"Annotated Sequence"']]
        mod_seq = mod_seq.strip('"').split(".")
        #print(pro, mod_seq)
        uAA = mod_seq[0][-2]
        dAA = mod_seq[2][1]
        seq = mod_seq[1].upper()
        try:
            full_seq = db_seq[pro]
        except:
            print("Protein not found!", pro, l)
            continue
        loc = full_seq.find(seq)
        #print(pro, uAA+'.'+seq+'.'+dAA, full_seq, loc)

        raw_fn =  elems[tags['"Spectrum File"']].strip('"').split('.')[0]
        scan_id = elems[tags['"First Scan"']].strip('"')
        charge = int(elems[tags['"Charge"']].strip('"'))
        score = float(elems[tags['"Percolator q-Value"']].strip('"'))

        #"m/z [Da]"	"MH+ [Da]"	"Theo. MH+ [Da]"
        #pep_mass = float(elems[tags["Calculated M/Z"]])
        mass = float(elems[tags['"MH+ [Da]"']].strip('"'))

        ts = raw_fn.split('_')
        com_fn = ""
        for t in ts[:-1]:
            com_fn = com_fn + t + "_"
        com_fn = com_fn[:-1]
        fn_ndx = ts[-1]

        rt = float(elems[tags['"RT [min]"']].strip('"'))

        current_label = None
        nterm_marker = "-"
        markers = []
        #check mod
        mods = elems[tags['"Modifications"']]
        if len(mods)>2:
            mods = [m.strip('"') for m in mods.split(';')]
            #print(mods)
        else:
            mods = []

        for assigned_mod in mods:
            if assigned_mod[:6] != "N-Term":
                tmps = assigned_mod[:-1].strip()
                loc = tmps.find('(')
                posA = tmps[:loc]
                modn = tmps[loc+1:]
            else:
                tmps = assigned_mod[:-1].strip()
                loc = tmps.find(')(')
                posA = "N-Term"
                modn = tmps[loc+2:]
            #skip N-term by now ...
            if posA != "N-Term":
                ndx = int(posA[1:])-1
                for m in mod_list:
                    if m[0] == modn:
                        current_label = tag
                        if m[2]!="-": markers.append((m[2], ndx))
                for m in norm_labels:
                    if m[0] == modn:
                        if m[2]!="-": markers.append((m[2], ndx))
            else:
                for m in nterm_labels:
                    current_label = tag
                    nterm_marker = m[2]

        if current_label == None and len(std_AA)>0:
            for aa in std_AA:
                if aa in seq:
                    current_label = tag
                    break

        if current_label != None:
            if pro not in ipi_lst.keys():
                ipi_lst[pro] = 1
            else:
                ipi_lst[pro] = ipi_lst[pro]+1
            gene = pro
            if pro[:7] == "Reverse": gene = "Reverse_" + gene
            print("seq:", seq)
            print("marks:", markers)
            all_seq = mark_seq(uAA, seq, dAA, markers)
            print("result:", all_seq)
            if nterm_marker!="-":
                all_seq = all_seq[:2] + nterm_marker + all_seq[2:]

            scan_key = gene + ":" + all_seq + ":" + str(charge) + ":" + fn_ndx
            print(scan_key)
            scan_lst.append(scan_key + " " + com_fn + " " + str(scan_id) + " " + current_label + "\n")
            cross_tab[scan_key].append((mass, scan_id, score))
    print("Done!")
    print("scan_list:", len(scan_lst))
    print("cross_table:", len(cross_tab.keys()))

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
    print(pro)
    elems = pro.split('\t')
    disc = "NO INFO"
    ipi = pro
    out_uni_lst.append(ipi)
    if elems[0][:7] == "Reverse":
        ipi = "Reverse_" + ipi
    #disc = disc.split("[REVIEWED]")[0]
    #label = disc.find("[NOT_REVIEWED]")
    label = 0
    symbol = "NONE"
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
print(out_uni_lst)
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

