#!/usr/bin/env python
import xml.etree.ElementTree as ET
from math import fabs,log
from collections import defaultdict
import sys, glob

ns = {'pep': 'http://regis-web.systemsbiology.net/pepXML'}
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

# Notes: nt refers to nterm of protein, n refers to nterm of any peptide

# Usage:
# pepXMLs2CimageTab.py name sc_cut L-mods H-mods other-mods light/ heavy/
# sc_cut -> 0.05

L_labels = []
H_labels = []
norm_labels = []
#supports multiple nterm modifications, but only one at a time
norm_nterm_labels = []

ipi_lst = defaultdict(None)
scan_lst = []
cross_tab = defaultdict(list)

sample_name = sys.argv[1]
sc_cut = float(sys.argv[2])

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
  #print markers
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

    fn_lst = glob.glob( pepxml_dir[tag]+"/"+sample_name+"*.pepXML" )
    print "Loading ..."
    print fn_lst
    for fn in fn_lst: #for each fraction
        raw_fn = fn.split('/')[-1]
        elems = raw_fn.split('_')
        com_fn = ""
        for elem in elems[:-1]:
            com_fn = com_fn + elem + "_"
        com_fn = com_fn[:-1]
        fn_ndx = elems[-1].split(".")[0]

        try:
            tree = ET.parse(fn)
        except:
            print "Open pepXML file", fn, "failed!"
            sys.exit()
        root = tree.getroot()
        for summary in root.findall('pep:msms_run_summary', ns):
            for query in summary.findall('pep:spectrum_query', ns):
                scan_id = int(query.attrib['start_scan'])
                if scan_id == 0: #for fixed mgf
                    scan_id = int(query.attrib['spectrum'].split('|')[-1])
                mass = float(query.attrib['precursor_neutral_mass'])
                rt = 0.0
                if 'retention_time_sec' in query.attrib.keys():
                    rt = float(query.attrib['retention_time_sec'])
                charge = int(query.attrib['assumed_charge'])
                for result in query.findall('pep:search_result', ns):
                    for hit in result.findall('pep:search_hit', ns):
                        #core code
                        score_map = {}
                        for sc in hit.findall('pep:search_score', ns ):
                            score_map[sc.attrib['name']]= float(sc.attrib['value'])
                        score = -log(score_map["expect"])
                        if score < sc_cut: continue

                        #pass filter
                        pro = hit.attrib['protein'].strip()
                        seq = hit.attrib['peptide']
                        pep_mass = float( hit.attrib['calc_neutral_pep_mass'] )
                        uAA = hit.attrib['peptide_prev_aa']
                        dAA = hit.attrib['peptide_next_aa']
                        num_tryp_ter = int(hit.attrib['num_tol_term'])
                        rank = int( hit.attrib['hit_rank'] )
                        tot_num_ions = int( hit.attrib['tot_num_ions'] )
                        num_matched_ions = int( hit.attrib['num_matched_ions'] )
                        delta_mass = float(hit.attrib['massdiff'])
                        num_miss_clv = int( hit.attrib['num_missed_cleavages'] )

                        #TBD: for SILAC all light is needed, but for ABPP only labeled light is required
                        current_label = None
                        #current_label = "light"
                        nterm_marker = "-"
                        markers = []
                        for mod_info in hit.findall('pep:modification_info', ns):
                            #check others
                            for mod in mod_info.findall('pep:mod_aminoacid_mass', ns):
                                ndx = int(mod.attrib['position'])-1
                                mod_mass = float(mod.attrib['mass']) #modified mass in pepxml
                                ideal_mass = mass_norm_AA[seq[ndx]]  #ideal AA mass
                                diff_mass = mod_mass - ideal_mass    #mass of modification
                                #check label mod
                                for m in mod_list:
                                    if m[0] == seq[ndx] and fabs(diff_mass-m[1]) <= mass_tol: #AA
                                        current_label = tag
                                        if m[2]!="-": markers.append((m[2],ndx))
                                #check normal mod
                                for m in norm_labels:
                                    if m[0] == seq[ndx] and fabs(diff_mass-m[1]) <= mass_tol: #AA
                                        if m[2]!="-": markers.append((m[2],ndx))
                            #check nterm
                            mod_nterm_mass = 0.0
                            if 'mod_nterm_mass' in mod_info.attrib.keys():
                                mod_nterm_mass = float(mod_info.attrib['mod_nterm_mass'])
                                if len(nterm_labels)>0:
                                    for nl in nterm_labels:
                                        if nl[0] == "nt" and uAA!="-": continue
                                        if fabs(mod_nterm_mass - mass_norm_AA["n"] - nl[1]) < mass_tol:
                                            current_label = tag
                                            nterm_marker = nl[2]
                                            #mark this nterm, only one mod allowed
                                            break
                                if len(norm_nterm_labels)>0:
                                    for nl in norm_nterm_labels:
                                        #protein N-term
                                        if nl[0] == "nt" and uAA!="-": continue
                                        #peptide N-term
                                        #print mod_nterm_mass, mass_norm_AA["n"], nl[1]
                                        #mod_nterm_mass
                                        #mass_norm_AA["n"]: natrual nterm massdiff, i.e., 1.0078
                                        #nl[1]: nterm label mass
                                        if fabs(mod_nterm_mass - mass_norm_AA["n"] - nl[1]) < mass_tol:
                                            nterm_marker = nl[2]
                                            #mark this nterm, only one mod allowed
                                            break

                        #after mod checking
                        if current_label == None and len(std_AA)>0:
                            #scan static mod like K,R of SILAC
                            #currentlly, nterm-only labeling are not supported
                            for aa in std_AA:
                                if aa in seq:
                                    current_label = tag
                                    break
                        #finish mod

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
                            #if scan_key not in cross_tab.keys():
                            #    cross_tab[scan_key] = []
                            #if fabs(mass-156.1263)<mass_tol:
                            #  print "DB:", mass, scan_id, score
                            cross_tab[scan_key].append((mass, scan_id, score))

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
    elems = pro.split('|')
    ipi = elems[1] #uniprot actually
    out_uni_lst.append(ipi)
    disc = elems[2]
    tmp = disc.split()
    gene = tmp[0]
    disc = disc[len(gene)+1:]
    symbol = get_symbol(disc)
    if pro[:7] == "Reverse":
        ipi = "Reverse_" + ipi
    disc = disc.split("[REVIEWED]")[0]
    label = disc.find("[NOT_REVIEWED]")
    if label>0:
         disc = disc[:label]
    elems = disc.split()
    disc = ""
    for elem in elems:
      if "=" in elem: break
      disc = disc + " " + elem
    disc = disc[:min(60,len(disc))]
    disc = disc.replace("'", "")
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

