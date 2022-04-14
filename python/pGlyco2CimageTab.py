import sys, os
from math import log
from Bio import SeqIO
from collections import defaultdict

map_atom_mass = {
    "C"  : 12.0,
    "H"  : 1.00783,
    "O"  : 15.99491,
    "N"  : 14.00307,
    "S"  : 31.97207,
    "P"  : 30.97376,
    "13C": 13.003355,
    "15N": 15.000109,
    "2H" : 2.014102
}

glc_DB = {}
dbfn = "../Glycan.txt"
lines = open(dbfn, 'r').readlines()
for l in lines:
    es = l.strip().split()
    name = es[0].split('=')[1]
    atoms = es[-1][:-1].split(')')
    glc_DB[name] = {}
    #print name, atoms
    for atom in atoms:
        atm, n = atom.split('(')
        n = int(n)
        glc_DB[name][atm] = n
#print glc_DB

seq_DB = {}
des_DB = {}
for record in SeqIO.parse("../UP000005640_9606.fasta", "fasta"):
    if record.id[:8] == "Reverse": continue
    pro = record.id
    seq_DB[pro] = record.seq
    des_DB[pro] = record.description
#print seq_DB.keys()

tab_DB = defaultdict(list)
lines = open("../tab.tmpl", 'r').readlines()
for l in lines:
    es = l.strip().split()
    name = es[0]
    markers = []
    for marker in es[1:]:
        atoms = []
        tmp = marker.split(':')
        mark = tmp[0]
        for atom in tmp[1:]:
            a, n = atom.split(',')
            n = int(n)
            atoms.append((a,n))
        markers.append( (mark, atoms) ) 
    tab_DB[name] = markers
#print tab_DB

L_labels = []
H_labels = []
norm_labels = {}
glc_label = {}
#light
for mod in sys.argv[1].split('|'):
    elems = mod.split(':')
    if len(elems) == 3:
        mod_label = elems[0]
        mod_mass = float(elems[1])
        tag = elems[2]
        L_labels.append((mod_label, mod_mass, tag))
        glc_label["light"] = mod_label
#heavy
for mod in sys.argv[2].split('|'):
    elems = mod.split(':')
    if len(elems) == 3:
        mod_label = elems[0]
        mod_mass = float(elems[1])
        tag = elems[2]
        L_labels.append((mod_label, mod_mass, tag))
        glc_label["heavy"] = mod_label
#normal
for mod in sys.argv[3].split('|'):
    elems = mod.split(':')
    if len(elems) == 3:
        mod_label = elems[0]
        mod_mass = float(elems[1])
        tag = elems[2]
        if tag != '-':
            norm_labels[mod_label] = tag

def mark_seq(uAA, seq, dAA, markers):
    ans = ""
    if len(markers)==0:
        ans = uAA + "." + seq + "." + dAA
    else:
        markers = sorted(markers, key=lambda x:x[1])
        ans = uAA + "."
        prev = 0
        for marker, ndx in markers:
            ans = ans + seq[prev:ndx+1] + marker
            prev = ndx + 1
        ans = ans + seq[prev:] + "." + dAA
    return ans

#for each glc
ipi_lst = defaultdict(lambda: defaultdict(None))
scan_lst = defaultdict(list)
cross_tab = defaultdict(lambda: defaultdict(list))

fn = {}
fn['light'] = sys.argv[4]+"/pGlycoDB-GP-FDR-Pro-Quant-Site.txt"
fn['heavy'] = sys.argv[5]+"/pGlycoDB-GP-FDR-Pro-Quant-Site.txt"
for tag in ["light", "heavy"]:
    lines = open(fn[tag], 'r').readlines()
    tags= {}
    for i, e in enumerate(lines[0].strip().split('\t')):
        tags[e] = i
    #print tags
    for l in lines[1:]:
        es = l.strip().split('\t')
        pro = es[tags["Proteins"]].split(';')[0]
        full_seq = seq_DB[pro]
        pep_seq = es[tags["Peptide"]].replace("J", "N")
        loc = full_seq.find(pep_seq)
        uAA = '-'
        if loc>0: uAA = full_seq[loc-1]
        dAA = '-'
        if loc+len(pep_seq)<len(full_seq): dAA = full_seq[loc+len(pep_seq)]
        gene = es[tags["Genes"]]
        ipi = pro.split('|')[1]
        raw_fn = es[tags["RawName"]]
        charge = int(es[tags["Charge"]])
        scan_id = es[tags["Scan"]]
        mods = es[tags["Mod"]]
        glc_site = int(es[tags["GlySite"]])
        current_label = glc_label[tag]
        label = es[tags["Glycan(H,N,A,F," + current_label + ")"]]
        label = "".join(label.split(' '))
        pep_mass = float(es[tags["PeptideMH"]])
        glc_mass = float(es[tags["GlyMass"]])
        tot_fdr = float(es[tags["TotalFDR"]])+1e-20
        tot_mass = pep_mass + glc_mass
        rt = float(es[tags["RT"]])

        ts = raw_fn.split('_')
        com_fn = "_".join(ts[:-1])
        fn_ndx = ts[-1]

        #quick check
        if label[-1] == "0": continue

        #normal
        markers = []
        markers.append(('*', glc_site-1))
        if len(mods)>2:
            mods = [m.strip('"') for m in mods.split(';') if len(m)>0]
        else:
            mods = []
        for assigned_mod in mods:
            tmp = assigned_mod.split(',')
            pos = int(tmp[0])-1
            mod = tmp[1]
            if mod in norm_labels.keys():
                markers.append((norm_labels[mod], pos))

        #gather
        pro = pro + '\t' + gene + '\t' + des_DB[pro]
        if pro not in ipi_lst[label].keys():
            ipi_lst[label][pro] = 1
        else:
            ipi_lst[label][pro] += 1
        all_seq = mark_seq(uAA, pep_seq, dAA, markers)
        scan_key = ipi + ":" + all_seq + ":" + str(charge) + ":" + fn_ndx
        scan_lst[label].append(scan_key + " " + com_fn + " " + str(scan_id) + " " + tag + "\n")
        cross_tab[label][scan_key].append((tot_mass, scan_id, -log(tot_fdr)))
        #print pro, ipi, uAA, pep_seq, dAA, loc, markers, all_seq

def gen_glc_atom( name, glc, tag ):
    #name LPG/HPG
    #glc atom list
    #00001
    return []

def get_str( marker ):
    str_mod = marker[0]
    tmp = []
    for atom, N in marker[1]:
        tmp.append(atom+","+str(N))
    return str_mod+":"+":".join(tmp)

for g in ipi_lst.keys():
    print "G:", g
    # generate dir
    try:
        os.makedirs(g)
    except: pass
    # generate tab.txt
    for name in tab_DB.keys():
        for m in tab_DB[name]:
            print get_str(m)
    break
    # generate cimage.params
    # generate ipi/all_scan/cross
    out_ipi_lst = [ p for p in ipi_lst[g].keys() if ipi_lst[g][p]>=1 ]
    out_uni_lst = [ p.split('\t')[0].split('|')[1] for p in ipi_lst[g].keys() if ipi_lst[g][p]>=1 ]
    with open(g+"/ipi_name.table", 'w') as ipiout:
        ipiout.write("name\n")
        for pro in out_ipi_lst:
            es = pro.split('\t')
            disc = es[2].replace("'", "")
            gene = es[1]
            ipi = es[0].split('|')[1]
            ipiout.write(ipi)
            ipiout.write("\t")
            ipiout.write(gene)
            ipiout.write(" ")
            ipiout.write(disc)
            ipiout.write("\n")
    with open(g+"/all_scan.table", 'w') as scanout:
        scanout.write("key run scan HL\n")
        for scan in scan_lst[g]:
            uni = scan.split(":")[0]
            if uni in out_uni_lst: 
                scanout.write(scan)
    with open(g+"/cross_scan.table", 'w') as crossout:
        crossout.write( "key mass %s\n" % (com_fn) )
        for scan_core in cross_tab[g].keys():
            uni = scan_core.split(":")[0]
            if uni not in out_uni_lst: continue
            rank = sorted( cross_tab[g][scan_core], key=lambda x: x[2] )
            neutral_mass = rank[-1][0]
            id_scan = rank[-1][1]
            crossout.write(scan_core+" "+str(neutral_mass)+" "+str(id_scan)+"\n")
