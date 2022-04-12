import sys
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
print glc_DB

glc_label = {}
glc_label["light"] = sys.argv[1]
glc_label["heavy"] = sys.argv[2]

glc_dict = defaultdict(list)
#for each glc
ipi_lst = defaultdict(list)
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
        pro = es[tags["Proteins"]]
        gene = es[tags["Genes"]]
        name = es[tags["RawName"]]
        charge = int(es[tags["Charge"]])
        scan = es[tags["Scan"]]
        mod = es[tags["Mod"]]
        label = es[tags["Glycan(H,N,A,F," + glc_label[tag] + ")"]]
        label = "".join(label.split(' '))
        glc_dict[label].append( scan )
        pep_mass = float(es[tags["PeptideMH"]])
        glc_mass = float(es[tags["GlyMass"]])
        tol_fdr = float(es[tags["TotalFDR"]])


        #print pro, charge, scan, mod, label
        #if len(mod)>5:
        #    break

for k in glc_dict.keys():
    print "G:", k, len(glc_dict[k])
