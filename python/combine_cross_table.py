import sys
from collections import defaultdict
import numpy as np

name = []
nf = len(sys.argv) - 1
cross = defaultdict(list)
for n, fn in enumerate(sys.argv[1:]):
    with open(fn, 'r') as fp:
        lines = fp.readlines()
        name.append(lines[0].strip().split()[-1])
        for l in lines[1:]:
            es = l.strip().split()
            key = es[0]
            mass = float(es[1])
            scan = int(es[2])
            cross[key].append((n, mass, scan))

print "key", "mass", " ".join(name)
for k in cross.keys():
    mass = []
    index = [0]*nf
    for n, m, s in cross[k]:
        mass.append(m)
        index[n] = s
    mass = np.median(mass)
    print k, "%6.4f" % mass, " ".join([str(s) for s in index])

