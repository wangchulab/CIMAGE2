#!/usr/bin/env python

# calculate peptide mass, either monoisotopic or AA-based avergaed mass
# mono -- absolute monoisotopic mass
# mz -- monoisotopic mass with C13 adjustment
# avg -- use average aa mass

import sys
import os
from sys import argv

if len(argv) != 3:
    print 'Usage: %s <peptide_sequence> <mono|mz|avg>'%argv[0]
    print 'calculate peptide mass based on SEQUEST'
    print 'mono -- monoisotopic mass'
    print 'mz -- monoisotopic mass + C13 adjustment'
    print 'avg -- sum of averaged aa mass'
    sys.exit()

from peptideParams import *

if argv[2] == 'mono' or argv[2] == 'mz':
    weightmap = AA_MonoMassMap
    dHplus = MonoHplus
    dH = MonoH
    dOH = MonoOH
elif argv[2] == 'avg':
    weightmap = AA_AvgMassMap
    dHplus = AvgHplus
    dH = AvgH
    dOH = AvgOH
else:
    print 'Unrecoganized option: %s'%argv[2]
    sys.exit()

if os.path.exists(argv[1]):
    sequences = open(argv[1]).read().split('\n')
    if sequences[-1] == '':
        sequences=sequences[:-1]
else:
    sequences=[argv[1],]

for pepseq in sequences:

    pepseq = pepseq.upper()
    mass = 0.0
    ncarbon = 0
    for c in pepseq:
        if c in weightmap:
            mass = mass + weightmap[c]
            ncarbon = ncarbon + AA_NumCarbonMap[c]
        else:
            print "unrecognized sequence", c
            sys.exit()

    mass = mass + dH + dOH
    C13_adjust = 0
    if (argv[2] == 'mz'):
        C13_adjust = int(ncarbon*percC13)
    mass = mass + deltaC12C13*C13_adjust
    print mass
