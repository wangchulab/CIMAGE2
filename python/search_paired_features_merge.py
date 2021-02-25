#!/usr/bin/env python
import sys
from rpy2.robjects.packages import importr
import rpy2.robjects as rob
from scipy.stats import pearsonr
import numpy as np
from math import log, fabs

tol = 0.001
dM0 = 1.003355
rt_gap = 0.2 #merge
max_Np = 4
max_R2 = 0.6
charges = [ 2, 3, 4, 5, 6 ]

class paired_feats:
    #for each pair, we save mz, mzApex, rt range, charge and ratio( fit chrom... )
    def __init__(self):
        self.rtStart = 0.0
        self.rtEnd = 0.0
        self.mz = 0.0
        self.mzLApex = 0.0
        self.mzHApex = 0.0
        self.rtApex = 0.0
        self.intApex = 0.0
        self.intR = 0.0
        self.intSum = 0.0
        self.corr = 0.0
        self.charge = 0
        self.Np = 0

    def output(self):
        return "MZ: %6.4f, charge: %d, L: %6.4f, H: %6.4f, RT: %6.4f, R: %6.4f, DB:( %3.2f %3.2f %d )" \
        % ( self.mz, self.charge, self.mzLApex, self.mzHApex, self.rtApex*60.0, self.intR, self.corr, self.intApex, self.Np )

def check_chromatograms_corr( intensL, intensH ):
    Ls = [0.0]
    Hs = [0.0]
    for l, h in zip(intensL[1], intensH[1]):
        if l+h > 0.01:
            Ls.append(l)
            Hs.append(h)

    Np = len(Ls)
    if Np > max_Np:
        corr, _ = pearsonr(Ls, Hs)
        ratio = np.sum(Ls)+0.01 / np.sum(Hs)+0.01
        if ratio > 100.0: ratio = 100.0
        if ratio < 0.01: ratio = 0.01
    else:
        corr = 0.0
        ratio = 0.0
    return corr, ratio, Np

def merge_pairs( closed_pairs ):
    merged_pairs = []
    rt1_s = -999.9
    rt1_e = -999.9
    for p in sorted(closed_pairs, key=lambda x: x.rtApex):
        rt0_s = rt1_s
        rt0_e = rt1_e
        rt1_s = p.rtStart
        rt1_e = p.rtEnd
        if rt0_s >= rt1_e:
            delta_rt = rt0_s - rt1_e
        elif rt1_s >= rt0_e:
            delta_rt = rt1_s - rt0_e
        elif rt0_s >= rt1_s and rt0_s <= rt1_e:
            delta_rt = 0.0
        elif rt1_s >= rt0_s and rt1_s <= rt0_e:
            delta_rt = 0.0

        #print delta_rt
        if delta_rt < rt_gap:
            #print("merge!")
            merged_pairs[-1].rtStart = min(rt0_s, rt1_s)
            merged_pairs[-1].rtEnd = max(rt1_e, rt1_e)
            I1 = merged_pairs[-1].intSum
            r1 = merged_pairs[-1].intR
            I2 = p.intSum
            r2 = p.intR
            merged_pairs[-1].intSum = (I1+I2)
            merged_pairs[-1].intR = (I1*r1+I2*r2)/(I1+I2)
            if merged_pairs[-1].intApex < p.intApex:
                merged_pairs[-1].intApex = p.intApex
                merged_pairs[-1].rtApex = p.rtApex
        else:
            #print(delta_rt)
            merged_pairs.append(p)
    return merged_pairs

def main():
    xcms = importr('xcms')

    fn_feat = sys.argv[1]
    fn_mzML = fn_feat.replace("features.tsv", "mzML")
    shift0 = float(sys.argv[2])

    #load mzML file
    fl = xcms.xcmsRaw(fn_mzML, profstep=0, includeMSn=False)

    lines = open(fn_feat, 'r').readlines()
    tags = { k:n for n,k in enumerate(lines[0].strip().split('\t')) }

    pairs = {}
    for c in charges:
        pairs[c] = []

    for n, l in enumerate(lines[1:]):
        es = l.strip().split('\t')

        charge = int(es[tags['charge']])
        if charge not in charges: continue

        rtStart = float(es[tags['rtStart']])
        rtEnd = float(es[tags['rtEnd']])
        mz = float(es[tags['mz']])
        mzApex = float(es[tags['mostAbundantMz']])
        rtApex = float(es[tags['rtApex']])
        intApex = float(es[tags['intensityApex']])
        intSum = float(es[tags['intensitySum']])

        mz0 = mzApex
        mz_p = mz0 - shift0 / charge
        mz_q = mz0 + shift0 / charge
        rt_range = rob.FloatVector([rtStart*60.0, rtEnd*60.0])

        mz_p_range = rob.FloatVector([mz_p-tol, mz_p+tol])
        EIC_p = xcms.rawEIC(fl, mz_p_range, rt_range)
        scan_p, intens_p = EIC_p.items()

        mz0_range = rob.FloatVector([mz0-tol, mz0+tol])
        EIC0 = xcms.rawEIC(fl, mz0_range, rt_range)
        scan0, intens0 = EIC0.items()

        mz_q_range = rob.FloatVector([mz_q-tol, mz_q+tol])
        EIC_q = xcms.rawEIC(fl, mz_q_range, rt_range)
        scan_q, intens_q = EIC_q.items()

        corr1, r1, Np1 = check_chromatograms_corr( intens_p, intens0 )
        if corr1 > max_R2:
            #print(corr1, r1)
            p = paired_feats()
            p.mz = mz
            p.mzLApex = mz_p
            p.mzHApex = mz0
            p.rtStart = rtStart
            p.rtEnd = rtEnd
            p.rtApex = rtApex
            p.intApex = intApex
            p.intR = r1
            p.intSum = intSum
            p.charge = charge
            p.corr = corr1
            p.Np = Np1
            pairs[charge].append(p)
            #print(p.output())

        corr2, r2, Np2 = check_chromatograms_corr( intens0, intens_q )
        if corr2 > max_R2:
            #print(corr2, r2)
            p = paired_feats()
            p.mz = mz
            p.mzLApex = mz0
            p.mzHApex = mz_q
            p.rtStart = rtStart
            p.rtEnd = rtEnd
            p.rtApex = rtApex
            p.intApex = intApex
            p.intR = r2
            p.intSum = intSum
            p.charge = charge
            p.corr = corr2
            p.Np = Np2
            pairs[charge].append(p)
            #print(p.output())

    #merge all pair
    for c in charges:
        #print("charge:", c, "Npair:", len(pairs[c]))
        closed_pairs = []
        sorted_pairs = sorted(pairs[c], key=lambda x: x.mzLApex)
        nlast = len(sorted_pairs)
        if nlast == 0: continue
        p0 = 0
        p1 = 1
        closed_pairs.append(sorted_pairs[p0])

        while p1 < nlast:
            mz0 = sorted_pairs[p0].mzLApex
            mz1 = sorted_pairs[p1].mzLApex
            mz_tol = tol
            if mz1 - mz0 > mz_tol:
                merged_pairs = merge_pairs(closed_pairs)
                for p in merged_pairs:
                    print( p.output() )

                closed_pairs = []
                p0 = p1
                p1 = p0 + 1
                closed_pairs.append(sorted_pairs[p0])
            else:
                closed_pairs.append(sorted_pairs[p1])
                p1 = p1 + 1

        merged_pairs = merge_pairs(closed_pairs)
        for p in merged_pairs:
            print( p.output() )

#MAIN
if __name__ == "__main__":
    main()

