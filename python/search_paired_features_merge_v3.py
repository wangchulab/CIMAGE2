#!/usr/bin/env python
# v3 = v2 + three data-driven optimizations (validated on YJ7 vs ion.tsv GT):
#
#   (1) max_Np 5 -> 3  [RECALL fix, the big one]
#       v2 requires Np>5 (>=6 co-eluting points), which drops sharp / low-abundance
#       peaks and caps recall at ~94% on several ratios. Np>3 (>=4) recovers recall
#       to >=95% everywhere. The adaptive corr threshold already demands corr>0.95
#       at Np=4, so short traces still need near-perfect co-elution -> few false adds.
#
#   (2) theoretical-ratio cut  [PRECISION, removes wrong-ratio false pairs]
#       a real pair's log2(L/H) sits at log2(nominal); random/self pairs are off.
#       pass nominal L:H ratio as argv[3]; keep only |log2(intR)-log2(nom)|<=RTOL.
#
#   (3) heavy-must-be-a-Dinosaur-feature  [PRECISION, removes self-pairs]
#       since the tag shift ~= 6x 13C, a light + its own 6th isotope mimics a pair
#       (and at high nominal ratio even passes the ratio cut). Dinosaur deisotopes,
#       so the light's isotope is NOT a separate feature: require the partner m/z to
#       match a feature -> self-pairs removed, real feature-feature pairs kept.
#
# Usage: search_paired_features_merge_v3.py <features.tsv> <shift0_Da> <nominal_L:H>
import sys
from rpy2.robjects.packages import importr
import rpy2.robjects as rob
import numpy as np
from math import log, fabs, log2

tol = 0.001
dM0 = 1.003355
rt_gap = 0.2
max_Np = 3          # (1) was 5
A = 10.0
charges = [ 2, 3, 4, 5, 6 ]
RTOL = 0.75         # (2) log2 ratio window
feat_ppm = 10.0     # (3) tolerance for "partner is a feature"

def pearson(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    xm = x - x.mean(); ym = y - y.mean()
    d = ((xm*xm).sum() * (ym*ym).sum()) ** 0.5
    return (xm*ym).sum()/d if d > 0 else 0.0

class paired_feats:
    def __init__(self):
        self.rtStart=self.rtEnd=self.mz=self.mzLApex=self.mzHApex=0.0
        self.rtApex=self.intApex=self.intR=self.intSum=self.corr=0.0
        self.charge=self.Np=0
    def output(self):
        return "MZ: %6.4f, charge: %d, L: %6.4f, H: %6.4f, RT: %6.4f, R: %6.4f, DB:( %3.2f %3.2f %d ), RT( %4.3f %4.3f )" \
        % ( self.mz, self.charge, self.mzLApex, self.mzHApex, self.rtApex*60.0, self.intR, self.corr, self.intApex, self.Np, self.rtStart*60.0, self.rtEnd*60.0 )

def check_chromatograms_corr( intensL, intensH ):
    Ls=[0.0]; Hs=[0.0]
    for l,h in zip(intensL[1], intensH[1]):
        if l*h > 0.01:
            Ls.append(l); Hs.append(h)
    Np=len(Ls)
    if Np > max_Np:
        corr = pearson(Ls, Hs)
        ratio = np.sum(Ls)/np.sum(Hs)
        ratio = min(100.0, max(0.01, ratio))
    else:
        corr, ratio = 0.0, 0.0
    return corr, ratio, Np

def merge_pairs( closed_pairs ):
    merged_pairs=[]; rt1_s=rt1_e=-999.9
    for p in sorted(closed_pairs, key=lambda x: x.rtApex):
        rt0_s, rt0_e = rt1_s, rt1_e
        rt1_s, rt1_e = p.rtStart, p.rtEnd
        if rt0_s >= rt1_e:   delta_rt = rt0_s - rt1_e
        elif rt1_s >= rt0_e: delta_rt = rt1_s - rt0_e
        else:                delta_rt = 0.0
        if merged_pairs and delta_rt < rt_gap:
            merged_pairs[-1].rtStart = min(rt0_s, rt1_s)
            merged_pairs[-1].rtEnd   = max(rt0_e, rt1_e)   # fixed (was max(rt1_e,rt1_e))
            I1=merged_pairs[-1].intSum; r1=merged_pairs[-1].intR
            I2=p.intSum; r2=p.intR
            merged_pairs[-1].intSum = I1+I2
            merged_pairs[-1].intR   = (I1*r1+I2*r2)/(I1+I2)
            if merged_pairs[-1].intApex < p.intApex:
                merged_pairs[-1].intApex=p.intApex; merged_pairs[-1].rtApex=p.rtApex
                merged_pairs[-1].mz=p.mz; merged_pairs[-1].mzLApex=p.mzLApex; merged_pairs[-1].mzHApex=p.mzHApex
        else:
            merged_pairs.append(p)
    return merged_pairs

def main():
    xcms = importr('xcms')
    fn_feat = sys.argv[1]
    fn_mzML = fn_feat.replace("features.tsv", "mzML")
    shift0  = float(sys.argv[2])
    nominal = float(sys.argv[3])          # (2) nominal L:H
    log2nom = log2(nominal)
    fl = xcms.xcmsRaw(fn_mzML, profstep=0, includeMSn=False)

    lines = open(fn_feat).readlines()
    tags = { k:n for n,k in enumerate(lines[0].strip().split('\t')) }

    # (3) build per-charge sorted feature m/z lists for the "partner is a feature" test
    feat_mz = { c:[] for c in charges }
    rows = []
    for l in lines[1:]:
        es = l.strip().split('\t'); c = int(es[tags['charge']])
        if c in charges:
            feat_mz[c].append(float(es[tags['mz']])); rows.append((c, es))
    for c in charges: feat_mz[c].sort()
    import bisect
    def is_feature(c, mz):
        a = feat_mz[c]; t = mz*feat_ppm/1e6
        i = bisect.bisect_left(a, mz-t); return i < len(a) and a[i] <= mz+t

    pairs = { c:[] for c in charges }
    for c, es in rows:
        charge = c
        rtStart=float(es[tags['rtStart']]); rtEnd=float(es[tags['rtEnd']])
        mz=float(es[tags['mz']]); mzApex=mz
        rtApex=float(es[tags['rtApex']]); intApex=float(es[tags['intensityApex']]); intSum=float(es[tags['intensitySum']])
        mz0=mz; mz_p=mz0-shift0/charge; mz_q=mz0+shift0/charge
        rt_range = rob.FloatVector([rtStart*60.0, rtEnd*60.0])
        EIC_p = xcms.rawEIC(fl, rob.FloatVector([mz_p-tol, mz_p+tol]), rt_range); _,intens_p = EIC_p.items()
        EIC0  = xcms.rawEIC(fl, rob.FloatVector([mz0-tol,  mz0+tol]),  rt_range); _,intens0 = EIC0.items()
        EIC_q = xcms.rawEIC(fl, rob.FloatVector([mz_q-tol, mz_q+tol]), rt_range); _,intens_q = EIC_q.items()

        # this feature as HEAVY (partner=mz_p below); keep if partner is a feature
        corr1,r1,Np1 = check_chromatograms_corr(intens_p, intens0)
        if corr1 > -Np1/A + 4.0/A + 0.95 and is_feature(charge, mz_p) \
           and fabs(log2(max(r1,1e-3)) - log2nom) <= RTOL:
            p=paired_feats(); p.mz=mz_p; p.mzLApex=mzApex-shift0/charge; p.mzHApex=mzApex
            p.rtStart=rtStart; p.rtEnd=rtEnd; p.rtApex=rtApex; p.intApex=intApex
            p.intR=r1; p.intSum=intSum; p.charge=charge; p.corr=corr1; p.Np=Np1
            pairs[charge].append(p)

        # this feature as LIGHT (partner=mz_q above)
        corr2,r2,Np2 = check_chromatograms_corr(intens0, intens_q)
        if corr2 > -Np2/A + 4.0/A + 0.95 and is_feature(charge, mz_q) \
           and fabs(log2(max(r2,1e-3)) - log2nom) <= RTOL:
            p=paired_feats(); p.mz=mz0; p.mzLApex=mzApex; p.mzHApex=mzApex+shift0/charge
            p.rtStart=rtStart; p.rtEnd=rtEnd; p.rtApex=rtApex; p.intApex=intApex
            p.intR=r2; p.intSum=intSum; p.charge=charge; p.corr=corr2; p.Np=Np2
            pairs[charge].append(p)

    for c in charges:
        sorted_pairs = sorted(pairs[c], key=lambda x: x.mz)
        nlast=len(sorted_pairs)
        if nlast==0: continue
        closed=[sorted_pairs[0]]; p0=0; p1=1
        while p1<nlast:
            if sorted_pairs[p1].mz - sorted_pairs[p0].mz > tol:
                for p in merge_pairs(closed): print(p.output())
                closed=[sorted_pairs[p1]]; p0=p1; p1=p0+1
            else:
                closed.append(sorted_pairs[p1]); p1+=1
        for p in merge_pairs(closed): print(p.output())

if __name__ == "__main__":
    main()
