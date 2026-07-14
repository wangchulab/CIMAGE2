#!/usr/bin/env python3
"""Pure-Python port of search_paired_features_merge_v2.py (no rpy2/xcms).
Builds EICs from the .ms1 file instead of xcms.rawEIC(mzML).

Per Dinosaur feature it extracts 3 chromatograms over the feature RT range
(mz-shift/z, mz, mz+shift/z), computes light/heavy Pearson corr + integrated
ratio, keeps a pair if corr passes an adaptive threshold, then merges pairs
overlapping in RT. Emits the same fields as the R script.

We cache (corr, ratio, Np, mz's, rt's) per feature-direction so downstream
threshold / ratio-cut sweeps are cheap.
"""
import sys, bisect, math, csv
import numpy as np

def pearson(x, y):
    x = np.asarray(x); y = np.asarray(y)
    xm = x - x.mean(); ym = y - y.mean()
    d = math.sqrt((xm*xm).sum() * (ym*ym).sum())
    return (xm*ym).sum() / d if d > 0 else 0.0

TOL   = 0.001       # EIC m/z half-window (Da), as in the R script
DM0   = 1.003355
CHARGES = {2, 3, 4, 5, 6}


def load_ms1(path, cache=True):
    """-> (RT[min], MZ[list of np arrays], IN[list of np arrays]). Cached to .pkl."""
    import os, pickle
    pkl = path + '.xiccache.pkl'
    if cache and os.path.exists(pkl) and os.path.getmtime(pkl) >= os.path.getmtime(path):
        with open(pkl, 'rb') as f:
            return pickle.load(f)
    rts = []; counts = []; mz_s = []; in_s = []
    cur_rt = 0.0; cnt = 0; have = False
    with open(path) as f:
        for line in f:
            c = line[0]
            if c == 'S':
                if have:
                    rts.append(cur_rt); counts.append(cnt)
                cnt = 0; cur_rt = 0.0; have = True
            elif c == 'I':
                if line[2] == 'R':  # RTime
                    cur_rt = float(line.rsplit('\t', 1)[1])
            elif '0' <= c <= '9':
                sp = line.split()
                mz_s.append(sp[0]); in_s.append(sp[1]); cnt += 1
        if have:
            rts.append(cur_rt); counts.append(cnt)
    allmz = np.asarray(mz_s, dtype=float)
    allin = np.asarray(in_s, dtype=float)
    MZ = []; IN = []; off = 0
    for n in counts:
        MZ.append(allmz[off:off+n]); IN.append(allin[off:off+n]); off += n
    order = sorted(range(len(rts)), key=lambda i: rts[i])
    RT = [rts[i] for i in order]; MZ = [MZ[i] for i in order]; IN = [IN[i] for i in order]
    if cache:
        with open(pkl, 'wb') as f:
            pickle.dump((RT, MZ, IN), f)
    return RT, MZ, IN


def load_mzml(path, cache=True):
    """Read MS1 scans directly from mzML (no rpy2/xcms). -> (RT[min], MZ, IN)."""
    import os, pickle
    pkl = path + '.xiccache.pkl'
    if cache and os.path.exists(pkl) and os.path.getmtime(pkl) >= os.path.getmtime(path):
        with open(pkl, 'rb') as f:
            return pickle.load(f)
    from pyteomics import mzml
    rts = []; MZ = []; IN = []
    for s in mzml.read(path):
        if s.get('ms level') != 1:
            continue
        rt = float(s['scanList']['scan'][0]['scan start time'])  # minutes
        mz = np.asarray(s['m/z array'], dtype=float)
        it = np.asarray(s['intensity array'], dtype=float)
        rts.append(rt); MZ.append(mz); IN.append(it)
    order = sorted(range(len(rts)), key=lambda i: rts[i])
    RT = [rts[i] for i in order]; MZ = [MZ[i] for i in order]; IN = [IN[i] for i in order]
    if cache:
        with open(pkl, 'wb') as f:
            pickle.dump((RT, MZ, IN), f)
    return RT, MZ, IN


def load_scans(path):
    """Dispatch: .mzML -> load_mzml, else .ms1 loader."""
    return load_mzml(path) if path.lower().endswith('mzml') else load_ms1(path)


def eic(RT, MZ, IN, rt_lo, rt_hi, mz, tol):
    """max intensity within [mz-tol,mz+tol] per scan in [rt_lo,rt_hi]."""
    lo = bisect.bisect_left(RT, rt_lo); hi = bisect.bisect_right(RT, rt_hi)
    out = []
    for i in range(lo, hi):
        a = MZ[i]
        if len(a) == 0:
            out.append(0.0); continue
        j0 = bisect.bisect_left(a, mz - tol); j1 = bisect.bisect_right(a, mz + tol)
        out.append(float(IN[i][j0:j1].max()) if j1 > j0 else 0.0)
    return out


def corr_ratio(Ls, Hs, max_Np):
    xs = []; ys = []
    for l, h in zip(Ls, Hs):
        if l * h > 0.01:
            xs.append(l); ys.append(h)
    Np = len(xs) + 1  # R script starts lists with a 0.0 seed
    if Np > max_Np:
        xs = [0.0] + xs; ys = [0.0] + ys
        corr = pearson(xs, ys)
        ratio = sum(xs) / sum(ys) if sum(ys) > 0 else 0
        ratio = min(100.0, max(0.01, ratio))
    else:
        corr, ratio = 0.0, 0.0
    return corr, ratio, Np


def candidates(features_path, ms1_path, shift0, tol=TOL, max_Np=5):
    """Return list of candidate pairs (before threshold): dicts with
    charge, L, H, rtStart, rtEnd, rtApex, intApex, intSum, corr, ratio, Np, dir."""
    RT, MZ, IN = load_scans(ms1_path)
    out = []
    with open(features_path) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            z = int(row['charge'])
            if z not in CHARGES:
                continue
            mz = float(row['mz']); mzApex = mz          # v2: use mono mz
            rtS = float(row['rtStart']); rtE = float(row['rtEnd'])
            rtA = float(row['rtApex']); iA = float(row['intensityApex'])
            iS = float(row['intensitySum'])
            sep = shift0 / z
            e0 = eic(RT, MZ, IN, rtS, rtE, mz, tol)
            ep = eic(RT, MZ, IN, rtS, rtE, mz - sep, tol)
            eq = eic(RT, MZ, IN, rtS, rtE, mz + sep, tol)
            # this feature as HEAVY (partner below): corr(ep, e0)
            c1, r1, n1 = corr_ratio(ep, e0, max_Np)
            if n1 > max_Np:
                out.append(dict(charge=z, L=mzApex - sep, H=mzApex, rtStart=rtS,
                                rtEnd=rtE, rtApex=rtA, intApex=iA, intSum=iS,
                                corr=c1, ratio=r1, Np=n1))
            # this feature as LIGHT (partner above): corr(e0, eq)
            c2, r2, n2 = corr_ratio(e0, eq, max_Np)
            if n2 > max_Np:
                out.append(dict(charge=z, L=mzApex, H=mzApex + sep, rtStart=rtS,
                                rtEnd=rtE, rtApex=rtA, intApex=iA, intSum=iS,
                                corr=c2, ratio=r2, Np=n2))
    return out


def adaptive_ok(corr, Np, A=10.0, base=0.95, k=4.0):
    return corr > -Np / A + k / A + base


def merge(pairs, rt_gap=0.2, mz_tol=TOL):
    """merge pairs overlapping in RT (within rt_gap) after grouping by charge+mz."""
    from collections import defaultdict
    bych = defaultdict(list)
    for p in pairs:
        bych[p['charge']].append(p)
    merged = []
    for z, ps in bych.items():
        ps.sort(key=lambda x: x['L'])
        i = 0; n = len(ps)
        while i < n:
            j = i + 1
            grp = [ps[i]]
            while j < n and ps[j]['L'] - ps[j-1]['L'] <= mz_tol:
                grp.append(ps[j]); j += 1
            # within an m/z cluster, merge by RT overlap
            grp.sort(key=lambda x: x['rtApex'])
            cur = None
            for p in grp:
                if cur is None:
                    cur = dict(p); merged.append(cur); continue
                overlap = not (p['rtStart'] > cur['rtEnd'] + rt_gap or
                               cur['rtStart'] > p['rtEnd'] + rt_gap)
                if overlap:
                    I1 = cur['intSum']; I2 = p['intSum']
                    cur['intR'] = (I1*cur['ratio'] + I2*p['ratio']) / (I1 + I2) \
                        if 'intR' not in cur else (cur['intR']*I1 + p['ratio']*I2)/(I1+I2)
                    cur['intSum'] = I1 + I2
                    cur['rtStart'] = min(cur['rtStart'], p['rtStart'])
                    cur['rtEnd'] = max(cur['rtEnd'], p['rtEnd'])
                    if p['intApex'] > cur['intApex']:
                        cur['intApex'] = p['intApex']; cur['rtApex'] = p['rtApex']
                        cur['L'] = p['L']; cur['H'] = p['H']
                else:
                    cur = dict(p); merged.append(cur)
            i = j
    for m in merged:
        m.setdefault('intR', m['ratio'])
    return merged


if __name__ == '__main__':
    feat, ms1, shift = sys.argv[1], sys.argv[2], float(sys.argv[3])
    cands = candidates(feat, ms1, shift)
    kept = [c for c in cands if adaptive_ok(c['corr'], c['Np'])]
    mg = merge(kept)
    print(f'candidates={len(cands)} kept(threshold)={len(kept)} merged={len(mg)}')
