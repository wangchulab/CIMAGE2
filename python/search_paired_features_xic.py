#!/usr/bin/env python3
"""Pure-Python v3: pair Dinosaur features via XIC read straight from mzML.

No rpy2 / R / xcms and no separate .ms1 -- mzML is decoded with pyteomics
(base64 + zlib + 64-bit floats). Same output format as the R scripts.

Recipe (v3, validated on YJ7 vs ion.tsv GT):
  (1) max_Np = 3  (>=4 co-eluting points; recovers recall the R-v2 Np>=6 lost)
  (2) theoretical-ratio cut  |log2(intR) - log2(nominal)| <= RTOL  (removes wrong-ratio false pairs)
  (3) heavy must be an independent Dinosaur feature  (removes 6x-13C self-pairs)

Usage:
  search_paired_features_xic.py <features.tsv> <shift0_Da> <nominal_L:H> [RTOL]
  (mzML is inferred by replacing 'features.tsv' -> 'mzML', like the R scripts)

Deps: numpy, pyteomics, psims  (pip install pyteomics psims)
"""
import sys, os, bisect, csv
from math import log2, fabs

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# reuse the validated port living in the envelope/ analysis dir if present,
# else fall back to a local copy of the primitives.
try:
    import xic_pair as XP
except ImportError:
    _here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(_here, '..', '..', 'envelope'))
    import xic_pair as XP

MAX_NP   = 3       # (1)  keep Np > MAX_NP  (>=4)
A        = 10.0
RTOL_DEF = 0.75    # (2)
FEAT_PPM = 10.0    # (3)
TOL      = 0.001   # EIC m/z half-window (Da)


def main():
    fn_feat = sys.argv[1]
    shift0  = float(sys.argv[2])
    nominal = float(sys.argv[3])
    rtol    = float(sys.argv[4]) if len(sys.argv) > 4 else RTOL_DEF
    log2nom = log2(nominal)
    fn_mzml = fn_feat.replace('features.tsv', 'mzML')

    # (3) per-charge sorted feature m/z, for the "partner is a feature" test
    feat = {}
    with open(fn_feat) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            z = int(r['charge'])
            feat.setdefault(z, []).append(float(r['mz']))
    for z in feat:
        feat[z].sort()

    def is_feature(z, mz):
        a = feat.get(z, []); t = mz * FEAT_PPM / 1e6
        i = bisect.bisect_left(a, mz - t)
        return i < len(a) and a[i] <= mz + t

    # candidates via XIC from mzML (Np > MAX_NP kept)
    cands = XP.candidates(fn_feat, fn_mzml, shift0, tol=TOL, max_Np=MAX_NP)

    kept = []
    for c in cands:
        if c['corr'] <= -c['Np'] / A + 4.0 / A + 0.95:     # adaptive corr gate
            continue
        if not is_feature(c['charge'], c['H']):            # (3) self-pair removal
            continue
        r = c['ratio']
        if fabs(log2(max(r, 1e-3)) - log2nom) > rtol:      # (2) ratio cut
            continue
        kept.append(c)

    for m in XP.merge(kept):
        print("MZ: %6.4f, charge: %d, L: %6.4f, H: %6.4f, RT: %6.4f, R: %6.4f, "
              "DB:( %3.2f %3.2f %d ), RT( %4.3f %4.3f )" % (
              m['L'], m['charge'], m['L'], m['H'], m['rtApex'] * 60.0, m['intR'],
              m['corr'], m['intApex'], m['Np'], m['rtStart'] * 60.0, m['rtEnd'] * 60.0))


if __name__ == '__main__':
    main()
