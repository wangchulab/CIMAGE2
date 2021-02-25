import sys
from math import fabs, log

ppm = 1e-6
tol = 4.0 * ppm #>2
rt_gap = 0.2 #min

charge = int(sys.argv[2])

lines = open(sys.argv[1], 'r').readlines()
tags = { k:n for n,k in enumerate(lines[0].strip().split('\t')) }

class feature:
    def __init__(self):
        self.rtStart = 0.0
        self.rtEnd = 0.0
        self.mz = 0.0
        self.rtApex = 0.0
        self.mzApex = 0.0
        self.intApex = 0.0
        self.intSum = 0.0
    
    def output(self):
        scaleH = self.intSum / (self.rtEnd - self.rtStart)
        return "(MZ: %6.4f, RT: %6.4f, Apex: %6.4f, scaleH: %6.4f)" \
        % (self.mz, self.rtApex*60.0, self.intApex, scaleH)

features = []
for n, l in enumerate(lines[1:]):
    es = l.strip().split('\t')
    c = int(es[tags['charge']])
    if c != charge: continue

    feat = feature()
    feat.rtStart = float(es[tags['rtStart']])
    feat.rtEnd = float(es[tags['rtEnd']])
    feat.mz = float(es[tags['mz']])
    feat.mzApex = float(es[tags['mostAbundantMz']])
    feat.rtApex = float(es[tags['rtApex']])
    feat.intApex = float(es[tags['intensityApex']])
    feat.intSum = float(es[tags['intensitySum']])
    
    #merge
    mz_tol = tol * feat.mz
    merged = False
    for ncan, candy in enumerate(features):
        d_mz = candy.mz - feat.mz
        d_mzA = candy.mzApex - feat.mzApex
        if fabs(d_mz) < mz_tol and fabs(d_mzA) < mz_tol:
            #split
            if candy.rtStart - feat.rtEnd > rt_gap:
                continue
            if feat.rtStart - candy.rtEnd > rt_gap:
                continue

            #merge
            #print "merge", d_mz, d_mzA
            candy.intSum += feat.intSum
            if feat.rtEnd > candy.rtEnd:
                candy.rtEnd = feat.rtEnd
            if feat.rtStart < candy.rtStart:
                candy.rtStart = feat.rtStart
            if feat.intApex > candy.intApex:
                candy.intApex = feat.intApex
                candy.mzApex = feat.mzApex
                candy.rtApex = feat.rtApex
                candy.mz = feat.mz
            merged = True
            break

    if not merged: features.append(feat)

#sorted by mz
sorted_features = sorted(features, key=lambda x: x.mz)
#for ncan, candy in enumerate(sorted_features):
#    print ncan, candy.output()
#print ncan+1

shift = float(sys.argv[3]) / charge
nlast = len(sorted_features)
p0 = 0
p1 = 0
while p1 < nlast:
    mz0 = sorted_features[p0].mz
    mz1 = sorted_features[p1].mz
    mz_tol = tol * mz1
    #print p0, p1, (mz1 - mz0)*charge
    if mz1 - mz0 > shift + mz_tol:
        p0 += 1
        continue
    if mz1 - mz0 < shift - mz_tol:
        p1 += 1
        continue

    #match
    delta_rt = fabs(sorted_features[p0].rtApex - sorted_features[p1].rtApex)
    #print delta_rt
    if delta_rt < rt_gap:
        log_fc = log(sorted_features[p0].intApex/sorted_features[p1].intApex)
        #print log_fc
        if fabs(log_fc) < 4.0:
            #check rt range
            rt0_s = sorted_features[p0].rtStart
            rt0_e = sorted_features[p0].rtEnd
            rt1_s = sorted_features[p1].rtStart
            rt1_e = sorted_features[p1].rtEnd
            status = True
            r = 0
            if rt0_s <= rt1_s and rt0_e >= rt1_e:
                #print "1 in 0"
                statuts = True
                r = 1.0
            elif rt1_s <= rt0_s and rt1_e >= rt0_e:
                #print "0 in 1"
                status = True
                r = 1.0
            elif rt0_e < rt1_s:
                status = False
            elif rt0_s > rt1_e:
                status = False
            elif rt0_s < rt1_s:
                inter = rt0_e - rt1_s
                r0 = inter / (rt0_e - rt0_s)
                r1 = inter / (rt1_e - rt1_s)
                r = max(r0, r1)
                if r > 0.2: status = True
            elif rt1_s < rt0_s:
                inter = rt1_e - rt0_s
                r0 = inter / (rt0_e - rt0_s)
                r1 = inter / (rt1_e - rt1_s)
                r = max(r0, r1)
                if r > 0.2: status = True
            else:
                print "FAIL"
                assert(False)

            #match!
            if status:
                print "match0:", charge, sorted_features[p0].output()
                print "match1:", charge, sorted_features[p1].output()
                print "dt:", delta_rt, "log_fc:", log_fc, "dm:", ((mz1-mz0)-shift)/shift
                print "--", rt0_s, rt0_e, rt1_s, rt1_e, r
    
    #do the smallest move
    if p0 >= p1-1: p1 += 1
    elif mz1 - sorted_features[p0+1].mz < shift - mz_tol:
        p1 += 1
    else:
        p0 += 1

