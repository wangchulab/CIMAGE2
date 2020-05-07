import numpy as np
from math import fabs
from collections import defaultdict
import sys

lines = open(sys.argv[1], 'r').readlines()
if len(sys.argv) > 2:
  gap = int(sys.argv[2])
else:
  gap = 6

rt_lst = []
mz_lst = []
it_lst = []

map_rt_it = defaultdict(float)

for l in lines:
  es = l.strip().split()
  rt = str("%4.3f" % float(es[0]))
  mz = float(es[1])
  it = float(es[2])
  rt_lst.append(rt)
  mz_lst.append(mz)
  it_lst.append(it)
  map_rt_it[rt] += it

best_rt = 0
best_it = 0
for rt in map_rt_it.keys():
  if map_rt_it[rt] > best_it:
    best_rt = float(rt)
    best_it = map_rt_it[rt]

tol = 0.001
peak_mz = []
peak_it = []
best_mz = 0
best_it = 0
for rt, mz, it in zip( rt_lst, mz_lst, it_lst ):
  if fabs( float(rt) - best_rt ) < tol:
    if it>best_it:
      best_mz = mz
      best_it = it
    peak_mz.append( mz )
    peak_it.append( it )
    #print mz, it

charge = int(round( 1.0/ (peak_mz[1]-peak_mz[0]) ))
mzL = np.min(peak_mz)-0.1/charge
mzH = np.max(peak_mz)+0.1/charge
test1 = best_mz - (1.0078 * gap)/charge
test2 = best_mz + (1.0078 * gap)/charge
#print charge, best_mz, best_it, test1, test2

def get_it( mzs, its, qmz ):
  for mz, it in zip( mzs, its ):
    if fabs( mz - qmz )<0.05:
      return it

if test1>mzL and test1<mzH:
  it_sec = get_it( peak_mz, peak_it, test1 )
  if it_sec > best_it/4.0:
    print sys.argv[1], charge, best_rt, test1, best_mz
if test2>mzL and test2<mzH:
  it_sec = get_it( peak_mz, peak_it, test2 )
  if it_sec > best_it/4.0:
    print sys.argv[1], charge, best_rt, best_mz, test2
