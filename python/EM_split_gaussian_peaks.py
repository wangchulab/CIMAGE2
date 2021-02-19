#!/usr/bin/env python
import sys
from math import fabs
import random
import numpy as np

tol_d = 0.05 #ACR:0.05 #HNE/IA:0.04
min_size = 50

d_lst = []
for l in open(sys.argv[1], 'r').readlines():
  try:
    d_lst.append(float(l))
  except:
    pass

x = np.array(d_lst)

#fast-K
nd = len(x)
clusters = []
seed = random.randint(0,nd)
clusters.append([seed])

for i in xrange(nd):
  min_d = 999
  min_c = -1
  for n,c in enumerate(clusters):
    d = fabs(x[i] - x[c[0]])
    if d < min_d:
      min_d = d
      min_c = n
  if min_d > tol_d:
    clusters.append([i])
  else:
    clusters[min_c].append(i)

print [len(c) for c in clusters]

gauss = []
for c in clusters:
  if len(c)>min_size:
    #fit cluster
    m = np.median([x[i] for i in c])
    gauss.append(m)

#print gauss
nc = len(gauss)
gauss = np.array(gauss)
sig = 0.004
eng = np.zeros([nc, nd]) #x^2/s^2
pc = np.zeros([nc, nd])  #prob of cluster

for step in xrange(20):
  #calc P for each x
  for n in xrange(nc):
    c = gauss[n]
    eng[n,:] = ((x - c)/sig)**2/2
  for i in xrange(nd):
    pc[:,i] = np.exp(-eng[:,i])
    sum_p = np.sum(pc[:,i])
    #print "pc:", pc[:,i]
    if sum_p > 1e-6:
      pc[:,i] /= sum_p
      #print "/sum_p", pc[:,i]
    else:
      pc[:,i] = 0.0
  
  #calc gaussian mean
  sum_p = np.zeros([nc])
  sum_px = np.zeros([nc])
  for i in xrange(nd):
    #print sum_p.shape, pc[:,i]
    sum_p += pc[:,i]
    sum_px += pc[:,i]*x[i]
  
  cs = sum_px/sum_p
  d2 = np.sum(((gauss-cs)/sig)**2)
  if d2<0.0001: break
  gauss = cs

#process results
print gauss
print sum_p
