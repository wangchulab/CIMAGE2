#!/usr/bin/env python
import sys
import ms
import peptideParams as pep
import numpy as np
from math import fabs, exp
from collections import defaultdict
#import msplotlib

ms1_fn = sys.argv[1]
tol = float(sys.argv[2])
shifts = []
for s in sys.argv[3:]:
  shifts.append(float(s))

# cheap version: [rt, [mz_list], [inten_list]]
LC_peaks = {}
candidates = {}
for c in [ 2, 3, 4, 5, 6 ]:
  LC_peaks[c] = []
  candidates[c] = []

saved_peaks = []
saved_envlist = []
saved_rt = []

for ms1 in ms.load_ms1( ms1_fn, tol ):
  print "scan sp=", ms1.scan, ms1.rt
  current_rt = ms1.rt
  #save it for reverse appending
  saved_peaks.append(ms1.peaks)
  saved_envlist.append(ms1.envelope_list)
  saved_rt.append(current_rt)

  #skip charge=1 by now
  for charge in [ 2, 3, 4, 5, 6 ]:
    bin_size = pep.deltaC12C13 / charge

    for env in ms1.envelope_list[charge]:
      moz_array = np.array([ ms1.peaks[p][0] for p in env ]) #array of m over z
      int_array = np.array([ ms1.peaks[p][1] for p in env ]) #array of intensity
      mono_moz = moz_array[ np.argmax( int_array ) ]         #m over z of the highest peak
      #print moz_array

      marked = False
      for ncan, candy in enumerate( candidates[charge] ):
        #check if we can put something into it
        mono_moz_check = candy[-1][1][ np.argmax(candy[-1][2]) ]
        dmz = mono_moz - mono_moz_check
        nshift = round( dmz/bin_size )
        delta = fabs(nshift*bin_size + mono_moz_check - mono_moz)
        #print mono_moz_check, nshift, delta
        if delta < 0.01 and fabs(nshift) <= 1:
          #save it
          candy.append((current_rt, moz_array, int_array))
          marked = True
          break
        
      if not marked:
        candidates[charge].append([(current_rt, moz_array, int_array)])
  
      #check if the candidates can be finished
      new_candidates = []
      for ncan, candy in enumerate( candidates[charge] ):
        if current_rt - candy[-1][0] > 0.1:
          #finish this LCpeak
          if len(candy) > 20:
            LC_peaks[charge].append(candy)
        else:
          new_candidates.append(candy)
      candidates[charge] = new_candidates

#final check
for ncan, candy in enumerate( candidates[charge] ):
  #finish this
  if len(candy) > 20:
    LC_peaks[charge].append(candy)

#do reverse appending
saved_peaks.reverse()
saved_envlist.reverse()
saved_rt.reverse()

for peaks, el, rt in zip(saved_peaks, saved_envlist, saved_rt):
  current_rt = rt
  print "rev_rt=", current_rt
  for charge in [ 2, 3, 4, 5, 6 ]:
    bin_size = pep.deltaC12C13 / charge
    for env in el[charge]:
      moz_array = np.array([ peaks[p][0] for p in env ])
      int_array = np.array([ peaks[p][1] for p in env ])
      mono_moz = moz_array[ np.argmax( int_array ) ]

      #print moz_array
      for candy in LC_peaks[charge]:
        first_rt = candy[0][0] 
        if first_rt > current_rt-0.001: continue
        if first_rt < current_rt-0.1: continue
        mono_moz_check = candy[0][1][ np.argmax(candy[0][2]) ]
        dmz = mono_moz - mono_moz_check
        nshift = round( dmz/bin_size )
        delta = fabs(nshift*bin_size + mono_moz_check - mono_moz)
        if delta < 0.01 and fabs(nshift) <= 1:
          candy.insert(0, (current_rt, moz_array, int_array))    
          break

env_ndx = 0
for c in [ 2, 3, 4, 5, 6 ]:
  print "charge=", c
  for LC_peak in LC_peaks[c]:
    env_ndx += 1
    fn = "env_" + str(env_ndx) + ".dat"
    fp = open(fn, 'w')
    for env in LC_peak:
      print len(env)   
      rt = env[0]    
      for (p, it) in zip( env[1], env[2] ):
        str_peak = "%f %f %f\n" % (rt, p, it)
        fp.write(str_peak)
    fp.close()
