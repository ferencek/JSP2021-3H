import os
import sys
import re

mX_min = 400
mX_max = 4000
mX_step = 200
mY_step = 200

failedMassPoints = []

content = os.listdir(sys.argv[1])

# loop over mass points
for mX in range(mX_min, mX_max + mX_step, mX_step):
    for mY in sorted(list(set([260,mX-140])) + range(300, mX-125, mY_step)):
        found = False
        for c in content:
            if ('HISTOGRAMS_TRSM_XToHY_6b_M3_%i_M2_%i*.root' % (mX, mY)) in c:
                found = True
                break
        if not found:
            failedMassPoints.append([mX, mY])

if len(failedMassPoints):
    print ('Missing output for the following %i mass points:' % len(failedMassPoints))
    print failedMassPoints
else:
    print ('No missing output')
