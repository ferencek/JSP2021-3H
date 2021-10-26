mX_min = 400
mX_max = 4000
mX_step = 200
mY_step = 200

counter = 0

for mX in range(mX_min, mX_max + mX_step, mX_step):
    for mY in sorted(list(set([260,mX-140])) + range(300, mX-125, mY_step)):
        print ("(mX, mY) = (%i, %i)" % (mX, mY))
        counter += 1

print ("Total number of mass points: %i" % counter)
