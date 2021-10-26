import ROOT
import math
import sys

from optparse import OptionParser

usage = "Usage: python %prog [options] \n\nExample: python %prog --mX 3000 --mY 300"

parser = OptionParser(usage=usage)

parser.add_option('--mX', type='float', action='store',
                  dest='mX',
                  help='X scalar mass')

parser.add_option('--mY', type='float', action='store',
                  dest='mY',
                  help='Y scalar mass')

(options, args) = parser.parse_args()

# make sure all necessary input parameters are provided
if not (options.mX and options.mY):
    print ('Mandatory parameters missing\n')
    parser.print_help()
    sys.exit(1)

mX = options.mX
mY = options.mY
mH = 125.

Delta = mX**2 - mY**2 - mH**2
pz = math.sqrt((Delta**2 - 4 * mH**2 * mY**2)/(4 * mX**2))

#print pz

tX = ROOT.TLorentzVector()
tX.SetPxPyPzE(0,0,0,mX)
#tX.Print()

tY = ROOT.TLorentzVector()
tY.SetPxPyPzE(0,0,pz,math.sqrt(mY**2 + pz**2))
#tY.Print()

print ("Momenta:\n")

tH1 = ROOT.TLorentzVector()
tH1.SetPxPyPzE(0,0,-pz,math.sqrt(mH**2 + pz**2))
print ("H1")
tH1.Print()

#tSumX = tH1 + tY
#tSumX.Print()

boostY = tY.BoostVector()
#tY.Boost(-boostY)
#tY.Print()

DeltaY = mY**2 - 2 * mH**2
pzH = math.sqrt((DeltaY**2 - 4 * mH**4)/(4 * mY**2))

tH2 = ROOT.TLorentzVector()
tH2.SetPxPyPzE(0,0,-pzH,math.sqrt(mH**2 + pzH**2))
tH2.Boost(boostY)
print ("H2")
tH2.Print()

tH3 = ROOT.TLorentzVector()
tH3.SetPxPyPzE(0,0,pzH,math.sqrt(mH**2 + pzH**2))
tH3.Boost(boostY)
print ("H2")
tH3.Print()

#tSumY = tH2 + tH3
#tSumY.Print()

print ("\nMasses:\n")

print ("m12 = %f GeV^2" % (tH1 + tH2).M2())
print ("m13 = %f GeV^2" % (tH1 + tH3).M2())
print ("m23 = %f GeV^2" % (tH2 + tH3).M2())
