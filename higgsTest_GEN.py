from DataFormats.FWLite import Events, Handle
from math import hypot
import ROOT
import sys

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat("nemruoi")
ROOT.gROOT.ForceStyle()

from optparse import OptionParser
parser = OptionParser()

parser.add_option('--maxEvents', type='int', action='store',
                  default=-1,
                  dest='maxEvents',
                  help='Number of events to run. -1 is all events')

parser.add_option('--reportEvery', type='int', action='store',
                  default=100,
                  dest='reportEvery',
                  help='Report every N events')

parser.add_option('--mX', type='int', action='store',
                  default=3000,
                  dest='mX',
                  help='X scalar mass')

parser.add_option('--mY', type='int', action='store',
                  default=300,
                  dest='mY',
                  help='Y scalar mass')

parser.add_option('--massPoint', action='store',
                  dest='massPoint',
                  help='Mass point')

parser.add_option("--withNu", action="store_true",
                  dest="withNu",
                  help="Include neutrinos in GenJets",
                  default=False)

(options, args) = parser.parse_args()


def DeltaPhi(v1, v2, c = 3.141592653589793):
    r = (v2 - v1) % (2.0 * c)
    if r < -c:
        r += 2.0 * c
    elif r > c:
        r -= 2.0 * c
    return abs(r)

# single input file for testing
ifile = "./data/TRSM_XToHY_6b_M3_2800_M2_700_GEN.root"

# open input file
events = Events(ifile)

# define collections to process
gpHandle = Handle ("std::vector<reco::GenParticle>")
jetHandle = Handle ("std::vector<reco::GenJet>")
gpLabel = ("genParticles")
jetLabel= ("ak8GenJetsNoNu")
if options.withNu:
    jetLabel = ("ak8GenJets")

# loop over events
for i,event in enumerate(events): 
	event.getByLabel(gpLabel, gpHandle)
	genparticles = gpHandle.product()
	event.getByLabel(jetLabel, jetHandle)
	jets = jetHandle.product()

	nDaughters = 0
	higgsList=[]
	higgscount=0
	for gp in genparticles:
		if not gp.pdgId()==25:
		    continue
		hasHiggsDaughter = False
		for d in range(gp.numberOfDaughters()):
		    if gp.daughter(d).pdgId()==25:
		        hasHiggsDaughter = True
		        break
		if hasHiggsDaughter:
		    nDaughters += 1
		    continue
	if i > 100:
		break

print(nDaughters)
