from DataFormats.FWLite import Events, Handle
from math import hypot
import ROOT

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

(options, args) = parser.parse_args()


# output histogram file
histo_filename = "HISTOGRAMS_MINIAOD_TRSM_XToHY_6b_M3_%i_M2_%i.root" % (options.mX, options.mY)
if options.massPoint:
    histo_filename = "HISTOGRAMS_MINIAOD_TRSM_XToHY_6b_%s.root" % options.massPoint
f = ROOT.TFile(histo_filename, "RECREATE")
f.cd()


# define histograms
h_PV_multiplicity = ROOT.TH1F("h_PV_multiplicity", "h_PV_multiplicity", 61,-0.5,60.5)


# input file
ifile = "/STORE/ferencek/TRSM_XToHY_6b/2017/13TeV/MINIAOD/TRSM_XToHY_6b_M3_%i_M2_%i_MINIAOD_0.root" % (options.mX, options.mY)
if options.massPoint:
    ifile = "/STORE/ferencek/TRSM_XToHY_6b/2017/13TeV/MINIAOD/TRSM_XToHY_6b_%s_MINIAOD_0.root" % options.massPoint

# open input file
events = Events(ifile)

# define collections to process
pvHandle = Handle ("vector<reco::Vertex>")
pvLabel = ("offlineSlimmedPrimaryVertices")

# loop over events
for i,event in enumerate(events):
    if options.maxEvents > 0 and (i+1) > options.maxEvents :
        break
    if i % options.reportEvery == 0 :
        print ('Event: %i' % (i+1) )
    event.getByLabel(pvLabel, pvHandle)
    primaryVertices = pvHandle.product()

    h_PV_multiplicity.Fill(primaryVertices.size())


f.Write()
f.Close()
