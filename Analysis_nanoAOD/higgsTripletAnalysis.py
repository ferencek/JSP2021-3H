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


# output histogram file
histo_filename = "HISTOGRAMS_TRSM_XToHY_6b_M3_%i_M2_%i.root" % (options.mX, options.mY)
if options.massPoint:
    histo_filename = "HISTOGRAMS_TRSM_XToHY_6b_%s.root" % options.massPoint
if options.withNu:
    histo_filename = histo_filename.replace(".root", "_WithNu.root")
f = ROOT.TFile(histo_filename, "RECREATE")
f.cd()


# define histograms
h_jetmass = ROOT.TH1F("h_jetmass", "h_jetmass", 300,0,300)
h_jetpt = ROOT.TH1F("h_jetpt", "h_jetpt", 250,0, 2500)
h_jetphi = ROOT.TH1F("h_jetphi", "h_jetphi", 150,-4,4)
h_jeteta = ROOT.TH1F("h_jeteta", "h_jeteta", 150,-4,4)

h_higgsmass = ROOT.TH1F("h_higgsmass", "h_higgsmass", 300,0,300)
h_higgspt = ROOT.TH1F("h_higgspt", "h_higgspt", 250,0, 2500)
h_higgsphi = ROOT.TH1F("h_higgsphi", "h_higgsphi", 150,-4,4)
h_higgseta = ROOT.TH1F("h_higgseta", "h_higgseta", 150,-4,4)

h_higgsmass_matched = ROOT.TH1F("h_higgsmass_matched", "h_higgsmass_matched", 300,0,300)
h_higgspt_matched = ROOT.TH1F("h_higgspt_matched", "h_higgspt_matched", 250,0, 2500)
h_higgsphi_matched = ROOT.TH1F("h_higgsphi_matched", "h_higgsphi_matched", 300,-4,4)
h_higgseta_matched = ROOT.TH1F("h_higgseta_matched", "h_higgseta_matched", 300,-4,4)

h_jetmass_matched = ROOT.TH1F("h_jetmass_matched", "h_jetmass_matched", 1000,0,1000)
h_jetpt_matched = ROOT.TH1F("h_jetpt_matched", "h_jetpt_matched", 250,0, 2500)
h_jetphi_matched = ROOT.TH1F("h_jetphi_matched", "h_jetphi_matched", 150,-4,4)
h_jeteta_matched = ROOT.TH1F("h_jeteta_matched", "h_jeteta_matched", 150,-4,4)

h_jet_pt_vs_higgs_pt = ROOT.TH2F("h_jet_pt_vs_higgs_pt", ";p^{Higgs}_{T} [GeV];p^{jet}_{T} [GeV]", 300,-100,2400,300,0,2000)
h_jet_mass_vs_higgs_pt= ROOT.TH2F("h_jet_mass_vs_higgs_pt", ";p^{Higgs}_{T} [GeV]; mass_jet [GeV]", 300,0,1200,300,100,1200)
h_DeltaR_vs_higgs_pt = ROOT.TH2F("h_DeltaR_y_vs_ptH", ";p^{Higgs}_{T} [GeV];#DeltaR_y",300,0,1500,300,0,0.2)
h_min_DR_vs_higgs_pt = ROOT.TH2F("h_min_DR_vs_higgs_pt", ";p^{Higgs}_{T} [GeV];#DeltaR_min",300,100,1500, 300,0,0.5)
h_max_DR_vs_higgs_pt = ROOT.TH2F("h_max_DR_vs_higgs_pt", ";p^{Higgs}_{T} [GeV];#DeltaR_max",300,100,1500,300,0,0.8)

h_higgs_pt_all= ROOT.TH1F("h_higgs_pt_all", ";p^{Higgs}_{T} [GeV]",300,0,2000)
h_DeltaR_bb_vs_higgspt = ROOT.TH2F("h_DeltaR_bb_vs_higgspt", ";p^{Higgs}_{T} [GeV];#DeltaR(b,b)",300,100,1500,300,0,1.5)
h_multiplicityN_higgs_candidates = ROOT.TH1F("h_multiplicityN_higgs_candidates", "h_multiplicityN_higgs_candidates", 5,-0.5,4.5)
h_multiplicityN_higgs_candidates_boosted = ROOT.TH1F("h_multiplicityN_higgs_candidates_boosted", "h_multiplicityN_higgs_candidates_boosted", 5,-0.5,4.5)


# input file
ifile = "/STORE/ferencek/TRSM_XToHY_6b/2017/13TeV/GEN/TRSM_XToHY_6b_M3_%i_M2_%i_GEN.root" % (options.mX, options.mY)
if options.massPoint:
    ifile = "/STORE/ferencek/TRSM_XToHY_6b/2017/13TeV/GEN/TRSM_XToHY_6b_%s_GEN.root" % options.massPoint

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
    if options.maxEvents > 0 and (i+1) > options.maxEvents :
        break
    if i % options.reportEvery == 0 :
        print ('Event: %i' % (i+1) )
    event.getByLabel(gpLabel, gpHandle)
    genparticles = gpHandle.product()
    event.getByLabel(jetLabel, jetHandle)
    jets = jetHandle.product()

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
            continue
        h_higgs_pt_all.Fill(gp.pt())
        if abs(gp.eta()) < 2:
            h_higgsmass.Fill(gp.mass())
            h_higgsphi.Fill(gp.phi())
            h_higgseta.Fill(gp.eta())
            h_higgspt.Fill(gp.pt())
            higgsList.append(gp)
            d1=gp.daughter(0)
            d2=gp.daughter(1)
            dphi=DeltaPhi(d1.phi(), d2.phi())
            dy=abs(d1.rapidity()-d2.rapidity())
            DeltaR = hypot(dphi, dy)
            h_DeltaR_bb_vs_higgspt.Fill(gp.pt(),DeltaR)
            if DeltaR < 0.8:
                higgscount +=1


    higgs_candidatesList=[]
    for jet in jets:
        h_jetmass.Fill(jet.mass())
        h_jetphi.Fill(jet.phi())
        h_jeteta.Fill(jet.eta())
        h_jetpt.Fill(jet.pt())
        for h in higgsList:
            dphi=DeltaPhi(jet.phi(),h.phi())
            dy=abs(jet.rapidity() - h.rapidity())
            DeltaR=hypot(dy, dphi)
            if DeltaR < 0.2:
                h_higgsmass_matched.Fill(h.mass())
                h_higgsphi_matched.Fill(h.phi())
                h_higgseta_matched.Fill(h.eta())
                h_higgspt_matched.Fill(h.pt())
                h_jetmass_matched.Fill(jet.mass())
                h_jetphi_matched.Fill(jet.phi())
                h_jeteta_matched.Fill(jet.eta())
                h_jetpt_matched.Fill(jet.pt())
                h_jet_pt_vs_higgs_pt.Fill(h.pt(),jet.pt())
                h_jet_mass_vs_higgs_pt.Fill(h.pt(),jet.mass())
                dphi1=DeltaPhi(jet.phi(), h.daughter(0).phi())
                dy1=abs(jet.rapidity() - h.daughter(0).rapidity())
                dphi2=DeltaPhi(jet.phi(), h.daughter(1).phi())
                dy2=abs(jet.rapidity() - h.daughter(1).rapidity())
                dR1=hypot(dy1, dphi1)
                dR2=hypot(dy2, dphi2)
                h_min_DR_vs_higgs_pt.Fill(h.pt(),min(dR1,dR2))
                h_max_DR_vs_higgs_pt.Fill(h.pt(),max(dR1,dR2))
                h_DeltaR_vs_higgs_pt.Fill(h.pt(),DeltaR)
        if (jet.pt() > 250 and abs(jet.eta()) < 2 and jet.mass() > 100 and jet.mass() < 150):
            higgs_candidatesList.append(jet)

    h_multiplicityN_higgs_candidates.Fill(len(higgs_candidatesList))
    h_multiplicityN_higgs_candidates_boosted.Fill(higgscount)
    


f.Write()
f.Close()
