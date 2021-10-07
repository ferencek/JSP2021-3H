import ROOT
from math import hypot
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

def findDaughters(n,partList):
    # loop over list of particle indices and check if any of them have 
    # the n-th particle as the mother, if yes then n-th particle
    # has m-th partcile as her daughter

    # return list of daughters' indicies 
    dIndex = []
    for m in partList:
        if n == e.GenPart_genPartIdxMother[m]:
            dIndex.append(m)

    return dIndex

# output histogram file
histo_filename = "./hist/nanoAOD_HISTOGRAMS_TRSM_XToHY_6b_M3_%i_M2_%i_FatJet.root" % (options.mX, options.mY)
if options.massPoint:
    histo_filename = "./hist/nanoAOD_HISTOGRAMS_TRSM_XToHY_6b_%s_FatJet.root" % options.massPoint
if options.withNu:
    histo_filename = histo_filename.replace(".root", "_WithNu.root")

# # single output file for testing
# histo_filename = "nanoAOD_HISTOGRAMS_TRSM_XToHY_6b_M3_2800_M2_700.root"

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
ifile = "/STORE/ferencek/TRSM_XToHY_6b/2017/13TeV/NANOAOD/TRSM_XToHY_6b_M3_%i_M2_%i_NANOAOD.root" % (options.mX, options.mY)
if options.massPoint:
    ifile = "/STORE/ferencek/TRSM_XToHY_6b/2017/13TeV/NANOAOD/TRSM_XToHY_6b_%s_NANOAOD.root" % options.massPoint

# # single input file for testing
# ifile = "./data/TRSM_XToHY_6b_M3_2800_M2_700_NANOAOD.root"

# open root input file directly 
evtFile = ROOT.TFile.Open(ifile)
events  = evtFile.Get("Events")

# for deleting previous printed line to have nice reportEvery 
CURSOR_UP_ONE = '\x1b[1A' 
ERASE_LINE = '\x1b[2K' 

nDaughters = 0
# loop over events
for i,e in enumerate(events):
    if options.maxEvents > 0 and (i+1) > options.maxEvents :
        break
    if i % options.reportEvery == 0:
        print('Event: %i' %(i+1))
        sys.stdout.write(CURSOR_UP_ONE) 
        sys.stdout.write(ERASE_LINE) 

    higgsList=[]
    higgscount=0

    # empty Lorentz4vector for calculating rapidity
    partVec1 = ROOT.TLorentzVector()
    partVec2 = ROOT.TLorentzVector()

    for n in range(e.nGenPart):
        # list of all particle indices
        allPart = [i for i in range(e.nGenPart)] 

        if not e.GenPart_pdgId[n] == 25:
            continue
        hasHiggsDaughter = False
        
        dIndex = findDaughters(n,allPart)
        for m in dIndex:
            if e.GenPart_pdgId[m]==25:
                hasHiggsDaughter = True
                break
        if hasHiggsDaughter: # at least one!
            nDaughters += 1
            continue

        h_higgs_pt_all.Fill(e.GenPart_pt[n])

        if abs(e.GenPart_eta[n]) < 2:
            h_higgsmass.Fill(e.GenPart_mass[n])
            h_higgsphi.Fill(e.GenPart_phi[n])
            h_higgseta.Fill(e.GenPart_eta[n])
            h_higgspt.Fill(e.GenPart_pt[n])
            higgsList.append(n)

            dIndex = findDaughters(n,allPart)

            partVec1.SetPtEtaPhiM(e.GenPart_pt[dIndex[0]],e.GenPart_eta[dIndex[0]],e.GenPart_phi[dIndex[0]],e.GenPart_mass[dIndex[0]])
            partVec2.SetPtEtaPhiM(e.GenPart_pt[dIndex[1]],e.GenPart_eta[dIndex[1]],e.GenPart_phi[dIndex[1]],e.GenPart_mass[dIndex[1]])
            
            dphi=DeltaPhi(e.GenPart_phi[dIndex[0]], e.GenPart_phi[dIndex[1]])
            dy=abs(partVec1.Rapidity()-partVec2.Rapidity())
            DeltaR = hypot(dphi, dy)

            h_DeltaR_bb_vs_higgspt.Fill(e.GenPart_pt[n],DeltaR)
            if DeltaR < 0.8:
                higgscount +=1

    # using FatJet - ak8 Jets made with visible genparticles
    higgs_candidatesList=[]

    # empty Lorentz4vector for calculating rapidity
    jetVec = ROOT.TLorentzVector()
    hVec   = ROOT.TLorentzVector()

    for j in range(e.nFatJet):
        h_jetmass.Fill(e.FatJet_mass[j])
        h_jetphi.Fill(e.FatJet_phi[j])
        h_jeteta.Fill(e.FatJet_eta[j])
        h_jetpt.Fill(e.FatJet_pt[j])

        for h in higgsList:    
            hVec.SetPtEtaPhiM(e.GenPart_pt[h],e.GenPart_eta[h],e.GenPart_phi[h],e.GenPart_mass[h])
            jetVec.SetPtEtaPhiM(e.FatJet_pt[j],e.FatJet_eta[j],e.FatJet_phi[j],e.FatJet_mass[j])
            
            dphi=DeltaPhi(e.FatJet_phi[j],e.GenPart_phi[h])
            dy=abs(jetVec.Rapidity() - hVec.Rapidity())
            DeltaR=hypot(dy, dphi)

            if DeltaR < 0.2:
                h_higgsmass_matched.Fill(e.GenPart_mass[h])
                h_higgsphi_matched.Fill(e.GenPart_phi[h])
                h_higgseta_matched.Fill(e.GenPart_eta[h])
                h_higgspt_matched.Fill(e.GenPart_pt[h])
                h_jetmass_matched.Fill(e.FatJet_mass[j])
                h_jetphi_matched.Fill(e.FatJet_phi[j])
                h_jeteta_matched.Fill(e.FatJet_eta[j])
                h_jetpt_matched.Fill(e.FatJet_pt[j])
                h_jet_pt_vs_higgs_pt.Fill(e.GenPart_pt[h],e.FatJet_pt[j])
                h_jet_mass_vs_higgs_pt.Fill(e.GenPart_pt[h],e.FatJet_mass[j])

                dIndex = findDaughters(h,allPart)
        
                partVec1.SetPtEtaPhiM(e.GenPart_pt[dIndex[0]],e.GenPart_eta[dIndex[0]],e.GenPart_phi[dIndex[0]],e.GenPart_mass[dIndex[0]])
                partVec2.SetPtEtaPhiM(e.GenPart_pt[dIndex[1]],e.GenPart_eta[dIndex[1]],e.GenPart_phi[dIndex[1]],e.GenPart_mass[dIndex[1]])
            
                dphi1=DeltaPhi(e.FatJet_phi[j], e.GenPart_phi[dIndex[0]])
                dy1=abs(jetVec.Rapidity() - partVec1.Rapidity())
                dphi2=DeltaPhi(e.FatJet_phi[j], e.GenPart_phi[dIndex[1]])
                dy2=abs(jetVec.Rapidity() - partVec2.Rapidity())

                dR1=hypot(dy1, dphi1)
                dR2=hypot(dy2, dphi2)
                h_min_DR_vs_higgs_pt.Fill(e.GenPart_pt[h],min(dR1,dR2))
                h_max_DR_vs_higgs_pt.Fill(e.GenPart_pt[h],max(dR1,dR2))
                h_DeltaR_vs_higgs_pt.Fill(e.GenPart_pt[h],DeltaR)

        if (e.FatJet_pt[j] > 250 and abs(e.FatJet_eta[j]) < 2 and e.FatJet_mass[j] > 100 and e.FatJet_mass[j] < 150):
            higgs_candidatesList.append(j)

    h_multiplicityN_higgs_candidates.Fill(len(higgs_candidatesList))
    h_multiplicityN_higgs_candidates_boosted.Fill(higgscount)
    
#print(nDaughters)

f.Write()
f.Close()
