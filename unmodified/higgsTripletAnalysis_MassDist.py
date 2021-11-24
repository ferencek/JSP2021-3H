from DataFormats.FWLite import Events, Handle
import ROOT


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

parser.add_option('--massPoint', action='store',
                  default='M3_2000_M2_300',
                  dest='massPoint',
                  help='Mass point')

(options, args) = parser.parse_args()


gpHandle = Handle ("std::vector<reco::GenParticle>")
gpLabel = ("genParticles")


f = ROOT.TFile("MassDistributions_TRSM_XToHY_6b_%s.root" % options.massPoint, "RECREATE")
f.cd()

h_m3H_distribution = ROOT.TH1F("h_m3H_distribution", ";m_{3H} [GeV]",2000,0,4000)
h_mH2H3_vs_mH1H2 = ROOT.TH2F("h_mH2H3_vs_mH1H2", ";m_{H1H2} [GeV];m_{H2H3} [GeV]",1000,0,2000,1000,0,2000)

ifile = "/STORE/ferencek/TRSM_XToHY_6b/2017/13TeV/GEN/TRSM_XToHY_6b_%s_GEN.root" % options.massPoint

events = Events(ifile)

for i,event in enumerate(events): 
    if options.maxEvents > 0 and (i+1) > options.maxEvents :
        break
    if i % options.reportEvery == 0 :
        print ('Event: %i' % (i+1) )
    event.getByLabel(gpLabel, gpHandle)
    genparticles = gpHandle.product()

    h=[]
    for gp in genparticles:
        if not gp.pdgId()==25:
            continue
        hasHiggsDaughter = False
        for d in range(gp.numberOfDaughters()):
            if gp.daughter(d).pdgId()==25:
                hasHiggsDaughter = True
                break
        if not hasHiggsDaughter:
            h.append(gp)

    #print('# of Higgs bosons found: %i' % len(h))

    h.sort(key=lambda x: x.phi())

    tH1 = ROOT.TLorentzVector()
    tH1.SetPtEtaPhiE(h[0].pt(),h[0].eta(),h[0].phi(),h[0].energy())

    tH2 = ROOT.TLorentzVector()
    tH2.SetPtEtaPhiE(h[1].pt(),h[1].eta(),h[1].phi(),h[1].energy())

    tH3 = ROOT.TLorentzVector()
    tH3.SetPtEtaPhiE(h[2].pt(),h[2].eta(),h[2].phi(),h[2].energy())

    h_m3H_distribution.Fill( (tH1+tH2+tH3).M() )
    h_mH2H3_vs_mH1H2.Fill( (tH1+tH2).M(), (tH2+tH3).M() )


f.Write()
f.Close()
