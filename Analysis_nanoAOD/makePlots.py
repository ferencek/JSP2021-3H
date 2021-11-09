import ROOT as r
import copy
import os

from optparse import OptionParser
parser = OptionParser()

parser.add_option("--withNu", action="store_true",
                  dest="withNu",
                  help="Include neutrinos in GenJets",
                  default=False)

parser.add_option('--omjet', '--outputMjet', action='store',
                  dest='condor_outputMjet',
                  help='Output folder of analysis in which the mjet histograms are stored')

parser.add_option('--omsoftdrop', '--outputMsoftdrop', action='store',
                  dest='condor_outputMsoftdrop',
                  help='Output folder of analysis in which the msoftdrop histograms are stored')


(options, args) = parser.parse_args()

#---------------------------------------------------------------------
# to run in the batch mode (to prevent canvases from popping up)
r.gROOT.SetBatch()

# set plot style
r.gROOT.SetStyle("Plain")
r.gStyle.SetPalette(57)

# suppress the statistics box
r.gStyle.SetOptStat(0)

# more detailed statistics box
#r.gStyle.SetOptStat("nemruoi")
#r.gStyle.SetOptStat(1111111)

# suppress the histogram title
r.gStyle.SetOptTitle(0)

r.gStyle.SetPadTickX(1)  # to get the tick marks on the opposite side of the frame
r.gStyle.SetPadTickY(1)  # to get the tick marks on the opposite side of the frame

# tweak margins
r.gStyle.SetPadTopMargin(0.1);
r.gStyle.SetPadBottomMargin(0.1);
r.gStyle.SetPadLeftMargin(0.10);
r.gStyle.SetPadRightMargin(0.15);

# tweak axis title offsets
r.gStyle.SetTitleOffset(1.25, "Y");
r.gStyle.SetTitleOffset(1.25, "Z");

# set nicer fonts
r.gStyle.SetTitleFont(42, "")
r.gStyle.SetTitleFont(42, "XYZ")
r.gStyle.SetLabelFont(42, "XYZ")
r.gStyle.SetTextFont(42)
r.gStyle.SetStatFont(42)
r.gROOT.ForceStyle()

#---------------------------------------------------------------------
def plot(graph, graphs_BP, name, plotBP=False):

    #############################
    # The "surface" drawing options ("surf", "col", "cont") impose their own internal coordinate system
    # which complicates overlay of additional objects on the plot and requires extra steps. The code
    # below is based on discussion and instructions provided in
    # https://root-forum.cern.ch/t/drawing-a-single-contour-over-a-tgraph2d-using-cont4/12061
    #############################

    c = r.TCanvas("c", "",1000,1000)
    c.cd()

    pad1 = r.TPad("pad1","",0,0,1,1)
    pad2 = r.TPad("pad2","",0.1,0.1,0.85,0.9)
    pad2.SetFillStyle(4000) # will be transparent

    pad1.Draw()
    pad1.cd()

    graph.SetMinimum(0) # has to be placed before calling TAxis methods (?!)
    graph.SetMaximum(1) # has to be placed before calling TAxis methods (?!)

    graph.GetXaxis().SetTickLength(0)
    graph.GetYaxis().SetTickLength(0)

    graph.Draw("cont4z")

    # get color palette
    pad1.Update()
    palette = graph.GetHistogram().GetListOfFunctions().FindObject("palette")

    xmin = graph.GetXaxis().GetXmin()
    xmax = graph.GetXaxis().GetXmax()
    ymin = graph.GetYaxis().GetXmin()
    ymax = graph.GetYaxis().GetXmax()

    graph.GetXaxis().SetLimits(graph.GetXaxis().GetXmin()/1000.,graph.GetXaxis().GetXmax()/1000.)
    graph.GetYaxis().SetLimits(graph.GetYaxis().GetXmin()/1000.,graph.GetYaxis().GetXmax()/1000.)

    #print xmin, xmax, ymin, ymax

    pad2.Range(xmin,ymin,xmax,ymax)
    pad2.Draw()
    pad2.cd()

    if plotBP:
        gr_list = []
        for gr_BP in graphs_BP:
            gr_list.append(r.TGraph())
            x = gr_BP.GetX()[0]
            y = gr_BP.GetY()[0]
            z = gr_BP.GetZ()[0]
            gr_list[-1].SetPoint(0,x,y)
            ci = palette.GetValueColor(z)
            gr_list[-1].SetMarkerStyle(8)
            gr_list[-1].SetMarkerSize(1.5)
            gr_list[-1].SetMarkerColor(ci)
            gr_list[-1].Draw("SAME P")
            gr_list.append(r.TGraph())
            gr_list[-1].SetPoint(0,x,y)
            gr_list[-1].SetMarkerStyle(4)
            gr_list[-1].SetMarkerSize(1.6)
            gr_list[-1].SetMarkerColor(r.kRed)
            gr_list[-1].Draw("SAME P")
    else:
        gr_dummy = r.TGraph()
        gr_dummy.SetPoint(0,1000,3000)
        gr_dummy.Draw("SAME *") # necessary to get the tick marks (unless the BP graphs are drawn)

    pline = r.TPolyLine()
    pline.SetPoint(0,xmin,xmin-125.)
    pline.SetPoint(1,xmax,xmax-125.)
    pline.SetPoint(2,xmin,xmax-125.)
    pline.SetPoint(3,xmin,xmin-125.)
    pline.SetFillColor(10)
    pline.SetFillStyle(1001)
    pline.SetLineColor(2)
    pline.SetLineWidth(2)
    pline.Draw("f")
    pline_hatched = copy.deepcopy(pline)
    pline_hatched.SetFillColor(r.kGray)
    #pline_hatched.SetFillStyle(3344)
    pline_hatched.SetFillStyle(3002)
    pline_hatched.Draw("f")

    pad2.RedrawAxis()

    c.SaveAs(name)

#---------------------------------------------------------------------
# regular mass points
gr_GenPart         = copy.deepcopy(r.TGraph2D())

gr_FatJet          = copy.deepcopy(r.TGraph2D())
gr_FatJetMatched   = copy.deepcopy(r.TGraph2D())
gr_FatJetUnmatched = copy.deepcopy(r.TGraph2D())

gr_FatJetDeepTag          = copy.deepcopy(r.TGraph2D())
gr_FatJetMatchedDeepTag   = copy.deepcopy(r.TGraph2D())
gr_FatJetUnmatchedDeepTag = copy.deepcopy(r.TGraph2D())

gr_FatJetParticleNet          = copy.deepcopy(r.TGraph2D())
gr_FatJetMatchedParticleNet   = copy.deepcopy(r.TGraph2D())
gr_FatJetUnmatchedParticleNet = copy.deepcopy(r.TGraph2D())

# softdrop mass 
gr_GenPart_msd         = copy.deepcopy(r.TGraph2D())
gr_FatJet_msd          = copy.deepcopy(r.TGraph2D())
gr_FatJetMatched_msd   = copy.deepcopy(r.TGraph2D())
gr_FatJetUnmatched_msd = copy.deepcopy(r.TGraph2D())

gr_FatJetDeepTag_msd          = copy.deepcopy(r.TGraph2D())
gr_FatJetMatchedDeepTag_msd   = copy.deepcopy(r.TGraph2D())
gr_FatJetUnmatchedDeepTag_msd = copy.deepcopy(r.TGraph2D())

gr_FatJetParticleNet_msd          = copy.deepcopy(r.TGraph2D())
gr_FatJetMatchedParticleNet_msd   = copy.deepcopy(r.TGraph2D())
gr_FatJetUnmatchedParticleNet_msd = copy.deepcopy(r.TGraph2D())

# ratios
gr_msdTOmjetRatioMatched   = copy.deepcopy(r.TGraph2D()) # matched softdrop / matched mjet
gr_msdTOmjetRatioUnmatched = copy.deepcopy(r.TGraph2D()) # unmatched softdrop / unmatched mjet

gr_msdTOmjetRatioMatchedDeepTag   = copy.deepcopy(r.TGraph2D()) # matched softdrop / matched mjet
gr_msdTOmjetRatioUnmatchedDeepTag = copy.deepcopy(r.TGraph2D()) # unmatched softdrop / unmatched mjet

gr_msdTOmjetRatioMatchedParticleNet   = copy.deepcopy(r.TGraph2D()) # matched softdrop / matched mjet
gr_msdTOmjetRatioUnmatchedParticleNet = copy.deepcopy(r.TGraph2D()) # unmatched softdrop / unmatched mjet

gr_matchedTOallRatio     = copy.deepcopy(r.TGraph2D()) # matched / all 
gr_matchedTOallRatio_msd = copy.deepcopy(r.TGraph2D()) # matched / all with softdrop mass

gr_matchedTOallRatioDeepTag     = copy.deepcopy(r.TGraph2D()) # matched / all 
gr_matchedTOallRatioDeepTag_msd = copy.deepcopy(r.TGraph2D()) # matched / all with softdrop mass

gr_matchedTOallRatioParticleNet     = copy.deepcopy(r.TGraph2D()) # matched / all 
gr_matchedTOallRatioParticleNet_msd = copy.deepcopy(r.TGraph2D()) # matched / all with softdrop mass


gr_GenPart.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates (GenPart)")
gr_FatJet.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates (FatJet)")
gr_FatJetMatched.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates matched to Higgs particle (FatJet)")
gr_FatJetUnmatched.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates unmatched to Higgs particle (FatJet)")

gr_FatJetDeepTag.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates (FatJet-DeepTag)")
gr_FatJetMatchedDeepTag.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates matched to Higgs particle (FatJet-DeepTag)")
gr_FatJetUnmatchedDeepTag.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates unmatched to Higgs particle (FatJet-DeepTag)")

gr_FatJetParticleNet.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates (FatJet-ParticleNet)")
gr_FatJetMatchedParticleNet.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates matched to Higgs particle (FatJet-ParticleNet)")
gr_FatJetUnmatchedParticleNet.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates unmatched to Higgs particle (FatJet-ParticleNet)")

gr_GenPart_msd.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates (GenPart)")
gr_FatJet_msd.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates (FatJet)")
gr_FatJetMatched_msd.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates matched to Higgs particle (FatJet)")
gr_FatJetUnmatched_msd.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates unmatched to Higgs particle (FatJet)")

gr_FatJetDeepTag_msd.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates (FatJet-DeepTag)")
gr_FatJetMatchedDeepTag_msd.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates matched to Higgs particle (FatJet-DeepTag)")
gr_FatJetUnmatchedDeepTag_msd.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates unmatched to Higgs particle (FatJet-DeepTag)")

gr_FatJetParticleNet_msd.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates (FatJet-ParticleNet)")
gr_FatJetMatchedParticleNet_msd.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates matched to Higgs particle (FatJet-ParticleNet)")
gr_FatJetUnmatchedParticleNet_msd.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Fraction of boosted Higgs boson candidates unmatched to Higgs particle (FatJet-ParticleNet)")

gr_msdTOmjetRatioMatched.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Softdrop mass over jet mass ratio (matched)")
gr_msdTOmjetRatioUnmatched.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Softdrop mass over jet mass ratio (unmatched)")

gr_msdTOmjetRatioMatchedDeepTag.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Softdrop mass over jet mass ratio (matched)")
gr_msdTOmjetRatioUnmatchedDeepTag.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Softdrop mass over jet mass ratio (unmatched)")

gr_msdTOmjetRatioMatchedParticleNet.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Softdrop mass over jet mass ratio (matched)")
gr_msdTOmjetRatioUnmatchedParticleNet.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Softdrop mass over jet mass ratio (unmatched)")

gr_matchedTOallRatio.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Ratio of matched over all (mjet)")
gr_matchedTOallRatio_msd.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Ratio of matched over all (msoftdrop)")

gr_matchedTOallRatioDeepTag.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Ratio of matched over all (DeepTag)(mjet)")
gr_matchedTOallRatioDeepTag_msd.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Ratio of matched over all (DeepTag)(msoftdrop)")

gr_matchedTOallRatioParticleNet.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Ratio of matched over all (ParticleNet)(mjet)")
gr_matchedTOallRatioParticleNet_msd.SetTitle(";m_{X} [TeV];m_{Y} [TeV];Ratio of matched over all (ParticleNet)(msoftdrop)")

boosted_higgs_graphsGenPart = [copy.deepcopy(r.TGraph2D()) for i in range(4)]
boosted_higgs_graphsFatJet  = [copy.deepcopy(r.TGraph2D()) for i in range(4)]
boosted_higgs_graphsGenPart_msd = [copy.deepcopy(r.TGraph2D()) for i in range(4)]
boosted_higgs_graphsFatJet_msd  = [copy.deepcopy(r.TGraph2D()) for i in range(4)]

for i in range(4):
    boosted_higgs_graphsGenPart[i].SetTitle(";m_{X} [TeV];m_{Y} [TeV];Event selection eff. (%i H cand.)"%i)
    boosted_higgs_graphsFatJet[i].SetTitle(";m_{X} [TeV];m_{Y} [TeV];Event selection eff. (%i H cand.)"%i)
    boosted_higgs_graphsGenPart_msd[i].SetTitle(";m_{X} [TeV];m_{Y} [TeV];Event selection eff. (%i H cand.)"%i)
    boosted_higgs_graphsFatJet_msd[i].SetTitle(";m_{X} [TeV];m_{Y} [TeV];Event selection eff. (%i H cand.)"%i)

mX_min  = 400
mX_max  = 4000
mX_step = 400
mY_step = 400

values = open('mass_point_values.txt', 'w')
header = '{:3s}  {:^4s}  {:^4s}  {:8s}  {:11s}  {:5s}  {:9s}  {:5s}  {:9s}  {:5s}  {:9s}\n'.format('idx','mX','mY','frac_GenPart','frac_FatJet','eff_1','eff_jet_1','eff_2','eff_jet_2','eff_3','eff_jet_3')
values.write(header)
values.write('-' * (len(header)-1) + '\n')

n = 0
for mX in range(mX_min, mX_max + mX_step, mX_step):
    for mY in sorted(list(set([260,mX-140])) + range(300, mX-125, mY_step)):
        values.write('{:>3d}  {:>4d}  {:>4d}'.format(n+1, mX, mY))

        readHist = "/users/ldulibic/nanoAODhiggs/Analysis_nanoAOD/"+options.condor_outputMjet+"/nanoAOD_HISTOGRAMS_TRSM_XToHY_6b_M3_%i_M2_%i_FatJet.root" % (mX,mY)
        if options.withNu:
            readHist = readHist.replace(".root", "_WithNu.root")
        f = r.TFile(readHist)

        readHist = "/users/ldulibic/nanoAODhiggs/Analysis_nanoAOD/"+options.condor_outputMsoftdrop+"/nanoAOD_HISTOGRAMS_TRSM_XToHY_6b_M3_%i_M2_%i_FatJet_msoftdrop.root" % (mX,mY)
        if options.withNu:
            readHist = readHist.replace(".root", "_WithNu.root")
        g = r.TFile(readHist)

        f.cd()
        pt_all  = f.Get('h_higgs_pt_all')   # all higgs in all events
        boosted = f.Get('h_HCands_boosted') # gen part boosted higgs
        candidates = f.Get("h_HCands")           # jet candidates for higgs
        matched    = f.Get('h_HCands_matched')   # jet candidates matched to higgs genpart
        unmatched  = f.Get('h_HCands_unmatched') # jet candidates not matched to higgs genpart 
        candidatesDeepTag = f.Get("h_HCands_deeptag")           
        matchedDeepTag    = f.Get('h_HCands_matched_deeptag')  
        unmatchedDeepTag  = f.Get('h_HCands_unmatched_deeptag') 
        candidatesParticleNet = f.Get("h_HCands_particlenet")           
        matchedParticleNet    = f.Get('h_HCands_matched_particlenet')   
        unmatchedParticleNet  = f.Get('h_HCands_unmatched_particlenet')  

        g.cd()
        pt_all_msd = g.Get('h_higgs_pt_all')    # all higgs in all events
        boosted_msd = g.Get('h_HCands_boosted') # gen part boosted higgs
        candidates_msd = g.Get("h_HCands")          # jet candidates for higgs
        matched_msd = g.Get('h_HCands_matched')     # jet candidates matched to higgs genpart
        unmatched_msd = g.Get('h_HCands_unmatched') # jet candidates not matched to higgs genpart
        candidatesDeepTag_msd = g.Get("h_HCands_deeptag")           
        matchedDeepTag_msd    = g.Get('h_HCands_matched_deeptag')  
        unmatchedDeepTag_msd  = g.Get('h_HCands_unmatched_deeptag') 
        candidatesParticleNet_msd = g.Get("h_HCands_particlenet")           
        matchedParticleNet_msd    = g.Get('h_HCands_matched_particlenet')   
        unmatchedParticleNet_msd  = g.Get('h_HCands_unmatched_particlenet')  

        nHCandsBoosted     = 0
        nHCandsBoosted_msd = 0
        for i in range(1,5):
            nHCandsBoosted     += boosted.GetBinContent(i+1)*i
            nHCandsBoosted_msd += boosted_msd.GetBinContent(i+1)*i
        frac_GenPart     = float(nHCandsBoosted) / pt_all.Integral(0, pt_all.GetNbinsX()+1)
        frac_GenPart_msd = float(nHCandsBoosted_msd) / pt_all_msd.Integral(0, pt_all_msd.GetNbinsX()+1)

        nHCands            = 0
        nHCandsDeepTag     = 0
        nHCandsParticleNet = 0
        nHCands_msd            = 0
        nHCandsDeepTag_msd     = 0
        nHCandsParticleNet_msd = 0
        for i in range(1,5):
            nHCands            += candidates.GetBinContent(i+1)*i
            nHCandsDeepTag     += candidatesDeepTag.GetBinContent(i+1)*i
            nHCandsParticleNet += candidatesParticleNet.GetBinContent(i+1)*i
            nHCands_msd            += candidates_msd.GetBinContent(i+1)*i
            nHCandsDeepTag_msd     += candidatesDeepTag_msd.GetBinContent(i+1)*i
            nHCandsParticleNet_msd += candidatesParticleNet_msd.GetBinContent(i+1)*i

        frac_FatJet            = float(nHCands) / pt_all.Integral(0, pt_all.GetNbinsX()+1)
        frac_FatJetDeepTag     = float(nHCandsDeepTag) / pt_all.Integral(0, pt_all.GetNbinsX()+1)
        frac_FatJetParticleNet = float(nHCandsParticleNet) / pt_all.Integral(0, pt_all.GetNbinsX()+1)

        frac_FatJet_msd            = float(nHCands_msd) / pt_all_msd.Integral(0, pt_all_msd.GetNbinsX()+1)
        frac_FatJetDeepTag_msd     = float(nHCandsDeepTag_msd) / pt_all_msd.Integral(0, pt_all_msd.GetNbinsX()+1)
        frac_FatJetParticleNet_msd = float(nHCandsParticleNet_msd) / pt_all_msd.Integral(0, pt_all_msd.GetNbinsX()+1)

        nHCandsMatched            = 0
        nHCandsMatchedDeepTag     = 0
        nHCandsMatchedParticleNet = 0
        nHCandsMatched_msd            = 0
        nHCandsMatchedDeepTag_msd     = 0
        nHCandsMatchedParticleNet_msd = 0
        for i in range(1,5):
            nHCandsMatched            += matched.GetBinContent(i+1)*i
            nHCandsMatchedDeepTag     += matchedDeepTag.GetBinContent(i+1)*i
            nHCandsMatchedParticleNet += matchedParticleNet.GetBinContent(i+1)*i
            nHCandsMatched_msd            += matched_msd.GetBinContent(i+1)*i
            nHCandsMatchedDeepTag_msd     += matchedDeepTag_msd.GetBinContent(i+1)*i
            nHCandsMatchedParticleNet_msd += matchedParticleNet_msd.GetBinContent(i+1)*i

        frac_FatJetMatched            = float(nHCandsMatched) / pt_all.Integral(0, pt_all.GetNbinsX()+1)
        frac_FatJetMatchedDeepTag     = float(nHCandsMatchedDeepTag) / pt_all.Integral(0, pt_all.GetNbinsX()+1)
        frac_FatJetMatchedParticleNet = float(nHCandsMatchedParticleNet) / pt_all.Integral(0, pt_all.GetNbinsX()+1)

        frac_FatJetMatched_msd            = float(nHCandsMatched_msd) / pt_all_msd.Integral(0, pt_all_msd.GetNbinsX()+1)
        frac_FatJetMatchedDeepTag_msd     = float(nHCandsMatchedDeepTag_msd) / pt_all_msd.Integral(0, pt_all_msd.GetNbinsX()+1)
        frac_FatJetMatchedParticleNet_msd = float(nHCandsMatchedParticleNet_msd) / pt_all_msd.Integral(0, pt_all_msd.GetNbinsX()+1)

        nHCandsUnmatched     = 0
        nHCandsUnmatchedDeepTag     = 0
        nHCandsUnmatchedParticleNet = 0
        nHCandsUnmatched_msd = 0
        nHCandsUnmatchedDeepTag_msd     = 0
        nHCandsUnmatchedParticleNet_msd = 0
        for i in range(1,5):
            nHCandsUnmatched     += unmatched.GetBinContent(i+1)*i
            nHCandsUnmatchedDeepTag     += unmatchedDeepTag.GetBinContent(i+1)*i
            nHCandsUnmatchedParticleNet += unmatchedParticleNet.GetBinContent(i+1)*i
            nHCandsUnmatched_msd += unmatched_msd.GetBinContent(i+1)*i
            nHCandsUnmatchedDeepTag_msd     += unmatchedDeepTag_msd.GetBinContent(i+1)*i
            nHCandsUnmatchedParticleNet_msd += unmatchedParticleNet_msd.GetBinContent(i+1)*i

        frac_FatJetUnmatched            = float(nHCandsUnmatched) / pt_all.Integral(0, pt_all.GetNbinsX()+1)
        frac_FatJetUnmatchedDeepTag     = float(nHCandsUnmatchedDeepTag) / pt_all.Integral(0, pt_all.GetNbinsX()+1)
        frac_FatJetUnmatchedParticleNet = float(nHCandsUnmatchedParticleNet) / pt_all.Integral(0, pt_all.GetNbinsX()+1)

        frac_FatJetUnmatched_msd            = float(nHCandsUnmatched_msd) / pt_all_msd.Integral(0, pt_all_msd.GetNbinsX()+1)
        frac_FatJetUnmatchedDeepTag_msd     = float(nHCandsUnmatchedDeepTag_msd) / pt_all_msd.Integral(0, pt_all_msd.GetNbinsX()+1)
        frac_FatJetUnmatchedParticleNet_msd = float(nHCandsUnmatchedParticleNet_msd) / pt_all_msd.Integral(0, pt_all_msd.GetNbinsX()+1)

        # print ("(mX, mY) = (%i, %i)" % (mX, mY))
        # print (frac_GenPart)
        # print (frac_FatJet)

        values.write('    {:.3f}      {:.3f}   '.format(frac_GenPart, frac_FatJet))

        # setting values for:
        # generator
        gr_GenPart.SetPoint(n,mX,mY,frac_GenPart)
        # fatjet 
        gr_FatJet.SetPoint(n,mX,mY,frac_FatJet)
        gr_FatJetDeepTag.SetPoint(n,mX,mY,frac_FatJetDeepTag)
        gr_FatJetParticleNet.SetPoint(n,mX,mY,frac_FatJetParticleNet)
        # fatjet matched
        gr_FatJetMatched.SetPoint(n,mX,mY,frac_FatJetMatched)
        gr_FatJetMatchedDeepTag.SetPoint(n,mX,mY,frac_FatJetMatchedDeepTag)
        gr_FatJetMatchedParticleNet.SetPoint(n,mX,mY,frac_FatJetMatchedParticleNet)
        # fatjet unmatched
        gr_FatJetUnmatched.SetPoint(n,mX,mY,frac_FatJetUnmatched)
        gr_FatJetUnmatchedDeepTag.SetPoint(n,mX,mY,frac_FatJetUnmatchedDeepTag)
        gr_FatJetUnmatchedParticleNet.SetPoint(n,mX,mY,frac_FatJetUnmatchedParticleNet)
        # all the same but with softdrop mass
        gr_GenPart_msd.SetPoint(n,mX,mY,frac_GenPart_msd)

        gr_FatJet_msd.SetPoint(n,mX,mY,frac_FatJet_msd)
        gr_FatJetDeepTag_msd.SetPoint(n,mX,mY,frac_FatJetDeepTag_msd)
        gr_FatJetParticleNet_msd.SetPoint(n,mX,mY,frac_FatJetParticleNet_msd)

        gr_FatJetMatched_msd.SetPoint(n,mX,mY,frac_FatJetMatched_msd)
        gr_FatJetMatchedDeepTag_msd.SetPoint(n,mX,mY,frac_FatJetMatchedDeepTag_msd)
        gr_FatJetMatchedParticleNet_msd.SetPoint(n,mX,mY,frac_FatJetMatchedParticleNet_msd)

        gr_FatJetUnmatched_msd.SetPoint(n,mX,mY,frac_FatJetUnmatched_msd)    
        gr_FatJetUnmatchedDeepTag_msd.SetPoint(n,mX,mY,frac_FatJetUnmatchedDeepTag_msd)     
        gr_FatJetUnmatchedParticleNet_msd.SetPoint(n,mX,mY,frac_FatJetUnmatchedParticleNet_msd)     

        # ratios
        gr_msdTOmjetRatioMatched.SetPoint(n,mX,mY,frac_FatJetMatched_msd/frac_FatJetMatched)
        gr_msdTOmjetRatioUnmatched.SetPoint(n,mX,mY,frac_FatJetUnmatched_msd/frac_FatJetUnmatched)   
        gr_matchedTOallRatio.SetPoint(n,mX,mY,frac_FatJetMatched/frac_FatJet)
        gr_matchedTOallRatio_msd.SetPoint(n,mX,mY,frac_FatJetMatched_msd/frac_FatJet_msd)

        gr_msdTOmjetRatioMatchedDeepTag.SetPoint(n,mX,mY,frac_FatJetMatchedDeepTag_msd/frac_FatJetMatchedDeepTag)
        gr_msdTOmjetRatioUnmatchedDeepTag.SetPoint(n,mX,mY,frac_FatJetUnmatchedDeepTag_msd/frac_FatJetUnmatchedDeepTag)   
        gr_matchedTOallRatioDeepTag.SetPoint(n,mX,mY,frac_FatJetMatchedDeepTag/frac_FatJetDeepTag)
        gr_matchedTOallRatioDeepTag_msd.SetPoint(n,mX,mY,frac_FatJetMatchedDeepTag_msd/frac_FatJetDeepTag_msd)

        gr_msdTOmjetRatioMatchedParticleNet.SetPoint(n,mX,mY,frac_FatJetMatchedParticleNet_msd/frac_FatJetMatchedParticleNet)
        gr_msdTOmjetRatioUnmatchedParticleNet.SetPoint(n,mX,mY,frac_FatJetUnmatchedParticleNet_msd/frac_FatJetUnmatchedParticleNet)   
        gr_matchedTOallRatioParticleNet.SetPoint(n,mX,mY,frac_FatJetMatchedParticleNet/frac_FatJetParticleNet)
        gr_matchedTOallRatioParticleNet_msd.SetPoint(n,mX,mY,frac_FatJetMatchedParticleNet_msd/frac_FatJetParticleNet_msd)

        for count in range(4):
            frac_GenPart     = boosted.GetBinContent(count+1) / boosted.Integral()
            frac_FatJet      = candidates.GetBinContent(count+1) / candidates.Integral()
            frac_GenPart_msd = boosted_msd.GetBinContent(count+1) / boosted_msd.Integral()
            frac_FatJet_msd  = candidates_msd.GetBinContent(count+1) / candidates_msd.Integral()

            boosted_higgs_graphsGenPart[count].SetPoint(n, mX, mY, frac_GenPart)
            boosted_higgs_graphsFatJet[count].SetPoint(n, mX, mY, frac_FatJet)
            boosted_higgs_graphsGenPart_msd[count].SetPoint(n, mX, mY, frac_GenPart)
            boosted_higgs_graphsFatJet_msd[count].SetPoint(n, mX, mY, frac_FatJet)

            if count > 0: values.write('  {:.3f}    {:.3f}  '.format(frac_GenPart, frac_FatJet))
        values.write('\n')
        n += 1

values.close()

#---------------------------------------------------------------------
# benchmark points
# points = [(1600, 500), (2000, 300), (2000, 800), (2500, 300)]
# suffix = ["BPb", "BPd", "BPe", "BPf"]

# gr_GenPart_BP = [copy.deepcopy(r.TGraph2D()) for point in points]
# gr_FatJet_BP = [copy.deepcopy(r.TGraph2D()) for point in points]
# boosted_higgs_graphsGenPart_BP = [[copy.deepcopy(r.TGraph2D()) for point in points] for i in range(4)]
# boosted_higgs_graphsFatJet_BP = [[copy.deepcopy(r.TGraph2D()) for point in points] for i in range(4)]


# for p in range(len(points)):
#     gr_GenPart_BP[p].SetTitle(";m_{X} [GeV];m_{Y} [GeV];Fraction of boosted Higgs boson candidates")
#     gr_FatJet_BP[p].SetTitle(";m_{X} [GeV];m_{Y} [GeV];Fraction of boosted Higgs boson candidates")
#     for i in range(4):
#         boosted_higgs_graphsGenPart_BP[i][p].SetTitle(";m_{X} [GeV];m_{Y} [GeV];Event selection eff. (%i H cand.)"%i)
#         boosted_higgs_graphsFatJet_BP[i][p].SetTitle(";m_{X} [GeV];m_{Y} [GeV];Event selection eff. (%i H cand.)"%i)

# for p, (mX,mY) in enumerate(points):
#         n = 0
#         f = r.TFile("/users/ldulibic/GEnH/Analysis/GEN_testAnalysisCondor_20211021105733/HISTOGRAMS_TRSM_XToHY_6b_M3_%i_M2_%i_%s.root" % (mX, mY, suffix[p]))
#         f.cd()
#         candidates = f.Get("h_HCands")
#         h2_b = f.Get('h_DeltaR_bb_vs_higgspt')
#         pt_all = f.Get('pt_all')
#         boosted = f.Get('h_HCands_boosted')
#         frac_GenPart = h2_b.Integral(0,h2_b.GetNbinsX()+1,0,h2_b.GetYaxis().FindBin(0.8)-1)/ pt_all.Integral(0,pt_all.GetNbinsX()+1)
#         nHCands=0
#         for i in range(1,5):
#             nHCands += candidates.GetBinContent(i+1)*i
#         frac_FatJet=float(nHCands)/pt_all.Integral(0,pt_all.GetNbinsX()+1)
#         print ("(mX, mY) = (%i, %i)" % (mX, mY))
#         print (frac_GenPart)
#         print (frac_FatJet)
#         gr_GenPart_BP[p].SetPoint(n,mX,mY,frac_GenPart)
#         gr_FatJet_BP[p].SetPoint(n,mX,mY,frac_FatJet)
#         for count in range(4):
#             frac_GenPart = boosted.GetBinContent(count+1) / boosted.Integral()
#             frac_FatJet = candidates.GetBinContent(count+1) / candidates.Integral()
#             boosted_higgs_graphsGenPart_BP[count][p].SetPoint(n, mX, mY, frac_GenPart)
#             boosted_higgs_graphsFatJet_BP[count][p].SetPoint(n, mX, mY, frac_FatJet)

#---------------------------------------------------------------------

# creating folder structure
try:
    os.mkdir('figs')
except OSError as error:
    pass
try:
    os.mkdir('figs/figsMjet')
except OSError as error:
    pass
try:
    os.mkdir('figs/figsMsoftdrop')
except OSError as error:
    pass
try:
    os.mkdir('figs/figsRatios')
except OSError as error:
    pass

# make plots
plot(gr_GenPart, [], "./figs/figsMjet/BoostedHiggsFraction_GenPart.pdf")
plot(gr_FatJet, [], "./figs/figsMjet/BoostedHiggsFraction_FatJet.pdf")
plot(gr_FatJetMatched,[],"./figs/figsMjet/BoostedHiggsFraction_FatJet_matched.pdf")
plot(gr_FatJetUnmatched,[],"./figs/figsMjet/BoostedHiggsFraction_FatJet_unmatched.pdf")
plot(gr_FatJetDeepTag, [], "./figs/figsMjet/BoostedHiggsFraction_FatJet_Deeptag.pdf")
plot(gr_FatJetMatchedDeepTag,[],"./figs/figsMjet/BoostedHiggsFraction_FatJet_DeepTag_matched.pdf")
plot(gr_FatJetUnmatchedDeepTag,[],"./figs/figsMjet/BoostedHiggsFraction_FatJet__DeepTag_unmatched.pdf")
plot(gr_FatJetParticleNet, [], "./figs/figsMjet/BoostedHiggsFraction_FatJet_ParticleNet.pdf")
plot(gr_FatJetMatchedParticleNet,[],"./figs/figsMjet/BoostedHiggsFraction_FatJet_ParticleNet_matched.pdf")
plot(gr_FatJetUnmatchedParticleNet,[],"./figs/figsMjet/BoostedHiggsFraction_FatJet_ParticleNet_unmatched.pdf")

plot(gr_GenPart_msd, [], "./figs/figsMsoftdrop/BoostedHiggsFraction_GenPart_msoftdrop.pdf")
plot(gr_FatJet_msd, [], "./figs/figsMsoftdrop/BoostedHiggsFraction_FatJet_msoftdrop.pdf")
plot(gr_FatJetMatched_msd,[],"./figs/figsMsoftdrop/BoostedHiggsFraction_FatJet_matched_msoftdrop.pdf")
plot(gr_FatJetUnmatched_msd,[],"./figs/figsMsoftdrop/BoostedHiggsFraction_FatJet_unmatched_msoftdrop.pdf")
plot(gr_FatJetDeepTag_msd, [], "./figs/figsMsoftdrop/BoostedHiggsFraction_FatJet_DeepTag_msoftdrop.pdf")
plot(gr_FatJetMatchedDeepTag_msd,[],"./figs/figsMsoftdrop/BoostedHiggsFraction_FatJet_DeepTag_matched_msoftdrop.pdf")
plot(gr_FatJetUnmatchedDeepTag_msd,[],"./figs/figsMsoftdrop/BoostedHiggsFraction_FatJet_DeepTag_unmatched_msoftdrop.pdf")
plot(gr_FatJetParticleNet_msd, [], "./figs/figsMsoftdrop/BoostedHiggsFraction_FatJet_Particlenet_msoftdrop.pdf")
plot(gr_FatJetMatchedParticleNet_msd,[],"./figs/figsMsoftdrop/BoostedHiggsFraction_FatJet_ParticleNet_matched_msoftdrop.pdf")
plot(gr_FatJetUnmatchedParticleNet_msd,[],"./figs/figsMsoftdrop/BoostedHiggsFraction_FatJet_ParticleNet_unmatched_msoftdrop.pdf")

plot(gr_msdTOmjetRatioMatched,[],"./figs/figsRatios/msdTOmjetRatioMatched.pdf")
plot(gr_msdTOmjetRatioUnmatched,[],"./figs/figsRatios/msdTOmjetRatioUnmatched.pdf")
plot(gr_matchedTOallRatio,[],"./figs/figsRatios/MatchedToAllRatio.pdf")
plot(gr_matchedTOallRatio_msd,[],"./figs/figsRatios/MatchedToAllRatio_msoftdrop.pdf")

plot(gr_msdTOmjetRatioMatchedDeepTag,[],"./figs/figsRatios/msdTOmjetRatioMatchedDeepTag.pdf")
plot(gr_msdTOmjetRatioUnmatchedDeepTag,[],"./figs/figsRatios/msdTOmjetRatioUnmatchedDeepTag.pdf")
plot(gr_matchedTOallRatioDeepTag,[],"./figs/figsRatios/MatchedToAllRatioDeepTag.pdf")
plot(gr_matchedTOallRatioDeepTag_msd,[],"./figs/figsRatios/MatchedToAllRatioDeepTag_msoftdrop.pdf")

plot(gr_msdTOmjetRatioMatchedParticleNet,[],"./figs/figsRatios/msdTOmjetRatioMatchedParticleNet.pdf")
plot(gr_msdTOmjetRatioUnmatchedParticleNet,[],"./figs/figsRatios/msdTOmjetRatioUnmatchedparticleNet.pdf")
plot(gr_matchedTOallRatioParticleNet,[],"./figs/figsRatios/MatchedToAllRatioParticleNet.pdf")
plot(gr_matchedTOallRatioParticleNet_msd,[],"./figs/figsRatios/MatchedToAllRatioParticleNet_msoftdrop.pdf")

# for i in range(4):
#     plot(boosted_higgs_graphsGenPart[i], [], "./figs/figsMjet/Event_Selection_eff_%i_GenPart.pdf"%i)
#     plot(boosted_higgs_graphsFatJet[i], [], "./figs/figsMjet/Event_Selection_eff_%i_FatJet.pdf"%i) 

#     plot(boosted_higgs_graphsGenPart_msd[i], [], "./figs/figsMsoftdrop/Event_Selection_eff_%i_GenPart_msoftdrop.pdf"%i)
#     plot(boosted_higgs_graphsFatJet_msd[i], [], "./figs/figsMsoftdrop/Event_Selection_eff_%i_FatJet_msoftdrop.pdf"%i)

#---------------------------------------------------------------------
