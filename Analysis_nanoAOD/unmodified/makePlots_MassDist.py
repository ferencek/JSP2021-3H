import ROOT as r

#---------------------------------------------------------------------
# to run in the batch mode (to prevent canvases from popping up)
r.gROOT.SetBatch()

# set plot style
r.gROOT.SetStyle("Plain")
r.gStyle.SetPalette(57)

# suppress the statistics box
r.gStyle.SetOptStat(0)

# more detailed statistics box
r.gStyle.SetOptStat("nemruoi")
#r.gStyle.SetOptStat(1111111)

# suppress the histogram title
#r.gStyle.SetOptTitle(0)

r.gStyle.SetPadTickX(1)  # to get the tick marks on the opposite side of the frame
r.gStyle.SetPadTickY(1)  # to get the tick marks on the opposite side of the frame

# tweak margins
#r.gStyle.SetPadTopMargin(0.05)
#r.gStyle.SetPadBottomMargin(0.13)
r.gStyle.SetPadLeftMargin(0.10)
r.gStyle.SetPadRightMargin(0.05)

# tweak axis title offsets
r.gStyle.SetTitleOffset(1.2, "Y")

# set nicer fonts
r.gStyle.SetTitleFont(42, "")
r.gStyle.SetTitleFont(42, "XYZ")
r.gStyle.SetLabelFont(42, "XYZ")
r.gStyle.SetTextFont(42)
r.gStyle.SetStatFont(42)
r.gROOT.ForceStyle()
#---------------------------------------------------------------------

massPoints = ['M3_1600_M2_500_BPb', 'M3_2000_M2_300_BPd', 'M3_1600_M2_500', 'M3_2000_M2_300', 'M3_2000_M2_800_BPe', 'M3_2500_M2_300_BPf']


def makePlot(f, m, histo, opt, rm):
    r.gStyle.SetPadRightMargin(rm)
    # create canvas
    c = r.TCanvas("c", "",1200,800)
    c.cd()

    h = f.Get(histo)
    h.SetTitle(m)
    h.Draw(opt)
    # need to update the gPad to be able to retrieve the stat boxes
    r.gPad.Update()

    statbox = h.GetListOfFunctions().FindObject('stats')

    x2 = statbox.GetX2NDC()
    y2 = statbox.GetY2NDC()
    x1 = statbox.GetX1NDC()
    y1 = statbox.GetY1NDC()

    statbox.SetX1NDC(x1-0.2*(x2-x1))
    statbox.SetX2NDC(x2-0.2*(x2-x1))
    statbox.SetY1NDC(y1)
    statbox.SetY2NDC(y2)

    r.gPad.Modified()
    r.gPad.Update()

    c.SaveAs("%s_TRSM_XToHY_6b_%s.png" % (histo, m))
    c.Close()
    del c


for m in massPoints:
    f = r.TFile.Open("MassDistributions_TRSM_XToHY_6b_%s.root" % m)
    makePlot(f, m, "h_m3H_distribution", "", 0.05)
    makePlot(f, m, "h_mH2H3_vs_mH1H2", "colz", 0.10)

