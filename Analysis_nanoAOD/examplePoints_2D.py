import ROOT as r
import os


r.gStyle.SetPadTickX(1)  # to get the tick marks on the opposite side of the frame
r.gStyle.SetPadTickY(1)
r.gStyle.SetOptTitle(0)
r.gStyle.SetOptStat(0)

# tweak margins
r.gStyle.SetPadLeftMargin(0.14);
r.gStyle.SetPadRightMargin(0.10);

# tweak axis title offsets
r.gStyle.SetTitleOffset(1.2, "X");
r.gStyle.SetTitleOffset(2.0, "Y");
#r.gStyle.SetTitleOffset(1.25, "Z");
r.gROOT.ForceStyle()


from optparse import OptionParser
parser = OptionParser()

parser.add_option("--withNu", action="store_true",
                  dest="withNu",
                  help="Include neutrinos in GenJets",
                  default=False)

(options, args) = parser.parse_args()

points_list = [(4000, 3700), (4000, 1700), (4000, 300), (2600, 1300), (2400, 300), (2000, 1700), (600, 260)]
hist_list = {
    "h_jet_pt_vs_higgs_pt": (0, 2500, 0, 3000, r"p_{T}^{Higgs}\left[GeV\right]", r"p_{T}^{jet}\left[GeV\right]", None)
}

h_dict = dict()
min_dict, max_dict = dict(), dict()


folder = "/users/ldulibic/GENhiggs/Analysis/GEN_testAnalysisCondor_20211021105733"
if options.withNu:
    folder = "/users/ldulibic/GENhiggs/Analysis/GEN_testAnalysisCondor_20211021105733/"

c = r.TCanvas('c', 'c', 1000, 1000)

for i, (mX, mY) in enumerate(points_list):
    f = r.TFile(os.path.join(folder, "HISTOGRAMS_TRSM_XToHY_6b_M3_%i_M2_%i%s.root" % (mX, mY, ("_WithNu" if options.withNu else ""))), "READ")
    print '{} {} loaded'.format(mX, mY)

    c.cd()
    for hist_name in hist_list:
        min_x, max_x, min_y, max_y, x_label, y_label, new_bin_count = hist_list[hist_name]
        hist = f.Get(hist_name)
        hist.GetXaxis().SetTitle(x_label)
        hist.GetYaxis().SetTitle(y_label)
        hist.Draw("colz")
        c.SaveAs(hist_name + '_{}_{}{}.pdf'.format(mX, mY, ("_WithNu" if options.withNu else "")))
