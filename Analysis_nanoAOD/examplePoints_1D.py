import ROOT as r


r.gStyle.SetOptStat(0)
r.gStyle.SetPadTickX(1)  # to get the tick marks on the opposite side of the frame
r.gStyle.SetPadTickY(1)
r.gStyle.SetOptTitle(0)


points_list = [(4000, 3700), (4000, 1700), (4000, 300), (2600, 1300), (2400, 300), (2000, 1700), (600, 260)]
hist_list = {
    "h_jet_pt_vs_higgs_pt": (0, 2500, 0, 3000, r"p_{T}^{Higgs}\left[GeV\right]", r"p_{T}^{jet}\left[GeV\right]", None),
    "h_jetpt_matched": (50, 1200, 0, 800, r"p_{T}^{jet}\left[GeV\right]", "Entries", None),
    "h_jetphi_matched": (-4, 4, 0, 2000, r"\Phi", "Entries", 5),
    "h_jetmass_matched": (0, 300, 0, 900, r"m^{jet}\left[GeV\right]", "Entries", None),
    "h_jeteta_matched": (-3, 3, 0, 2000, r"\eta^{jet}", "Entries", 2),
    "h_jetpt": (0, 500, 0, 6500, r"p_{T}^{jet}\left[GeV\right]", "Entries", 5),
    "h_jetphi": (-4, 4, 0, 7000, r"\Phi^{jet}\left[rad\right]", "Entries", 5),
    "h_jetmass": (0, 300, 0, 10000, r"m^{jet}\left[GeV\right]", "Entries", None),
    "h_jeteta": (-10, 10, 0, 6000, r"\eta^{jet}", "Entries", 5),
    "h_DeltaR_y_vs_ptH": (0, 1800, 0, 0.5, r"p_{T}^{Higgs}\left[GeV\right]", "\Delta R(b,b)", None),
    "h_multiplicityN_higgs_candidates_boosted": (-1, 4, 0, 14000, r"N_{candidatesBoosted}", "Entries", None),
    "h_multiplicityN_higgs_candidates":  (-1, 4, 0, 14000, r"N_{candidates}", "Entries", None),
    "h_higgsmass": (100, 150, 0, 10000, r"m_{Higgs}\left[GeV\right]", "Entries", None),
    "h_higgspt_matched": (0, 2500, 0, 800, r"p_{T}^{Higgs}\left[GeV\right]", "Entries", None),
}

h_dict = dict()
min_dict, max_dict = dict(), dict()
f_List=[]

colors=[1, 2, 3, 4, r.kOrange+1, 6, 7]
c = r.TCanvas('c', 'c', 1000, 600)
for i, (mX, mY) in enumerate(points_list):
    f = r.TFile("/users/ldulibic/GENhiggs/Analysis/GEN_testAnalysisCondor_20211021105733/HISTOGRAMS_TRSM_XToHY_6b_M3_%i_M2_%i.root" % (mX, mY), "READ")
    print("/users/ldulibic/GENhiggs/Analysis/GEN_testAnalysisCondor_20211021105733/HISTOGRAMS_TRSM_XToHY_6b_M3_%i_M2_%i.root" % (mX, mY))
    f_List.append(f)
    print '{} {} loaded'.format(mX, mY)
    for hist_name in hist_list:
        print '\t hist {} loaded'.format(hist_name)
        hist = f.Get(hist_name)
        print(hist_name)
        h_list = h_dict.get(hist_name, [])
        h_list.append(hist)
        h_dict[hist_name] = h_list
        max_val = hist.GetMaximum()
        min_val = hist.GetMinimum()
        min_dict[hist_name] = min(min_val, min_dict.get(hist_name, float('inf')))
        max_dict[hist_name] = max(max_val, max_dict.get(hist_name, float('-inf')))
        hist.SetLineStyle(i+1)
        hist.SetLineColor(colors[i])
        #hist.Draw("HIST" if i == 0 else "HIST SAMES")
        #print '\t hist {} drawn'.format(hist_name)
        #c.Modified()
        #c.Update()
    #f.Close()
    print '{} {} closed'.format(mX, mY)

#print len(f_List)
#print len(h_dict)

c.cd()
for hist_index, hist in enumerate(h_dict, 1):
    #print hist, max_dict[hist]
    min_x, max_x, min_y, max_y, x_label, y_label, new_bin_count = hist_list[hist]
    #c.cd(hist_index)
    list_of_hists = h_dict[hist]
    for h,histogram in enumerate(list_of_hists):
        if new_bin_count != None:
            histogram.Rebin(new_bin_count)
        #if hist != 'h_DeltaR_y_vs_ptH':
        #    histogram.SetMaximum(max_dict[hist]*2)
        #    histogram.SetMaximum(max_y)
        #    histogram.SetMinimum(min_dict[hist])
        #    histogram.SetMinimum(
        #histogram.SetMaximum(max_y)
        #histogram.SetMinimum(min_y)
        print min_y, max_y
        if hist != 'h_DeltaR_y_vs_ptH':
            histogram.SetMaximum(max_y)
            histogram.SetMinimum(min_y)
        #histogram.SetAxisRange(min_x, max_x, 'X')
        #histogram.SetAxisRange(min_y, max_y, 'Y')
        histogram.GetXaxis().SetTitle(x_label)
        histogram.GetYaxis().SetTitle(y_label)
        #if h == 0:
        #    histogram.SetTitle("(mX, mY) = ({}, {}) [GeV]".format(*points_list[h]))
        #else:
        histogram.SetTitle("({}, {})".format(*points_list[h]))
        #if hist in ["h_jet_pt_vs_higgs_pt", "h_DeltaR_y_vs_ptH"]:
           # histogram.Draw("colz" if h == 0 else "HIST SAME")
        histogram.Draw("HIST" if h == 0 else "HIST SAME")
        if hist in ["h_jetmass_matched"]:
            histogram.GetXaxis().SetRange(0,600)
        if hist in ["h_higgsmass"]:
            histogram.SetMinimum(0)
            histogram.SetMaximum(14000)
        if hist in ["h_jetmass", "h_jetpt"]:
            histogram.SetMinimum(0)
            histogram.SetMaximum(40000)
        #if hist in ["h_jetpt_matched"]:
            #histogram.SetMinimum(0)
            #histogram.SetMaximum(1000)
    #c.SetLogy(1 if hist == 'h_jetpt' else 0)
    legend = c.BuildLegend(0.65, 0.65, 0.9, 0.9)
    legend.SetHeader("(m_{X}, m_{Y}) [GeV]", "C")
    legend.GetListOfPrimitives().First().SetTextFont(62)
    list_of_hists[0].SetTitle(hist)
    print(hist)
    c.SaveAs(hist + '.pdf')

