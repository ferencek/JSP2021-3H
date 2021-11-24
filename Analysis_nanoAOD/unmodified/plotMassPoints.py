import ROOT as r

#---------------------------------------------------------------------
# to run in the batch mode (to prevent canvases from popping up)
r.gROOT.SetBatch()

# set plot style
r.gROOT.SetStyle("Plain")

# suppress the statistics box
r.gStyle.SetOptStat(0)

# more detailed statistics box
#r.gStyle.SetOptStat(1111111)

# suppress the histogram title
r.gStyle.SetOptTitle(0)

r.gStyle.SetPadTickX(1)  # to get the tick marks on the opposite side of the frame
r.gStyle.SetPadTickY(1)  # to get the tick marks on the opposite side of the frame

# tweak margins
r.gStyle.SetPadTopMargin(0.02);
r.gStyle.SetPadBottomMargin(0.1);
r.gStyle.SetPadLeftMargin(0.1);
r.gStyle.SetPadRightMargin(0.02);

# tweak axis title offsets
r.gStyle.SetTitleOffset(1.25, "Y");

# set nicer fonts
r.gStyle.SetTitleFont(42, "")
r.gStyle.SetTitleFont(42, "XYZ")
r.gStyle.SetLabelFont(42, "XYZ")
r.gStyle.SetTextFont(42)
r.gStyle.SetStatFont(42)
r.gROOT.ForceStyle()
#---------------------------------------------------------------------
plotBP=False
plotCP=True

# create canvas
c = r.TCanvas("c", "",1000,1000)
c.cd()

# create mass plane
bkg = r.TH2F("bkg",";m_{X} [TeV];m_{Y} [TeV]",105,0,4200/1000.,105,0,4200/1000.)
bkg.Draw()


gr_BP = r.TGraph()

n = 0
for (mX, mY) in [(1600, 500), (2000, 300), (2000, 800),  (2500, 300)]:
        gr_BP.SetPoint(n,mX/1000.,mY/1000.)
        n += 1

gr_BP.SetMarkerStyle(4)
gr_BP.SetMarkerColor(r.kBlue)

if plotBP: gr_BP.Draw("SAME P")


gr_CP = r.TGraph()

n = 0
for (mX, mY) in [(4000, 3700), (4000, 1700), (4000, 300),(2600, 1300), (2400, 300), (2000, 1700), (600, 260)]:
        gr_CP.SetPoint(n,mX/1000.,mY/1000.)
        n += 1

gr_CP.SetMarkerStyle(4)
gr_CP.SetMarkerColor(r.kBlue)

if plotCP: gr_CP.Draw("SAME P")


mX_min = 400
mX_max = 4000
mX_step = 200
mY_step = 200

gr = r.TGraph()

n = 0
for mX in range(mX_min, mX_max + mX_step, mX_step):
    for mY in sorted(list(set([260,mX-140])) + range(300, mX-125, mY_step)):
        print ("(mX, mY) = (%i, %i)" % (mX, mY))
        gr.SetPoint(n,mX/1000.,mY/1000.)
        n += 1

gr.SetMarkerStyle(20)
gr.SetMarkerColor(r.kRed)
gr.SetMarkerSize(0.6)

gr.Draw("SAME P")

pline = r.TPolyLine()
pline.SetPoint(0,0,0)
pline.SetPoint(1,4200/1000.,0)
pline.SetPoint(2,4200/1000.,2*125/1000.)
pline.SetPoint(3,3*125/1000.,2*125/1000.)
pline.SetPoint(4,4200/1000.,4200/1000.-125/1000.)
pline.SetPoint(5,4200/1000.,4200/1000.)
pline.SetPoint(6,0,4200/1000.)
pline.SetPoint(7,0,0)
pline.SetFillColor(r.kGray)
#pline.SetFillStyle(3344)
pline.SetFillStyle(3002)
pline.Draw("f")

r.gPad.RedrawAxis()

c.SaveAs("mass_points.pdf")
