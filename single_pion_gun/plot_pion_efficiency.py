import os
import logging
from ROOT import TH1D, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLatex, TLine, TMath, TEfficiency, TF1, TH1, TEfficiency, TGraphErrors
from ROOT import kBlack, kBlue, kRed, kYellow, kGreen, kGray, kOrange
from ROOT import gStyle, gPad
from ROOT import gROOT
from ROOT import TStyle
from optparse import OptionParser
import itertools
from math import *
from array import array
import numpy as np

def draw_title(title):
    t3 = TLatex()
    t3.SetTextFont(42)
    t3.SetTextColor(1)
    t3.SetTextSize(0.035)
    t3.SetTextAlign(12)
    t3.SetNDC()
    t3.DrawLatex(0.24, 0.94, title)

    t4 = TLatex()
    t4.SetTextFont(42)
    t4.SetTextColor(1)
    t4.SetTextSize(0.035)
    t4.SetTextAlign(12)
    t4.SetNDC()
    t4.DrawLatex(0.82, 0.94, '#sqrt{s} = 10 TeV')

def init_legend(x1, y1, x2, y2):
    leg = TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    return leg

# Options
parser = OptionParser()
parser.add_option("-i", "--inFile",   dest='inFile',
                  default="ntup_tracks.root", help="Name of the ROOT file")
parser.add_option("-o", "--outFolder",   dest='outFolder',
                  default="/data", help="Name of the output folder")
(options, args) = parser.parse_args()

gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

# Load files
fFile = TFile(options.inFile, "READ")

tree = fFile.Get("pion_tree")

arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 10., 20., 30., 40., 50., 75.,
                  100., 125., 150., 175., 200., 225., 250., 300., 350.,
                  400., 450., 500., 600., 700., 800., 900., 1000.))

h_truth_E = TH1D('truth_E', 'truth_E', len(arrBins_E)-1, arrBins_E)
h_reco_E = TH1D('reco_E', 'reco_E', len(arrBins_E)-1, arrBins_E)

h_frame = TH1D('framec3', 'framec3', len(arrBins_E)-1, arrBins_E)

h_resolution = TH2D('resolution_E', 'resolution_E', len(arrBins_E)-1, arrBins_E, 100, -400, 400)
h_total_resolution = TH1D('total_resolution','total_resolution', 150, -1.5, 1.5)

for entry in tree:
    h_truth_E.Fill(entry.E_truth)
    h_reco_E.Fill(entry.E)
    h_resolution.Fill(entry.E_truth, entry.E-entry.E_truth)
    h_total_resolution.Fill((entry.E-entry.E_truth)/entry.E_truth)

######################################################
#Reco vs True energy comparison plot
c1 = TCanvas("", "", 800, 600)

h_truth_E.SetLineColor(kBlue+1)
h_truth_E.SetLineWidth(2)
h_truth_E.SetTitle("")
h_truth_E.GetYaxis().SetTitle("Number of Events")
h_truth_E.GetYaxis().SetTitleOffset(1.7)
h_truth_E.GetXaxis().SetTitleOffset(1.3)
h_truth_E.GetXaxis().SetTitle("Pion Energy [GeV]")
h_truth_E.Draw("HIST")

h_reco_E.SetLineColor(kRed+1)
h_reco_E.SetLineWidth(2)
h_reco_E.Draw("SAME HIST")

leg = init_legend(.6, .64, .9, .78)
leg.AddEntry(h_truth_E, "True Pion Energy", "l")
leg.AddEntry(h_reco_E, "Reco Pion Energy", "l")
leg.Draw()

draw_title("Pion Energy Reconstruction Comparison")

file_name = options.outFolder + options.inFile.split('.')[0] + "_energy_plot.pdf"
c1.SaveAs(file_name)


######################################################
#Efficiency plot adapted from https://gist.github.com/rmanzoni/d31db4d0c692fc5a9580
'''
gStyle.SetStatX(0.90)                
gStyle.SetStatY(0.45)                
gStyle.SetStatW(0.25)                
gStyle.SetStatH(0.15)   

TH1.SetDefaultSumw2()

#bins = np.array(arrBins_E)
#print(bins)

c1 = TCanvas("", "", 800, 600)

eff = TEfficiency(h_reco_E, h_truth_E)
eff.Draw()

file_name = options.outFolder + options.inFile.split('.')[0] + "_eff_plot.pdf"
c1.SaveAs(file_name)'''

######################################################
#Efficiency plot

#Total efficiency plot
ctotal = TCanvas("", "", 800, 600)
h_total_resolution.SetTitle("")
h_total_resolution.GetYaxis().SetTitle("Number of Events")
h_total_resolution.GetYaxis().SetTitleOffset(1.7)
h_total_resolution.GetXaxis().SetTitleOffset(1.3)
h_total_resolution.GetXaxis().SetTitle("Energy Resolution  #frac{E_{reco} - E_{true}}{E_{true}}")

h_total_resolution.Draw("HIST")
draw_title("Pion Energy Resolution With Mismatched Jets")
ctotal.SaveAs(options.outFolder + options.inFile.split('.')[0] + "_resolution.pdf")




#Total efficiency plot with fit
cfit = TCanvas("", "", 800, 600)

#h_total_resolution.Draw("E")
res_err = TGraphErrors(h_total_resolution)
res_err.GetXaxis().SetRangeUser(-0.9,0.9)
res_err.SetMarkerStyle(8)
res_err.SetMarkerSize(0.6)
res_err.SetMarkerColorAlpha(kBlack, 1)
res_err.SetLineColor(kBlack)

res_err.SetTitle("")
res_err.GetYaxis().SetTitle("Number of Events")
res_err.GetYaxis().SetTitleOffset(1.7)
res_err.GetXaxis().SetTitleOffset(1.3)
res_err.GetXaxis().SetTitle("Energy Resolution  #frac{E_{reco} - E_{true}}{E_{true}}")


gaussFit = TF1("gaussfit", "gaus", -.9,.9)
#gaussFit.FixParameter(1,0)
res_err.Fit(gaussFit, "", "", -.1,.1)


constant = gaussFit.GetParameter(0)
mean = gaussFit.GetParameter(1)
sigma = gaussFit.GetParameter(2)
chi2 = gaussFit.GetChisquare()
ndof = gaussFit.GetNDF()

res_err.Draw("AP")

leg = init_legend(.58, .74, .9, .89)
leg.AddEntry(res_err, "R = 0.4 anti-kt jets", "P")
leg.Draw()

draw_title("Single Pion Reconstruction Energy Resolution") 
#{:.3f}".format(gaussFit.GetParameter(1)))

t2 = TLatex()
t2.SetTextFont(42)
t2.SetTextColor(1)
t2.SetTextSize(0.035)
t2.SetTextAlign(12)
t2.SetNDC()
t2.DrawLatex(.66, 0.77, "\sigma = %.3f"%(sigma))
t2.DrawLatex(.66, 0.72, "\mu = %.3f"%(mean))
t2.DrawLatex(.66, 0.67, "Constant = %.1f"%(constant))
t2.DrawLatex(.66, 0.62, "\chi^{2}/ndof = %.1f / %d = %.1f"%(chi2, ndof, chi2/ndof))


cfit.SaveAs(options.outFolder + options.inFile.split('.')[0] + "_resfit.pdf")





#Binned resolution plot
h_reso_E = TH1D('reso_E', 'reso_E', len(arrBins_E)-1, arrBins_E)
for bin in range(0, len(arrBins_E)-1):
    h_my_proj = h_resolution.ProjectionY("_py", bin, bin+1)
    gaussFit = TF1("gaussfit", "gaus")

    h_my_proj.Fit(gaussFit, "E", "", -75, 75)
    sigma = gaussFit.GetParameter(2)
    sigma_err = gaussFit.GetParError(2)
    h_reso_E.SetBinContent(bin+1, sigma/h_reso_E.GetBinCenter(bin+1))
    h_reso_E.SetBinError(bin+1, sigma_err/h_reso_E.GetBinCenter(bin+1))

    cdebug = TCanvas("", "", 800, 600)
    # h_resolutionAlt.Draw("COLZ")

    h_my_proj.SetTitle("")
    h_my_proj.SetLineWidth(2)
    h_my_proj.SetLineColor(kBlack)
    h_my_proj.GetYaxis().SetTitle("Number of Events")
    h_my_proj.GetYaxis().SetTitleOffset(1.7)
    h_my_proj.GetXaxis().SetTitleOffset(1.3)
    h_my_proj.GetXaxis().SetTitle("Energy Resolution  E_{reco} - E_{true}")

    h_my_proj.Draw("E")

    leg = init_legend(.58, .74, .9, .89)
    leg.AddEntry(h_my_proj, "R = 0.4 anti-kt jets", "LE")
    leg.Draw()

    draw_title("Single Pion Energy Resolution: E_{true} \in  [" + str(int(arrBins_E[bin])) + ", " + str(int(arrBins_E[bin+1])) + "]") 

    constant = gaussFit.GetParameter(0)
    mean = gaussFit.GetParameter(1)
    sigma = gaussFit.GetParameter(2)
    mean_err = gaussFit.GetParError(1)
    sigma_err = gaussFit.GetParError(2)

    t2 = TLatex()
    t2.SetTextFont(42)
    t2.SetTextColor(1)
    t2.SetTextSize(0.035)
    t2.SetTextAlign(12)
    t2.SetNDC()
    t2.DrawLatex(.66, 0.77, "#sigma_{E} = %.3f \pm %.3f"%(sigma,sigma_err))
    t2.DrawLatex(.66, 0.72, "\mu = %.3f \pm %.3f"%(mean,mean_err))

    cdebug.SaveAs(options.outFolder + "subfits/" + options.inFile.split('.')[0] + "_" + str(arrBins_E[bin])+ "_" + str(arrBins_E[bin+1]) + "_debug.pdf")


c2 = TCanvas("", "", 800, 600)
h_reso_E.SetTitle(" ")
h_reso_E.SetLineColor(kBlack)
h_reso_E.SetLineWidth(2)
# h_reso_E.SetMaximum(0.035)
# h_reso_E.SetMinimum(0.015)
h_reso_E.GetYaxis().SetTitle(
    "Pion Energy Resolution   #sigma_{E} / E")
h_reso_E.GetYaxis().SetTitleOffset(1.5)
h_reso_E.GetYaxis().SetRangeUser(0,.25)
h_reso_E.GetXaxis().SetTitleOffset(1.2)
h_reso_E.GetXaxis().SetRangeUser(20., 1000.)
h_reso_E.GetXaxis().SetTitle("Pion Energy [GeV]")
h_reso_E.Draw("E0")

leg = init_legend(.58, .74, .9, .89)
leg.AddEntry(h_reso_E, "R = 0.4 anti-kt jets", "LE")
leg.Draw()

draw_title("Reconstructed Pion Energy Resolution")


c2.SaveAs(options.outFolder + options.inFile.split('.')[0] + "_reso_vs_E.pdf")


fFile.Close()