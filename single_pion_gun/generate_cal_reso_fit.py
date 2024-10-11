# -*- coding: utf-8 -*-


import os
import logging
from ROOT import TH1D, TFile, TTree, TColor, TCanvas, TLegend, TLatex, TLine, TF1
from ROOT import kBlack, kBlue, kRed, kYellow, kGray
from ROOT import gStyle, gPad
from ROOT import gROOT
from ROOT import TStyle
from optparse import OptionParser
import itertools
from math import fabs


def check_output_directory(output_path):
    # checks if output directory exists; if not, mkdir
    if not os.path.exists(str(output_path)):
        os.makedirs(output_path)


# Options
parser = OptionParser()
parser.add_option("-i", "--inFile",   dest='inFile',
                  default="histos_cal_resolution_v1.root", help="Name of the ROOT histo file")
parser.add_option("-c", "--caption", dest='caption', 
                  default = "MuColl v1", help="geometry of simulated detector")
parser.add_option("-o", "--outPrefix",   dest='outPrefix',
                  default="histos_cal_reso_", help="Name of the output histo file")
(options, args) = parser.parse_args()

# Init the logger for the script and various modules
fFile = TFile(options.inFile, "READ")
caption = options.caption
outPrefix = options.outPrefix

gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
    
#fit the gaussian to the plot of the difference between reconstructed and generator energies
diff_with_fit = fFile.Get("clustered_cal_resolution")
gaussFit = TF1("gaussfit", "gausn", -0.9, 0.9)
diff_with_fit.Fit(gaussFit, "E")

#set plot parameters
c_diff = TCanvas("", "", 800, 600)
    
diff_with_fit.SetLineColor(kBlue+1)
diff_with_fit.SetLineWidth(2)
diff_with_fit.SetTitle("")
diff_with_fit.GetYaxis().SetTitle("Number of Events")
diff_with_fit.GetYaxis().SetTitleOffset(1.7)
diff_with_fit.GetXaxis().SetNdivisions(10)
diff_with_fit.GetXaxis().SetLabelSize(0.04)
diff_with_fit.GetXaxis().SetTitleOffset(1.3)
diff_with_fit.GetXaxis().SetTitle("Resolution")
diff_with_fit.Draw("HIST")
diff_with_fit.Draw("pe")

gPad.RedrawAxis()

constant = gaussFit.GetParameter(0)
mean = gaussFit.GetParameter(1)
sigma = gaussFit.GetParameter(2)
chi2 = gaussFit.GetChisquare()
ndof = gaussFit.GetNDF()

t2 = TLatex()
t2.SetTextFont(42)
t2.SetTextColor(1)
t2.SetTextSize(0.035)
t2.SetTextAlign(12)
t2.SetNDC()
t2.DrawLatex(.62, 0.8, "Sigma = %.3f"%(sigma))
t2.DrawLatex(.62, 0.75, "Mean = %.3f"%(mean))
t2.DrawLatex(.62, 0.7, "Constant = %.1f"%(constant))
t2.DrawLatex(.62, 0.65, "chi2/ndof = %.1f / %d = %.1f"%(chi2, ndof, chi2/ndof))


t3 = TLatex() 
t3.SetTextFont(42)
t3.SetTextColor(1)
t3.SetTextSize(0.035)
t3.SetTextAlign(12)
t3.SetNDC()
t3.DrawLatex(0.24, 0.94, "(Reconstructed Energy - True Energy)/True Energy")

t4 = TLatex()
t4.SetTextFont(42)
t4.SetTextColor(1)
t4.SetTextSize(0.035)
t4.SetTextAlign(12)
t4.SetNDC()
t4.DrawLatex(0.82, 0.94, '#sqrt{s} = 10 TeV')

file_name = outPrefix + "resolution_with_fit.pdf"
c_diff.SaveAs(file_name)

fFile.Close()
