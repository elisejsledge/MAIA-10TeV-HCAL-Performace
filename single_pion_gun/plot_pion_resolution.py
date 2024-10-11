import os
import logging
from ROOT import TH1D, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLatex, TLine, TMath, TEfficiency, TF1
from ROOT import kBlack, kBlue, kRed, kYellow, kGreen, kGray, kOrange
from ROOT import gStyle, gPad
from ROOT import gROOT
from ROOT import TStyle
from optparse import OptionParser
import itertools
from math import *
from array import array


def check_output_directory(output_path):
    # checks if output directory exists; if not, mkdir
    if not os.path.exists(str(output_path)):
        os.makedirs(output_path)


# Options
parser = OptionParser()
parser.add_option("-i", "--inFile",   dest='inFile',
                  default="ntup_tracks.root", help="Name of the ROOT file")
parser.add_option("-o", "--outFolder",   dest='outFolder',
                  default="/data", help="Name of the output folder")
(options, args) = parser.parse_args()

# Load files
fFile = TFile(options.inFile, "READ")

gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

tree = fFile.Get("pion_tree")

arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 10., 15., 20., 25., 50.,
                  100., 250., 500., 1000.))

c0 = TCanvas("", "", 800, 600)
h_frame = TH1D('framec3', 'framec3', len(arrBins_E)-1, arrBins_E)

h_resolution = TH2D('resolution_E', 'resolution_E', len(
    arrBins_E)-1, arrBins_E, 100, -400, 400)

for entry in tree:
    h_resolution.Fill(entry.E_truth, entry.E-entry.E_truth)

h_reso_E = TH1D('reso_E', 'reso_E', len(arrBins_E)-1, arrBins_E)
for bin in range(0, len(arrBins_E)-1):
    h_my_proj = h_resolution.ProjectionY("_py", bin, bin+1)
    gaussFit = TF1("gaussfit", "gaus")
    h_my_proj.Fit(gaussFit, "E")
    sigma = gaussFit.GetParameter(2)
    sigma_err = gaussFit.GetParError(2)
    h_reso_E.SetBinContent(bin+1, sigma/h_reso_E.GetBinCenter(bin+1))
    h_reso_E.SetBinError(bin+1, sigma_err/h_reso_E.GetBinCenter(bin+1))

    if bin == 2:
        cdebug = TCanvas("", "", 800, 600)
        # h_resolutionAlt.Draw("COLZ")
        h_my_proj.Draw()
        cdebug.SaveAs("debug.pdf")

c2 = TCanvas("", "", 800, 600)
h_reso_E.SetTitle(" ")
h_reso_E.SetLineColor(kBlack)
h_reso_E.SetLineWidth(2)
# h_reso_E.SetMaximum(0.035)
# h_reso_E.SetMinimum(0.015)
h_reso_E.GetYaxis().SetTitle(
    "Pion Energy resolution #sigma_{E} / E")
h_reso_E.GetYaxis().SetTitleOffset(1.5)
h_reso_E.GetXaxis().SetTitleOffset(1.2)
h_reso_E.GetXaxis().SetRangeUser(20., 1000.)
h_reso_E.GetXaxis().SetTitle("Pion Energy [GeV]")
h_reso_E.Draw("E0")

t2_3 = TLatex()
t2_3.SetTextFont(42)
t2_3.SetTextColor(1)
t2_3.SetTextSize(0.035)
t2_3.SetTextAlign(12)
t2_3.SetNDC()
t2_3.DrawLatex(.5, 0.36, 'Pion particle gun')
t2_3.DrawLatex(.5, 0.3, 'Uniform in E (10, 1000) GeV and #theta')

t3 = TLatex()
t3.SetTextFont(42)
t3.SetTextColor(1)
t3.SetTextSize(0.035)
t3.SetTextAlign(12)
t3.SetNDC()
t3.DrawLatex(0.15, 0.94, 'Background hits overlay in [-0.5 10] ns range')

t4 = TLatex()
t4.SetTextFont(42)
t4.SetTextColor(1)
t4.SetTextSize(0.035)
t4.SetTextAlign(12)
t4.SetNDC()
t4.DrawLatex(0.77, 0.94, '#sqrt{s} = 10 TeV')

c2.SaveAs("pion_reso_vs_E.pdf")

fFile.Close()
