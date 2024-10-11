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

parser = OptionParser()
parser.add_option("-i", "--inFile",   dest='inFile',
                  default="pion_jet_study_0_50.root", help="Name of the ROOT histo file")
parser.add_option("-c", "--caption", dest='caption', 
                  default = "MuColl v0A", help="geometry of simulated detector")
parser.add_option("-o", "--outPrefix",   dest='outPrefix',
                  default="pion_jet_study_0_50", help="Name of the output file(s)")
(options, args) = parser.parse_args()

fFile = TFile(options.inFile, "READ")
caption = options.caption

gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

root_histos_list = ["pion_energy", "pion_pT", "pion_pz", "pion_charge",
                    "pion_momentum_theta", "pion_momentum_phi", "pion_end_theta", 
                    "pion_end_phi", "pion_theta_diff", "pion_phi_diff",
                    "jet_pT", "jet_E", "jet_num", 
                    "jet_const_num", "jet_const_pt", "jet_const_E",
                ]

title_list = ["Energy of Generator Pion", "Transverse Momentum of Generator Pion", 
              "z-direction Momentum of Generator Pion", "Charge of Generator Pion", 
              "Polar Angle of Generator Pion Momentum", "Azimuthal Angle of Generator Pion Momentum", 
              "Polar angle of Generator Pion End Point", "Azimuthal Angle of Generator Pion End Point",
              "Theta Difference", "Phi Difference", "Jet Transverse Momentum", "Jet Energy",
              "Number of Jets per Generator Pion", "Number of Constituent Particles per Jet",
              "Jet Constituent Transverse Momentum", "Jet Constituent Energy",
            ]

X_axis_list = ["Energy [GeV]", "Transverse Momentum [GeV]", "z-direction Momentum [GeV]", 
               "charge", "Polar Angle (radians)", "Azimuthal Angle (radians)", "Polar Angle (radians)", 
               "Azimuthal Angle (radians)", "Polar Angle (radians)", "Azimuthal Angle (radians)",
               "Energy [GeV]", "Energy [GeV]", "Number of Jets", "Number of Jet Constituents",
               "Energy [GeV]",
            ]

Y_axis_list = ["Number of Events", "Number of Events", "Number of Events", "Number of Events",
               "Number of Events", "Number of Events", "Number of Events", "Number of Events",
               "Number of Events", "Number of Events", "Number of Events", "Number of Events",
               "Number of Events", "Number of Events", "Number of Events", "Number of Events",
            ]

TH2D_histos_list = ["pion_jet_pt", "pion_pt_const_pt", 
                    "pion_pt_jet_num", "pion_pt_const_num",
                    ]

TH2D_title_list = ["Jet pT vs Pion pT", "Jet Constituent pT vs Pion pT", 
                   "Number of Jets vs Pion pT", "Number of Jet Constituents vs Pion pT"
                ]

TH2D_X_axis_list = [ "Pion pT [GeV]", "Pion pT [GeV]","Pion pT [GeV]","Pion pT [GeV]"]

TH2D_Y_axis_list = ["Jet pT [GeV]", "Jet Constituent pT [GeV]", "Number of Jets", "Number of Jet Consitutents"]

plots_jet_list = ['pion_jet_reco_eff', 'pion_jet_reco_resolution']
plots_pfo_list = ['pion_pfo_reco_eff', 'pion_pfo_reco_resolution']

plots_title_list = ["Reconstructed Pion Energy Efficiency", "Reconstructed Pion Energy Resolution"]

plots_X_axis_list = ["Energy [GeV]", "Resolution"]

plots_Y_axis_list = ["Efficiency", "Number of Events"]


for i in range(len(root_histos_list[:-2])):
    #get the file
    h_barrel = fFile.Get(root_histos_list[i])

    #set plot parameters
    c1 = TCanvas("", "", 800, 600)
    if(i>9): c1.SetLogy()
    
    h_barrel.SetLineColor(kBlue+1)
    h_barrel.SetLineWidth(2)
    h_barrel.SetTitle("")
    h_barrel.GetYaxis().SetTitle(Y_axis_list[i])
    h_barrel.GetYaxis().SetTitleOffset(1.7)
    h_barrel.GetXaxis().SetNdivisions(10)
    h_barrel.GetXaxis().SetLabelSize(0.04)
    h_barrel.GetXaxis().SetTitleOffset(1.3)
    h_barrel.GetXaxis().SetTitle(X_axis_list[i])
    h_barrel.Draw("HIST")

    gPad.RedrawAxis()

    leg = TLegend(.6, .64, .9, .78)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    #leg.AddEntry(h_barrel, caption, "L")
    leg.Draw()

    t2_3 = TLatex()
    t2_3.SetTextFont(42)
    t2_3.SetTextColor(1)
    t2_3.SetTextSize(0.035)
    t2_3.SetTextAlign(12)
    t2_3.SetNDC()
    #t2_3.DrawLatex(.62, 0.82, 'B_{solenoid} = 5 T')

    t3 = TLatex()
    t3.SetTextFont(42)
    t3.SetTextColor(1)
    t3.SetTextSize(0.035)
    t3.SetTextAlign(12)
    t3.SetNDC()
    t3.DrawLatex(0.24, 0.94, title_list[i])

    t4 = TLatex()
    t4.SetTextFont(42)
    t4.SetTextColor(1)
    t4.SetTextSize(0.035)
    t4.SetTextAlign(12)
    t4.SetNDC()
    t4.DrawLatex(0.82, 0.94, '#sqrt{s} = 10 TeV')

    file_name = options.outPrefix + title_list[i] + ".pdf"
    c1.SaveAs(file_name)

for i in range(len(TH2D_histos_list)):
    #get the file
    h_barrel = fFile.Get(TH2D_histos_list[i])

    #set plot parameters
    c1 = TCanvas("", "", 800, 600)

    h_barrel.SetTitle("")
    h_barrel.GetYaxis().SetTitle(TH2D_Y_axis_list[i])
    h_barrel.GetYaxis().SetTitleOffset(1.7)
    h_barrel.GetXaxis().SetNdivisions(10)
    h_barrel.GetXaxis().SetLabelSize(0.04)
    h_barrel.GetXaxis().SetTitleOffset(1.3)
    h_barrel.GetXaxis().SetTitle(TH2D_X_axis_list[i])
    pt_min = float(options.outPrefix.split("_")[-2])
    pt_max = float(options.outPrefix.split("_")[-1])
    h_barrel.GetXaxis().SetRangeUser(pt_min, pt_max)
    if(i<2): h_barrel.GetYaxis().SetRangeUser(pt_min, pt_max)
    h_barrel.Draw("COLZ")

    gPad.RedrawAxis()

    leg = TLegend(.6, .64, .9, .78)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    #leg.AddEntry(h_barrel, caption, "L")
    leg.Draw()

    t2_3 = TLatex()
    t2_3.SetTextFont(42)
    t2_3.SetTextColor(1)
    t2_3.SetTextSize(0.035)
    t2_3.SetTextAlign(12)
    t2_3.SetNDC()
    #t2_3.DrawLatex(.62, 0.82, 'B_{solenoid} = 5 T')

    t3 = TLatex()
    t3.SetTextFont(42)
    t3.SetTextColor(1)
    t3.SetTextSize(0.035)
    t3.SetTextAlign(12)
    t3.SetNDC()
    t3.DrawLatex(0.24, 0.94, TH2D_title_list[i])

    t4 = TLatex()
    t4.SetTextFont(42)
    t4.SetTextColor(1)
    t4.SetTextSize(0.035)
    t4.SetTextAlign(12)
    t4.SetNDC()
    t4.DrawLatex(0.82, 0.94, '#sqrt{s} = 10 TeV')

    file_name = options.outPrefix + TH2D_title_list[i] + ".pdf"
    c1.SaveAs(file_name)
    

for i in range(len(plots_jet_list)):
    h_barrel = fFile.Get(plots_jet_list[i])
    h_pfo = fFile.Get(plots_pfo_list[i])
    #h_barrel.Sumw2()

    #set plot parameters
    c1 = TCanvas("", "", 800, 600)

    h_barrel.SetTitle("")
    h_barrel.GetYaxis().SetTitle(plots_Y_axis_list[i])
    h_barrel.GetYaxis().SetTitleOffset(1.7)
    h_barrel.GetXaxis().SetNdivisions(10)
    h_barrel.GetXaxis().SetLabelSize(0.04)
    h_barrel.GetXaxis().SetTitleOffset(1.3)
    h_barrel.GetXaxis().SetTitle(plots_X_axis_list[i])

    if(plots_title_list[i] == "Reconstructed Pion Energy Resolution"):
        h_barrel.Rebin(2)
        h_pfo.Rebin(2)

    if(plots_title_list[i] == "Reconstructed Pion Energy Efficiency"):
        pt_min = float(options.outPrefix.split("_")[-2])
        pt_max = float(options.outPrefix.split("_")[-1])
        if(int(pt_min) == 50):
            h_barrel.Scale(0.25)
            h_pfo.Scale(0.25)
            h_barrel.Rebin(4)
            h_pfo.Rebin(4)
        elif (int(pt_min) == 250):
            print("hi")
            h_barrel.Scale(0.1)
            h_pfo.Scale(0.1)
            h_barrel.Rebin(10)
            h_pfo.Rebin(10)
        h_barrel.GetXaxis().SetRangeUser(pt_min, pt_max)

    h_barrel.SetMarkerStyle(20)
    h_barrel.SetMarkerColorAlpha(kBlue, 0.85)
    h_barrel.Draw("HIST P")

    h_pfo.SetMarkerStyle(20)
    h_pfo.SetMarkerColorAlpha(kRed, 0.85)
    h_pfo.Draw("SAME HIST P")

    gPad.RedrawAxis()

    if(plots_title_list[i] == "Reconstructed Pion Energy Efficiency"):
        leg = TLegend(.6, .34, .9, .48)
    if(plots_title_list[i] == "Reconstructed Pion Energy Resolution"):
        leg = TLegend(.6, .64, .9, .78)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.AddEntry(h_barrel, "Jet Reconstruction", "p")
    leg.AddEntry(h_pfo, "PFO Reconstruction", "p")
    leg.Draw()

    t2_3 = TLatex()
    t2_3.SetTextFont(42)
    t2_3.SetTextColor(1)
    t2_3.SetTextSize(0.035)
    t2_3.SetTextAlign(12)
    t2_3.SetNDC()
    #t2_3.DrawLatex(.62, 0.82, 'B_{solenoid} = 5 T')

    t3 = TLatex()
    t3.SetTextFont(42)
    t3.SetTextColor(1)
    t3.SetTextSize(0.035)
    t3.SetTextAlign(12)
    t3.SetNDC()
    t3.DrawLatex(0.24, 0.94, plots_title_list[i])

    t4 = TLatex()
    t4.SetTextFont(42)
    t4.SetTextColor(1)
    t4.SetTextSize(0.035)
    t4.SetTextAlign(12)
    t4.SetNDC()
    t4.DrawLatex(0.82, 0.94, '#sqrt{s} = 10 TeV')

    file_name = options.outPrefix + plots_title_list[i] + ".pdf"
    c1.SaveAs(file_name)




fFile.Close()