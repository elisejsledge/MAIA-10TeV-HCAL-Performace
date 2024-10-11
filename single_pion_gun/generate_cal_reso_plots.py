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
                    "pion_end_phi", "pion_theta_diff", "pion_phi_diff", "ecal_clustered_energy", 
                    "hcal_clustered_energy", "ecal_endcap_clustered_energy", 
                    "hcal_endcap_clustered_energy", "clustered_cal_energy", 
                    "clustered_cal_resolution"]
title_list = ["Energy of Generator Pion", "Transverse Momentum of Generator Pion", 
              "z-direction Momentum of Generator Pion", "Charge of Generator Pion", 
              "Polar Angle of Generator Pion Momentum", "Azimuthal Angle of Generator Pion Momentum", 
              "Polar angle of Generator Pion End Point", "Azimuthal Angle of Generator Pion End Point",
              "Theta Difference", "Phi Difference", 
              "ECAL Barrel Clustered Energy", "HCAL Barrel Clustered Energy", 
              "ECAL Endcap Clustered Energy", "HCAL Endcap Clustered Energy", 
              "Calorimeter Clustered Energy", "Calorimeter Resolution"]
X_axis_list = ["Energy [GeV]", "Transverse Momentum [GeV]", "z-direction Momentum [GeV]", 
               "charge", "Polar Angle (radians)", "Azimuthal Angle (radians)", "Polar Angle (radians)", 
               "Azimuthal Angle (radians)", "Polar Angle (radians)", "Azimuthal Angle (radians)", "Energy [GeV]", 
               "Energy [GeV]", "Energy [GeV]", "Energy [GeV]", 
               "Energy [GeV]", "Resolution"]
Y_axis_list = ["Number of Events", "Number of Events", "Number of Events", "Number of Events",
               "Number of Events", "Number of Events", "Number of Events", "Number of Events",
               "Number of Events", "Number of Events", "Number of Events", "Number of Events", 
               "Number of Events", "Number of Events", "Number of Events", "Number of Events"]

for i in range(len(root_histos_list)):
    #get the file
    h_barrel = fFile.Get(root_histos_list[i])
    
    #set plot parameters
    c1 = TCanvas("", "", 800, 600)
    
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
    leg.AddEntry(h_barrel, caption, "L")
    leg.Draw()

    t2_3 = TLatex()
    t2_3.SetTextFont(42)
    t2_3.SetTextColor(1)
    t2_3.SetTextSize(0.035)
    t2_3.SetTextAlign(12)
    t2_3.SetNDC()
    t2_3.DrawLatex(.62, 0.82, 'B_{solenoid} = 5 T')

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

fFile.Close()