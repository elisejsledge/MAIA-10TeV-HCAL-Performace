import pyLCIO
import glob
import ctypes
import math
from optparse import OptionParser
import os
import logging
from ROOT import TH1D, TH2F, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLatex, TLine, TMath, TEfficiency, TF1, TH1, TEfficiency, TGraphErrors, THStack
from ROOT import kBlack, kBlue, kRed, kYellow, kGreen, kGray, kOrange
from ROOT import gStyle, gPad
from ROOT import gROOT
from ROOT import TStyle
from ROOT import TLorentzVector
from ROOT import TH3D, THStack
from optparse import OptionParser
import itertools
from math import *
from array import array
import numpy as np

print("Imports Succesful!")

parser = OptionParser()
parser.add_option("-i", "--inFile",   dest='inFile',
                  default="neutron_study_0_50_v2_ntup_pfoPFO.root", help="Name of the ROOT file")
parser.add_option("-o", "--outFolder",   dest='outFolder',
                  default="transfer_plots_0_50/", help="Name of the output folder")
parser.add_option("-e", "--energyMax", dest='eMax',
                  default="50", help="Maximum energy bin")
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
tree = fFile.Get("jet_tree")

arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_pT = array('d', (0., 10., 15., 20., 25., 30., 35., 40., 45., 55.,))

if(options.eMax == "50"):
    print("Max energy is 50 GeV")
if(options.eMax == "250"):
    print("Max energy is 250 GeV")
    arrBins_pT = array('d', (0., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., 100, 125, 150, 200, 265))

if(options.eMax == "1000"):
    print("Max energy is 1000 GeV")
    arrBins_pT = array('d', (0., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., 100, 125, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1050))


'''creating the transfer function x,y = jet-pT, true-pT
plotting points of the mean of the jet pT in that bin with error bars from the std'''

#1D
h_theta_percent_error = TH1D('theta_percent_error','theta_percent_error', 50, -0.5, 0.5) #resolution for the anti-kt jet

#2D
h_true_pT_vs_jet_pT = TH2D("true_pT_vs_jet_pT", "true_pT_vs_jet_pT", len(arrBins_pT)-1, arrBins_pT, len(arrBins_pT)-1, arrBins_pT)

#3D
h_true_pT_vs_theta_vs_jet_pT = TH3D("true_pT_vs_theta_vs_jet_pT", "true_pT_vs_theta_vs_jet_pT", len(arrBins_pT)-1, arrBins_pT, len(arrBins_theta)-1, arrBins_theta, len(arrBins_pT)-1, arrBins_pT)

histos_list = [h_theta_percent_error, h_true_pT_vs_jet_pT, h_true_pT_vs_theta_vs_jet_pT]

for entry in tree:
    if (0.175 < entry.neutron_theta_angle < 2.96): #exclude nozzle neutrons
        h_theta_percent_error.Fill((entry.neutron_theta_angle - entry.jet_theta)/entry.neutron_theta_angle)
        h_true_pT_vs_theta_vs_jet_pT.Fill(entry.neutron_pT, entry.neutron_theta_angle, entry.jet_pT)
        h_true_pT_vs_jet_pT.Fill(entry.neutron_pT, entry.jet_pT)

#Get the transfer function analysis
subfit_root_file = TFile(options.outFolder + "fits.root", 'RECREATE')

for t_bin in range(0, len(arrBins_theta) - 1):
    transfer_func_title = "transfer_function_{}".format(arrBins_theta[t_bin])
    h_transfer_function = TH1D(transfer_func_title, transfer_func_title, len(arrBins_pT)-1, arrBins_pT)

    for e_bin in range(0, len(arrBins_pT)-1):
        h_pT_proj = h_true_pT_vs_theta_vs_jet_pT.ProjectionX("_pX_{}_{}".format(arrBins_pT[e_bin], arrBins_theta[t_bin]), e_bin, e_bin+1, t_bin, t_bin+1)
        gaussFit = TF1("gaussFit", "gaus")
        h_pT_proj.Fit(gaussFit, "Q0E")

        h_pT_proj.Write()

        transfer_val = gaussFit.GetParameter(1)
        transfer_val_err = gaussFit.GetParameter(2)

        h_transfer_function.SetBinContent(e_bin+1, transfer_val)
        h_transfer_function.SetBinError(e_bin+1, transfer_val_err)
    
    h_transfer_function.Write()

subfit_root_file.Close()



#outputs
output_file = TFile(options.outFolder + "transfer_histos.root", 'RECREATE')
for histo in histos_list:
    histo.Write()
output_file.Close()
fFile.Close()