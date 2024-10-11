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

'''
Goals:
1) Use the new calibration constant and generate a 2D calibration matrix which is basically the response in E and Theta

To-Do:
1) Replace k_calc with correct value
2) See if jet energy scaling is done properly
3) do with sim level and reco level
4) do with subtracting off ECAL component?
4) see if we need to insure h_sim_E>0
5) save subfigs
6) rescaling jet energy with k
'''

#Define some constants
k_old = 0.0231348530678/0.0004825
k_calc = 49 #Replace this later
k_scaling = k_calc/k_old

# Options
parser = OptionParser()
parser.add_option("-i", "--inFile",   dest='inFile',
                  default="neutron_hcal_calibration_const_0_50.root", help="Name of the ROOT file")
parser.add_option("-j", "--jetInFile",   dest='jetInFile',
                  default="neutron_study_0_50_v2_ntup_pfoPFO.root", help="Name of the ROOT file w/ jet")
parser.add_option("-o", "--outFolder",   dest='outFolder',
                  default="calib_plots_0_50/", help="Name of the output folder")
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
jetFile = TFile(options.jetInFile, "READ")
neutronTree = jetFile.Get("neutron_tree")

fFile = TFile(options.inFile, "READ")
tree = fFile.Get("calorimeter_tree")

arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 10., 15., 20., 25., 30., 35., 40., 45., 50.))

if(options.eMax == "50"):
    print("Max energy is 50 GeV")
if(options.eMax == "250"):
    print("Max energy is 250 GeV")
    arrBins_E = array('d', (0., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., 100, 125, 150, 200, 250))

if(options.eMax == "1000"):
    print("Max energy is 1000 GeV")
    arrBins_E = array('d', (0., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., 100, 125, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000))


#3D Histo
nres_bins = 150
res_max = 1.5
res_bins = np.linspace(-1.5, 1.5, 150)
h_sim_error_vs_E_Theta = TH3D("h_sim_error_vs_E_Theta", "h_sim_error_vs_E_Theta", len(arrBins_E)-1, arrBins_E, len(arrBins_theta)-1, arrBins_theta, len(res_bins) - 1, res_bins)
h_rec_error_vs_E_Theta = TH3D("h_rec_error_vs_E_Theta", "h_rec_error_vs_E_Theta", len(arrBins_E)-1, arrBins_E, len(arrBins_theta)-1, arrBins_theta, len(res_bins) - 1, res_bins)
h_jet_error_vs_E_Theta = TH3D("h_jet_error_vs_E_Theta", "h_jet_error_vs_E_Theta", len(arrBins_E)-1, arrBins_E, len(arrBins_theta)-1, arrBins_theta, len(res_bins) - 1, res_bins)


#2D analysis histos
h_sim_response_vs_E_Theta = TH2D("h_sim_response_vs_E_Theta", "h_sim_response_vs_E_Theta", len(arrBins_E)-1, arrBins_E, len(arrBins_theta)-1, arrBins_theta)
h_rec_response_vs_E_Theta = TH2D("h_rec_response_vs_E_Theta", "h_rec_response_vs_E_Theta", len(arrBins_E)-1, arrBins_E, len(arrBins_theta)-1, arrBins_theta)
h_jet_response_vs_E_Theta = TH2D("h_jet_response_vs_E_Theta", "h_jet_response_vs_E_Theta", len(arrBins_E)-1, arrBins_E, len(arrBins_theta)-1, arrBins_theta)

h_jet_statistics_vs_E_Theta = TH2D("h_jet_statistics_vs_E_Theta", "h_jet_statistics_vs_E_Theta", len(arrBins_E)-1, arrBins_E, len(arrBins_theta)-1, arrBins_theta)


histos_list = [h_sim_error_vs_E_Theta, h_rec_error_vs_E_Theta, h_jet_error_vs_E_Theta, h_sim_response_vs_E_Theta, h_rec_response_vs_E_Theta, h_jet_response_vs_E_Theta, h_jet_statistics_vs_E_Theta]

#ensure hcal_sim_E > 0 maybe
for entry in tree:
    #if entry.HCAL_sim_total
    ecal_E = entry.ECAL_rec_total
    hcal_sim_E = entry.HCAL_sim_total * k_calc
    hcal_rec_E_corrected = entry.HCAL_rec_total * k_scaling

    if (0.175 < entry.theta_truth < 2.96): #exclude nozzle neutrons
        h_sim_error_vs_E_Theta.Fill(entry.E_truth, entry.theta_truth, (entry.E_truth - (hcal_sim_E + ecal_E))/entry.E_truth)
        h_rec_error_vs_E_Theta.Fill(entry.E_truth, entry.theta_truth, (entry.E_truth - (hcal_rec_E_corrected + ecal_E))/entry.E_truth)

for entry in neutronTree:
    jet_E = entry.E #need to rescale this later
    if (0.175 < entry.theta_truth < 2.96 and jet_E > 0): #exclude nozzle neutrons
        h_jet_error_vs_E_Theta.Fill(entry.E_truth, entry.theta_truth, (entry.E_truth - jet_E)/entry.E_truth)

sim_calib_mat = np.full([len(arrBins_E)-1, len(arrBins_theta)-1], -1, float)
rec_calib_mat = np.full([len(arrBins_E)-1, len(arrBins_theta)-1], -1, float)
jet_calib_mat = np.full([len(arrBins_E)-1, len(arrBins_theta)-1], -1, float)

subfit_root_file = TFile(options.outFolder + "fits.root", 'RECREATE')

for e_bin in range(0, len(arrBins_E) - 1):
    for t_bin in range(0, len(arrBins_theta)-1):
        h_sim_proj = h_sim_error_vs_E_Theta.ProjectionZ("_sim_py_{}_{}".format(arrBins_E[e_bin], arrBins_theta[t_bin]), e_bin, e_bin+1, t_bin, t_bin+1)
        h_rec_proj = h_rec_error_vs_E_Theta.ProjectionZ("_rec_py_{}_{}".format(arrBins_E[e_bin], arrBins_theta[t_bin]), e_bin, e_bin+1, t_bin, t_bin+1)
        h_jet_proj = h_jet_error_vs_E_Theta.ProjectionZ("_jet_py_{}_{}".format(arrBins_E[e_bin], arrBins_theta[t_bin]), e_bin, e_bin+1, t_bin, t_bin+1)

        h_sim_proj.Write()
        h_rec_proj.Write()
        h_jet_proj.Write()

        simGaussFit = TF1("simGaussfit", "gaus")
        recGaussFit = TF1("recGaussfit", "gaus")
        jetGaussFit = TF1("jetGaussfit", "gaus")

        h_sim_proj.Fit(simGaussFit, "Q0E")
        h_rec_proj.Fit(recGaussFit, "Q0E")
        h_jet_proj.Fit(jetGaussFit, "Q0E")

        #calibration const calculation
        sim_const = simGaussFit.GetParameter(1)
        rec_const = recGaussFit.GetParameter(1)
        jet_const = jetGaussFit.GetParameter(1)
        sim_const_err = simGaussFit.GetParameter(2)
        rec_const_err = recGaussFit.GetParameter(2)
        jet_const_err = jetGaussFit.GetParameter(2)

        #set values
        h_sim_response_vs_E_Theta.SetBinContent(e_bin+1, t_bin+1, sim_const)
        h_sim_response_vs_E_Theta.SetBinError(e_bin+1, t_bin+1, sim_const_err)

        h_rec_response_vs_E_Theta.SetBinContent(e_bin+1, t_bin+1, rec_const)
        h_rec_response_vs_E_Theta.SetBinError(e_bin+1, t_bin+1, rec_const_err)

        h_jet_response_vs_E_Theta.SetBinContent(e_bin+1, t_bin+1, jet_const)
        h_jet_response_vs_E_Theta.SetBinError(e_bin+1, t_bin+1, jet_const_err)

        #get binning stats
        h_jet_statistics_vs_E_Theta.SetBinContent(e_bin+1, t_bin+1, h_jet_proj.GetEntries())

        #print("sim_const: {}, rec_const: {}".format(sim_const, rec_const))
        sim_calib_mat[e_bin, t_bin] = sim_const
        rec_calib_mat[e_bin, t_bin] = rec_const
        jet_calib_mat[e_bin, t_bin] = jet_const



        #save subfigs
subfit_root_file.Close()

#print(np.array(h_sim_response_vs_E_Theta.GetArray()))
#print(h_rec_response_vs_E_Theta.GetArray())

#sim_calib_mat = h_sim_response_vs_E_Theta.GetArray()
print("SIM MAT -------------------")
print(sim_calib_mat)

print("REC MAT -------------------")
print(rec_calib_mat)

print("JET MAT -------------------")
print(jet_calib_mat)

output_file = TFile(options.outFolder + "histos.root", 'RECREATE')
for histo in histos_list:
    histo.Write()

output_file.Close()
fFile.Close()
jetFile.Close()