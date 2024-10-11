import os
import logging
from ROOT import TH1D, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLatex, TLine, TMath, TEfficiency, TF1, TH1, TEfficiency, TGraphErrors, THStack
from ROOT import kBlack, kBlue, kRed, kYellow, kGreen, kGray, kOrange
from ROOT import gStyle, gPad
from ROOT import gROOT
from ROOT import TStyle
from optparse import OptionParser
import itertools
from math import *
from array import array
import numpy as np

print("Imports Succesful!")

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
                  default="neutron_study_0_50_v2_ntup_pfoPFO.root", help="Name of the ROOT file")
parser.add_option("-c", "--calibFile",   dest='calibFile',
                  default="calib_plots_0_50/histos.root", help="Name of the ROOT calibration file")
parser.add_option("-o", "--outFolder",   dest='outFolder',
                  default="v2Calibrated/reso_resp_eff_0_50/", help="Name of the output folder")
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

tree = fFile.Get("neutron_tree")

calibrationFile = TFile(options.calibFile, "READ")
h_hcalCalib_E_theta = calibrationFile.h_jet_response_vs_E_Theta
'''calibrationFile.Get("h_jet_response_vs_E_theta")
print(calibrationFile.GetListOfKeys())
quit()'''
#h_hcalCalib_E_theta = TFile.Open(options.calibFile).GetPrimitive("h_jet_response_vs_E_theta")


arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50.))


if(options.eMax == "50"):
    print("Max energy is 50 GeV")
if(options.eMax == "250"):
    print("Max energy is 250 GeV")
    arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., 100, 125, 150, 200, 250))

if(options.eMax == "1000"):
    print("Max energy is 1000 GeV")
    arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., 100, 125, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000))


#Reco/truth histograms
h_truth_E = TH1D('truth_E', 'truth_E', len(arrBins_E)-1, arrBins_E)
h_truth_theta = TH1D('truth_theta', 'truth_theta', len(arrBins_theta)-1, arrBins_theta) #Rads: true generator theta
h_jet_E = TH1D('jet_E', 'jet_E', len(arrBins_E)-1, arrBins_E)

h_frame = TH1D('framec3', 'framec3', len(arrBins_E)-1, arrBins_E)

#analysis 2D histos
nres_bins = 150
res_max = 1.5
h_jet_error_vs_E = TH2D('jet_error_vs_E','jet_error_vs_E', len(arrBins_E)-1, arrBins_E, nres_bins, -res_max, res_max)
h_jet_error_corrected_vs_E = TH2D('jet_error_corrected_vs_E','jet_error_corrected_vs_E', len(arrBins_E)-1, arrBins_E, nres_bins, -res_max, res_max)

#theta resolution
h_jet_error_vs_theta = TH2D('jet_error_vs_theta','jet_error_vs_theta', len(arrBins_theta)-1, arrBins_theta, nres_bins, -res_max, res_max) #jet theta resolution 2d hist
h_jet_error_corrected_vs_theta = TH2D('jet_error_corrected_vs_theta','jet_error_corrected_vs_theta', len(arrBins_theta)-1, arrBins_theta, nres_bins, -res_max, res_max) #jet theta resolution 2d hist

#analysis 1D histograms
#reso vs E
h_jet_resolution_vs_E = TH1D('jet_resolution_vs_E','jet_resolution_vs_E', len(arrBins_E)-1, arrBins_E) #resolution for the anti-kt jet
h_jet_resolution_corrected_vs_E = TH1D('jet_resolution_corrected_vs_E','jet_resolution_corrected_vs_E', len(arrBins_E)-1, arrBins_E) #resolution for the anti-kt jet
#reso vs theta
h_jet_resolution_vs_theta = TH1D('jet_resolution_vs_theta','jet_resolution_vs_theta', len(arrBins_theta)-1, arrBins_theta) #resolution for the anti-kt jet
h_jet_resolution_corrected_vs_theta = TH1D('jet_resolution_corrected_vs_theta','jet_resolution_corrected_vs_theta', len(arrBins_theta)-1, arrBins_theta) #resolution for the anti-kt jet
#resp vs E
h_jet_response_vs_E = TH1D('jet_response_vs_E','jet_response_vs_E', len(arrBins_E)-1, arrBins_E) #resolution for the anti-kt jet
h_jet_response_corrected_vs_E = TH1D('jet_response_corrected_vs_E','jet_response_corrected_vs_E', len(arrBins_E)-1, arrBins_E) #resolution for the anti-kt jet
#resp vs theta
h_jet_response_vs_theta = TH1D('jet_response_vs_theta','jet_response_vs_theta', len(arrBins_theta)-1, arrBins_theta) #resolution for the anti-kt jet
h_jet_response_corrected_vs_theta = TH1D('jet_response_corrected_vs_theta','jet_response_corrected_vs_theta', len(arrBins_theta)-1, arrBins_theta) #resolution for the anti-kt jet


#efficiency histos
h_truth_E_containing_jet = TH1D('truth_E_containing_jet', 'truth_E_containing_jet', len(arrBins_E)-1, arrBins_E)
h_truth_theta_containing_jet = TH1D('truth_theta_containing_jet','truth_theta_containing_jet', len(arrBins_theta)-1, arrBins_theta)

histos_list = [h_truth_E, h_jet_E,
               h_jet_error_vs_E, h_jet_error_corrected_vs_E, h_jet_error_vs_theta, h_jet_error_corrected_vs_theta,
               h_truth_E_containing_jet, h_truth_theta_containing_jet, h_truth_theta,
               h_jet_resolution_vs_E, h_jet_resolution_corrected_vs_E, h_jet_response_vs_E, h_jet_response_corrected_vs_E,
               h_jet_resolution_vs_theta, h_jet_resolution_corrected_vs_theta, h_jet_response_vs_theta, h_jet_response_corrected_vs_theta
            ]

for entry in tree:
    h_truth_E.Fill(entry.E_truth) #true energy
    h_jet_E.Fill(entry.E) #jet energy
    h_truth_theta.Fill(entry.theta_truth)

    #if(0.175 < entry.theta_truth < 2.96):
    if(entry.theta > 0 and 0.175 < entry.theta_truth < 2.96):
        hcal_calib_value = h_hcalCalib_E_theta.GetBinContent(h_hcalCalib_E_theta.FindBin(entry.E_truth, entry.theta_truth))

        #for jet efficiencies
        h_truth_E_containing_jet.Fill(entry.E_truth)
        h_truth_theta_containing_jet.Fill(entry.theta_truth)

        #Energy error vs E
        h_jet_error_vs_E.Fill(entry.E_truth, (entry.E - entry.E_truth)/entry.E_truth)
        h_jet_error_corrected_vs_E.Fill(entry.E_truth, (entry.E*(1+hcal_calib_value) - entry.E_truth)/entry.E_truth)

        #Energy error vs Theta
        h_jet_error_vs_theta.Fill(entry.theta_truth, (entry.E - entry.E_truth)/entry.E_truth)
        h_jet_error_corrected_vs_theta.Fill(entry.theta_truth, (entry.E*(1+hcal_calib_value) - entry.E_truth)/entry.E_truth)
    

h_jet_efficiency_vs_E = TEfficiency(h_truth_E_containing_jet, h_truth_E)
h_jet_efficiency_vs_theta = TEfficiency(h_truth_theta_containing_jet, h_truth_theta)

h_jet_efficiency_vs_E.SetName("jet_efficiency_vs_E")
h_jet_efficiency_vs_theta.SetName("jet_efficiency_vs_theta")

histos_list.append(h_jet_efficiency_vs_E)
histos_list.append(h_jet_efficiency_vs_theta)


#energy binned analysis
for bin in range(1, len(arrBins_E)-1):
    h_jet_proj = h_jet_error_vs_E.ProjectionY("_py", bin, bin+1)
    h_jet_corrected_proj = h_jet_error_corrected_vs_E.ProjectionY("_pyC", bin, bin+1)

    jetGaussFit = TF1("gaussfit", "gaus")
    jetCorrectedGaussFit = TF1("gaussfit", "gaus")

    h_jet_proj.Fit(jetGaussFit, "E")
    h_jet_corrected_proj.Fit(jetCorrectedGaussFit, "E")

    #calo resolution
    jet_sigma = jetGaussFit.GetParameter(2)
    jet_sigma_err = jetGaussFit.GetParError(2)
    jet_corrected_sigma = jetCorrectedGaussFit.GetParameter(2)
    jet_corrected_sigma_err = jetCorrectedGaussFit.GetParError(2)

    h_jet_resolution_vs_E.SetBinContent(bin+1, jet_sigma)
    h_jet_resolution_vs_E.SetBinError(bin+1, jet_sigma_err)
    h_jet_resolution_corrected_vs_E.SetBinContent(bin+1, jet_corrected_sigma)
    h_jet_resolution_corrected_vs_E.SetBinError(bin+1, jet_corrected_sigma_err)

    #calo response
    jet_mu = jetGaussFit.GetParameter(1)
    jet_mu_err = jetGaussFit.GetParError(1)
    jet_corrected_mu = jetCorrectedGaussFit.GetParameter(1)
    jet_corrected_mu_err = jetCorrectedGaussFit.GetParError(1)

    h_jet_response_vs_E.SetBinContent(bin+1, jet_mu)
    h_jet_response_vs_E.SetBinError(bin+1, jet_mu_err)
    h_jet_response_corrected_vs_E.SetBinContent(bin+1, jet_corrected_mu)
    h_jet_response_corrected_vs_E.SetBinError(bin+1, jet_corrected_mu_err)

    #save subfits
    cdebug = TCanvas("", "", 800, 600)
    cdebug.SetLogy()
    h_jet_corrected_proj.SetTitle("")
    h_jet_corrected_proj.GetYaxis().SetTitle("Number of Events")
    h_jet_corrected_proj.GetYaxis().SetTitleOffset(1.7)
    h_jet_corrected_proj.GetXaxis().SetTitleOffset(1.3)
    h_jet_corrected_proj.GetXaxis().SetTitle("Energy Reconstruction Difference  [E_{reco} - E_{true}]")

    h_jet_corrected_proj.SetLineColor(kBlack)
    h_jet_corrected_proj.SetLineWidth(2)
    h_jet_proj.SetLineColor(kGreen + 3)
    h_jet_proj.SetLineWidth(2)

    h_jet_corrected_proj.Draw("E")
    h_jet_proj.Draw("SAME E")

    leg = init_legend(.58, .74, .9, .89)
    leg.AddEntry(h_jet_proj, "R = 0.4 anti-kt jets", "LE")
    leg.AddEntry(h_jet_corrected_proj, "Jets Corrected", "LE")
    leg.Draw()

    draw_title("Single Neutron Energy Resolution: E_{true} \in  [" + str(int(arrBins_E[bin])) + ", " + str(int(arrBins_E[bin+1])) + "]") 
    cdebug.SaveAs(options.outFolder + "subfits/energyFit" + "_" + str(arrBins_E[bin])+ "_" + str(arrBins_E[bin+1]) + "_debug.pdf")


#theta binned analysis
for bin in range(0, len(arrBins_theta)-1):
    h_jet_proj = h_jet_error_vs_theta.ProjectionY("_py", bin, bin+1)
    h_jet_corrected_proj = h_jet_error_corrected_vs_theta.ProjectionY("_pyC", bin, bin+1)

    jetGaussFit = TF1("gaussfit", "gaus")
    jetCorrectedGaussFit = TF1("gaussfit", "gaus")

    h_jet_proj.Fit(jetGaussFit, "E")
    h_jet_corrected_proj.Fit(jetCorrectedGaussFit, "E")

    #calo resolution
    jet_sigma = jetGaussFit.GetParameter(2)
    jet_sigma_err = jetGaussFit.GetParError(2)
    jet_corrected_sigma = jetCorrectedGaussFit.GetParameter(2)
    jet_corrected_sigma_err = jetCorrectedGaussFit.GetParError(2)

    h_jet_resolution_vs_theta.SetBinContent(bin+1, jet_sigma)
    h_jet_resolution_vs_theta.SetBinError(bin+1, jet_sigma_err)
    h_jet_resolution_corrected_vs_theta.SetBinContent(bin+1, jet_corrected_sigma)
    h_jet_resolution_corrected_vs_theta.SetBinError(bin+1, jet_corrected_sigma_err)

    #calo response
    jet_mu = jetGaussFit.GetParameter(1)
    jet_mu_err = jetGaussFit.GetParError(1)
    jet_corrected_mu = jetCorrectedGaussFit.GetParameter(1)
    jet_corrected_mu_err = jetCorrectedGaussFit.GetParError(1)

    h_jet_response_vs_theta.SetBinContent(bin+1, jet_mu)
    h_jet_response_vs_theta.SetBinError(bin+1, jet_mu_err)
    h_jet_response_corrected_vs_theta.SetBinContent(bin+1, jet_corrected_mu)
    h_jet_response_corrected_vs_theta.SetBinError(bin+1, jet_corrected_mu_err)

    #Save subfits
    cdebug = TCanvas("", "", 800, 600)
    cdebug.SetLogy()
    h_jet_corrected_proj.SetTitle("")
    h_jet_corrected_proj.GetYaxis().SetTitle("Number of Events")
    h_jet_corrected_proj.GetYaxis().SetTitleOffset(1.7)
    h_jet_corrected_proj.GetXaxis().SetTitleOffset(1.3)
    h_jet_corrected_proj.GetXaxis().SetTitle("Energy Reconstruction Difference  [E_{reco} - E_{true}]")

    h_jet_corrected_proj.SetLineColor(kBlack)
    h_jet_corrected_proj.SetLineWidth(2)
    h_jet_proj.SetLineColor(kGreen + 3)
    h_jet_proj.SetLineWidth(2)

    h_jet_corrected_proj.Draw("E")
    h_jet_proj.Draw("SAME E")

    leg = init_legend(.58, .74, .9, .89)
    leg.AddEntry(h_jet_proj, "All R = 0.4 anti-kt jets", "LE")
    leg.AddEntry(h_jet_corrected_proj, "Jets Corrected", "LE")
    leg.Draw()

    draw_title("Single Neutron Energy Resolution: #theta_{true} \in  [" + str(int(arrBins_theta[bin])) + ", " + str(int(arrBins_theta[bin+1])) + "]") 
    cdebug.SaveAs(options.outFolder + "subfits/thetaFit" + "_" + str(arrBins_theta[bin])+ "_" + str(arrBins_theta[bin+1]) + "_debug.pdf")


output_file = TFile(options.outFolder + "histos.root", 'RECREATE')
for histo in histos_list:
    histo.Write()
output_file.Close()
fFile.Close()
calibrationFile.Close()
