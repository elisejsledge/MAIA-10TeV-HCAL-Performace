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
                  default="neutron_jet_study_test_ntup_pfoPFO.root", help="Name of the ROOT file")
parser.add_option("-o", "--outFolder",   dest='outFolder',
                  default="reso_resp_eff_plots_0_50/", help="Name of the output folder")
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

arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50.))

'''if(options.inFile == "neutron_study_50_250_ntup_pfoPFO.root"):
    print("The infile is "+ "neutron_study_50_250_ntup_pfoPFO.root")
    arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., 100, 125, 150, 200, 250))

if(options.inFile == "neutron_study_250_1000_ntup_pfoPFO.root"):
    print("The infile is "+ "neutron_study_250_1000_ntup_pfoPFO.root")
    arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., 100, 125, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000))

if(options.inFile == "neutron_study_total.root"):
    print("The infile is "+ "neutron_jet_study_total.root")
    arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., 100, 125, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000))
'''
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
h_cone_E = TH1D('cone_E', 'cone_E', len(arrBins_E)-1, arrBins_E)

h_frame = TH1D('framec3', 'framec3', len(arrBins_E)-1, arrBins_E)

#analysis 2D histos
nres_bins = 150
res_max = 1.5
h_jet_error_vs_E = TH2D('jet_error_vs_E','jet_error_vs_E', len(arrBins_E)-1, arrBins_E, nres_bins, -res_max, res_max)
h_cone_error_vs_E = TH2D('cone_error_vs_E','cone_error_vs_E', len(arrBins_E)-1, arrBins_E, nres_bins, -res_max, res_max)
#theta resolution
h_jet_error_vs_theta = TH2D('jet_error_vs_theta','jet_error_vs_theta', len(arrBins_theta)-1, arrBins_theta, nres_bins, -res_max, res_max) #jet theta resolution 2d hist
h_cone_error_vs_theta = TH2D('cone_error_vs_theta','cone_error_vs_theta', len(arrBins_theta)-1, arrBins_theta, nres_bins, -res_max, res_max) #cone theta resolution 2d hist

#analysis 1D histograms
#reso vs E
h_jet_resolution_vs_E = TH1D('jet_resolution_vs_E','jet_resolution_vs_E', len(arrBins_E)-1, arrBins_E) #resolution for the anti-kt jet
h_cone_resolution_vs_E = TH1D('cone_resolution_vs_E','cone_resolution_vs_E', len(arrBins_E)-1, arrBins_E) #resolution for cone clustering
#reso vs theta
h_jet_resolution_vs_theta = TH1D('jet_resolution_vs_theta','jet_resolution_vs_theta', len(arrBins_theta)-1, arrBins_theta) #resolution for the anti-kt jet
h_cone_resolution_vs_theta = TH1D('cone_resolution_vs_theta','cone_resolution_vs_theta', len(arrBins_theta)-1, arrBins_theta) #resolution for cone clustering
#resp vs E
h_jet_response_vs_E = TH1D('jet_response_vs_E','jet_response_vs_E', len(arrBins_E)-1, arrBins_E) #resolution for the anti-kt jet
h_cone_response_vs_E = TH1D('cone_response_vs_E','cone_response_vs_E', len(arrBins_E)-1, arrBins_E) #resolution for cone clustering
#resp vs theta
h_jet_response_vs_theta = TH1D('jet_response_vs_theta','jet_response_vs_theta', len(arrBins_theta)-1, arrBins_theta) #resolution for the anti-kt jet
h_cone_response_vs_theta = TH1D('cone_response_vs_theta','cone_response_vs_theta', len(arrBins_theta)-1, arrBins_theta) #resolution for cone clustering


#efficiency histos
h_truth_E_containing_jet = TH1D('truth_E_containing_jet', 'truth_E_containing_jet', len(arrBins_E)-1, arrBins_E)
h_truth_theta_containing_jet = TH1D('truth_theta_containing_jet','truth_theta_containing_jet', len(arrBins_theta)-1, arrBins_theta)

histos_list = [h_truth_E, h_jet_E, h_cone_E,
               h_jet_error_vs_E, h_cone_error_vs_E, h_jet_error_vs_theta, h_cone_error_vs_theta,
               h_truth_E_containing_jet, h_truth_theta_containing_jet,
               h_jet_resolution_vs_E, h_cone_resolution_vs_E, h_jet_response_vs_E, h_cone_response_vs_E,
               h_jet_resolution_vs_theta, h_cone_resolution_vs_theta, h_jet_response_vs_theta, h_cone_response_vs_theta
            ]

for entry in tree:
    h_truth_E.Fill(entry.E_truth) #true energy
    h_jet_E.Fill(entry.E) #jet energy
    h_truth_theta.Fill(entry.theta_truth)

    #cone energy
    cone_energy = entry.E_ECal_Barrel + entry.E_HCal_Barrel + entry.E_ECal_Endcap + entry.E_HCal_Endcap
    h_cone_E.Fill(cone_energy)

    if(entry.theta >= 0):
        #for jet efficiencies
        h_truth_E_containing_jet.Fill(entry.E_truth)
        h_truth_theta_containing_jet.Fill(entry.theta_truth)

        #Energy error vs E
        h_jet_error_vs_E.Fill(entry.E_truth, (entry.E - entry.E_truth)/entry.E_truth)
        h_cone_error_vs_E.Fill(entry.E_truth, (cone_energy - entry.E_truth)/entry.E_truth)

        #Energy error vs Theta
        h_jet_error_vs_theta.Fill(entry.theta_truth, (entry.E - entry.E_truth)/entry.E_truth)
        h_cone_error_vs_theta.Fill(entry.theta_truth, (cone_energy - entry.E_truth)/entry.E_truth)
    

h_jet_efficiency_vs_E = TEfficiency(h_truth_E_containing_jet, h_truth_E)
h_jet_efficiency_vs_theta = TEfficiency(h_truth_theta_containing_jet, h_truth_theta)

h_jet_efficiency_vs_E.SetName("jet_efficiency_vs_E")
h_jet_efficiency_vs_theta.SetName("jet_efficiency_vs_theta")

histos_list.append(h_jet_efficiency_vs_E)
histos_list.append(h_jet_efficiency_vs_theta)


#energy binned analysis
for bin in range(1, len(arrBins_E)-1):
    h_jet_proj = h_jet_error_vs_E.ProjectionY("_py", bin, bin+1)
    h_cone_proj = h_cone_error_vs_E.ProjectionY("_py", bin, bin+1)

    jetGaussFit = TF1("gaussfit", "gaus")
    coneGaussFit = TF1("gaussfit", "gaus")

    h_jet_proj.Fit(jetGaussFit, "E")
    h_cone_proj.Fit(coneGaussFit, "E")

    #calo resolution
    cone_sigma = coneGaussFit.GetParameter(2)
    cone_sigma_err = coneGaussFit.GetParError(2)
    jet_sigma = jetGaussFit.GetParameter(2)
    jet_sigma_err = jetGaussFit.GetParError(2)

    h_cone_resolution_vs_E.SetBinContent(bin+1, cone_sigma)
    h_cone_resolution_vs_E.SetBinError(bin+1, cone_sigma_err)
    h_jet_resolution_vs_E.SetBinContent(bin+1, jet_sigma)
    h_jet_resolution_vs_E.SetBinError(bin+1, jet_sigma_err)

    #calo response
    cone_mu = coneGaussFit.GetParameter(1)
    cone_mu_err = coneGaussFit.GetParError(1)
    jet_mu = jetGaussFit.GetParameter(1)
    jet_mu_err = jetGaussFit.GetParError(1)

    h_cone_response_vs_E.SetBinContent(bin+1, cone_mu)
    h_cone_response_vs_E.SetBinError(bin+1, cone_mu_err)
    h_jet_response_vs_E.SetBinContent(bin+1, jet_mu)
    h_jet_response_vs_E.SetBinError(bin+1, jet_mu_err)

    #save subfits
    cdebug = TCanvas("", "", 800, 600)
    cdebug.SetLogy()
    h_cone_proj.SetTitle("")
    h_cone_proj.GetYaxis().SetTitle("Number of Events")
    h_cone_proj.GetYaxis().SetTitleOffset(1.7)
    h_cone_proj.GetXaxis().SetTitleOffset(1.3)
    h_cone_proj.GetXaxis().SetTitle("Energy Reconstruction Difference  [E_{reco} - E_{true}]")

    h_cone_proj.SetLineColor(kBlack)
    h_cone_proj.SetLineWidth(2)
    h_jet_proj.SetLineColor(kGreen + 3)
    h_jet_proj.SetLineWidth(2)

    h_cone_proj.Draw("E")
    h_jet_proj.Draw("SAME E")

    leg = init_legend(.58, .74, .9, .89)
    leg.AddEntry(h_jet_proj, "R = 0.4 anti-kt jets", "LE")
    leg.AddEntry(h_cone_proj, "Cone Clustering", "LE")
    leg.Draw()

    draw_title("Single Neutron Energy Resolution: E_{true} \in  [" + str(int(arrBins_E[bin])) + ", " + str(int(arrBins_E[bin+1])) + "]") 
    cdebug.SaveAs(options.outFolder + "subfits/energyFit" + "_" + str(arrBins_E[bin])+ "_" + str(arrBins_E[bin+1]) + "_debug.pdf")


#theta binned analysis
for bin in range(0, len(arrBins_theta)-1):
    h_jet_proj = h_jet_error_vs_theta.ProjectionY("_py", bin, bin+1)
    h_cone_proj = h_cone_error_vs_theta.ProjectionY("_py", bin, bin+1)

    jetGaussFit = TF1("gaussfit", "gaus")
    coneGaussFit = TF1("gaussfit", "gaus")

    h_jet_proj.Fit(jetGaussFit, "E")
    h_cone_proj.Fit(coneGaussFit, "E")

    #calo resolution
    cone_sigma = coneGaussFit.GetParameter(2)
    cone_sigma_err = coneGaussFit.GetParError(2)
    jet_sigma = jetGaussFit.GetParameter(2)
    jet_sigma_err = jetGaussFit.GetParError(2)

    h_cone_resolution_vs_theta.SetBinContent(bin+1, cone_sigma)
    h_cone_resolution_vs_theta.SetBinError(bin+1, cone_sigma_err)
    h_jet_resolution_vs_theta.SetBinContent(bin+1, jet_sigma)
    h_jet_resolution_vs_theta.SetBinError(bin+1, jet_sigma_err)

    #calo response
    cone_mu = coneGaussFit.GetParameter(1)
    cone_mu_err = coneGaussFit.GetParError(1)
    jet_mu = jetGaussFit.GetParameter(1)
    jet_mu_err = jetGaussFit.GetParError(1)

    h_cone_response_vs_theta.SetBinContent(bin+1, cone_mu)
    h_cone_response_vs_theta.SetBinError(bin+1, cone_mu_err)
    h_jet_response_vs_theta.SetBinContent(bin+1, jet_mu)
    h_jet_response_vs_theta.SetBinError(bin+1, jet_mu_err)

    #Save subfits
    cdebug = TCanvas("", "", 800, 600)
    cdebug.SetLogy()
    h_cone_proj.SetTitle("")
    h_cone_proj.GetYaxis().SetTitle("Number of Events")
    h_cone_proj.GetYaxis().SetTitleOffset(1.7)
    h_cone_proj.GetXaxis().SetTitleOffset(1.3)
    h_cone_proj.GetXaxis().SetTitle("Energy Reconstruction Difference  [E_{reco} - E_{true}]")

    h_cone_proj.SetLineColor(kBlack)
    h_cone_proj.SetLineWidth(2)
    h_jet_proj.SetLineColor(kGreen + 3)
    h_jet_proj.SetLineWidth(2)

    h_cone_proj.Draw("E")
    h_jet_proj.Draw("SAME E")

    leg = init_legend(.58, .74, .9, .89)
    leg.AddEntry(h_jet_proj, "All R = 0.4 anti-kt jets", "LE")
    leg.AddEntry(h_cone_proj, "Cone Clustering", "LE")
    leg.Draw()

    draw_title("Single Neutron Energy Resolution: #theta_{true} \in  [" + str(int(arrBins_theta[bin])) + ", " + str(int(arrBins_theta[bin+1])) + "]") 
    cdebug.SaveAs(options.outFolder + "subfits/thetaFit" + "_" + str(arrBins_theta[bin])+ "_" + str(arrBins_theta[bin+1]) + "_debug.pdf")



'''
#make plots###############################
#jet/cone resolution vs E###############################
c1 = TCanvas("", "", 800, 600)
h_jet_resolution_vs_E.SetTitle(" ")

h_jet_resolution_vs_E.SetLineColor(kGreen + 3)
h_jet_resolution_vs_E.SetLineWidth(2)
h_cone_resolution_vs_E.SetLineColor(kBlack)
h_cone_resolution_vs_E.SetLineWidth(2)

h_jet_resolution_vs_E.GetYaxis().SetTitle(
    "Neutron Energy Resolution   #sigma_{E} / E")
h_jet_resolution_vs_E.GetYaxis().SetTitleOffset(1.5)
h_jet_resolution_vs_E.GetYaxis().SetRangeUser(-0,.3)
h_jet_resolution_vs_E.GetXaxis().SetTitleOffset(1.2)
#h_jet_resolution_vs_E.GetXaxis().SetRangeUser(0., 250.)
h_jet_resolution_vs_E.GetXaxis().SetTitle("Neutron Energy [GeV]")

h_jet_resolution_vs_E.Draw("E0")
h_cone_resolution_vs_E.Draw("SAME E0")

leg = init_legend(.58, .74, .9, .89)
leg.AddEntry(h_jet_resolution_vs_E, "All R = 0.4 anti-kt jets", "LE")
leg.AddEntry(h_cone_resolution_vs_E, "Cone Clustering R = 0.4", "LE")
leg.Draw()

draw_title("Reconstructed Neutron Energy Resolution")
c1.SaveAs(options.outFolder + "reso_vs_E.pdf")

#jet/cone response vs E###############################
c2 = TCanvas("", "", 800, 600)
h_jet_response_vs_E.SetTitle(" ")
h_jet_response_vs_E.SetLineColor(kGreen + 3)
h_jet_response_vs_E.SetLineWidth(2)
h_cone_response_vs_E.SetLineColor(kBlack)
h_cone_response_vs_E.SetLineWidth(2)


h_jet_response_vs_E.GetYaxis().SetTitle(
    "Neutron Energy Response   #mu_{E} / E")
h_jet_response_vs_E.GetYaxis().SetTitleOffset(1.5)
h_jet_response_vs_E.GetYaxis().SetRangeUser(-0.6,0.2)
h_jet_response_vs_E.GetXaxis().SetTitleOffset(1.2)
#h_jet_response_vs_E.GetXaxis().SetRangeUser(0., 250.)
h_jet_response_vs_E.GetXaxis().SetTitle("Neutron Energy [GeV]")

h_jet_response_vs_E.Draw("E0")
h_cone_response_vs_E.Draw("SAME E0")

h_jet_response_vs_E.GetYaxis().SetRangeUser(-1,1)

leg = init_legend(.58, .74, .9, .89)
leg.AddEntry(h_jet_response_vs_E, "All R = 0.4 anti-kt jets", "LE")
leg.AddEntry(h_cone_response_vs_E, "Cone Clustering R = 0.4", "LE")
leg.Draw()

draw_title("Reconstructed Neutron Energy Response")


c2.SaveAs(options.outFolder + "resp_vs_E.pdf")

#jet/cone resolution vs theta###############################
c3 = TCanvas("", "", 800, 600)
h_jet_resolution_vs_theta.SetTitle(" ")

h_jet_resolution_vs_theta.SetLineColor(kGreen +3)
h_jet_resolution_vs_theta.SetLineWidth(2)
h_cone_resolution_vs_theta.SetLineColor(kBlack)
h_cone_resolution_vs_theta.SetLineWidth(2)

h_jet_resolution_vs_theta.GetYaxis().SetTitle(
    "Neutron Energy Resolution   #sigma_{E} / E")
h_jet_resolution_vs_theta.GetYaxis().SetTitleOffset(1.5)
h_jet_resolution_vs_theta.GetYaxis().SetRangeUser(-0,.25)
#gPad.SetLogy()
h_jet_resolution_vs_theta.GetXaxis().SetTitleOffset(1.2)
#h_jet_reso_theta.GetXaxis().SetRangeUser(-3.5, 3.)
h_jet_resolution_vs_theta.GetXaxis().SetTitle("True Neutron Theta [radians]")

h_jet_resolution_vs_theta.Draw("E0")
h_cone_resolution_vs_theta.Draw("SAME E0")

leg = init_legend(.58, .74, .9, .89)
leg.AddEntry(h_jet_resolution_vs_theta, "All R = 0.4 anti-kt jets", "LE")
leg.AddEntry(h_cone_resolution_vs_theta, "Cone Clustering R = 0.4", "LE")
leg.Draw()

draw_title("Reconstructed Neutron Energy Resolution")

c3.SaveAs(options.outFolder + "reso_vs_theta.pdf")


#jet/cone response vs theta###############################
c4 = TCanvas("", "", 800, 600)
h_jet_response_vs_theta.SetTitle(" ")

h_jet_response_vs_theta.SetLineColor(kGreen +3)
h_jet_response_vs_theta.SetLineWidth(2)
h_cone_response_vs_theta.SetLineColor(kBlack)
h_cone_response_vs_theta.SetLineWidth(2)

h_jet_response_vs_theta.GetYaxis().SetTitle(
    "Neutron Energy Response   #mu_{E} / E (using E = 25 avg)")
h_jet_response_vs_theta.GetYaxis().SetTitleOffset(1.5)
h_jet_response_vs_theta.GetYaxis().SetRangeUser(-0.3,0.3)
#gPad.SetLogy()
h_jet_response_vs_theta.GetXaxis().SetTitleOffset(1.2)
#h_jet_reso_theta.GetXaxis().SetRangeUser(-3.5, 3.)
h_jet_response_vs_theta.GetXaxis().SetTitle("True Neutron Theta [radians]")

h_jet_response_vs_theta.Draw("E0")
h_cone_response_vs_theta.Draw("SAME E0")

leg = init_legend(.58, .74, .9, .89)
leg.AddEntry(h_jet_response_vs_theta, "All R = 0.4 anti-kt jets", "LE")
leg.AddEntry(h_cone_response_vs_theta, "Cone Clustering R = 0.4", "LE")
leg.Draw()

draw_title("Reconstructed Neutron Energy Response")

c4.SaveAs(options.outFolder + "resp_vs_theta.pdf")


#jet efficiency vs E###############################
c5 = TCanvas("", "", 800, 600)
h_jet_efficiency_vs_E.SetTitle("; True Neutron Energy [GeV]; Reconstruction Efficiency [%]")
h_jet_efficiency_vs_E.Draw()

draw_title("R = 0.4 Anti-kt Jet Reconstruction Efficiency")
c5.SaveAs(options.outFolder + "eff_vs_E.pdf")


#jet efficiency vs theta###############################
c6 = TCanvas("", "", 800, 600)
h_jet_efficiency_vs_theta.SetTitle("; True Neutron Theta [radians]; Reconstruction Efficiency [%]")
h_jet_efficiency_vs_theta.Draw()

draw_title("R = 0.4 Anti-kt Jet Reconstruction Efficiency")
c6.SaveAs(options.outFolder + "eff_vs_theta.pdf")



'''




output_file = TFile(options.outFolder + "histos.root", 'RECREATE')
for histo in histos_list:
    histo.Write()
output_file.Close()
fFile.Close()