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
                  default="plots/", help="Name of the output folder")
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
#arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., 100, 125, 150, 200, 250))

h_truth_E = TH1D('truth_E', 'truth_E', len(arrBins_E)-1, arrBins_E)
h_reco_E = TH1D('reco_E', 'reco_E', len(arrBins_E)-1, arrBins_E)
h_barrel_E = TH1D('barrel_E', 'barrel_E', len(arrBins_E)-1, arrBins_E) #energy deposit in barrel
h_endcap_E = TH1D('endcap_E', 'endcap_E', len(arrBins_E)-1, arrBins_E) #energy deposit in endcap
h_total_depo_E = TH1D('total_depo_E', 'total_depo_E', len(arrBins_E)-1, arrBins_E) #energy deposit in all calorimeters

h_frame = TH1D('framec3', 'framec3', len(arrBins_E)-1, arrBins_E)

#h_resolution = TH2D('resolution_E', 'resolution_E', len(arrBins_E)-1, arrBins_E, 100, -400, 400)
h_resolution = TH2D('resolution_E', 'resolution_E', len(arrBins_E)-1, arrBins_E, 400, -50, 50)
h_barrel_resolution = TH2D('barrel_resolution_E', 'barrel_resolution_E', len(arrBins_E)-1, arrBins_E, 400, -50, 50)
h_endcap_resolution = TH2D('endcap_resolution_E', 'endcap_resolution_E', len(arrBins_E)-1, arrBins_E, 400, -50, 50)
h_total_depo_resolution = TH2D('total_depo_resolution_E', 'total_depo_resolution_E', len(arrBins_E)-1, arrBins_E, 400, -50, 50)

h_barrel_reco_resolution = TH2D('barrel_reco_resolution', 'barrel_reco_resolution', len(arrBins_E)-1, arrBins_E, 400, -50, 50)
h_endcap_reco_resolution = TH2D('endcap_reco_resolution', 'endcap_reco_resolution', len(arrBins_E)-1, arrBins_E, 400, -50, 50)

h_total_resolution = TH1D('total_resolution','total_resolution', 150, -1.5, 1.5) #resolution for the anti-kt jet
h_total_barrel_resolution = TH1D('total_barrel_resolution','total_barrel_resolution', 150, -1.5, 1.5) #resolution for the barrel calos
h_total_endcap_resolution = TH1D('total_endcap_resolution','total_endcap_resolution', 150, -1.5, 1.5) #resolution for the encap calos
h_total_total_depo_resolution = TH1D('total_total_depo_resolution','total_total_depo_resolution', 150, -1.5, 1.5) #resolution for the encap + barrel calos

h_total_barrel_reco_resolution = TH1D('total_barrel_reco_resolution','total_barrel_reco_resolution', 150, -1.5, 1.5) #resolution for the barrel calos using anti-kt
h_total_endcap_reco_resolution = TH1D('total_endcap_reco_resolution','total_endcap_reco_resolution', 150, -1.5, 1.5) #resolution for the endcap calos using anti-kt

h_barrel_frac = TH1D('barrel_frac', 'barrel_frac', 10, 0, 1) #fraction of energy deposited in the barrel
h_endcap_frac = TH1D('endcap_frac', 'endcap_frac', 10, 0, 1) #fraction of depo in endcap

#theta resolution
h_resolution_theta = TH2D('resolution_theta','resolution_theta', len(arrBins_theta)-1, arrBins_theta, 400, -50, 50) #jet theta resolution 2d hist
h_resolution_cone_theta = TH2D('resolution_cone_theta','resolution_cone_theta', len(arrBins_theta)-1, arrBins_theta, 400, -50, 50) #cone theta resolution 2d hist

for entry in tree:
    h_truth_E.Fill(entry.E_truth)
    h_reco_E.Fill(entry.E)

    #if(1 < entry.theta and entry.theta < 2 ): h_resolution.Fill(entry.E_truth, entry.E-entry.E_truth)

    if(entry.theta >= 0): 
        h_resolution.Fill(entry.E_truth, entry.E-entry.E_truth)
        h_resolution_theta.Fill(entry.theta_truth, entry.E-entry.E_truth)
    h_total_resolution.Fill((entry.E-entry.E_truth)/entry.E_truth)

    #my additions
    clearance_fraction = 0.95
    barrel_energy = entry.E_ECal_Barrel + entry.E_HCal_Barrel
    endcap_energy = entry.E_ECal_Endcap + entry.E_HCal_Endcap
    h_barrel_E.Fill(barrel_energy)
    h_endcap_E.Fill(endcap_energy)

    #barrel_frac = barrel_energy/(barrel_energy + endcap_energy)
    #endcap_frac = endcap_energy/(barrel_energy + endcap_energy)
    #h_barrel_frac.Fill(barrel_frac)
    #h_endcap_frac.Fill(endcap_frac)

    if(entry.theta >= 0 and (barrel_energy + endcap_energy) > 0):
        if(barrel_energy/(barrel_energy + endcap_energy) > clearance_fraction): #most energy in the barrel
            h_barrel_resolution.Fill(entry.E_truth, barrel_energy - entry.E_truth)
            h_total_barrel_resolution.Fill((barrel_energy - entry.E_truth)/entry.E_truth)
            h_barrel_reco_resolution.Fill(entry.E_truth, entry.E-entry.E_truth)
            h_total_barrel_reco_resolution.Fill((entry.E-entry.E_truth)/entry.E_truth)

        elif(endcap_energy/(barrel_energy + endcap_energy) > clearance_fraction): #most energy in the endcap
            h_endcap_resolution.Fill(entry.E_truth, endcap_energy - entry.E_truth)
            h_total_endcap_resolution.Fill((endcap_energy - entry.E_truth)/entry.E_truth)
            h_endcap_reco_resolution.Fill(entry.E_truth, entry.E-entry.E_truth)
            h_total_endcap_reco_resolution.Fill((entry.E-entry.E_truth)/entry.E_truth)

    #fill all calo histograms
    h_total_depo_E.Fill(barrel_energy + endcap_energy)
    h_total_depo_resolution.Fill(entry.E_truth, barrel_energy + endcap_energy - entry.E_truth)
    h_total_total_depo_resolution.Fill((endcap_energy + barrel_energy - entry.E_truth)/entry.E_truth)
    h_resolution_cone_theta.Fill(entry.theta_truth, barrel_energy + endcap_energy - entry.E_truth)


######################################################
#Reco vs True energy comparison plot
c1 = TCanvas("", "", 800, 600)

h_truth_E.SetLineColor(kBlue+1)
h_truth_E.SetLineWidth(2)
h_truth_E.SetTitle("")
h_truth_E.GetYaxis().SetTitle("Number of Events")
h_truth_E.GetYaxis().SetTitleOffset(1.7)
h_truth_E.GetXaxis().SetTitleOffset(1.3)
h_truth_E.GetXaxis().SetTitle("Neutron Energy [GeV]")
h_truth_E.Draw("HIST")

h_reco_E.SetLineColor(kRed+1)
h_reco_E.SetLineWidth(2)
h_reco_E.Draw("SAME HIST")

#my additions
h_total_depo_E.SetLineColor(kGreen+3)
h_total_depo_E.SetLineWidth(2)
h_total_depo_E.Draw("SAME HIST")

'''
h_barrel_E.SetLineColor(kOrange+3)
h_barrel_E.SetLineWidth(2)
h_barrel_E.Draw("SAME HIST")

h_endcap_E.SetLineColor(kGray+3)
h_endcap_E.SetLineWidth(2)
h_endcap_E.Draw("SAME HIST")'''

#end additions

#c1.SetLogy()
#h_truth_E.GetYaxis().SetRangeUser(0,12000)

leg = init_legend(.6, .64, .9, .78)
leg.AddEntry(h_truth_E, "True Neutron Energy", "l")
leg.AddEntry(h_reco_E, "Anti-kT Neutron Energy", "l")
leg.AddEntry(h_total_depo_E, "All Calorimeter Hits", "l")
#leg.AddEntry(h_barrel_E, "Barrel Energy", "l")
#leg.AddEntry(h_endcap_E, "Endcap Energy", "l")
leg.Draw()

draw_title("Neutron Energy Reconstruction Comparison")

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
h_total_resolution.GetXaxis().SetTitle("Error  [#frac{E_{reco} - E_{true}}{E_{true}}]")

h_total_resolution.Draw("HIST")
draw_title("Neutron Energy Reco Error With Mismatched Jets")
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
res_err.GetXaxis().SetTitle("Error  [#frac{E_{reco} - E_{true}}{E_{true}}]")


gaussFit = TF1("gaussfit", "gaus", -.9,.9)
#gaussFit.FixParameter(1,0)
res_err.Fit(gaussFit, "", "", -.9,.9)


constant = gaussFit.GetParameter(0)
mean = gaussFit.GetParameter(1)
sigma = gaussFit.GetParameter(2)
chi2 = gaussFit.GetChisquare()
ndof = gaussFit.GetNDF()

res_err.Draw("AP")
res_err.GetYaxis().SetRangeUser(0,2000)

leg = init_legend(.58, .74, .9, .89)
leg.AddEntry(res_err, "R = 0.4 anti-kt jets", "P")
leg.Draw()

draw_title("Single Neutron Reconstruction Energy Error") 
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


#barrel vs endcap vs sum resolution plots
#Total efficiency plot
ctotal = TCanvas("", "", 800, 600)
h_total_total_depo_resolution.SetTitle("")
h_total_total_depo_resolution.GetYaxis().SetTitle("Number of Events")
h_total_total_depo_resolution.GetYaxis().SetTitleOffset(1.7)
h_total_total_depo_resolution.GetXaxis().SetTitleOffset(1.3)
h_total_total_depo_resolution.GetXaxis().SetTitle("Error  [#frac{E_{reco} - E_{true}}{E_{true}}]")

h_total_total_depo_resolution.GetXaxis().SetRangeUser(-1,1)

h_total_barrel_resolution.SetLineColor(kRed)
h_total_endcap_resolution.SetLineColor(kGreen+3)

h_total_barrel_reco_resolution.SetLineColor(kBlack)
h_total_endcap_reco_resolution.SetLineColor(kOrange + 1)

h_total_total_depo_resolution.Draw("HIST")
h_total_barrel_resolution.Draw("SAME HIST")
h_total_endcap_resolution.Draw("SAME HIST")
h_total_barrel_reco_resolution.Draw("SAME HIST")
h_total_endcap_reco_resolution.Draw("SAME HIST")

leg = init_legend(.58, .74, .9, .89)
leg.AddEntry(h_total_total_depo_resolution, "Total Detector Error", "l")
leg.AddEntry(h_total_barrel_resolution, "Total Barrel Error", "l")
leg.AddEntry(h_total_endcap_resolution, "Total Endcap Error", "l")
leg.AddEntry(h_total_barrel_reco_resolution, "Anti-kt Barrel Error", "l")
leg.AddEntry(h_total_endcap_reco_resolution, "Anti-kt Endcap Error", "l")
leg.Draw()

draw_title("Energy Reconstruction Error Comparison")
ctotal.SaveAs(options.outFolder + options.inFile.split('.')[0] + "_calo_resolutions.pdf")


#binned resolution plots
h_barrel_reso_E = TH1D('barrel_reso_E ', 'barrel_reso_E ', len(arrBins_E)-1, arrBins_E) #barrel resolution
h_barrel_resp_E = TH1D('barrel_resp_E ', 'barrel_resp_E ', len(arrBins_E)-1, arrBins_E) #barrel response
h_endcap_reso_E = TH1D('endcap_reso_E ', 'endcap_reso_E ', len(arrBins_E)-1, arrBins_E) #endcap resolution
h_endcap_resp_E = TH1D('endcap_resp_E ', 'endcap_resp_E ', len(arrBins_E)-1, arrBins_E) #endcap response
h_total_depo_reso_E = TH1D('total_depo_reso_E ', 'total_depo_reso_E ', len(arrBins_E)-1, arrBins_E) #all hits resolution
h_total_depo_resp_E = TH1D('total_depo_resp_E ', 'total_depo_resp_E ', len(arrBins_E)-1, arrBins_E) #all hits response

h_reso_E = TH1D('reso_E', 'reso_E', len(arrBins_E)-1, arrBins_E) #all anti-kt jets resolution
h_resp_E = TH1D('resp_E', 'resp_E', len(arrBins_E)-1, arrBins_E) #all anti-ket jets response
for bin in range(0, len(arrBins_E)-1):
    h_barrel_proj = h_barrel_reco_resolution.ProjectionY("_py", bin, bin+1)
    h_endcap_proj = h_endcap_reco_resolution.ProjectionY("_py", bin, bin+1)
    h_total_depo_proj = h_total_depo_resolution.ProjectionY("_py", bin, bin+1)
    h_my_proj = h_resolution.ProjectionY("_py", bin, bin+1)

    barrelGaussFit = TF1("gaussfit", "gaus")
    endcapGaussFit = TF1("gaussfit", "gaus")
    totalDepoGaussFit = TF1("gaussfit", "gaus")
    gaussFit = TF1("gaussfit", "gaus")

    h_barrel_proj.Fit(barrelGaussFit, "E")
    h_endcap_proj.Fit(endcapGaussFit, "E")
    h_total_depo_proj.Fit(totalDepoGaussFit, "E")
    h_my_proj.Fit(gaussFit, "E")

    #calorimeter resolution
    barrel_sigma = barrelGaussFit.GetParameter(2)
    barrel_sigma_err = barrelGaussFit.GetParError(2)
    endcap_sigma = endcapGaussFit.GetParameter(2)
    endcap_sigma_err = endcapGaussFit.GetParError(2)
    total_depo_sigma = totalDepoGaussFit.GetParameter(2)
    total_depo_sigma_err = totalDepoGaussFit.GetParError(2)
    sigma = gaussFit.GetParameter(2)
    sigma_err = gaussFit.GetParError(2)

    h_barrel_reso_E.SetBinContent(bin+1, barrel_sigma/h_barrel_reso_E.GetBinCenter(bin+1))
    h_barrel_reso_E.SetBinError(bin+1, barrel_sigma_err/h_barrel_reso_E.GetBinCenter(bin+1))
    h_endcap_reso_E.SetBinContent(bin+1, endcap_sigma/h_endcap_reso_E.GetBinCenter(bin+1))
    h_endcap_reso_E.SetBinError(bin+1, endcap_sigma_err/h_endcap_reso_E.GetBinCenter(bin+1))
    h_total_depo_reso_E.SetBinContent(bin+1, total_depo_sigma/h_total_depo_reso_E.GetBinCenter(bin+1))
    h_total_depo_reso_E.SetBinError(bin+1, total_depo_sigma_err/h_total_depo_reso_E.GetBinCenter(bin+1))
    h_reso_E.SetBinContent(bin+1, sigma/h_reso_E.GetBinCenter(bin+1))
    h_reso_E.SetBinError(bin+1, sigma_err/h_reso_E.GetBinCenter(bin+1))

    #calorimeter response
    barrel_mu = barrelGaussFit.GetParameter(1)
    barrel_mu_err = barrelGaussFit.GetParError(1)
    endcap_mu = endcapGaussFit.GetParameter(1)
    endcap_mu_err = endcapGaussFit.GetParError(1)
    total_depo_mu = totalDepoGaussFit.GetParameter(1)
    total_depo_mu_err = totalDepoGaussFit.GetParError(1)
    mu = gaussFit.GetParameter(1)
    mu_err = gaussFit.GetParError(1)

    h_barrel_resp_E.SetBinContent(bin+1, barrel_mu/h_barrel_resp_E.GetBinCenter(bin+1))
    h_barrel_resp_E.SetBinError(bin+1, barrel_mu_err/h_barrel_resp_E.GetBinCenter(bin+1))
    h_endcap_resp_E.SetBinContent(bin+1, endcap_mu/h_endcap_resp_E.GetBinCenter(bin+1))
    h_endcap_resp_E.SetBinError(bin+1, endcap_mu_err/h_endcap_resp_E.GetBinCenter(bin+1))
    h_total_depo_resp_E.SetBinContent(bin+1, total_depo_mu/h_total_depo_resp_E.GetBinCenter(bin+1))
    h_total_depo_resp_E.SetBinError(bin+1, total_depo_mu_err/h_total_depo_resp_E.GetBinCenter(bin+1))
    h_resp_E.SetBinContent(bin+1, mu/h_resp_E.GetBinCenter(bin+1))
    h_resp_E.SetBinError(bin+1, mu_err/h_resp_E.GetBinCenter(bin+1))

    #save subfits
    cdebug = TCanvas("", "", 800, 600)
    cdebug.SetLogy()
    h_total_depo_proj.SetTitle("")
    h_total_depo_proj.GetYaxis().SetTitle("Number of Events")
    h_total_depo_proj.GetYaxis().SetTitleOffset(1.7)
    h_total_depo_proj.GetXaxis().SetTitleOffset(1.3)
    h_total_depo_proj.GetXaxis().SetTitle("Energy Reconstruction Difference  [E_{reco} - E_{true}]")

    h_barrel_proj.SetLineColor(kRed)
    h_barrel_proj.SetLineWidth(2)
    h_endcap_proj.SetLineColor(kBlue)
    h_endcap_proj.SetLineWidth(2)
    h_total_depo_proj.SetLineColor(kBlack)
    h_total_depo_proj.SetLineWidth(2)
    h_my_proj.SetLineColor(kGreen + 3)
    h_my_proj.SetLineWidth(2)

    h_total_depo_proj.Draw("E")
    h_barrel_proj.Draw("SAME E")
    h_endcap_proj.Draw("SAME E")
    h_my_proj.Draw("SAME E")


    leg = init_legend(.58, .74, .9, .89)
    leg.AddEntry(h_barrel_proj, "Barrel R = 0.4 anti-kt jets", "LE")
    leg.AddEntry(h_endcap_proj, "Endcap R = 0.4 anti-kt jets", "LE")
    leg.AddEntry(h_my_proj, "All R = 0.4 anti-kt jets", "LE")
    leg.AddEntry(h_total_depo_proj, "Total Calo Hit Energy", "LE")
    leg.Draw()

    draw_title("Single Neutron Energy Resolution: E_{true} \in  [" + str(int(arrBins_E[bin])) + ", " + str(int(arrBins_E[bin+1])) + "]") 
    cdebug.SaveAs(options.outFolder + "/subfits/" + options.inFile.split('.')[0] + "_" + str(arrBins_E[bin])+ "_" + str(arrBins_E[bin+1]) + "_debug.pdf")


#binned resolution plot over theta
h_cone_reso_theta = TH1D('cone_reso_theta ', 'cone_reso_theta ', len(arrBins_theta)-1, arrBins_theta) #all hits resolution
h_cone_resp_theta = TH1D('cone_resp_theta ', 'cone_resp_theta ', len(arrBins_theta)-1, arrBins_theta) #all hits response

h_jet_reso_theta = TH1D('jet_reso_theta', 'jet_reso_theta', len(arrBins_theta)-1, arrBins_theta) #all anti-kt jets resolution
h_jet_resp_theta = TH1D('jet_resp_theta', 'jet_resp_theta', len(arrBins_theta)-1, arrBins_theta) #all anti-ket jets response
for bin in range(0, len(arrBins_theta)-1):
    h_cone_proj = h_resolution_cone_theta.ProjectionY("_py", bin, bin+1)
    h_jet_proj = h_resolution_theta.ProjectionY("_py", bin, bin+1)

    coneGaussFit = TF1("gaussfit", "gaus")
    jetGaussFit = TF1("gaussfit", "gaus")

    h_cone_proj.Fit(coneGaussFit, "E")
    h_jet_proj.Fit(jetGaussFit, "E")

    #calo resolution
    cone_sigma = coneGaussFit.GetParameter(2)
    cone_sigma_err = coneGaussFit.GetParError(2)
    jet_sigma = jetGaussFit.GetParameter(2)
    jet_sigma_err = jetGaussFit.GetParError(2)

    h_cone_reso_theta.SetBinContent(bin+1, cone_sigma/25)
    h_cone_reso_theta.SetBinError(bin+1, cone_sigma_err/25)
    h_jet_reso_theta.SetBinContent(bin+1, jet_sigma/25)
    h_jet_reso_theta.SetBinError(bin+1, jet_sigma_err/25)

    #calo response
    cone_mu = coneGaussFit.GetParameter(1)
    cone_mu_err = coneGaussFit.GetParError(1)
    jet_mu = jetGaussFit.GetParameter(1)
    jet_mu_err = jetGaussFit.GetParError(1)

    h_cone_resp_theta.SetBinContent(bin+1, cone_mu/25)
    h_cone_resp_theta.SetBinError(bin+1, cone_mu_err/25)
    h_jet_resp_theta.SetBinContent(bin+1, jet_mu/25)
    h_jet_resp_theta.SetBinError(bin+1, jet_mu_err/25)

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
    cdebug.SaveAs(options.outFolder + "/subfits/" + options.inFile.split('.')[0] + "_theta_" + str(arrBins_theta[bin])+ "_" + str(arrBins_theta[bin+1]) + "_debug.pdf")






'''
#Binned resolution plot
h_reso_E = TH1D('reso_E', 'reso_E', len(arrBins_E)-1, arrBins_E)
h_resp_E = TH1D('resp_E', 'resp_E', len(arrBins_E)-1, arrBins_E)
for bin in range(0, len(arrBins_E)-1):
    h_my_proj = h_resolution.ProjectionY("_py", bin, bin+1)
    gaussFit = TF1("gaussfit", "gaus")

    #gaussFit.SetParLimits(1,-arrBins_E[bin]/2,0)

    h_my_proj.Fit(gaussFit, "E")# "", -40, 5)

    #calorimeter resolution
    sigma = gaussFit.GetParameter(2)
    sigma_err = gaussFit.GetParError(2)
    #print("Bin: {} | Reso: {}".format(bin+1, sigma/h_reso_E.GetBinCenter(bin+1)))
    h_reso_E.SetBinContent(bin+1, sigma/h_reso_E.GetBinCenter(bin+1))
    h_reso_E.SetBinError(bin+1, sigma_err/h_reso_E.GetBinCenter(bin+1))

    #calorimeter response
    mu = gaussFit.GetParameter(1)
    mu_err = gaussFit.GetParError(1)
    h_resp_E.SetBinContent(bin+1, mu/h_resp_E.GetBinCenter(bin+1))
    h_resp_E.SetBinError(bin+1, mu_err/h_resp_E.GetBinCenter(bin+1))

    cdebug = TCanvas("", "", 800, 600)
    # h_resolutionAlt.Draw("COLZ")

    h_my_proj.SetTitle("")
    h_my_proj.SetLineWidth(2)
    h_my_proj.SetLineColor(kBlack)
    h_my_proj.GetYaxis().SetTitle("Number of Events")
    h_my_proj.GetYaxis().SetTitleOffset(1.7)
    h_my_proj.GetXaxis().SetTitleOffset(1.3)
    h_my_proj.GetXaxis().SetTitle("Energy Reconstruction Difference  [E_{reco} - E_{true}]")

    h_my_proj.Draw("E")

    leg = init_legend(.58, .74, .9, .89)
    leg.AddEntry(h_my_proj, "R = 0.4 anti-kt jets", "LE")
    leg.Draw()

    draw_title("Single Neutron Energy Resolution: E_{true} \in  [" + str(int(arrBins_E[bin])) + ", " + str(int(arrBins_E[bin+1])) + "]") 

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

    cdebug.SaveAs(options.outFolder + "/subfits/" + options.inFile.split('.')[0] + "_" + str(arrBins_E[bin])+ "_" + str(arrBins_E[bin+1]) + "_debug.pdf")'''



c2 = TCanvas("", "", 800, 600)
h_reso_E.SetTitle(" ")

h_reso_E.SetLineColor(kGreen + 3)
h_reso_E.SetLineWidth(2)
h_barrel_reso_E.SetLineColor(kRed)
h_barrel_reso_E.SetLineWidth(2)
h_endcap_reso_E.SetLineColor(kBlue)
h_endcap_reso_E.SetLineWidth(2)
h_total_depo_reso_E.SetLineColor(kBlack)
h_total_depo_reso_E.SetLineWidth(2)


# h_reso_E.SetMaximum(0.035)
# h_reso_E.SetMinimum(0.015)
h_reso_E.GetYaxis().SetTitle(
    "Neutron Energy Resolution   #sigma_{E} / E")
h_reso_E.GetYaxis().SetTitleOffset(1.5)
h_reso_E.GetYaxis().SetRangeUser(-0,.5)
#gPad.SetLogy()
h_reso_E.GetXaxis().SetTitleOffset(1.2)
h_reso_E.GetXaxis().SetRangeUser(0., 250.)
h_reso_E.GetXaxis().SetTitle("Neutron Energy [GeV]")

h_reso_E.Draw("E0")
h_barrel_reso_E.Draw("SAME E0")
h_endcap_reso_E.Draw("SAME E0")
h_total_depo_reso_E.Draw("SAME E0")

leg = init_legend(.58, .74, .9, .89)
leg.AddEntry(h_reso_E, "All R = 0.4 anti-kt jets", "LE")
leg.AddEntry(h_barrel_reso_E, "Barrel jets", "LE")
leg.AddEntry(h_endcap_reso_E, "Endcap jets", "LE")
leg.AddEntry(h_total_depo_reso_E, "Cone Clustering R = 0.1", "LE")
leg.Draw()

draw_title("Reconstructed Neutron Energy Resolution")


c2.SaveAs(options.outFolder + options.inFile.split('.')[0] + "_reso_vs_E.pdf")


c3 = TCanvas("", "", 800, 600)
h_resp_E.SetTitle(" ")
h_resp_E.SetLineColor(kGreen + 3)
h_resp_E.SetLineWidth(2)
h_barrel_resp_E.SetLineColor(kRed)
h_barrel_resp_E.SetLineWidth(2)
h_endcap_resp_E.SetLineColor(kBlue)
h_endcap_resp_E.SetLineWidth(2)
h_total_depo_resp_E.SetLineColor(kBlack)
h_total_depo_resp_E.SetLineWidth(2)


h_resp_E.GetYaxis().SetTitle(
    "Neutron Energy Response   #mu_{E} / E")
h_resp_E.GetYaxis().SetTitleOffset(1.5)
h_resp_E.GetYaxis().SetRangeUser(0,5)
h_resp_E.GetXaxis().SetTitleOffset(1.2)
h_resp_E.GetXaxis().SetRangeUser(0., 250.)
h_resp_E.GetXaxis().SetTitle("Neutron Energy [GeV]")

h_resp_E.Draw("E0")
h_barrel_resp_E.Draw("SAME E0")
h_endcap_resp_E.Draw("SAME E0")
h_total_depo_resp_E.Draw("SAME E0")

h_resp_E.GetYaxis().SetRangeUser(-1,1)

leg = init_legend(.58, .74, .9, .89)
leg.AddEntry(h_resp_E, "All R = 0.4 anti-kt jets", "LE")
leg.AddEntry(h_barrel_resp_E, "Barrel jets", "LE")
leg.AddEntry(h_endcap_resp_E, "Endcap jets", "LE")
leg.AddEntry(h_total_depo_resp_E, "Cone Clustering R = 0.1", "LE")
leg.Draw()

draw_title("Reconstructed Neutron Energy Response")


c3.SaveAs(options.outFolder + options.inFile.split('.')[0] + "_resp_vs_E.pdf")


c4 = TCanvas("", "", 800, 600)
h_barrel_frac.SetTitle("")
h_barrel_frac.SetLineColor(kRed)
h_barrel_frac.SetLineWidth(2)
h_endcap_frac.SetLineColor(kBlue)
h_endcap_frac.SetLineWidth(2)

h_barrel_frac.GetYaxis().SetTitle("Events")
h_barrel_frac.GetYaxis().SetTitleOffset(1.5)
h_barrel_frac.GetYaxis().SetRangeUser(0,5)
h_barrel_frac.GetXaxis().SetTitleOffset(1.2)
h_barrel_frac.GetXaxis().SetRangeUser(0., 100.)
h_barrel_frac.GetXaxis().SetTitle("Percentage of Total Energy Deposition")

h_barrel_frac.Draw("HIST")
#h_endcap_frac.Draw("SAME HIST")

leg = init_legend(.58, .74, .9, .89)
leg.AddEntry(h_barrel_frac, "Barrel Deposition Fraction", "L")
#leg.AddEntry(h_endcap_frac, "Endcap Deposition Fraction", "L")
leg.Draw()

draw_title("depo fracs")

c4.SaveAs(options.outFolder + options.inFile.split('.')[0] + "_depo_fracs.pdf")


#theta resolution
c5 = TCanvas("", "", 800, 600)
h_jet_reso_theta.SetTitle(" ")

h_jet_reso_theta.SetLineColor(kGreen +3)
h_jet_reso_theta.SetLineWidth(2)
h_cone_reso_theta.SetLineColor(kBlack)
h_cone_reso_theta.SetLineWidth(2)

h_jet_reso_theta.GetYaxis().SetTitle(
    "Neutron Energy Resolution   #sigma_{E} / E")
h_jet_reso_theta.GetYaxis().SetTitleOffset(1.5)
h_jet_reso_theta.GetYaxis().SetRangeUser(-0,.5)
#gPad.SetLogy()
h_jet_reso_theta.GetXaxis().SetTitleOffset(1.2)
#h_jet_reso_theta.GetXaxis().SetRangeUser(-3.5, 3.)
h_jet_reso_theta.GetXaxis().SetTitle("True Neutron Theta [radians]")

h_jet_reso_theta.Draw("E0")
h_cone_reso_theta.Draw("SAME E0")

leg = init_legend(.58, .74, .9, .89)
leg.AddEntry(h_jet_reso_theta, "All R = 0.4 anti-kt jets", "LE")
leg.AddEntry(h_cone_reso_theta, "Cone Clustering R = 0.1", "LE")
leg.Draw()

draw_title("Reconstructed Neutron Energy Resolution")

c5.SaveAs(options.outFolder + options.inFile.split('.')[0] + "_reso_vs_theta.pdf")

#theta response
c6 = TCanvas("", "", 800, 600)
h_jet_resp_theta.SetTitle(" ")

h_jet_resp_theta.SetLineColor(kGreen +3)
h_jet_resp_theta.SetLineWidth(2)
h_cone_resp_theta.SetLineColor(kBlack)
h_cone_resp_theta.SetLineWidth(2)

h_jet_resp_theta.GetYaxis().SetTitle(
    "Neutron Energy Response   #mu_{E} / E (using E = 25 avg)")
h_jet_resp_theta.GetYaxis().SetTitleOffset(1.5)
h_jet_resp_theta.GetYaxis().SetRangeUser(-1,1)
#gPad.SetLogy()
h_jet_resp_theta.GetXaxis().SetTitleOffset(1.2)
#h_jet_reso_theta.GetXaxis().SetRangeUser(-3.5, 3.)
h_jet_resp_theta.GetXaxis().SetTitle("True Neutron Theta [radians]")

h_jet_resp_theta.Draw("E0")
h_cone_resp_theta.Draw("SAME E0")

leg = init_legend(.58, .74, .9, .89)
leg.AddEntry(h_jet_resp_theta, "All R = 0.4 anti-kt jets", "LE")
leg.AddEntry(h_cone_resp_theta, "Cone Clustering R = 0.1", "LE")
leg.Draw()

draw_title("Reconstructed Neutron Energy Response")

c6.SaveAs(options.outFolder + options.inFile.split('.')[0] + "_resp_vs_theta.pdf")


fFile.Close()