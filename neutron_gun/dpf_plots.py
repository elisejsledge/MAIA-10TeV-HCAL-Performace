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

    '''t4 = TLatex()
    t4.SetTextFont(42)
    t4.SetTextColor(1)
    t4.SetTextSize(0.035)
    t4.SetTextAlign(12)
    t4.SetNDC()
    t4.DrawLatex(0.82, 0.94, '#sqrt{s} = 10 TeV')'''

def init_legend(x1, y1, x2, y2):
    leg = TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    return leg

gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., 100, 125, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000))



#v1File = TFile("v1Files/reso_resp_eff_total/histos.root", "READ")
v2File = TFile("v2Files/reso_resp_eff_total/histos.root", "READ")
#TeV3File = TFile("3TeVFiles/reso_resp_eff_total/histos.root", "READ") #3TeV Samples

##########efficiency vs energy
#h_v1_eff_E = v1File.Get("jet_efficiency_vs_E")
h_v2_eff_E = v2File.Get("jet_efficiency_vs_E")
#h_3TeV_eff_E = TeV3File.Get("jet_efficiency_vs_E")


c1 = TCanvas("", "", 800, 600)
h_v2_eff_E.SetTitle("; True Neutron Energy [GeV]; Reconstruction Efficiency [%]")

#h_v2_eff_E.SetLineColor(kBlue)
#h_v1_eff_E.SetLineColor(kRed)
#h_3TeV_eff_E.SetLineColor(kGreen +2)

h_v2_eff_E.Draw()

#h_v1_eff_E.Draw("SAME")
#h_3TeV_eff_E.Draw("SAME")


#gPad.SetLogy()
gPad.SetLogx()

gPad.Update()
h_v2_eff_E.GetPaintedGraph().GetXaxis().SetRangeUser(10,1000)
h_v2_eff_E.GetPaintedGraph().GetYaxis().SetTitleOffset(1.5)
h_v2_eff_E.GetPaintedGraph().GetXaxis().SetTitleOffset(1.4)

leg = init_legend(.65, .25, .97, .40)
#leg.AddEntry(h_v1_eff_E, "10 TeV v1 jets", "LE")
leg.AddEntry(h_v2_eff_E, "Anti-k_{T} R = 0.4", "LE")
#leg.AddEntry(h_3TeV_eff_E, "3 TeV jets", "LE")
leg.Draw()

draw_title("Neutron Reconstruction Efficiency")
c1.SaveAs("v2Files/plots/eff_vs_E.pdf")



##########efficiency vs theta
#h_v1_eff_theta = v1File.Get("jet_efficiency_vs_theta")
h_v2_eff_theta = v2File.Get("jet_efficiency_vs_theta")
#h_3TeV_eff_theta = TeV3File.Get("jet_efficiency_vs_theta")

c2 = TCanvas("", "", 800, 600)
h_v2_eff_theta.SetTitle("; True Neutron Theta [radians]; Reconstruction Efficiency [%]")

#h_v2_eff_theta.SetLineColor(kBlue)
#h_v1_eff_theta.SetLineColor(kRed)
#h_3TeV_eff_theta.SetLineColor(kGreen +2)

h_v2_eff_theta.Draw()
#h_v1_eff_theta.Draw("SAME")
#h_3TeV_eff_theta.Draw("SAME")



gPad.Update()
h_v2_eff_theta.GetPaintedGraph().GetYaxis().SetRangeUser(0.75,1)
h_v2_eff_theta.GetPaintedGraph().GetYaxis().SetTitleOffset(1.5)
h_v2_eff_theta.GetPaintedGraph().GetXaxis().SetTitleOffset(1.4)

leg = init_legend(.65, .20, .97, .35)
#leg.AddEntry(h_v1_eff_theta, "10 TeV v1 jets", "LE")
leg.AddEntry(h_v2_eff_theta, "Anti-k_{T} R = 0.4", "LE")
#leg.AddEntry(h_3TeV_eff_theta, "3 TeV jets", "LE")
leg.Draw()

draw_title("Neutron Reconstruction Efficiency")
c2.SaveAs("v2Files/plots/eff_vs_theta.pdf")


##########resolution vs energy
#h_v1_reso_E = v1File.Get("jet_resolution_vs_E")
h_v2_reso_E = v2File.Get("jet_resolution_vs_E")
#h_3TeV_reso_E = TeV3File.Get("jet_resolution_vs_E")

h_v2_reso_E.SetLineColor(kBlack)
#h_v1_reso_E.SetLineColor(kRed)
#h_3TeV_reso_E.SetLineColor(kGreen +2)

c3 = TCanvas("", "", 800, 600)

h_v2_reso_E.SetTitle(" ")

h_v2_reso_E.GetYaxis().SetTitle(
    "Neutron Energy Resolution   #sigma_{E} / E")
h_v2_reso_E.GetYaxis().SetTitleOffset(1.5)
h_v2_reso_E.GetYaxis().SetRangeUser(-0,.2)
h_v2_reso_E.GetXaxis().SetTitleOffset(1.2)
#h_jet_resolution_vs_E.GetXaxis().SetRangeUser(0., 250.)
h_v2_reso_E.GetXaxis().SetTitle("True Neutron Energy [GeV]")

gPad.SetLogx()

h_v2_reso_E.Draw()
#h_v1_reso_E.Draw("SAME")
#h_3TeV_reso_E.Draw("SAME")

leg = init_legend(.65, .74, .9, .89)
#leg.AddEntry(h_v1_reso_E, "10 TeV v1 jets", "LE")
leg.AddEntry(h_v2_reso_E, "Anti-k_{T} R = 0.4", "LE")
#leg.AddEntry(h_3TeV_reso_E, "3 TeV jets", "LE")
leg.Draw()

draw_title("Neutron Reconstruction Energy Resolution")
c3.SaveAs("v2Files/plots/reso_vs_E.pdf")



##########resolution vs theta
#h_v1_reso_theta = v1File.Get("jet_resolution_vs_theta")
h_v2_reso_theta = v2File.Get("jet_resolution_vs_theta")
#h_3TeV_reso_theta = TeV3File.Get("jet_resolution_vs_theta")

h_v2_reso_theta.SetLineColor(kBlack)
#h_v1_reso_theta.SetLineColor(kRed)
#h_3TeV_reso_theta.SetLineColor(kGreen +2)

c4 = TCanvas("", "", 800, 600)

h_v2_reso_theta.SetTitle(" ")

h_v2_reso_theta.GetYaxis().SetTitle(
    "Neutron Energy Resolution   #sigma_{E} / E")
h_v2_reso_theta.GetYaxis().SetTitleOffset(1.5)
h_v2_reso_theta.GetYaxis().SetRangeUser(-0,.1)
h_v2_reso_theta.GetXaxis().SetTitleOffset(1.2)
#h_jet_resolution_vs_E.GetXaxis().SetRangeUser(0., 250.)
h_v2_reso_theta.GetXaxis().SetTitle("True Neutron Theta [radians]")

h_v2_reso_theta.Draw()
#h_v1_reso_theta.Draw("SAME")
#h_3TeV_reso_theta.Draw("SAME")

#leg = init_legend(.7, .74, .9, .89)
leg = init_legend(.65, .20, .97, .35)
#leg.AddEntry(h_v1_reso_theta, "10 TeV v1 jets", "LE")
leg.AddEntry(h_v2_reso_theta, "Anti-k_{T} R = 0.4", "LE")
#leg.AddEntry(h_3TeV_reso_theta, "3 TeV jets", "LE")
leg.Draw()

draw_title("Neutron Reconstruction Energy Resolution")
c4.SaveAs("v2Files/plots/reso_vs_theta.pdf")


##########response vs energy
#h_v1_resp_E = v1File.Get("jet_response_vs_E")
h_v2_resp_E = v2File.Get("jet_response_vs_E")
#h_3TeV_resp_E = TeV3File.Get("jet_response_vs_E")

h_v2_resp_E.SetLineColor(kBlack)
#h_v1_resp_E.SetLineColor(kRed)
#h_3TeV_resp_E.SetLineColor(kGreen +2)

c5 = TCanvas("", "", 800, 600)

h_v2_resp_E.SetTitle(" ")

h_v2_resp_E.GetYaxis().SetTitle(
    "Neutron Energy Response   #mu_{E} / E")
h_v2_resp_E.GetYaxis().SetTitleOffset(1.5)
h_v2_resp_E.GetYaxis().SetRangeUser(-0.2,0.025)
h_v2_resp_E.GetXaxis().SetTitleOffset(1.2)
h_v2_resp_E.GetXaxis().SetTitle("True Neutron Energy [GeV]")

gPad.SetLogx()

h_v2_resp_E.Draw()
#h_v1_resp_E.Draw("SAME")
#h_3TeV_resp_E.Draw("SAME")

leg = init_legend(.65, .20, .97, .35)
#leg.AddEntry(h_v1_resp_E, "10 TeV v1 jets", "LE")
leg.AddEntry(h_v2_resp_E, "Anti-k_{T} R = 0.4", "LE")
#leg.AddEntry(h_3TeV_resp_E, "3 TeV jets", "LE")
leg.Draw()

draw_title("Neutron Reconstruction Energy Response")
c5.SaveAs("v2Files/plots/resp_vs_E.pdf")


##########response vs theta
#h_v1_resp_theta = v1File.Get("jet_response_vs_theta")
h_v2_resp_theta = v2File.Get("jet_response_vs_theta")
#h_3TeV_resp_theta = TeV3File.Get("jet_response_vs_theta")

h_v2_resp_theta.SetLineColor(kBlack)
'''h_v1_resp_theta.SetLineColor(kRed)
h_3TeV_resp_theta.SetLineColor(kGreen +2)'''

c6 = TCanvas("", "", 800, 600)

h_v2_resp_theta.SetTitle(" ")

h_v2_resp_theta.GetYaxis().SetTitle(
    "Neutron Energy Resolution   #mu_{E} / E")
h_v2_resp_theta.GetYaxis().SetTitleOffset(1.5)
#h_v2_resp_theta.GetYaxis().SetRangeUser(-0,.3)
h_v2_resp_theta.GetXaxis().SetTitleOffset(1.2)
#h_jet_resolution_vs_E.GetXaxis().SetRangeUser(0., 250.)
h_v2_resp_theta.GetYaxis().SetRangeUser(-0.1, 0.025)
h_v2_resp_theta.GetXaxis().SetTitle("True Neutron Theta [radians]")

h_v2_resp_theta.Draw()
#h_v1_resp_theta.Draw("SAME")
#h_3TeV_resp_theta.Draw("SAME")

#leg = init_legend(.7, .5, .9, .65)
leg = init_legend(.65, .74, .9, .89)
#leg.AddEntry(h_v1_resp_theta, "10 TeV v1 jets", "LE")
leg.AddEntry(h_v2_resp_theta, "Anti-k_{T} R = 0.4", "LE")
#leg.AddEntry(h_3TeV_resp_theta, "3 TeV jets", "LE")
leg.Draw()

draw_title("Neutron Reconstruction Energy Resolution")
c6.SaveAs("v2Files/plots/resp_vs_theta.pdf")





















#v1File.Close()
v2File.Close()
#TeV3File.Close()
