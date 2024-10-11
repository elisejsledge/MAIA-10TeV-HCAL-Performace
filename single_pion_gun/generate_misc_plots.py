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
import sys



fFile = TFile(sys.argv[1], "READ")

gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

clus_ene = fFile.Get("clustered_cal_energy")
ecal_ene = fFile.Get("ecal_clustered_energy")
hcal_ene = fFile.Get("hcal_clustered_energy")
ecap_ene = fFile.Get("ecal_endcap_clustered_energy")
hcap_ene = fFile.Get("hcal_endcap_clustered_energy")

c1 = TCanvas("", "", 800, 600)
c1.SetLogy()
clus_ene.SetTitle("Clustered Calorimeter Energy Response")
clus_ene.GetXaxis().SetTitle("Energy Deposition [GeV]")
clus_ene.GetYaxis().SetTitle("Number of Events")
clus_ene.GetXaxis().SetRangeUser(0,float(sys.argv[2]))
clus_ene.SetLineColor(kBlack)
clus_ene.Draw()

leg = TLegend(.6, .64, .9, .78)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.AddEntry(clus_ene,"Total", "L")

ecal_ene.SetLineColor(kBlue-3)
ecal_ene.Draw("SAME")
leg.AddEntry(ecal_ene,"ECAL Barrel", "L")

hcal_ene.SetLineColor(kRed-3)
hcal_ene.Draw("SAME")
leg.AddEntry(hcal_ene,"HCAL Barrel", "L")

ecap_ene.SetLineColor(kBlue-3)
ecap_ene.SetLineStyle(2)
ecap_ene.Draw("SAME")
leg.AddEntry(ecap_ene,"ECAL Endcap", "L")

hcap_ene.SetLineColor(kRed-3)
hcap_ene.SetLineStyle(2)
hcap_ene.Draw("SAME")
leg.AddEntry(hcap_ene,"HCAL Endcap", "L")

leg.Draw()
file_name = sys.argv[1].split(".")[0] + "_cal_response.pdf"
c1.SaveAs(file_name)

