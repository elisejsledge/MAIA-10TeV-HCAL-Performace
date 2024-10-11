import pyLCIO
import glob
import ctypes
import math
from optparse import OptionParser
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

print("Imports succesful!")

exec(open("./plotHelper.py").read())
print("Plot helper opened")

# ############## SETUP #############################

# Prevent ROOT from drawing while you're running -- good for slow remote servers
# Instead, save files and view them with an sftp client like Fetch (feel free to ask me for my UTK license)
ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1
obj_type = "ne"
append = "240417" #idk what this is
do_ECAL = False
do_HCAL = True
ECAL_calib_file = "calib/theta_prof_240415.root" #make sure this isn't used later

e_calibration_mip = 0.0001575
e_calibration_mip_to_reco = 0.00641222630095
e_sampling_scaling = e_calibration_mip_to_reco/e_calibration_mip

h_calibration_mip = 0.0004825
h_calibration_mip_to_reco = 0.0231348530678
h_sampling_scaling = h_calibration_mip_to_reco/h_calibration_mip


# Set up things for each object
settings = {
        "labelname": {  "ne": "Neutron",
                        "pi": "Pion",
                        "ph": "Photon",
                        "mu": "Muon",
                        "el": "Electron"},
        "plotdir":{ "ne": "neutrons",
                    "pi": "pions",
                    "ph": "photons",
                    "mu": "muons",
                    "el": "electrons"},
        "pdgid":  { "ne": [2112],
                    "pi": [211, 111],
                    "ph": [22],
                    "mu": [13],
                    "el": [11]},
        "mass":   { "ne": 0.940,
                    "pi": 0.135,
                    "ph": 0,
                    "mu": 0.106,
                    "el": 0.000511}
}
print("Running on", settings["labelname"][obj_type])

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

def getCollection(event, name):
    if name in event.getCollectionNames():
        return event.getCollection(name)
    return []

#Adapted from Tova Holmes, Federico Meloni, and Junjia Zhang
#########################
# parameters
parser = OptionParser()
parser.add_option('-i', '--inFileDir', help='--inFileDir /data/fmeloni/DataMuC_MuColl_v1/reco/',
                  type=str, default='/data/fmeloni/DataMuC_MuColl10_v0A/v2/reco/')
parser.add_option('-p', '--inFilePrefix', help = '--inFilePrefix neutronGun_E_0_50', 
                  type=str, default='neutronGun_E_0_50')
parser.add_option('-o', '--outFile', help='--output file name', 
                  type=str, default = 'neutron_calibration_test_0_50')
parser.add_option('-v', action="store_true", help="-v if verbose", default=False)

(options, args) = parser.parse_args()

gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

# Define good particle
def isGood(tlv):
    if abs(tlv.Eta()) < 2.436:
        return True
    return False

# Perform matching between two TLVs
def isMatched(tlv1, tlv2, req_pt = True):
    if tlv1.DeltaR(tlv2) > 0.1: return False
    if req_pt:
        drelpt = abs(tlv1.Perp()-tlv2.Perp())/tlv2.Perp()
        if drelpt > 0.1*tlv2.Perp()/100: return False # Require 10% at 100, 20% at 200, ...
    return True

def getClusterEta(cluster):
    theta = cluster.getITheta()
    return -1*math.ln(math.tan(theta/2))


calibration_tree = TTree("calibration_tree", "calibration_tree")
E_truth = array('d', [0]) #true neutron energy
E_ECAL_barrel_sim = array('d', [0]) #sim energy in ecal barrel
E_ECAL_endcap_sim = array('d', [0]) #sim energy in ecal endcap
E_ECAL_sim_total = array('d', [0]) #total sim energy in ecal 

E_HCAL_barrel_sim = array('d', [0]) #sim energy in hcal barrel
E_HCAL_endcap_sim = array('d', [0]) #sim energy in hcal endcap
E_HCAL_sim_total = array('d', [0]) #total sim energy in ecal

#E_ECAL_barrel_digi = array('d', [0]) #digi energy in ecal barrel
#E_ECAL_endcap_digi = array('d', [0]) #digi energy in ecal endcap
E_ECAL_digi_total = array('d', [0]) #total digi energy in ecal 

#E_HCAL_barrel_digi = array('d', [0]) #digi energy in hcal barrel
#E_HCAL_endcap_digi = array('d', [0]) #digi energy in hcal endcap
E_HCAL_digi_total = array('d', [0]) #total digi energy in hcal

#E_ECAL_barrel_mip_to_digi = array('d', [0]) #converted digi energy in ecal barrel
#E_ECAL_endcap_mip_to_digi = array('d', [0]) #converted digi energy in ecal endcap
E_ECAL_mip_to_digi_total = array('d', [0]) #converted total digi energy in ecal 

#E_HCAL_barrel_mip_to_digi = array('d', [0]) #converted digi energy in hcal barrel
#E_HCAL_endcap_mip_to_digi = array('d', [0]) #converted digi energy in hcal endcap
E_HCAL_mip_to_digi = array('d', [0]) #total converted digi energy in hcal

E_ECAL_total_corrected = array('d', [0]) #corrected ECAL energy
E_HCAL_total_corrected = array('d', [0]) #corrected HCAL energy
E_ECAL_HCAL_total_corrected = array('d', [0]) #corrected energy of all calos

calibration_tree.Branch("E_truth",  E_truth,  'var/D')
calibration_tree.Branch("E_ECAL_barrel_sim",  E_ECAL_barrel_sim,  'var/D')
calibration_tree.Branch("E_ECAL_endcap_sim",  E_ECAL_endcap_sim,  'var/D')
calibration_tree.Branch("E_ECAL_sim_total",  E_ECAL_sim_total,  'var/D')

calibration_tree.Branch("E_HCAL_barrel_sim",  E_HCAL_barrel_sim,  'var/D')
calibration_tree.Branch("E_HCAL_endcap_sim",  E_HCAL_endcap_sim,  'var/D')
calibration_tree.Branch("E_HCAL_sim_total",  E_HCAL_sim_total,  'var/D')

#calibration_tree.Branch("E_ECAL_barrel_digi",  E_ECAL_barrel_digi,  'var/D')
#calibration_tree.Branch("E_ECAL_endcap_digi",  E_ECAL_endcap_digi,  'var/D')
calibration_tree.Branch("E_ECAL_digi_total",  E_ECAL_digi_total,  'var/D')

#calibration_tree.Branch("E_HCAL_barrel_digi",  E_HCAL_barrel_digi,  'var/D')
#calibration_tree.Branch("E_HCAL_endcap_digi",  E_HCAL_endcap_digi,  'var/D')
calibration_tree.Branch("E_HCAL_digi_total",  E_HCAL_digi_total,  'var/D')

#calibration_tree.Branch("E_ECAL_barrel_mip_to_digi",  E_ECAL_barrel_mip_to_digi,  'var/D')
#calibration_tree.Branch("E_ECAL_endcap_mip_to_digi",  E_ECAL_endcap_mip_to_digi,  'var/D')
calibration_tree.Branch("E_ECAL_mip_to_digi_total",  E_ECAL_mip_to_digi_total,  'var/D')

#calibration_tree.Branch("E_HCAL_barrel_mip_to_digi",  E_HCAL_barrel_mip_to_digi,  'var/D')
#calibration_tree.Branch("E_HCAL_endcap_mip_to_digi",  E_HCAL_endcap_mip_to_digi,  'var/D')
calibration_tree.Branch("E_HCAL_mip_to_digi",  E_HCAL_mip_to_digi,  'var/D')

calibration_tree.Branch("E_ECAL_total_corrected",  E_ECAL_total_corrected,  'var/D')
calibration_tree.Branch("E_HCAL_total_corrected",  E_HCAL_total_corrected,  'var/D')
calibration_tree.Branch("E_ECAL_HCAL_total_corrected",  E_ECAL_HCAL_total_corrected,  'var/D')



# ############## CREATE EMPTY HISTOGRAM OBJECTS  #############################
# Set up histograms
# This is an algorithmic way of making a bunch of histograms and storing them in a dictionary
variables = {}
#variables["z"] = {"nbins":500, "xmin": -10000, "xmax": 10000, "title": "z [mm]"}
variables["scale"] = {"nbins":80000, "xmin": 0, "xmax": 16000, "title": "true-over-sim ratio"}
#variables["pt"] =  {"nbins": 30, "xmin": 0, "xmax": 3000,   "title": "p_{T} [GeV]"}
#variables["E"] =   {"nbins": 50, "xmin": 0, "xmax": 1000,   "title": "E [GeV]"}
variables["eta"] = {"nbins": 48, "xmin": -2.4, "xmax": 2.4,     "title": "#eta"}
variables["theta"] = {"nbins": 63, "xmin": 0, "xmax": 3.15,     "title": "#theta"}
#variables["phi"] = {"nbins": 30, "xmin": -3.5, "xmax": 3.5, "title": "#phi"}
#variables["n"] =   {"nbins": 20, "xmin": 0, "xmax": 20,     "title": "n"}
hists = {}

# Initialize all the 2D histograms: the each of the above variables at each level vs the mcp value
hists2d = {}
for var in variables:
    if var=="scale": continue
    hists2d[var] = ROOT.TH2F(var, var, variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"], variables["scale"]["nbins"], variables["scale"]["xmin"], variables["scale"]["xmax"])

# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
# Loop over events
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "ECalBarrelCollection", "ECalEndcapCollection", "HCalBarrelCollection", "HCalEndcapCollection",
                                "EcalBarrelCollectionDigi", "EcalEndcapCollectionDigi", "HcalBarrelCollectionDigi", "HcalEndcapCollectionDigi"])

#create an array of the slcio files to process
to_process = []
file_dir = options.inFileDir + options.inFilePrefix
#print(file_dir)
if os.path.isdir(file_dir):
    for file_name in os.listdir(file_dir): 
        to_process.append(file_dir + "/" + file_name)
        #print(to_process[-1])

i = 0 #count number of files processed
for file in to_process:
    i+=1
    #if (i==11): break
    reader.open(file)
    #print(file)
    for ievt, event in enumerate(reader):
        if options.v: print("Processing event " + str(ievt))
        ###alex's code tests
        '''ecals_b = getCollection(event, "EcalBarrelCollectionDigi")
        ecals_e = getCollection(event, "EcalEndcapCollectionDigi")
        hcals_b = getCollection(event, "HcalBarrelCollectionDigi")
        hcals_e = getCollection(event, "HcalEndcapCollectionDigi")
        print(f"Event {ievt}")
        print(f"len(ecals_b): {len(ecals_b)}")
        print(f"len(ecals_e): {len(ecals_e)}")
        print(f"len(hcals_b): {len(hcals_b)}")
        print(f"len(hcals_e): {len(hcals_e)}")
        break'''

        # Get the collections we care about
        mcpCollection = event.getCollection("MCParticle")
        ecal_simCollection_b = event.getCollection("ECalBarrelCollection")
        ecal_simCollection_e = event.getCollection("ECalEndcapCollection")
        hcal_simCollection_b = event.getCollection("HCalBarrelCollection")
        hcal_simCollection_e = event.getCollection("HCalEndcapCollection")


        try: ecal_digCollection_b = event.getCollection("EcalBarrelCollectionDigi")
        except: ecal_digCollection_b = None
        try:ecal_digCollection_e = event.getCollection("EcalEndcapCollectionDigi")
        except: ecal_digCollection_e = None

        try: hcal_digCollection_b = event.getCollection("HcalBarrelCollectionDigi")
        except: hcal_digCollection_b = None
        try: hcal_digCollection_e = event.getCollection("HcalEndcapCollectionDigi")
        except: hcal_digCollection_e = None


        '''
        try: 
            ecal_recCollection_b = event.getCollection("EcalBarrelCollectionRec")
            print("EcalBarrelCollectionRec not empty!")
        except: ecal_digCollection_b = None
        try: 
            ecal_recCollection_e = event.getCollection("EcalEndcapCollectionRec")
            print("EcalEndcapCollectionRec not empty!")
        except: ecal_digCollection_e = None

        try: 
            hcal_recCollection_b = event.getCollection("HcalBarrelCollectionRec")
            print("HcalBarrelCollectionRec not empty!")
        except: hcal_digCollection_b = None
        try: 
            hcal_digCollection_e = event.getCollection("HcalEndcapCollectionRec")
            print("HcalEndcapCollectionRec not empty!")
        except: hcal_digCollection_e = None'''

        # Make counter variables
        n_mcp_ob = 0
        has_mcp_ob = False
        my_mcp_ob = None


        # Loop over the truth objects and fill histograms
        for mcp in mcpCollection:
            mcp_tlv = getTLV(mcp)

            if abs(mcp.getPDG()) in settings['pdgid'][obj_type] and mcp.getGeneratorStatus()==1 and isGood(mcp_tlv):
                has_mcp_ob = True
                n_mcp_ob += 1
                my_mcp_ob = mcp_tlv

        # Figure out total sim energy
        e_values_b = [x.getEnergy() for x in ecal_simCollection_b]
        e_values_e = [x.getEnergy() for x in ecal_simCollection_e]
        e_sim_E = sum(e_values_b) + sum(e_values_e)

        '''
        e_values_digi_b = [x.getEnergy() for x in ecal_digCollection_b]
        e_values_digi_e = [x.getEnergy() for x in ecal_digCollection_e]
        e_digi_E = sum(e_values_digi_b) + sum(e_values_digi_e)'''

        #calibration_mip_to_reco = 0.00641222630095
        # Loop over digi hits and sum
        e_dig_E = 0
        if ecal_digCollection_b:
            #print("n barrel digi hits:", len(digCollection_b))
            #print(type(digCollection_b.data()))
            #tmp_E = sum(digCollection_b.data())
            for dig in ecal_digCollection_b: e_dig_E += dig.getEnergy()*e_calibration_mip_to_reco
        if ecal_digCollection_e:
            for dig in ecal_digCollection_e: e_dig_E += dig.getEnergy()*e_calibration_mip_to_reco

        # Loop over digi hits and sum
        h_dig_E = 0
        if hcal_digCollection_b:
            #print("n barrel digi hits:", len(digCollection_b))
            #print(type(digCollection_b.data()))
            #tmp_E = sum(digCollection_b.data())
            for dig in hcal_digCollection_b: h_dig_E += dig.getEnergy()*h_calibration_mip_to_reco
        if hcal_digCollection_e:
            for dig in hcal_digCollection_e: h_dig_E += dig.getEnergy()*h_calibration_mip_to_reco


        h_values_b = [x.getEnergy() for x in hcal_simCollection_b]
        h_values_e = [x.getEnergy() for x in hcal_simCollection_e]
        h_sim_E = sum(h_values_b) + sum(h_values_e)
        #print("e_Sim_E: {} | h_sim_E: {}".format(e_sim_E, h_sim_E))

        # Only make plots for events with isGood mcps
        if has_mcp_ob and h_sim_E > 0:
            # Calibrate the ecal energy and subtract off
            #ecal_calib_value = h_calib.GetBinContent(h_calib.FindBin(my_mcp_ob.Theta()))
            #estimated_hcal_energy = my_mcp_ob.E() - e_sim_E*e_sampling_scaling
            estimated_hcal_energy = my_mcp_ob.E() - e_dig_E
            #*ecal_calib_value
            #print(ecal_calib_value, my_mcp_ob.E(), estimated_hcal_energy)

            #print("filling mcp_ob_pt with", my_mcp_ob.Perp())
            hists2d["eta"].Fill(my_mcp_ob.Eta(), estimated_hcal_energy/h_sim_E)
            #hists2d["eta"].Fill(my_mcp_ob.Eta(), (estimated_hcal_energy - h_sim_E)/estimated_hcal_energy)
            #hists2d["eta"].Fill(my_mcp_ob.Eta(), (estimated_hcal_energy - h_sim_E)/estimated_hcal_energy)
            hists2d["theta"].Fill(my_mcp_ob.Theta(), estimated_hcal_energy/h_sim_E)
            #hists2d["E"].Fill(my_mcp_ob.E())

            E_truth[0] = my_mcp_ob.E()
            E_ECAL_barrel_sim[0] = sum(e_values_b)
            E_ECAL_endcap_sim[0] = sum(e_values_e)
            E_ECAL_sim_total[0] = e_sim_E

            E_HCAL_barrel_sim[0] = sum(h_values_b)
            E_HCAL_endcap_sim[0] = sum(h_values_e)
            E_HCAL_sim_total[0] = h_sim_E

            E_ECAL_digi_total[0] = e_dig_E/e_calibration_mip_to_reco
            E_HCAL_digi_total[0] = h_dig_E/h_calibration_mip_to_reco
            E_ECAL_mip_to_digi_total[0] = e_dig_E
            E_HCAL_mip_to_digi[0] = h_dig_E

            E_ECAL_total_corrected[0] = e_sim_E*e_sampling_scaling
            E_HCAL_total_corrected[0] = h_sim_E*h_sampling_scaling
            E_ECAL_HCAL_total_corrected[0] = E_ECAL_total_corrected[0] + E_HCAL_total_corrected[0]

        calibration_tree.Fill()
    reader.close()
    if(i % 10 == 0): print("Processed {} Files".format(i))

output_file = TFile(options.outFile + "_digi.root", 'RECREATE')
for histo in hists2d:
    hists2d[histo].Write()

    #c = TCanvas("cprof", "cprof")
    h_prof = hists2d[histo].ProfileX("_pfx", 1, -1, "s")
    h_prof.Write()

calibration_tree.Write()
        

#print("WEEEEE!")
