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
from optparse import OptionParser
import itertools
from math import *
from array import array
import numpy as np

print("Imports succesful!")
ROOT.gROOT.SetBatch()

'''
Goal: find the digitization constant for the HCAL outside of the solonoid shadow
Forward region: 0.2 < theta < 0.5 & 2.6 < theta < 2.9
Barrel region: 1.25 < theta < 1.75

k = (Integral (HCal*Rec) - Integral (Ecal* Rec))/ Integral(HCal*Collection) # i think this is wrong
    for neutrons outside of the Solenoid shadow (i.e. ~ tetha <0.5)

k = (E_true - Integral (Ecal* Rec))/ Integral(HCal*Collection)

["MCParticle"] True particle collections

Look at the following collections
["ECalBarrelCollection"]: calculated neutron energies
["EcalBarrelCollectionDigi"]: 
["EcalBarrelCollectionRec"] 
    for both HCAL and ECAL
'''


#calibration constants
e_calibration_mip = 0.0001575
e_calibration_mip_to_reco = 0.00641222630095
e_sampling_scaling = e_calibration_mip_to_reco/e_calibration_mip

#hcal calibration constant guesses
h_calibration_mip = 0.0004825
h_calibration_mip_to_reco = 0.0231348530678
h_sampling_scaling = h_calibration_mip_to_reco/h_calibration_mip



# SOME IMPORTANT FUNCTIONS
def getTLV(obj):
    obj_p = obj.getMomentum()
    obj_e = obj.getEnergy()
    obj_tlv = TLorentzVector()
    #obj_tlv.SetPxPyPzE(obj_p.x, obj_p.y, obj_p.z, obj_e)
    obj_tlv.SetPxPyPzE(obj_p[0], obj_p[1], obj_p[2], obj_e)
    return obj_tlv

def IterativeGaussianFit(h, nSig = 1.5, fit_settings = "Q0"):
    #Adapted from Gillian Kopp: https://github.com/gk199/Run3-HCAL-LLP-Analysis/blob/main/MiniTuplePlotter/PlotFunctions.h#L103
    mean = h.GetMean()
    sigma = h.GetRMS()
    dm = 999
    ds = 999
    epsilon = 0.0000001 #gotta look if these vars are useful
    maxIter = 100
    iteration = 0

    fit_g = TF1("fit_g", "gaus")
    h.Fit(fit_g, "Q0", "", mean - (nSig * sigma), mean + (nSig * sigma))
    mean = fit_g.GetParameter(1)
    sigma = fit_g.GetParameter(2)

    while ((dm > epsilon) or (ds > epsilon)):
        h.Fit(fit_g, fit_settings, "", mean - (nSig * sigma), mean + (nSig * sigma))

        prevMean = mean
        prevSigma = sigma

        mean = fit_g.GetParameter(1)
        sigma = fit_g.GetParameter(2)

        dm = abs((mean - prevMean) / prevMean)
        ds = abs((sigma - prevSigma) / prevSigma)

        iteration += 1
        if(iteration > maxIter):
            break

    print("iter: {}, dm: {}, ds: {}, mean: {}, sigma: {}".format(iteration, dm, ds, mean, sigma))
    return fit_g



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



#create tree of for calorimeter values
calorimeter_tree = TTree("calorimeter_tree", "calorimeter_tree")
#Truth values
E_truth = array('d', [0]) #true neutron energy
theta_truth = array('d', [0]) #true neutron theta
#ECAL
ECAL_sim_total = array('d', [0]) #total sim energy in ecal 
ECAL_digi_total = array('d', [0]) #total digi energy in ecal 
ECAL_rec_total = array('d', [0]) #total rec energy in ecal 
#HCAL
HCAL_sim_total = array('d', [0]) #total sim energy in hcal 
HCAL_digi_total = array('d', [0]) #total digi energy in hcal 
HCAL_rec_total = array('d', [0]) #total rec energy in hcal 

#create tree branches
calorimeter_tree.Branch("E_truth",  E_truth,  'var/D')
calorimeter_tree.Branch("theta_truth",  theta_truth,  'var/D')
calorimeter_tree.Branch("ECAL_sim_total",  ECAL_sim_total,  'var/D')
calorimeter_tree.Branch("ECAL_digi_total",  ECAL_digi_total,  'var/D')
calorimeter_tree.Branch("ECAL_rec_total",  ECAL_rec_total,  'var/D')
calorimeter_tree.Branch("HCAL_sim_total",  HCAL_sim_total,  'var/D')
calorimeter_tree.Branch("HCAL_digi_total",  HCAL_digi_total,  'var/D')
calorimeter_tree.Branch("HCAL_rec_total",  HCAL_rec_total,  'var/D')


#create histos
'''
ideas:
estimated hcal energy
'''
# ############## CREATE EMPTY HISTOGRAM OBJECTS  #############################
# Set up histograms
# This is an algorithmic way of making a bunch of histograms and storing them in a dictionary
variables = {}
#variables["z"] = {"nbins":500, "xmin": -10000, "xmax": 10000, "title": "z [mm]"}
variables["scale"] = {"nbins":500, "xmin": 0, "xmax": 16000, "title": "true-over-sim ratio"}
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
    hists2d[var] = TH2F(var, var, variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"], variables["scale"]["nbins"], variables["scale"]["xmin"], variables["scale"]["xmax"])


#initialize 1D histos

h_hcal_cap_const = TH1D('hcal_cap_const', 'hcal_cap_const', 800, 0, 200)
h_hcal_barrel_const = TH1D('hcal_barrel_const', 'hcal_barrel_const', 800, 0, 200)

histos_list = [h_hcal_cap_const, h_hcal_barrel_const]


# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
# Loop over events
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "ECalBarrelCollection", "ECalEndcapCollection", "HCalBarrelCollection", "HCalEndcapCollection",
                                "EcalBarrelCollectionDigi", "EcalEndcapCollectionDigi", "HcalBarrelCollectionDigi", "HcalEndcapCollectionDigi",
                                "EcalBarrelCollectionRec", "EcalEndcapCollectionRec", "HcalBarrelCollectionRec", "HcalEndcapCollectionRec"])


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
    #if (i==51): break
    reader.open(file)
    for ievt, event in enumerate(reader):
        if options.v: print("Processing event " + str(ievt))

        #Get relevant collections
        mcpCollection = event.getCollection("MCParticle")

        ecal_barrel = event.getCollection("ECalBarrelCollection") #sim ecal barrel hits
        ecal_cap = event.getCollection("ECalEndcapCollection") #sim ecal endcap hits

        hcal_barrel = event.getCollection("HCalBarrelCollection") #sim hcal barrel hits
        hcal_cap = event.getCollection("HCalEndcapCollection") #sim hcal endcap hits

        #get DIGI and REC ecal collections
        try: ecal_barrel_digi = event.getCollection("EcalBarrelCollectionDigi")
        except: ecal_barrel_digi = None
        try: ecal_barrel_rec = event.getCollection("EcalBarrelCollectionRec")
        except: ecal_barrel_rec = None

        try: ecal_cap_digi = event.getCollection("EcalEndcapCollectionDigi")
        except: ecal_cap_digi = None
        try: ecal_cap_rec = event.getCollection("EcalEndcapCollectionRec")
        except: ecal_cap_rec = None

        #get DIGI and REC hcal collections
        try: hcal_barrel_digi = event.getCollection("HcalBarrelCollectionDigi")
        except: hcal_barrel_digi = None
        try: hcal_barrel_rec = event.getCollection("HcalBarrelCollectionRec")
        except: hcal_barrel_rec = None

        try: hcal_cap_digi = event.getCollection("HcalEndcapCollectionDigi")
        except: hcal_cap_digi = None
        try: hcal_cap_rec = event.getCollection("HcalEndcapCollectionRec")
        except: hcal_cap_rec = None

        neutron_tlv = getTLV(mcpCollection[0])
        if options.v: print("True Energy: {} | Theta: {} | PDG: {}".format(neutron_tlv.Energy(), neutron_tlv.Theta(), mcpCollection[0].getPDG()))
        #if (neutron_tlv.Theta() < 0.5 or neutron_tlv.Theta() > 2.6):
            #true_E = neutron_tlv.Energy()
            #print("len Barrel: {}, Digi: {}, Rec: {}".format(len(ecal_barrel), len(ecal_barrel_digi), len(ecal_barrel_rec)))
        
        #get all 3 different energy types for ECAL
        ecal_sim_E = 0
        #n_ecal_sim_hits = 0
        if ecal_barrel:
            if options.v: print("num ecal barrel sim hits: {}".format(len(ecal_barrel)))
            for sim in ecal_barrel: ecal_sim_E += sim.getEnergy()#multiply by sample scaling later
        if ecal_cap: 
            if options.v: print("num ecal endcap sim hits: {}".format(len(ecal_cap)))
            for sim in ecal_cap: ecal_sim_E += sim.getEnergy()#multiply by sample scaling later
        
        ecal_digi_E = 0
        if ecal_barrel_digi:
            if options.v: print("num ecal barrel digi hits: {}".format(len(ecal_barrel_digi)))
            for sim in ecal_barrel_digi: ecal_digi_E += sim.getEnergy()#multiply by mip to reco later
        if ecal_cap_digi: 
            if options.v: print("num ecal endcap digi hits: {}".format(len(ecal_cap_digi)))
            for sim in ecal_cap_digi: ecal_digi_E += sim.getEnergy()#multiply by mip to reco later
        
        ecal_rec_E = 0
        if ecal_barrel_rec:
            if options.v: print("num ecal barrel rec hits: {}".format(len(ecal_barrel_rec)))
            for sim in ecal_barrel_rec: ecal_rec_E += sim.getEnergy()#multiply by mip to reco later
        if ecal_cap_rec: 
            if options.v: print("num ecal endcap rec hits: {}".format(len(ecal_cap_rec)))
            for sim in ecal_cap_rec: ecal_rec_E += sim.getEnergy()#multiply by mip to reco later

        #get all 3 different energy types for HCAL
        hcal_sim_E = 0
        if hcal_barrel:
            if options.v: print("num hcal barrel sim hits: {}".format(len(hcal_barrel)))
            for sim in hcal_barrel: hcal_sim_E += sim.getEnergy()
        if hcal_cap: 
            if options.v: print("num hcal endcap sim hits: {}".format(len(hcal_cap)))
            for sim in hcal_cap: hcal_sim_E += sim.getEnergy()
        
        hcal_digi_E = 0
        if hcal_barrel_digi:
            if options.v: print("num hcal barrel digi hits: {}".format(len(hcal_barrel_digi)))
            for sim in hcal_barrel_digi: hcal_digi_E += sim.getEnergy()
        if hcal_cap_digi: 
            if options.v: print("num hcal endcap digi hits: {}".format(len(hcal_cap_digi)))
            for sim in hcal_cap_digi: hcal_digi_E += sim.getEnergy()
        
        hcal_rec_E = 0
        if hcal_barrel_rec:
            if options.v: print("num hcal barrel rec hits: {}".format(len(hcal_barrel_rec)))
            for sim in hcal_barrel_rec: hcal_rec_E += sim.getEnergy()
        if hcal_cap_rec: 
            if options.v: print("num hcal endcap rec hits: {}".format(len(hcal_cap_rec)))
            for sim in hcal_cap_rec: hcal_rec_E += sim.getEnergy()


        #fill tree

        E_truth[0] = neutron_tlv.Energy()
        theta_truth[0] = neutron_tlv.Theta()

        ECAL_sim_total[0] = ecal_sim_E
        ECAL_digi_total[0] = ecal_digi_E
        ECAL_rec_total[0] = ecal_rec_E

        HCAL_sim_total[0]= hcal_sim_E
        HCAL_digi_total[0] = hcal_digi_E
        HCAL_rec_total[0] = hcal_rec_E

        calorimeter_tree.Fill()

        #make 2d histos

        if hcal_sim_E > 0:
            k_calc_fede = (hcal_rec_E - ecal_rec_E)/hcal_sim_E
            k_calc_me = (neutron_tlv.Energy() - ecal_rec_E)/hcal_sim_E
            hists2d["eta"].Fill(neutron_tlv.Eta(), k_calc_me)
            hists2d["theta"].Fill(neutron_tlv.Theta(), k_calc_me)
            if 0.2 < neutron_tlv.Theta() < 0.5 or 2.6 < neutron_tlv.Theta() < 3.0: #in endcap region
                h_hcal_cap_const.Fill(k_calc_me)
            
            if 1.25 < neutron_tlv.Theta() < 1.75: #in barrel region
                h_hcal_barrel_const.Fill(k_calc_me)



        if options.v and hcal_sim_E > 0: 
            print("ECAL sim: {} | digi: {} | rec: {}".format(ecal_sim_E, ecal_digi_E, ecal_rec_E))
            print("ECAL CORRECTED sim: {} | digi: {} | rec: {}".format(ecal_sim_E*e_sampling_scaling, ecal_digi_E*e_calibration_mip_to_reco, ecal_rec_E))

            print("HCAL sim: {} | digi: {} | rec: {}".format(hcal_sim_E, hcal_digi_E, hcal_rec_E))
            print("HCAL CORRECTED guess sim: {} | digi: {} | rec: {}".format(hcal_sim_E*h_sampling_scaling, hcal_digi_E*h_calibration_mip_to_reco, hcal_rec_E))

            print("ECAL k:", e_sampling_scaling)
            print("HCAL k guess:", h_sampling_scaling)
            print("HCAL k fede calc: {} | HCAL k me calc: {}".format(k_calc_fede, k_calc_me))
            print("\n------------------------------------\n")
        '''
        next steps:
        1) make 2d digitization constant plots then profile them
        2) compare k guess and k calc by profiling the values
        3) make hists of all of the values that I have gotten out of this 
        4) Create plots with prettyness and good titles
        5) Compare my plots to 3TeV and 1.5 TeV to make sure I know its what I want/am looking for
        '''


    reader.close()
    if(i % 10 == 0): print("Processed {} Files".format(i))

output_file = TFile(options.outFile + ".root", 'RECREATE')
for histo in hists2d:
    hists2d[histo].Write()

    #c = TCanvas("cprof", "cprof")
    h_prof = hists2d[histo].ProfileX("_pfx", 1, -1, "s")
    h_prof.Write()

for histo in histos_list:
    nSig = 1.5
    fit_settings = "Q"
    gausIterFit = IterativeGaussianFit(histo, nSig, fit_settings)
    #gausIterFit.SetLineColor(kBlack) #this doesn't work to change the color

    #make sure to draw both
    histo.Fit("gaus", "+")

    histo.Write()

calorimeter_tree.Write()

'''
Note of to-dos:
'''