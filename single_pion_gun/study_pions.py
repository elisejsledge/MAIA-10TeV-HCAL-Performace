from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TMath, TTree, TVector3
import math
from optparse import OptionParser
from array import array
import os
import fnmatch


#Adapted from Federico Meloni and Junjia Zhang
#########################
# parameters
parser = OptionParser()
parser.add_option('-i', '--inFileDir', help='--inFileDir /data/reco/',
                  type=str, default='/data/reco/')
parser.add_option('-p', '--inFilePrefix', help = '--inFilePrefix pionGun_pT_0_50', 
                  type=str, default='pionGun_pT_0_50')
parser.add_option('-o', '--outFile', help='--output file name', 
                  type=str, default = 'pion_jet_study_test')
parser.add_option('-v', action="store_true", help="-v if verbose", default=False)

(options, args) = parser.parse_args()

#Set Histo bins
arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 50., 100., 250., 500., 1000., 2500., 5000.))


# declare histograms
h_truth_E = TH1D('truth_E', 'truth_E', len(arrBins_E)-1, arrBins_E) #GeV: true generator energy
h_truth_pT = TH1D('truth_pT', 'truth_pT', len(arrBins_E)-1, arrBins_E) #GeV: true generator pT
h_truth_theta = TH1D('truth_theta', 'truth_theta', len(arrBins_theta)-1, arrBins_theta) #Rads: true generator theta
h_matched_theta = TH1D('matched_theta', 'matched_theta', len(arrBins_theta)-1, arrBins_theta)
h_matched_E = TH1D('matched_E', 'matched_E', len(arrBins_E)-1, arrBins_E)

h_Npfo = TH1D('Npfo', "Npfo", 1000, 0, 1000) #Number of PFO objets
h_EMpfo_E = TH1D('EMpfo_E', 'EMpfo_E', len(arrBins_E)-1, arrBins_E)
h_HADpfo_E = TH1D('HADpfo_E', 'HADpfo_E', len(arrBins_E)-1, arrBins_E)
h_matchedEMpfo_E = TH1D('matchedEMpfo_E', 'matchedEMpfo_E', len(arrBins_E)-1, arrBins_E)
h_matchedEMpfo_theta = TH1D('matchedEMpfo_theta', 'matchedEMpfo_theta',
                            len(arrBins_theta)-1, arrBins_theta)
h_matchedHADpfo_E = TH1D('matchedHADpfo_E', 'matchedHADpfo_E', len(arrBins_E)-1, arrBins_E)
h_matchedHADpfo_theta = TH1D('matchedHADpfo_theta', 'matchedHADpfo_theta',
                             len(arrBins_theta)-1, arrBins_theta)
h_pfo_type = TH1D('pfo_type', "pfo_type", 3000, 0, 3000) #pfo object PDG number
h_jet_const_type = TH1D('jet_const_type', "jet_const_type", 3000, 0, 3000) #pfo object PDG number

h_deltaEM_E = TH1D('deltaEM_E', 'deltaEM_E', 250, -1000, 1000)
h_deltaHAD_E = TH1D('deltaHAD_E', 'deltaHAD_E', 250, -1000, 1000)
h_delta_E_sumE = TH1D('delta_E_sumE', 'delta_E_sumE', 250, -1000, 1000)

'''
# Low-level sim hit distributions
h_ECAL_simhit_E = TH1D('ECAL_simhit_E', 'ECAL_simhit_E', 100, 0, 20)  # GeV
h_ECAL_simhit_layer = TH1D(
    'ECAL_simhit_layer', 'ECAL_simhit_layer', 100, 0, 100)
h_ECAL_simhit_layer_ele = TH1D(
    'ECAL_simhit_layer_ele', 'ECAL_simhit_layer_ele', 100, 0, 100)
h_ECAL_simhit_layer_gamma = TH1D(
    'ECAL_simhit_layer_gamma', 'ECAL_simhit_layer_gamma', 100, 0, 100)
h_ECAL_simhit_layer_other = TH1D(
    'ECAL_simhit_layer_other', 'ECAL_simhit_layer_other', 100, 0, 100)

h_HCAL_simhit_E = TH1D('HCAL_simhit_E', 'HCAL_simhit_E', 100, 0, 20)  # GeV
h_HCAL_simhit_layer = TH1D(
    'HCAL_simhit_layer', 'HCAL_simhit_layer', 100, 0, 100)
'''

'''
# Low-level digitised hit distributions
h_ECAL_hit_time = TH1D('ECAL_hit_time', 'ECAL_hit_time', 100, -10, 10)  # ns
h_ECAL_hit_E = TH1D('ECAL_hit_E', 'ECAL_hit_E', 100, 0, 20)  # GeV
h_ECAL_hit_R = TH1D('ECAL_hit_R', 'ECAL_hit_R', 100, 1700, 4000)  # m
h_ECAL_hit_layer = TH1D('ECAL_hit_layer', 'ECAL_hit_layer', 100, 0, 100)

h_HCAL_hit_time = TH1D('HCAL_hit_time', 'HCAL_hit_time', 100, -10, 10)  # ns
h_HCAL_hit_E = TH1D('HCAL_hit_E', 'HCAL_hit_E', 100, 0, 20)  # GeV
h_HCAL_hit_R = TH1D('HCAL_hit_R', 'HCAL_hit_R', 100, 1700, 4000)  # m
h_HCAL_hit_layer = TH1D('HCAL_hit_layer', 'HCAL_hit_layer', 100, 0, 100)

# Aggregated energy info
h_sumE = TH1D('sumE', 'sumE', 120, 0, 6000)  # GeV
h_ECAL_sumE = TH1D('ECAL_sumE', 'ECAL_sumE', 120, 0, 6000)  # GeV
h_HCAL_sumE = TH1D('HCAL_sumE', 'HCAL_sumE', 120, 0, 6000)  # GeV
'''
h_EMfrac = TH1D('EMfrac', 'EMfrac', 100, 0, 1)  # GeV
h_EMfrac_PFO = TH1D('EMfrac_PFO', 'EMfrac_PFO', 100, 0, 1)  # GeV


# Histo list for writing to outputs
histos_list = [h_truth_E, h_truth_theta, h_matched_theta,
               h_EMpfo_E, h_matched_E,
               h_HADpfo_E,
               h_matchedEMpfo_E, h_matchedEMpfo_theta,
               h_matchedHADpfo_E, h_matchedHADpfo_theta,
               h_deltaEM_E, h_deltaHAD_E, #h_delta_E_sumE,
               h_Npfo, h_pfo_type,
               h_jet_const_type, h_truth_pT,
               #h_ECAL_hit_time, h_ECAL_hit_E, h_ECAL_hit_R,
               #h_HCAL_hit_time, h_HCAL_hit_E, h_HCAL_hit_R,
               #h_ECAL_simhit_E, h_HCAL_simhit_E,
               #h_ECAL_sumE, h_HCAL_sumE,
               h_EMfrac, h_EMfrac_PFO,
               #h_ECAL_hit_layer, h_HCAL_hit_layer,
               #h_ECAL_simhit_layer, h_ECAL_simhit_layer_ele, h_ECAL_simhit_layer_gamma, h_ECAL_simhit_layer_other,
               #h_HCAL_simhit_layer
               ]

for histo in histos_list:
    histo.SetDirectory(0)

####################################
pion_tree = TTree("pion_tree", "pion_tree")
E = array('d', [0]) #reconstructed pion energy
pT = array('d', [0]) #reconstructed pion pT
phi = array('d', [0]) #reconstructed pion phi
theta = array('d', [0]) #reconstructed pion theta
E_truth = array('d', [0]) #true pion energy
pT_truth = array('d', [0]) #true pion pT
phi_truth = array('d', [0]) #true pion phi
theta_truth = array('d', [0]) #true pion theta
pion_tree.Branch("E",  E,  'var/D')
pion_tree.Branch("pT",  pT,  'var/D')
pion_tree.Branch("phi", phi, 'var/D')
pion_tree.Branch("theta", theta, 'var/D')
pion_tree.Branch("E_truth",  E_truth,  'var/D')
pion_tree.Branch("pT_truth",  pT_truth,  'var/D')
pion_tree.Branch("phi_truth", phi_truth, 'var/D')
pion_tree.Branch("theta_truth", theta_truth, 'var/D')


#create an array of the slcio files to process
to_process = []
file_dir = options.inFileDir + options.inFilePrefix
if os.path.isdir(file_dir):
    for file_name in os.listdir(file_dir): 
        to_process.append(file_dir + "/" + file_name)

i = 0 #count number of files processed
for file in to_process:
    i+=1
    #create reader and open LCIO file
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    for ievt, event in enumerate(reader):
        if options.v:
            print("Processing event " + str(ievt))
        
        #Fill the truth level histos 
        mcpCollection = event.getCollection('MCParticle')
        h_truth_E.Fill(mcpCollection[0].getEnergy()) #true energy (first particle is generator)
        dp3 = mcpCollection[0].getMomentum()
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcpCollection[0].getEnergy())
        h_truth_theta.Fill(tlv.Theta())
        h_truth_pT.Fill(tlv.Perp())

        E_truth[0] = mcpCollection[0].getEnergy()
        pT_truth[0] = tlv.Perp()
        phi_truth[0] = tlv.Phi()
        theta_truth[0] = tlv.Theta()

        # Fill the reco-level histos
        pfoCollection = event.getCollection('PandoraPFOs')
        h_Npfo.Fill(len(pfoCollection))

        # Match true pfo with closest reco PFO in deltaR
        matchedEM_E = -1.
        matchedEM_theta = -1.
        matchedHAD_E = -1.
        matchedHAD_theta = -1.
        allEM_E = 0.
        allHAD_E = 0.

        minDREM = 999999.
        minDRHAD = 999999.

        for pfo in pfoCollection:
            h_pfo_type.Fill(abs(pfo.getType())) #fill PFO object PDG number
            dp3 = pfo.getMomentum()
            tlv_pfo = TLorentzVector()
            tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getEnergy())

            if abs(pfo.getType()) == 22: #pfo is a photon
                allEM_E = allEM_E + pfo.getEnergy() #add photon energy to total EM energy
            elif abs(pfo.getType()) == 2112: #pfo is a neutron
                allHAD_E = allHAD_E + pfo.getEnergy()
            
            dR = tlv_pfo.DeltaR(tlv) #DeltaR from pfo to inital photon

            #if dR < minDREM and abs(pfo.getType()) == 211: #closest reconstructed pion #I FEEL LIKE THIS SHOULD BE TYPE 22
            if dR < minDREM and abs(pfo.getType()) == 22:
                minDREM = dR
                matchedEM_E = pfo.getEnergy()
                matchedEM_theta = tlv_pfo.Theta()
            
            if dR < minDRHAD and abs(pfo.getType()) == 2112: #closest reconstructed neutron
                minDRHAD = dR
                matchedHAD_E = pfo.getEnergy()
                matchedHAD_theta = tlv_pfo.Theta()
        
        if(options.v):
            print("Total EM_E: {}, Total HAD_E: {}, Matched EM_E: {}, Matched HAD_E: {}".format(allEM_E, allHAD_E, matchedEM_E, matchedHAD_E))

        h_EMpfo_E.Fill(allEM_E)
        h_HADpfo_E.Fill(allHAD_E)

        if matchedEM_E > 0:
            h_matchedEMpfo_E.Fill(matchedEM_E)
            h_matchedEMpfo_theta.Fill(matchedEM_theta)
            h_deltaEM_E.Fill(matchedEM_E-mcpCollection[0].getEnergy())
        if matchedHAD_E > 0:
            h_matchedHADpfo_E.Fill(matchedHAD_E)
            h_matchedHADpfo_theta.Fill(matchedHAD_theta)
            h_deltaHAD_E.Fill(matchedHAD_E-mcpCollection[0].getEnergy())

        if allHAD_E+allEM_E > 0:
            h_EMfrac_PFO.Fill(allEM_E/(allHAD_E+allEM_E))


        #look at the anti-kt jet
        jetCollection = event.getCollection('JetOut')
        minDRjet = 999999.
        E[0] = -1
        pT[0] = -1
        phi[0] = -4
        theta[0] = -1

        for jet in jetCollection:
            dp3 = jet.getMomentum()
            tlv_pfo = TLorentzVector()
            tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

            dR = tlv_pfo.DeltaR(tlv)

            if dR < minDRjet:
                minDRjet = dR
                E[0] =jet.getEnergy()
                pT[0] = tlv.Perp()
                phi[0] = tlv_pfo.Phi()
                theta[0] = tlv_pfo.Theta()
            
            for constituent in jet.getParticles():
                h_jet_const_type.Fill(abs(constituent.getType()))
            
        
        if minDRjet<0.2:
            h_matched_theta.Fill(tlv.Theta())
            h_matched_E.Fill(tlv.Energy())


        pion_tree.Fill()
    
    if(i % 100 == 0): print("Processed {} Files".format(i))

    #close reader
    reader.close()

# write histograms
output_file = TFile(options.outFile + "_ntup_pfoPFO.root", 'RECREATE')
for histo in histos_list:
    histo.Write()
pion_tree.Write()
output_file.Close()
