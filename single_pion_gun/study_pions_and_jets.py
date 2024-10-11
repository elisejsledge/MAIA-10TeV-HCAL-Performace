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

#########################
# define functions
def get_clustered_energy(event, collection_name, pion_theta, pion_phi):
    #calculate clustered energy in one of the calorimeters
    clustered_energy = 0
    hitscollection = event.getCollection(collection_name)

    for ihit, hit in enumerate(hitscollection):
        #get theta and phi of the hit position
        position = hit.getPosition()
        position_vec = TVector3()
        position_vec.SetXYZ(position[0], position[1], position[2])
        position_theta = position_vec.Theta()
        position_phi = position_vec.Phi()

        #check if the hit is in a cone of 0.1 in theta/phi around the generator level pion
        if (math.fabs(position_theta - pion_theta) < 0.1 and math.fabs(position_phi - pion_phi) < 0.1):
            clustered_energy += hit.getEnergy()

    return clustered_energy

def collection_is_not_empty(event, collection_name):
    # check if a calorimeter has a hit in it
    collection = event.getCollection(collection_name)
    count = 0
    for ihit, hit in enumerate(collection):
        count += 1
        break
    if count == 1: return True #check if at least 1 hit
    else: return False


#set histogram params
pt_min = 0
pt_max = 1000.
pt_bins = 250
pt_buffer = 250

#energy and transverse momentum of generator level pion
pion_energy = TH1D('pion_energy', 'pion_gen_energy', pt_bins, pt_min, pt_max+2*pt_buffer) #energy of generator pion
pion_pT = TH1D('pion_pT', 'pion_gen_pT', 2*pt_bins, pt_min, pt_max) #transverse momentum of generator pion
pion_pz = TH1D('pion_pz', 'pion_gen_pz', pt_bins, -pt_max, pt_max) #z-direction momentum of generator pion
pion_charge = TH1D('pion_charge', 'pion_gen_charge', 4, -2, 2) #charge of the pion (-1, 0, 1)

#theta and phi of the generator level pion
pion_momentum_theta = TH1D('pion_momentum_theta', 'pion_gen_momentum_theta', 100, 0, math.pi) #angle in radians
pion_momentum_phi = TH1D('pion_momentum_phi', 'pion_gen_momentum_phi', 100, -math.pi, math.pi) #angle in radians
pion_end_theta = TH1D('pion_end_theta', 'pion_gen_end_theta', 100, 0, math.pi) #angle in radians
pion_end_phi = TH1D('pion_end_phi', 'pion_gen_end_phi', 100, -math.pi, math.pi) #angle in radians
pion_theta_diff = TH1D('pion_theta_diff', 'pion_theta_diff', 100, -1, 1) #angle in radians
pion_phi_diff = TH1D('pion_phi_diff', 'pion_phi_diff', 100, -1, 1) #angle in radians

#jet histograms
jet_pT = TH1D('jet_pT','jet_pT', 2*pt_bins, pt_min, pt_max) #GeV
jet_E = TH1D('jet_E','jet_E', pt_bins, pt_min, pt_max+2*pt_buffer) #GeV
jet_num = TH1D('jet_num','jet_num', 15, 0, 15) #number of jets per initial pion
jet_const_num = TH1D('jet_const_num','jet_const_num', 30, 0, 30) #number of constituent particles per jet
jet_const_pT = TH1D('jet_const_pT','jet_const_pT', 2*pt_bins, pt_min, pt_max) #GeV
jet_const_E = TH1D('jet_const_E','jet_const_E', pt_bins, pt_min, pt_max+2*pt_buffer) #GeV
const_type = TH1D('const_type', "const_type", 3000, 0, 3000) #PDG of jet constituent

#PFO histograms
Npfo = TH1D('Npfo', "Npfo", 1000, 0, 1000) #number of PFOs
pfo_type = TH1D('pfo_type', "pfo_type", 3000, 0, 3000) #PDG of jet constituent

#jet reco histograms
pion_jet_reco_pt = TH1D('pion_jet_reco_pt','pion_jet_reco_pt', 2*pt_bins, pt_min, pt_max) #GeV
pion_jet_reco_eff = TH1D('pion_jet_reco_eff', 'pion_jet_reco_eff', 2*pt_bins, pt_min, pt_max) #error
pion_jet_reco_resolution = TH1D('pion_jet_reco_resolution', 'pion_jet_reco_resolution', 100, -1, 1)

#pfo reco histograms
pion_pfo_reco_pt = TH1D('pion_pfo_reco_pt','pion_pfo_reco_pt', 2*pt_bins, pt_min, pt_max) #GeV
pion_pfo_reco_eff = TH1D('pion_pfo_reco_eff', 'pion_pfo_reco_eff', 2*pt_bins, pt_min, pt_max) #error
pion_pfo_reco_resolution = TH1D('pion_pfo_reco_resolution', 'pion_pfo_reco_resolution', 100, -1, 1)

#clustered calorimeter energy
clustered_cal_energy = TH1D('clustered_cal_energy', 'clustered_cal_energy', pt_bins, pt_min, pt_max+2*pt_buffer) #GeV
#resolution defined as (clustered_energy - gen_energy)/gen_energy
clustered_cal_resolution = TH1D('clustered_cal_resolution', 'clustered_cal_resolution', 100, -1, 1)

#2D Histograms
pion_jet_pt = TH2D('pion_jet_pt','pion_jet_pt', 2*pt_bins, pt_min, pt_max, 4*pt_bins, 0, pt_max)
pion_pt_const_pt = TH2D('pion_pt_const_pt','pion_pt_const_pt', 2*pt_bins, pt_min, pt_max, 4*pt_bins, 0, pt_max)
pion_pt_jet_num = TH2D('pion_pt_jet_num','pion_pt_jet_num', 2*pt_bins, pt_min, pt_max, 15, 0, 15)
pion_pt_const_num = TH2D('pion_pt_const_num','pion_pt_const_num', 2*pt_bins, pt_min, pt_max, 30, 0, 30)


histos_list = [pion_energy, pion_pT, pion_pz, pion_charge,
               pion_momentum_theta, pion_momentum_phi,
               pion_end_theta, pion_end_phi, pion_theta_diff, pion_phi_diff,
               jet_pT, jet_E, jet_num,
               jet_const_num, jet_const_pT, jet_const_E, const_type,
               pfo_type, pion_pfo_reco_pt, pion_pfo_reco_eff, pion_pfo_reco_resolution,
               pion_jet_reco_pt, pion_jet_reco_eff, pion_jet_reco_resolution, 
               clustered_cal_energy, clustered_cal_resolution,
               pion_jet_pt, pion_pt_const_pt, pion_pt_jet_num, pion_pt_const_num,
            ]


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

        mcpCollection = event.getCollection('MCParticle')
        initial_energy = mcpCollection[0].getEnergy() #keep initial energy to compare to jets later
        pion_energy.Fill(mcpCollection[0].getEnergy()) #first particle is the pion
        pion_charge.Fill(mcpCollection[0].getCharge())

        #find pion endpoint coordinates
        endpoint = mcpCollection[0].getEndpoint()
        endpoint_vec = TVector3()
        endpoint_vec.SetXYZ(endpoint[0], endpoint[1], endpoint[2])
        pion_end_theta.Fill(endpoint_vec.Theta())
        pion_end_phi.Fill(endpoint_vec.Phi())

        #find pion initial diretion
        dp3 = mcpCollection[0].getMomentum()
        pion_pz.Fill(dp3[2])
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcpCollection[0].getEnergy())
        pion_trans_mom = tlv.Perp()
        pion_pT.Fill(pion_trans_mom)
        pion_momentum_theta.Fill(tlv.Theta())
        pion_momentum_phi.Fill(tlv.Phi())

        #find angular separation of pion initial and final coords
        pion_theta_diff.Fill(tlv.Theta( )- endpoint_vec.Theta())
        pion_phi_diff.Fill(tlv.Phi() - endpoint_vec.Phi())

        #for part in mcpCollection:
            #print("PDG: {}".format(part.getPDG()))

        #iterate through pre-reconstructed jet objets
        jetCollection = event.getCollection('JetOut')
        jet_num.Fill(len(jetCollection))
        pion_pt_jet_num.Fill(pion_trans_mom, len(jetCollection))

        ################################
        #Junjia Analysis
        clustered_energy = 0
        ecal_energy = 0
        hcal_energy = 0
        ecal_endcap_energy = 0
        hcal_endcap_energy = 0

        pion_theta = tlv.Theta()
        pion_phi = tlv.Phi()

        if collection_is_not_empty(event, "ECalBarrelCollection"):
            ecal_energy = get_clustered_energy(event, "EcalBarrelCollectionRec", pion_theta, pion_phi)
        if collection_is_not_empty(event, "HCalBarrelCollection"):
            hcal_energy = get_clustered_energy(event, "HcalBarrelsCollectionRec", pion_theta, pion_phi)
        if collection_is_not_empty(event, "ECalEndcapCollection"):
            ecal_endcap_energy = get_clustered_energy(event, "EcalEndcapCollectionRec", pion_theta, pion_phi)
        if collection_is_not_empty(event, "HCalEndcapCollection"):
            hcal_endcap_energy = get_clustered_energy(event, "HcalEndcapsCollectionRec", pion_theta, pion_phi)

        #fill in the individual calorimeter clustered energy histograms
        #ecal_clustered_energy.Fill(ecal_energy)
        #hcal_clustered_energy.Fill(hcal_energy)
        #ecal_endcap_clustered_energy.Fill(ecal_endcap_energy)
        #hcal_endcap_clustered_energy.Fill(hcal_endcap_energy)

        clustered_energy = ecal_energy + hcal_energy + ecal_endcap_energy + hcal_endcap_energy

        if(options.v):
            print("ECAL energy: {energy} GeV".format(energy=ecal_energy))
            print("HCAL energy: {energy} GeV".format(energy=hcal_energy))
            print("ECAP energy: {energy} GeV".format(energy=ecal_endcap_energy))
            print("HCAP energy: {energy} GeV".format(energy=hcal_endcap_energy))
            print("Cluster energy: {energy} GeV".format(energy=clustered_energy))
        
        #fill in the histogram
        clustered_cal_energy.Fill(clustered_energy)
        resolution = (clustered_energy - initial_energy)/initial_energy
        clustered_cal_resolution.Fill(resolution)

        #################################################
        #Jet Analysis
        jet_energies = []
        for jet in jetCollection:
            #get jet momentum and energy
            dp3 = jet.getMomentum()
            tlv = TLorentzVector()
            tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

            jet_E.Fill(tlv.Energy())
            jet_pT.Fill(tlv.Perp())
            jet_energies.append(tlv.Energy())

            pion_jet_pt.Fill(pion_trans_mom, tlv.Perp())

            #iterate through jet constituent particles
            jet_const_num.Fill(len(jet.getParticles()))
            pion_pt_const_num.Fill(pion_trans_mom,len(jet.getParticles()))

            for constituent in jet.getParticles():
                dp3 = constituent.getMomentum()
                tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

                jet_const_E.Fill(tlv.Energy())
                jet_const_pT.Fill(tlv.Perp())
                const_type.Fill(abs(constituent.getType()))
                pion_pt_const_pt.Fill(pion_trans_mom, tlv.Perp())

                if(abs(constituent.getType()) == 211): #check if constituent is a pion
                    pion_jet_reco_pt.Fill(tlv.Perp()) #get reconstructed energy
                    pion_jet_reco_resolution.Fill((tlv.Energy() - initial_energy)/initial_energy) #get energy reconsonstruction resolution

                if(options.v): print("Constituent Mass: {} | Energy: {} | Type: {} ".format(constituent.getMass(), constituent.getEnergy(), constituent.getType()))
        
        pfoCollection = event.getCollection("PandoraPFOs")
        Npfo.Fill(len(pfoCollection))


        #################################
        #PFO Analysis
        #go through PFOs
        for pfo in pfoCollection:
            pfo_type.Fill(abs(pfo.getType()))
            tlv = TLorentzVector()
            tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getEnergy())

            if(abs(pfo.getType()) == 211): #check if PFO is a pion
                pion_pfo_reco_pt.Fill(tlv.Perp())
                pion_pfo_reco_resolution.Fill((tlv.Energy() - initial_energy)/initial_energy)

    if(i % 100 == 0): print("Processed {} Files".format(i))
    
    #close reader
    reader.close()

pion_jet_reco_eff.Divide(pion_jet_reco_pt,pion_pT)
pion_pfo_reco_eff.Divide(pion_pfo_reco_pt, pion_pT)

#save histograms
output_file = TFile(options.outFile + ".root", 'RECREATE')
for histo in histos_list:
    histo.Write()


output_file.Close()