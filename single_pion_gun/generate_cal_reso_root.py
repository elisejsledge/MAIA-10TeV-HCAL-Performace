from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TFile, TLorentzVector, TMath, TVector3
import math
from optparse import OptionParser
import os
import fnmatch
import sys

#ADAPTED FROM JUNJIA ZHANG
#########################
# parameters

parser = OptionParser()
parser.add_option('-i', '--inFileDir', help='--inFileDir /data/reco',
                  type=str, default='/data/reco')
parser.add_option('-p', '--inFilePrefix', help = '--inFilePrefix pionGun_pT_0_50', 
                  type=str, default='pionGun_pT_0_50')
parser.add_option('-o', '--outFile', help='--output file name', 
                  type=str, default = '/result/histos_gen_calo_0_50')
parser.add_option('-v', action="store_true", help="-v if verbose", default=False)
(options, args) = parser.parse_args()

Bfield = 5 #Tesla

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

#########################
# declare histograms

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

#individual clustered energy in calorimeter
ecal_clustered_energy = TH1D('ecal_clustered_energy', 'ecal_clustered_energy', pt_bins, pt_min, pt_max+2*pt_buffer) #GeV
hcal_clustered_energy = TH1D('hcal_clustered_energy', 'hcal_clustered_energy', pt_bins, pt_min, pt_max+2*pt_buffer) #GeV
ecal_endcap_clustered_energy = TH1D('ecal_endcap_clustered_energy', 'ecal_endcap_clustered_energy', pt_bins, pt_min, pt_max+2*pt_buffer) #GeV
hcal_endcap_clustered_energy = TH1D('hcal_endcap_clustered_energy', 'hcal_endcap_clustered_energy', pt_bins, pt_min, pt_max+2*pt_buffer) #GeV

#clustered calorimeter energy
clustered_cal_energy = TH1D('clustered_cal_energy', 'clustered_cal_energy', pt_bins, pt_min, pt_max+2*pt_buffer) #GeV
#resolution defined as (clustered_energy - gen_energy)/gen_energy
clustered_cal_resolution = TH1D('clustered_cal_resolution', 'clustered_cal_resolution', 100, -1, 1)


#########################
#loop over input files
print("Input file format: " + options.inFileDir + "/" + options.inFilePrefix + "/" + options.inFilePrefix + "_reco_" + "<num>")
for i in range(1000):
    file_label = i*10
    file_name = options.inFileDir + "/" + options.inFilePrefix + "/" + options.inFilePrefix + "_reco_" + str(file_label) + ".slcio"
    if(options.v): print(file_name)
    elif(i%10 == 0 and i>0): print("Completed " + str(file_label) + " files")
    
    # create a reader and open an LCIO file
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file_name)

    #loop over all events in file
    for ievt, event in enumerate(reader):
        if(options.v): print("New event, #" + str(ievt+file_label))

        pion_theta = 0
        pion_phi = 0
        has_pion = False
        generator_energy = 0

        #find the generator-level pion
        mcpCollection = event.getCollection('MCParticle')
        for part in mcpCollection:
            #the generator-level pion should have PDG number matching that of a pion and no daughters
            if (math.fabs(part.getPDG()) == 211 or part.getPDG() == 111) and len(part.getParents()) == 0:
                #get the generator-level pion endpiont
                endpoint = part.getEndpoint()
                endpoint_vec = TVector3()
                endpoint_vec.SetXYZ(endpoint[0], endpoint[1], endpoint[2])

                if endpoint_vec.Mag() > 0:
                    #sanity check, every event has a generator level pion with endpoint not at the origin
                    if (has_pion == True):
                        print("Warning: found more than one generator-level pions")
                    has_pion = True

                    #get the energy, pT, pz and charge of the generator level pion
                    generator_energy = part.getEnergy()
                    momentum = part.getMomentum()
                    pt = math.sqrt(momentum[0]*momentum[0] + momentum[1]*momentum[1])
                    pz = momentum[2]
                    pion_energy.Fill(generator_energy)
                    pion_pT.Fill(pt)
                    pion_pz.Fill(pz)
                    pion_charge.Fill(part.getCharge())
                    if(options.v): print(generator_energy)

                    #get the theta and phi of the momentum of generator level pion
                    momentum_vec = TVector3()
                    momentum_vec.SetXYZ(momentum[0], momentum[1], momentum[2])
                    pion_momentum_theta.Fill(momentum_vec.Theta())
                    pion_momentum_phi.Fill(momentum_vec.Phi())

                    #get the theta and phi of the generator level pion's endpoint
                    pion_theta = endpoint_vec.Theta()
                    pion_phi = endpoint_vec.Phi()
                    pion_end_theta.Fill(pion_theta)
                    pion_end_phi.Fill(pion_phi)

                    #the difference between pion momentum angles and pion endpoint angles
                    theta_diff = momentum_vec.Theta() - pion_theta
                    phi_diff = momentum_vec.Phi() - pion_phi
                    pion_theta_diff.Fill(theta_diff)
                    pion_phi_diff.Fill(phi_diff)
        
        if has_pion:
            clustered_energy = 0
            ecal_energy = 0
            hcal_energy = 0
            ecal_endcap_energy = 0
            hcal_endcap_energy = 0

            if collection_is_not_empty(event, "ECalBarrelCollection"):
                ecal_energy = get_clustered_energy(event, "EcalBarrelCollectionRec", pion_theta, pion_phi)
            if collection_is_not_empty(event, "HCalBarrelCollection"):
                hcal_energy = get_clustered_energy(event, "HcalBarrelsCollectionRec", pion_theta, pion_phi)
            if collection_is_not_empty(event, "ECalEndcapCollection"):
                ecal_endcap_energy = get_clustered_energy(event, "EcalEndcapCollectionRec", pion_theta, pion_phi)
            if collection_is_not_empty(event, "HCalEndcapCollection"):
                hcal_endcap_energy = get_clustered_energy(event, "HcalEndcapsCollectionRec", pion_theta, pion_phi)

            #fill in the individual calorimeter clustered energy histograms
            ecal_clustered_energy.Fill(ecal_energy)
            hcal_clustered_energy.Fill(hcal_energy)
            ecal_endcap_clustered_energy.Fill(ecal_endcap_energy)
            hcal_endcap_clustered_energy.Fill(hcal_endcap_energy)

            clustered_energy = ecal_energy + hcal_energy + ecal_endcap_energy + hcal_endcap_energy

            if(options.v):
                print("ECAL energy: {energy} GeV".format(energy=ecal_energy))
                print("HCAL energy: {energy} GeV".format(energy=hcal_energy))
                print("ECAP energy: {energy} GeV".format(energy=ecal_endcap_energy))
                print("HCAP energy: {energy} GeV".format(energy=hcal_endcap_energy))
                print("Cluster energy: {energy} GeV".format(energy=clustered_energy))
            
            #fill in the histogram
            clustered_cal_energy.Fill(clustered_energy)
            resolution = (clustered_energy - generator_energy)/generator_energy
            clustered_cal_resolution.Fill(resolution)
            
#########################
#write histograms
output_file_name = options.outFile + '.root'
output_file = TFile(output_file_name, 'RECREATE')

#energy and transverse momentum of generator level pion
pion_energy.Write()
pion_pT.Write()
pion_pz.Write()
pion_charge.Write()

#theta and phi of the generator level pion
pion_momentum_theta.Write()
pion_momentum_phi.Write()
pion_end_theta.Write()
pion_end_phi.Write()
pion_theta_diff.Write()
pion_phi_diff.Write()

#individual clustered energy in calorimeter
ecal_clustered_energy.Write()
hcal_clustered_energy.Write()
ecal_endcap_clustered_energy.Write()
hcal_endcap_clustered_energy.Write()

#clustered calorimeter energy
clustered_cal_energy.Write()
clustered_cal_resolution.Write()