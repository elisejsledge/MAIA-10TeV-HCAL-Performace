from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TFile, TLorentzVector, TMath, TTree
from math import *
from optparse import OptionParser
from array import array
import os
import fnmatch

#Adapted from Federico Meloni and Junjia Zhang

#########################
# parameters

parser = OptionParser()
parser.add_option('-i', '--inFileDir', help='--inFileDir /data/reco',
                  type=str, default='/data/reco')
parser.add_option('-p', '--inFilePrefix', help = '--inFilePrefix pionGun_pT_0_50', 
                  type=str, default='pionGun_pT_0_50')
parser.add_option('-o', '--outFile', help='--output file name', 
                  type=str, default = '/result/pion_jet_study_test')
parser.add_option('-v', action="store_true", help="-v if verbose", default=False)
(options, args) = parser.parse_args()


#########################
#loop over input files
print("Input file format: " + options.inFileDir + "/" + options.inFilePrefix + "/" + options.inFilePrefix + "_reco_" + "<num>")
for i in range(1000):
    file_label = i*10
    file_name = options.inFileDir + "/" + options.inFilePrefix + "/" + options.inFilePrefix + "_reco_" + str(file_label) + ".slcio"
    if(options.v): print(file_name)
    elif(i%10 == 0 and i>0): print("Completed " + str(file_label) + " files")

    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file_name)

    #loop over events in file
    for ievt, event in enumerate(reader):
        if ievt % 1 == 0:
            print("Processing event " + str(ievt))
        
        #look at the jets (which are Reconstructed Particle arrarys)
        jetCollection = event.getCollection('JetOut')
        print("Number of jets: {}".format(len(jetCollection)))
        for jet in jetCollection:
            dp3 = jet.getMomentum()
            tlv = TLorentzVector()
            tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

            print("pT: {} GeV || # of jet particles: {}".format(tlv.Perp(),len(jet.getParticles())))
            
            print("___________________________________\nJet Constituent Info:")
            for constituent in jet.getParticles():
                ids = constituent.getParticleIDs()
                print("PID: {} || Energy: {} GeV || Charge: {} || len(tracks): {}".format(ids, constituent.getEnergy(), constituent.getCharge(), len(constituent.getTracks())))


        
    break