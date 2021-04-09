import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os

class photonRECOSFProducer(Module):
  def __init__( self , year ):
    self.year = year
    self.reco_input1 = "egamma_RECO_low.root"
    self.reco_input2 = "egamma_RECO_high.root"
    self.SF_location_path = "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/data/year%s/" %(os.environ['CMSSW_BASE'], self.year)
#    print 'SF location:', self.SF_location_path

  def beginJob(self):
    print 'begin to set Electron and Photon RECO SF --->>>'
    print 'start to open SF root file --->>>'
    self.reco_th2f_low = ROOT.TH2F()
    self.reco_th2f_high = ROOT.TH2F()
    self.file_reco_th2f_low = ROOT.TFile.Open(self.SF_location_path+self.reco_input1)
    self.file_reco_th2f_high = ROOT.TFile.Open(self.SF_location_path+self.reco_input2)
    self.file_reco_th2f_low.GetObject('EGamma_SF2D', self.reco_th2f_low)
    self.file_reco_th2f_high.GetObject('EGamma_SF2D', self.reco_th2f_high)
    print 'open SF files successfully --->>>'

  def endJob(self):
    print 'close SF root file --->>>'
    self.file_reco_th2f_low.Close()
    self.file_reco_th2f_high.Close()
    print 'finish setting Electron and Photon RECO and SF --->>>'
    
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch('Photon_RECO_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_RECO_SFerr','F', lenVar='nPhoton')
  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    
    photons = Collection(event, "Photon")
    if not (len(photons)>0): pass
    Photon_RECO_SF = []
    Photon_RECO_SFerr = []
    
    for ipho in range(0, len(photons)):
      if photons[ipho].pt < 20: 
        Photon_RECO_SF.append(self.reco_th2f_low.GetBinContent(self.reco_th2f_low.FindBin(photons[ipho].eta, photons[ipho].pt)))
        Photon_RECO_SFerr.append(self.reco_th2f_low.GetBinError(self.reco_th2f_low.FindBin(photons[ipho].eta, photons[ipho].pt)))
      if photons[ipho].pt >= 20 and photons[ipho].pt < 500:
        Photon_RECO_SF.append(self.reco_th2f_high.GetBinContent(self.reco_th2f_high.FindBin(photons[ipho].eta, photons[ipho].pt)))
        Photon_RECO_SFerr.append(self.reco_th2f_high.GetBinError(self.reco_th2f_high.FindBin(photons[ipho].eta, photons[ipho].pt)))
      if photons[ipho].pt >= 500:
        Photon_RECO_SF.append(self.reco_th2f_high.GetBinContent(self.reco_th2f_high.FindBin(photons[ipho].eta, 499)))
        Photon_RECO_SFerr.append(self.reco_th2f_high.GetBinError(self.reco_th2f_high.FindBin(photons[ipho].eta, 499)))
    self.out.fillBranch('Photon_RECO_SF', Photon_RECO_SF)
    self.out.fillBranch('Photon_RECO_SFerr', Photon_RECO_SFerr)

    return True

photonRECOSF2016 = lambda: photonRECOSFProducer("2016")
photonRECOSF2017 = lambda: photonRECOSFProducer("2017")
photonRECOSF2018 = lambda: photonRECOSFProducer("2018")
