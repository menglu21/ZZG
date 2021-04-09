import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os

class photonIDSFProducer(Module):
  def __init__( self , year ):
    self.year = year
    self.id_loose = "photon_Loose.root"
    self.id_medium = "photon_Medium.root"
    self.id_tight = "photon_Tight.root"
    self.id_mva80 = "photon_MVA80.root"
    self.id_mva90 = "photon_MVA90.root"
    self.eleveto_CSEV = "eleveto_CSEV.root"
    self.eleveto_PV = "eleveto_PV.root"
    self.SF_location_path = "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/data/year%s/" %(os.environ['CMSSW_BASE'], self.year)
    print 'SF location:', self.SF_location_path

  def beginJob(self):
    print 'begin to set Photon ID SF --->>>'
    print 'start to open SF root file --->>>'
    # init the TH*F
    self.id_loose_th2f= ROOT.TH2F()
    self.id_medium_th2f= ROOT.TH2F()
    self.id_tight_th2f= ROOT.TH2F()
    self.id_mva80_th2f= ROOT.TH2F()
    self.id_mva90_th2f= ROOT.TH2F()
    self.dirLoose_CSEV = ROOT.TDirectoryFile()
    self.dirMedium_CSEV = ROOT.TDirectoryFile()
    self.dirTight_CSEV = ROOT.TDirectoryFile()
    self.dirMVA_CSEV = ROOT.TDirectoryFile()
    self.eleveto_loose_th1f_CSEV = ROOT.TH1F()
    self.eleveto_medium_th1f_CSEV = ROOT.TH1F()
    self.eleveto_tight_th1f_CSEV = ROOT.TH1F()
    self.eleveto_mva_th1f_CSEV = ROOT.TH1F()
    self.dirLoose_PV = ROOT.TDirectoryFile()
    self.dirMedium_PV = ROOT.TDirectoryFile()
    self.dirTight_PV = ROOT.TDirectoryFile()
    self.dirMVA_PV = ROOT.TDirectoryFile()
    self.eleveto_loose_th1f_PV = ROOT.TH1F()
    self.eleveto_medium_th1f_PV = ROOT.TH1F()
    self.eleveto_tight_th1f_PV = ROOT.TH1F()
    self.eleveto_mva_th1f_PV = ROOT.TH1F()
    #Open the SF root file
    self.file_id_loose= ROOT.TFile.Open(self.SF_location_path+self.id_loose)
    self.file_id_medium= ROOT.TFile.Open(self.SF_location_path+self.id_medium)
    self.file_id_tight= ROOT.TFile.Open(self.SF_location_path+self.id_tight)
    self.file_id_mva80= ROOT.TFile.Open(self.SF_location_path+self.id_mva80)
    self.file_id_mva90= ROOT.TFile.Open(self.SF_location_path+self.id_mva90)
    self.file_eleveto_CSEV= ROOT.TFile.Open(self.SF_location_path+self.eleveto_CSEV)
    self.file_eleveto_PV= ROOT.TFile.Open(self.SF_location_path+self.eleveto_PV)
    #access to the TH*F
    self.file_id_loose.GetObject('EGamma_SF2D', self.id_loose_th2f)
    self.file_id_medium.GetObject('EGamma_SF2D', self.id_medium_th2f)
    self.file_id_tight.GetObject('EGamma_SF2D', self.id_tight_th2f)
    self.file_id_mva80.GetObject('EGamma_SF2D', self.id_mva80_th2f)
    self.file_id_mva90.GetObject('EGamma_SF2D', self.id_mva90_th2f)

    self.file_eleveto_CSEV.GetObject('LooseID', self.dirLoose_CSEV)
    self.file_eleveto_CSEV.GetObject('MediumID', self.dirMedium_CSEV)
    self.file_eleveto_CSEV.GetObject('TightID', self.dirTight_CSEV)
    self.file_eleveto_CSEV.GetObject('MVAID', self.dirMVA_CSEV)
    self.dirLoose_CSEV.GetObject('SF_CSEV_LooseID', self.eleveto_loose_th1f_CSEV)
    self.dirMedium_CSEV.GetObject('SF_CSEV_MediumID', self.eleveto_medium_th1f_CSEV)
    self.dirTight_CSEV.GetObject('SF_CSEV_TightID', self.eleveto_tight_th1f_CSEV)
    self.dirMVA_CSEV.GetObject('SF_CSEV_MVAID', self.eleveto_mva_th1f_CSEV)

    self.file_eleveto_PV.GetObject('LooseID', self.dirLoose_PV)
    self.file_eleveto_PV.GetObject('MediumID', self.dirMedium_PV)
    self.file_eleveto_PV.GetObject('TightID', self.dirTight_PV)
    self.file_eleveto_PV.GetObject('MVAID', self.dirMVA_PV)
    self.dirLoose_PV.GetObject('SF_HasPix_LooseID', self.eleveto_loose_th1f_PV)
    self.dirMedium_PV.GetObject('SF_HasPix_MediumID', self.eleveto_medium_th1f_PV)
    self.dirTight_PV.GetObject('SF_HasPix_TightID', self.eleveto_tight_th1f_PV)
    self.dirMVA_PV.GetObject('SF_HasPix_MVAID', self.eleveto_mva_th1f_PV)
    print 'open SF files successfully --->>>'

  def endJob(self):
    print 'close SF root file --->>>'
    self.file_id_loose.Close()
    self.file_id_medium.Close()
    self.file_id_tight.Close()
    self.file_id_mva80.Close()
    self.file_id_mva90.Close()
    self.file_eleveto_CSEV.Close()
    self.file_eleveto_PV.Close()
    print 'finish setting Photon ID SF --->>>'
    
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch('Photon_CutBased_LooseID_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_CutBased_LooseID_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_CutBased_MediumID_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_CutBased_MediumID_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_CutBased_TightID_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_CutBased_TightID_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_MVAFall17V2_WP80_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_MVAFall17V2_WP80_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_MVAFall17V2_WP90_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_MVAFall17V2_WP90_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_CSEV_LooseID_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_CSEV_LooseID_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_CSEV_MediumID_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_CSEV_MediumID_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_CSEV_TightID_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_CSEV_TightID_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_CSEV_MVAID_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_CSEV_MVAID_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_PV_LooseID_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_PV_LooseID_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_PV_MediumID_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_PV_MediumID_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_PV_TightID_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_PV_TightID_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_PV_MVAID_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_eleVeto_PV_MVAID_SFerr','F', lenVar='nPhoton')
  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    
    photons = Collection(event, "Photon")
    if not (len(photons)>0): pass
    Photon_CutBased_LooseID_SF = []
    Photon_CutBased_LooseID_SFerr = []
    Photon_CutBased_MediumID_SF = []
    Photon_CutBased_MediumID_SFerr = []
    Photon_CutBased_TightID_SF = []
    Photon_CutBased_TightID_SFerr = []
    Photon_MVAFall17V2_WP80_SF = []
    Photon_MVAFall17V2_WP80_SFerr = []
    Photon_MVAFall17V2_WP90_SF = []
    Photon_MVAFall17V2_WP90_SFerr = []
    Photon_eleVeto_CSEV_LooseID_SF = []
    Photon_eleVeto_CSEV_LooseID_SFerr = []
    Photon_eleVeto_CSEV_MediumID_SF = []
    Photon_eleVeto_CSEV_MediumID_SFerr = []
    Photon_eleVeto_CSEV_TightID_SF = []
    Photon_eleVeto_CSEV_TightID_SFerr = []
    Photon_eleVeto_CSEV_MVAID_SF = []
    Photon_eleVeto_CSEV_MVAID_SFerr = []
    Photon_eleVeto_PV_LooseID_SF = []
    Photon_eleVeto_PV_LooseID_SFerr = []
    Photon_eleVeto_PV_MediumID_SF = []
    Photon_eleVeto_PV_MediumID_SFerr = []
    Photon_eleVeto_PV_TightID_SF = []
    Photon_eleVeto_PV_TightID_SFerr = []
    Photon_eleVeto_PV_MVAID_SF = []
    Photon_eleVeto_PV_MVAID_SFerr = []
    
    for ipho in range(0, len(photons)):
      if photons[ipho].isScEtaEB:
        if photons[ipho].r9 > 0.96:
          Photon_eleVeto_CSEV_LooseID_SF.append(self.eleveto_loose_th1f_CSEV.GetBinContent(2))
          Photon_eleVeto_CSEV_LooseID_SFerr.append(self.eleveto_loose_th1f_CSEV.GetBinError(2))
          Photon_eleVeto_CSEV_MediumID_SF.append(self.eleveto_medium_th1f_CSEV.GetBinContent(2))
          Photon_eleVeto_CSEV_MediumID_SFerr.append(self.eleveto_medium_th1f_CSEV.GetBinError(2))
          Photon_eleVeto_CSEV_TightID_SF.append(self.eleveto_tight_th1f_CSEV.GetBinContent(2))
          Photon_eleVeto_CSEV_TightID_SFerr.append(self.eleveto_tight_th1f_CSEV.GetBinError(2))
          Photon_eleVeto_CSEV_MVAID_SF.append(self.eleveto_mva_th1f_CSEV.GetBinContent(2))
          Photon_eleVeto_CSEV_MVAID_SFerr.append(self.eleveto_mva_th1f_CSEV.GetBinError(2))
          Photon_eleVeto_PV_LooseID_SF.append(self.eleveto_loose_th1f_PV.GetBinContent(2))
          Photon_eleVeto_PV_LooseID_SFerr.append(self.eleveto_loose_th1f_PV.GetBinError(2))
          Photon_eleVeto_PV_MediumID_SF.append(self.eleveto_medium_th1f_PV.GetBinContent(2))
          Photon_eleVeto_PV_MediumID_SFerr.append(self.eleveto_medium_th1f_PV.GetBinError(2))
          Photon_eleVeto_PV_TightID_SF.append(self.eleveto_tight_th1f_PV.GetBinContent(2))
          Photon_eleVeto_PV_TightID_SFerr.append(self.eleveto_tight_th1f_PV.GetBinError(2))
          Photon_eleVeto_PV_MVAID_SF.append(self.eleveto_mva_th1f_PV.GetBinContent(2))
          Photon_eleVeto_PV_MVAID_SFerr.append(self.eleveto_mva_th1f_PV.GetBinError(2))
        else:
          Photon_eleVeto_CSEV_LooseID_SF.append(self.eleveto_loose_th1f_CSEV.GetBinContent(3))
          Photon_eleVeto_CSEV_LooseID_SFerr.append(self.eleveto_loose_th1f_CSEV.GetBinError(3))
          Photon_eleVeto_CSEV_MediumID_SF.append(self.eleveto_medium_th1f_CSEV.GetBinContent(3))
          Photon_eleVeto_CSEV_MediumID_SFerr.append(self.eleveto_medium_th1f_CSEV.GetBinError(3))
          Photon_eleVeto_CSEV_TightID_SF.append(self.eleveto_tight_th1f_CSEV.GetBinContent(3))
          Photon_eleVeto_CSEV_TightID_SFerr.append(self.eleveto_tight_th1f_CSEV.GetBinError(3))
          Photon_eleVeto_CSEV_MVAID_SF.append(self.eleveto_mva_th1f_CSEV.GetBinContent(3))
          Photon_eleVeto_CSEV_MVAID_SFerr.append(self.eleveto_mva_th1f_CSEV.GetBinError(3))
          Photon_eleVeto_PV_LooseID_SF.append(self.eleveto_loose_th1f_PV.GetBinContent(3))
          Photon_eleVeto_PV_LooseID_SFerr.append(self.eleveto_loose_th1f_PV.GetBinError(3))
          Photon_eleVeto_PV_MediumID_SF.append(self.eleveto_medium_th1f_PV.GetBinContent(3))
          Photon_eleVeto_PV_MediumID_SFerr.append(self.eleveto_medium_th1f_PV.GetBinError(3))
          Photon_eleVeto_PV_TightID_SF.append(self.eleveto_tight_th1f_PV.GetBinContent(3))
          Photon_eleVeto_PV_TightID_SFerr.append(self.eleveto_tight_th1f_PV.GetBinError(3))
          Photon_eleVeto_PV_MVAID_SF.append(self.eleveto_mva_th1f_PV.GetBinContent(3))
          Photon_eleVeto_PV_MVAID_SFerr.append(self.eleveto_mva_th1f_PV.GetBinError(3))
      if photons[ipho].isScEtaEE:
        if photons[ipho].r9 > 0.96:
          Photon_eleVeto_CSEV_LooseID_SF.append(self.eleveto_loose_th1f_CSEV.GetBinContent(5))
          Photon_eleVeto_CSEV_LooseID_SFerr.append(self.eleveto_loose_th1f_CSEV.GetBinError(5))
          Photon_eleVeto_CSEV_MediumID_SF.append(self.eleveto_medium_th1f_CSEV.GetBinContent(5))
          Photon_eleVeto_CSEV_MediumID_SFerr.append(self.eleveto_medium_th1f_CSEV.GetBinError(5))
          Photon_eleVeto_CSEV_TightID_SF.append(self.eleveto_tight_th1f_CSEV.GetBinContent(5))
          Photon_eleVeto_CSEV_TightID_SFerr.append(self.eleveto_tight_th1f_CSEV.GetBinError(5))
          Photon_eleVeto_CSEV_MVAID_SF.append(self.eleveto_mva_th1f_CSEV.GetBinContent(5))
          Photon_eleVeto_CSEV_MVAID_SFerr.append(self.eleveto_mva_th1f_CSEV.GetBinError(5))
          Photon_eleVeto_PV_LooseID_SF.append(self.eleveto_loose_th1f_PV.GetBinContent(5))
          Photon_eleVeto_PV_LooseID_SFerr.append(self.eleveto_loose_th1f_PV.GetBinError(5))
          Photon_eleVeto_PV_MediumID_SF.append(self.eleveto_medium_th1f_PV.GetBinContent(5))
          Photon_eleVeto_PV_MediumID_SFerr.append(self.eleveto_medium_th1f_PV.GetBinError(5))
          Photon_eleVeto_PV_TightID_SF.append(self.eleveto_tight_th1f_PV.GetBinContent(5))
          Photon_eleVeto_PV_TightID_SFerr.append(self.eleveto_tight_th1f_PV.GetBinError(5))
          Photon_eleVeto_PV_MVAID_SF.append(self.eleveto_mva_th1f_PV.GetBinContent(5))
          Photon_eleVeto_PV_MVAID_SFerr.append(self.eleveto_mva_th1f_PV.GetBinError(5))
        else:
          Photon_eleVeto_CSEV_LooseID_SF.append(self.eleveto_loose_th1f_CSEV.GetBinContent(6))
          Photon_eleVeto_CSEV_LooseID_SFerr.append(self.eleveto_loose_th1f_CSEV.GetBinError(6))
          Photon_eleVeto_CSEV_MediumID_SF.append(self.eleveto_medium_th1f_CSEV.GetBinContent(6))
          Photon_eleVeto_CSEV_MediumID_SFerr.append(self.eleveto_medium_th1f_CSEV.GetBinError(6))
          Photon_eleVeto_CSEV_TightID_SF.append(self.eleveto_tight_th1f_CSEV.GetBinContent(6))
          Photon_eleVeto_CSEV_TightID_SFerr.append(self.eleveto_tight_th1f_CSEV.GetBinError(6))
          Photon_eleVeto_CSEV_MVAID_SF.append(self.eleveto_mva_th1f_CSEV.GetBinContent(6))
          Photon_eleVeto_CSEV_MVAID_SFerr.append(self.eleveto_mva_th1f_CSEV.GetBinError(6))
          Photon_eleVeto_PV_LooseID_SF.append(self.eleveto_loose_th1f_PV.GetBinContent(6))
          Photon_eleVeto_PV_LooseID_SFerr.append(self.eleveto_loose_th1f_PV.GetBinError(6))
          Photon_eleVeto_PV_MediumID_SF.append(self.eleveto_medium_th1f_PV.GetBinContent(6))
          Photon_eleVeto_PV_MediumID_SFerr.append(self.eleveto_medium_th1f_PV.GetBinError(6))
          Photon_eleVeto_PV_TightID_SF.append(self.eleveto_tight_th1f_PV.GetBinContent(6))
          Photon_eleVeto_PV_TightID_SFerr.append(self.eleveto_tight_th1f_PV.GetBinError(6))
          Photon_eleVeto_PV_MVAID_SF.append(self.eleveto_mva_th1f_PV.GetBinContent(6))
          Photon_eleVeto_PV_MVAID_SFerr.append(self.eleveto_mva_th1f_PV.GetBinError(6))

      if photons[ipho].pt < 500: 
        Photon_CutBased_LooseID_SF.append(self.id_loose_th2f.GetBinContent(self.id_loose_th2f.FindBin(photons[ipho].eta, photons[ipho].pt)))
        Photon_CutBased_LooseID_SFerr.append(self.id_loose_th2f.GetBinError(self.id_loose_th2f.FindBin(photons[ipho].eta, photons[ipho].pt)))
        Photon_CutBased_MediumID_SF.append(self.id_medium_th2f.GetBinContent(self.id_medium_th2f.FindBin(photons[ipho].eta, photons[ipho].pt)))
        Photon_CutBased_MediumID_SFerr.append(self.id_medium_th2f.GetBinError(self.id_medium_th2f.FindBin(photons[ipho].eta, photons[ipho].pt)))
        Photon_CutBased_TightID_SF.append(self.id_tight_th2f.GetBinContent(self.id_tight_th2f.FindBin(photons[ipho].eta, photons[ipho].pt)))
        Photon_CutBased_TightID_SFerr.append(self.id_tight_th2f.GetBinError(self.id_tight_th2f.FindBin(photons[ipho].eta, photons[ipho].pt)))
        Photon_MVAFall17V2_WP80_SF.append(self.id_mva80_th2f.GetBinContent(self.id_mva80_th2f.FindBin(photons[ipho].eta, photons[ipho].pt)))
        Photon_MVAFall17V2_WP80_SFerr.append(self.id_mva80_th2f.GetBinError(self.id_mva80_th2f.FindBin(photons[ipho].eta, photons[ipho].pt)))
        Photon_MVAFall17V2_WP90_SF.append(self.id_mva90_th2f.GetBinContent(self.id_mva90_th2f.FindBin(photons[ipho].eta, photons[ipho].pt)))
        Photon_MVAFall17V2_WP90_SFerr.append(self.id_mva90_th2f.GetBinError(self.id_mva90_th2f.FindBin(photons[ipho].eta, photons[ipho].pt)))
      else: 
        Photon_CutBased_LooseID_SF.append(self.id_loose_th2f.GetBinContent(self.id_loose_th2f.FindBin(photons[ipho].eta, 499)))
        Photon_CutBased_LooseID_SFerr.append(self.id_loose_th2f.GetBinError(self.id_loose_th2f.FindBin(photons[ipho].eta, 499)))
        Photon_CutBased_MediumID_SF.append(self.id_medium_th2f.GetBinContent(self.id_medium_th2f.FindBin(photons[ipho].eta, 499)))
        Photon_CutBased_MediumID_SFerr.append(self.id_medium_th2f.GetBinError(self.id_medium_th2f.FindBin(photons[ipho].eta, 499)))
        Photon_CutBased_TightID_SF.append(self.id_tight_th2f.GetBinContent(self.id_tight_th2f.FindBin(photons[ipho].eta, 499)))
        Photon_CutBased_TightID_SFerr.append(self.id_tight_th2f.GetBinError(self.id_tight_th2f.FindBin(photons[ipho].eta, 499)))
        Photon_MVAFall17V2_WP80_SF.append(self.id_mva80_th2f.GetBinContent(self.id_mva80_th2f.FindBin(photons[ipho].eta, 499)))
        Photon_MVAFall17V2_WP80_SFerr.append(self.id_mva80_th2f.GetBinError(self.id_mva80_th2f.FindBin(photons[ipho].eta, 499)))
        Photon_MVAFall17V2_WP90_SF.append(self.id_mva90_th2f.GetBinContent(self.id_mva90_th2f.FindBin(photons[ipho].eta, 499)))
        Photon_MVAFall17V2_WP90_SFerr.append(self.id_mva90_th2f.GetBinError(self.id_mva90_th2f.FindBin(photons[ipho].eta, 499)))


    self.out.fillBranch('Photon_CutBased_LooseID_SF', Photon_CutBased_LooseID_SF)
    self.out.fillBranch('Photon_CutBased_LooseID_SFerr', Photon_CutBased_LooseID_SFerr)
    self.out.fillBranch('Photon_CutBased_MediumID_SF', Photon_CutBased_MediumID_SF)
    self.out.fillBranch('Photon_CutBased_MediumID_SFerr', Photon_CutBased_MediumID_SFerr)
    self.out.fillBranch('Photon_CutBased_TightID_SF', Photon_CutBased_TightID_SF)
    self.out.fillBranch('Photon_CutBased_TightID_SFerr', Photon_CutBased_TightID_SFerr)
    self.out.fillBranch('Photon_MVAFall17V2_WP80_SF', Photon_MVAFall17V2_WP80_SF)
    self.out.fillBranch('Photon_MVAFall17V2_WP80_SFerr', Photon_MVAFall17V2_WP80_SFerr)
    self.out.fillBranch('Photon_MVAFall17V2_WP90_SF', Photon_MVAFall17V2_WP90_SF)
    self.out.fillBranch('Photon_MVAFall17V2_WP90_SFerr', Photon_MVAFall17V2_WP90_SFerr)
    self.out.fillBranch('Photon_eleVeto_CSEV_LooseID_SF',Photon_eleVeto_CSEV_LooseID_SF)
    self.out.fillBranch('Photon_eleVeto_CSEV_LooseID_SFerr',Photon_eleVeto_CSEV_LooseID_SFerr)
    self.out.fillBranch('Photon_eleVeto_CSEV_MediumID_SF',Photon_eleVeto_CSEV_MediumID_SF)
    self.out.fillBranch('Photon_eleVeto_CSEV_MediumID_SFerr',Photon_eleVeto_CSEV_MediumID_SFerr)
    self.out.fillBranch('Photon_eleVeto_CSEV_TightID_SF',Photon_eleVeto_CSEV_TightID_SF)
    self.out.fillBranch('Photon_eleVeto_CSEV_TightID_SFerr',Photon_eleVeto_CSEV_TightID_SFerr)
    self.out.fillBranch('Photon_eleVeto_CSEV_MVAID_SF',Photon_eleVeto_CSEV_MVAID_SF)
    self.out.fillBranch('Photon_eleVeto_CSEV_MVAID_SFerr',Photon_eleVeto_CSEV_MVAID_SFerr)
    self.out.fillBranch('Photon_eleVeto_PV_LooseID_SF',Photon_eleVeto_PV_LooseID_SF)
    self.out.fillBranch('Photon_eleVeto_PV_LooseID_SFerr',Photon_eleVeto_PV_LooseID_SFerr)
    self.out.fillBranch('Photon_eleVeto_PV_MediumID_SF',Photon_eleVeto_PV_MediumID_SF)
    self.out.fillBranch('Photon_eleVeto_PV_MediumID_SFerr',Photon_eleVeto_PV_MediumID_SFerr)
    self.out.fillBranch('Photon_eleVeto_PV_TightID_SF',Photon_eleVeto_PV_TightID_SF)
    self.out.fillBranch('Photon_eleVeto_PV_TightID_SFerr',Photon_eleVeto_PV_TightID_SFerr)
    self.out.fillBranch('Photon_eleVeto_PV_MVAID_SF',Photon_eleVeto_PV_MVAID_SF)
    self.out.fillBranch('Photon_eleVeto_PV_MVAID_SFerr',Photon_eleVeto_PV_MVAID_SFerr)

    return True

photonIDSF2016 = lambda: photonIDSFProducer("2016")
photonIDSF2017 = lambda: photonIDSFProducer("2017")
photonIDSF2018 = lambda: photonIDSFProducer("2018")
