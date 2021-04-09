import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os

class ZZGProducer(Module):
  def __init__( self , year ):
    self.year = year
  def beginJob(self):
    pass
  def endJob(self):
    pass
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch("channel", "I")
    self.out.branch("Z1_l1_pt", "F")
    self.out.branch("Z1_l1_eta", "F")
    self.out.branch("Z1_l1_phi", "F")
    self.out.branch("Z1_l1_mass", "F")
    self.out.branch("Z1_l1_pdgId", "I")
    self.out.branch("Z1_l2_pt", "F")
    self.out.branch("Z1_l2_eta", "F")
    self.out.branch("Z1_l2_phi", "F")
    self.out.branch("Z1_l2_mass", "F")
    self.out.branch("Z1_l2_pdgId", "I")
    self.out.branch("Z2_l1_pt", "F")
    self.out.branch("Z2_l1_eta", "F")
    self.out.branch("Z2_l1_phi", "F")
    self.out.branch("Z2_l1_mass", "F")
    self.out.branch("Z2_l1_pdgId", "I")
    self.out.branch("Z2_l2_pt", "F")
    self.out.branch("Z2_l2_eta", "F")
    self.out.branch("Z2_l2_phi", "F")
    self.out.branch("Z2_l2_mass", "F")
    self.out.branch("Z2_l2_pdgId", "I")
    self.out.branch("Z1_pt", "F")
    self.out.branch("Z1_eta", "F")
    self.out.branch("Z1_phi", "F")
    self.out.branch("Z1_mass", "F")
    self.out.branch("Z2_pt", "F")
    self.out.branch("Z2_eta", "F")
    self.out.branch("Z2_phi", "F")
    self.out.branch("Z2_mass", "F")
    self.out.branch("photon_pt", "F")
    self.out.branch("photon_eta", "F")
    self.out.branch("photon_phi", "F")
    self.out.branch("photon_mass", "F")
    self.out.branch("photon_isScEtaEB","B")
    self.out.branch("photon_isScEtaEE","B")
    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))
  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    
    # PV selection
    if (event.PV_npvsGood<1): return False

    # trigger selection

    # total number of ele+muon
    if ((event.nMuon + event.nElectron) < 4): return False

    # Muon selection: tight cut-based ID + tight PF iso
    muons = Collection(event, 'Muon')
    muon_v4_temp=TLorentzVector()
    tightMuons = []
    tightMuons_pdgid = []
    for imu in range(0, event.nMuon):
      if (muons[imu].tightId and muons[imu].pfRelIso04_all<0.15 and abs(muons[imu].eta)<2.4 and event.Muon_corrected_pt[imu]>5):
        muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
        tightMuons.append(muon_v4_temp.Clone())
        tightMuons_pdgid.append(muons[imu].pdgId)
    
    # electron selection: tight cut-based ID + impact parameter cut
    electrons = Collection(event, 'Electron')
    electron_v4_temp=TLorentzVector()
    tightElectrons = []
    tightElectrons_pdgid = []
    for iele in range(0, event.nElectron):
      if (electrons[iele].cutBased==4 and electrons[iele].tightCharge==2 and ((abs(electrons[iele].eta+electrons[iele].deltaEtaSC) <1.4442 and abs(electrons[iele].dxy)<0.05 and abs(electrons[iele].dz<0.1)) or (abs(electrons[iele].eta + electrons[iele].deltaEtaSC)>1.566 and abs(electrons[iele].eta + electrons[iele].deltaEtaSC)<2.4 and abs(electrons[iele].dxy)<0.1 and abs(electrons[iele].dz)<0.2)) and electrons[iele].pt>7):
        electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
        tightElectrons.append(electron_v4_temp.Clone())
        tightElectrons_pdgid.append(electrons[iele].pdgId)

    # require at elase 4 tight leptons
    if not (len(tightElectrons) + len(tightMuons) ==4):
      return False

    tightLeptons = tightMuons + tightElectrons
    tightLeptons.sort(key=lambda x: x.Pt(), reverse=True)
    # require 0mu-4ele or 2mu-2ele or 4mu-0ele
    if (len(tightElectrons)==1 or len(tightElectrons)==3 or len(tightMuons)==1 or len(tightMuons)==3):
      return False

    # 4 ele case: at least one ele with pt>20, at least 2 eles with pt>12
    if len(tightMuons)==0:
      self.out.fillBranch("channel",1)
      if not (tightElectrons[0].Pt()>20 and tightElectrons[1].Pt()>12):
        return False
      if not (tightElectrons_pdgid[0]+tightElectrons_pdgid[1]+tightElectrons_pdgid[2]+tightElectrons_pdgid[3] ==0):
        return False

    # 2mu-2ele case: at least one lep with pt>20, at least 2 leps with pt>12
    if len(tightMuons)==2:
      self.out.fillBranch("channel",2)
      if not (tightLeptons[0].Pt()>20 and tightLeptons[1].Pt()>12):
        return False
      if not (tightMuons_pdgid[0]+tightMuons_pdgid[1]==0):
        return False
      if not (tightElectrons_pdgid[0]+tightElectrons_pdgid[1]==0):
        return False

    # 4 muon case: at least one muon with pt>20, at least 2 muon with pt>10
    if len(tightMuons)==4:
      self.out.fillBranch("channel",3)
      if not (tightMuons[0].Pt()>20 and tightMuons[1].Pt()>10):
        return False
      if not (tightMuons_pdgid[0]+tightMuons_pdgid[1]+tightMuons_pdgid[2]+tightMuons_pdgid[3]==0):
        return False

    # require DeltaR(l1,l2)>0.02, following ZZto4L
    if not (tightLeptons[0].DeltaR(tightLeptons[1])>0.02 and tightLeptons[0].DeltaR(tightLeptons[2])>0.02 and tightLeptons[0].DeltaR(tightLeptons[3])>0.02 and tightLeptons[1].DeltaR(tightLeptons[2])>0.02 and tightLeptons[1].DeltaR(tightLeptons[3])>0.02 and tightLeptons[2].DeltaR(tightLeptons[3])>0.02):
      return False

    # require mll>4 regardless the flavor and charge
    if not ((tightLeptons[0]+tightLeptons[1]).M()>4 and (tightLeptons[0]+tightLeptons[2]).M()>4 and (tightLeptons[0]+tightLeptons[3]).M()>4 and (tightLeptons[1]+tightLeptons[2]).M()>4 and (tightLeptons[1]+tightLeptons[3]).M()>4 and (tightLeptons[2]+tightLeptons[3]).M()>4):
      return False

    # require DeltaR(ele, mu)>0.05 to remove spurious ghost leptons formed from ambiguities in track reconstruction, following ZZto4L
    if len(tightMuons)==2:
      if not (tightMuons[0].DeltaR(tightElectrons[0])>0.05 and tightMuons[0].DeltaR(tightElectrons[1])>0.05 and tightMuons[1].DeltaR(tightElectrons[0])>0.05 and tightMuons[1].DeltaR(tightElectrons[1])>0.05):
        return False

    Z1=TLorentzVector()
    Z2=TLorentzVector()
    Z1_lepton1=TLorentzVector()
    Z1_lepton2=TLorentzVector()
    Z2_lepton1=TLorentzVector()
    Z2_lepton2=TLorentzVector()

    # construct two Z boson, the leading Z is the one closer to Z mass
    # both Z mass need between (4, 120) GeV

    if len(tightMuons)==2:
      if not ((tightMuons[0]+tightMuons[1]).M()>4 and (tightMuons[0]+tightMuons[1]).M()<120 and (tightElectrons[0]+tightElectrons[1]).M()>4 and (tightElectrons[0]+tightElectrons[1]).M()<120):
        return False
      if abs((tightMuons[0]+tightMuons[1]).M()-91.1876) < abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876):
        Z1=tightMuons[0]+tightMuons[1]
        Z2=tightElectrons[0]+tightElectrons[1]
        Z1_lepton1 = tightMuons[0]
        Z1_lepton2 = tightMuons[1]
        Z2_lepton1 = tightElectrons[0]
        Z2_lepton2 = tightElectrons[1]
        self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[0])
        self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[1])
        self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
        self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[1])
      else:
        Z1=tightElectrons[0]+tightElectrons[1]
        Z2=tightMuons[0]+tightMuons[1]
        Z1_lepton1 = tightElectrons[0]
        Z1_lepton2 = tightElectrons[1]
        Z2_lepton1 = tightMuons[0]
        Z2_lepton2 = tightMuons[1]
        self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[0])
        self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[1])
        self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
        self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[1])
        
    pass_Z1_basic_1stpair = False
    pass_Z2_basic_1stpair = False
    pass_Z1_basic_2ndpair = False
    pass_Z2_basic_2ndpair = False

    # 4 muons case
    if len(tightElectrons)==0:
      if tightMuons_pdgid[0]+tightMuons_pdgid[1] == 0:
        if tightMuons_pdgid[0]+tightMuons_pdgid[2] == 0:
          pass_Z1_basic_1stpair = (tightMuons[0]+tightMuons[1]).M()>4 and (tightMuons[0]+tightMuons[1]).M()<120
          pass_Z2_basic_1stpair = (tightMuons[2]+tightMuons[3]).M()>4 and (tightMuons[2]+tightMuons[3]).M()<120
          pass_Z1_basic_2ndpair = (tightMuons[0]+tightMuons[2]).M()>4 and (tightMuons[0]+tightMuons[2]).M()<120
          pass_Z2_basic_2ndpair = (tightMuons[1]+tightMuons[3]).M()>4 and (tightMuons[1]+tightMuons[3]).M()<120
          if pass_Z1_basic_1stpair and pass_Z2_basic_1stpair:
            if pass_Z1_basic_2ndpair and pass_Z2_basic_2ndpair:
              if min(abs((tightMuons[0]+tightMuons[1]).M()-91.1876), abs((tightMuons[2]+tightMuons[3]).M()-91.1876)) < min(abs((tightMuons[0]+tightMuons[2]).M()-91.1876), abs((tightMuons[1]+tightMuons[3]).M()-91.1876)):
                if abs((tightMuons[0]+tightMuons[1]).M()-91.1876) < abs((tightMuons[2]+tightMuons[3]).M()-91.1876):
                  Z1 = tightMuons[0]+tightMuons[1]
                  Z2 = tightMuons[2]+tightMuons[3]
                  Z1_lepton1 = tightMuons[0]
                  Z1_lepton2 = tightMuons[1]
                  Z2_lepton1 = tightMuons[2]
                  Z2_lepton2 = tightMuons[3]
                  self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[0])
                  self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[1])
                  self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[2])
                  self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[3])
                else:
                  Z1 = tightMuons[2]+tightMuons[3]
                  Z2 = tightMuons[0]+tightMuons[1]
                  Z1_lepton1 = tightMuons[2]
                  Z1_lepton2 = tightMuons[3]
                  Z2_lepton1 = tightMuons[0]
                  Z2_lepton2 = tightMuons[1]
                  self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[2])
                  self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[3])
                  self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
                  self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[1])
              else:

                if abs((tightMuons[0]+tightMuons[2]).M()-91.1876) < abs((tightMuons[1]+tightMuons[3]).M()-91.1876):
                  Z1 = tightMuons[0]+tightMuons[2]
                  Z2 = tightMuons[1]+tightMuons[3]
                  Z1_lepton1 = tightMuons[0]
                  Z1_lepton2 = tightMuons[2]
                  Z2_lepton1 = tightMuons[1]
                  Z2_lepton2 = tightMuons[3]
                  self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[0])
                  self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[2])
                  self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[1])
                  self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[3])
                else:
                  Z1 = tightMuons[1]+tightMuons[3]
                  Z2 = tightMuons[0]+tightMuons[2]
                  Z1_lepton1 = tightMuons[1]
                  Z1_lepton2 = tightMuons[3]
                  Z2_lepton1 = tightMuons[0]
                  Z2_lepton2 = tightMuons[2]
                  self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[1])
                  self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[3])
                  self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
                  self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[2])

            else:
              if abs((tightMuons[0]+tightMuons[1]).M()-91.1876) < abs((tightMuons[2]+tightMuons[3]).M()-91.1876):
                Z1 = tightMuons[0]+tightMuons[1]
                Z2 = tightMuons[2]+tightMuons[3]
                Z1_lepton1 = tightMuons[0]
                Z1_lepton2 = tightMuons[1]
                Z2_lepton1 = tightMuons[2]
                Z2_lepton2 = tightMuons[3]
                self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[0])
                self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[1])
                self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[2])
                self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[3])
              else:
                Z1 = tightMuons[2]+tightMuons[3]
                Z2 = tightMuons[0]+tightMuons[1]
                Z1_lepton1 = tightMuons[2]
                Z1_lepton2 = tightMuons[3]
                Z2_lepton1 = tightMuons[0]
                Z2_lepton2 = tightMuons[1]
                self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[2])
                self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[3])
                self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[1])
          else:
            if pass_Z1_basic_2ndpair and pass_Z2_basic_2ndpair:
              if abs((tightMuons[0]+tightMuons[2]).M()-91.1876) < abs((tightMuons[1]+tightMuons[3]).M()-91.1876):
                Z1 = tightMuons[0]+tightMuons[2]
                Z2 = tightMuons[1]+tightMuons[3]
                Z1_lepton1 = tightMuons[0]
                Z1_lepton2 = tightMuons[2]
                Z2_lepton1 = tightMuons[1]
                Z2_lepton2 = tightMuons[3]
                self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[0])
                self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[2])
                self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[1])
                self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[3])
              else:
                Z1 = tightMuons[1]+tightMuons[3]
                Z2 = tightMuons[0]+tightMuons[2]
                Z1_lepton1 = tightMuons[1]
                Z1_lepton2 = tightMuons[3]
                Z2_lepton1 = tightMuons[0]
                Z2_lepton2 = tightMuons[2]
                self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[1])
                self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[3])
                self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[2])

        else:
          pass_Z1_basic_1stpair = (tightMuons[0]+tightMuons[1]).M()>4 and (tightMuons[0]+tightMuons[1]).M()<120
          pass_Z2_basic_1stpair = (tightMuons[2]+tightMuons[3]).M()>4 and (tightMuons[2]+tightMuons[3]).M()<120
          pass_Z1_basic_2ndpair = (tightMuons[0]+tightMuons[3]).M()>4 and (tightMuons[0]+tightMuons[3]).M()<120
          pass_Z2_basic_2ndpair = (tightMuons[1]+tightMuons[2]).M()>4 and (tightMuons[1]+tightMuons[2]).M()<120

          if pass_Z1_basic_1stpair and pass_Z2_basic_1stpair:
            if pass_Z1_basic_2ndpair and pass_Z2_basic_2ndpair:
              if min(abs((tightMuons[0]+tightMuons[1]).M()-91.1876), abs((tightMuons[2]+tightMuons[3]).M()-91.1876)) < min(abs((tightMuons[0]+tightMuons[3]).M()-91.1876), abs((tightMuons[1]+tightMuons[2]).M()-91.1876)):
                if abs((tightMuons[0]+tightMuons[1]).M()-91.1876) < abs((tightMuons[2]+tightMuons[3]).M()-91.1876):
                  Z1 = tightMuons[0]+tightMuons[1]
                  Z2 = tightMuons[2]+tightMuons[3]
                  Z1_lepton1 = tightMuons[0]
                  Z1_lepton2 = tightMuons[1]
                  Z2_lepton1 = tightMuons[2]
                  Z2_lepton2 = tightMuons[3]
                  self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[0])
                  self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[1])
                  self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[2])
                  self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[3])
                else:
                  Z1 = tightMuons[2]+tightMuons[3]
                  Z2 = tightMuons[0]+tightMuons[1]
                  Z1_lepton1 = tightMuons[2]
                  Z1_lepton2 = tightMuons[3]
                  Z2_lepton1 = tightMuons[0]
                  Z2_lepton2 = tightMuons[1]
                  self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[2])
                  self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[3])
                  self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
                  self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[1])
              else:

                if abs((tightMuons[0]+tightMuons[3]).M()-91.1876) < abs((tightMuons[1]+tightMuons[2]).M()-91.1876):
                  Z1 = tightMuons[0]+tightMuons[3]
                  Z2 = tightMuons[1]+tightMuons[2]
                  Z1_lepton1 = tightMuons[0]
                  Z1_lepton2 = tightMuons[3]
                  Z2_lepton1 = tightMuons[1]
                  Z2_lepton2 = tightMuons[2]
                  self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[0])
                  self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[3])
                  self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[1])
                  self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[2])
                else:
                  Z1 = tightMuons[1]+tightMuons[2]
                  Z2 = tightMuons[0]+tightMuons[3]
                  Z1_lepton1 = tightMuons[1]
                  Z1_lepton2 = tightMuons[2]
                  Z2_lepton1 = tightMuons[0]
                  Z2_lepton2 = tightMuons[3]
                  self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[1])
                  self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[2])
                  self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
                  self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[3])

            else:
              if abs((tightMuons[0]+tightMuons[1]).M()-91.1876) < abs((tightMuons[2]+tightMuons[3]).M()-91.1876):
                Z1 = tightMuons[0]+tightMuons[1]
                Z2 = tightMuons[2]+tightMuons[3]
                Z1_lepton1 = tightMuons[0]
                Z1_lepton2 = tightMuons[1]
                Z2_lepton1 = tightMuons[2]
                Z2_lepton2 = tightMuons[3]
                self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[0])
                self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[1])
                self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[2])
                self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[3])
              else:
                Z1 = tightMuons[2]+tightMuons[3]
                Z2 = tightMuons[0]+tightMuons[1]
                Z1_lepton1 = tightMuons[2]
                Z1_lepton2 = tightMuons[3]
                Z2_lepton1 = tightMuons[0]
                Z2_lepton2 = tightMuons[1]
                self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[2])
                self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[3])
                self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[1])
          else:
            if pass_Z1_basic_2ndpair and pass_Z2_basic_2ndpair:
              if abs((tightMuons[0]+tightMuons[3]).M()-91.1876) < abs((tightMuons[1]+tightMuons[2]).M()-91.1876):
                Z1 = tightMuons[0]+tightMuons[3]
                Z2 = tightMuons[1]+tightMuons[2]
                Z1_lepton1 = tightMuons[0]
                Z1_lepton2 = tightMuons[3]
                Z2_lepton1 = tightMuons[1]
                Z2_lepton2 = tightMuons[2]
                self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[0])
                self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[3])
                self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[1])
                self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[2])
              else:
                Z1 = tightMuons[1]+tightMuons[2]
                Z2 = tightMuons[0]+tightMuons[3]
                Z1_lepton1 = tightMuons[1]
                Z1_lepton2 = tightMuons[2]
                Z2_lepton1 = tightMuons[0]
                Z2_lepton2 = tightMuons[3]
                self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[1])
                self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[2])
                self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[3])

      else: 
        pass_Z1_basic_1stpair = (tightMuons[0]+tightMuons[2]).M()>4 and (tightMuons[0]+tightMuons[2]).M()<120
        pass_Z2_basic_1stpair = (tightMuons[1]+tightMuons[3]).M()>4 and (tightMuons[1]+tightMuons[3]).M()<120
        pass_Z1_basic_2ndpair = (tightMuons[0]+tightMuons[3]).M()>4 and (tightMuons[0]+tightMuons[3]).M()<120
        pass_Z2_basic_2ndpair = (tightMuons[1]+tightMuons[2]).M()>4 and (tightMuons[1]+tightMuons[2]).M()<120

        if pass_Z1_basic_1stpair and pass_Z2_basic_1stpair:
          if pass_Z1_basic_2ndpair and pass_Z2_basic_2ndpair:
            if min(abs((tightMuons[0]+tightMuons[2]).M()-91.1876), abs((tightMuons[1]+tightMuons[3]).M()-91.1876)) < min(abs((tightMuons[0]+tightMuons[3]).M()-91.1876), abs((tightMuons[1]+tightMuons[2]).M()-91.1876)):
              if abs((tightMuons[0]+tightMuons[2]).M()-91.1876) < abs((tightMuons[1]+tightMuons[3]).M()-91.1876):
                Z1 = tightMuons[0]+tightMuons[2]
                Z2 = tightMuons[1]+tightMuons[3]
                Z1_lepton1 = tightMuons[0]
                Z1_lepton2 = tightMuons[2]
                Z2_lepton1 = tightMuons[1]
                Z2_lepton2 = tightMuons[3]
                self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[0])
                self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[2])
                self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[1])
                self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[3])
              else:
                Z1 = tightMuons[1]+tightMuons[3]
                Z2 = tightMuons[0]+tightMuons[2]
                Z1_lepton1 = tightMuons[1]
                Z1_lepton2 = tightMuons[3]
                Z2_lepton1 = tightMuons[0]
                Z2_lepton2 = tightMuons[2]
                self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[1])
                self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[3])
                self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[2])
            else:

              if abs((tightMuons[0]+tightMuons[3]).M()-91.1876) < abs((tightMuons[1]+tightMuons[2]).M()-91.1876):
                Z1 = tightMuons[0]+tightMuons[3]
                Z2 = tightMuons[1]+tightMuons[2]
                Z1_lepton1 = tightMuons[0]
                Z1_lepton2 = tightMuons[3]
                Z2_lepton1 = tightMuons[1]
                Z2_lepton2 = tightMuons[2]
                self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[0])
                self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[3])
                self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[1])
                self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[2])
              else:
                Z1 = tightMuons[1]+tightMuons[2]
                Z2 = tightMuons[0]+tightMuons[3]
                Z1_lepton1 = tightMuons[1]
                Z1_lepton2 = tightMuons[2]
                Z2_lepton1 = tightMuons[0]
                Z2_lepton2 = tightMuons[3]
                self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[1])
                self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[2])
                self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[3])

          else:
            if abs((tightMuons[0]+tightMuons[2]).M()-91.1876) < abs((tightMuons[1]+tightMuons[3]).M()-91.1876):
              Z1 = tightMuons[0]+tightMuons[2]
              Z2 = tightMuons[1]+tightMuons[3]
              Z1_lepton1 = tightMuons[0]
              Z1_lepton2 = tightMuons[2]
              Z2_lepton1 = tightMuons[1]
              Z2_lepton2 = tightMuons[3]
              self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[0])
              self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[2])
              self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[1])
              self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[3])
            else:
              Z1 = tightMuons[1]+tightMuons[3]
              Z2 = tightMuons[0]+tightMuons[2]
              Z1_lepton1 = tightMuons[1]
              Z1_lepton2 = tightMuons[3]
              Z2_lepton1 = tightMuons[0]
              Z2_lepton2 = tightMuons[2]
              self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[1])
              self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[3])
              self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
              self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[2])
        else:
          if pass_Z1_basic_2ndpair and pass_Z2_basic_2ndpair:
            if abs((tightMuons[0]+tightMuons[3]).M()-91.1876) < abs((tightMuons[1]+tightMuons[2]).M()-91.1876):
              Z1 = tightMuons[0]+tightMuons[3]
              Z2 = tightMuons[1]+tightMuons[2]
              Z1_lepton1 = tightMuons[0]
              Z1_lepton2 = tightMuons[3]
              Z2_lepton1 = tightMuons[1]
              Z2_lepton2 = tightMuons[2]
              self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[0])
              self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[3])
              self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[1])
              self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[2])
            else:
              Z1 = tightMuons[1]+tightMuons[2]
              Z2 = tightMuons[0]+tightMuons[3]
              Z1_lepton1 = tightMuons[1]
              Z1_lepton2 = tightMuons[2]
              Z2_lepton1 = tightMuons[0]
              Z2_lepton2 = tightMuons[3]
              self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[1])
              self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[2])
              self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
              self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[3])

    # 4 electrons case
    if len(tightMuons)==0:
      if tightElectrons_pdgid[0]+tightElectrons_pdgid[1] == 0:
        if tightElectrons_pdgid[0]+tightElectrons_pdgid[2] == 0:
          pass_Z1_basic_1stpair = (tightElectrons[0]+tightElectrons[1]).M()>4 and (tightElectrons[0]+tightElectrons[1]).M()<120
          pass_Z2_basic_1stpair = (tightElectrons[2]+tightElectrons[3]).M()>4 and (tightElectrons[2]+tightElectrons[3]).M()<120
          pass_Z1_basic_2ndpair = (tightElectrons[0]+tightElectrons[2]).M()>4 and (tightElectrons[0]+tightElectrons[2]).M()<120
          pass_Z2_basic_2ndpair = (tightElectrons[1]+tightElectrons[3]).M()>4 and (tightElectrons[1]+tightElectrons[3]).M()<120
          if pass_Z1_basic_1stpair and pass_Z2_basic_1stpair:
            if pass_Z1_basic_2ndpair and pass_Z2_basic_2ndpair:
              if min(abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876), abs((tightElectrons[2]+tightElectrons[3]).M()-91.1876)) < min(abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876), abs((tightElectrons[1]+tightElectrons[3]).M()-91.1876)):
                if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876) < abs((tightElectrons[2]+tightElectrons[3]).M()-91.1876):
                  Z1 = tightElectrons[0]+tightElectrons[1]
                  Z2 = tightElectrons[2]+tightElectrons[3]
                  Z1_lepton1 = tightElectrons[0]
                  Z1_lepton2 = tightElectrons[1]
                  Z2_lepton1 = tightElectrons[2]
                  Z2_lepton2 = tightElectrons[3]
                  self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[0])
                  self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[1])
                  self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[2])
                  self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[3])
                else:
                  Z1 = tightElectrons[2]+tightElectrons[3]
                  Z2 = tightElectrons[0]+tightElectrons[1]
                  Z1_lepton1 = tightElectrons[2]
                  Z1_lepton2 = tightElectrons[3]
                  Z2_lepton1 = tightElectrons[0]
                  Z2_lepton2 = tightElectrons[1]
                  self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[2])
                  self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[3])
                  self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
                  self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[1])
              else:

                if abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876) < abs((tightElectrons[1]+tightElectrons[3]).M()-91.1876):
                  Z1 = tightElectrons[0]+tightElectrons[2]
                  Z2 = tightElectrons[1]+tightElectrons[3]
                  Z1_lepton1 = tightElectrons[0]
                  Z1_lepton2 = tightElectrons[2]
                  Z2_lepton1 = tightElectrons[1]
                  Z2_lepton2 = tightElectrons[3]
                  self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[0])
                  self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[2])
                  self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[1])
                  self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[3])
                else:
                  Z1 = tightElectrons[1]+tightElectrons[3]
                  Z2 = tightElectrons[0]+tightElectrons[2]
                  Z1_lepton1 = tightElectrons[1]
                  Z1_lepton2 = tightElectrons[3]
                  Z2_lepton1 = tightElectrons[0]
                  Z2_lepton2 = tightElectrons[2]
                  self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[1])
                  self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[3])
                  self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
                  self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[2])

            else:
              if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876) < abs((tightElectrons[2]+tightElectrons[3]).M()-91.1876):
                Z1 = tightElectrons[0]+tightElectrons[1]
                Z2 = tightElectrons[2]+tightElectrons[3]
                Z1_lepton1 = tightElectrons[0]
                Z1_lepton2 = tightElectrons[1]
                Z2_lepton1 = tightElectrons[2]
                Z2_lepton2 = tightElectrons[3]
                self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[0])
                self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[1])
                self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[2])
                self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[3])
              else:
                Z1 = tightElectrons[2]+tightElectrons[3]
                Z2 = tightElectrons[0]+tightElectrons[1]
                Z1_lepton1 = tightElectrons[2]
                Z1_lepton2 = tightElectrons[3]
                Z2_lepton1 = tightElectrons[0]
                Z2_lepton2 = tightElectrons[1]
                self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[2])
                self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[3])
                self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[1])
          else:
            if pass_Z1_basic_2ndpair and pass_Z2_basic_2ndpair:
              if abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876) < abs((tightElectrons[1]+tightElectrons[3]).M()-91.1876):
                Z1 = tightElectrons[0]+tightElectrons[2]
                Z2 = tightElectrons[1]+tightElectrons[3]
                Z1_lepton1 = tightElectrons[0]
                Z1_lepton2 = tightElectrons[2]
                Z2_lepton1 = tightElectrons[1]
                Z2_lepton2 = tightElectrons[3]
                self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[0])
                self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[2])
                self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[1])
                self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[3])
              else:
                Z1 = tightElectrons[1]+tightElectrons[3]
                Z2 = tightElectrons[0]+tightElectrons[2]
                Z1_lepton1 = tightElectrons[1]
                Z1_lepton2 = tightElectrons[3]
                Z2_lepton1 = tightElectrons[0]
                Z2_lepton2 = tightElectrons[2]
                self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[1])
                self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[3])
                self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[2])

        else:
          pass_Z1_basic_1stpair = (tightElectrons[0]+tightElectrons[1]).M()>4 and (tightElectrons[0]+tightElectrons[1]).M()<120
          pass_Z2_basic_1stpair = (tightElectrons[2]+tightElectrons[3]).M()>4 and (tightElectrons[2]+tightElectrons[3]).M()<120
          pass_Z1_basic_2ndpair = (tightElectrons[0]+tightElectrons[3]).M()>4 and (tightElectrons[0]+tightElectrons[3]).M()<120
          pass_Z2_basic_2ndpair = (tightElectrons[1]+tightElectrons[2]).M()>4 and (tightElectrons[1]+tightElectrons[2]).M()<120

          if pass_Z1_basic_1stpair and pass_Z2_basic_1stpair:
            if pass_Z1_basic_2ndpair and pass_Z2_basic_2ndpair:
              if min(abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876), abs((tightElectrons[2]+tightElectrons[3]).M()-91.1876)) < min(abs((tightElectrons[0]+tightElectrons[3]).M()-91.1876), abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876)):
                if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876) < abs((tightElectrons[2]+tightElectrons[3]).M()-91.1876):
                  Z1 = tightElectrons[0]+tightElectrons[1]
                  Z2 = tightElectrons[2]+tightElectrons[3]
                  Z1_lepton1 = tightElectrons[0]
                  Z1_lepton2 = tightElectrons[1]
                  Z2_lepton1 = tightElectrons[2]
                  Z2_lepton2 = tightElectrons[3]
                  self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[0])
                  self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[1])
                  self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[2])
                  self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[3])
                else:
                  Z1 = tightElectrons[2]+tightElectrons[3]
                  Z2 = tightElectrons[0]+tightElectrons[1]
                  Z1_lepton1 = tightElectrons[2]
                  Z1_lepton2 = tightElectrons[3]
                  Z2_lepton1 = tightElectrons[0]
                  Z2_lepton2 = tightElectrons[1]
                  self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[2])
                  self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[3])
                  self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
                  self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[1])
              else:

                if abs((tightElectrons[0]+tightElectrons[3]).M()-91.1876) < abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876):
                  Z1 = tightElectrons[0]+tightElectrons[3]
                  Z2 = tightElectrons[1]+tightElectrons[2]
                  Z1_lepton1 = tightElectrons[0]
                  Z1_lepton2 = tightElectrons[3]
                  Z2_lepton1 = tightElectrons[1]
                  Z2_lepton2 = tightElectrons[2]
                  self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[0])
                  self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[3])
                  self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[1])
                  self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[2])
                else:
                  Z1 = tightElectrons[1]+tightElectrons[2]
                  Z2 = tightElectrons[0]+tightElectrons[3]
                  Z1_lepton1 = tightElectrons[1]
                  Z1_lepton2 = tightElectrons[2]
                  Z2_lepton1 = tightElectrons[0]
                  Z2_lepton2 = tightElectrons[3]
                  self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[1])
                  self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[2])
                  self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
                  self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[3])

            else:
              if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876) < abs((tightElectrons[2]+tightElectrons[3]).M()-91.1876):
                Z1 = tightElectrons[0]+tightElectrons[1]
                Z2 = tightElectrons[2]+tightElectrons[3]
                Z1_lepton1 = tightElectrons[0]
                Z1_lepton2 = tightElectrons[1]
                Z2_lepton1 = tightElectrons[2]
                Z2_lepton2 = tightElectrons[3]
                self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[0])
                self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[1])
                self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[2])
                self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[3])
              else:
                Z1 = tightElectrons[2]+tightElectrons[3]
                Z2 = tightElectrons[0]+tightElectrons[1]
                Z1_lepton1 = tightElectrons[2]
                Z1_lepton2 = tightElectrons[3]
                Z2_lepton1 = tightElectrons[0]
                Z2_lepton2 = tightElectrons[1]
                self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[2])
                self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[3])
                self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[1])
          else:
            if pass_Z1_basic_2ndpair and pass_Z2_basic_2ndpair:
              if abs((tightElectrons[0]+tightElectrons[3]).M()-91.1876) < abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876):
                Z1 = tightElectrons[0]+tightElectrons[3]
                Z2 = tightElectrons[1]+tightElectrons[2]
                Z1_lepton1 = tightElectrons[0]
                Z1_lepton2 = tightElectrons[3]
                Z2_lepton1 = tightElectrons[1]
                Z2_lepton2 = tightElectrons[2]
                self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[0])
                self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[3])
                self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[1])
                self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[2])
              else:
                Z1 = tightElectrons[1]+tightElectrons[2]
                Z2 = tightElectrons[0]+tightElectrons[3]
                Z1_lepton1 = tightElectrons[1]
                Z1_lepton2 = tightElectrons[2]
                Z2_lepton1 = tightElectrons[0]
                Z2_lepton2 = tightElectrons[3]
                self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[1])
                self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[2])
                self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[3])

      else: 
        pass_Z1_basic_1stpair = (tightElectrons[0]+tightElectrons[2]).M()>4 and (tightElectrons[0]+tightElectrons[2]).M()<120
        pass_Z2_basic_1stpair = (tightElectrons[1]+tightElectrons[3]).M()>4 and (tightElectrons[1]+tightElectrons[3]).M()<120
        pass_Z1_basic_2ndpair = (tightElectrons[0]+tightElectrons[3]).M()>4 and (tightElectrons[0]+tightElectrons[3]).M()<120
        pass_Z2_basic_2ndpair = (tightElectrons[1]+tightElectrons[2]).M()>4 and (tightElectrons[1]+tightElectrons[2]).M()<120

        if pass_Z1_basic_1stpair and pass_Z2_basic_1stpair:
          if pass_Z1_basic_2ndpair and pass_Z2_basic_2ndpair:
            if min(abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876), abs((tightElectrons[1]+tightElectrons[3]).M()-91.1876)) < min(abs((tightElectrons[0]+tightElectrons[3]).M()-91.1876), abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876)):
              if abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876) < abs((tightElectrons[1]+tightElectrons[3]).M()-91.1876):
                Z1 = tightElectrons[0]+tightElectrons[2]
                Z2 = tightElectrons[1]+tightElectrons[3]
                Z1_lepton1 = tightElectrons[0]
                Z1_lepton2 = tightElectrons[2]
                Z2_lepton1 = tightElectrons[1]
                Z2_lepton2 = tightElectrons[3]
                self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[0])
                self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[2])
                self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[1])
                self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[3])
              else:
                Z1 = tightElectrons[1]+tightElectrons[3]
                Z2 = tightElectrons[0]+tightElectrons[2]
                Z1_lepton1 = tightElectrons[1]
                Z1_lepton2 = tightElectrons[3]
                Z2_lepton1 = tightElectrons[0]
                Z2_lepton2 = tightElectrons[2]
                self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[1])
                self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[3])
                self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[2])
            else:

              if abs((tightElectrons[0]+tightElectrons[3]).M()-91.1876) < abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876):
                Z1 = tightElectrons[0]+tightElectrons[3]
                Z2 = tightElectrons[1]+tightElectrons[2]
                Z1_lepton1 = tightElectrons[0]
                Z1_lepton2 = tightElectrons[3]
                Z2_lepton1 = tightElectrons[1]
                Z2_lepton2 = tightElectrons[2]
                self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[0])
                self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[3])
                self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[1])
                self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[2])
              else:
                Z1 = tightElectrons[1]+tightElectrons[2]
                Z2 = tightElectrons[0]+tightElectrons[3]
                Z1_lepton1 = tightElectrons[1]
                Z1_lepton2 = tightElectrons[2]
                Z2_lepton1 = tightElectrons[0]
                Z2_lepton2 = tightElectrons[3]
                self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[1])
                self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[2])
                self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[3])

          else:
            if abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876) < abs((tightElectrons[1]+tightElectrons[3]).M()-91.1876):
              Z1 = tightElectrons[0]+tightElectrons[2]
              Z2 = tightElectrons[1]+tightElectrons[3]
              Z1_lepton1 = tightElectrons[0]
              Z1_lepton2 = tightElectrons[2]
              Z2_lepton1 = tightElectrons[1]
              Z2_lepton2 = tightElectrons[3]
              self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[0])
              self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[2])
              self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[1])
              self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[3])
            else:
              Z1 = tightElectrons[1]+tightElectrons[3]
              Z2 = tightElectrons[0]+tightElectrons[2]
              Z1_lepton1 = tightElectrons[1]
              Z1_lepton2 = tightElectrons[3]
              Z2_lepton1 = tightElectrons[0]
              Z2_lepton2 = tightElectrons[2]
              self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[1])
              self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[3])
              self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
              self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[2])
        else:
          if pass_Z1_basic_2ndpair and pass_Z2_basic_2ndpair:
            if abs((tightElectrons[0]+tightElectrons[3]).M()-91.1876) < abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876):
              Z1 = tightElectrons[0]+tightElectrons[3]
              Z2 = tightElectrons[1]+tightElectrons[2]
              Z1_lepton1 = tightElectrons[0]
              Z1_lepton2 = tightElectrons[3]
              Z2_lepton1 = tightElectrons[1]
              Z2_lepton2 = tightElectrons[2]
              self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[0])
              self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[3])
              self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[1])
              self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[2])
            else:
              Z1 = tightElectrons[1]+tightElectrons[2]
              Z2 = tightElectrons[0]+tightElectrons[3]
              Z1_lepton1 = tightElectrons[1]
              Z1_lepton2 = tightElectrons[2]
              Z2_lepton1 = tightElectrons[0]
              Z2_lepton2 = tightElectrons[3]
              self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[1])
              self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[2])
              self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
              self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[3])

    if not (Z1.M()>40):
      return False

    self.out.fillBranch("Z1_l1_pt", Z1_lepton1.Pt())
    self.out.fillBranch("Z1_l1_eta", Z1_lepton1.Eta())
    self.out.fillBranch("Z1_l1_phi", Z1_lepton1.Phi())
    self.out.fillBranch("Z1_l1_mass", Z1_lepton1.M())
    self.out.fillBranch("Z1_l2_pt", Z1_lepton2.Pt())
    self.out.fillBranch("Z1_l2_eta", Z1_lepton2.Eta())
    self.out.fillBranch("Z1_l2_phi", Z1_lepton2.Phi())
    self.out.fillBranch("Z1_l2_mass", Z1_lepton2.M())
    self.out.fillBranch("Z2_l1_pt", Z2_lepton1.Pt())
    self.out.fillBranch("Z2_l1_eta", Z2_lepton1.Eta())
    self.out.fillBranch("Z2_l1_phi", Z2_lepton1.Phi())
    self.out.fillBranch("Z2_l1_mass", Z2_lepton1.M())
    self.out.fillBranch("Z2_l2_pt", Z2_lepton2.Pt())
    self.out.fillBranch("Z2_l2_eta", Z2_lepton2.Eta())
    self.out.fillBranch("Z2_l2_phi", Z2_lepton2.Phi())
    self.out.fillBranch("Z2_l2_mass", Z2_lepton2.M())
    self.out.fillBranch("Z1_pt",Z1.Pt())
    self.out.fillBranch("Z1_eta",Z1.Eta())
    self.out.fillBranch("Z1_phi",Z1.Phi())
    self.out.fillBranch("Z1_mass",Z1.M())
    self.out.fillBranch("Z2_pt",Z2.Pt())
    self.out.fillBranch("Z2_eta",Z2.Eta())
    self.out.fillBranch("Z2_phi",Z2.Phi())
    self.out.fillBranch("Z2_mass",Z2.M())

    # photon selection
    # (MinPtCut,PhoSCEtaMultiRangeCut,PhoSingleTowerHadOverEmCut,PhoFull5x5SigmaIEtaIEtaCut,PhoGenericRhoPtScaledCut,PhoGenericRhoPtScaledCut,PhoGenericRhoPtScaledCut), 2 bits per cut, correspond to 3 IDs, loose, medium, tight
    # Ptcut is the right most bit
    # for each cut, loose:01, medium:1X, tight:11

    mask_medium_full = (1<<1) | (1<<3) | (1 << 5) | (1 << 7) | (1 << 9) | (1 << 11) | (1 << 13)

    photons = Collection(event, 'Photon')
    photon_v4temp=TLorentzVector()
    mediumPhotons = []
    mediumPhotons_i = []

    for ipho in range(0, event.nPhoton):
      if photons[ipho].pt < 10:
        continue
    
      if not (photons[ipho].isScEtaEE or photons[ipho].isScEtaEB):
        continue

      photon_bitmap = photons[ipho].vidNestedWPBitmap & mask_medium_full
 
      # cut-based medium ID
      if not (photon_bitmap==mask_medium_full):
        continue

      # pixel veto for electron mis-ID
      if photons[ipho].pixelSeed:
        continue

      # require DeltaR(photon, lepton)>0.3
      photon_v4temp.SetPtEtaPhiM(photons[ipho].pt,photons[ipho].eta,photons[ipho].phi,photons[ipho].mass)
      if not (photon_v4temp.DeltaR(Z1_lepton1)>0.3 and photon_v4temp.DeltaR(Z2_lepton1)>0.3 and photon_v4temp.DeltaR(Z1_lepton2)>0.3 and photon_v4temp.DeltaR(Z2_lepton2)>0.3):
        continue

      mediumPhotons.append(photon_v4temp.Clone())
      mediumPhotons_i.append(ipho)
    if not len(mediumPhotons)>0:
      return False

    self.out.fillBranch("photon_pt", mediumPhotons[0].Pt())
    self.out.fillBranch("photon_eta", mediumPhotons[0].Eta())
    self.out.fillBranch("photon_phi", mediumPhotons[0].Phi())
    self.out.fillBranch("photon_mass", mediumPhotons[0].M())
    self.out.fillBranch("photon_isScEtaEB", photons[mediumPhotons_i[0]].isScEtaEB)
    self.out.fillBranch("photon_isScEtaEE", photons[mediumPhotons_i[0]].isScEtaEE)
    
      
    return True

ZZG2016 = lambda: ZZGProducer("2016")
ZZG2017 = lambda: ZZGProducer("2017")
ZZG2018 = lambda: ZZGProducer("2018")
