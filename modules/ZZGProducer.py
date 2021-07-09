import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os
import numpy as np

class ZZGProducer(Module):
  def __init__( self , year ):
    self.year = year
  def beginJob(self):
    pass
  def endJob(self):
    pass
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch("HLT_passEle32WPTight", "I")
    self.out.branch("WZ_region", "I")
    self.out.branch("WZ_zl1_id", "I")
    self.out.branch("WZ_zl2_id", "I")
    self.out.branch("WZ_wl_id", "I")
    self.out.branch("ZZ_Z1_id", "I")
    self.out.branch("ZZ_Z2_id", "I")
    self.out.branch("ZZ_Z1_l1_id", "I")
    self.out.branch("ZZ_Z1_l2_id", "I")
    self.out.branch("ZZ_Z2_l1_id", "I")
    self.out.branch("ZZ_Z2_l2_id", "I")
    self.out.branch("ZZ_SR", "I")
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
    self.out.branch("matched_photon_pt", "F")
    self.out.branch("matched_photon_eta", "F")
    self.out.branch("matched_photon_phi", "F")
    self.out.branch("matched_photon_mass", "F")
    self.out.branch("matched_photon_isScEtaEB","B")
    self.out.branch("matched_photon_isScEtaEE","B")
    self.out.branch("unmatched_photon_pt", "F")
    self.out.branch("unmatched_photon_eta", "F")
    self.out.branch("unmatched_photon_phi", "F")
    self.out.branch("unmatched_photon_mass", "F")
    self.out.branch("unmatched_photon_isScEtaEB","B")
    self.out.branch("unmatched_photon_isScEtaEE","B")
    self.out.branch("tightJets_nob_CSVloose_id","I",lenVar="nJet")
    self.out.branch("tightJets_b_CSVloose_id","I",lenVar="nJet")
    self.out.branch("tightJets_nob_CSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightJets_b_CSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightJets_nob_CSVtight_id","I",lenVar="nJet")
    self.out.branch("tightJets_b_CSVtight_id","I",lenVar="nJet")
    self.out.branch("tightJets_nob_DeepCSVloose_id","I",lenVar="nJet")
    self.out.branch("tightJets_b_DeepCSVloose_id","I",lenVar="nJet")
    self.out.branch("tightJets_nob_DeepCSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightJets_b_DeepCSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightJets_nob_DeepCSVtight_id","I",lenVar="nJet")
    self.out.branch("tightJets_b_DeepCSVtight_id","I",lenVar="nJet")
    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))
  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    
    # PV selection
    if (event.PV_npvsGood<1): return False

    # trigger selection
    # special action for 2017 single ele HLT, https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2017DataRecommendations#Single_Electron_Triggers
    HLT_passEle32WPTight=0
    if self.year=="2017":
      trgobjs=Collection(event, 'TrigObj')
      if event.HLT_Ele32_WPTight_Gsf_L1DoubleEG==1:
	for iobj in range(0,event.nTrigObj):
	  if trgobjs[iobj].id==11 and (trgobjs[iobj].filterBits & (1<<10))== (1<<10):
	    HLT_passEle32WPTight=1

    self.out.fillBranch("HLT_passEle32WPTight",HLT_passEle32WPTight)

    # total number of ele+muon, currently require at least 3 leptons (WZ/ZZ/ttZ regions will be used)
    if ((event.nMuon + event.nElectron) < 3): return False

    # Muon selection: tight cut-based ID + tight PF iso, or loose cut-based ID + loose PF iso
    muons = Collection(event, 'Muon')
    muon_v4_temp=TLorentzVector()
    tightMuons = []
    tightMuons_pdgid = []
    tightMuons_id = []
    additional_looseMuons = []
    additional_looseMuons_pdgid = []
    additional_looseMuons_id = []
    for imu in range(0, event.nMuon):
      if (muons[imu].tightId):
        if (muons[imu].pfRelIso04_all<0.15 and abs(muons[imu].eta)<2.4 and event.Muon_corrected_pt[imu]>5):
          muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
          tightMuons.append(muon_v4_temp.Clone())
          tightMuons_pdgid.append(muons[imu].pdgId)
          tightMuons_id.append(imu)
      elif (muons[imu].looseId):
        if (muons[imu].pfRelIso04_all<0.25 and abs(muons[imu].eta)<2.4 and event.Muon_corrected_pt[imu]>5):
          muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
          additional_looseMuons.append(muon_v4_temp.Clone())
          additional_looseMuons_pdgid.append(muons[imu].pdgId)
          additional_looseMuons_id.append(imu)

    # electron selection: tight (veto) cut-based ID + impact parameter cut
    electrons = Collection(event, 'Electron')
    electron_v4_temp=TLorentzVector()
    tightElectrons = []
    tightElectrons_pdgid = []
    tightElectrons_id = []
    additional_vetoElectrons = []
    additional_vetoElectrons_pdgid = []
    additional_vetoElectrons_id = []
    for iele in range(0, event.nElectron):
      if (electrons[iele].cutBased==4):
        if (electrons[iele].tightCharge==2 and ((abs(electrons[iele].eta+electrons[iele].deltaEtaSC) <1.4442 and abs(electrons[iele].dxy)<0.05 and abs(electrons[iele].dz<0.1)) or (abs(electrons[iele].eta + electrons[iele].deltaEtaSC)>1.566 and abs(electrons[iele].eta + electrons[iele].deltaEtaSC)<2.4 and abs(electrons[iele].dxy)<0.1 and abs(electrons[iele].dz)<0.2)) and electrons[iele].pt>7):
          electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
          tightElectrons.append(electron_v4_temp.Clone())
          tightElectrons_pdgid.append(electrons[iele].pdgId)
          tightElectrons_id.append(iele)
      elif (electrons[iele].cutBased==1):
        if (((abs(electrons[iele].eta+electrons[iele].deltaEtaSC) <1.4442 and abs(electrons[iele].dxy)<0.05 and abs(electrons[iele].dz<0.1)) or (abs(electrons[iele].eta + electrons[iele].deltaEtaSC)>1.566 and abs(electrons[iele].eta + electrons[iele].deltaEtaSC)<2.4 and abs(electrons[iele].dxy)<0.1 and abs(electrons[iele].dz)<0.2)) and electrons[iele].pt>7):
          electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
          additional_vetoElectrons.append(electron_v4_temp.Clone())
          additional_vetoElectrons_pdgid.append(electrons[iele].pdgId)
          additional_vetoElectrons_id.append(iele)

    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
    # tight ak4 jets, 2016 (111=7), 2017/2018 (110=6), medium B-tag WP
    # DeepCSV=(nanoaod btagDeepB) loose: 0.1355, medium: 0.4506, tight: 0.7738
    # DeepFlavor=(nanoaod btagDeepFlavB) loose: 0.0532, medium: 0.3040, tight: 0.7476

    # c-jet tag is based on two-D cuts, medium DeepJet WP:
    # CvsL=btagDeepFlavCvL: 0.085, CvsB=btagDeepFlavCvB: 0.34
    # c-tag not available in NANOAOD yet
    jets = Collection(event, 'Jet')
    tightJets_nob_CSVloose_id = []
    tightJets_b_CSVloose_id = []
    tightJets_nob_CSVmedium_id = []
    tightJets_b_CSVmedium_id = []
    tightJets_nob_CSVtight_id = []
    tightJets_b_CSVtight_id = []

    tightJets_nob_DeepCSVloose_id = []
    tightJets_b_DeepCSVloose_id = []
    tightJets_nob_DeepCSVmedium_id = []
    tightJets_b_DeepCSVmedium_id = []
    tightJets_nob_DeepCSVtight_id = []
    tightJets_b_DeepCSVtight_id = []

    for ijet in range(0, event.nJet):
      if self.year=="2016":
        if jets[ijet].jetId==7 and jets[ijet].pt>30 and abs(jets[ijet].eta)<4.7: 
          if (jets[ijet].btagCSVV2 > 0.9693):
            tightJets_b_CSVtight_id.append(ijet)
          else:tightJets_nob_CSVtight_id.append(ijet)
          if (jets[ijet].btagCSVV2 > 0.8838):
            tightJets_b_CSVmedium_id.append(ijet)
          else:tightJets_nob_CSVmedium_id.append(ijet)
          if (jets[ijet].btagCSVV2 > 0.5803):
            tightJets_b_CSVloose_id.append(ijet)
          else:tightJets_nob_CSVloose_id.append(ijet)
          if (jets[ijet].btagDeepB > 0.8001):
            tightJets_b_DeepCSVtight_id.append(ijet)
          else:tightJets_nob_DeepCSVtight_id.append(ijet)
          if (jets[ijet].btagDeepB > 0.4941):
            tightJets_b_DeepCSVmedium_id.append(ijet)
          else:tightJets_nob_DeepCSVmedium_id.append(ijet)
          if (jets[ijet].btagDeepB > 0.1522):
            tightJets_b_DeepCSVloose_id.append(ijet)
          else:tightJets_nob_DeepCSVloose_id.append(ijet)

      elif (self.year=="2017" or self.year=="2018"):
	if jets[ijet].jetId==6 and jets[ijet].pt>30 and abs(jets[ijet].eta)<4.7:
          if (jets[ijet].btagCSVV2 > 0.9693):
            tightJets_b_CSVtight_id.append(ijet)
          else:tightJets_nob_CSVtight_id.append(ijet)
          if (jets[ijet].btagCSVV2 > 0.8838):
            tightJets_b_CSVmedium_id.append(ijet)
          else:tightJets_nob_CSVmedium_id.append(ijet)
          if (jets[ijet].btagCSVV2 > 0.5803):
            tightJets_b_CSVloose_id.append(ijet)
          else:tightJets_nob_CSVloose_id.append(ijet)
          if (jets[ijet].btagDeepB > 0.7738):
            tightJets_b_DeepCSVtight_id.append(ijet)
          else:tightJets_nob_DeepCSVtight_id.append(ijet)
          if (jets[ijet].btagDeepB > 0.4506):
            tightJets_b_DeepCSVmedium_id.append(ijet)
          else:tightJets_nob_DeepCSVmedium_id.append(ijet)
          if (jets[ijet].btagDeepB > 0.1355):
            tightJets_b_DeepCSVloose_id.append(ijet)
          else:tightJets_nob_DeepCSVloose_id.append(ijet)

    tightJets_b_CSVtight_id.extend(np.zeros(event.nJet-len(tightJets_b_CSVtight_id),int)-1)
    tightJets_nob_CSVtight_id.extend(np.zeros(event.nJet-len(tightJets_nob_CSVtight_id),int)-1)
    tightJets_b_CSVmedium_id.extend(np.zeros(event.nJet-len(tightJets_b_CSVmedium_id),int)-1)
    tightJets_nob_CSVmedium_id.extend(np.zeros(event.nJet-len(tightJets_nob_CSVmedium_id),int)-1)
    tightJets_b_CSVloose_id.extend(np.zeros(event.nJet-len(tightJets_b_CSVloose_id),int)-1)
    tightJets_nob_CSVloose_id.extend(np.zeros(event.nJet-len(tightJets_nob_CSVloose_id),int)-1)
    tightJets_b_DeepCSVtight_id.extend(np.zeros(event.nJet-len(tightJets_b_DeepCSVtight_id),int)-1)
    tightJets_nob_DeepCSVtight_id.extend(np.zeros(event.nJet-len(tightJets_nob_DeepCSVtight_id),int)-1)
    tightJets_b_DeepCSVmedium_id.extend(np.zeros(event.nJet-len(tightJets_b_DeepCSVmedium_id),int)-1)
    tightJets_nob_DeepCSVmedium_id.extend(np.zeros(event.nJet-len(tightJets_nob_DeepCSVmedium_id),int)-1)
    tightJets_b_DeepCSVloose_id.extend(np.zeros(event.nJet-len(tightJets_b_DeepCSVloose_id),int)-1)
    tightJets_nob_DeepCSVloose_id.extend(np.zeros(event.nJet-len(tightJets_nob_DeepCSVloose_id),int)-1)
    
    self.out.fillBranch("tightJets_nob_CSVloose_id",tightJets_nob_CSVloose_id)
    self.out.fillBranch("tightJets_b_CSVloose_id",tightJets_b_CSVloose_id)
    self.out.fillBranch("tightJets_nob_CSVmedium_id",tightJets_nob_CSVmedium_id)
    self.out.fillBranch("tightJets_b_CSVmedium_id",tightJets_b_CSVmedium_id)
    self.out.fillBranch("tightJets_nob_CSVtight_id",tightJets_nob_CSVtight_id)
    self.out.fillBranch("tightJets_b_CSVtight_id",tightJets_b_CSVtight_id)
    self.out.fillBranch("tightJets_nob_DeepCSVloose_id",tightJets_nob_DeepCSVloose_id)
    self.out.fillBranch("tightJets_b_DeepCSVloose_id",tightJets_b_DeepCSVloose_id)
    self.out.fillBranch("tightJets_nob_DeepCSVmedium_id",tightJets_nob_DeepCSVmedium_id)
    self.out.fillBranch("tightJets_b_DeepCSVmedium_id",tightJets_b_DeepCSVmedium_id)
    self.out.fillBranch("tightJets_nob_DeepCSVtight_id",tightJets_nob_DeepCSVtight_id)
    self.out.fillBranch("tightJets_b_DeepCSVtight_id",tightJets_b_DeepCSVtight_id)

    # tight leptons and additional loose leptons collection
    tightLeptons = tightMuons + tightElectrons
    tightLeptons.sort(key=lambda x: x.Pt(), reverse=True)
    looseLeptons = additional_looseMuons + additional_vetoElectrons
    looseLeptons.sort(key=lambda x: x.Pt(), reverse=True)

    # WW     WWW     WW ZZZZZZZZZ   region: only 3 tight leptons(pt>20), no b-jet, mll>4, |Z-91.1876|<15
    # WW     WWW     WW       ZZ            MET>30
    #  WW   WW WW   WW      ZZ
    #   WW WW   WW WW     ZZ
    #    WWW     WWW    ZZZZZZZZZ

    #WZ region lepton number selections
    WZ_nl=False
    #WZ region b-jet selection
    WZ_nb=False
    #WZ region lepton kinematics selctions
    WZ_leptons=False
    #WZ region MET selection
    WZ_MET=False
    #WZ region tag, 0: fail to pass the WZ selection, 1:3 muon, 2:2muon, 3:1muon, 4:0 muon
    WZ_region=0
    WZ_zl1_id=-1
    WZ_zl2_id=-1
    WZ_wl_id=-1

    if len(tightLeptons)==3 and tightLeptons[2].Pt()>20 and (len(looseLeptons)==0 or (len(looseLeptons)>0 and looseLeptons[0].Pt()<15)):
      WZ_nl=True
    if WZ_nl and tightJets_b_DeepCSVmedium_id[0]==-1:
      WZ_nb=True

    if WZ_nb and (tightLeptons[0]+tightLeptons[1]).M()>4 and (tightLeptons[2]+tightLeptons[1]).M()>4 and (tightLeptons[0]+tightLeptons[2]).M()>4:
      WZ_leptons=True
    
    if WZ_leptons and event.MET_pt>30:
      WZ_MET=True

    if WZ_MET:
      # 3 muons case
      if len(tightElectrons)==0 and abs(tightMuons_pdgid[0]+tightMuons_pdgid[1]+tightMuons_pdgid[2])==13:
	#two combination 0+2 or 1+2
        if (tightMuons_pdgid[0]-tightMuons_pdgid[1])==0:
          if abs((tightMuons[0]+tightMuons[2]).M()-91.1876)<abs((tightMuons[1]+tightMuons[2]).M()-91.1876) and abs((tightMuons[0]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[0]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[1]
	  if abs((tightMuons[0]+tightMuons[2]).M()-91.1876)>abs((tightMuons[1]+tightMuons[2]).M()-91.1876) and abs((tightMuons[1]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[1]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[0]
	#two combination 0+1 or 1+2
	elif (tightMuons_pdgid[0]-tightMuons_pdgid[2])==0:
	  if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<abs((tightMuons[1]+tightMuons[2]).M()-91.1876) and abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[0]
            WZ_zl2_id=tightMuons_id[1]
            WZ_wl_id=tightMuons_id[2]
	  if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)>abs((tightMuons[1]+tightMuons[2]).M()-91.1876) and abs((tightMuons[1]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[1]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[0]
	#two combination 0+1 or 0+2
	else:
	  if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<abs((tightMuons[0]+tightMuons[2]).M()-91.1876) and abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[0]
            WZ_zl2_id=tightMuons_id[1]
            WZ_wl_id=tightMuons_id[2]
	  if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)>abs((tightMuons[0]+tightMuons[2]).M()-91.1876) and abs((tightMuons[0]+tightMuons[2]).M()-91.1876)<15:
	    WZ_region=1
            WZ_zl1_id=tightMuons_id[0]
            WZ_zl2_id=tightMuons_id[2]
            WZ_wl_id=tightMuons_id[1]

      # 2 muons case
      if len(tightElectrons)==1 and (tightMuons_pdgid[0]-tightMuons_pdgid[1])==0:
	if abs((tightMuons[0]+tightMuons[1]).M()-91.1876)<15:
	  WZ_region=2
	  WZ_zl1_id=tightMuons_id[0]
	  WZ_zl2_id=tightMuons_id[1]
	  WZ_wl_id=tightElectrons_id[0]

      # 1 muon case
      if len(tightElectrons)==2 and (tightElectrons_pdgid[0]-tightElectrons_pdgid[1])==0:
	if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<15:
	  WZ_region=3
	  WZ_zl1_id=tightElectrons_id[0]
	  WZ_zl2_id=tightElectrons_id[1]
	  WZ_wl_id=tightMuons_id[0]

      # 0 muon case
      if len(tightElectrons)==3 and abs(tightElectrons_pdgid[0]+tightElectrons_pdgid[1]+tightElectrons_pdgid[2])==11:
	#two combination 0+2 or 1+2
        if (tightElectrons_pdgid[0]-tightElectrons_pdgid[1])==0:
          if abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876)<abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[0]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[1]
	  if abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876)>abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[1]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[0]
	#two combination 0+1 or 1+2
	elif (tightElectrons_pdgid[0]-tightElectrons_pdgid[2])==0:
	  if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[0]
            WZ_zl2_id=tightElectrons_id[1]
            WZ_wl_id=tightElectrons_id[2]
	  if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)>abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[1]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[1]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[0]
	#two combination 0+1 or 0+2
	else:
	  if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[0]
            WZ_zl2_id=tightElectrons_id[1]
            WZ_wl_id=tightElectrons_id[2]
	  if abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876)>abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876) and abs((tightElectrons[0]+tightElectrons[2]).M()-91.1876)<15:
	    WZ_region=4
            WZ_zl1_id=tightElectrons_id[0]
            WZ_zl2_id=tightElectrons_id[2]
            WZ_wl_id=tightElectrons_id[1]
	
    self.out.fillBranch("WZ_region", WZ_region)
    self.out.fillBranch("WZ_zl1_id", WZ_zl1_id)
    self.out.fillBranch("WZ_zl2_id", WZ_zl2_id)
    self.out.fillBranch("WZ_wl_id", WZ_wl_id)
    

    #  ZZZZZZZZZ  ZZZZZZZZZ  region: only 4 tight leptons, |Z-91.1876|<15
    #        ZZ         ZZ           
    #      ZZ         ZZ
    #    ZZ         ZZ
    #  ZZZZZZZZZ  ZZZZZZZZZ

    # and 

    #   tt        tt       ZZZZZZZZZ  region: 4 tight leptons, 1 Z, 2 bjets
    # tttttttt  tttttttt        ZZ           
    #   tt        tt           ZZ
    #   tt  tt    tt  tt      ZZ
    #   tttttt    tttttt   ZZZZZZZZZ

    ZZ_nl=False
    # ZZ_region: 1: 4ele, 2: 2mu+2ele, 3: 4mu
    ZZ_region=0
    ZZ_drll=False
    ZZ_mll=False
    ZZ_dremu=True

#    ttZ_bb=False
#    ttZ_nl=False
#    # ttZ_region: 1: 4ele, 2: 2mu+2ele, 3: 4mu, 4: 3ele+1mu, 5: 1ele+3mu
#    ttZ_region=0
#    ttZ_drll=False
#    ttZ_mll=False

#    # at lease 2 jets and at least 1 b jet in ttZ region
#    if tightJets_b_DeepCSVmedium_id[0]>-1 and tightJets_nob_DeepCSVmedium_id[0]>-1:
#      ttZ_1b=True

    if (len(tightElectrons) + len(tightMuons) ==4):
      ZZ_nl=True

    # 4 ele case: at least one ele with pt>20, at least 2 eles with pt>12
    if ZZ_nl and len(tightMuons)==0 and tightElectrons[0].Pt()>20 and tightElectrons[1].Pt()>12 and (tightElectrons_pdgid[0]+tightElectrons_pdgid[1]+tightElectrons_pdgid[2]+tightElectrons_pdgid[3] ==0):
      ZZ_region=1

    # 2mu-2ele case: at least one lep with pt>20, at least 2 leps with pt>12
    if ZZ_nl and len(tightMuons)==2 and tightLeptons[0].Pt()>20 and tightLeptons[1].Pt()>12 and (tightMuons_pdgid[0]+tightMuons_pdgid[1]==0) and (tightElectrons_pdgid[0]+tightElectrons_pdgid[1]==0):
      ZZ_region=2

    # 4 muon case: at least one muon with pt>20, at least 2 muon with pt>10
    if ZZ_nl and len(tightMuons)==4 and (tightMuons[0].Pt()>20 and tightMuons[1].Pt()>10) and (tightMuons_pdgid[0]+tightMuons_pdgid[1]+tightMuons_pdgid[2]+tightMuons_pdgid[3]==0):
      ZZ_region=3

    # require DeltaR(l1,l2)>0.02, following ZZto4L
    if ZZ_nl and (tightLeptons[0].DeltaR(tightLeptons[1])>0.02 and tightLeptons[0].DeltaR(tightLeptons[2])>0.02 and tightLeptons[0].DeltaR(tightLeptons[3])>0.02 and tightLeptons[1].DeltaR(tightLeptons[2])>0.02 and tightLeptons[1].DeltaR(tightLeptons[3])>0.02 and tightLeptons[2].DeltaR(tightLeptons[3])>0.02):
      ZZ_drll=True

    # require mll>4 regardless the flavor and charge
    if ZZ_nl and ((tightLeptons[0]+tightLeptons[1]).M()>4 and (tightLeptons[0]+tightLeptons[2]).M()>4 and (tightLeptons[0]+tightLeptons[3]).M()>4 and (tightLeptons[1]+tightLeptons[2]).M()>4 and (tightLeptons[1]+tightLeptons[3]).M()>4 and (tightLeptons[2]+tightLeptons[3]).M()>4):
      ZZ_mll=True

    # require DeltaR(ele, mu)>0.05 to remove spurious ghost leptons formed from ambiguities in track reconstruction, following ZZto4L
    if ZZ_nl and len(tightMuons)==2:
      if not (tightMuons[0].DeltaR(tightElectrons[0])>0.05 and tightMuons[0].DeltaR(tightElectrons[1])>0.05 and tightMuons[1].DeltaR(tightElectrons[0])>0.05 and tightMuons[1].DeltaR(tightElectrons[1])>0.05):
        ZZ_dremu=False

    Z1=TLorentzVector()
    Z2=TLorentzVector()
    Z1_lepton1=TLorentzVector()
    Z1_lepton2=TLorentzVector()
    Z2_lepton1=TLorentzVector()
    Z2_lepton2=TLorentzVector()

    # construct two Z boson, the leading Z is the one closer to Z mass
    # both Z mass need between (4, 120) GeV
    ZZ_2Zmass=False
    #ZZ_Z1_id: 26 means 2 muons, 22 means 2 eles
    ZZ_Z1_id=-99
    ZZ_Z2_id=-99
    ZZ_Z1_l1_id=-1
    ZZ_Z1_l2_id=-1
    ZZ_Z2_l1_id=-1
    ZZ_Z2_l2_id=-1
    ZZ_SR=-1

    # 2 muons and 2 eles case
    if ZZ_nl and ZZ_region>0 and len(tightMuons)==2:
      if ((tightMuons[0]+tightMuons[1]).M()>4 and (tightMuons[0]+tightMuons[1]).M()<120 and (tightElectrons[0]+tightElectrons[1]).M()>4 and (tightElectrons[0]+tightElectrons[1]).M()<120):
        ZZ_2Zmass=True
      if abs((tightMuons[0]+tightMuons[1]).M()-91.1876) < abs((tightElectrons[0]+tightElectrons[1]).M()-91.1876):
        Z1=tightMuons[0]+tightMuons[1]
        Z2=tightElectrons[0]+tightElectrons[1]
        Z1_lepton1 = tightMuons[0]
        Z1_lepton2 = tightMuons[1]
        Z2_lepton1 = tightElectrons[0]
        Z2_lepton2 = tightElectrons[1]
	ZZ_Z1_id = 26
	ZZ_Z2_id = 22
	ZZ_Z1_l1_id = tightMuons_id[0]
	ZZ_Z1_l2_id = tightMuons_id[1]
	ZZ_Z2_l1_id = tightElectrons_id[0]
	ZZ_Z2_l2_id = tightElectrons_id[1]
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
	ZZ_Z1_id = 22
	ZZ_Z2_id = 26
	ZZ_Z1_l1_id = tightElectrons_id[0]
	ZZ_Z1_l2_id = tightElectrons_id[1]
	ZZ_Z2_l1_id = tightMuons_id[0]
	ZZ_Z2_l2_id = tightMuons_id[1]
        self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[0])
        self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[1])
        self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
        self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[1])
        
    pass_Z1_basic_1stpair = False
    pass_Z2_basic_1stpair = False
    pass_Z1_basic_2ndpair = False
    pass_Z2_basic_2ndpair = False

    # 4 muons case
    if ZZ_nl and ZZ_region>0 and len(tightElectrons)==0:
      ZZ_Z1_id = 26
      ZZ_Z2_id = 26
      if tightMuons_pdgid[0]+tightMuons_pdgid[1] == 0:
        if tightMuons_pdgid[0]+tightMuons_pdgid[2] == 0:
          pass_Z1_basic_1stpair = (tightMuons[0]+tightMuons[1]).M()>4 and (tightMuons[0]+tightMuons[1]).M()<120
          pass_Z2_basic_1stpair = (tightMuons[2]+tightMuons[3]).M()>4 and (tightMuons[2]+tightMuons[3]).M()<120
          pass_Z1_basic_2ndpair = (tightMuons[0]+tightMuons[2]).M()>4 and (tightMuons[0]+tightMuons[2]).M()<120
          pass_Z2_basic_2ndpair = (tightMuons[1]+tightMuons[3]).M()>4 and (tightMuons[1]+tightMuons[3]).M()<120
	  ZZ_2Zmass=True
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
	          ZZ_Z1_l1_id = tightMuons_id[0]
	          ZZ_Z1_l2_id = tightMuons_id[1]
	          ZZ_Z2_l1_id = tightMuons_id[2]
	          ZZ_Z2_l2_id = tightMuons_id[3]
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
	          ZZ_Z1_l1_id = tightMuons_id[2]
	          ZZ_Z1_l2_id = tightMuons_id[3]
	          ZZ_Z2_l1_id = tightMuons_id[0]
	          ZZ_Z2_l2_id = tightMuons_id[1]
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
	          ZZ_Z1_l1_id = tightMuons_id[0]
	          ZZ_Z1_l2_id = tightMuons_id[2]
	          ZZ_Z2_l1_id = tightMuons_id[1]
	          ZZ_Z2_l2_id = tightMuons_id[3]
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
	          ZZ_Z1_l1_id = tightMuons_id[1]
	          ZZ_Z1_l2_id = tightMuons_id[3]
	          ZZ_Z2_l1_id = tightMuons_id[0]
	          ZZ_Z2_l2_id = tightMuons_id[2]
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
	        ZZ_Z1_l1_id = tightMuons_id[0]
	        ZZ_Z1_l2_id = tightMuons_id[1]
	        ZZ_Z2_l1_id = tightMuons_id[2]
	        ZZ_Z2_l2_id = tightMuons_id[3]
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
	        ZZ_Z1_l1_id = tightMuons_id[2]
	        ZZ_Z1_l2_id = tightMuons_id[3]
	        ZZ_Z2_l1_id = tightMuons_id[0]
	        ZZ_Z2_l2_id = tightMuons_id[1]
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
	        ZZ_Z1_l1_id = tightMuons_id[0]
	        ZZ_Z1_l2_id = tightMuons_id[2]
	        ZZ_Z2_l1_id = tightMuons_id[1]
	        ZZ_Z2_l2_id = tightMuons_id[3]
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
	        ZZ_Z1_l1_id = tightMuons_id[1]
	        ZZ_Z1_l2_id = tightMuons_id[3]
	        ZZ_Z2_l1_id = tightMuons_id[0]
	        ZZ_Z2_l2_id = tightMuons_id[2]
                self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[1])
                self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[3])
                self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[2])

        else:
          pass_Z1_basic_1stpair = (tightMuons[0]+tightMuons[1]).M()>4 and (tightMuons[0]+tightMuons[1]).M()<120
          pass_Z2_basic_1stpair = (tightMuons[2]+tightMuons[3]).M()>4 and (tightMuons[2]+tightMuons[3]).M()<120
          pass_Z1_basic_2ndpair = (tightMuons[0]+tightMuons[3]).M()>4 and (tightMuons[0]+tightMuons[3]).M()<120
          pass_Z2_basic_2ndpair = (tightMuons[1]+tightMuons[2]).M()>4 and (tightMuons[1]+tightMuons[2]).M()<120
	  ZZ_2Zmass=True
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
	          ZZ_Z1_l1_id = tightMuons_id[0]
	          ZZ_Z1_l2_id = tightMuons_id[1]
	          ZZ_Z2_l1_id = tightMuons_id[2]
	          ZZ_Z2_l2_id = tightMuons_id[3]
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
	          ZZ_Z1_l1_id = tightMuons_id[2]
	          ZZ_Z1_l2_id = tightMuons_id[3]
	          ZZ_Z2_l1_id = tightMuons_id[0]
	          ZZ_Z2_l2_id = tightMuons_id[1]
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
	          ZZ_Z1_l1_id = tightMuons_id[0]
	          ZZ_Z1_l2_id = tightMuons_id[3]
	          ZZ_Z2_l1_id = tightMuons_id[1]
	          ZZ_Z2_l2_id = tightMuons_id[2]
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
	          ZZ_Z1_l1_id = tightMuons_id[1]
	          ZZ_Z1_l2_id = tightMuons_id[2]
	          ZZ_Z2_l1_id = tightMuons_id[0]
	          ZZ_Z2_l2_id = tightMuons_id[3]
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
	        ZZ_Z1_l1_id = tightMuons_id[0]
	        ZZ_Z1_l2_id = tightMuons_id[1]
	        ZZ_Z2_l1_id = tightMuons_id[2]
	        ZZ_Z2_l2_id = tightMuons_id[3]
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
	        ZZ_Z1_l1_id = tightMuons_id[2]
	        ZZ_Z1_l2_id = tightMuons_id[3]
	        ZZ_Z2_l1_id = tightMuons_id[0]
	        ZZ_Z2_l2_id = tightMuons_id[1]
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
	        ZZ_Z1_l1_id = tightMuons_id[0]
	        ZZ_Z1_l2_id = tightMuons_id[3]
	        ZZ_Z2_l1_id = tightMuons_id[1]
	        ZZ_Z2_l2_id = tightMuons_id[2]
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
	        ZZ_Z1_l1_id = tightMuons_id[1]
	        ZZ_Z1_l2_id = tightMuons_id[2]
	        ZZ_Z2_l1_id = tightMuons_id[0]
	        ZZ_Z2_l2_id = tightMuons_id[3]
                self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[1])
                self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[2])
                self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[3])

      else: 
        pass_Z1_basic_1stpair = (tightMuons[0]+tightMuons[2]).M()>4 and (tightMuons[0]+tightMuons[2]).M()<120
        pass_Z2_basic_1stpair = (tightMuons[1]+tightMuons[3]).M()>4 and (tightMuons[1]+tightMuons[3]).M()<120
        pass_Z1_basic_2ndpair = (tightMuons[0]+tightMuons[3]).M()>4 and (tightMuons[0]+tightMuons[3]).M()<120
        pass_Z2_basic_2ndpair = (tightMuons[1]+tightMuons[2]).M()>4 and (tightMuons[1]+tightMuons[2]).M()<120
	ZZ_2Zmass=True
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
	        ZZ_Z1_l1_id = tightMuons_id[0]
	        ZZ_Z1_l2_id = tightMuons_id[2]
	        ZZ_Z2_l1_id = tightMuons_id[1]
	        ZZ_Z2_l2_id = tightMuons_id[3]
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
	        ZZ_Z1_l1_id = tightMuons_id[1]
	        ZZ_Z1_l2_id = tightMuons_id[3]
	        ZZ_Z2_l1_id = tightMuons_id[0]
	        ZZ_Z2_l2_id = tightMuons_id[2]
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
	        ZZ_Z1_l1_id = tightMuons_id[0]
	        ZZ_Z1_l2_id = tightMuons_id[3]
	        ZZ_Z2_l1_id = tightMuons_id[1]
	        ZZ_Z2_l2_id = tightMuons_id[2]
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
	        ZZ_Z1_l1_id = tightMuons_id[1]
	        ZZ_Z1_l2_id = tightMuons_id[2]
	        ZZ_Z2_l1_id = tightMuons_id[0]
	        ZZ_Z2_l2_id = tightMuons_id[3]
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
	      ZZ_Z1_l1_id = tightMuons_id[0]
	      ZZ_Z1_l2_id = tightMuons_id[2]
	      ZZ_Z2_l1_id = tightMuons_id[1]
	      ZZ_Z2_l2_id = tightMuons_id[3]
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
	      ZZ_Z1_l1_id = tightMuons_id[1]
	      ZZ_Z1_l2_id = tightMuons_id[3]
	      ZZ_Z2_l1_id = tightMuons_id[0]
	      ZZ_Z2_l2_id = tightMuons_id[2]
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
	      ZZ_Z1_l1_id = tightMuons_id[0]
	      ZZ_Z1_l2_id = tightMuons_id[3]
	      ZZ_Z2_l1_id = tightMuons_id[1]
	      ZZ_Z2_l2_id = tightMuons_id[2]
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
	      ZZ_Z1_l1_id = tightMuons_id[1]
	      ZZ_Z1_l2_id = tightMuons_id[2]
	      ZZ_Z2_l1_id = tightMuons_id[0]
	      ZZ_Z2_l2_id = tightMuons_id[3]
              self.out.fillBranch("Z1_l1_pdgId", tightMuons_pdgid[1])
              self.out.fillBranch("Z1_l2_pdgId", tightMuons_pdgid[2])
              self.out.fillBranch("Z2_l1_pdgId", tightMuons_pdgid[0])
              self.out.fillBranch("Z2_l2_pdgId", tightMuons_pdgid[3])

    # 4 electrons case
    if ZZ_nl and ZZ_region>0 and len(tightMuons)==0:
      ZZ_Z1_id = 22
      ZZ_Z2_id = 22
      if tightElectrons_pdgid[0]+tightElectrons_pdgid[1] == 0:
        if tightElectrons_pdgid[0]+tightElectrons_pdgid[2] == 0:
          pass_Z1_basic_1stpair = (tightElectrons[0]+tightElectrons[1]).M()>4 and (tightElectrons[0]+tightElectrons[1]).M()<120
          pass_Z2_basic_1stpair = (tightElectrons[2]+tightElectrons[3]).M()>4 and (tightElectrons[2]+tightElectrons[3]).M()<120
          pass_Z1_basic_2ndpair = (tightElectrons[0]+tightElectrons[2]).M()>4 and (tightElectrons[0]+tightElectrons[2]).M()<120
          pass_Z2_basic_2ndpair = (tightElectrons[1]+tightElectrons[3]).M()>4 and (tightElectrons[1]+tightElectrons[3]).M()<120
	  ZZ_2Zmass=True
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
	          ZZ_Z1_l1_id = tightElectrons_id[0]
	          ZZ_Z1_l2_id = tightElectrons_id[1]
	          ZZ_Z2_l1_id = tightElectrons_id[2]
	          ZZ_Z2_l2_id = tightElectrons_id[3]
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
	          ZZ_Z1_l1_id = tightElectrons_id[2]
	          ZZ_Z1_l2_id = tightElectrons_id[3]
	          ZZ_Z2_l1_id = tightElectrons_id[0]
	          ZZ_Z2_l2_id = tightElectrons_id[1]
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
	          ZZ_Z1_l1_id = tightElectrons_id[0]
	          ZZ_Z1_l2_id = tightElectrons_id[2]
	          ZZ_Z2_l1_id = tightElectrons_id[1]
	          ZZ_Z2_l2_id = tightElectrons_id[3]
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
	          ZZ_Z1_l1_id = tightElectrons_id[1]
	          ZZ_Z1_l2_id = tightElectrons_id[3]
	          ZZ_Z2_l1_id = tightElectrons_id[0]
	          ZZ_Z2_l2_id = tightElectrons_id[2]
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
	        ZZ_Z1_l1_id = tightElectrons_id[0]
	        ZZ_Z1_l2_id = tightElectrons_id[1]
	        ZZ_Z2_l1_id = tightElectrons_id[2]
	        ZZ_Z2_l2_id = tightElectrons_id[3]
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
	        ZZ_Z1_l1_id = tightElectrons_id[2]
	        ZZ_Z1_l2_id = tightElectrons_id[3]
	        ZZ_Z2_l1_id = tightElectrons_id[0]
	        ZZ_Z2_l2_id = tightElectrons_id[1]
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
	        ZZ_Z1_l1_id = tightElectrons_id[0]
	        ZZ_Z1_l2_id = tightElectrons_id[2]
	        ZZ_Z2_l1_id = tightElectrons_id[1]
	        ZZ_Z2_l2_id = tightElectrons_id[3]
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
	        ZZ_Z1_l1_id = tightElectrons_id[1]
	        ZZ_Z1_l2_id = tightElectrons_id[3]
	        ZZ_Z2_l1_id = tightElectrons_id[0]
	        ZZ_Z2_l2_id = tightElectrons_id[2]
                self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[1])
                self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[3])
                self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[2])

        else:
          pass_Z1_basic_1stpair = (tightElectrons[0]+tightElectrons[1]).M()>4 and (tightElectrons[0]+tightElectrons[1]).M()<120
          pass_Z2_basic_1stpair = (tightElectrons[2]+tightElectrons[3]).M()>4 and (tightElectrons[2]+tightElectrons[3]).M()<120
          pass_Z1_basic_2ndpair = (tightElectrons[0]+tightElectrons[3]).M()>4 and (tightElectrons[0]+tightElectrons[3]).M()<120
          pass_Z2_basic_2ndpair = (tightElectrons[1]+tightElectrons[2]).M()>4 and (tightElectrons[1]+tightElectrons[2]).M()<120
	  ZZ_2Zmass=True
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
	          ZZ_Z1_l1_id = tightElectrons_id[0]
	          ZZ_Z1_l2_id = tightElectrons_id[1]
	          ZZ_Z2_l1_id = tightElectrons_id[2]
	          ZZ_Z2_l2_id = tightElectrons_id[3]
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
	          ZZ_Z1_l1_id = tightElectrons_id[2]
	          ZZ_Z1_l2_id = tightElectrons_id[3]
	          ZZ_Z2_l1_id = tightElectrons_id[0]
	          ZZ_Z2_l2_id = tightElectrons_id[1]
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
	          ZZ_Z1_l1_id = tightElectrons_id[0]
	          ZZ_Z1_l2_id = tightElectrons_id[3]
	          ZZ_Z2_l1_id = tightElectrons_id[1]
	          ZZ_Z2_l2_id = tightElectrons_id[2]
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
	          ZZ_Z1_l1_id = tightElectrons_id[1]
	          ZZ_Z1_l2_id = tightElectrons_id[2]
	          ZZ_Z2_l1_id = tightElectrons_id[0]
	          ZZ_Z2_l2_id = tightElectrons_id[3]
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
	        ZZ_Z1_l1_id = tightElectrons_id[0]
	        ZZ_Z1_l2_id = tightElectrons_id[1]
	        ZZ_Z2_l1_id = tightElectrons_id[2]
	        ZZ_Z2_l2_id = tightElectrons_id[3]
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
	        ZZ_Z1_l1_id = tightElectrons_id[2]
	        ZZ_Z1_l2_id = tightElectrons_id[3]
	        ZZ_Z2_l1_id = tightElectrons_id[0]
	        ZZ_Z2_l2_id = tightElectrons_id[1]
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
	        ZZ_Z1_l1_id = tightElectrons_id[0]
	        ZZ_Z1_l2_id = tightElectrons_id[3]
	        ZZ_Z2_l1_id = tightElectrons_id[1]
	        ZZ_Z2_l2_id = tightElectrons_id[2]
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
	        ZZ_Z1_l1_id = tightElectrons_id[1]
	        ZZ_Z1_l2_id = tightElectrons_id[2]
	        ZZ_Z2_l1_id = tightElectrons_id[0]
	        ZZ_Z2_l2_id = tightElectrons_id[3]
                self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[1])
                self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[2])
                self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
                self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[3])

      else: 
        pass_Z1_basic_1stpair = (tightElectrons[0]+tightElectrons[2]).M()>4 and (tightElectrons[0]+tightElectrons[2]).M()<120
        pass_Z2_basic_1stpair = (tightElectrons[1]+tightElectrons[3]).M()>4 and (tightElectrons[1]+tightElectrons[3]).M()<120
        pass_Z1_basic_2ndpair = (tightElectrons[0]+tightElectrons[3]).M()>4 and (tightElectrons[0]+tightElectrons[3]).M()<120
        pass_Z2_basic_2ndpair = (tightElectrons[1]+tightElectrons[2]).M()>4 and (tightElectrons[1]+tightElectrons[2]).M()<120
	ZZ_2Zmass=True
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
	        ZZ_Z1_l1_id = tightElectrons_id[0]
	        ZZ_Z1_l2_id = tightElectrons_id[2]
	        ZZ_Z2_l1_id = tightElectrons_id[1]
	        ZZ_Z2_l2_id = tightElectrons_id[3]
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
	        ZZ_Z1_l1_id = tightElectrons_id[1]
	        ZZ_Z1_l2_id = tightElectrons_id[3]
	        ZZ_Z2_l1_id = tightElectrons_id[0]
	        ZZ_Z2_l2_id = tightElectrons_id[2]
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
	        ZZ_Z1_l1_id = tightElectrons_id[0]
	        ZZ_Z1_l2_id = tightElectrons_id[3]
	        ZZ_Z2_l1_id = tightElectrons_id[1]
	        ZZ_Z2_l2_id = tightElectrons_id[2]
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
	        ZZ_Z1_l1_id = tightElectrons_id[1]
	        ZZ_Z1_l2_id = tightElectrons_id[2]
	        ZZ_Z2_l1_id = tightElectrons_id[0]
	        ZZ_Z2_l2_id = tightElectrons_id[3]
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
	      ZZ_Z1_l1_id = tightElectrons_id[0]
	      ZZ_Z1_l2_id = tightElectrons_id[2]
	      ZZ_Z2_l1_id = tightElectrons_id[1]
	      ZZ_Z2_l2_id = tightElectrons_id[3]
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
	      ZZ_Z1_l1_id = tightElectrons_id[1]
	      ZZ_Z1_l2_id = tightElectrons_id[3]
	      ZZ_Z2_l1_id = tightElectrons_id[0]
	      ZZ_Z2_l2_id = tightElectrons_id[2]
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
	      ZZ_Z1_l1_id = tightElectrons_id[0]
	      ZZ_Z1_l2_id = tightElectrons_id[3]
	      ZZ_Z2_l1_id = tightElectrons_id[1]
	      ZZ_Z2_l2_id = tightElectrons_id[2]
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
	      ZZ_Z1_l1_id = tightElectrons_id[1]
	      ZZ_Z1_l2_id = tightElectrons_id[2]
	      ZZ_Z2_l1_id = tightElectrons_id[0]
	      ZZ_Z2_l2_id = tightElectrons_id[3]
              self.out.fillBranch("Z1_l1_pdgId", tightElectrons_pdgid[1])
              self.out.fillBranch("Z1_l2_pdgId", tightElectrons_pdgid[2])
              self.out.fillBranch("Z2_l1_pdgId", tightElectrons_pdgid[0])
              self.out.fillBranch("Z2_l2_pdgId", tightElectrons_pdgid[3])

    if ZZ_nl and ZZ_region>0 and ZZ_drll and ZZ_mll and ZZ_dremu and ZZ_2Zmass and Z1.M()>60 and Z1.M()<120 and Z2.M()>60 and Z2.M()<120:
      ZZ_SR=1

    if not (WZ_region>0 or ZZ_SR>0):
      return False

    self.out.fillBranch("ZZ_Z1_id",ZZ_Z1_id)
    self.out.fillBranch("ZZ_Z2_id",ZZ_Z2_id)
    self.out.fillBranch("ZZ_Z1_l1_id",ZZ_Z1_l1_id)
    self.out.fillBranch("ZZ_Z1_l2_id",ZZ_Z1_l2_id)
    self.out.fillBranch("ZZ_Z2_l1_id",ZZ_Z2_l1_id)
    self.out.fillBranch("ZZ_Z2_l2_id",ZZ_Z2_l2_id)
    self.out.fillBranch("ZZ_SR",ZZ_SR)
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
    # https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/PhotonIdentification/python/Identification/cutBasedPhotonID_tools.py#L333-L339
    # https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/PhotonIdentification/python/Identification/cutBasedPhotonID_Fall17_94X_V2_cff.py#L54
    # (MinPtCut,PhoSCEtaMultiRangeCut,PhoSingleTowerHadOverEmCut,PhoFull5x5SigmaIEtaIEtaCut,PhoGenericRhoPtScaledCut,PhoGenericRhoPtScaledCut,PhoGenericRhoPtScaledCut), 2 bits per cut, correspond to 3 IDs, loose, medium, tight
    # Ptcut is the right most bit
    # for each cut, loose:01, medium:1X, tight:11

    mask_medium_full = (1<<1) | (1<<3) | (1 << 5) | (1 << 7) | (1 << 9) | (1 << 11) | (1 << 13)

    photons = Collection(event, 'Photon')
    photon_v4temp=TLorentzVector()
    mediumPhotons_matched = []
    mediumPhotons_unmatched = []
    mediumPhotons_matched_id = []
    mediumPhotons_unmatched_id = []

    # medium ID photon pass full medium cuts
    print 'npho', event.nPhoton
    for ipho in range(0, event.nPhoton):
      if photons[ipho].pt < 15:
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
      print 'HAHAHA', photons[ipho].pt

      # require DeltaR(photon, lepton)>0.5
      photon_v4temp.SetPtEtaPhiM(photons[ipho].pt,photons[ipho].eta,photons[ipho].phi,photons[ipho].mass)
      if not (photon_v4temp.DeltaR(Z1_lepton1)>0.5 and photon_v4temp.DeltaR(Z2_lepton1)>0.5 and photon_v4temp.DeltaR(Z1_lepton2)>0.5 and photon_v4temp.DeltaR(Z2_lepton2)>0.5):
        continue

      # MC photon match to gen-level prompt photon
      if self.is_mc:
	if (photons[ipho].genPartFlav==1):
	  mediumPhotons_matched.append(photon_v4temp.Clone())
	  mediumPhotons_matched_id.append(ipho)
	else:
	  mediumPhotons_unmatched.append(photon_v4temp.Clone())
          mediumPhotons_unmatched_id.append(ipho)
      # use the same name for data
      else:
	mediumPhotons_matched.append(photon_v4temp.Clone())
	mediumPhotons_matched_id.append(ipho)

    if len(mediumPhotons_matched_id)==0:
      self.out.fillBranch("matched_photon_pt", -99)
      self.out.fillBranch("matched_photon_eta", -99)
      self.out.fillBranch("matched_photon_phi", -99)
      self.out.fillBranch("matched_photon_mass", -99)
      self.out.fillBranch("matched_photon_isScEtaEB", -99)
      self.out.fillBranch("matched_photon_isScEtaEE", -99)
    if len(mediumPhotons_matched_id)>0:
      self.out.fillBranch("matched_photon_pt", mediumPhotons_matched[0].Pt())
      self.out.fillBranch("matched_photon_eta", mediumPhotons_matched[0].Eta())
      self.out.fillBranch("matched_photon_phi", mediumPhotons_matched[0].Phi())
      self.out.fillBranch("matched_photon_mass", mediumPhotons_matched[0].M())
      self.out.fillBranch("matched_photon_isScEtaEB", photons[mediumPhotons_matched_id[0]].isScEtaEB)
      self.out.fillBranch("matched_photon_isScEtaEE", photons[mediumPhotons_matched_id[0]].isScEtaEE)

    if len(mediumPhotons_unmatched_id)==0:
      self.out.fillBranch("unmatched_photon_pt", -99)
      self.out.fillBranch("unmatched_photon_eta", -99)
      self.out.fillBranch("unmatched_photon_phi", -99)
      self.out.fillBranch("unmatched_photon_mass", -99)
      self.out.fillBranch("unmatched_photon_isScEtaEB", -99)
      self.out.fillBranch("unmatched_photon_isScEtaEE", -99)
    if len(mediumPhotons_unmatched_id)>0:
      self.out.fillBranch("unmatched_photon_pt", mediumPhotons_unmatched[0].Pt())
      self.out.fillBranch("unmatched_photon_eta", mediumPhotons_unmatched[0].Eta())
      self.out.fillBranch("unmatched_photon_phi", mediumPhotons_unmatched[0].Phi())
      self.out.fillBranch("unmatched_photon_mass", mediumPhotons_unmatched[0].M())
      self.out.fillBranch("unmatched_photon_isScEtaEB", photons[mediumPhotons_unmatched_id[0]].isScEtaEB)
      self.out.fillBranch("unmatched_photon_isScEtaEE", photons[mediumPhotons_unmatched_id[0]].isScEtaEE)

#    mask_sieie = (1<<1) | (1<<3) | (1 << 5) | (1 << 9) | (1 << 11) | (1 << 13)
#    mask_sieie_chiso = (1<<1) | (1<<3) | (1 << 5) | (1 << 11) | (1 << 13)
#
#    # store photon info for template fit
#    for ipho in range(0, event.nPhoton):
#      if photons[ipho].pt < 15:
#        continue
#      if not (photons[ipho].isScEtaEE or photons[ipho].isScEtaEB):
#        continue
#      photon_nosieie_bitmap = photons[ipho].vidNestedWPBitmap & mask_sieie
#      if not (photon_nosieie_bitmap==mask_sieie):
#        continue
#      # require DeltaR(photon, lepton)>0.5
#      photon_v4temp.SetPtEtaPhiM(photons[ipho].pt,photons[ipho].eta,photons[ipho].phi,photons[ipho].mass)
#      if not (photon_v4temp.DeltaR(Z1_lepton1)>0.5 and photon_v4temp.DeltaR(Z2_lepton1)>0.5 and photon_v4temp.DeltaR(Z1_lepton2)>0.5 and photon_v4temp.DeltaR(Z2_lepton2)>0.5):
#        continue
      
    return True

ZZG2016 = lambda: ZZGProducer("2016")
ZZG2017 = lambda: ZZGProducer("2017")
ZZG2018 = lambda: ZZGProducer("2018")
