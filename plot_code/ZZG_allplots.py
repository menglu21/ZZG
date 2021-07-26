import ROOT
import time
import os
import math
from math import sqrt
import plot_ZZGregion

#ZZG_header_path = os.path.join("aa.h")
#ROOT.gInterpreter.Declare('#include "{}"'.format(ZZG_header_path))

# the EnableImplicitMT option should only use in cluster, at lxplus, it will make the code slower(my experience)
#ROOT.ROOT.EnableImplicitMT()

def overunder_flowbin(h1):
  h1.SetBinContent(1,h1.GetBinContent(0)+h1.GetBinContent(1))
  h1.SetBinError(1,sqrt(h1.GetBinError(0)*h1.GetBinError(0)+h1.GetBinError(1)*h1.GetBinError(1)))
  h1.SetBinContent(h1.GetNbinsX(),h1.GetBinContent(h1.GetNbinsX())+h1.GetBinContent(h1.GetNbinsX()+1))
  h1.SetBinError(h1.GetNbinsX(),sqrt(h1.GetBinError(h1.GetNbinsX())*h1.GetBinError(h1.GetNbinsX())+h1.GetBinError(h1.GetNbinsX()+1)*h1.GetBinError(h1.GetNbinsX()+1)))
  return h1

def get_mcEventnumber(filename):
  print 'opening file ', filename
  ftemp=ROOT.TFile.Open(path+filename)
  htemp=ftemp.Get('nEventsGenWeighted')
  return htemp.GetBinContent(1)

def all_trigger(df):
  all_trigger = df.Filter("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_DoubleEle33_CaloIdL_MW || HLT_Ele35_WPTight_Gsf || HLT_Ele38_WPTight_Gsf || HLT_Ele40_WPTight_Gsf || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_TripleMu_10_5_5_DZ || HLT_TripleMu_12_10_5 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_TripleMu_10_5_5_DZ || HLT_TripleMu_12_10_5 || HLT_IsoMu27")
  return all_trigger

def for_diele_trigger(df):
  ditri_ele_trigger = df.Filter("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_DoubleEle33_CaloIdL_MW")
  return ditri_ele_trigger

def for_singleele_trigger(df):
  sin_ele_trigger = df.Filter("(HLT_passEle32WPTight || HLT_Ele35_WPTight_Gsf || HLT_Ele38_WPTight_Gsf || HLT_Ele40_WPTight_Gsf) && !(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_DoubleEle33_CaloIdL_MW)")
  return sin_ele_trigger

def for_dimuon_trigger(df):
  ditri_mu_trigger = df.Filter("(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_TripleMu_10_5_5_DZ || HLT_TripleMu_12_10_5) && !(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_DoubleEle33_CaloIdL_MW || HLT_Ele35_WPTight_Gsf || HLT_Ele38_WPTight_Gsf || HLT_Ele40_WPTight_Gsf || HLT_passEle32WPTight)")
  return ditri_mu_trigger

def for_cross_trigger(df):
  x_trigger = df.Filter("(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_DiMu9_Ele9_CaloIdL_TrackIdL || HLT_Mu8_DiEle12_CaloIdL_TrackIdL) && !(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_TripleMu_10_5_5_DZ || HLT_TripleMu_12_10_5 || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_DoubleEle33_CaloIdL_MW || HLT_Ele35_WPTight_Gsf || HLT_Ele38_WPTight_Gsf || HLT_Ele40_WPTight_Gsf || HLT_passEle32WPTight)")
  return x_trigger

def for_singlemuon_trigger(df):
  single_mu_trigger = df.Filter("HLT_IsoMu27 && !(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_DiMu9_Ele9_CaloIdL_TrackIdL || HLT_Mu8_DiEle12_CaloIdL_TrackIdL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_TripleMu_10_5_5_DZ || HLT_TripleMu_12_10_5 || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_DoubleEle33_CaloIdL_MW || HLT_Ele35_WPTight_Gsf || HLT_Ele38_WPTight_Gsf || HLT_Ele40_WPTight_Gsf || HLT_passEle32WPTight)")
  return single_mu_trigger

#path='/eos/cms/store/user/melu/ZZG2017/'
path='/eos/user/m/melu/ZZG2017/'

doubleMu_names = ROOT.std.vector('string')()
for f in ["DoubleMuonB.root","DoubleMuonC.root","DoubleMuonD.root","DoubleMuonE.root","DoubleMuonF.root"]:
  doubleMu_names.push_back(path+f)

singleMu_names = ROOT.std.vector('string')()
for f in ["SingleMuonB.root","SingleMuonC.root","SingleMuonD.root","SingleMuonE.root","SingleMuonF.root"]:
  singleMu_names.push_back(path+f)

doubleEle_names = ROOT.std.vector('string')()
for f in ["DoubleEGB.root","DoubleEGC.root","DoubleEGD.root","DoubleEGE.root","DoubleEGF.root"]:
  doubleEle_names.push_back(path+f)

singleEle_names = ROOT.std.vector('string')()
for f in ["SingleEGB.root","SingleEGC.root","SingleEGD.root","SingleEGE.root","SingleEGF.root"]:
  singleEle_names.push_back(path+f)

muonEle_names = ROOT.std.vector('string')()
for f in ["MuonEGB.root","MuonEGC.root","MuonEGD.root","MuonEGE.root","MuonEGF.root"]:
  muonEle_names.push_back(path+f)

signal_list = ['ZZG.root']
zz_list = ['ZZ.root']
ggzz_list = ['ggZZ_4e.root','ggZZ_4mu.root','ggZZ_4tau.root','ggZZ_2e2mu.root','ggZZ_2e2tau.root','ggZZ_2mu2tau.root']
top_list = ['TTZ.root','TTG.root']
vvv_list = ['WWZ.root','WZG.root','WZZ.root','ZZZ.root']
vv_list = ['WZTo3L.root','WZTo2L.root','ZG.root']

#histograms name
hists_name = ['Z1_l1_pt','Z1_l1_eta','Z1_l1_phi','Z1_l2_pt','Z1_l2_eta','Z1_l2_phi','Z2_l1_pt','Z2_l1_eta','Z2_l1_phi','Z2_l2_pt','Z2_l2_eta','Z2_l2_phi','Z1_mass','Z1_pt','Z1_eta','Z1_phi','Z2_mass','Z2_pt','Z2_eta','Z2_phi', 'photon_pt','photon_eta','photon_phi']

#histograms bins
histos_bins = {
hists_name[0]:10,
hists_name[1]:10,
hists_name[2]:10,
hists_name[3]:10,
hists_name[4]:10,
hists_name[5]:10,
hists_name[6]:10,
hists_name[7]:10,
hists_name[8]:10,
hists_name[9]:10,
hists_name[10]:10,
hists_name[11]:10,
hists_name[12]:20,
hists_name[13]:10,
hists_name[14]:10,
hists_name[15]:10,
hists_name[16]:20,
hists_name[17]:10,
hists_name[18]:10,
hists_name[19]:10,
hists_name[20]:10,
hists_name[21]:10,
hists_name[22]:10,
}

#low edge
histos_bins_low = {
hists_name[0]:0,
hists_name[1]:-3,
hists_name[2]:-4,
hists_name[3]:0,
hists_name[4]:-3,
hists_name[5]:-4,
hists_name[6]:0,
hists_name[7]:-3,
hists_name[8]:-4,
hists_name[9]:0,
hists_name[10]:-3,
hists_name[11]:-4,
hists_name[12]:60,
hists_name[13]:0,
hists_name[14]:-5,
hists_name[15]:-4,
hists_name[16]:60,
hists_name[17]:0,
hists_name[18]:-5,
hists_name[19]:-4,
hists_name[20]:0,
hists_name[21]:-2.5,
hists_name[22]:-4,
}

#high edge
histos_bins_high = {
hists_name[0]:200,
hists_name[1]:3,
hists_name[2]:4,
hists_name[3]:100,
hists_name[4]:3,
hists_name[5]:4,
hists_name[6]:200,
hists_name[7]:3,
hists_name[8]:4,
hists_name[9]:100,
hists_name[10]:3,
hists_name[11]:4,
hists_name[12]:120,
hists_name[13]:200,
hists_name[14]:5,
hists_name[15]:4,
hists_name[16]:120,
hists_name[17]:160,
hists_name[18]:5,
hists_name[19]:4,
hists_name[20]:100,
hists_name[21]:2.5,
hists_name[22]:4,
}

def ZZG_Analysis():

  histos = []

  lumi = 41480.

  ZZG_xs = 0.02814
  ZZG_ev = get_mcEventnumber(signal_list[0])

  ZZ_xs = 1.256*1.1
  ZZ_ev = get_mcEventnumber(zz_list[0])

  ggZZ_4l_xs = 0.00159*1.7
  ggZZ_4e_ev = get_mcEventnumber(ggzz_list[0])
  ggZZ_4mu_ev = get_mcEventnumber(ggzz_list[1])
  ggZZ_4tau_ev = get_mcEventnumber(ggzz_list[2])

  ggZZ_2l2l_xs = 0.00319*1.7
  ggZZ_2e2mu_ev = get_mcEventnumber(ggzz_list[3])
  ggZZ_2e2tau_ev = get_mcEventnumber(ggzz_list[4])
  ggZZ_2mu2tau_ev = get_mcEventnumber(ggzz_list[5])
  
  TTG_xs = 0.632
  TTG_ev = get_mcEventnumber(top_list[1])

  TTZ_xs = 0.253
  TTZ_ev = get_mcEventnumber(top_list[0])

  WZTo2L_xs = 12.178
  WZTo2L_ev = get_mcEventnumber(vv_list[1])

  WZTo3L_xs = 4.4297
  WZTo3L_ev = get_mcEventnumber(vv_list[0])

  ZG_xs = 123.9
  ZG_ev = get_mcEventnumber(vv_list[2])

  WZG_xs = 0.0384
  WZG_ev = get_mcEventnumber(vvv_list[1])

  WZZ_xs = 0.05565
  WZZ_ev = get_mcEventnumber(vvv_list[2])

  WWZ_xs = 0.1651
  WWZ_ev = get_mcEventnumber(vvv_list[0])

  ZZZ_xs = 0.01398
  ZZZ_ev = get_mcEventnumber(vvv_list[3])

  # define the filters here, 1:4e, 2:2e2m, 3:4m
  #filters="ZZ_region==2 && Z1_mass>60 && Z1_mass<120 && Z2_mass>60 && Z2_mass<120 && Muon_sip3d[0]<4&&Muon_sip3d[1]<4&&Electron_sip3d[0]<4&&Electron_sip3d[1]<4 && nPhoton>0"
  #filters="ZZ_region==1 && Z1_mass>60 && Z1_mass<120 && Z2_mass>60 && Z2_mass<120 && Electron_sip3d[0]<4&&Electron_sip3d[1]<4&&Electron_sip3d[2]<4&&Electron_sip3d[3]<4 && nPhoton>0"
  filters="ZZ_region==3 && Z1_mass>60 && Z1_mass<120 && Z2_mass>60 && Z2_mass<120 && Muon_sip3d[0]<4&&Muon_sip3d[1]<4&&Muon_sip3d[2]<4&&Muon_sip3d[3]<4 && nPhoton>0"

  photon_filter1="loosePhotons_matched_id[0]>-1 && "
  photon_filter2="loosePhotons_unmatched_id[0]>-1"

  df_ZZG_tree = ROOT.RDataFrame("Events",path+signal_list[0])
  df_ZZG_tree = df_ZZG_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_ZZG = df_ZZG_tree.Filter(filters)
  df_ZZG_photon = df_ZZG.Filter(photon_filter1)
  df_ZZG_trigger = all_trigger(df_ZZG_photon)
  df_ZZG_trigger = df_ZZG_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_ZZG_trigger = df_ZZG_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_ZZG_trigger = df_ZZG_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_ZZG_histos=[]
  for i in hists_name:
    df_ZZG_histo = df_ZZG_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_ZZG_histos.append(df_ZZG_histo)

  df_ZZ_tree = ROOT.RDataFrame("Events",path+zz_list[0])
  df_ZZ_tree = df_ZZ_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_ZZ = df_ZZ_tree.Filter(filters)
  df_ZZ_photon = df_ZZ.Filter(photon_filter2)
  df_ZZ_trigger = all_trigger(df_ZZ_photon)
  df_ZZ_trigger = df_ZZ_trigger.Define('photon_pt','Photon_pt[loosePhotons_unmatched_id[0]]')
  df_ZZ_trigger = df_ZZ_trigger.Define('photon_eta','Photon_eta[loosePhotons_unmatched_id[0]]')
  df_ZZ_trigger = df_ZZ_trigger.Define('photon_phi','Photon_phi[loosePhotons_unmatched_id[0]]')
  df_ZZ_histos=[]
  for i in hists_name:
    df_ZZ_histo = df_ZZ_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_ZZ_histos.append(df_ZZ_histo)

  df_ggZZ_4e_tree = ROOT.RDataFrame("Events",path+ggzz_list[0])
  df_ggZZ_4e_tree = df_ggZZ_4e_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_ggZZ_4e = df_ggZZ_4e_tree.Filter(filters)
  df_ggZZ_4e_photon = df_ggZZ_4e.Filter(photon_filter1)
  df_ggZZ_4e_trigger = all_trigger(df_ggZZ_4e_photon)
  df_ggZZ_4e_trigger = df_ggZZ_4e_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_ggZZ_4e_trigger = df_ggZZ_4e_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_ggZZ_4e_trigger = df_ggZZ_4e_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_ggZZ_4e_histos=[]
  for i in hists_name:
    df_ggZZ_4e_histo = df_ggZZ_4e_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_ggZZ_4e_histos.append(df_ggZZ_4e_histo)

  df_ggZZ_4mu_tree = ROOT.RDataFrame("Events",path+ggzz_list[1])
  df_ggZZ_4mu_tree = df_ggZZ_4mu_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_ggZZ_4mu = df_ggZZ_4mu_tree.Filter(filters)
  df_ggZZ_4mu_photon = df_ggZZ_4mu.Filter(photon_filter1)
  df_ggZZ_4mu_trigger = all_trigger(df_ggZZ_4mu_photon)
  df_ggZZ_4mu_trigger = df_ggZZ_4mu_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_ggZZ_4mu_trigger = df_ggZZ_4mu_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_ggZZ_4mu_trigger = df_ggZZ_4mu_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_ggZZ_4mu_histos=[]
  for i in hists_name:
    df_ggZZ_4mu_histo = df_ggZZ_4mu_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_ggZZ_4mu_histos.append(df_ggZZ_4mu_histo)

  df_ggZZ_4tau_tree = ROOT.RDataFrame("Events",path+ggzz_list[2])
  df_ggZZ_4tau_tree = df_ggZZ_4tau_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_ggZZ_4tau = df_ggZZ_4tau_tree.Filter(filters)
  df_ggZZ_4tau_photon = df_ggZZ_4tau.Filter(photon_filter1)
  df_ggZZ_4tau_trigger = all_trigger(df_ggZZ_4tau_photon)
  df_ggZZ_4tau_trigger = df_ggZZ_4tau_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_ggZZ_4tau_trigger = df_ggZZ_4tau_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_ggZZ_4tau_trigger = df_ggZZ_4tau_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_ggZZ_4tau_histos=[]
  for i in hists_name:
    df_ggZZ_4tau_histo = df_ggZZ_4tau_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_ggZZ_4tau_histos.append(df_ggZZ_4tau_histo)

  df_ggZZ_2e2mu_tree = ROOT.RDataFrame("Events",path+ggzz_list[3])
  df_ggZZ_2e2mu_tree = df_ggZZ_2e2mu_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_ggZZ_2e2mu = df_ggZZ_2e2mu_tree.Filter(filters)
  df_ggZZ_2e2mu_photon = df_ggZZ_2e2mu.Filter(photon_filter1)
  df_ggZZ_2e2mu_trigger = all_trigger(df_ggZZ_2e2mu_photon)
  df_ggZZ_2e2mu_trigger = df_ggZZ_2e2mu_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_ggZZ_2e2mu_trigger = df_ggZZ_2e2mu_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_ggZZ_2e2mu_trigger = df_ggZZ_2e2mu_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_ggZZ_2e2mu_histos=[]
  for i in hists_name:
    df_ggZZ_2e2mu_histo = df_ggZZ_2e2mu_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_ggZZ_2e2mu_histos.append(df_ggZZ_2e2mu_histo)

  df_ggZZ_2e2tau_tree = ROOT.RDataFrame("Events",path+ggzz_list[4])
  df_ggZZ_2e2tau_tree = df_ggZZ_2e2tau_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_ggZZ_2e2tau = df_ggZZ_2e2tau_tree.Filter(filters)
  df_ggZZ_2e2tau_photon = df_ggZZ_2e2tau.Filter(photon_filter1)
  df_ggZZ_2e2tau_trigger = all_trigger(df_ggZZ_2e2tau_photon)
  df_ggZZ_2e2tau_trigger = df_ggZZ_2e2tau_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_ggZZ_2e2tau_trigger = df_ggZZ_2e2tau_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_ggZZ_2e2tau_trigger = df_ggZZ_2e2tau_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_ggZZ_2e2tau_histos=[]
  for i in hists_name:
    df_ggZZ_2e2tau_histo = df_ggZZ_2e2tau_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_ggZZ_2e2tau_histos.append(df_ggZZ_2e2tau_histo)

  df_ggZZ_2mu2tau_tree = ROOT.RDataFrame("Events",path+ggzz_list[5])
  df_ggZZ_2mu2tau_tree = df_ggZZ_2mu2tau_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_ggZZ_2mu2tau = df_ggZZ_2mu2tau_tree.Filter(filters)
  df_ggZZ_2mu2tau_photon = df_ggZZ_2mu2tau.Filter(photon_filter1)
  df_ggZZ_2mu2tau_trigger = all_trigger(df_ggZZ_2mu2tau_photon)
  df_ggZZ_2mu2tau_trigger = df_ggZZ_2mu2tau_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_ggZZ_2mu2tau_trigger = df_ggZZ_2mu2tau_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_ggZZ_2mu2tau_trigger = df_ggZZ_2mu2tau_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_ggZZ_2mu2tau_histos=[]
  for i in hists_name:
    df_ggZZ_2mu2tau_histo = df_ggZZ_2mu2tau_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_ggZZ_2mu2tau_histos.append(df_ggZZ_2mu2tau_histo)

  df_TTZ_tree = ROOT.RDataFrame("Events",path+top_list[0])
  df_TTZ_tree = df_TTZ_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_TTZ = df_TTZ_tree.Filter(filters)
  df_TTZ_photon = df_TTZ.Filter(photon_filter1)
  df_TTZ_trigger = all_trigger(df_TTZ_photon)
  df_TTZ_trigger = df_TTZ_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_TTZ_trigger = df_TTZ_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_TTZ_trigger = df_TTZ_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_TTZ_histos=[]
  for i in hists_name:
    df_TTZ_histo = df_TTZ_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_TTZ_histos.append(df_TTZ_histo)

  df_TTG_tree = ROOT.RDataFrame("Events",path+top_list[1])
  df_TTG_tree = df_TTG_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_TTG = df_TTG_tree.Filter(filters)
  df_TTG_photon = df_TTG.Filter(photon_filter1)
  df_TTG_trigger = all_trigger(df_TTG_photon)
  df_TTG_trigger = df_TTG_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_TTG_trigger = df_TTG_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_TTG_trigger = df_TTG_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_TTG_histos=[]
  for i in hists_name:
    df_TTG_histo = df_TTG_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_TTG_histos.append(df_TTG_histo)

  df_WWZ_tree = ROOT.RDataFrame("Events",path+vvv_list[0])
  df_WWZ_tree = df_WWZ_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_WWZ = df_WWZ_tree.Filter(filters)
  df_WWZ_photon = df_WWZ.Filter(photon_filter1)
  df_WWZ_trigger = all_trigger(df_WWZ_photon)
  df_WWZ_trigger = df_WWZ_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_WWZ_trigger = df_WWZ_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_WWZ_trigger = df_WWZ_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_WWZ_histos=[]
  for i in hists_name:
    df_WWZ_histo = df_WWZ_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_WWZ_histos.append(df_WWZ_histo)

  df_WZG_tree = ROOT.RDataFrame("Events",path+vvv_list[1])
  df_WZG_tree = df_WZG_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_WZG = df_WZG_tree.Filter(filters)
  df_WZG_photon = df_WZG.Filter(photon_filter1)
  df_WZG_trigger = all_trigger(df_WZG_photon)
  df_WZG_trigger = df_WZG_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_WZG_trigger = df_WZG_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_WZG_trigger = df_WZG_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_WZG_histos=[]
  for i in hists_name:
    df_WZG_histo = df_WZG_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_WZG_histos.append(df_WZG_histo)

  df_WZZ_tree = ROOT.RDataFrame("Events",path+vvv_list[2])
  df_WZZ_tree = df_WZZ_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_WZZ = df_WZZ_tree.Filter(filters)
  df_WZZ_photon = df_WZZ.Filter(photon_filter1)
  df_WZZ_trigger = all_trigger(df_WZZ_photon)
  df_WZZ_trigger = df_WZZ_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_WZZ_trigger = df_WZZ_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_WZZ_trigger = df_WZZ_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_WZZ_histos=[]
  for i in hists_name:
    df_WZZ_histo = df_WZZ_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_WZZ_histos.append(df_WZZ_histo)

  df_ZZZ_tree = ROOT.RDataFrame("Events",path+vvv_list[3])
  df_ZZZ_tree = df_ZZZ_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_ZZZ = df_ZZZ_tree.Filter(filters)
  df_ZZZ_photon = df_ZZZ.Filter(photon_filter1)
  df_ZZZ_trigger = all_trigger(df_ZZZ_photon)
  df_ZZZ_trigger = df_ZZZ_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_ZZZ_trigger = df_ZZZ_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_ZZZ_trigger = df_ZZZ_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_ZZZ_histos=[]
  for i in hists_name:
    df_ZZZ_histo = df_ZZZ_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_ZZZ_histos.append(df_ZZZ_histo)

  df_WZTo3L_tree = ROOT.RDataFrame("Events",path+vv_list[0])
  df_WZTo3L_tree = df_WZTo3L_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_WZTo3L = df_WZTo3L_tree.Filter(filters)
  df_WZTo3L_photon = df_WZTo3L.Filter(photon_filter2)
  df_WZTo3L_trigger = all_trigger(df_WZTo3L_photon)
  df_WZTo3L_trigger = df_WZTo3L_trigger.Define('photon_pt','Photon_pt[loosePhotons_unmatched_id[0]]')
  df_WZTo3L_trigger = df_WZTo3L_trigger.Define('photon_eta','Photon_eta[loosePhotons_unmatched_id[0]]')
  df_WZTo3L_trigger = df_WZTo3L_trigger.Define('photon_phi','Photon_phi[loosePhotons_unmatched_id[0]]')
  df_WZTo3L_histos=[]
  for i in hists_name:
    df_WZTo3L_histo = df_WZTo3L_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_WZTo3L_histos.append(df_WZTo3L_histo)

  df_WZTo2L_tree = ROOT.RDataFrame("Events",path+vv_list[1])
  df_WZTo2L_tree = df_WZTo2L_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_WZTo2L = df_WZTo2L_tree.Filter(filters)
  df_WZTo2L_photon = df_WZTo2L.Filter(photon_filter2)
  df_WZTo2L_trigger = all_trigger(df_WZTo2L_photon)
  df_WZTo2L_trigger = df_WZTo2L_trigger.Define('photon_pt','Photon_pt[loosePhotons_unmatched_id[0]]')
  df_WZTo2L_trigger = df_WZTo2L_trigger.Define('photon_eta','Photon_eta[loosePhotons_unmatched_id[0]]')
  df_WZTo2L_trigger = df_WZTo2L_trigger.Define('photon_phi','Photon_phi[loosePhotons_unmatched_id[0]]')
  df_WZTo2L_histos=[]
  for i in hists_name:
    df_WZTo2L_histo = df_WZTo2L_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_WZTo2L_histos.append(df_WZTo2L_histo)

  df_ZG_tree = ROOT.RDataFrame("Events",path+vv_list[2])
  df_ZG_tree = df_ZG_tree.Define("genweight","puWeight*genWeight/abs(genWeight)")
  df_ZG = df_ZG_tree.Filter(filters)
  df_ZG_photon = df_ZG.Filter(photon_filter1)
  df_ZG_trigger = all_trigger(df_ZG_photon)
  df_ZG_trigger = df_ZG_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_ZG_trigger = df_ZG_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_ZG_trigger = df_ZG_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_ZG_histos=[]
  for i in hists_name:
    df_ZG_histo = df_ZG_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i,'genweight')
    df_ZG_histos.append(df_ZG_histo)

  df_DoubleMu_tree = ROOT.RDataFrame("Events", doubleMu_names)
  df_DoubleMu = df_DoubleMu_tree.Filter(filters)
  df_DoubleMu_photon = df_DoubleMu.Filter(photon_filter1)
  df_DoubleMu_trigger = for_dimuon_trigger(df_DoubleMu_photon) 
  df_DoubleMu_trigger = df_DoubleMu_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_DoubleMu_trigger = df_DoubleMu_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_DoubleMu_trigger = df_DoubleMu_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_DoubleMu_histos=[]
  for i in hists_name:
    df_DoubleMu_histo = df_DoubleMu_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i)
    df_DoubleMu_histos.append(df_DoubleMu_histo)

  df_SingleMu_tree = ROOT.RDataFrame("Events", singleMu_names)
  df_SingleMu = df_SingleMu_tree.Filter(filters)
  df_SingleMu_photon = df_SingleMu.Filter(photon_filter1)
  df_SingleMu_trigger = for_singlemuon_trigger(df_SingleMu_photon) 
  df_SingleMu_trigger = df_SingleMu_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_SingleMu_trigger = df_SingleMu_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_SingleMu_trigger = df_SingleMu_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_SingleMu_histos=[]
  for i in hists_name:
    df_SingleMu_histo = df_SingleMu_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i)
    df_SingleMu_histos.append(df_SingleMu_histo)

  df_DoubleEle_tree = ROOT.RDataFrame("Events", doubleEle_names)
  df_DoubleEle = df_DoubleEle_tree.Filter(filters)
  df_DoubleEle_photon = df_DoubleEle.Filter(photon_filter1)
  df_DoubleEle_trigger = for_diele_trigger(df_DoubleEle_photon)
  df_DoubleEle_trigger = df_DoubleEle_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_DoubleEle_trigger = df_DoubleEle_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_DoubleEle_trigger = df_DoubleEle_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_DoubleEle_histos=[]
  for i in hists_name:
    df_DoubleEle_histo = df_DoubleEle_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i)
    df_DoubleEle_histos.append(df_DoubleEle_histo)

  df_SingleEle_tree = ROOT.RDataFrame("Events", singleEle_names)
  df_SingleEle = df_SingleEle_tree.Filter(filters)
  df_SingleEle_photon = df_SingleEle.Filter(photon_filter1)
  df_SingleEle_trigger = for_singleele_trigger(df_SingleEle_photon)
  df_SingleEle_trigger = df_SingleEle_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_SingleEle_trigger = df_SingleEle_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_SingleEle_trigger = df_SingleEle_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_SingleEle_histos=[]
  for i in hists_name:
    df_SingleEle_histo = df_SingleEle_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i)
    df_SingleEle_histos.append(df_SingleEle_histo)

  df_MuonEle_tree = ROOT.RDataFrame("Events", muonEle_names)
  df_MuonEle = df_MuonEle_tree.Filter(filters)
  df_MuonEle_photon = df_MuonEle.Filter(photon_filter1)
  df_MuonEle_trigger = for_cross_trigger(df_MuonEle_photon)
  df_MuonEle_trigger = df_MuonEle_trigger.Define('photon_pt','Photon_pt[loosePhotons_matched_id[0]]')
  df_MuonEle_trigger = df_MuonEle_trigger.Define('photon_eta','Photon_eta[loosePhotons_matched_id[0]]')
  df_MuonEle_trigger = df_MuonEle_trigger.Define('photon_phi','Photon_phi[loosePhotons_matched_id[0]]')
  df_MuonEle_histos=[]
  for i in hists_name:
    df_MuonEle_histo = df_MuonEle_trigger.Histo1D((i,'',histos_bins[i],histos_bins_low[i],histos_bins_high[i]), i)
    df_MuonEle_histos.append(df_MuonEle_histo)

  for ij in range(0,len(hists_name)):
    df_ZZG_histos[ij].Draw()
    df_ZZ_histos[ij].Draw()
    df_ggZZ_4e_histos[ij].Draw()
    df_ggZZ_4mu_histos[ij].Draw()
    df_ggZZ_4tau_histos[ij].Draw()
    df_ggZZ_2e2mu_histos[ij].Draw()
    df_ggZZ_2e2tau_histos[ij].Draw()
    df_ggZZ_2mu2tau_histos[ij].Draw()
    df_TTZ_histos[ij].Draw()
    df_TTG_histos[ij].Draw()
    df_WWZ_histos[ij].Draw()
    df_WZG_histos[ij].Draw()
    df_WZZ_histos[ij].Draw()
    df_ZZZ_histos[ij].Draw()
    df_WZTo3L_histos[ij].Draw()
    df_WZTo2L_histos[ij].Draw()
    df_ZG_histos[ij].Draw()
    df_DoubleMu_histos[ij].Draw()
    df_SingleMu_histos[ij].Draw()
    df_DoubleEle_histos[ij].Draw()
    df_SingleEle_histos[ij].Draw()
    df_MuonEle_histos[ij].Draw()

# ROOT version 6.14 don;t have function "ROOT.RDF.RunGraphs"
#  ROOT.RDF.RunGraphs({df_ZZG_histo, df_ZZ_histo, df_ggZZ_4e_histo,df_ggZZ_4mu_histo, df_ggZZ_4tau_histo, df_ggZZ_2e2mu_histo,df_ggZZ_2e2tau_histo, df_ggZZ_2mu2tau_histo, df_TTZ_histo,df_TTG_histo, df_WWZ_histo, df_WZG_histo,df_WZZ_histo, df_ZZZ_histo, df_WZTo3L_histo,df_WZTo2L_histo, df_ZG_histo})

    h_ZZG = df_ZZG_histos[ij].GetValue()
    h_ZZ = df_ZZ_histos[ij].GetValue()
    h_ggZZ_4e = df_ggZZ_4e_histos[ij].GetValue()
    h_ggZZ_4mu = df_ggZZ_4mu_histos[ij].GetValue()
    h_ggZZ_4tau = df_ggZZ_4tau_histos[ij].GetValue()
    h_ggZZ_2e2mu = df_ggZZ_2e2mu_histos[ij].GetValue()
    h_ggZZ_2e2tau = df_ggZZ_2e2tau_histos[ij].GetValue()
    h_ggZZ_2mu2tau = df_ggZZ_2mu2tau_histos[ij].GetValue()
    h_TTZ = df_TTZ_histos[ij].GetValue()
    h_TTG = df_TTG_histos[ij].GetValue()
    h_WWZ = df_WWZ_histos[ij].GetValue()
    h_WZG = df_WZG_histos[ij].GetValue()
    h_WZZ = df_WZZ_histos[ij].GetValue()
    h_ZZZ = df_ZZZ_histos[ij].GetValue()
    h_WZTo3L = df_WZTo3L_histos[ij].GetValue()
    h_WZTo2L = df_WZTo2L_histos[ij].GetValue()
    h_ZG = df_ZG_histos[ij].GetValue()
    h_DoubleMu = df_DoubleMu_histos[ij].GetValue()
    h_SingleMu = df_SingleMu_histos[ij].GetValue()
    h_DoubleEle = df_DoubleEle_histos[ij].GetValue()
    h_SingleEle = df_SingleEle_histos[ij].GetValue()
    h_MuonEle = df_MuonEle_histos[ij].GetValue()

    h_ZZG.Scale(ZZG_xs/ZZG_ev)
    h_ZZ.Scale(ZZ_xs/ZZ_ev)
    h_ggZZ_4e.Scale(ggZZ_4l_xs/ggZZ_4e_ev)
    h_ggZZ_4mu.Scale(ggZZ_4l_xs/ggZZ_4mu_ev)
    h_ggZZ_4tau.Scale(ggZZ_4l_xs/ggZZ_4tau_ev)
    h_ggZZ_2e2mu.Scale(ggZZ_2l2l_xs/ggZZ_2e2mu_ev)
    h_ggZZ_2e2tau.Scale(ggZZ_2l2l_xs/ggZZ_2e2tau_ev)
    h_ggZZ_2mu2tau.Scale(ggZZ_2l2l_xs/ggZZ_2mu2tau_ev)
    h_TTZ.Scale(TTZ_xs/TTZ_ev)
    h_TTG.Scale(TTG_xs/TTG_ev)
    h_WWZ.Scale(WWZ_xs/WWZ_ev)
    h_WZG.Scale(WZG_xs/WZG_ev)
    h_WZZ.Scale(WZZ_xs/WZZ_ev)
    h_ZZZ.Scale(ZZZ_xs/ZZZ_ev)
    h_WZTo3L.Scale(WZTo3L_xs/WZTo3L_ev)
    h_WZTo2L.Scale(WZTo3L_xs/WZTo3L_ev)
    h_ZG.Scale(ZG_xs/ZG_ev)

    histos.append(h_ZZG.Clone())
    histos.append(h_ZZ.Clone())
    histos.append(h_ggZZ_4e.Clone())
    histos.append(h_ggZZ_4mu.Clone())
    histos.append(h_ggZZ_4tau.Clone())
    histos.append(h_ggZZ_2e2mu.Clone())
    histos.append(h_ggZZ_2e2tau.Clone())
    histos.append(h_ggZZ_2mu2tau.Clone())
    histos.append(h_TTZ.Clone())
    histos.append(h_TTG.Clone())
    histos.append(h_WWZ.Clone())
    histos.append(h_WZG.Clone())
    histos.append(h_WZZ.Clone())
    histos.append(h_ZZZ.Clone())
    histos.append(h_WZTo3L.Clone())
    histos.append(h_WZTo2L.Clone())
    histos.append(h_ZG.Clone())
    histos.append(h_DoubleMu.Clone()) 
    histos.append(h_SingleMu.Clone())
    histos.append(h_DoubleEle.Clone())
    histos.append(h_SingleEle.Clone())
    histos.append(h_MuonEle.Clone())

    for i in range(0,22):
      histos[i]=overunder_flowbin(histos[i])

    c1 = plot_ZZGregion.draw_plots(histos, 1, hists_name[ij])
    del histos[:]
 
if __name__ == "__main__":
  start = time.time()
  start1 = time.clock() 
  ZZG_Analysis()
  end = time.time()
  end1 = time.clock()
  print "wall time:", end-start
  print "process time:", end1-start1
