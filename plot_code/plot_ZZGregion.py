import ROOT
import numpy as np

import CMSTDRStyle
CMSTDRStyle.setTDRStyle().cd()
import CMSstyle
from array import array

def set_axis(the_histo, coordinate, title, is_energy):

	if coordinate == 'x':
		axis = the_histo.GetXaxis()
	elif coordinate == 'y':
		axis = the_histo.GetYaxis()
	else:
		raise ValueError('x and y axis only')

	axis.SetLabelFont(42)
	axis.SetLabelOffset(0.015)
	axis.SetLabelSize(0.03)
	axis.SetNdivisions(505)
	axis.SetTitleFont(42)
	axis.SetTitleOffset(1.15)
	axis.SetTitleSize(0.04)
	if (coordinate == "y"):axis.SetTitleOffset(1.2)
	if is_energy:
		axis.SetTitle(title+' [GeV]')
	else:
		axis.SetTitle(title) 

def draw_plots(hist_array =[], draw_data=0, x_name=''):

	lumi=41480.
	ZZG = hist_array[0].Clone()
	ZZG.SetFillColor(ROOT.kRed)
	ZZG.Scale(lumi)

	ZG = hist_array[16].Clone()
	ZG.SetFillColor(ROOT.kYellow-4)
	ZG.Scale(lumi)

	Top = hist_array[8].Clone()
	Top.SetFillColor(ROOT.kGray)
	Top.Add(hist_array[9])
	Top.Scale(lumi)

	VVV = hist_array[10].Clone()
	VVV.SetFillColor(ROOT.kSpring-9)
	VVV.Add(hist_array[11])
	VVV.Add(hist_array[12])
	VVV.Add(hist_array[13])
	VVV.Scale(lumi)

	WZ = hist_array[14].Clone()
	WZ.SetFillColor(ROOT.kCyan-7)
	WZ.Add(hist_array[15])
	WZ.Scale(lumi)

	ggZZ = hist_array[2].Clone()
	ggZZ.SetFillColor(ROOT.kBlue-7)
	ggZZ.Add(hist_array[3])
	ggZZ.Add(hist_array[4])
	ggZZ.Add(hist_array[5])
	ggZZ.Add(hist_array[6])
	ggZZ.Add(hist_array[7])
	ggZZ.Scale(lumi)

	ZZ = hist_array[1].Clone()
	ZZ.SetFillColor(ROOT.kViolet-4)
	ZZ.Scale(lumi)

	Data = hist_array[17].Clone()
	Data.Add(hist_array[18])
	Data.Add(hist_array[19])
	Data.Add(hist_array[20])
	Data.Add(hist_array[21])
	if not draw_data: Data.Reset('ICE')
	Data.SetMarkerStyle(20)
	Data.SetMarkerSize(0.85)
	Data.SetMarkerColor(1)
	Data.SetLineWidth(1)

	h_stack = ROOT.THStack()
	h_stack.Add(ZG)
	h_stack.Add(Top)
	h_stack.Add(VVV)
	h_stack.Add(WZ)
	h_stack.Add(ggZZ)
	h_stack.Add(ZZ)
	h_stack.Add(ZZG)
	max_yields = 0
	Nbins=h_stack.GetStack().Last().GetNbinsX()
	for i in range(1,Nbins+1):
		max_yields_temp = h_stack.GetStack().Last().GetBinContent(i)
		if max_yields_temp>max_yields:max_yields=max_yields_temp

	max_yields_data = 0
	for i in range(1,Nbins+1):
		max_yields_data_temp = Data.GetBinContent(i)
		if max_yields_data_temp>max_yields_data:max_yields_data=max_yields_data_temp

	h_stack.SetMaximum(max(max_yields, max_yields_data)*1.8)

	##MC error
	h_error = h_stack.GetStack().Last()
	h_error.SetBinErrorOption(ROOT.TH1.kPoisson);
	binsize = h_error.GetSize()-2;
	x = [];
	y = [];
	xerror_l = [];
	xerror_r = [];
	yerror_u = [];
	yerror_d = [];
	for i in range(0,binsize):
		x.append(h_error.GetBinCenter(i+1))
		y.append(h_error.GetBinContent(i+1))
		xerror_l.append(0.5*h_error.GetBinWidth(i+1))
		xerror_r.append(0.5*h_error.GetBinWidth(i+1))
		yerror_u.append(h_error.GetBinErrorUp(i+1))
		yerror_d.append(h_error.GetBinErrorLow(i+1))
	gr = ROOT.TGraphAsymmErrors(len(x), np.array(x), np.array(y),np.array(xerror_l),np.array(xerror_r), np.array(yerror_d), np.array(yerror_u))

	ZZG_yield =round(ZZG.Integral(),1)
	ZG_yield =round(ZG.Integral(),1)
	Top_yield =round(Top.Integral(),1)
	VVV_yield =round(VVV.Integral(),1)
	WZ_yield = round(WZ.Integral(),1)
	ggZZ_yield = round(ggZZ.Integral(),1)
	ZZ_yield = round(ZZ.Integral(),1)
	Data_yield = round(Data.Integral())

	c = ROOT.TCanvas()
	pad = ROOT.TPad()
	pad.Draw()
	h_stack.Draw('HIST')
	Data.Draw("SAME pe")

	gr.SetFillColor(1)
	gr.SetFillStyle(3005)
	gr.Draw("SAME 2")
	if 'lepspt' in x_name:set_axis(h_stack,'x', 'pt of four leptons', True)
	if 'lepseta' in x_name:set_axis(h_stack,'x', '#eta of four leptons', False)
	if 'lepsphi' in x_name:set_axis(h_stack,'x', 'phi of four leptons', False)
	if 'Z1_l1_pt' in x_name:set_axis(h_stack,'x', 'pt of Z1 leading lepton', True)
	if 'Z1_l1_eta' in x_name:set_axis(h_stack,'x', '#eta of Z1 leading lepton', False)
	if 'Z1_l1_phi' in x_name:set_axis(h_stack,'x', 'phi of Z1 leading lepton', False)
	if 'Z1_l2_pt' in x_name:set_axis(h_stack,'x', 'pt of Z1 subleading lepton', True)
	if 'Z1_l2_eta' in x_name:set_axis(h_stack,'x', '#eta of Z1 subleading lepton', False)
	if 'Z1_l2_phi' in x_name:set_axis(h_stack,'x', 'phi of Z1 subleading lepton', False)
	if 'Z2_l1_pt' in x_name:set_axis(h_stack,'x', 'pt of Z2 leading lepton', True)
	if 'Z2_l1_eta' in x_name:set_axis(h_stack,'x', '#eta of Z2 leading lepton', False)
	if 'Z2_l1_phi' in x_name:set_axis(h_stack,'x', 'phi of Z2 leading lepton', False)
	if 'Z2_l2_pt' in x_name:set_axis(h_stack,'x', 'pt of Z2 subleading lepton', True)
	if 'Z2_l2_eta' in x_name:set_axis(h_stack,'x', '#eta of Z2 subleading lepton', False)
	if 'Z2_l2_phi' in x_name:set_axis(h_stack,'x', 'phi of Z2 subleading lepton', False)
	if 'photon_pt' in x_name:set_axis(h_stack,'x', 'photon pt', True)
	if 'photon_eta' in x_name:set_axis(h_stack,'x', 'photon eta', False)
	if 'photon_phi' in x_name:set_axis(h_stack,'x', 'photon phi', False)
	if 'Zmass' in x_name:set_axis(h_stack,'x', 'Two Z mass', True)
	if 'Zpt' in x_name:set_axis(h_stack,'x', 'Two Z pt', True)
	if 'Zeta' in x_name:set_axis(h_stack,'x', 'Two Z #eta', False)
	if 'Zphi' in x_name:set_axis(h_stack,'x', 'Two Z phi', False)
	if 'Z1_mass' in x_name:set_axis(h_stack,'x', 'Leading Z mass', True)
	if 'Z1_pt' in x_name:set_axis(h_stack,'x', 'Leading Z pt', True)
	if 'Z1_eta' in x_name:set_axis(h_stack,'x', 'Leading Z #eta', False)
	if 'Z1_phi' in x_name:set_axis(h_stack,'x', 'Leading Z phi', False)
	if 'Z2_mass' in x_name:set_axis(h_stack,'x', 'Subleading Z mass', True)
	if 'Z2_pt' in x_name:set_axis(h_stack,'x', 'Subleading Z pt', True)
	if 'Z2_eta' in x_name:set_axis(h_stack,'x', 'Subleading Z #eta', False)
	if 'Z2_phi' in x_name:set_axis(h_stack,'x', 'Subleading Z phi', False)
	if 'Z1Gmass' in x_name:set_axis(h_stack,'x', 'Z1G mass', True)
	if 'Z1Gpt' in x_name:set_axis(h_stack,'x', 'Z1G pt', True)
	if 'Z2Gmass' in x_name:set_axis(h_stack,'x', 'Z2G mass', True)
	if 'Z2Gpt' in x_name:set_axis(h_stack,'x', 'Z2G pt', True)
	if 'ZZGmass' in x_name:set_axis(h_stack,'x', 'ZZG mass', True)
	if 'ZZGpt' in x_name:set_axis(h_stack,'x', 'ZZG pt', True)
	if 'ZZGeta' in x_name:set_axis(h_stack,'x', 'ZZG #eta', False)
	if 'ZZmass' in x_name:set_axis(h_stack,'x', 'ZZ mass', True)
	if 'ZZpt' in x_name:set_axis(h_stack,'x', 'ZZ pt', True)
	if 'ZZeta' in x_name:set_axis(h_stack,'x', 'ZZ #eta', False)
	# WZ region
	if 'wl_pt' in x_name:set_axis(h_stack,'x', 'W lepton pt', True)
	if 'wl_eta' in x_name:set_axis(h_stack,'x', 'W lepton #eta', False)
	if 'zl1_pt' in x_name:set_axis(h_stack,'x', 'Z Leading lepton pt', True)
	if 'zl1_eta' in x_name:set_axis(h_stack,'x', 'Z Leading lepton #eta', False)
	if 'zl2_pt' in x_name:set_axis(h_stack,'x', 'Z Subleading lepton pt', True)
	if 'zl2_eta' in x_name:set_axis(h_stack,'x', 'Z Subleading lepton #eta', False)
	if 'met' in x_name:set_axis(h_stack,'x', 'MET', True)
	if 'zmass' in x_name:set_axis(h_stack,'x', 'Z mass', True)
	
	set_axis(h_stack,'y', 'Event/Bin', False)

	CMSstyle.SetStyle(pad)

	##legend
	leg1 = ROOT.TLegend(0.72, 0.78, 0.94, 0.88)
	leg2 = ROOT.TLegend(0.48, 0.78, 0.7, 0.88)
	leg3 = ROOT.TLegend(0.24, 0.78, 0.46, 0.88)
	leg1.SetMargin(0.4)
	leg2.SetMargin(0.4)
	leg3.SetMargin(0.4)
	leg1.AddEntry(ZG,'ZG ['+str(ZG_yield)+']','f')
	leg1.AddEntry(Top,'Top ['+str(Top_yield)+']','f')
	leg1.AddEntry(VVV,'VVV ['+str(VVV_yield)+']','f')
	leg2.AddEntry(WZ,'WZ ['+str(WZ_yield)+']','f')
	leg2.AddEntry(ggZZ,'ggZZ ['+str(ggZZ_yield)+']','f')
	leg2.AddEntry(ZZ,'ZZ ['+str(ZZ_yield)+']','f')
	leg3.AddEntry(ZZG,'ZZG ['+str(ZZG_yield)+']','f')
	leg3.AddEntry(gr,'Stat. unc','f')
	leg3.AddEntry(Data,'Data ['+str(Data_yield)+']','pe')
	leg1.SetFillColor(ROOT.kWhite)
	leg1.Draw('same')
	leg2.SetFillColor(ROOT.kWhite);
	leg2.Draw('same');
	leg3.SetFillColor(ROOT.kWhite);
	leg3.Draw('same');

	c.Update()
	c.SaveAs(x_name+'.pdf')
	c.SaveAs(x_name+'.png')
	return c
	pad.Close()
	del hist_array
