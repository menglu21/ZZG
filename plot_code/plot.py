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

	lumi=41500.
#	DY = hist_array[0].Clone()
#	DY.SetFillColor(ROOT.kYellow-4)
#	DY.Scale(lumi)
	ZG = hist_array[6].Clone()
	ZG.SetFillColor(ROOT.kYellow-4)
	ZG.Scale(lumi)

	Top = hist_array[2].Clone()
	Top.SetFillColor(ROOT.kGray)
	Top.Add(hist_array[8])
	Top.Scale(lumi)

	VVV = hist_array[3].Clone()
	VVV.SetFillColor(ROOT.kSpring-9)
	VVV.Add(hist_array[4])
	VVV.Add(hist_array[9])
	VVV.Scale(lumi)

	WZG = hist_array[5].Clone()
	WZG.SetFillColor(ROOT.kCyan-7)
	WZG.Scale(lumi)

	ggZZ = hist_array[7].Clone()
	ggZZ.SetFillColor(ROOT.kBlue-7)
	ggZZ.Scale(lumi)

	ZZ = hist_array[1].Clone()
	ZZ.SetFillColor(ROOT.kViolet-4)
	ZZ.Scale(lumi)

	ZZG = hist_array[0].Clone()
	ZZG.SetFillColor(ROOT.kRed-4)
	ZZG.Scale(lumi)

	Data = hist_array[10].Clone()
	Data.Add(hist_array[11])
	Data.Add(hist_array[12])
	Data.Add(hist_array[13])
	Data.Add(hist_array[14])
	if not draw_data: Data.Reset('ICE')
	Data.SetMarkerStyle(20)
	Data.SetMarkerSize(0.85)
	Data.SetMarkerColor(1)
	Data.SetLineWidth(1)

	h_stack = ROOT.THStack()
	h_stack.Add(ZG)
	h_stack.Add(Top)
	h_stack.Add(VVV)
	h_stack.Add(WZG)
	h_stack.Add(ggZZ)
	h_stack.Add(ZZ)
	h_stack.Add(ZZG)
	max_yields = 0
	Nbins=h_stack.GetStack().Last().GetNbinsX()
	for i in range(1,Nbins+1):
		max_yields_temp = h_stack.GetStack().Last().GetBinContent(i)
		if max_yields_temp>max_yields:max_yields=max_yields_temp
	h_stack.SetMaximum(max_yields*1.8)

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

	ZG_yield =round(ZG.Integral(),1)
	Top_yield =round(Top.Integral(),1)
	VVV_yield =round(VVV.Integral(),1)
	WZG_yield = round(WZG.Integral(),1)
	ggZZ_yield = round(ggZZ.Integral(),1)
	ZZ_yield = round(ZZ.Integral(),1)
	ZZG_yield = round(ZZG.Integral(),1)
	Data_yield = round(Data.Integral())

	c = ROOT.TCanvas()
	pad = ROOT.TPad()
	pad.Draw()
	h_stack.Draw('HIST')
	Data.Draw("SAME pe")

	gr.SetFillColor(1)
	gr.SetFillStyle(3005)
	gr.Draw("SAME 2")
	if 'Z1l1pt' in x_name:set_axis(h_stack,'x', 'pt of Z1 leading lepton', True)
	if 'Z1l1eta' in x_name:set_axis(h_stack,'x', '#eta of Z1 leading lepton', False)
	if 'Z1l1phi' in x_name:set_axis(h_stack,'x', 'phi of Z1 leading lepton', False)
	if 'Z1l2pt' in x_name:set_axis(h_stack,'x', 'pt of Z1 subleading lepton', True)
	if 'Z1l2eta' in x_name:set_axis(h_stack,'x', '#eta of Z1 subleading lepton', False)
	if 'Z1l2phi' in x_name:set_axis(h_stack,'x', 'phi of Z1 subleading lepton', False)
	if 'Z2l1pt' in x_name:set_axis(h_stack,'x', 'pt of Z2 leading lepton', True)
	if 'Z2l1eta' in x_name:set_axis(h_stack,'x', '#eta of Z2 leading lepton', False)
	if 'Z2l1phi' in x_name:set_axis(h_stack,'x', 'phi of Z2 leading lepton', False)
	if 'Z2l2pt' in x_name:set_axis(h_stack,'x', 'pt of Z2 subleading lepton', True)
	if 'Z2l2eta' in x_name:set_axis(h_stack,'x', '#eta of Z2 subleading lepton', False)
	if 'Z2l2phi' in x_name:set_axis(h_stack,'x', 'phi of Z2 subleading lepton', False)
	if 'photonpt' in x_name:set_axis(h_stack,'x', 'photon pt', True)
	if 'photoneta' in x_name:set_axis(h_stack,'x', 'photon eta', False)
	if 'photonphi' in x_name:set_axis(h_stack,'x', 'photon phi', False)
	if 'Z1mass' in x_name:set_axis(h_stack,'x', 'Leading Z mass', True)
	if 'Z1pt' in x_name:set_axis(h_stack,'x', 'Leading Z pt', True)
	if 'Z1eta' in x_name:set_axis(h_stack,'x', 'Leading Z #eta', False)
	if 'Z1phi' in x_name:set_axis(h_stack,'x', 'Leading Z phi', False)
	if 'Z2mass' in x_name:set_axis(h_stack,'x', 'Subleading Z mass', True)
	if 'Z2pt' in x_name:set_axis(h_stack,'x', 'Subleading Z pt', True)
	if 'Z2eta' in x_name:set_axis(h_stack,'x', 'Subleading Z #eta', False)
	if 'Z2phi' in x_name:set_axis(h_stack,'x', 'Subleading Z phi', False)
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
	leg2.AddEntry(WZG,'WZG ['+str(WZG_yield)+']','f')
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
