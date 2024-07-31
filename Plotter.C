#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "THStack.h"
#include "TLegend.h"
#include "TF1.h"
#include "TF2.h"
#include "TLine.h"
#include "TMath.h"
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

void ReadFill(TH1D *&h, std::string filename, int init_bin){
	std::ifstream file(filename);
	std::string line;
	int num_bins = 5;
	for (int i = 0; i < init_bin; i++){ std::getline(file, line); }
	for (int i = 1; i < num_bins+1; i++){
		std::getline(file, line);
		double cov = std::stod(line);
		h->SetBinContent(i, cov);
	}
	file.close();
}

double Mode(std::string filename, int init_bin, int num_univ){
	std::ifstream file(filename);
	std::cout << filename << std::endl;
	std::string line;
	std::vector<double> values;
	while (std::getline(file, line)){ values.push_back(std::stod(line)); }//read file to vector
	file.close();
	TH1D* htemp = new TH1D("htemp","htemp",100,0,1.2e6);
	for (int i = init_bin; i < num_univ; i+=7){
		htemp->Fill(values[i]);
	}
	values.clear();
	double mode = htemp->GetBinCenter(htemp->GetMaximumBin());
	std::cout << mode << std::endl;
	htemp->Reset();
	return mode;
}

void Plotter(){
	//Start canvas and histogram	
	TCanvas *c1 = new TCanvas("TruevRecoTrackComp","Simple Reco Muon E 1 Track",900,600);


	gROOT->SetBatch(false);
	gStyle->SetOptStat(1000111);
//	gStyle->SetOptStat(0);
	int num_universes = 100;
	int num_bins = 7;

	THStack* hs = new THStack("hs","Flux CV v Univ (mode)");

	TH1D* hflux = new TH1D("hflux","Flux weights",100,0,100);
	TH1D* hxsec = new TH1D("hxsec","XSec weights",100,0,100);
	TH1D* hcv   = new TH1D("hcv","CV weights",5,0,5);
	TH1D* hexps = new TH1D("hexps","hexps",5,0,5);
	TH1D* hhcur = new TH1D("hhcur","hexps",5,0,5);
	TH1D* hnuci = new TH1D("hnuci","hexps",5,0,5);
	TH1D* hnucq = new TH1D("hnucq","hexps",5,0,5);
	TH1D* hnuct = new TH1D("hnuct","hexps",5,0,5);
	TH1D* hpine = new TH1D("hpine","hexps",5,0,5);
	TH1D* hpqex = new TH1D("hpqex","hexps",5,0,5);
	TH1D* hptot = new TH1D("hptot","hexps",5,0,5);
	TH1D* hkmns = new TH1D("hkmns","hexps",5,0,5);
	TH1D* hkpls = new TH1D("hkpls","hexps",5,0,5);
	TH1D* hkzro = new TH1D("hkzro","hexps",5,0,5);
	TH1D* hpimi = new TH1D("hpimi","hexps",5,0,5);
	TH1D* hpipl = new TH1D("hpipl","hexps",5,0,5);

	TH1D* hall = new TH1D("hall","hall",5,0,5);

//	TH1D* huniv[num_universes];
//	for (int a = 0; a < num_universes; a++){ huniv[a] = new TH1D("a","a",5,0,5); }

	TH2D* hEvAngle = new TH2D("hEvAngle","Mock Reco E vs Azimuthal Angle CCInc",100,0,2000,100,0,90);

	//Fill histograms
/*	for (int j = 0; j < num_universes; j++){
		int bin = 7*j + 2;
		ReadFill(huniv[j],"univ/weight_All_UBGenie_univ.txt",bin);
		huniv[j]->SetLineColor(kRed);
		hs->Add(huniv[j]);
	}
*/
	for (int j = 2; j < num_bins; j++){
		hexps->SetBinContent(j-1, Mode("univ/weight_expskin_FluxUnisim_univ.txt", j, 1000));
		hhcur->SetBinContent(j-1, Mode("univ/weight_horncurrent_FluxUnisim_univ.txt",j,1000));
		hnuci->SetBinContent(j-1, Mode("univ/weight_nucleoninexsec_FluxUnisim_univ.txt", j, 1000));
		hnucq->SetBinContent(j-1, Mode("univ/weight_nucleonqexsec_FluxUnisim_univ.txt",j,1000));
		hnuct->SetBinContent(j-1, Mode("univ/weight_nucleontotxsec_FluxUnisim_univ.txt", j, 1000));
		hpine->SetBinContent(j-1, Mode("univ/weight_pioninexsec_FluxUnisim_univ.txt",j,1000));
		hpqex->SetBinContent(j-1, Mode("univ/weight_pionqexsec_FluxUnisim_univ.txt", j, 1000));
		hptot->SetBinContent(j-1, Mode("univ/weight_piontotxsec_FluxUnisim_univ.txt",j,1000));
		hkmns->SetBinContent(j-1, Mode("univ/weight_kminus_PrimaryHadronNormalization_univ.txt", j, 1000));
		hkpls->SetBinContent(j-1, Mode("univ/weight_kplus_PrimaryHadronFeynmanScaling_univ.txt",j,1000));
		hkzro->SetBinContent(j-1, Mode("univ/weight_kzero_PrimaryHadronSanfordWang_univ.txt", j, 1000));
		hpimi->SetBinContent(j-1, Mode("univ/weight_piminus_PrimaryHadronSWCentralSplineVariation_univ.txt",j,1000));
		hpipl->SetBinContent(j-1, Mode("univ/weight_piplus_PrimaryHadronSWCentralSplineVariation_univ.txt", j, 1000));
	}
	hexps->SetLineColor(kRed);
	hhcur->SetLineColor(kGreen);
	hnuci->SetLineColor(kYellow);
	hnucq->SetLineColor(kOrange);
	hnuct->SetLineColor(kBlack);
	hpine->SetLineColor(kCyan);
	hpqex->SetLineColor(kMagenta);
	hptot->SetLineColor(kGray);
	hkmns->SetLineColor(kSpring);
	hkpls->SetLineColor(kTeal);
	hkzro->SetLineColor(kAzure);
	hpimi->SetLineColor(kViolet);
	hpipl->SetLineColor(kPink);
	
	hs->Add(hexps);
	hs->Add(hhcur);
	hs->Add(hnuci);
	hs->Add(hnucq);
	hs->Add(hnuct);
	hs->Add(hpine);
	hs->Add(hpqex);
	hs->Add(hptot);
	hs->Add(hkmns);
	hs->Add(hkpls);
	hs->Add(hkzro);
	hs->Add(hpimi);
	hs->Add(hpipl);

/*	ReadFill(hcv,"cv/weight_All_UBGenie_cv.txt");
	ReadFill(hcv,"cv/weight_AxFFCCQEshape_UBGenie_cv.txt");
	ReadFill(hcv,"cv/weight_DecayAngMEC_UBGenie_cv.txt");
	ReadFill(hcv,"cv/weight_NormCCCOH_UBGenie_cv.txt");
	ReadFill(hcv,"cv/weight_Norm_NCCOH_UBGenie_cv.txt");
	ReadFill(hcv,"cv/weight_RPA_CCQE_UBGenie_cv.txt");
	ReadFill(hcv,"cv/weight_ThetaDelta2NRad_UBGenie_cv.txt");
	ReadFill(hcv,"cv/weight_Theta_Delta2Npi_UBGenie_cv.txt");
	ReadFill(hcv,"cv/weight_VecFFCCQEshape_UBGenie_cv.txt");
	ReadFill(hcv,"cv/weight_XSecShape_CCMEC_UBGenie_cv.txt");

	ReadFill(hcv,"cv/weight_expskin_FluxUnisim_cv.txt");
	ReadFill(hcv,"cv/weight_horncurrent_FluxUnisim_cv.txt");
	ReadFill(hcv,"cv/weight_nucleoninexsec_FluxUnisim_cv.txt");
	ReadFill(hcv,"cv/weight_nucleonqexsec_FluxUnisim_cv.txt");
	ReadFill(hcv,"cv/weight_nucleontotxsec_FluxUnisim_cv.txt");
	ReadFill(hcv,"cv/weight_pioninexsec_FluxUnisim_cv.txt");
	ReadFill(hcv,"cv/weight_pionqexsec_FluxUnisim_cv.txt");
	ReadFill(hcv,"cv/weight_piontotxsec_FluxUnisim_cv.txt");
	ReadFill(hcv,"cv/weight_kminus_PrimaryHadronNormalization_cv.txt");
	ReadFill(hcv,"cv/weight_kplus_PrimaryHadronFeynmanScaling_cv.txt");
	ReadFill(hcv,"cv/weight_kzero_PrimaryHadronSanfordWang_cv.txt");
	ReadFill(hcv,"cv/weight_piminus_PrimaryHadronSWCentralSplineVariation_cv.txt");
	ReadFill(hcv,"cv/weight_piplus_PrimaryHadronSWCentralSplineVariation_cv.txt");
*/
	ReadFill(hcv,"cv/weight_expskin_FluxUnisim_cv.txt",2);
	hcv->SetLineColor(kBlue);
//	hcv->SetLineWidth(6);

	hs->Add(hcv);

	hs->Draw("nostack");
	hs->GetXaxis()->SetTitle("bin index");
	hs->GetYaxis()->SetTitle("count");

//	gPad->SetLogx();
//	gPad->SetLogy();

	c1->Update();

	auto legend = new TLegend(0.6,0.4,0.9,0.9);
	legend->SetHeader("Legend","C");
	legend->AddEntry("hcv","CV","l");
	legend->AddEntry("hexps","Expskin","l");
	legend->AddEntry("hhcur","Horncurrent","l");
	legend->AddEntry("hnuci","Nuc Ine","l");
	legend->AddEntry("hnucq","Nuc QE","l");
	legend->AddEntry("hnuct","Nuc Tot","l");
	legend->AddEntry("hpine","Pion Ine","l");
	legend->AddEntry("hpqex","Pion QE","l");
	legend->AddEntry("hptot","Pion Tot","l");
	legend->AddEntry("hkmns","KMinus","l");
	legend->AddEntry("hkpls","KPlus","l");
	legend->AddEntry("hkzro","KZero","l");
	legend->AddEntry("hpimi","PiMinus","l");
	legend->AddEntry("hpipl","PiPlus","l");
	legend->Draw();

}
