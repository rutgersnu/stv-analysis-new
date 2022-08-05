// Authors: Serenity Engel, Andy Mastbaum

{

TFile* f = TFile::Open("/uboone/data/users/mastbaum/stv-analysis-ru3/rustv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root");
TTree* stv_tree = (TTree*) f->Get("stv_tree");

// Efficiency = (selected signal events) / (all signal events)

TCanvas* c2 = new TCanvas;

stv_tree->Draw("p3_mu.Mag()>>h(20,0,2)", "mc_is_cc0pi_signal", "goff");
TH1F* h = (TH1F*) gDirectory->Get("h");

stv_tree->Draw("p3_mu.Mag()>>h2(20,0,2)", "mc_is_cc0pi_signal&sel_CC0pi", "goff");
TH1F* h2 = (TH1F*) gDirectory->Get("h2");

TEfficiency* eff = new TEfficiency(*h2, *h);
eff->SetTitle(";Muon momentum p_{#mu} (GeV/c);Efficiency");
eff->Draw();
c2->SaveAs("cc0pi_eff.pdf");

// Purity = (selected signal events) / (all selected events)

TCanvas* c3 = new TCanvas;

stv_tree->Draw("p3_mu.Mag()>>h3(20,0,2)", "sel_CC0pi", "goff");
TH1F* h3 = (TH1F*) gDirectory->Get("h3");

stv_tree->Draw("p3_mu.Mag()>>h4(20,0,2)", "mc_is_cc0pi_signal&sel_CC0pi", "goff");
TH1F* h4 = (TH1F*) gDirectory->Get("h4");

TEfficiency* pur = new TEfficiency(*h4, *h3);
pur->SetTitle(";Muon momentum p_{#mu} (GeV/c);Purity");
pur->Draw();
c2->SaveAs("cc0pi_pur.pdf");

}

