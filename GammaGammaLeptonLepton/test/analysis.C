#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "Canvas.h"
#include "TH1.h"
#include "TLorentzVector.h"

#include <iostream>

using namespace std;

double Wgg, Wmiss, delta_eta, deltaR, delta_phi, mreco, Wlep;

bool paper_cuts() {
    return Wgg > 200
           && Wmiss > 186
           && deltaR < 3.
           && delta_eta < 2.1;
};

bool my_cuts() {
    return 100 < mreco && mreco < 250;
}

void analysis(
        const char *signalfile = "computed.root",
        const char *wwfile = "computed_ww.root",
    const char *dyfile = "computed_dy.root") {

    THStack * bgh = new THStack("bgh"," stacked");
    THStack * sih = new THStack("sih"," stacked");

    TH1F *wwh = new TH1F("h_0","h_0",300,0,250);
    TH1F *dyh = new TH1F("h_1","h_1",300,0,250);
    TH1F *lm1mu = new TH1F("h_2","h_2",300,0,250);

    double lm1mur_crossx = 0.2396;
    double ww_crossx = 0.7358;
    double dy_crossx = 0.7358;

    TFile file1(signalfile);
    TFile file2(wwfile);
    TFile file3(dyfile);
    TTree *tree1 = (TTree *) file1.Get("computed");
    TTree *tree2 = (TTree *) file2.Get("computed");
    TTree *tree3 = (TTree *) file3.Get("computed");

//    tree1->SetBranchAddress("pair_dphi", &delta_phi);
//    tree1->SetBranchAddress("Wgg", &Wgg);
//    tree1->SetBranchAddress("Wmiss", &Wmiss);
//    tree1->SetBranchAddress("pair_eta_diff", &delta_eta);
//    tree1->SetBranchAddress("deltaR", &deltaR);
//    tree1->SetBranchAddress("mreco", &mreco);
    tree1->SetBranchAddress("Wlep", &Wlep);
    int pass1 = 0;
    int pass2 = 0;

    int N1 = tree1->GetEntriesFast();
    for (int i = 0; i < N1; i++) {
        tree1->GetEntry(i);
        lm1mu->Fill(Wlep);
        if(paper_cuts() && my_cuts()) pass1++;
    }

//    tree2->SetBranchAddress("pair_dphi", &delta_phi);
//    tree2->SetBranchAddress("Wgg", &Wgg);
//    tree2->SetBranchAddress("Wmiss", &Wmiss);
//    tree2->SetBranchAddress("pair_eta_diff", &delta_eta);
//    tree2->SetBranchAddress("deltaR", &deltaR);
//    tree2->SetBranchAddress("mreco", &mreco);
    tree2->SetBranchAddress("Wlep", &Wlep);
    int N2 = tree2->GetEntriesFast();
    for (int i = 0; i < N2; i++) {
        tree2->GetEntry(i);
        wwh->Fill(Wlep);
        if(paper_cuts() && my_cuts()) pass2++;
    }

//    tree3->SetBranchAddress("pair_dphi", &delta_phi);
//    tree3->SetBranchAddress("Wgg", &Wgg);
//    tree3->SetBranchAddress("Wmiss", &Wmiss);
//    tree3->SetBranchAddress("pair_eta_diff", &delta_eta);
//    tree3->SetBranchAddress("deltaR", &deltaR);
//    tree3->SetBranchAddress("mreco", &mreco);
    tree3->SetBranchAddress("Wlep", &Wlep);
    for (int i = 0; i < tree3->GetEntriesFast(); i++) {
        tree3->GetEntry(i);
        dyh->Fill(Wlep);
        if(paper_cuts() && my_cuts()) pass2++;
    }

    double signal_crossx = 0.2396;
    double bg_crossx = 0.7358;
    double crossx_ratio = signal_crossx/bg_crossx;
    double detector_acceptance = (double)N1/(double)N2;
    double analysis_ratio = (double)pass1/(double)pass2;

    cout << "Signal events passing analysis cuts: " << pass1 << endl;
    cout << "BG events passing analysis cuts: " << pass2 << endl;
    cout << "Detector ratio: " << detector_acceptance << endl;
    cout << "Analysis improvement: " << analysis_ratio << endl;
    cout << "S/B ratio: " << crossx_ratio * detector_acceptance * analysis_ratio << endl;
    cout << "Post-analysis Signal cross-section: " << signal_crossx * (double)pass1/(double)1000 << endl;

    // Plot stacked
    lm1mu->SetFillColor(kRed);
    sih->Add(lm1mu);
    wwh->SetFillColor(kBlue);
    bgh->Add(wwh);
    dyh->SetFillColor(kGreen);
    bgh->Add(dyh);
    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->AddEntry("lm1mu","GG -> muR+muR- (LM1)","l");
    legend->AddEntry("ww","WW -> muR+muR-","l");
    legend->AddEntry("dy","Drell-Yan -> muR+muR-","l");

    TCanvas * c1 = new TCanvas("c1","stacked hists",700,900);
    bgh->Draw();
    sih->Draw("same");
    legend->Draw();
    c1->Update();


//    With the cuts in the paper
//    Signal events passing analysis cuts: 489
//    BG events passing analysis cuts: 79
//    Detector ratio: 1.02421
//    Analysis improvement: 6.18987
//    S/B ratio: 2.06442
//    Post-analysis Signal cross-section: 0.117164
}