#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "Canvas.h"
#include "TH1.h"
#include "TLorentzVector.h"

#include <iostream>

#include <time.h>

clock_t start = clock(), diff;

void start_time() {
    start = clock();
}

void end_time() {
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds", msec/1000, msec%1000);
}

using namespace std;

double Wgg, Wmiss, delta_eta, deltaR, delta_phi, mreco, Wlep;

bool paper_cuts() {
    return Wgg > 200
           && Wmiss > 186
           && deltaR < 3.
           && delta_eta < 2.1;
};

bool my_cuts() {
    return true;
}

void analysis(
        const char *justplot = "n",
        const char *signal_r = "computed_slr2.root",
        const char *signal_l = "computed_sll3.root",
        const char *wwfile = "computed_ww.root") {

    THStack * sh1 = new THStack("sh1"," Without cuts");
    THStack * sh2 = new THStack("sh2"," With analysis cuts");

//    TH1D *lm1e = new TH1D("h_2","h_1",10,0,250);
//    TH1D *wwh2 = new TH1D("h_02","h_02",10,0,250);
//    TH1D *lm1mu2 = new TH1D("h_12","h_12",10,0,250);
//    TH1D *lm1e2 = new TH1D("h_22","h_12",10,0,250);
//    TH1D *dyh = new TH1D("h_2","h_1",100,0,250);

    double lm1mur_crossx = 0.263;
    double lm1mul_crossx = 0.055;
    double ww_crossx = 0.7583;
//    double lm1er_crossx = 0.3324;
//    double dy_crossx = 0;

    TH1D *wwh = new TH1D("h_ww","h_0",20,0,250);
    TH1D *wwh2 = new TH1D("h_ww2","h_0",20,0,250);
    TH1D *lm1mur = new TH1D("h_mur","h_1",20,0,250);
    TH1D *lm1mul = new TH1D("h_mul","h_1",20,0,250);
    TH1D *lm1mur2 = new TH1D("h_mur2","h_1",20,0,250);
    TH1D *lm1mul2 = new TH1D("h_mul2","h_1",20,0,250);

    TFile file1(signal_r);
    TFile file2(wwfile);
    TFile file3(signal_l);
    TTree *tree1 = (TTree *) file1.Get("computed");
    TTree *tree2 = (TTree *) file2.Get("computed");
    TTree *tree3 = (TTree *) file3.Get("computed");

    start_time();

    tree1->SetBranchAddress("pair_dphi", &delta_phi);
    tree1->SetBranchAddress("Wgg", &Wgg);
    tree1->SetBranchAddress("Wgg", &Wgg);
    tree1->SetBranchAddress("Wmiss", &Wmiss);
    tree1->SetBranchAddress("pair_eta_diff", &delta_eta);
    tree1->SetBranchAddress("deltaR", &deltaR);
    tree1->SetBranchAddress("mreco", &mreco);
    tree1->SetBranchAddress("Wlep", &Wlep);
    int pass1 = 0;
    int pass2 = 0;

    int N1 = tree1->GetEntriesFast();
    for (int i = 0; i < N1; i++) {
        tree1->GetEntry(i);
        lm1mur->Fill(Wlep);
        lm1mur2->Fill(Wlep);
        if(paper_cuts() && my_cuts()) {
            pass1++;
            lm1mur2->Fill(Wlep);
        }
    }

    end_time();

    tree2->SetBranchAddress("pair_dphi", &delta_phi);
    tree2->SetBranchAddress("Wgg", &Wgg);
    tree2->SetBranchAddress("Wmiss", &Wmiss);
    tree2->SetBranchAddress("pair_eta_diff", &delta_eta);
    tree2->SetBranchAddress("deltaR", &deltaR);
    tree2->SetBranchAddress("mreco", &mreco);
    tree2->SetBranchAddress("Wlep", &Wlep);
    int N2 = tree2->GetEntriesFast();
    for (int i = 0; i < N2; i++) {
        tree2->GetEntry(i);
        wwh->Fill(Wlep);
        if(paper_cuts() && my_cuts()) {
            pass2++;
            wwh2->Fill(Wlep);
        }
    }

    tree3->SetBranchAddress("pair_dphi", &delta_phi);
    tree3->SetBranchAddress("Wgg", &Wgg);
    tree3->SetBranchAddress("Wmiss", &Wmiss);
    tree3->SetBranchAddress("pair_eta_diff", &delta_eta);
    tree3->SetBranchAddress("deltaR", &deltaR);
    tree3->SetBranchAddress("mreco", &mreco);
    tree3->SetBranchAddress("Wlep", &Wlep);
    int N3 = tree3->GetEntriesFast();
    for (int i = 0; i < N3; i++) {
        tree3->GetEntry(i);
        lm1mul->Fill(Wlep);
        if(paper_cuts() && my_cuts()) {
            pass1++;
            lm1mul2->Fill(Wlep);
        }
    }

    double signal_crossx = lm1mur_crossx + lm1mul_crossx;
    double bg_crossx = ww_crossx;
    double crossx_ratio = signal_crossx/bg_crossx;
    double detector_acceptance = (double)N1/(double)N2;
    double analysis_ratio = (double)pass1/(double)pass2;

    cout << endl;
    cout << "Signal events passing analysis cuts: " << pass1 << endl;
    cout << "BG events passing analysis cuts: " << pass2 << endl;
    cout << "Detector ratio: " << detector_acceptance << endl;
    cout << "Analysis improvement: " << analysis_ratio << endl;
    cout << "S/B ratio: " << crossx_ratio * detector_acceptance * analysis_ratio << endl;
    cout << "Post-analysis signal cross-section: " << signal_crossx * (double)pass1/(double)2000 << endl;
    cout << endl;

    // Plot stacked
    wwh->Scale(bg_crossx/N2);
//    wwh->SetFillColor(kBlue);
    wwh->SetLineColor(kRed);
    wwh->SetLineWidth(2);
    lm1mur->Scale(lm1mur_crossx/N1);
    lm1mur->SetFillColor(kBlue);
    lm1mur->SetLineColor(kBlack);
    lm1mur->SetLineWidth(1);
    sh1->Add(lm1mur);
    lm1mul->Scale(lm1mul_crossx/N1);
    lm1mul->SetFillColor(kGreen);
    lm1mul->SetLineColor(kBlack);
    lm1mul->SetLineWidth(1);
    sh1->Add(lm1mul);
    wwh2->Scale(bg_crossx/N2);
//    wwh2->SetFillColor(kBlue);
    wwh2->SetLineColor(kRed);
    wwh2->SetLineWidth(2);
    lm1mur2->Scale(lm1mur_crossx/N1);
    lm1mur2->SetFillColor(kBlue);
    lm1mur2->SetLineColor(kBlack);
    lm1mur2->SetLineWidth(1);
    sh2->Add(lm1mur2);
    lm1mul2->Scale(lm1mul_crossx/N1);
    lm1mul2->SetFillColor(kGreen);
    lm1mul2->SetLineColor(kBlack);
    lm1mul2->SetLineWidth(1);
    sh2->Add(lm1mul2);

    TCanvas * c1 = new TCanvas("c1","stacked hists",700,900);
    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->AddEntry(wwh,"WW -> l+l-","fl");
    legend->AddEntry(lm1mur,"GG -> muR+muR- -> l+l- (LM1)","f");
    legend->AddEntry(lm1mul,"GG -> muL+muL- -> l+l- (LM1)","f");
//    legend->AddEntry("dy","Drell-Yan -> muR+muR-","l");
//    sh2->Draw();
    c1->Divide(1, 3);
    c1->cd(1);
    sh1->Draw("HIST");
    wwh->Draw("HIST SAME");
    c1->cd(2);
    sh2->Draw("HIST");
    wwh2->Draw("HIST SAME");
    c1->cd(3);
    legend->Draw();
    c1->Update();
    sh1->GetXaxis()->SetTitle("Wlep [GeV]");
    sh2->GetXaxis()->SetTitle("Wlep [GeV]");
    sh1->SetMaximum(0.08);
    c1->Modified();
//    c1->Print("analysis.pdf");


//    With the cuts in the paper
//    Signal events passing analysis cuts: 489
//    BG events passing analysis cuts: 79
//    Detector ratio: 1.02421
//    Analysis improvement: 6.18987
//    S/B ratio: 2.06442
//    Post-analysis Signal cross-section: 0.117164
}