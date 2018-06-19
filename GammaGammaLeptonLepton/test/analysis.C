#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "Canvas.h"
#include "TH1.h"
#include "TLorentzVector.h"

#include <iostream>

using namespace std;

double Wgg, Wmiss, delta_eta, deltaR, delta_phi;

bool paper_cuts() {
    return Wgg > 200
           && Wmiss > 186
           && deltaR < 3.
           && delta_eta < 2.1;
};

void analysis(
        const char *signalfile = "computed.root",
        const char *bgfile = "computed_ww.root") {
    TFile file1(signalfile);
    TFile file2(bgfile);
    TTree *tree1 = (TTree *) file1.Get("computed");
    TTree *tree2 = (TTree *) file2.Get("computed");

    tree1->SetBranchAddress("pair_dphi", &delta_phi);
    tree1->SetBranchAddress("Wgg", &Wgg);
    tree1->SetBranchAddress("Wmiss", &Wmiss);
    tree1->SetBranchAddress("pair_eta_diff", &delta_eta);
    tree1->SetBranchAddress("deltaR", &deltaR);
    int pass1 = 0;
    int pass2 = 0;

    int N1 = tree1->GetEntriesFast();
    for (int i = 0; i < N1; i++) {
        tree1->GetEntry(i);
        if(paper_cuts()) pass1++;
    }

    tree2->SetBranchAddress("pair_dphi", &delta_phi);
    tree2->SetBranchAddress("Wgg", &Wgg);
    tree2->SetBranchAddress("Wmiss", &Wmiss);
    tree2->SetBranchAddress("pair_eta_diff", &delta_eta);
    tree2->SetBranchAddress("deltaR", &deltaR);
    int N2 = tree2->GetEntriesFast();
    for (int i = 0; i < N2; i++) {
        tree2->GetEntry(i);
        if(paper_cuts()) pass2++;
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


//    With the cuts in the paper
//    Signal events passing analysis cuts: 489
//    BG events passing analysis cuts: 79
//    Detector ratio: 1.02421
//    Analysis improvement: 6.18987
//    S/B ratio: 2.06442
//    Post-analysis Signal cross-section: 0.117164
}