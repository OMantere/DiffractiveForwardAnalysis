//#include "DiffarctiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"
#include "Canvas.h"
#include "TH1.h"
#include "TLorentzVector.h"

#include <iostream>

double pair_mass,delta_phi, Wgg, Wmiss, delta_eta, mreco;
double Emiss, ptMiss, ptTot;
using namespace std;
/*
void draw_canvas( TH1D h, const char* histo_name, const char* histo_text) {
  Canvas c(histo_name, histo_text);
  h.Draw();
  c.Prettify(&h);
  c.Save( "pdf" );
}
*/

void lataaPuutTesti( 
        const char *LM6 = "computedLM6.root",
        const char *LM2 = "computedLM2.root"){
    TFile LM2f(LM2);
    TFile LM6f(LM6);


    TTree *tree1 = (TTree *) LM2f.Get("computed");
    TTree *tree2 = (TTree *) LM6f.Get("computed");


    tree1->SetBranchAddress("pair_dphi", &delta_phi);
    tree1->SetBranchAddress("Wgg",&Wgg);
    tree1->SetBranchAddress("mreco",&mreco);
    tree1->SetBranchAddress("pair_mass",&pair_mass);
    tree1->SetBranchAddress("Emiss", &Emiss);
    tree1->SetBranchAddress("ptMiss",&ptMiss);
    tree1->SetBranchAddress("ptTot",&ptTot);
    int pass1 = 0;
    
    int N1 = tree1->GetEntriesFast();
    

    tree2->SetBranchAddress("mreco",&mreco);
    tree2->SetBranchAddress("pair_mass",&pair_mass);
    tree2->SetBranchAddress("Emiss", &Emiss);
    tree2->SetBranchAddress("ptMiss",&ptMiss);
    tree2->SetBranchAddress("ptTot",&ptTot);
    int N2 = tree2-> GetEntriesFast();

    int N = 2;
    
    TH1D *hmasses[100000];

    for(int i = 0; i < N; i++){
        hmasses[i] = new TH1D("a","b",70,0,1000);
    }

    //TH1D hmass1("Mreco","mreco",70,0,400);
            
    for(int i=0;i<N1;i++){
        tree1 -> GetEntry(i);
       // cout << "EVENT  " << i << endl;
        hmasses[0]->Fill(pair_mass);
    }

//    hmasses[0] = hmass1;
    
    TCanvas *c = new TCanvas("Durr durr");
    TH1D hmass2("testi2","mreco",70,0,400);
    
    for(int i = 0; i < N2; i++){
        tree2 -> GetEntry(i);
//        cout << "LM6 Mass: " << mreco << endl;

        hmasses[1]->Fill(pair_mass);

    }

  //  hmasses[1] = hmass2;
    //Canvas c("testi","hurrdurr");
    //TCanvas c("Testikanvas");
    //c->Divide(2,2,0,0);
    c->cd(1);
    //hmass1.SetName("LM2");
    hmasses[0]->Draw();
    hmasses[0]->SetLineColor(kRed);
// c.Update();
   // hmass2.SetName("LM6");
    hmasses[1]->Draw("same");
   // hmass2.SetLineColor(kBlue);
    c->Update(); 
    c->Print("asdf.png"); 
    
    return;
}
