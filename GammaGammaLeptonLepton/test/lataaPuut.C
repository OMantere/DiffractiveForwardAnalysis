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

void lataaPuut( 
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

    TH1D hmass1("Mreco","mreco",70,0,400);
    TH1D hME1("ME","Emiss",70,0,1000);
    TH1D hptMiss1("ptMiss","ptMiss",70,0,300);
    TH1D hptTot1("ptTot","ptTot",10,0,10);
            
    for(int i=0;i<N1;i++){
        tree1 -> GetEntry(i);
       // cout << "EVENT  " << i << endl;
        cout << "Mass:  " << mreco << ", ME: " << Emiss << ", ptMiss: " << ptMiss << ", ptTot: " << ptTot<< endl;
        hmass1.Fill(pair_mass);
        hME1.Fill(Emiss);
        hptMiss1.Fill(ptMiss);
        hptTot1.Fill(ptTot);
    }

  
    TH1D hmass2("testi2","mreco",70,0,400);
    TH1D hME2("Me","Emiss",70,0,1000);
    TH1D hptMiss2("ptMiss","ptMiss",70,0,300);
    TH1D hptTot2("ptTot","ptTot",10,0,10);
    
    for(int i = 0; i < N2; i++){
        tree2 -> GetEntry(i);
     //   cout << "LM6 Mass: " << mreco << endl;

        hmass2.Fill(pair_mass);
        hME2.Fill(Emiss);
        hptMiss2.Fill(ptMiss);
        hptTot2.Fill(ptTot);

    }

    //Canvas c("testi","hurrdurr");
    TCanvas c("Testikanvas");
    c.Divide(2,2,0,0);
    c.cd(1);
    hmass1.SetName("LM2");
    hmass1.Draw();
    hmass1.SetLineColor(kRed);
    hmass1.SetLineColor(kRed);
   // c.Update();
    hmass2.SetName("LM6");
    hmass2.Draw("same");
    hmass2.SetMarkerColor(kBlue);
    hmass2.SetLineColor(kBlue);
    

    c.cd(2);
    

    hME1.SetName("LM2");
    hME1.Draw();
    hME1.SetLineColor(kRed);
    hME1.SetLineColor(kRed);
   // c.Update();
    hME2.SetName("LM6");
    hME2.Draw("same");
    hME2.SetMarkerColor(kBlue);
    hME2.SetLineColor(kBlue);

    c.cd(3);


    hptMiss1.SetName("LM2");
    hptMiss1.Draw();
    hptMiss1.SetLineColor(kRed);
    hptMiss1.SetLineColor(kRed);
    hptMiss2.SetName("LM6");
    hptMiss2.Draw("same");
    hptMiss2.SetMarkerColor(kBlue);
    hptMiss2.SetLineColor(kBlue);

    c.cd(4);
    hptTot1.SetName("LM2");
    hptTot1.Draw();
    hptTot1.SetLineColor(kRed);
    hptTot1.SetLineColor(kRed);
    hptTot2.SetName("LM6");
    hptTot2.Draw("same");
    hptTot2.SetMarkerColor(kBlue);
    hptTot2.SetLineColor(kBlue);

    c.Print("Testi.png");
    
}
