#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "Canvas.h"
#include "TH1.h"
#include "TLorentzVector.h"

#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <vector>


using namespace std;

std::ifstream infile("fileList.txt");

std::string line;

int Colorlist[] = {600, 632, 416, 880, 432, 860, 900};



// Function for loading the filelist
vector<string> fileloader(string fnam="fileList.txt"){

    
    std::string line;
    vector<string> filelist;
    while(getline(infile, line)){
        if(line[0]!='#'){

           filelist.push_back(line);
        }
    }
    return filelist;

}

char *str2char (string str){ //Helper function for turning strings into char *:s
    char * out = new char[str.size() + 1];
    std::copy(str.begin(),str.end(),out);
    out[str.size()] = '\0';

    return out;
}

void rescaleY(TH1D *h, double scale){
    if(h->GetMaximum()<scale){
        h->SetMaximum(scale*1.1);
    }
    return;
}
// Function for drawing the histograms on each iteration
TCanvas *draw_histos(TCanvas *c,
        TH1D *hmasses[], TH1D *hMEs[],
        TH1D *hptMisses[], TH1D *hptTots[], 
        TH1D *hWgg[], string foldstr,
        string data_name, int id, double xsec)
{
    
    Double_t Wgg, pair_mass, Emiss, ptMiss, ptTot;
    
    char * name = str2char(data_name);
    char * fold = str2char(foldstr);
    cout << name << endl;
    TFile File(fold);
    TTree *tree = (TTree *) File.Get("computed");
 
    //Set the desired parameters here

    tree->SetBranchAddress("pair_mass",&pair_mass);
    tree->SetBranchAddress("Emiss", &Emiss);
    tree->SetBranchAddress("ptMiss",&ptMiss);
    tree->SetBranchAddress("ptTot",&ptTot);
    tree->SetBranchAddress("Wgg",&Wgg);

    int N = tree->GetEntriesFast();


    for(int j=0;j<N;j++){
        // Fill the historgams
        tree->GetEntry(j);
        hmasses[id]->Fill(pair_mass);
        hMEs[id]->Fill(Emiss);
        hptMisses[id]->Fill(ptMiss);
        hptTots[id]->Fill(ptTot);
        hWgg[id]->Fill(Wgg);
    }
    // Select the color for current iteration
    int color = Colorlist[id];
    
    //Plot the values
//hmasses[id]->Sumw2();    

    double norm = xsec/N;
    cout<<norm<<endl;
    hmasses[id]->Scale(norm);
    hMEs[id]->Scale(norm);
    hptMisses[id]->Scale(norm);
    hptTots[id]->Scale(norm);
    hWgg[id]->Scale(norm);
    //hmasses[id]->SetMaximum(1);

    rescaleY(hmasses[0], hmasses[id]->GetMaximum());
    rescaleY(hMEs[0], hMEs[id]->GetMaximum());
    rescaleY(hptMisses[0], hptMisses[id]->GetMaximum());
    rescaleY(hptTots[0], hptTots[id]->GetMaximum());
    rescaleY(hWgg[0],hWgg[id]->GetMaximum());

    //gStyle->SetErrorY(0.0001);
    c->cd(1);

    hmasses[id]->Draw("HIST SAME");
    hmasses[id]->SetLineColor(color);
    hmasses[id]->SetXTitle("[GeV]");
    c->cd(2);
    
    hMEs[id]->SetName(name);
    hMEs[id]->SetLineColor(color);
    hMEs[id]->Draw("HIST SAME");
    hMEs[id]->SetXTitle("[GeV]");
  //  hMEs[id]->SetFillColor(color);
    c->cd(3);

    hptMisses[id]->SetName(name);
    hptMisses[id]->SetLineColor(color);
    hptMisses[id]->Draw("HIST SAME");
    hptMisses[id]->SetXTitle("[GeV]");
//    hptMisses[id]->SetFillColor(color);
    c->cd(4);

    hptTots[id]->SetName(name);
    hptTots[id]->SetLineColor(color);
 //   hptTots[id]->SetFillColor(color);
    hptTots[id]->Draw("HIST SAME");
    hptTots[id]->SetXTitle("[GeV]");

    c->cd(5);

    hWgg[id]->SetName(name);
    hWgg[id]->SetLineColor(color);
    hWgg[id]->Draw("HIST SAME");
    hWgg[id]->SetXTitle("[GeV]");

    delete[] name;
    delete[] fold;
    
    c->Update();
    
    return c;
}

void combinedHisto(){

    vector<string> flist = fileloader();
 
    // Create the canvas and initialize a bunch of stuff
    
    TCanvas *c = new TCanvas("Testi");
    c->Divide(3,2);


    int N = flist.size();
    TH1D *hmasses[N];
    TH1D *hMEs[N];
    TH1D *hptMisses[N];
    TH1D *hptTots[N];
    TH1D *hWgg[N];

    vector<string> names;
    vector<double> xsecs;
    vector<string> outputs;
    auto* legend = new TLegend(0.2,0.2,0.8,0.8);
    // The iteration loop
    // 1. Obtain the parameters from the filelist
    // 2. Initialize the histograms
    // 3. Run the plotting function

    for(int i = 0; i < N; i++){
       vector<string> words;
        
       istringstream iss(flist[i]);
       copy(istream_iterator<string>(iss),
               istream_iterator<string>(),
               back_inserter(words));
       string name = words[0];
       names.push_back(name);

       double xsec;
       std::string::size_type sz;
       xsec = std::stof(words[1], &sz);
       xsecs.push_back(xsec);

       string output = words[2];
       outputs.push_back(output);

       hmasses[i] = new TH1D("Pair mass","Pair mass",40,0,1000);       
       hMEs[i] = new TH1D("Missing Energy","Missing Energy",10, 0, 1000);
       hptMisses[i] = new TH1D("Missing PT","Missing PT",10,0,400);
       hptTots[i] = new TH1D("PT protons","PT of protons",10,0,10);
       hWgg[i] = new TH1D("Wgg", "Wgg",30,0,1500);
       c = draw_histos(c,hmasses,hMEs,hptMisses,hptTots,hWgg, output, name, i, xsec);
        
       char *name_ =str2char(name);
       legend->AddEntry(hmasses[i],name_);

       

    }
    c->cd(6); 
    legend->Draw();
    c->Print("c.pdf");
    //c->BuildLegend();

    return;

}

