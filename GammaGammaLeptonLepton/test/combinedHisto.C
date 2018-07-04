#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "Canvas.h"
#include "TH1.h"
#include "TLorentzVector.h"

#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <vector>
#include <algorithm>

using namespace std;

std::ifstream infile("fileList.txt");

std::string line;

int Colorlist[] = {600, 632, 416, 880, 432, 860, 900};
std::string datafold = "/afs/cern.ch/user/k/karjas/private/CMSSW/dataFold/";

double d1,d2,d3,d4,d5;

// SELECT PARAMETERS HERE!

// Currently supported options
// Wgg: Invariant mass of the photons
// mass: Invariant mass of the lepton pair
// ptTot: Sum of the transverse momenta of the protons
// ptMiss: pt of neutralinos (pt_photons - pt_leptons)
// ptl1: pt of lepton 1
// ptl2: pt of lepton 2
// MET_noProt: pt of leptons
vector<std::string> paramList = {"ptl1","ptl2","Wgg", "Emiss","MET_noProt"};
TTree *tree;

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!CUTS!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// 0 corresponds to 1st parameter in paramlist etc
bool cuts(double *vals[]){
//cout << *vals[4] <<endl;
//    return true;

        return *vals[0] < 60
        && *vals[1] < 60
        && *vals[2] > -1 
        && *vals[3] > 210

        && *vals[4] < 70;


};

// Currently supports only 5 histograms 
// Function for loading the filelist
vector<string> fileloader(){
    std::string infile = datafold+"filelist.txt";
    
    std::ifstream input(infile, std::ios::binary | ios::in);

    std::string line;
    vector<string> filelist;
    while(std::getline(input, line)){
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

void updateHisto(TH1D **h,int id, char *name,double norm){
    int color = Colorlist[id];

    h[id]->Scale(norm);
    
    rescaleY(h[0], h[id]->GetMaximum()); 
    h[id]->SetName(name);
    h[id]->Draw("HIST SAME");
    h[id]->SetLineColor(color);
    h[id]->SetXTitle("[GeV]");
    h[id]->Rebin(); 
    return; 
} 




// Function for drawing the histograms on each iteration
TCanvas *draw_histos(TCanvas *c,
        TH1D *h1[], TH1D *h2[],
        TH1D *h3[], TH1D *h4[], 
        TH1D *h5[], string foldstr,
        string data_name, int id, double xsec)
{
    
    int NParam = paramList.size();
    char * name = str2char(data_name);
    char * fold = str2char(foldstr);
    cout << name << endl;
    TFile File(fold);
    TTree *tree = (TTree *) File.Get("computed");

    double *p[5]={&d1,&d2,&d3,&d4,&d5}; // Pointer array

    TH1D **histolist[5]={h1,h2,h3,h4,h5}; // Histogram array
    
    for(int hID=0; hID<5; hID++){
        char *param = str2char(paramList[hID]);
        if(std::find(paramList.begin(), paramList.end(), param) != paramList.end()){
            cout << param << " found!" << endl;
            tree->SetBranchAddress(param, p[hID]);
        }
    }
   
    int N = tree->GetEntriesFast();


    int pass = 0;

    for(int j=0;j<N;j++){
        // Fill the historgams
        tree->GetEntry(j);
        if(cuts(p)){
            pass += 1;
            for(int k = 0; k < 5; k++){
                histolist[k][id]->Fill(*p[k]);
            
            }
        }
    }
    // Select the color for current iteration
    int color = Colorlist[id];
 
    cout << "In " << data_name << " " << N - pass<<" out of " << N  << " (" << (N*1.-pass*1.)/N*100. << "%) events were cut" << endl;

    double norm = xsec/N;
    for(int j=0;j<5;j++){
        c->cd(j+1);
        updateHisto(histolist[j],id,name,norm);
    }
    c->Update();
    return c;
}

void combinedHisto(){

    vector<string> flist = fileloader();
 
    // Create the canvas and initialize a bunch of stuff
    
    TCanvas *c = new TCanvas("Testi");
    c->Divide(3,2);

    int N = flist.size();
    int NParam = paramList.size(); 

    TH1D *h1[N];
    TH1D *h2[N];
    TH1D *h3[N];
    TH1D *h4[N];
    TH1D *h5[N];

    vector<string> names;
    vector<double> xsecs;
    vector<double> errors;
    vector<string> outputs;

    auto* legend = new TLegend(0.2,0.2,0.8,0.8);
    // The iteration loop
    // 1. Obtain the parameters from the filelist
    // 2. Initialize the histograms
    // 3. Run the plotting function
    cout<<N<<endl;
    for(int i = 0; i < N; i++){
       
       vector<string> words; // Iterator for the line
        
       istringstream iss(flist[i]);

       copy(istream_iterator<string>(iss), // Loads the data from the line to an vector
               istream_iterator<string>(),
               back_inserter(words));
       string name = words[0];
       names.push_back(name);

       double xsec;
       std::string::size_type sz;
       xsec = std::stof(words[1], &sz);
       xsecs.push_back(xsec);

       double error;
       error = std::stof(words[2], &sz);
       errors.push_back(error);

       string output = words[3]; 

       std::string out = datafold+"Computed/" + output;
       outputs.push_back(out);

       cout << out << endl; 
       
       char* name_ = str2char(name);
       h1[i] = new TH1D(paramList[0].c_str(),paramList[0].c_str(),40,0,0);
       h2[i] = new TH1D(paramList[1].c_str(),paramList[1].c_str(),40,0,0);
       h3[i] = new TH1D(paramList[2].c_str(),paramList[2].c_str(),40,0,0);
       h4[i] = new TH1D(paramList[3].c_str(),paramList[3].c_str(),40,0,0);
       h5[i] = new TH1D(paramList[4].c_str(),paramList[4].c_str(),40,0,0);

       cout<<"Check"<<endl;
       c = draw_histos(c,h1,h2,h3,h4,h5,  out, name, i, xsec);
     
       
       cout<<name_<<endl;

       legend->AddEntry(h1[i],name_);

       

    }
    c->cd(6); 
    legend->Draw();
    //c->BuildLegend();

    return;

}


