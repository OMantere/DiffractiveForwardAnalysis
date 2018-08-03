#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "Canvas.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "Parameter.h"

#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <limits>

using namespace std;

std::string line;

double d1,d2,d3,d4,d5;
int Colorlist[] = {600, 632, 416, 880, 432, 860, 900};
std::string datafold = "/afs/cern.ch/work/k/karjas/private/dataFold/";
double plm, pslm, rapiditysl, rapidityl1l2;

bool plot1 = true;
bool plot2 = false;
bool plot3 = false;
// SELECT PARAMETERS HERE!

// Currently supported options
// Wgg: Invariant mass of the photons
// mass: Invariant mass of the lepton pair
// ptTot: Sum of the transverse momenta of the protons
// ptMiss: pt of neutralinos (pt_photons - pt_leptons)
// ptl1: pt of lepton 1
// ptl2: pt of lepton 2
// MET_noProt: pt of leptons
vector<std::string> paramList = {"trackiso", "ecaliso","nPrimTracks","Wlep","ptl1", "ptl2", "deltaR","ptl1l2","ptTot","eta1", "eta2","xi11", "xi21", "invdetA","xi3", "ptMiss", "pair_mass","xim", "xip"};

string param1="nPrimTracks";
string param2="trackiso";

TTree *tree;

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!CUTS!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// 0 corresponds to 1st parameter in paramlist etc



// Currently supports only 5 histograms 
// Function for loading the filelist
vector<string> fileloader(){
    std::string infile = datafold+"filelist2.txt";
    
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

void setParamRange(Parameter *pList[], string name, double low, double high){
    for(int i = 0; i < paramList.size(); i++){
        if(pList[i]->pNam == name){
            //cout << name << endl;
            pList[i]->setRange(low,high);
            return;
        }
    }
    cout << "Parameter '" << name << "' not found"<<endl;
    return;
}

////////////////////////////////////
//SET RANGES HERE!!!!!
//ASDASDAS
///////////////////////////////////

void setRanges(Parameter *pList[]){
//    setParamRange(pList, "eta1", -0.6, 3);
//    setParamRange(pList, "ptl1", 10, 2000);
//    setParamRange(pList, "ptl2", 10, 2000);
//    setParamRange(pList, "xip", 0.03, 0.15);
//    setParamRange(pList, "xim", 0.03, 0.15);
//    setParamRange(pList, "ptl1l2", 40, 100);
//    setParamRange(pList, "deltaR",0, 3);
//    setParamRange(pList, "Wgg", 0,2000);
//    setParamRange(pList, "nPrimTracks", -1, 3);

//    setParamRange(pList, "trackiso", 0, 0.1);

//    setParamRange(pList, "invdetA", -1,1);
//    setParamRange(pList, "mtaus2", 2000, 22500);
    
//    setParamRange(pList, "xi21", -5e-15, 5e-15);
//    setParamRange(pList, "xi3", -100,100);
//    setParamRange(pList, "pair_mass", 0, 70);
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

bool cuts(Parameter *pList[]){
//cout << *vals[4] <<endl;
    
    bool pass = true;

    for(int i=0; i<paramList.size(); i++){
        double *val = pList[i]->pAdd;
        if(!pList[i]->cut(*val)){
            pass = false;
        }
    }

    return pass;


};
void updateHisto(TH1D **h,int id, char *name,double norm){
    int color = Colorlist[id];

    h[id]->Scale(norm*50);
    
    rescaleY(h[0], h[id]->GetMaximum()); 
    
    auto ax0 = h[0]->GetXaxis();
    auto axis = h[id]->GetXaxis();
    axis->SetLimits(ax0->GetXmin(),ax0->GetXmax());
    h[id]->SetName(name);
    h[id]->Draw("HIST SAME");
    //h[id]->SetFillColor(color);
    h[id]->SetLineColor(color);
   // h[id]->SetXTitle("[GeV]");
    return; 
}

///////////////////////////////////
//////////////////////////////////

bool mysort(std::pair<double,double> p1, std::pair<double,double> p2) {return p1.second<p2.second; }

void scanParameter(TCanvas *c, TH1D* h1[], TH1D* h2[],
        string foldstr1, string foldstr2, 
        string data_name1, string data_name2, 
        double Norm, TLegend *legend){


    char * name1 = str2char(data_name1);
    char * name2 = str2char(data_name2);
    
    char * fold1 = str2char(foldstr1);
    char * fold2 = str2char(foldstr2);
    //cout << name << endl;
    TFile File1(fold1);
    TTree *tree1 = (TTree *) File1.Get("computed");

    TFile File2(fold2);
    TTree *tree2 = (TTree *) File2.Get("computed");
    
    char* Param1 = str2char(param1);
    char* Paramcut = str2char(param2);
    int NParam = paramList.size();

    Parameter *pp1;
    Parameter *pp2;

    if(std::find(paramList.begin(), paramList.end(), Param1) != paramList.end()){
        tree1->SetBranchAddress(Param1, &d1);
        tree2->SetBranchAddress(Param1, &d1);
        pp1 = new Parameter(Param1, &d1);
    }
    if(std::find(paramList.begin(), paramList.end(), Paramcut) != paramList.end()){
        tree1->SetBranchAddress(Paramcut, &d2);
        tree2->SetBranchAddress(Paramcut, &d2);
        pp2 = new Parameter(Paramcut, &d2);
    }


    int N1 = tree1->GetEntriesFast();
    int N2 = tree2->GetEntriesFast();

    vector<std::pair<double, double>> points1;
    vector<std::pair<double, double>> points2;
    
//    cout << N1 << endl;

    for(int j = 0; j<N1; j++){
        tree1->GetEntry(j);
        std::pair<double, double> point1 = std::make_pair(d1,d2);
        points1.push_back(point1);
        cout << d2<< endl;
    }

    for(int j = 0; j<N2; j++){
        tree2->GetEntry(j);
        std::pair<double, double> point2 = std::make_pair(d1,d2);
        points2.push_back(point2);
    }
    
    vector<std::pair<double, double>> sortd1  = points1;
    std::sort (sortd1.begin(), sortd1.end(), mysort);
    
    vector<std::pair<double, double>> sortd2  = points2;
    std::sort (sortd2.begin(), sortd2.end(), mysort);

    //cout << sortd[0].second<<", " << sortd[N].second << endl;
    //cout << "min: " << mind2 << ", max: " << maxd2 << endl;

    int nSteps = 100;

    int step = N1/(nSteps + 1);
    int cutoff = 0;
    

    double SB[nSteps];
    double cutoffvals[nSteps];

    int idx = 0;
    for(int k = 0; k < nSteps ; k++){
        cutoff = k*step;
        double cutoffVal = sortd1[cutoff].second;
        cutoffvals[k] = cutoffVal;


        int cutoff2 = 0;
        
        while(sortd2[cutoff2].second < cutoffVal){
            cutoff2 += 1;
        }
    
        if(k%10 == 0 && idx < 5){
            
            for(int ii = cutoff; ii < N1 ; ii++){
                h1[idx]->Fill(sortd1[ii].first);
            }
            for(int jj = cutoff2; jj < N2; jj++){
                h2[idx]->Fill(sortd2[jj].first);
                //h3->Fill(sortd2[jj].first);
            }
            
            c->cd(1);
            updateHisto(h1, idx, name1, 1);


            c->cd(3);
            updateHisto(h2, idx, name2, 1);
            //h2[idx]->Draw("SAME");
            string entry = param2 + ">" + std::to_string(sortd1[cutoff].second);
            legend->AddEntry(h1[idx],str2char(entry),"l");
            idx += 1;
        }
        SB[k] = (1.0*(N1-cutoff)/(N2-cutoff2)*Norm);

        
        

    }
    
    TGraph* gr = new TGraph(nSteps, cutoffvals, SB);

    gr->SetTitle();
    gr->GetYaxis()->SetTitle("S/B");
    
    gr->GetXaxis()->SetTitle(str2char(param2+" cutoff"));

    c->cd(2);
    gr->Draw();


    return;
}
void draw_2D_histos(TCanvas *c, TH2D *H1,
        TH2D *H2, string foldstr,
        string data_name, int id, double xsec){

    int NParam = paramList.size();
    char * name = str2char(data_name);
    char * fold = str2char(foldstr);
    //cout << name << endl;
    TFile File(fold);
    TTree *tree = (TTree *) File.Get("computed");

    tree->SetBranchAddress("pair_mass", &plm);
    tree->SetBranchAddress("slMass", &pslm);
    tree->SetBranchAddress("rapiditysl", &rapiditysl);
    tree->SetBranchAddress("rapidityl1l2", &rapidityl1l2);
    int N = tree->GetEntriesFast();

    for(int j = 0; j<N; j++){
        tree->GetEntry(j);
        H1->Fill(plm, pslm);
        H2->Fill(rapiditysl, rapidityl1l2);
    }

    if(true){
        
        int color = Colorlist[id];
        //cout << color << endl;

        c->cd(id+1); 
        H1->Draw("COLZ"); 
        H1->SetTitle(str2char(data_name));
        H1->SetYTitle("slepton pair mass");
        H1->SetXTitle("lepton pair mass");

        c->cd(id+4);

        H2->Draw("COLZ");
        H2->SetTitle(str2char(data_name));
        H2->SetYTitle("slepton rapidity");
        H2->SetXTitle("lepton rapidity");
    }

    c->Update();
    return;
}
// Function for drawing the histograms on each iteration
/*
TCanvas *draw_histos(TCanvas *c,
        TH1D *h1[], TH1D *h2[],
        TH1D *h3[], TH1D *h4[], 
        TH1D *h5[], string foldstr,
        string data_name, int id, double xsec)
*/
TCanvas *draw_histos(TCanvas *c,
        TH1D **histolist[] , Parameter *pList[],  string foldstr,
        string data_name, int id, double norm)
{
    

    int NParam = paramList.size();
    char * name = str2char(data_name);
    char * fold = str2char(foldstr);
    //cout << name << endl;
    TFile File(fold);
    TTree *tree = (TTree *) File.Get("computed");

    
    double *p[NParam]; // Pointer array 
    
    double pp[NParam];
    for(int x=0;x<NParam;x++){
        p[x]=&pp[x];
    }
    
    for(int hID=0; hID<NParam; hID++){
        char *param = str2char(paramList[hID]);
        if(std::find(paramList.begin(), paramList.end(), param) != paramList.end()){
            tree->SetBranchAddress(param, p[hID]);
            pList[hID]->changeAddr(p[hID]); 
        }
    }
 

    int N = tree->GetEntriesFast();
    int pass = 0;
    cout << N << endl;
    for(int j=0;j<N;j++){
        // Fill the historgams
        tree->GetEntry(j);
        if(cuts(pList)){
            pass += 1;
            for(int k = 0; k < 5; k++){
                histolist[k][id]->Fill(*p[k]);
            
            }
        }
    }
 
    cout << "In " << data_name << " " << N - pass<<" out of " << N  << " (" << (N*1.-pass*1.)/N*100. << "%) events were cut" << endl;


    //double norm = 1; //xsec/10000; // N is incorrect, the correct value needs to be included in the filelist.txt (the number of generated events from madgraph)
    for(int j=0;j<5;j++){
        c->cd(j+1)->SetLogy();
        updateHisto(histolist[j],id,name,norm);
    }
    
    c->Update();
    return c;
}

void saveJSON(string runName = "",vector<string> datasets={""},
        Parameter *pList[]={}, vector<vector<int>> cutPoints={}){
    string out="{\n";
    
    out += "\"name\": \""+runName+"\",\n";
    out += "\"datasets\": [";
    int ii = 0;
    for(ii; ii<datasets.size();ii++){
        out += "\""+datasets[ii]+"\", ";
    }
    if(ii>0){
        out.pop_back();
        out.pop_back();
    }
    out += "],\n";
    
    out += "\"cuts\": [\n";
    for(ii=0; ii<paramList.size();ii++){

        out+="\t{ \"name\": \""+pList[ii]->pNam+"\", \"range\": [\""+std::to_string(pList[ii]->low)+"\",\""+std::to_string(pList[ii]->high)+"\"]},\n";    
    }

    if(ii>0){
        out.pop_back();
        out.pop_back();
    }
    out += "\n\t],\n";

    out += "\"pointsCut\": [\n";
    int jj = 0;
    for(ii=0; ii<datasets.size();ii++){
        for(jj = 0; jj < paramList.size();jj++){
            int points = cutPoints[ii][jj];
            if(points>0){
                out += "\t{ \"dataset\": \""+datasets[ii]+"\",\n";
                out += "\t  \"parameter\": \""+pList[jj]->pNam+"\",\n";
                out += "\t  \"cutPoints\": \""+std::to_string(points)+"\"},\n";
            }
        }
    }

    if(ii>0&&jj>0){
        out.pop_back();
        out.pop_back();
    }

    out += "\n\t]\n";

    out += "}\n";
    std::ofstream outfile;
    outfile.open("logfile.txt");
    //cout << out << endl;
    outfile << out;

    return;
}

//////
//MAIN FUNCTION
//////


void combinedHisto(){

    vector<string> flist = fileloader();
 
    // Create the canvas and initialize a bunch of stuff
    int N = flist.size();
    int NParam = paramList.size(); 


    TCanvas *c;
    TCanvas *c2;
    TCanvas *c3;
    
    TH1D *h1[N];
    TH1D *h2[N];
    TH1D *h3[N];
    TH1D *h4[N];
    TH1D *h5[N];
    TH1D **histolist[5]={h1,h2,h3,h4,h5};

    TH2D *H1[N];
    TH2D *H2[N];

    TH1D *hscans[2][5];
    TH2D *hscan2d[2];

    Parameter *pList[NParam];
    TLegend* legend;

    if(plot1){
        c = new TCanvas("Testi");
        legend = new TLegend (0.2, 0.2, 0.8, 0.8);
        c->Divide(3,2);
    }
    if(plot2){
        c2 = new TCanvas("2DHistos");
        c2->Divide(3,2);
    }
    if(plot3){
        c3 = new TCanvas("scanplots");
        c3->Divide(2,2);
        legend = new TLegend (0.2, 0.2, 0.8, 0.8);
    }


    vector<string> names;
    vector<double> xsecs;
    vector<double> errors;
    vector<string> outputs;
    vector<int> events;

    
    for(int hID=0; hID<NParam; hID++){
        char *param = str2char(paramList[hID]);
        if(std::find(paramList.begin(), paramList.end(), param) != paramList.end()){
            pList[hID] = new Parameter(param,NULL);
        }
    }
    setRanges(pList);
    

    vector<vector<int>> cutPoints;
    //auto* legend = new TLegend(0.2,0.2,0.8,0.8);
    // The iteration loop
    // 1. Obtain the parameters from the filelist
    // 2. Initialize the histograms
    // 3. Run the plotting function
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

       int eventcount = std::stoi(words[4], &sz);


       events.push_back(eventcount);
 
       double norm = xsec/eventcount * 40;
       cout << norm << endl;

       char* name_ = str2char(name);
       if(plot1){
           cout << "Durr" << endl;
           histolist[0][i] = new TH1D(paramList[0].c_str(),paramList[0].c_str(),50,0,0);
           histolist[1][i] = new TH1D(paramList[1].c_str(),paramList[1].c_str(),50,0,0);
           histolist[2][i] = new TH1D(paramList[2].c_str(),paramList[2].c_str(),50,0,0);
           histolist[3][i] = new TH1D(paramList[3].c_str(),paramList[3].c_str(),50,0,0);
           histolist[4][i] = new TH1D(paramList[4].c_str(),paramList[4].c_str(),50,0,0);
           cout << "Check" << endl;
           c = draw_histos(c,histolist, pList, out, name, i, norm);
           legend->AddEntry(h1[i],name_);
       }

       

       if(plot2){
           draw_2D_histos(c2, H1[i], H2[i], out, name, i, xsec);
           H1[i] = new TH2D("", "", 30, 0, 0, 30, 0, 0);
           H2[i] = new TH2D("", "", 30, 0, 0, 30, 0, 0);
       }
       if(plot3 && i < 2){

            hscans[i][0] = new TH1D("","",50,0,0);
            hscans[i][1] = new TH1D("","",50,0,0);
            hscans[i][2] = new TH1D("","",50,0,0);
            hscans[i][3] = new TH1D("","",50,0,0);
            hscans[i][4] = new TH1D("","",50,0,0);
       }

       //cout<<name_<<endl;


       vector<int> tempCuts;
       for(int ii=0; ii < NParam; ii++){
           tempCuts.push_back(pList[ii]->printCuts());
       }
       cutPoints.push_back(tempCuts);  

    }
    
    if(plot3){
        
        scanParameter(c3, hscans[0], hscans[1], outputs[0], outputs[1], names[0], names[1], (xsecs[0]/xsecs[1])*events[1]/events[0], legend);
        hscans[0][0]->SetTitle(str2char(param1));
        hscans[1][0]->SetTitle(str2char(param1));
        //legend->AddEntry(hscan2d[i], name_, "p");
        c3->cd(4);
        legend->Draw();
        c3->Update();    
    }
    
    if(plot1){
        c->cd(6); 
        legend->Draw();
        std::string foutName;
        for(int i = 0; i < paramList.size();i++){
            foutName += paramList[i]+"_";
        }
        foutName += ".pdf";
//c->Print(str2char(foutName));
    
        saveJSON("testi",names,pList, cutPoints);
    }
    return;

}


