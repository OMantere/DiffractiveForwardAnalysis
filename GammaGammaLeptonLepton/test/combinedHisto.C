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
int nPads = 0;
int ColorlistDef[] = {600, 632, 416, 880, 432, 860, 900};
vector<int> Colorlist;
std::string datafold = "/afs/cern.ch/work/k/karjas/private/dataFold/";
double plm, pslm, rapiditysl, rapidityl1l2;

bool plot1 = true;
bool plot2 = false;
bool plot3 = false;
bool plot4 = true;

bool logScale = true; 
int legendPad = 0;

// SELECT PARAMETERS HERE

vector<std::string> paramList = {"ptl1","ptl1l2", "pair_aco","mt2", "Wlep", "nPrimTracks","ptl2", "ecaliso","deltaEta", "pair_dphi","deltaR","ximProt","xipProt","Wgg", "Emiss", "trackiso",  "ptTot","eta1", "eta2","xi11", "xi21", "invdetA","xi3", "ptMiss", "pair_mass"};



vector<int> lowbounds = {0,0,0,0,0,0};
vector<double> upbounds = {360, 260, 1, 100, 400, 100};
    vector<string> units = {"[GeV]","[GeV]","" ,"[GeV]","[GeV]",""};


string param1="xipProt";
string param2="ximProt";

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

void setParamRange(Parameter *pList[], string name, double low, double high, bool incl){

    for(int i = 0; i < paramList.size(); i++){
        if(pList[i]->pNam == name){
            //cout << name << endl;
            pList[i]->setRange(low,high, incl);
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
//    setParamRange(pList, "trackiso", -1, 0.1);
    setParamRange(pList, "eta1", -2.4, 2.4,true);
    setParamRange(pList, "eta2", -2.4, 2.4, true);
    setParamRange(pList, "ptl2", 7, 1000, true);
    setParamRange(pList, "Wlep", 10,76, true);
    setParamRange(pList, "nPrimTracks",2, 2, true);
    setParamRange(pList, "Wgg", 236, 1000, true);
    setParamRange(pList, "pair_aco", 0.3, 1, true);
    setParamRange(pList, "ptl1l2", 17, 55, true);
    setParamRange(pList, "ptl1", 10, 1000, true);
    setParamRange(pList, "mt2", 10, 42, true);
    setParamRange(pList, "Emiss", 196, 4000, true);
//    setParamRange(pList, "trackiso", 0, 0.1);
    setParamRange(pList, "xipProt", 0.015, 0.15,true);
    setParamRange(pList, "ximProt", 0.015, 0.15,true);
//    setParamRange(pList, "deltaR",0, 3);
    
//    setParamRange(pList, "pair_dphi", -2,2);
    
}

char *str2char (string str){ //Helper function for turning strings into char *:s
    char * out = new char[str.size() + 1];
    std::copy(str.begin(),str.end(),out);
    out[str.size()] = '\0';

    return out;
}

void rescaleY(TH1D *h, TH1D *h2, int log){
 
    double scale = h2->GetMaximum();

    double ss = 1.1;
    if(log == 1) ss = 5;

    if(h->GetMaximum()<=scale){
        h->SetMaximum(scale*ss);
    }
    

    h2->SetMinimum(0.1);

    return;
}

bool cuts(Parameter *pList[], vector<int> idx){
   
    vector<int>  idx_;

    if(idx[0]==-1){
        for(int i = 0; i<paramList.size();i++){
            idx_.push_back(i);
        }
    } else idx_=idx;



    bool pass = true;

    for(int i=0; i < idx_.size(); i++){
        double *val = pList[idx_[i]]->pAdd;
        if(!pList[idx_[i]]->cut(*val)){
            pass = false;
        }
    }
    return pass;
};

TH1D *updateHisto(TH1D **h,int id, char *name,double norm, int log){
    int color = Colorlist[id];

    h[id]->Scale(norm*50);
    double scale = 1.1;

    rescaleY(h[0], h[id], log); 
    h[id]->SetStats(kFALSE);
    auto ax0 = h[0]->GetXaxis();
    auto axis = h[id]->GetXaxis();
    axis->SetLimits(ax0->GetXmin(),ax0->GetXmax());
    //h[id]->SetName(name);
    //h[id]->SetTitle(name);
   // h[id]->Draw("HIST SAME");
    //h[id]->SetFillColor(color);
    //h[id]->SetLineColor(1);
   // h[id]->SetXTitle("[GeV]");
    return h[id]; 
}
void setHistoBins(TH1D **h, int N, double min, double max){

    for(int id = 0; id < N; id++){
        auto axis = h[id]->GetXaxis();
        axis->SetLimits(min,max); 
    }
    return;
}

///////////////////////////////////
//////////////////////////////////

bool mysort(std::pair<double,double> p1, std::pair<double,double> p2) {return p1.second<p2.second; }

void scanParameter(TCanvas *c, TMultiGraph *gra, int id,
        string foldstr1, string foldstr2, 
        string data_name1, string data_name2, 
        double Norm, TLegend *legend){


    char * name1 = str2char(data_name1);
    char * name2 = str2char(data_name2);
    
    char * fold1 = str2char(foldstr1);
    char * fold2 = str2char(foldstr2);
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
        tree2->SetBranchAddress(Param1, &d2);
        pp1 = new Parameter(Param1, &d1);
    }
    if(std::find(paramList.begin(), paramList.end(), Paramcut) != paramList.end()){
        tree1->SetBranchAddress(Paramcut, &d3);
        tree2->SetBranchAddress(Paramcut, &d4);
        pp2 = new Parameter(Paramcut, &d2);
    }

    int N1 = tree1->GetEntriesFast();
    int N2 = tree2->GetEntriesFast();

    vector<std::pair<double, double>> points1;
    vector<std::pair<double, double>> points2;
    
    double minX1, minX2, maxX1, maxX2;

    minX1 = minX2 = INFINITY;
    maxX1 = maxX2 = -minX1;

    for(int j = 0; j<N1; j++){
        tree1->GetEntry(j);
        std::pair<double, double> point1 = std::make_pair(d1,d3);
        points1.push_back(point1);
        if(d1 < minX1) minX1 = d1;
        if(d1 > maxX1) maxX1 = d1;
    }

    for(int j2 = 0; j2<N2; j2++){
        tree2->GetEntry(j2);
        std::pair<double, double> point2 = std::make_pair(d2,d4);
        points2.push_back(point2);
        if(d2 < minX2) minX2 = d2;
        if(d2 > maxX2) maxX2 = d2;
    }
    
    vector<std::pair<double, double>> sortd1  = points1;
    std::sort (sortd1.begin(), sortd1.end(), mysort);
    
    vector<std::pair<double, double>> sortd2  = points2;
    std::sort (sortd2.begin(), sortd2.end(), mysort);
    cout << sortd1[0].first<<", "<<sortd1[N1-1].first<<endl;


    //setHistoBins(h1, 5, minX1, maxX1);
    //setHistoBins(h2, 5, minX2, maxX2);
    //setHistoBins(h1, 5, 0,0.1);
    //setHistoBins(h2, 5, 0,0.1);    
    
    int nSteps = 31;

    int step = N1/(nSteps + 1);
    int cutoff = 0;
    

    double SB[nSteps];
    double cutoffvals[nSteps];

    int idx = 0;
    for(int k = 0; k < nSteps ; k++){
        cout << k << endl;
        cutoff = k*step;
        //cout << cutoff << " Cutoff"  << endl;
        double cutoffVal = sortd1[cutoff].second;
        //cout << cutoffVal << endl;
        cutoffvals[k] = cutoffVal;


        int cutoff2 = 0;
        //cout << "Cutoffval" << sortd2[cutoff2] << endl;

        
        while(sortd2[cutoff2].second < cutoffVal){
            cutoff2 += 1;
        }
   /* 
        if(k%10 == 0 && idx < 5){
            
            for(int ii = cutoff; ii < N1 ; ii++){
                h1[idx]->Fill(sortd1[ii].first);
            }
            for(int jj = cutoff2; jj < N2; jj++){
                h2[idx]->Fill(sortd2[jj].first);
                //h3->Fill(sortd2[jj].first);
            }
            
            c->cd(1);
            h1[idx]->SetLineColor(ColorlistDef[idx]);
            updateHisto(h1, idx, name1, 1, c->cd(1)->GetLogy())->Draw("HIST SAME");
            h1[idx]->SetTitle(str2char(data_name1));


            c->cd(3);
            h2[idx]->SetLineColor(ColorlistDef[idx]);
            h2[idx] = updateHisto(h2, idx, name2, 1, c->cd(3)->GetLogy());
            h2[idx]->Draw("HIST SAME");
            h2[idx]->SetTitle(str2char(data_name2));
            string entry = param2 + ">" + std::to_string(sortd1[cutoff].second);
    //        legend->AddEntry(h1[idx],str2char(entry),"l");
            idx += 1;
        }*/
        cout<<(N1-cutoff)<<", "<<(N2-cutoff2)<<endl;
        SB[k] = (1.0/(N1-cutoff)*(N2-cutoff2)*Norm);

        
        

    }
    
    TGraph *gr = new TGraph(nSteps, cutoffvals, SB);
    //graphs[id] = gr;
    gr->SetTitle();

    gra->Add(gr);
    gr->SetMinimum(0);
    
    gr->GetXaxis()->SetTitle(str2char(param2+" cutoff"));
    gr->SetLineColor(ColorlistDef[id]);
    gr->GetYaxis()->SetTitle("S/B");
    legend->AddEntry(gr, str2char(data_name2),"l");

   // c->cd(1);


    gr->Draw("");
    c->Update();

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

    tree->SetBranchAddress("ptl1", &plm);
    tree->SetBranchAddress("ptl2", &pslm);
    tree->SetBranchAddress("nPrimTracks", &rapiditysl);
    tree->SetBranchAddress("trackiso", &rapidityl1l2);
    int N = tree->GetEntriesFast();

    for(int j = 0; j<N; j++){
        tree->GetEntry(j);
//H1->Fill(plm, pslm);
        H2->Fill(rapiditysl, rapidityl1l2);
    }

    if(true){
        
        int color = Colorlist[id];
        //cout << color << endl;
/*
        c->cd(id+1); 
        H1->Draw("COLZ"); 
        H1->SetTitle(str2char(data_name));
        H1->SetYTitle("slepton pair mass");
        H1->SetXTitle("lepton pair mass");
*/
        c->cd(id+1);

        H2->Draw("COLZ");
        H2->SetTitle(str2char(data_name));
        H2->SetYTitle("nPrimTracks");
        H2->SetXTitle("trackiso");
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

int countPads(TCanvas *c){

    int counter = 0;
    while(c->GetPad(counter)){
        counter += 1;
        c->cd(counter);
    }

    return counter - 1;

}

Canvas *draw_histos(Canvas *c,
        TH1D **histolist[] , Parameter *pList[],  string foldstr,
        string data_name, int id, double norm)
{
    

    int NParam = paramList.size();
    char * name = str2char(data_name);
    char * fold = str2char(foldstr);
    //cout << name << endl;
    TFile File(fold);
    TTree *tree = (TTree *) File.Get("computed");

    int nPads = countPads(c);

    if(nPads == 1) legendPad = 0;

    
    double *p[NParam]; // Pointer array 
    
    double pp[NParam];
    for(int x=0;x<NParam;x++){
        p[x]=&pp[x];
    }
    
    int N = tree->GetEntriesFast();

    for(int hID=0; hID<NParam; hID++){
        char *param = str2char(paramList[hID]);
        if(std::find(paramList.begin(), paramList.end(), param) != paramList.end()){
            tree->SetBranchAddress(param, p[hID]);
            pList[hID]->changeAddr(p[hID]); 
            pList[hID]->setnPoints(N);
        }
    }
 

    Parameter *pL1[2]={pList[0],pList[1]};
    Parameter *pL2[2]={pList[2],pList[3]};
    Parameter *pL3[2]={pList[4],pList[5]};

    vector<int> s1 = {0,3};
    vector<int> s2 = {0,1,3,4};
    vector<int> s3 = {0,1,2,3,4,5};
    vector<int> s4 = {-1};

    int pass1 = 0;
    int pass2 = 0;
    int pass3 = 0;

    int pass = 0;
    for(int j=0;j<N;j++){
        // Fill the historgams
        tree->GetEntry(j);
        /*if(cuts(pList)){
            pass += 1;
            
            for(int k = 0; k < nPads-legendPad; k++){
                if(id>0){
                    auto ax0 = histolist[0][id]->GetXaxis();
                    histolist[k][id]->GetXaxis()->SetLimits(ax0->GetXmin(), ax0->GetXmax());
                }
                histolist[k][id]->Fill(*p[k]);
            
            }
        }*/
        if(!plot4){
        if(cuts(pList, s4)){
            pass1 += 1;
            histolist[0][id]->Fill(*p[0]);
            histolist[3][id]->Fill(*p[3]);
        }
        if(cuts(pList, s4)){
            pass2 += 1;
            histolist[1][id]->Fill(*p[1]);
            histolist[4][id]->Fill(*p[4]);
        }
        if(cuts(pList,s4)){
            pass3+=1;
            histolist[2][id]->Fill(*p[2]);
            histolist[5][id]->Fill(*p[5]);
        }
    
        }else{
            if(cuts(pList,s4)){
                histolist[0][id]->Fill(*p[0]);
                pass3+=1;
            }
        }
    }
 
    cout << "In " << data_name << " " << N - pass3<<" out of " << N  << " (" << (N*1.-pass3*1.)/N*100. << "%) events were cut" << endl;


    for(int j=0;j<nPads - legendPad;j++){
        
        if(logScale)c->cd(j+1)->SetLogy();
        else c->cd(j+1);

        histolist[j][id]->GetYaxis()->SetTitle(str2char(""));
        
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
        if(pList[ii]->low!=-INFINITY){
        out+="\t{ \"name\": \""+pList[ii]->pNam+"\", \"range\": [\""+std::to_string(pList[ii]->low)+"\",\""+std::to_string(pList[ii]->high)+"\"]},\n";
        }
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
    outfile.open("logfile.txt",std::ios_base::app);
    outfile << out;

    return;
}

int setColor(string str){
    int out = 1;

    int L = str.size();
    std::size_t npos = str.npos;
    std::size_t foundDY = str.find("DY");
    if(foundDY!=npos){
        out = 800 -5; //800 = kOrange
    }
    std::size_t foundLM1 = str.find("LM1");
    if(foundLM1!=npos){
        out = 820-5; //820 = kSpring
    }
    std::size_t foundttbar = str.find("ttbar");
    if(foundttbar!=npos){
        out = 860-2; //860 = kAzure
    }
    std::size_t foundWW = str.find("WW");
    if(foundWW!=npos){
        out = 600 -2;
    }
    std::size_t founddiboson = str.find("WWjets");
    if(founddiboson!=npos){
        out = 600+3;
    }
    std::size_t foundelastic = str.find("elastic");
    if(foundelastic!=npos){
        out = 616-2;
    }
    std::size_t foundBG = str.find("BG");
    if(foundBG!=npos){
        out = 632-2;
    }
    std::size_t foundLM6 = str.find("LM6");
    if(foundLM6!=npos){
        out = 880-2;
    }
    std::size_t foundm35 = str.find("m35");
    if(foundm35!=npos){
        out = 920;
    }
    std::size_t foundm5 = str.find("m5");
    if(foundm5!=npos){
        out = 922;
    }
    return out;
}

//////
//MAIN FUNCTION
//////


void combinedHisto(){

    vector<string> flist = fileloader();
 
    // Create the canvas and initialize a bunch of stuff
    int N = flist.size();
    int NParam = paramList.size(); 

    Canvas *c;
    TCanvas *c2;
    TCanvas *c3;
    
    TH1D *h1[N];
    TH1D *h2[N];
    TH1D *h3[N];
    TH1D *h4[N];
    TH1D *h5[N];
    TH1D *h6[N];
    TH1D **histolist[6]={h1,h2,h3,h4,h5,h6};


    TH2D *H1[N];
    TH2D *H2[N];

    TH1D *hscans1[5];
    TH1D *hscans2[5];

    TMultiGraph *graphs = new TMultiGraph();

    Parameter *pList[NParam];
    TLegend* legend;

    if(plot1){
        c = new Canvas("Testi");
        if(!plot4) c->Divide(2,3);
        else c->Divide(1,1);
        c->SetLegendX1(0.62);
        nPads = countPads(c);
        cout << nPads << endl;
        //c->AddLegendEntry((TObject*)0, "L","");
        //
    }

    if(nPads == 1) legendPad = 0;

    if(plot2){
        c2 = new TCanvas("2DHistos");
        c2->Divide(2);
    }
    if(plot3){
        cout << "How about here?"<<endl;
        c3 = new TCanvas("scanplots");
       // c3->Divide(2,2);

        legend = new TLegend (0.75, 0.75, 0.9, 0.9);
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
    
    bool bg = true;

    int nPlots = nPads - legendPad;

    THStack *stacks[nPlots];
    TAxis *ax0;

    for(int j = 0; j < nPlots; j++){

        stacks[j] = new THStack(str2char(paramList[j]),"");
        stacks[j]->SetMinimum(0.1);

    }

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

       Colorlist.push_back(setColor(name));

       if(Colorlist[i]==815) bg = false;

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

       char* name_ = str2char(name);
       if(plot1){
           for(int j = 0; j < nPads - legendPad; j++){
               histolist[j][i] = new TH1D((paramList[j]+std::to_string(i)).c_str(),"",50,lowbounds[j],upbounds[j]);
          }

           c = draw_histos(c,histolist, pList, out, name, i, norm);
           for(int j = 0; j < nPads - legendPad; j++){
               
               if(bg){
                   histolist[j][i]->SetFillColor(Colorlist[i]);
                   stacks[j]->Add(updateHisto(histolist[j],i,name_,norm, c->cd(j+1)->GetLogy())); 
                   c->cd(j+1);
                   stacks[j]->Draw("HIST");
                   //cout<<stacks[j]->GetXaxis()<<endl;
                   stacks[j]->GetXaxis()->SetTitle(str2char(paramList[j]+" "+units[j]));
               }else{ 
                   histolist[j][i]->SetLineColor(Colorlist[i]);
                   histolist[j][i]->SetLineWidth(2);
                   c->Prettify(histolist[j][i]);
                   if(i>0){
                   //ax0[j] = stacks[j]->GetXaxis();
                   //cout<<ax0[j]->GetXmin()<<endl;
                   histolist[j][i]->GetXaxis()->SetLimits(stacks[j]->GetXaxis()->GetXmin(), stacks[j]->GetXaxis()->GetXmax());
                   }
                   updateHisto(histolist[j],i,name_,norm,c->cd(j+1)->GetLogy())->Draw("HIST SAME");
                   if(i==0){
                       histolist[j][i]->GetXaxis()->SetTitle(str2char(paramList[j]));
                   }

               }


           }
           //if(legendPad==1)c->AddLegendEntry(h1[i],name_, "f");
       }

       if(i>0 && plot3){
           cout << "Found this branch!" << endl;

           scanParameter(c3, graphs, i, outputs[0], outputs[i], names[0], names[i],1, legend); 
           c3->Update();

       }
       //c->Update();
/*
       if(legendPad==1){
           c->cd(6);
           TLegend *legend = c->GetLegend();
           legend->Draw();
       }
*/
       

       if(plot2){
           H1[i] = new TH2D("", "", 6, 0, 0, 5, 0, 4);
           H2[i] = new TH2D("", "", 30, 0, 100, 5, 0, 4);
           draw_2D_histos(c2, H1[i], H2[i], out, name, i, xsec);
       }



       vector<int> tempCuts;
       for(int ii=0; ii < NParam; ii++){
           tempCuts.push_back(pList[ii]->printCuts());
       }
       cutPoints.push_back(tempCuts);  

    }
    
    cout << "Got this far..?" << endl;
    if(plot3){
        graphs->Draw();
        legend->Draw();
/*
        hscans1[0] = new TH1D("","",50,0,0);
        hscans1[1] = new TH1D("","",50,0,0);
        hscans1[2] = new TH1D("","",50,0,0);
        hscans1[3] = new TH1D("","",50,0,0);
        hscans1[4] = new TH1D("","",50,0,0);
        
        hscans2[0] = new TH1D("","",50,0,0);
        hscans2[1] = new TH1D("","",50,0,0);
        hscans2[2] = new TH1D("","",50,0,0);
        hscans2[3] = new TH1D("","",50,0,0);
        hscans2[4] = new TH1D("","",50,0,0);
*/
        //cout << "So far so good" << endl;
        //scanParameter(c3, hscans1, hscans2, outputs[0], outputs[1], names[0], names[1], (xsecs[0]/xsecs[1])*events[1]/events[0]);
        return;
        //hscans1[0]->SetTitle(str2char(param1));
        //hscans2[0]->SetTitle(str2char(param1));
        //legend->AddEntry(hscan2d[i], name_, "p");
        //c3->cd(4);
        //c3->Update();    
    }

    TLine *l1, *l2;
    if(plot1 && plot4){
        l1 = new TLine(0.03,0,0.03,10e9);
        l1->SetLineColor(1);
        l1->SetLineWidth(2);
        l1->Draw();
        l2 = new TLine(0.15,0,0.15,10e9);
        l2->SetLineColor(1);
        l2->SetLineWidth(2);
        l2->Draw();
        c->cd(1)->RedrawAxis();
        c->Update();
    }

    
    /*if(!plot4){
        c->cd(6); 
        TLegend *legend = c->GetLegend();
        legend->Draw();
        std::string foutName;
        for(int i = 0; i < paramList.size();i++){
            foutName += paramList[i]+"_";
        }
        foutName += ".pdf";
//c->Print(str2char(foutName));
    
       // saveJSON("testi",names,pList, cutPoints);
   */ 
    return;

}


