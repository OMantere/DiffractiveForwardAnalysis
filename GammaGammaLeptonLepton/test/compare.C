#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "Canvas.h"
#include "TH1.h"
#include "TLorentzVector.h"

#include <iostream>

using namespace std;

const char* display_name(const char* name) {
    if(!strcmp("computed.root", name))
        return "Signal";
    if(!strcmp("computed_ww.root", name))
        return "W+W- Background";
    else
        return name;
}

void compare(
        const char* filename1 = "computed.root",
        const char* filename2 = "computed_pu.root",
        const char* identifier = "Wmiss",
        const char* output_name = "comparison.pdf",
        int bins = -1,
        float minn = -1,
        float maxx = -1
)
{
    double var1, var2;
    TFile file1(filename1);
    TFile file2(filename2);
    TTree* tree1 = (TTree*)file1.Get("computed");
    TTree* tree2 = (TTree*)file2.Get("computed");

    tree1->SetBranchAddress(identifier, &var1);
    tree2->SetBranchAddress(identifier, &var2);

    TCanvas c("Canvas for comparing of quantities");

    int N1 = tree1->GetEntriesFast();
    int N2 = tree2->GetEntriesFast();
    int N = N1;
    if(N1 != N2) {
        cout << "Warning: Your files have a different number of entries: " << N << ", " << tree2->GetEntriesFast() << endl;
        N = min(N1, N2);
    }

    double get_min, get_max;
    for(int i = 0; i < N; i++) {
        tree1->GetEntry(i);
        tree2->GetEntry(i);
        if(var1 < get_min) get_min = var1;
        if(var2 < get_min) get_min = var2;
        if(var1 > get_max) get_max = var1;
        if(var2 > get_max) get_max = var2;
    }

    if(minn == -1) minn = get_min;
    if(maxx == -1) minn = get_max;
    if(bins == -1) bins = N/3;

    TH1D h1("compare1", identifier, bins, minn, maxx);
    TH1D h2("compare2", identifier, bins, minn, maxx);
    h1.SetMaximum(80);
    h2.SetMaximum(80);
    for(int i = 0; i < N; i++) {
        tree1->GetEntry(i);
        tree2->GetEntry(i);
        h1.Fill(var1);
        h2.Fill(var2);
    }
    h1.Draw();
    h1.SetLineColor(kRed);
    h1.GetXaxis()->SetTitle("GeV");
    h1.GetYaxis()->SetTitle("Events");
    h2.Draw("same");
    h2.GetXaxis()->SetTitle("GeV");
    h2.GetYaxis()->SetTitle("Events");
    h2.SetLineColor(kBlue);
    c.Update();
    auto legend = new TLegend(0.1,0.7,0.48,0.9);
//    legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry("compare1",display_name(filename1),"l");
    legend->AddEntry("compare2",display_name(filename2),"l");
    legend->Draw();
    c.Update();
    c.Print(output_name);
}