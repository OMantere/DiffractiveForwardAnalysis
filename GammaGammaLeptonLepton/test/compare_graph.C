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

void compare_graph(
        const char* filename1 = "computed.root",
        const char* filename2 = "computed_pu.root",
        const char* identifier1 = "Wgg",
        const char* identifier2 = "Wmiss",
        const char* output_name = "comparison_graph.pdf"
)
{
    double var1a, var2a, var1b, var2b;
    TFile file1(filename1);
    TFile file2(filename2);
    TTree* tree1 = (TTree*)file1.Get("computed");
    TTree* tree2 = (TTree*)file2.Get("computed");

    tree1->SetBranchAddress(identifier1, &var1a);
    tree1->SetBranchAddress(identifier2, &var1b);
    tree2->SetBranchAddress(identifier1, &var2a);
    tree2->SetBranchAddress(identifier2, &var2b);

    TCanvas c("Canvas for comparing of quantities");

    int N1 = tree1->GetEntriesFast();
    int N2 = tree2->GetEntriesFast();
    int N = N1;
    if(N1 != N2) {
        cout << "Warning: Your files have a different number of entries: " << N << ", " << tree2->GetEntriesFast() << endl;
        N = min(N1, N2);
    }
    double* arr1a = new double[N];
    double* arr2a = new double[N];
    double* arr1b = new double[N];
    double* arr2b = new double[N];
    for(int i = 0; i < N; i++) {
        tree1->GetEntry(i);
        tree2->GetEntry(i);
        arr1a[i] = var1a;
        arr2a[i] = var2a;
        arr1b[i] = var1b;
        arr2b[i] = var2b;
    }
    TGraph h1(N, arr1a, arr1b);
    h1.SetName("compare1");
    h1.SetTitle(((string)identifier1 + "-" + (string)identifier2).c_str());
    h1.Draw("ap");
    h1.SetMarkerColor(kRed);
    h1.SetLineColor(kRed);
    h1.SetMarkerSize(5);
    c.Update();
    TGraph h2(N, arr2a, arr2b);
    h2.SetName("compare2");
    h2.SetTitle(((string)identifier1 + "-" + (string)identifier2).c_str());
    h2.Draw("p");
    h2.SetMarkerColor(kBlue);
    h2.SetLineColor(kBlue);
    h2.SetMarkerSize(5);
    c.Update();
    auto legend = new TLegend(0.1,0.7,0.48,0.9);
//    legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry("compare1",display_name(filename1),"l");
    legend->AddEntry("compare2",display_name(filename2),"l");
    legend->Draw();
    c.Update();
    c.Print(output_name);
}
