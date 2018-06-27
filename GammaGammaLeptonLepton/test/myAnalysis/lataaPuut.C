//#include "DiffarctiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"
#include "Canvas.h"
#include "TH1.h"
#include "TLorentzVector.h"

#include <iostream>

using namespace std;

void lataaPuut( const char *file1 = "../computedLM6.root"){
    TFile f1(file1);

    TTree *tree1 = (TTree *) f1.Get("computed");

}
