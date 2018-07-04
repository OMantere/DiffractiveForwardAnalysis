#include <cfloat>
#include <time.h>
#include "TH1.h"

clock_t start = clock(), diff;

void start_time() {
    start = clock();
}

void end_time() {
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds", msec/1000, msec%1000);
}

class VariableCut {
public:
    VariableCut(string var, double low, double up, double plow, double pup) : variable(var), low_bound(low), up_bound(up), low_plot(plow), up_plot(pup) {};
    VariableCut(string var, double low, double up) : variable(var), low_bound(low), up_bound(up), low_plot(0), up_plot(1000) {};
    VariableCut(string var, double low) : variable(var), low_bound(low), up_bound(DBL_MAX), low_plot(0), up_plot(1000) {};
    VariableCut(string var) : variable(var), low_bound(DBL_MIN), up_bound(DBL_MAX), low_plot(0), up_plot(1000) {};
    string variable;
    double low_bound;
    double up_bound;
    double low_plot;
    double up_plot;

    bool pass(double val) {
        return val < up_bound && val >= low_bound;
    }
};

int max_events = 10000;

vector<TH1F> load_files(vector<string> files, vector<VariableCut> cuts) {
    int n_vars = cuts.size();
    int n_files = files.size();
    vector<double> val(n_vars);
    vector<TH1F> histos;
    const char* c_names[n_vars];
    for(int i = 0; i < n_vars; i++) {
        c_names[i] = cuts[i].variable.c_str();
    }
    for(int i = 0; i < n_vars; i++) {
        histos.push_back(TH1F(c_names[i], c_names[i], 20, cuts[i].low_plot, cuts[i].up_plot));
    }
    for(int i = 0; i < n_files; i++) {
        TFile f(files[i].c_str());
        TTree* t = (TTree*)f.Get("computed");
        t->SetCacheSize(3010241024);
        for(int k = 0; k < n_vars; k++) {
            t->SetBranchAddress(c_names[k], &val[k]);
        }
        int N = t->GetEntriesFast();
        for(int j = 0; j < N; j++) {
            t->GetEntry(j);
            for(int k = 0; k < n_vars; k++) {
                histos[i].Fill(val[i]);
                if(!cuts[i].pass(val[i]))
                    continue;
            }
        }
    }
    return histos;
}

void l1_analysis() {
    vector<string> signal_files = {
            "computed_slr2.root",
            "computed_sll3.root",
    };

    vector<string> bg_files = {
            "computed_ww.root",
            "computed_dy.root",
            "computed_ggll_elastic.root"
    };

    vector<VariableCut> l1_cuts = {
            VariableCut("pt1", 20),
            VariableCut("pt2", 20),
    };

    start_time();
    vector<TH1F> histos = load_files(signal_files, l1_cuts);
    end_time();


}