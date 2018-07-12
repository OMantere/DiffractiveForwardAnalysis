#include <cfloat>
#include <time.h>
#include "TH1.h"
#include <string>
#include <iomanip> // setprecision
#include <sstream> // stringstream

clock_t start = clock(), diff;

void start_time() {
    start = clock();
}

void end_time() {
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds\n", msec / 1000, msec % 1000);
}

string to_str(double val, int n, bool scient = false) {
    stringstream stream;
    if (scient)
        stream << scientific << setprecision(n) << val;
    else
        stream << setprecision(n) << val;
    return stream.str();
}

class VariableCut {
public:
    VariableCut(string var, double low, double up, double plow, double pup, bool plt, bool inverted) : variable(var),
                                                                                                       low_bound(low),
                                                                                                       up_bound(up),
                                                                                                       low_plot(plow),
                                                                                                       up_plot(pup),
                                                                                                       plot(plt),
                                                                                                       inverted(
                                                                                                               inverted) {};

    VariableCut(string var, double low, double up, double plow, double pup, bool plt) : variable(var), low_bound(low),
                                                                                        up_bound(up), low_plot(plow),
                                                                                        up_plot(pup), plot(plt) {};

    VariableCut(string var, double low, double up, double plow, double pup) : variable(var), low_bound(low),
                                                                              up_bound(up), low_plot(plow),
                                                                              up_plot(pup), plot(true) {};

    VariableCut(string var, double plow, double pup) : variable(var), low_plot(plow), up_plot(pup), low_bound(DBL_MIN),
                                                       up_bound(DBL_MAX), plot(true) {};

    VariableCut(string var) : variable(var), low_bound(DBL_MIN), up_bound(DBL_MAX), low_plot(0), up_plot(1000),
                              plot(true) {};
    string variable;
    double low_bound;
    double up_bound;
    double low_plot;
    double up_plot;
    bool plot;
    bool inverted = false;

    bool pass(double val) {
        if(inverted)
            return val > up_bound || val < low_bound;
        else
            return val < up_bound && val > low_bound;
    }

    double plot_range() {
        return up_plot - low_plot;
    }
};

class Sample {
public:
    Sample(string file, string leg, string short_leg, double cx, int c, bool line = false) : filename(file),
                                                                                             legend(leg),
                                                                                             short_legend(short_leg),
                                                                                             crossx(cx), color(c),
                                                                                             line(line) {};

    Sample(string file, string leg, double cx, int c, bool line = false) : filename(file), legend(leg),
                                                                           short_legend(leg), crossx(cx), color(c),
                                                                           line(line) {};
    string filename;
    string legend;
    string short_legend;
    double crossx;
    int color;
    bool line;
    double factor = 1;

    void scale(double fac) {
        factor = fac;
    }
};

void plot_text(string text) {
    TText *xlabel = new TText();
    xlabel->SetNDC();
    xlabel->SetTextFont(1);
    xlabel->SetTextColor(1);
    xlabel->SetTextSize(0.03);
    xlabel->SetTextAlign(22);
    xlabel->SetTextAngle(0);
    xlabel->DrawText(0.5, 0.04, text.c_str());
}

double lhc_lumi = 50;
double elastic_crossx = 20.0781 * 1000;
double ww_crossx = 0.7583;
double dy_crossx = 1141 * 1000;
double lm1mul_crossx = 0.055;
double lm1mur_crossx = 0.263;
double lm1mur_x10_80_crossx_old = 1.7297;
double lm1mur_x10_80_crossx = 0.2472;
double lm1mur_lm6_crossx = 0.07169;
double lm1mur_lm6_x10_130_crossx = 0.07316;
double signal_crossx = lm1mur_crossx;

vector <string> display_names;
vector <string> files;
vector<int> colors;
vector<double> crossx;

Sample *elastic = new Sample("computed_elastic.root", "Elastic #gamma#gamma -> #mu_{R}^{+}#mu_{R}^{-}", elastic_crossx,
                             4);
Sample *slr2 = new Sample("computed_slr2.root", "LM1 #mu_{R}", lm1mur_crossx, 12);
Sample *slr2_x10_80_old = new Sample("computed_slr2_x10_80_old.root", "LM1 #mu_{R} with m(#chi_{1}^{0}) = 80 GeV)",
                                     "LM1 #mu_{R}",
                                     lm1mur_x10_80_crossx_old, 2);
Sample *slr2_x10_80 = new Sample("computed_slr2_x10_80.root", "LM1 #mu_{R} with m(#chi_{1}^{0}) = 80 GeV)",
                                 "LM1 #mu_{R}",
                                 lm1mur_x10_80_crossx, 6);
Sample *slr2_lm6 = new Sample("computed_slr2_lm6.root", "LM6 #mu_{R}", lm1mur_lm6_crossx, 7);
Sample *slr2_lm6_x10_130 = new Sample("computed_slr2_lm6_x10_130.root", "LM6 #mu_{R} with m(#chi_{1}^{0}) = 130 GeV)",
                                      "LM6 #mu_{R}",
                                      lm1mur_lm6_x10_130_crossx, 10);
Sample *sll3 = new Sample("computed_sll3.root", "LM1 #mu_{L}", lm1mul_crossx, 11);
Sample *ww = new Sample("computed_ww.root", "W^{+}W^{-} -> #mu^{+}#mu^{-}", ww_crossx, 9);
Sample *dy = new Sample("computed_dy_10k.root", "Drell-Yan -> #mu^{+}#mu^{-}", "DY", dy_crossx, 8);

vector <VariableCut> l1_cuts = {
        VariableCut("Pt", 2, 10000, 0, 50, false),
        VariableCut("pair_eta_diff", -100, 3, 0, 5, false),
        VariableCut("Pt", -100, 50, 0, 150, false),
        VariableCut("Et", -100, 70, 0, 150, false),
        VariableCut("pt1", -100, 10000, 0, 150),
        VariableCut("pt2", -100, 10000, 0, 150),
        VariableCut("Mt", -100, 10000, 0, 150),
//        VariableCut("pt1", -100, 50, 0, 500),
//        VariableCut("pt2", -100, 50, 0, 500, false),
//        VariableCut("Emiss", 200, 10000, 0, 1000),
//        VariableCut("Emiss", 200, 10000, 0, 500),
//        VariableCut("pt2", 50, 10000, 0, 500)
        VariableCut("Wlep", -1e9, 1e9, 0, 150, false),
};

vector <VariableCut> kristian_lm1_cuts = {
        VariableCut("pair_eta_diff", -100, 3, 0, 5, false),
        VariableCut("Pt", -100, 50, 0, 150, false),
        VariableCut("pt1", -100, 50, 0, 150, false),
        VariableCut("pt2", -100, 40, 0, 150, false),
        VariableCut("Emiss", 200, 10000000, 150, 500),
        VariableCut("Wlep", -10000, 100000, 0, 100, false)
};

vector <VariableCut> laurent_cuts = {
        VariableCut("Pt", 2, 1e9, 0, 10),
        VariableCut("Pt", -100, 10000, 0, 150),
};

vector <VariableCut> xip_cuts = {
        VariableCut("xip", 0.03, 0.15, 0, 0.2),
        VariableCut("xip", -1, 1, 0, 0.2)
};

vector <VariableCut> dy_test = {
        VariableCut("xip", 0.03, 0.15, 0, 0.2),
        VariableCut("Wlep", 80, 100, 50, 150, true, true),
        VariableCut("pair_aco", 0.10, 2, 0, 1),
        VariableCut("Wlep", -10, 2000, 0, 1,false)
};

//vector <Sample> all_samples = {
////        dy,
////        ww,
////        elastic,
//        slr2,
//        slr2_lm6,
//        slr2_lm6_x10_130,
//        slr2_x10_80,
//        sll3
//};
//
//vector <Sample> signal_samples = {
//        dy,
//        ww,
////        elastic,
//        slr2,
//        slr2_lm6,
//        sll3
//};

// Change
vector <VariableCut> use_cuts = dy_test;
const char *pdfname = "dy_reduction.pdf";
vector<Sample *> use_samples = {dy, slr2, slr2_lm6};

int max_events = 10000;
int n_bins = 80;
vector<int> n_events;
int n_vars;
int n_files;
vector <vector<int>> passed;
map<string, double> val_map;
vector <vector<double>> pass_frac;
vector <vector<TH1F *>> histos;
int n_plots = 0;

vector <vector<TH1F *>> load_files(vector <string> files, vector <VariableCut> cuts) {
    n_vars = cuts.size();
    n_files = files.size();
    vector <vector<TH1F *>> histos;
    const char *c_names[n_vars];
    for (int i = 0; i < n_vars; i++) {
        const char *c_name = cuts[i].variable.c_str();
        c_names[i] = c_name;
        val_map.insert(pair<string, double>(cuts[i].variable, 0.0));
        if (cuts[i].plot)
            n_plots++;
    }
    for (int i = 0; i < n_files; i++) {
        vector < TH1F * > histo_v;
        vector<int> pass_v;
        for (int j = 0; j < n_vars; j++) {
            char *histo_name = strdup((cuts[j].variable + "_" + files[i]).c_str());
            TH1F *histo = new TH1F(histo_name, histo_name, n_bins, cuts[j].low_plot, cuts[j].up_plot);
            histo_v.push_back(histo);
            pass_v.push_back(0);
        }
        histos.push_back(histo_v);
        passed.push_back(pass_v);
    }
    for (int i = 0; i < n_files; i++) {
        TFile f(files[i].c_str());
        TTree *t = (TTree *) f.Get("computed");
        for (int k = 0; k < n_vars; k++) {
            t->SetBranchAddress(c_names[k], &val_map.find(cuts[k].variable)->second);
        }
        int N = t->GetEntriesFast();
        n_events.push_back(N);
        for (int j = 0; j < N; j++) {
            t->GetEntry(j);
            for (int k = 0; k < n_vars; k++) {
                histos[i][k]->Fill(val_map.find(cuts[k].variable)->second);
                passed[i][k]++;
                if (!cuts[k].pass(val_map.find(cuts[k].variable)->second))
                    break;
            }
        }
    }
    return histos;
}

void plot_hist_legend(int j) {
    auto hist_legend = new TLegend(0.65, 0.5, 0.89, 0.87);
    hist_legend->SetHeader("After cut", "C");
    hist_legend->SetTextSize(0.07);
    hist_legend->SetBorderSize(0);
    for (int k = 0; k < n_files; k++) {
        double pass_fraction = pass_frac[k][j + 1];
        if (j == n_vars - 1)
            pass_fraction = pass_frac[k][j];
        char *str = strdup((use_samples[k]->short_legend + ": " + to_str(pass_fraction * 100, 3) + " %").c_str());
        if (use_samples[k]->line)
            hist_legend->AddEntry(histos[k][0], str, "l");
        else
            hist_legend->AddEntry(histos[k][0], str, "f");
    }
    hist_legend->Draw("same");
}

void l1_analysis() {
    // Scale samples
//    dy->scale(5 / 2.7e+06);
    dy->line = true;
    slr2->line = true;
    slr2_lm6->line = true;

    for (Sample *sample : use_samples) {
        crossx.push_back(sample->crossx);
        colors.push_back(sample->color);
        files.push_back(sample->filename);
        display_names.push_back(sample->legend);
    }

    start_time();
    histos = load_files(files, use_cuts);
//    vector<TH1F> bg_histos = load_files(bg_files, use_cuts);
    end_time();

    vector < THStack * > sh;
    for (auto cut : use_cuts) {
        const char *c_name = cut.variable.c_str();
        sh.push_back(new THStack(c_name, c_name));
    }
//    for(int j = 0; j < n_vars; j++) {
//        char* canv_name = strdup(("c" + to_string(j)).c_str());
//        canvas.push_back(new TCanvas(canv_name,"stacked hists",700,900));
//    }
    TCanvas *canvas = new TCanvas("name", "stacked hists", 750, 1200);
    auto legend = new TLegend(0.15, 0.15, 0.85, 0.85);
    vector<double> h_max(n_vars, 0.0);
    cout << endl;
    cout << endl;
    for (int i = 0; i < n_files; i++) {
        pass_frac.push_back(vector<double>(0));
        for (int j = 0; j < n_vars; j++) {
            histos[i][j]->Scale(lhc_lumi * crossx[i] / n_events[i] * use_samples[i]->factor);
            if (use_samples[i]->line) {
                histos[i][j]->SetLineColor(colors[i]);
                histos[i][j]->SetLineWidth(3);
                double this_h_max = histos[i][j]->GetMaximum();
                if (h_max[j] < this_h_max)
                    h_max[j] += this_h_max;
            } else {
                histos[i][j]->SetLineColor(kBlack);
                histos[i][j]->SetLineWidth(1);
                histos[i][j]->SetFillColor(colors[i]);
                sh[j]->Add(histos[i][j]);
                h_max[j] += histos[i][j]->GetMaximum();
            }
            pass_frac.back().push_back(((double) passed[i][j] / (double) n_events[i]));
        }
        cout << "Sample " + files[i] << " reduced by: " << to_string((1.0 - pass_frac[i][n_vars - 1]) * 100);
        cout << " %, total expected events: " << to_str(lhc_lumi * pass_frac[i][n_vars - 1] * use_samples[i]->crossx, 2)
             << endl;
        if (use_samples[i]->line)
            legend->AddEntry(histos[i][0], display_names[i].c_str(), "l");
        else
            legend->AddEntry(histos[i][0], display_names[i].c_str(), "f");
    }
    cout << endl;
    cout << endl;

    canvas->Divide(1, n_plots + 1);
    canvas->cd(1);
    legend->SetTextSize(0.07);
    legend->Draw();
    int i_plt = 0;
    for (int j = 0; j < n_vars; j++) {
        if (!use_cuts[j].plot)
            continue;
        canvas->cd(i_plt + 2);
        int k = 0;
        if (sh[j]->GetNhists() == 0) {
            histos[0][j]->SetStats(0);
            histos[0][j]->SetMaximum(h_max[j] + h_max[j] * 0.1);
            histos[0][j]->SetTitle(use_cuts[j].variable.c_str());
            histos[0][j]->Draw("HIST");
            k = 1;
        } else {
            sh[j]->SetMaximum(h_max[j] + h_max[j] * 0.1);
            sh[j]->Draw("HIST");
        }
        while (k < n_files) {
            if (use_samples[k]->line)
                histos[k][j]->Draw("HIST SAME");
            k++;
        }
        gStyle->SetTitleFontSize(0.07);
        canvas->Update();
//        plot_text("Events ")
        int y_max = 1e9;
        TLine *low_line = new TLine(use_cuts[j].low_bound, 0, use_cuts[j].low_bound, y_max);
        TLine *up_line = new TLine(use_cuts[j].up_bound, 0, use_cuts[j].up_bound, y_max);
        low_line->SetLineWidth(2);
        up_line->SetLineWidth(2);
        low_line->Draw("same");
        up_line->Draw("same");
        plot_hist_legend(j);
        canvas->Update();
        char *x_axis = strdup((use_cuts[j].variable + " [GeV]").c_str());
        char *y_axis = strdup(("Events / " + to_str(use_cuts[j].plot_range() / n_bins, 1) + " GeV").c_str());
        if (sh[j]->GetNhists() == 0) {
            histos[0][j]->GetXaxis()->SetTitle(x_axis);
            histos[0][j]->GetYaxis()->SetTitle(y_axis);
        } else {
            sh[j]->GetXaxis()->SetTitle(x_axis);
            sh[j]->GetYaxis()->SetTitle(y_axis);
        }
        canvas->Modified();
        i_plt++;
    }
    canvas->Print(pdfname);
}