#include <time.h>
#include "TH1.h"
#include <string>
#include <iomanip> // setprecision
#include <sstream> // stringstream


double lhc_lumi = 50;
double elastic_crossx = 20.0781 * 1000;
double ww_crossx = 0.7583;
double dy_crossx = 1141 * 1000;
double cms_dy_crossx1 = 1975 * 1000;
double cms_dy_crossx2 = 19.32 * 1000;
double cms_dy_crossx3 = 2.731 * 1000;
double lm1mul_crossx = 0.055;
double lm1mur_crossx = 0.263;
double lm1mur_x10_80_crossx_old = 1.7297;
double lm1mur_x10_80_crossx = 0.2472;
double lm1mur_lm6_crossx = 0.07169;
double lm1mur_lm6_x10_130_crossx = 0.07316;
double ppwz_crossx = 4.429 * 1.109 * 1000;
double ppww_crossx = (118.7 - 3.974) * pow(0.10862, 2) * 9 * 1000;
double dyjets_crossx = 1921.8 * 1000 * 801439.0 / 4999851.0;
double ttbar_crossx = 831.76 * pow(0.10862, 2) * 1000 * 1430956 / 5500000.0;
double ppzz_crossx = 0.564 * 1000;
double signal_crossx = lm1mur_crossx;
double totem_ctpps_fake_rate = 9.28 / 1000;


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

    VariableCut(string var, double plow, double pup) : variable(var), low_plot(plow), up_plot(pup), low_bound(-DBL_MAX),
                                                       up_bound(DBL_MAX), plot(true) {};

    VariableCut(string var) : variable(var), low_bound(-DBL_MAX), up_bound(DBL_MAX), low_plot(0), up_plot(1000),
                              plot(true) {};
    string variable;
    double low_bound;
    double up_bound;
    double low_plot;
    double up_plot;
    bool plot;
    bool inverted = false;

    bool pass(double val) {
        if (inverted)
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
    Sample(string file, string leg, string short_leg, double cx, bool line = false) : filename(file),
                                                                                      legend(leg),
                                                                                      short_legend(short_leg),
                                                                                      crossx(cx),
                                                                                      line(line) {};

    Sample(string file, string leg, double cx, bool line = false) : filename(file), legend(leg),
                                                                    short_legend(leg), crossx(cx),
                                                                    line(line) {};
    string filename;
    string legend;
    string short_legend;
    double crossx;
    bool line = false;
    double factor = 1;
    bool background = true;

    void init() {
        line = false;
        factor = 1;
    }

    void scale(double fac) {
        factor = fac;
    }

    double getFactor() {
        if(background)
            return factor * totem_ctpps_fake_rate;
        else
            return factor;
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

vector <string> display_names;
vector <string> files;
vector<int> colors = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15};
vector<double> crossx;

Sample *elastic = new Sample("computed_elastic.root", "Elastic #gamma#gamma -> #mu_{R}^{+}#mu_{R}^{-}", elastic_crossx);
Sample *slr2 = new Sample("computed_slr2.root", "LM1 #mu_{R}", lm1mur_crossx);
Sample *slr2full = new Sample("computed_lm1full.root", "LM1 #mu_{R} fullsim", "LM1 #mu_{R}", lm1mur_crossx);
Sample *slr2full_lm6 = new Sample("computed_slr2full_lm6.root", "LM1 #mu_{R} fullsim", "LM1 #mu_{R}", lm1mur_crossx);
Sample *slr2_x10_80_old = new Sample("computed_slr2_x10_80_old.root", "LM1 #mu_{R} with m(#chi_{1}^{0}) = 80 GeV)",
                                     "LM1 #mu_{R}",
                                     lm1mur_x10_80_crossx_old);
Sample *slr2_x10_80 = new Sample("computed_slr2_x10_80.root", "LM1 #mu_{R} with m(#chi_{1}^{0}) = 80 GeV)",
                                 "LM1 #mu_{R}",
                                 lm1mur_x10_80_crossx);
Sample *slr2_lm6 = new Sample("computed_slr2_lm6.root", "LM6 #mu_{R}", lm1mur_lm6_crossx);
Sample *slr2_lm6_x10_130 = new Sample("computed_slr2_lm6_x10_130.root", "LM6 #mu_{R} with m(#chi_{1}^{0}) = 130 GeV)",
                                      "LM6 #mu_{R}",
                                      lm1mur_lm6_x10_130_crossx);
Sample *sll3 = new Sample("computed_sll3.root", "LM1 #mu_{L}", lm1mul_crossx);
Sample *aaww = new Sample("computed_ww.root", "#gamma#gamma -> W^{+}W^{-}", "#gamma#gamma -> W^{+}W^{-}", ww_crossx);
//Sample *dy = new Sample("computed_dy_10k.root", "Drell-Yan -> #mu^{+}#mu^{-}", "DY", dy_crossx);
Sample *dy = new Sample("computed_dy_1M.root", "Drell-Yan -> #mu^{+}#mu^{-}", "DY1", cms_dy_crossx1);
Sample *dy2 = new Sample("computed_dy2.root", "Drell-Yan -> #mu^{+}#mu^{-}", "DY2", cms_dy_crossx2);
Sample *dy3 = new Sample("computed_dy3.root", "Drell-Yan -> #mu^{+}#mu^{-}", "DY3", cms_dy_crossx3);
Sample *wz = new Sample("computed_ppwz.root", "pp -> WZ", "WZ", ppwz_crossx);
Sample *zz = new Sample("computed_ppzz.root", "pp -> ZZ", "ZZ", ppzz_crossx);
Sample *ww = new Sample("computed_ppww.root", "pp -> W^{+}W^{-}", "W^{+}W^{-}", ppww_crossx);
Sample *ttbar = new Sample("computed_ttbar.root", "pp -> t#bar{t}", "t#bar{t}", ttbar_crossx);
Sample *dyjets = new Sample("computed_dyjets_big.root", "DY+jets", dyjets_crossx);

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
//        VariableCut("xip", 0.03, 0.15, 0, 0.2, false),
        VariableCut("Wlep", 80, 100, 50, 150, true, true),
        VariableCut("pair_aco", 0.25, 2, 0, 1),
        VariableCut("Wlep", 0, 250)
};

vector <VariableCut> slr2_dy_test = {
        VariableCut("Wlep", 0, 75, 0, 150),
        VariableCut("pt2", 5, DBL_MAX, 0, 60),
        VariableCut("mt2_100", 0, 122, 100, 200),
        VariableCut("Pt", 0, 100),
};

vector <VariableCut> slr2_ttbar_test = {
        VariableCut("mt2_100", 0, 122, 100, 200),
//        VariableCut("Wlep", 0, 75, 0, 150),
        VariableCut("pt1", 10, 50, 0, 100),
        VariableCut("pt2", 6, 40, 0, 60),
        VariableCut("Pt", 0, 150),
//        VariableCut("Pt", 0, 150),
};

vector <VariableCut> slr2_ww_test = {
        VariableCut("mt2_100", 0, 122, 100, 200),
        VariableCut("pt1", 10, 50, 0, 100),
        VariableCut("pt2", 6, 40, 0, 60),
        VariableCut("Pt", 0, 44, 0, 100),
        VariableCut("Wlep", 0, 150),
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

long long max_events = 100000;
int n_bins = 100;
vector<long int> n_events;
int n_vars;
int n_files;
vector <vector<long int>> passed;
map<string, double> val_map;
vector <vector<double>> pass_frac;
vector <vector<TH1F *>> histos;
int n_plots = 0;

void plot_hist_legend(int j, vector<Sample *> samples) {
    auto hist_legend = new TLegend(0.65, 0.5, 0.89, 0.87);
    hist_legend->SetHeader("After cut", "C");
    hist_legend->SetTextSize(0.07);
    hist_legend->SetBorderSize(0);
    for (int k = 0; k < n_files; k++) {
        double pass_fraction = pass_frac[k][j + 1];
        if (j == n_vars - 1)
            pass_fraction = pass_frac[k][j];
        char *str = strdup((samples[k]->short_legend + ": " + to_str(pass_fraction * 100, 3) + " %").c_str());
        if (samples[k]->line)
            hist_legend->AddEntry(histos[k][0], str, "l");
        else
            hist_legend->AddEntry(histos[k][0], str, "f");
    }
    hist_legend->Draw("same");
}

//// TODO
vector <string> dimension_vars = {
        "asdf"
};

//char* x_axis_label(string var) {
//    return "TODO".c_str();
//}

class Analysis {
public:
    Analysis(vector<Sample *> samples, vector<double> scales, vector <VariableCut> cuts, vector<bool> use_line,
             string name, bool propagate = true) : samples(samples), scales(scales), cuts(cuts), use_line(use_line),
                                                   name(name), propagate(propagate) {};
    vector<Sample *> samples;
    vector <VariableCut>cuts;
    vector<double> scales;
    vector<bool> use_line;
    string name;
    bool propagate;

    vector <vector<TH1F *>> load_files() {
        for (int i = 0; i < samples.size(); i++)
            files.push_back(samples[i]->filename);
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
            vector<long int> pass_v;
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
            const char* file_name = files[i].c_str();
            cout << "Reading " << file_name << "..." << endl;
            TFile f(file_name);
            TTree *t = (TTree *) f.Get("computed");
            for (int k = 0; k < n_vars; k++) {
                t->SetBranchAddress(c_names[k], &val_map.find(cuts[k].variable)->second);
            }
            long int N = std::min(t->GetEntriesFast(), max_events);
            n_events.push_back(N);
            for (int j = 0; j < N; j++) {
                t->GetEntry(j);
                for (int k = 0; k < n_vars; k++) {
                    histos[i][k]->Fill(val_map.find(cuts[k].variable)->second);
                    if (propagate) {
                        passed[i][k]++;
                        if (!cuts[k].pass(val_map.find(cuts[k].variable)->second))
                            break;
                    } else {
                        if(k == n_vars - 1)
                            break;
                        if (!cuts[k].pass(val_map.find(cuts[k].variable)->second))
                            passed[i+1][k]++;
                    }
                }
            }
        }
        return histos;
    }

    void run() {
        if (scales.size() > samples.size()) {
            cout << "There are more scale factors than samples" << endl;
            return;
        }
        if (use_line.size() > samples.size()) {
            cout << "There are more line settings than samples" << endl;
            return;
        }
        cuts.push_back(VariableCut("Pt", -DBL_MAX, DBL_MAX, 0, 1, false));
        for (Sample *sample : samples)
            sample->init();
        for (int i = 0; i < scales.size(); i++) {
            samples[i]->scale(scales[i]);
        }
        for (int i = 0; i < use_line.size(); i++) {
            samples[i]->line = use_line[i];
        }
        for (Sample *sample : samples) {
            crossx.push_back(sample->crossx);
            display_names.push_back(sample->legend);
        }

        start_time();
        histos = load_files();
        //    vector<TH1F> bg_histos = load_files(bg_files, cuts);
        end_time();

        vector < THStack * > sh;
        for (auto cut : cuts) {
            const char *c_name = cut.variable.c_str();
            sh.push_back(new THStack(c_name, c_name));
        }
        //    for(int j = 0; j < n_vars; j++) {
        //        char* canv_name = strdup(("c" + to_string(j)).c_str());
        //        canvas.push_back(new TCanvas(canv_name,"stacked hists",700,900));
        //    }
        TCanvas *canvas = new TCanvas(name.c_str(), "stacked hists", 750, 1200);
        auto legend = new TLegend(0.15, 0.15, 0.85, 0.85);
        vector<double> h_max(n_vars, 0.0);
        vector<double> expected_events(n_files);
        double bg_event_sum = 0;
        double signal_event_sum = 0;
        cout << endl;
        cout << endl;
        for (int i = 0; i < n_files; i++) {
            pass_frac.push_back(vector<double>(0));
            for (int j = 0; j < n_vars; j++) {
                histos[i][j]->Scale(lhc_lumi * crossx[i] / n_events[i] * samples[i]->getFactor());
                if (samples[i]->line) {
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
            if(i == 0) {
                expected_events.push_back(lhc_lumi * pass_frac[i][n_vars - 1] * samples[i]->crossx);
                signal_event_sum += expected_events.back();
            }
            else {
                expected_events.push_back(lhc_lumi * pass_frac[i][n_vars - 1] * samples[i]->crossx * totem_ctpps_fake_rate);
                bg_event_sum += expected_events.back();
            }
            cout << " %, total expected events: " << to_str(expected_events.back(), 2) << endl;
            if (samples[i]->line)
                legend->AddEntry(histos[i][0], display_names[i].c_str(), "l");
            else
                legend->AddEntry(histos[i][0], display_names[i].c_str(), "f");
        }
        cout << endl;
        cout << "S/B ratio: " << signal_event_sum/bg_event_sum * 100 << " %" << endl;
        cout << endl;
        cout << endl;

        canvas->Divide(1, n_plots + 1);
        canvas->cd(1);
        legend->SetTextSize(0.07);
        legend->Draw();
        int i_plt = 0;
        for (int j = 0; j < n_vars; j++) {
            if (!cuts[j].plot)
                continue;
            canvas->cd(i_plt + 2);
            int k = 0;
            if (sh[j]->GetNhists() == 0) {
                histos[0][j]->SetStats(0);
                histos[0][j]->SetMaximum(h_max[j] + h_max[j] * 0.1);
                histos[0][j]->SetTitle(cuts[j].variable.c_str());
                histos[0][j]->Draw("HIST");
                k = 1;
            } else {
                sh[j]->SetMaximum(h_max[j] + h_max[j] * 0.1);
                sh[j]->Draw("HIST");
            }
            while (k < n_files) {
                if (samples[k]->line)
                    histos[k][j]->Draw("HIST SAME");
                k++;
            }
            gStyle->SetTitleFontSize(0.07);
            canvas->Update();
            //        plot_text("Events ")
            int y_max = 1e9;
            TLine *low_line = new TLine(cuts[j].low_bound, 0, cuts[j].low_bound, y_max);
            TLine *up_line = new TLine(cuts[j].up_bound, 0, cuts[j].up_bound, y_max);
            low_line->SetLineWidth(2);
            up_line->SetLineWidth(2);
            low_line->Draw("same");
            up_line->Draw("same");
            plot_hist_legend(j, samples);
            canvas->Update();
            char *x_axis = strdup(cuts[j].variable.c_str());
            char *y_axis = strdup(("Events / " + to_str(cuts[j].plot_range() / n_bins, 1) + " GeV").c_str());
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
        canvas->Print((name + ".pdf").c_str());
    }
};


// Analysis definitions
Analysis *both_ww_analysis = new Analysis({slr2, ww, aaww}, {4.7 * 1e3 / 5, 1, 4.7 * 1e3 / 5}, slr2_ww_test, {true},
                                          "both_ww_reduction");

//2.5e+05
//Sample computed_ppzz.root reduced by: 0.000000 %, total expected events: 2.8e+04
//Sample computed_ttbar.root reduced by: 0.000000 %, total expected events: 4.4e+06
//Sample computed_dyjets.root reduced by: 0.000000 %, total expected events: 2.9e+08
vector<Sample *> all_samples = {slr2full, ttbar, dyjets};
vector<double> all_scales = {
        1.0/20.0,
//    13.0/38.0,
        13 / 6.1e+05,
        13 / 2.5e+05,
        13 / 2.8e+04,
        13 / 4.4e+06,
        13 / 2.9e+08,
};
vector<bool> all_lines = {true, true, true, true, true, true, true};
Analysis *all_analysis_new = new Analysis(
        all_samples,
        all_scales,
        {
                VariableCut("HTlep", 0, 400),
                VariableCut("pt1", 0, 100),
        },
        all_lines,
        "all_analysis_new"
);

Analysis *all_analysis_eta = new Analysis(
        all_samples,
        all_scales,
        {
                VariableCut("eta1", -2.5, 2.5),
                VariableCut("eta2", -2.5, 2.5),
        },
        all_lines,
        "all_analysis_eta"
);

Analysis *all_analysis_new2 = new Analysis(
        all_samples,
        all_scales,
        {
            VariableCut("Mt", 0, 300),
            VariableCut("mt2_100", 100, 150),
        },
        all_lines,
        "all_analysis_new2"
);


/// July 25th: Cuts don't seem to affect the XiP distribution for signal
Analysis *all_analysis_extra = new Analysis(
        all_samples,
        {},
        {
            VariableCut("pt1", 10, 50, 0, 150),
            VariableCut("pt2", 10, 50, 0, 150),
            VariableCut("Wlep", -DBL_MAX, 75, 0, 150),
            VariableCut("extratracks2", -DBL_MAX, 3, 0, 100),
            VariableCut("Pt", 0, 150),
        },
        {},
        "all_analysis_extra"
);

/// 4th August:
/// This is cool, bump global PT cut from 10 -> 12 results in 5 % -> 11 % S/B ratio... <3
/// With GGLL acceptance (kind of hand wavy right now) for BGs we get 60 %... :D
double pt_cut2 = 12.0;
Analysis *all_analysis_kristian = new Analysis(
        all_samples,
        {},
        {
            VariableCut("closest_extra", 0.01, DBL_MAX, 0, 0.2, false),
            VariableCut("pt1", pt_cut2, 200, 0, 150, false),
            VariableCut("Wlep", 0, 75, 0, 150, false),
            VariableCut("pair_dphi", -2.5, 2.5, -3.14, 3.14, false),
//            VariableCut("closest_extra", 0.01, DBL_MAX, 0, 0.2),
            VariableCut("pt2", pt_cut2, 200, 0, 150, false),
            VariableCut("extratracks2", -DBL_MAX, 3, 0, 100, false),
//            VariableCut("xiP1", 0.02, 0.15, 0, 1),
//            VariableCut("xiP2", 0.02, 0.15, 0, 1),
            VariableCut("xiP1", 0.02, 0.15, 0,1),
            VariableCut("xiP2", 0.02, 0.15, 0,1)
        },
        {true},
        "all_analysis_kristian"
);

Analysis *test_vars = new Analysis(
        all_samples,
        {},
        {
            VariableCut("closest_extra", 0.01, DBL_MAX, 0, 0.2),
            VariableCut("extratracks2", -DBL_MAX, 3, 0, 100),
        },
        {true},
        "test_vars"
);

double pt_cut = 4.0;
Analysis *pt_analysis = new Analysis(
        all_samples,
        {},
        {
            VariableCut("pt1", pt_cut, DBL_MAX, 0, 100),
            VariableCut("pt2", pt_cut, DBL_MAX, 0, 100),
            VariableCut("eta2", -3, 3),
        },
        {true},
        "pt_analysis"
);

//VariableCut("mt2_100", 0, 122, 100, 200),
//        VariableCut("pt1", 10, 50, 0, 100),
//        VariableCut("pt2", 6, 40, 0, 60),
//        VariableCut("Pt", 0, 44, 0, 100),
//        VariableCut("Wlep", 0, 150),


void l1_analysis() {
    slr2->background = false;
    slr2full->background = false;
    slr2full_lm6->background = false;
    slr2_lm6->background = false;
//    all_analysis_new->run();
//    all_analysis_new2->run();
//    all_analysis_eta->run();
    all_analysis_kristian->run();
//    test_vars->run();
//    pt_analysis->run();
}

