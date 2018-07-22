#include <math>

using namespace std;

void merge() {
    vector<string> input = {"dy_1M", "dy2", "dy3"};
    double cms_dy_crossx1 = 1975 * 1000;
    double cms_dy_crossx2 = 19.32 * 1000;
    double cms_dy_crossx3 = 2.731 * 1000;
    vector<double> crossx = {cms_dy_crossx1, cms_dy_crossx2, cms_dy_crossx3};
    double ref_crossx = crossx[0];

    const char* output_filename = "computed_dy_cms.root";

    TFile output_file(output_filename, "recreate");
    auto output_tree = new TTree("computed", "Computed quantities from the ROOT analysis");

    TFile reference_file(input[0]);
    auto reference_tree = dynamic_cast<TTree *>( file.Get("computed"));
    int n_ref = reference_tree->GetEntriesFast();

    for(int i = 1; i < input.size(); i++) {
        int n_events = (int)floor(crossx[i] / ref_crossx * n_ref);
    }
}