#include <cfloat>

void load_files() {
    vector<string> signal_files = {
            "computed_slr2.root",
            "computed_sll3.root"
    };

    vector<string> bg_files = {
            "computed_ww.root",
            "computed_dy.root",
            "computed_ggll_elastic.root"
    };

    vector<TH1D> hist

    vector<VariableCut> variables = {
            VariableCut("Wlep", 20, 40)
            VariableCut("pt", 20, 40)
    };

    for(auto file : bg_files) {
        cout << "hello" << endl;
    }
    for(auto t : variables) {
        cout << t.a << t.b << endl;
    }
}