#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "Canvas.h"
#include "TH1.h"
#include "TLorentzVector.h"

#include <iostream>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>

#include <time.h>

clock_t start = clock(), diff;

void start_time() {
    start = clock();
}

void end_time() {
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
}

using namespace std;

void draw_canvas( TH1D h, const char* histo_name, const char* histo_text) {
  Canvas c(histo_name, histo_text);
  h.Draw();
  c.Prettify(&h);
  c.Save( "pdf" );
}

random_device rd;
mt19937 e2(rd());

double smeared(double q) {
  // This sucks try to come up with something better
  normal_distribution<double> smear(0, max(1.5, q/100));
  return q + smear(e2);
}

char delimiter = ',';
int max_events = 10000;
int max_vars = 50;
int n_vars = 0;
int n_rows = 0;
vector<vector<double>> var_v(max_events, vector<double>(max_vars));
std::map<string, int> var_map;
vector<std::pair<string, int>> keypairs;
vector<char*> keys;

void store_var(string name, double val) {
  if(var_map.count(name) == 0) {
    var_map.insert(std::pair<string,int>(name,n_vars));
    n_vars++;
  }
  var_v[n_rows][var_map.find(name)->second] = val;
}

void write_csv(const char* csv_filename) {
  ofstream ofs(csv_filename, ofstream::out); 
  for(int j = 0; j < n_vars; j++) {
    ofs << keys[j] << delimiter << "\n";
  }
  for(int i = 0; i < n_rows; i++) {
    for(int j = 0; j < n_vars; j++) {
        ofs << var_v[i][j] << delimiter << "\n";
    }
  }
  ofs.close();
}

void reader(const char* c_name = "elastic")
{
  string name(c_name);
  string cfilename = "../ggll_"+name+".root";
  string cfilename2 = "computed_"+name+".root";
  string ccsv_filename= name+".root";
  const char* filename = cfilename.c_str();
  const char* filename2 = cfilename2.c_str();
  const char* csv_filename=ccsv_filename.c_str();
  TFile file( filename);
  auto tree = dynamic_cast<TTree*>( file.Get( "ggll_aod/ntp1") );
  int N = tree->GetEntriesFast();
  ggll::AnalysisEvent evt;

  string type = "Muon";
//  if(filename.find("e.root") != string::npos)
//    type = "Electron";

  if(type == "Muon")
    evt.load( tree, ggll::DiMuon, true );
  else
    evt.load( tree, ggll::DiElectron, true );

//  auto newtree1 = new TTree("ggll", "The AnalysisEvent from GGLL");
//  if(type == "Electron")
//    evt.attach(newtree1, ggll::DiElectron, true);
//  else
//    evt.attach(newtree1, ggll::DiMuon, true);

//  TH1D h_pair_mass( "h_pair_mass", "m(#mu^{+}#mu^{-})\\Events\\GeV?.2f", 200, 0., 1000. );
//  TH1D h_pair_dphi( "h_pair_dphi", "m(#mu^{+}#mu^{-})", 200, -1. * M_PI, 1. * M_PI );
//  TH1D h_extratracks( "h_extratracks", "m(#mu^{+}#mu^{-})", 200, 0.,  10);
//  TH1D h_kvc_z( "h_kvc_z", "m(#mu^{+}#mu^{-})", 200, 0., 15. );
//  TH1D h_susy_mass( "h_susy_mass", "m(#mu^{+}#mu^{-})", 200, 0., 1000. );

  double xim;
  double xip;
  double neutralino_mass;
  double Wmiss;
  double Wgg;
  double mreco;
  double susymass;
  double Wlep;
  double pair_aco, pho_aco, slep_aco;

  double* WmissA = new double[N];
  double* WggA = new double[N];

//  TTree *tree = new TTree("analysis", "The analysis tree");
//  tree->Branch("xip", &xip);
//  tree->Branch("xim",&xim);
//  tree->Branch("neutralino_mass",&neutralino_mass);

  cout << "Entries: " << N << endl;
  for ( unsigned long long i = 0; i < tree->GetEntriesFast(); ++i ) {
    tree->GetEntry( i );
    //cout << ">>> event " << i << endl;
    /*cout << "- fired triggers:" << endl;
    for ( const auto& hlt : *evt.HLT_Name ) {
      cout << "  *) " << hlt << endl;
    }*/

    TLorentzVector P1(0, 0, sqrt(pow(6500, 2) - pow(0.938, 2)), 6500); // P1, P2: Beam protons of energy 6500
    TLorentzVector P2(0, 0,-sqrt(pow(6500, 2) - pow(0.938, 2)), 6500);
    TLorentzVector p1, p2; // p1, p2: measured protons
    TLorentzVector com(0, 0, 0, 13.e3);
    p1.SetPtEtaPhiE(evt.GenProtCand_pt[0], evt.GenProtCand_eta[0], evt.GenProtCand_phi[0], evt.GenProtCand_e[0]);
    p2.SetPtEtaPhiE(evt.GenProtCand_pt[1], evt.GenProtCand_eta[1], evt.GenProtCand_phi[1], evt.GenProtCand_e[1]);


    for ( unsigned int j = 0; j < evt.nPair; ++j ) {

      // Lets do acceptance cuts
//      double pt_acc = 7;
//      double eta_acc = 2.5;
//      if(evt.MuonCand_eta[l1] > eta_acc || evt.MuonCand_eta[l2] > eta_acc || evt.MuonCand_pt[l1] < pt_acc || evt.MuonCand_pt[l2] < pt_acc)
//        continue;

      //if ( evt.Pair_extratracks0p5mm[j] != 0 ) continue;
      if ( fabs( evt.KalmanVertexCand_z[j] ) > 15. ) continue;
      //if ( 1.-fabs( evt.Pair_dphi[j] )/M_PI > 0.009 ) continue;
      //if ( evt.Pair_mass[j] < 110. ) continue;

      set<string> bgs = {"ww", "elastic", "dy"};

      const unsigned int l1 = evt.Pair_lepton1[j], l2 = evt.Pair_lepton2[j];
      double El1, El2;
      TLorentzVector pl1, pl2, pg1, pg2, pl1g, pl2g, psl1, psl2;
      if(type == "Muon") {
        El1 = evt.MuonCand_e[l1];
        El2 = evt.MuonCand_e[l2];
        if(bgs.find(name) == bgs.end()) {
          xip = ( evt.GenRMuonCand_pt[l1]*exp( +evt.GenRMuonCand_eta[l1] ) + evt.GenRMuonCand_pt[l2]*exp( +evt.GenRMuonCand_eta[l2] ) ) / 13.e3;
          xim = ( evt.GenRMuonCand_pt[l1]*exp( -evt.GenRMuonCand_eta[l1] ) + evt.GenRMuonCand_pt[l2]*exp( -evt.GenRMuonCand_eta[l2] ) ) / 13.e3;  
        } else {
          xip = ( evt.MuonCand_pt[l1]*exp( +evt.MuonCand_eta[l1] ) + evt.MuonCand_pt[l2]*exp( +evt.MuonCand_eta[l2] ) ) / 13.e3;
          xim = ( evt.MuonCand_pt[l1]*exp( -evt.MuonCand_eta[l1] ) + evt.MuonCand_pt[l2]*exp( -evt.MuonCand_eta[l2] ) ) / 13.e3; 
        }
        pl1g.SetPtEtaPhiE(evt.GenMuonCand_pt[0], evt.GenMuonCand_eta[0], evt.GenMuonCand_phi[0], evt.GenMuonCand_e[0]);
        pl2g.SetPtEtaPhiE(evt.GenMuonCand_pt[1], evt.GenMuonCand_eta[1], evt.GenMuonCand_phi[1], evt.GenMuonCand_e[1]);
        pl1.SetPtEtaPhiE(evt.MuonCand_pt[l1], evt.MuonCand_eta[l1], evt.MuonCand_phi[l1], evt.MuonCand_e[l1]);
        pl2.SetPtEtaPhiE(evt.MuonCand_pt[l2], evt.MuonCand_eta[l2], evt.MuonCand_phi[l2], evt.MuonCand_e[l2]);
        psl1.SetPtEtaPhiE(evt.GenRMuonCand_pt[l1], evt.GenRMuonCand_eta[l1], evt.GenRMuonCand_phi[l1], evt.GenRMuonCand_e[l1]);
        psl2.SetPtEtaPhiE(evt.GenRMuonCand_pt[l2], evt.GenRMuonCand_eta[l2], evt.GenRMuonCand_phi[l2], evt.GenRMuonCand_e[l2]);
      } else {
        El1 = evt.EleCand_e[l1];
        El2 = evt.EleCand_e[l2];
        if(bgs.find(name) == bgs.end()) {
          xip = ( evt.GenREleCand_pt[l1]*exp( +evt.GenREleCand_eta[l1] ) + evt.GenREleCand_pt[l2]*exp( +evt.GenREleCand_eta[l2] ) ) / 13.e3;
          xim = ( evt.GenREleCand_pt[l1]*exp( -evt.GenREleCand_eta[l1] ) + evt.GenREleCand_pt[l2]*exp( -evt.GenREleCand_eta[l2] ) ) / 13.e3;
        } else {
          xip = ( evt.EleCand_pt[l1]*exp( +evt.EleCand_eta[l1] ) + evt.EleCand_pt[l2]*exp( +evt.EleCand_eta[l2] ) ) / 13.e3;
          xim = ( evt.EleCand_pt[l1]*exp( -evt.EleCand_eta[l1] ) + evt.EleCand_pt[l2]*exp( -evt.EleCand_eta[l2] ) ) / 13.e3; 
        }
        pl1g.SetPtEtaPhiE(evt.GenEleCand_pt[0], evt.GenEleCand_eta[0], evt.GenEleCand_phi[0], evt.GenEleCand_e[0]);
        pl2g.SetPtEtaPhiE(evt.GenEleCand_pt[1], evt.GenEleCand_eta[1], evt.GenEleCand_phi[1], evt.GenEleCand_e[1]);
        psl1.SetPtEtaPhiE(evt.GenREleCand_pt[l1], evt.GenREleCand_eta[l1], evt.GenREleCand_phi[l1], evt.GenREleCand_e[l1]);
        psl2.SetPtEtaPhiE(evt.GenREleCand_pt[l2], evt.GenREleCand_eta[l2], evt.GenREleCand_phi[l2], evt.GenREleCand_e[l2]);
      }

      // We don't have working photons for now
//      pg1.SetPtEtaPhiE(evt.GenPhotCand_pt[0], evt.GenPhotCand_eta[0], evt.GenPhotCand_phi[0], evt.GenPhotCand_e[0]);
//      pg2.SetPtEtaPhiE(evt.GenPhotCand_pt[1], evt.GenPhotCand_eta[1], evt.GenPhotCand_phi[1], evt.GenPhotCand_e[1]);
      pg1 = P1 - p1;
      pg2 = P2 - p2;
      Wgg = 2 * sqrt(pg1.E()*pg2.E());
      const TLorentzVector lep_pair = pl1 + pl2;
      const TLorentzVector slep_pair = psl1 + psl2;
      const TLorentzVector gen_lep_pair = pl1g + pl2g;
      Wlep = lep_pair.M();
      double Wgenlep = gen_lep_pair.M();
      const TLorentzVector p_miss = pg1+pg2-pl1-pl2;
      Wmiss = p_miss.M();
      WmissA[i] = Wmiss;
      WggA[i] = Wgg;
      double Emiss = p_miss.E();
      double Chi10mass = 97;
      double mumass = 105.658 * 1e-3;
      double pair_dphi = evt.Pair_dphi[j];
      double extratracks = evt.Pair_extratracks0p5mm[j];
      double kvc_z = evt.KalmanVertexCand_z[j];
      double pair_mass = evt.Pair_mass[j];
      double pho_dphi = fabs( pg1.Phi()-pg2.Phi() );
      double closest_extra = evt.ClosestExtraTrack_vtxdxyz[j];
      double closest_hp_extra = evt.ClosestHighPurityExtraTrack_vtxdxyz[j];
      while(pho_dphi < -M_PI) pho_dphi += 2 * M_PI;
      while(pho_dphi > M_PI) pho_dphi -= 2 * M_PI;
      pho_aco = 1.-fabs( pho_dphi )/M_PI;
      slep_aco = 1.-fabs( evt.GenSLRPair_dphi )/M_PI;
      pair_aco = 1.-fabs( evt.Pair_dphi[j] )/M_PI;
      mreco = 0.5*sqrt(pow(Wgg, 2) - sqrt(pow(Wmiss,2) - 4*pow(Chi10mass, 2)) + sqrt(pow(Wlep, 2) - 4*pow(mumass, 2)));
      double mreco2 = 2 * mreco;

      double ptMiss = sqrt(p_miss.Px()*p_miss.Px() + p_miss.Py()*p_miss.Py()); // Missing transfers momentum
      double ptTot = sqrt(p1.Px()*p1.Px() + p1.Py()*p1.Py()) + sqrt(p2.Px()*p2.Px() + p2.Px()*p2.Px()); // Total Pt, calculated from diffracted protons
      double pair_eta_diff = fabs(evt.MuonCand_eta[l1] - evt.MuonCand_eta[l2]);
      double eta1 = pl1.Eta();
      double eta2 = pl2.Eta();
      double pt1 = pl1.Pt();
      double pt2 = pl2.Pt();
      double deltaR = sqrt(pow(pl2.Eta() - pl1.Eta(), 2) + pow(pl2.Phi() - pl1.Phi(), 2));
      double Et = lep_pair.Et();


      store_var("Wgg", Wgg);
      store_var("Wmiss", Wmiss);
      store_var("Emiss", Emiss);
      store_var("Wlep", Wlep);
      store_var("Wgenlep", Wgenlep);
      store_var("mreco", mreco);
      store_var("mreco2", mreco2);
      store_var("pair_mass", pair_mass);
      store_var("extratracks", extratracks);
      store_var("pair_dphi", pair_dphi);
      store_var("kvc_z", kvc_z);
      store_var("slep_dphi", evt.GenSLRPair_dphi);
      store_var("pho_dphi", pho_dphi);
      store_var("slep_aco", slep_aco);
      store_var("pair_aco", pair_aco);
      store_var("pho_aco", pho_aco);
      store_var("closest_extra", closest_extra);
      store_var("closest_hp_extra", closest_hp_extra);
      store_var("pair_eta_diff", pair_eta_diff);
      store_var("eta1", eta1);
      store_var("eta2", eta2);
      store_var("pt1", pt1);
      store_var("pt2", pt2);
      store_var("deltaR", deltaR);
      store_var("ptMiss", ptMiss);
      store_var("ptTot", ptTot);
      store_var("Et", Et);
      store_var("xip", xip);
      store_var("xim", xim);
      store_var("Pt", lep_pair.Pt());
      store_var("Mt", sqrt(lep_pair.E()*lep_pair.E()-lep_pair.Pz()*lep_pair.Pz()));
      if(bgs.find(name) == bgs.end())
        store_var("pair_y", slep_pair.Rapidity());
      else
        store_var("pair_y", lep_pair.Rapidity());
      n_rows++;
//      h_pair_mass.Fill( evt.Pair_mass[j] );
//      h_extratracks.Fill( evt.Pair_extratracks0p5mm[j] );
//      h_pair_dphi.Fill( evt.Pair_dphi[j] );
//      h_kvc_z.Fill( evt.KalmanVertexCand_z[j] );

//      newtree1->Fill();

      //cout << evt.Pair_extratracks2mm[j] << endl;

      /*auto kvc = TVector3( evt.KalmanVertexCand_x[j], evt.KalmanVertexCand_y[j], evt.KalmanVertexCand_z[j] );
      double min_dist = 999.;
      int min_dist_vtx = -1;
      for ( unsigned int k = 0; k < evt.nPrimVertexCand; ++k ) {
        auto pvc = TVector3( evt.PrimVertexCand_x[k], evt.PrimVertexCand_y[k], evt.PrimVertexCand_z[k] );
        const double dist = ( kvc-pvc ).Mag();
        if ( dist < min_dist ) {
          min_dist = dist;
          min_dist_vtx = k;
        }
      }*/
    }
  }

  // Scatter plots
//  TH2F h2("h2", "Wmiss and Wgg", 100, WggA, 100, WmissA);
//  TCanvas c("Wmiss-Wgg");
//  TGraph h2(N, WggA, WmissA);
//  h2.Draw("ap");
//  c.Print("missgg.pdf");

  file.Close();
  TFile file2(filename2, "recreate");
  auto tree2 = new TTree("computed", "Computed quantities from the ROOT analysis");

  for(map<string,int>::iterator it = var_map.begin(); it != var_map.end(); ++it) {
    keypairs.push_back(*it);
  }
  std::sort(keypairs.begin(), keypairs.end(), [](auto &left, auto &right) {
      return left.second < right.second;
  });
  vector<double> reg(n_vars);
  for(int j = 0; j < n_vars; j++) {
    keys.push_back(strdup(keypairs[j].first.c_str()));
    tree2->Branch(keys[j], &reg[j]);
  }

  write_csv(csv_filename);

  start_time();
  for(int i = 0; i < n_rows; i++) {
    for(int j = 0; j < n_vars; j++) {
      reg[j] = var_v[i][j];
    }
    tree2->Fill();
  }
  end_time();

//  newtree1->Write();
  tree2->Write();
  file2.Close();

  // Friend tree
//  tree->AddFriend("computed",filename2);

//  draw_canvas(h_pair_mass, "h_pair_mass", "text 1");
//  draw_canvas(h_pair_dphi, "h_pair_dphi", "text 2");
//  draw_canvas(h_extratracks, "h_extratracks", "text 3");
//  draw_canvas(h_kvc_z, "h_kvc_z", "text 4");
//  draw_canvas(h_susy_mass, "h_susy_mass", "text 5");
}
