#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "Canvas.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "mt2.h"

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

double rapidity(const TLorentzVector v){
  double rap = 1.0/2*log((v.E()+v.Pz())/(v.E()-v.Pz()));
  return rap;
}

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
int max_events = 7500000;
int max_protons = 100000;
int max_vars = 75;
int n_vars = 0;
int n_rows = 0;
double e_mass, mumass, Chi10mass;

string type;

default_random_engine generator;
uniform_int_distribution<int> proton_random;

vector <vector<double>> var_v(max_events, vector<double>(max_vars));
std::map<string, int> var_map;
vector <std::pair<string, int>> keypairs;
vector<char *> keys;
vector <TLorentzVector> l_protons(max_protons);
vector <TLorentzVector> r_protons(max_protons);
int total_acc = 0;
int total_pairs = 0;
double low_xi_acc = 6500 * (1 - 0.15);
double high_xi_acc = 6500 * (1 - 0.02);
void store_var(string name, double val) {
  if(var_map.count(name) == 0) {
    var_map.insert(std::pair<string,int>(name,n_vars));
    n_vars++;
  }
  var_v[n_rows][var_map.find(name)->second] = val;
}


void read_protons(const char *proton_filename = "acc_events.csv") {
    ifstream ifs(proton_filename);
    string line;
    long long i = -1;
    double pt, eta, phi, e;
    bool left_proton = false;
    bool right_proton = false;
    while (ifs >> line) {
        if (i == max_protons) break;
        if (line != "event" && (!left_proton || !right_proton)) {
            pt = stod(line);
            ifs >> eta;
            ifs >> phi;
            ifs >> e;
            if (eta < 0) {
                l_protons[i].SetPtEtaPhiE(pt, eta, phi, e);
                left_proton = true;
            } else {
                r_protons[i].SetPtEtaPhiE(pt, eta, phi, e);
                right_proton = true;
            }
        } else if (line == "event") {
            left_proton = right_proton = false;
            i++;
        }
    }
    cout << "Added " << (i + 1) << " proton events from file " << proton_filename << endl;
    proton_random = uniform_int_distribution<int>(0, i);
    ifs.close();
}


double MT2(TLorentzVector l1, TLorentzVector l2) {
  double visible_mass;
  if(type == "Muon")
      visible_mass = mumass;
  else
      visible_mass = e_mass;
  TLorentzVector pair = l1+l2;
  asymm_mt2_lester_bisect::disableCopyrightMessage();
  return asymm_mt2_lester_bisect::get_mT2(
      visible_mass, l1.Px(), l1.Py(),
      visible_mass, l2.Px(), l2.Py(),
      -pair.Px(), -pair.Py(),
      Chi10mass, Chi10mass);
}


pair <TLorentzVector, TLorentzVector> get_diffraction_protons() {
    int i = proton_random(generator);
    return make_pair(l_protons[i], r_protons[i]);
}



void write_csv(const char *csv_filename) {
    ofstream ofs(csv_filename, ofstream::out); 
    for (int j = 0; j < n_vars; j++) {
        ofs << keys[j] << delimiter << "\n";
  }
  for(int i = 0; i < n_rows; i++) {
    for(int j = 0; j < n_vars; j++) {
        ofs << var_v[i][j] << delimiter << "\n";
    }
  }
  ofs.close();
}

void reader(const char* c_name = "elastic", const char* fout = "computed")
{
  string name(c_name);
  string cfilename = name;
  string name2(fout);
  string cfilename2 = name2;
  string ccsv_filename= name+".root";
  const char* filename = cfilename.c_str();
  const char* filename2 = cfilename2.c_str();
  const char* csv_filename=ccsv_filename.c_str();
  TFile file( filename);
  auto tree = dynamic_cast<TTree*>( file.Get( "ggll_aod/ntp1") );
  int N = tree->GetEntriesFast();
  ggll::AnalysisEvent evt;

  string type = "Muon";
  if(type == "Muon")
    evt.load( tree, ggll::DiMuon, true );
  else
    evt.load( tree, ggll::DiElectron, true );

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
 
  

  cout << "Entries: " << N << endl;
  start_time();
  if (false)
      read_protons();  
  if(N>500000) N = 500000;
 

  for (unsigned long long i = 0; i < N; i++){
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
  if (false) {
      pair <TLorentzVector, TLorentzVector> protons = get_diffraction_protons();
      p2 = protons.first;
      p1 = protons.second;
  } else {
        p1.SetPtEtaPhiE(evt.GenProtCand_pt[0], evt.GenProtCand_eta[0], evt.GenProtCand_phi[0], evt.GenProtCand_e[0]);
        p2.SetPtEtaPhiE(evt.GenProtCand_pt[1], evt.GenProtCand_eta[1], evt.GenProtCand_phi[1], evt.GenProtCand_e[1]);
    }

    for ( unsigned int j = 0; j < evt.nPair; ++j ) {
      if ( fabs( evt.KalmanVertexCand_z[j] ) > 15. ) continue;
      //cout << j << endl;

      const unsigned int l1 = evt.Pair_lepton1[j], l2 = evt.Pair_lepton2[j];
      double El1, El2;
      TLorentzVector pl1, pl2, pg1, pg2, pl1g, pl2g;
      if(type == "Muon") {
        El1 = evt.MuonCand_e[l1];
        El2 = evt.MuonCand_e[l2];
        xip = ( evt.MuonCand_pt[l1]*exp( +evt.MuonCand_eta[l1] ) + evt.MuonCand_pt[l2]*exp( +evt.MuonCand_eta[l2] ) ) / 13.e3;
        xim = ( evt.MuonCand_pt[l1]*exp( -evt.MuonCand_eta[l1] ) + evt.MuonCand_pt[l2]*exp( -evt.MuonCand_eta[l2] ) ) / 13.e3;
        pl1g.SetPtEtaPhiE(evt.GenMuonCand_pt[0], evt.GenMuonCand_eta[0], evt.GenMuonCand_phi[0], evt.GenMuonCand_e[0]);
        pl2g.SetPtEtaPhiE(evt.GenMuonCand_pt[1], evt.GenMuonCand_eta[1], evt.GenMuonCand_phi[1], evt.GenMuonCand_e[1]);
        pl1.SetPtEtaPhiE(evt.MuonCand_pt[l1], evt.MuonCand_eta[l1], evt.MuonCand_phi[l1], evt.MuonCand_e[l1]);
        pl2.SetPtEtaPhiE(evt.MuonCand_pt[l2], evt.MuonCand_eta[l2], evt.MuonCand_phi[l2], evt.MuonCand_e[l2]);
      } else {
        El1 = evt.EleCand_e[l1];
        El2 = evt.EleCand_e[l2];
       // xip = ( evt.EleCand_pt[l1]*exp( +evt.EleCand_eta[l1] ) + evt.EleCand_pt[l2]*exp( +evt.EleCand_eta[l2] ) ) / 13.e3;
       // xim = ( evt.EleCand_pt[l1]*exp( -evt.EleCand_eta[l1] ) + evt.EleCand_pt[l2]*exp( -evt.EleCand_eta[l2] ) ) / 13.e3;
        pl1g.SetPtEtaPhiE(evt.GenEleCand_pt[0], evt.GenEleCand_eta[0], evt.GenEleCand_phi[0], evt.GenEleCand_e[0]);
        pl2g.SetPtEtaPhiE(evt.GenEleCand_pt[1], evt.GenEleCand_eta[1], evt.GenEleCand_phi[1], evt.GenEleCand_e[1]);
        pl1.SetPtEtaPhiE(evt.EleCand_pt[l1], evt.EleCand_eta[l1], evt.EleCand_phi[l1], evt.EleCand_e[l1]);
        pl2.SetPtEtaPhiE(evt.EleCand_pt[l2], evt.EleCand_eta[l2], evt.EleCand_phi[l2], evt.EleCand_e[l2]);
      }

      // We don't have working photons for now
//      pg1.SetPtEtaPhiE(evt.GenPhotCand_pt[0], evt.GenPhotCand_eta[0], evt.GenPhotCand_phi[0], evt.GenPhotCand_e[0]);
//      pg2.SetPtEtaPhiE(evt.GenPhotCand_pt[1], evt.GenPhotCand_eta[1], evt.GenPhotCand_phi[1], evt.GenPhotCand_e[1]);
      double xipProt = (P1.Pz() - p1.Pz())/P1.Pz();
      double ximProt = (P2.Pz() - p2.Pz())/P2.Pz();
      pg1 = P1 - p1;
      pg2 = P2 - p2;
      Wgg = 2 * sqrt(pg1.E()*pg2.E());
      const TLorentzVector lep_pair = pl1 + pl2;
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


      int nPrimVertex = evt.nPrimVertexCand;
//      cout<<nPrimVertex<<endl;

      double dd = 10000;
      
      int nTracks2 = 0;

      for(int k = 0; k<nPrimVertex; k++){
          double dx = evt.PrimVertexCand_x[k] - evt.KalmanVertexCand_x[k];
          double dy = evt.PrimVertexCand_y[k] - evt.KalmanVertexCand_y[k];
          double dz = evt.PrimVertexCand_z[k] - evt.KalmanVertexCand_z[k];

          double distance = sqrt(dx*dx + dy*dy + dz*dz);

          if(dd > distance){
              dd = distance;
              nTracks2 = evt.PrimVertexCand_tracks[k];
          }

      }

      //cout << dd << endl;

      store_var("nPrimTracks", nTracks2);
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
      store_var("ptMiss", ptMiss);
      store_var("ptTot", ptTot);
      store_var("xip", xip);
      store_var("xim", xim);
      store_var("xipProt", xipProt);
      store_var("ximProt", ximProt);

      store_var("minDist", dd);

      double deltaR = sqrt(pow(evt.MuonCand_eta[l2] - evt.MuonCand_eta[l1], 2) + pow(evt.MuonCand_phi[l2] - evt.MuonCand_phi[l1], 2));

      double deltaEta = evt.MuonCand_eta[l2] - evt.MuonCand_eta[l1];

      double ptl1l2 = sqrt(lep_pair.Px()*lep_pair.Px() + lep_pair.Py()*lep_pair.Py());
      double ptl1l2Check = lep_pair.Pt();
      double ptl1 = sqrt(pl1.Px()*pl1.Px()+pl1.Py()*pl1.Py());
      double ptl2 = sqrt(pl2.Px()*pl2.Px()+pl2.Py()*pl2.Py());

      double etal1l2 = lep_pair.Eta();
      double eta1 = evt.MuonCand_eta[l1];
      double eta2 = evt.MuonCand_eta[l2];


      double rapidityl1l2 = rapidity(lep_pair);
      double rapidityl1 = rapidity(pl1);
      double rapidityl2 = rapidity(pl2);
      double rapiditysl = rapidity(pg1+pg2);

      double l1Phi = pl1.Phi();
      double l2Phi = pl2.Phi();

      double a = 1 - abs(pair_dphi)/M_PI; // a = pair_aco
      double slMass = (pg1 + pg2).M();

      double detA = (pl1.Px()*pl2.Py() - pl2.Px()*pl1.Py());

      double xi1 = 1/detA*(lep_pair.Px()*pl2.Py()-lep_pair.Py()*pl2.Px());//Atlas paper, page 10
      double xi2 = 1/detA*(lep_pair.Py()*pl1.Px()-lep_pair.Px()*pl1.Py());

      double mtaus2 = 2*pl1*pl2*(1+xi1)*(1+xi2);
      double xi11 = xi1-1;
      double xi21 = xi2-1;

      double xiPrime = lep_pair.Px()*pl2.Py()-lep_pair.Py()*pl1.Px();

      double xi3 = 1/detA*(lep_pair.Px()*pl2.Py()-lep_pair.Py()*pl1.Px());
      double xi4 = 1/detA*(lep_pair.Px()*pl1.Py()-lep_pair.Py()*pl2.Px());
      double xi5 = lep_pair.Px()*pl1.Py()-lep_pair.Py()*pl1.Px();


      //cout << xi1 << ", " << xi2 << endl;
      store_var("deltaR", deltaR);
      store_var("deltaEta", deltaEta);
      store_var("ptl1l2", ptl1l2);
      store_var("ptl1", ptl1);
      store_var("ptl2", ptl2);
      store_var("etal1l2", etal1l2);
      store_var("eta1", eta1);
      store_var("eta2", eta2);
      store_var("rapidityl1l2", rapidityl1l2);
      store_var("rapidityl1", rapidityl1);
      store_var("rapidityl2", rapidityl2);
      store_var("rapiditysl", rapiditysl);
      store_var("a", a);
      store_var("slMass",  slMass);
      store_var("xi1", xi1);
      store_var("xi2", xi2);
      store_var("xi3", xi3);
      store_var("xi4", xi4);
      store_var("xi5", xi5);
      store_var("xi11", xi11);
      store_var("xi21", xi21);
      store_var("invdetA", 1/detA);
      store_var("xiPrime", xiPrime);
      store_var("mtaus2", mtaus2);
      store_var("l1Phi", l1Phi);
      store_var("l2Phi", l2Phi);

    // Isolation parameters
      double trackiso = evt.MuonCand_trackiso[j];
      double ecaliso = evt.MuonCand_ecaliso[j];
      double hcaliso = evt.MuonCand_hcaliso[j];

      //cout<<trackiso<<endl;

      store_var("trackiso", trackiso);
      store_var("ecaliso", ecaliso);
      store_var("hcaliso", hcaliso);


      double mt2 = MT2(pl1, pl2);

      store_var("mt2", mt2);


      n_rows++;
    }
  }

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

  tree2->Write();
  file2.Close();
}
