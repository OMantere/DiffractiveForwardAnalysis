#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AnalysisEvent.h"

#include "Canvas.h"
#include "TH1.h"
#include "TLorentzVector.h"

#include <iostream>

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

double rapidity(const TLorentzVector v){
    double rap = 1.0/2*log((v.E()+v.Pz())/(v.E()-v.Pz()));
    return rap;
}

void reader( const char* filename = "/afs/cern.ch/user/k/karjas/private/CMSSW/dataFold/GammaGammaOutput/wwllbgKristianoutput.root", const char* filename2 = "computedwwllbgKristian.root")
{
  TFile file( filename);
  auto tree = dynamic_cast<TTree*>( file.Get( "ggll_aod/ntp1") );
  int N = tree->GetEntriesFast();
  ggll::AnalysisEvent evt;
  evt.load( tree, ggll::DiMuon, true );

  TFile file2(filename2, "recreate");
  auto tree2 = new TTree("computed", "Computed quantities from the ROOT analysis");
  auto newtree1 = new TTree("ggll", "The AnalysisEvent from GGLL");
  evt.attach(newtree1, ggll::DiMuon, true);

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

 /* TTree *tree = new TTree("analysis", "The analysis tree");
  tree->Branch("xip", &xip);
  tree->Branch("xim",&xim);
  tree->Branch("neutralino_mass",&neutralino_mass);

  */
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
 //   cout << "Proton count in this event: " << evt.nGenProtCand << endl;
    p1.SetPtEtaPhiE(evt.GenProtCand_pt[0], evt.GenProtCand_eta[0], evt.GenProtCand_phi[0], evt.GenProtCand_e[0]);
    p2.SetPtEtaPhiE(evt.GenProtCand_pt[1], evt.GenProtCand_eta[1], evt.GenProtCand_phi[1], evt.GenProtCand_e[1]);
 //   cout << "Proton 1 pz: " << p1.Z() << endl;
 //   cout << "Proton 2 pz: " << p2.Z() << endl;
 //   cout << "New proton pair mass: " << (p1+p2).M() << endl;

    //cout << "- dilepton pairs:" << endl;
    for ( unsigned int j = 0; j < evt.nPair; ++j ) {
      //cout << "  *) pair " << j << " has invariant mass = " << evt.Pair_mass[j] << endl;

      //if ( evt.Pair_extratracks0p5mm[j] != 0 ) continue;
      if ( fabs( evt.KalmanVertexCand_z[j] ) > 15. ) continue;
      //if ( 1.-fabs( evt.Pair_dphi[j] )/M_PI > 0.009 ) continue;
      //if ( evt.Pair_mass[j] < 110. ) continue;


      const unsigned int l1 = evt.Pair_lepton1[j], l2 = evt.Pair_lepton2[j];
      const double El1 = evt.MuonCand_e[l1];
      const double El2 = evt.MuonCand_e[l2];

      xip = (P1.Pz() - p1.Pz())/P1.Pz();

      //xip = ( evt.MuonCand_pt[l1]*exp( +evt.MuonCand_eta[l1] ) + evt.MuonCand_pt[l2]*exp( +evt.MuonCand_eta[l2] ) ) / 13.e3;
      xim = ( evt.MuonCand_pt[l1]*exp( -evt.MuonCand_eta[l1] ) + evt.MuonCand_pt[l2]*exp( -evt.MuonCand_eta[l2] ) ) / 13.e3;
  //    cout << "central system xip/xim = " << xip << " / " << xim << endl;

      TLorentzVector pl1, pl2, pg1, pg2, pl1g, pl2g;
      pl1g.SetPtEtaPhiE(evt.GenMuonCand_pt[0], evt.GenMuonCand_eta[0], evt.GenMuonCand_phi[0], evt.GenMuonCand_e[0]);
      pl2g.SetPtEtaPhiE(evt.GenMuonCand_pt[1], evt.GenMuonCand_eta[1], evt.GenMuonCand_phi[1], evt.GenMuonCand_e[1]);
      pl1.SetPtEtaPhiE(evt.MuonCand_pt[l1], evt.MuonCand_eta[l1], evt.MuonCand_phi[l1], evt.MuonCand_e[l1]);
      pl2.SetPtEtaPhiE(evt.MuonCand_pt[l2], evt.MuonCand_eta[l2], evt.MuonCand_phi[l2], evt.MuonCand_e[l2]);
      double pair_eta_diff = fabs(evt.MuonCand_eta[l1] - evt.MuonCand_eta[l2]);
      // We don't have working photons for now
//      pg1.SetPtEtaPhiE(evt.GenPhotCand_pt[0], evt.GenPhotCand_eta[0], evt.GenPhotCand_phi[0], evt.GenPhotCand_e[0]);
//      pg2.SetPtEtaPhiE(evt.GenPhotCand_pt[1], evt.GenPhotCand_eta[1], evt.GenPhotCand_phi[1], evt.GenPhotCand_e[1]);
      pg1 = P1 - p1;
      pg2 = P2 - p2;
     // cout << "Gamma 2 pZ: " << pg2.Z() << endl;
     // cout << "Gamma 2 E: " << pg2.E() << endl;
     // cout << "GenLep 1 E: " << pl1g.E() << endl;
     // cout << "GenLep 2 E: " << pl2g.E() << endl;
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
      //cout << "ptMiss: " << ptMiss << endl;
      
      double deltaR = sqrt(pow(evt.MuonCand_eta[l2] - evt.MuonCand_eta[l1], 2) + pow(evt.MuonCand_phi[l2] - evt.MuonCand_phi[l1], 2)); 


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
     
      double a = 1 - pair_dphi/M_PI;

      double slMass = (pg1 + pg2).M();

      tree2->Branch("Wgg", &Wgg);
      tree2->Branch("Wmiss", &Wmiss);
      tree2->Branch("Emiss", &Emiss);
      tree2->Branch("Wlep", &Wlep);
      tree2->Branch("Wgenlep", &Wgenlep);
      tree2->Branch("mreco", &mreco); // Reconstructed mass
      tree2->Branch("mreco2", &mreco2); // Reconstructed mass of the pair
      tree2->Branch("pair_mass", &pair_mass); // Mass of the pair
      tree2->Branch("extratracks", &extratracks); // Number of tracks within 5mm of the main jet
      tree2->Branch("pair_dphi", &pair_dphi); // Delta phi of the 
      tree2->Branch("kvc_z", &kvc_z);
      tree2->Branch("slep_dphi", &evt.GenSLRPair_dphi);
      tree2->Branch("pho_dphi", &pho_dphi);
      tree2->Branch("slep_aco", &slep_aco);
      tree2->Branch("pair_aco", &pair_aco);
      tree2->Branch("pho_aco", &pho_aco);
      tree2->Branch("closest_extra", &closest_extra);
      tree2->Branch("closest_hp_extra", &closest_hp_extra);
      tree2->Branch("pair_eta_diff", &pair_eta_diff);
      tree2->Branch("deltaR", &deltaR);
      tree2->Branch("ptMiss", &ptMiss);
      tree2->Branch("ptTot", &ptTot);
      tree2->Branch("ptl1l2", &ptl1l2);
      tree2->Branch("ptl1l2Check", &ptl1l2Check);
      tree2->Branch("ptl1", &ptl1);
      tree2->Branch("ptl2", &ptl2);
      tree2->Branch("xip", &xip);
      tree2->Branch("xim", &xim);
    
      tree2->Branch("etal1l2",&etal1l2);
      tree2->Branch("eta1", &eta1);
      tree2->Branch("eta2", &eta2);
      tree2->Branch("rapidityl1l2",&rapidityl1l2);
      tree2->Branch("rapidityl1",&rapidityl1);
      tree2->Branch("rapidityl2",&rapidityl2);
      tree2->Branch("rapiditysl",&rapiditysl);
      tree2->Branch("a",&a);
      tree2->Branch("slMass", &slMass);
      //      h_pair_mass.Fill( evt.Pair_mass[j] );
//      h_extratracks.Fill( evt.Pair_extratracks0p5mm[j] );
//      h_pair_dphi.Fill( evt.Pair_dphi[j] );
//      h_kvc_z.Fill( evt.KalmanVertexCand_z[j] );

      tree2->Fill();
      newtree1->Fill();

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


  file2.cd();
  newtree1->Write();
  tree2->Write();
  file.Close();
  file2.Close();

  // Friend tree
//  tree->AddFriend("computed",filename2);

//  draw_canvas(h_pair_mass, "h_pair_mass", "text 1");
//  draw_canvas(h_pair_dphi, "h_pair_dphi", "text 2");
//  draw_canvas(h_extratracks, "h_extratracks", "text 3");
//  draw_canvas(h_kvc_z, "h_kvc_z", "text 4");
//  draw_canvas(h_susy_mass, "h_susy_mass", "text 5");
}
