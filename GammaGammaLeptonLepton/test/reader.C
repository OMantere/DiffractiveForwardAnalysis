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

void reader( const char* filename = "../ggll.root" )
{
  auto file = TFile::Open( filename );
  auto tree = dynamic_cast<TTree*>( file->Get( "ggll_aod/ntp1" ) );
  ggll::AnalysisEvent evt;
  evt.load( tree, ggll::DiMuon, true );

//  TH1D h_pair_mass( "h_pair_mass", "m(#mu^{+}#mu^{-})\\Events\\GeV?.2f", 200, 0., 1000. );
//  TH1D h_pair_dphi( "h_pair_dphi", "m(#mu^{+}#mu^{-})", 200, -1. * M_PI, 1. * M_PI );
//  TH1D h_extratracks( "h_extratracks", "m(#mu^{+}#mu^{-})", 200, 0.,  10);
//  TH1D h_kvc_z( "h_kvc_z", "m(#mu^{+}#mu^{-})", 200, 0., 15. );
//  TH1D h_susy_mass( "h_susy_mass", "m(#mu^{+}#mu^{-})", 200, 0., 1000. );

  double xim;
  double xip;
  double neutralino_mass;
  double minpairmass;

//  TTree *tree = new TTree("analysis", "The analysis tree");
//  tree->Branch("xip", &xip);
//  tree->Branch("xim",&xim);
//  tree->Branch("neutralino_mass",&neutralino_mass);

  cout << "Entries: " << tree->GetEntriesFast() << endl;
  for ( unsigned long long i = 0; i < tree->GetEntriesFast(); ++i ) {
    tree->GetEntry( i );
    //cout << ">>> event " << i << endl;
    /*cout << "- fired triggers:" << endl;
    for ( const auto& hlt : *evt.HLT_Name ) {
      cout << "  *) " << hlt << endl;
    }*/

    TLorentzVector p1, p2;
    TLorentzVector com(0, 0, 0, 13.e3);
    cout << "Proton count in this event: " << evt.nGenProtCand << endl;
    p1.SetPtEtaPhiE(evt.GenProtCand_pt[0], evt.GenProtCand_eta[0], evt.GenProtCand_phi[0], evt.GenProtCand_e[0]);
    p2.SetPtEtaPhiE( evt.GenProtCand_pt[1], evt.GenProtCand_eta[1], evt.GenProtCand_phi[1], evt.GenProtCand_e[1] );
    cout << "Proton 1 E" << evt.GenProtCand_e[0] << endl;
    cout << "New proton pair mass: " << (p1+p2).M() << endl;

    //cout << "- dilepton pairs:" << endl;
    for ( unsigned int j = 0; j < evt.nPair; ++j ) {
      //cout << "  *) pair " << j << " has invariant mass = " << evt.Pair_mass[j] << endl;

      //if ( evt.Pair_extratracks0p5mm[j] != 0 ) continue;
      //if ( fabs( evt.KalmanVertexCand_z[j] ) > 15. ) continue;
      //if ( 1.-fabs( evt.Pair_dphi[j] )/M_PI > 0.009 ) continue;
      //if ( evt.Pair_mass[j] < 110. ) continue;

      cout << "extratracks: " << evt.Pair_extratracks0p5mm[j] << endl;
      cout << "KalmanVertexCand_z: " << fabs( evt.KalmanVertexCand_z[j] ) << endl;
      cout << "Pair_dphi: " << 1.-fabs( evt.Pair_dphi[j] )/M_PI << endl;
      cout << "Pair_mass: " << evt.Pair_mass[j] << endl;
      cout << "CANDIDATE!!!" << endl;

      const unsigned int l1 = evt.Pair_lepton1[j], l2 = evt.Pair_lepton2[j];
      xip = ( evt.MuonCand_pt[l1]*exp( +evt.MuonCand_eta[l1] ) + evt.MuonCand_pt[l2]*exp( +evt.MuonCand_eta[l2] ) ) / 13.e3;
      xim = ( evt.MuonCand_pt[l1]*exp( -evt.MuonCand_eta[l1] ) + evt.MuonCand_pt[l2]*exp( -evt.MuonCand_eta[l2] ) ) / 13.e3;
      cout << "central system xip/xim = " << xip << " / " << xim << endl;

      h_pair_mass.Fill( evt.Pair_mass[j] );
      h_extratracks.Fill( evt.Pair_extratracks0p5mm[j] );
      h_pair_dphi.Fill( evt.Pair_dphi[j] );
      h_kvc_z.Fill( evt.KalmanVertexCand_z[j] );

      TLorentzVector lep1, lep2;
      lep1.SetPtEtaPhiE(evt.GenMuonCand_pt[l1], evt.GenMuonCand_eta[l1], evt.GenMuonCand_phi[l1], evt.GenMuonCand_e[l1]);
      lep2.SetPtEtaPhiE( evt.GenMuonCand_pt[l2], evt.GenMuonCand_eta[l2], evt.GenMuonCand_phi[l2], evt.GenMuonCand_e[l2] );

      double susymass = (com - (p1+p2) - (lep1+lep2)).M();
      h_susy_mass.Fill( susymass );

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
  draw_canvas(h_pair_mass, "h_pair_mass", "text 1");
  draw_canvas(h_pair_dphi, "h_pair_dphi", "text 2");
  draw_canvas(h_extratracks, "h_extratracks", "text 3");
  draw_canvas(h_kvc_z, "h_kvc_z", "text 4");
  draw_canvas(h_susy_mass, "h_susy_mass", "text 5");
}
