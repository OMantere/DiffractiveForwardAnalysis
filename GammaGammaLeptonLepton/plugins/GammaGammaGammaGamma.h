// -*- C++ -*-
//
// Package:    DiffractiveForwardAnalysis/GammaGammaLeptonLepton
// Class:      GammaGammaGammaGamma
// 
/**\class GammaGammaGammaGamma GammaGammaGammaGamma.cc DiffractiveForwardAnalysis/GammaGammaLeptonLepton/plugins/GammaGammaGammaGamma.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Forthomme
//         Created:  Sun, 29 May 2016 19:05:16 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// CT-PPS objects
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/VertexFinder.h"
//
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "TTree.h"
#include "TLorentzVector.h"

#define MAX_PHOTONS 250
#define MAX_PHOTON_PAIRS 50
#define MAX_PROTONS 100
#define MAX_PROTON_PAIRS 10

class GammaGammaGammaGamma : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
 public:
  explicit GammaGammaGammaGamma(const edm::ParameterSet&);
  ~GammaGammaGammaGamma();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  void clearTree();

  bool passPhotonId(const edm::Ptr< pat::Photon >& photon_ref) const;

  edm::EDGetTokenT< edm::View<pat::Photon> > fPhotonToken;
  edm::EDGetTokenT< std::vector<reco::Conversion> > fConversionToken;
  edm::EDGetTokenT< edm::View<pat::Electron> > fElectronToken;
  edm::EDGetTokenT< std::vector<reco::GsfElectron> > fRecoElectronToken;
  edm::EDGetTokenT< edm::View<pat::PackedCandidate> > fPFCandidateToken;
  edm::EDGetTokenT< edm::View<reco::Vertex> > fVertexToken;
  edm::EDGetTokenT< reco::BeamSpot > fBeamSpotToken;

  edm::EDGetTokenT< edm::ValueMap<bool> > fPhotonMediumIdBoolMapToken;
  edm::EDGetTokenT< edm::ValueMap<vid::CutFlowResult> > fPhotonMediumIdFullInfoMapToken;
  edm::EDGetTokenT< edm::ValueMap<float> > fPhotonMVAIdValuesMapToken;
  edm::EDGetTokenT< edm::ValueMap<int> > fPhotonMVAIdCategoriesMapToken;

  edm::EDGetTokenT< edm::DetSetVector<TotemRPLocalTrack> > fProtonToken;
    
  // The first map simply has pass/fail for each particle
  edm::Handle< edm::ValueMap<bool> > fPhotonMediumIdDecisions;
  // The second map has the full info about the cut flow
  //edm::Handle< edm::ValueMap<vid::CutFlowResult> > fPhotonMediumIdCutflowData;

  //
  double fPhotonMinPt;
  bool fFetchProtons;
    
  TTree* fTree;

  //VertexFinder fVertexFinder;
  //
  int aRunId;
  int aLumiSection;
  int aEventNum;

  int aNumPhotons;
  double aPhotonPt[MAX_PHOTONS];
  double aPhotonEta[MAX_PHOTONS];
  double aPhotonPhi[MAX_PHOTONS];
  double aPhotonVtxX[MAX_PHOTONS];
  double aPhotonVtxY[MAX_PHOTONS];
  double aPhotonVtxZ[MAX_PHOTONS];
  double aPhotonIdMVAvalue[MAX_PHOTONS];
  int aPhotonIdMVAcategory[MAX_PHOTONS];
  double aPhotonIsoTrack[MAX_PHOTONS];
  double aPhotonIsoECAL[MAX_PHOTONS];
  double aPhotonIsoHCAL[MAX_PHOTONS];
  double aPhotonIsoCalo[MAX_PHOTONS];
  double aPhotonR9[MAX_PHOTONS];
  int aPhotonElectronVeto[MAX_PHOTONS];
  int aPhotonConverted[MAX_PHOTONS];

  int aNumPhotonPairs;
  double aPhotonPairVertexDist[MAX_PHOTON_PAIRS];
  double aPhotonPairMass[MAX_PHOTON_PAIRS];
  double aPhotonPairPt[MAX_PHOTON_PAIRS];
  double aPhotonPairDpt[MAX_PHOTON_PAIRS];
  double aPhotonPairDphi[MAX_PHOTON_PAIRS];

  int aNumProtons;
  double aProtonX[MAX_PROTONS];
  double aProtonZ[MAX_PROTONS];
  double aProtonY[MAX_PROTONS];
  double aProtonXsigma[MAX_PROTONS];
  double aProtonYsigma[MAX_PROTONS];
  int aProtonArm[MAX_PROTONS];
  int aProtonSide[MAX_PROTONS];
};

//define this as a plug-in
DEFINE_FWK_MODULE(GammaGammaGammaGamma);