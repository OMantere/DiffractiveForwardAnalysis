import FWCore.ParameterSet.Config as cms

ggll_aod = cms.EDAnalyzer(
    'GammaGammaLL',
    SqrtS = cms.double(8000.),
    HLTMenuLabel = cms.string('HLT'),
    TriggerResults = cms.InputTag('TriggerResults', '', 'HLT'),
    LeptonsType = cms.InputTag('electron', 'muon'),
    maxExtraTracks = cms.untracked.uint32(10000),
    isoValInputTags = cms.VInputTag(
        cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
        cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
        cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')
    ),
    pileupInfo = cms.InputTag('addPileupInfo'),
    beamSpotInputTag = cms.InputTag('offlineBeamSpot'),
    RecoVertexLabel = cms.InputTag('offlinePrimaryVertices'),
    conversionsInputTag = cms.InputTag('allConversions'),
    #eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
    #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
    #eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    #eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
    #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
    #eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
    rhoIsoInputTag = cms.InputTag('kt6PFJetsForIsolation', 'rho'),
    JetCollectionLabel = cms.InputTag('selectedPatJetsPFlow'),
    MetLabel = cms.InputTag('pfMet'),
    RunOnMC = cms.untracked.bool(True),
    MCAcceptPtCut = cms.untracked.double(0.),
    MCAcceptEtaCut = cms.untracked.double(-1.),
    GenParticlesCollectionLabel = cms.InputTag('genParticles'),
    PrintCandidates = cms.untracked.bool(False),
)

ggll = ggll_aod.clone() ## for backward-compatibility

ggll_miniaod = ggll_aod.clone(
    RecoVertexLabel = cms.InputTag('offlineSlimmedPrimaryVertices'),
    GlobalMuonCollectionLabel = cms.untracked.InputTag("slimmedMuons"), # PAT
    GlobalEleCollectionLabel = cms.untracked.InputTag("slimmedElectrons"), # PAT
    JetCollectionLabel = cms.InputTag('slimmedJetsAK8'),
    MetLabel = cms.InputTag('slimmedMETs'), # PAT
    GenParticlesCollectionLabel = cms.InputTag('prunedGenParticles'),
    PFLabel = cms.untracked.InputTag('packedPFCandidates'),
    pileupInfo = cms.InputTag('slimmedAddPileupInfo'),
)

