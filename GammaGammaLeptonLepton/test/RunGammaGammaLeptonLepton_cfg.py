import FWCore.ParameterSet.Config as cms

runOnMC = True
useAOD = True # AOD or MiniAOD?

process = cms.Process('NoSplit')

#########################
#    General options    #
#########################


#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
# process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#########################
#      Input files      #
#########################

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'/store/data/Run2016G/DoubleEG/AOD/23Sep2016-v1/100000/0042DBD3-BA8E-E611-919E-002481ACDAA8.root',
# 'file:/afs/cern.ch/user/j/jmantere/private/cms/CMSSW_9_2_3/src/%s.root' % in_file,
#   'file:/afs/cern.ch/user/j/jmantere/private/cms/CMSSW_9_4_0/src/final.root',
#   '/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM',
#   '/store/mc/RunIIFall17DRPremix/TTTo2L2Nu_TuneCP5down_PSweights_13TeV-powheg-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/00000/02705333-A03B-E811-8232-1CC1DE1D1E3C.root',
  '/store/mc/RunIIFall17DRPremix/TTTo2L2Nu_TuneCP5up_PSweights_13TeV-powheg-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/560BE6FE-763A-E811-B8D6-A4BF0112BE24.root',
# '/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v1/AODSIM',
#   '/store/mc/RunIIFall17DRPremix/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/90000/FEAA6CA7-8645-E811-90BE-FA163EC11CAA.root',
# 'file:/afs/cern.ch/user/k/karjas/private/CMSSW/dataFold/Events/wwllbg.root'
    ),
    #firstEvent = cms.untracked.uint32(0)
)

#########################
#        Triggers       #
#########################

process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults", "", "RECO")
process.hltFilter.HLTPaths = cms.vstring(
    'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL*'
#    'HLT_DoubleEle33_CaloIdL_MW_v*',
#    'HLT_Ele27_HighEta_Ele20_Mass55_v*',
#    'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*',
)

#########################
#      Preskimming      #
#########################
process.load("Configuration.StandardSequences.GeometryDB_cff") ## FIXME need to ensure that this is the good one
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")


#########################
#     PAT-ification     #
#########################
## Look at https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#Core_Tools for more information

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('PATuple.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_offline*PrimaryVertices*_*_*',
        'keep *_selectedPatMuons*_*_*',
        'keep *_*lectron*_*_*',
        'keep *_selectedPatElectrons*_*_*',
        'keep *_selectedPat*Photons*_*_*',
        'keep *_selectedPatJets*_*_*',
        'keep *_*MET*_*_*',
        'keep *_*particleFlow*_*_*',
    ),
)
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patAlgosToolsTask = getPatAlgosToolsTask(process)

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
patAlgosToolsTask.add(process.patCandidatesTask)

process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
patAlgosToolsTask.add(process.selectedPatCandidatesTask)

from PhysicsTools.PatAlgos.tools.coreTools import runOnData
if not runOnMC:
    runOnData( process )

#########################
#      Electron ID      #
#########################

#from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDElectronIdProducer, setupVIDElectronSelection, setupAllVIDIdsInModule, DataFormat
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

switchOnVIDElectronIdProducer(process, DataFormat.AOD)
#setupAllVIDIdsInModule(process, 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff', setupVIDElectronSelection)
setupAllVIDIdsInModule(process, 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff', setupVIDElectronSelection)

#########################
#       Photon ID       #
#########################

switchOnVIDPhotonIdProducer(process, DataFormat.AOD)
setupAllVIDIdsInModule(process, 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff', setupVIDPhotonSelection)

#########################
#       Analysis        #
#########################

process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.GammaGammaLL_cfi")

process.ggll_aod.triggersList = process.hltFilter.HLTPaths
process.ggll_aod.leptonsType = cms.string('Muon')
# process.ggll_aod.leptonsType = cms.string('ElectronMuon')
# process.ggll_aod.leptonsType = cms.string('Electron')
process.ggll_aod.runOnMC = cms.bool(runOnMC)
process.ggll_aod.fetchProtons = cms.bool(False)

# E/gamma identification
process.ggll_aod.eleIdLabels = cms.PSet(
    mediumLabel = cms.InputTag('mvaEleID-Spring16-GeneralPurpose-V1-wp90'),
       tightLabel = cms.InputTag('mvaEleID-Spring16-GeneralPurpose-V1-wp80'),
)
process.ggll_aod.phoIdLabels = cms.PSet(
    mediumLabel = cms.InputTag('mvaPhoID-Spring16-nonTrig-V1-wp90'),
    tightLabel = cms.InputTag('mvaPhoID-Spring16-nonTrig-V1-wp80'),
)
#process.ggll_aod.eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90")
#process.ggll_aod.eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80")
#process.ggll_aod.phoMediumIdMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp90")
#process.ggll_aod.phoTightIdMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp80")

# prepare the output file
process.TFileService = cms.Service('TFileService',
    # fileName = cms.string('%s.root' % out_file),
    fileName = cms.string('ggll.root'),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(
    # process.hltFilter*
    process.egmPhotonIDSequence*
    process.egmGsfElectronIDSequence*
    process.ggll_aod
)

#process.outpath = cms.EndPath(process.out, patAlgosToolsTask)
process.outpath = cms.EndPath(patAlgosToolsTask)

