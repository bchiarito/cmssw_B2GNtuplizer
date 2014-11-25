import FWCore.ParameterSet.Config as cms

process = cms.Process("JETCOL")
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'file:/eos/uscms/store/user/bchiari1/miniaod/ver1/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.00.root')
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.OUT = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('miniaod_plus.root'),
    outputCommands = cms.untracked.vstring(['keep *'])
)
process.endpath = cms.EndPath(process.OUT)

from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
from RecoJets.JetProducers.caTopTaggers_cff import cmsTopTagPFJetsCHS
from RecoJets.JetProducers.caTopTaggers_cff import caTopTagInfos

process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
process.ca8PFJetsCHS = ak5PFJets.clone(
					src = 'chs',
					rParam = cms.double(0.8),
					jetPtMin = cms.double(3.0),
					jetAlgorithm = cms.string("CambridgeAachen"))
process.ca8PFJetsCHSPruned = ak5PFJets.clone(
					src = 'chs',
					rParam = cms.double(0.8),
					jetPtMin = cms.double(3.0),
					jetAlgorithm = cms.string("CambridgeAachen"),
					usePruning = cms.bool(True),
					useExplicitGhosts = cms.bool(True),
					jetCollInstanceName = cms.string("SubJets"),
					writeCompound = cms.bool(True),
					nFilt = cms.int32(2), # number of subjets to decluster the fat jet into
					zcut = cms.double(0.1),
					rcut_factor = cms.double(0.5))

process.cmsTopTagPFJetsCHS = cmsTopTagPFJetsCHS.clone(src = cms.InputTag('chs'))
process.caTopTagInfos = caTopTagInfos.clone(src = cms.InputTag("cmsTopTagPFJetsCHS"))

process.selectedca8PFJetsCHS = cms.EDFilter('PFJetSelector',
					src = cms.InputTag('ca8PFJetsCHS'),
					cut = cms.string('pt > 25 && abs(eta) < 2.4'))

process.ca8Njettiness = Njettiness.clone(src = 'ca8PFJetsCHS')
process.selectedca8Njettiness = Njettiness.clone(src = 'selectedca8PFJetsCHS')

process.options = cms.untracked.PSet( 
        wantSummary = cms.untracked.bool(True), 
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.extraJetCols = cms.Path(
	process.chs *
	process.ca8PFJetsCHS *
	process.ca8PFJetsCHSPruned *
	process.ca8Njettiness *
	process.cmsTopTagPFJetsCHS *
	process.caTopTagInfos *
	process.selectedca8PFJetsCHS *
	process.selectedca8Njettiness
)	
