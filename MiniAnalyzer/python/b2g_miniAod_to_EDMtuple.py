import FWCore.ParameterSet.Config as cms

process = cms.Process("B2G")

# Input file
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
#'file:/eos/uscms/store/user/bchiari1/miniaod/ver1/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.00.root'
'/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E873348E-BC70-E411-BFA8-0025907B4FD6.root'
))
# Number of events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

# Message Service
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Output file
process.out = cms.OutputModule("PoolOutputModule",
		   fileName = cms.untracked.string("ntuple.root"),
		   SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p')),
		   outputCommands = cms.untracked.vstring('drop *','keep *_TriggerResults_*_HLT','keep *_patTrigger_*_*','keep *_selectedPatTrigger_*_*','keep *_b2g_*_*'))
process.outpath = cms.EndPath(process.out)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.out.dropMetaData = cms.untracked.string("DROPPED")

# Add extra jet collections to miniAOD (ca jets)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
from RecoJets.JetProducers.caTopTaggers_cff import cmsTopTagPFJetsCHS
from RecoJets.JetProducers.caTopTaggers_cff import caTopTagInfos

process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV()>0"))
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
					nFilt = cms.int32(2),
					zcut = cms.double(0.1),
					rcut_factor = cms.double(0.5))
process.cmsTopTagPFJetsCHS = cmsTopTagPFJetsCHS.clone(src = cms.InputTag('chs'))
process.caTopTagInfos = caTopTagInfos.clone(src = cms.InputTag("cmsTopTagPFJetsCHS"))
process.ca8Njettiness = Njettiness.clone(src = 'ca8PFJetsCHS')

process.selectedca8PFJetsCHS = cms.EDFilter('PFJetSelector',
					src = cms.InputTag('ca8PFJetsCHS'),
					cut = cms.string('pt > 25 && abs(eta) < 2.4'))
process.selectedca8Njettiness = Njettiness.clone(src = 'selectedca8PFJetsCHS')

process.extraJetCols = cms.Sequence(
	process.chs *
	process.ca8PFJetsCHS *
	process.ca8PFJetsCHSPruned *
	process.ca8Njettiness *
	process.cmsTopTagPFJetsCHS *
	process.caTopTagInfos *
	process.selectedca8PFJetsCHS *
	process.selectedca8Njettiness
)

# Ntuplizer
process.trigger = cms.EDAnalyzer('b2g_miniAodAnalyzer_trigger',
					bits		= cms.InputTag("TriggerResults"),
	    				prescales 	= cms.InputTag("patTrigger"),
					objects		= cms.InputTag("selectedPatTrigger")
)
process.b2g = cms.EDFilter('b2g_miniAodAnalyzer_general',
					vertexToken	= cms.InputTag("offlineSlimmedPrimaryVertices"),
	    				metToken 	= cms.InputTag("slimmedMETs"),
					electronToken	= cms.InputTag("slimmedElectrons"),
		    			muonToken	= cms.InputTag("slimmedMuons"),
		    			tauToken	= cms.InputTag("slimmedTaus"),
					photonToken	= cms.InputTag("slimmedPhotons"),
					ak4slimmedToken	= cms.InputTag("slimmedJets"),
					ak8slimmedToken	= cms.InputTag("slimmedJetsAK8"),
					pfcandsToken	= cms.InputTag("packedPFCandidates"),
					chscandsToken	= cms.InputTag("chs"),
					ca8Token 	= cms.InputTag("selectedca8PFJetsCHS"),
					ca8tau1Token 	= cms.InputTag("selectedca8Njettiness", "tau1"),
					ca8tau2Token 	= cms.InputTag("selectedca8Njettiness", "tau2"),
					ca8tau3Token 	= cms.InputTag("selectedca8Njettiness", "tau3"),
					ca8prunedToken 	= cms.InputTag("ca8PFJetsCHSPruned"),
					ca8subjetsToken	= cms.InputTag("ca8PFJetsCHSPruned", "SubJets"),
					ca8TopTagToken 	= cms.InputTag("cmsTopTagPFJetsCHS"),
					ca8TopTagInfoToken = cms.InputTag("caTopTagInfos")
)

# Path
process.p = cms.Path(
	process.extraJetCols *
	process.trigger *
	process.b2g
)
print "---+++---+++---+++---"
