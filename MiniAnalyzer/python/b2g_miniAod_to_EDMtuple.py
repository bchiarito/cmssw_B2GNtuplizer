import FWCore.ParameterSet.Config as cms

process = cms.Process("B2G")

# Input file
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
#'file:/eos/uscms/store/user/bchiari1/miniaod/ver1/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.00.root'
'/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E873348E-BC70-E411-BFA8-0025907B4FD6.root'
))
# Number of events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# Message Service
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

# Output file
process.out = cms.OutputModule("PoolOutputModule",
		   fileName = cms.untracked.string("ntuple.root"),
		   SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p')),
		   outputCommands = cms.untracked.vstring('drop *',
		     # Below passes on trigger collections from miniAOD
                     #'keep *_TriggerResults_*_HLT','keep *_patTrigger_*_*','keep *_selectedPatTrigger_*_*',
                     'keep *_*_*_B2G'
))
process.outpath = cms.EndPath(process.out)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.out.dropMetaData = cms.untracked.string("DROPPED")

# Add extra jet collections to miniAOD (ca jets)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
from RecoJets.JetProducers.caTopTaggers_cff import cmsTopTagPFJetsCHS
from RecoJets.JetProducers.caTopTaggers_cff import caTopTagInfos

process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV()>2"))
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
process.trigger = cms.EDFilter('b2g_miniAodAnalyzer_trigger',
					printAll	= cms.bool(False), # Dump all trigger information
					isData		= cms.bool(True), # If false will not save HLT Objects
					bits		= cms.InputTag("TriggerResults", "", "HLT"),
	    				prescales 	= cms.InputTag("patTrigger"),
					objects		= cms.InputTag("selectedPatTrigger"),
					useTriggerList  = cms.bool(False), # Only save information for triggers in the list
					triggerList	= cms.vstring("HLT_Ele25WP60_SC4_Mass55_v1",
								      "HLT_Mu40_v1",
								     )
)
process.general = cms.EDFilter('b2g_miniAodAnalyzer_general',
					vertecies	= cms.InputTag("offlineSlimmedPrimaryVertices"),
	    				met 		= cms.InputTag("slimmedMETs"),
					electrons	= cms.InputTag("slimmedElectrons"),
		    			muons		= cms.InputTag("slimmedMuons"),
		    			taus		= cms.InputTag("slimmedTaus"),
					photons		= cms.InputTag("slimmedPhotons"),
					ak4slimmed	= cms.InputTag("slimmedJets"),
					ak8slimmed	= cms.InputTag("slimmedJetsAK8"),
					ak8grommedMasses = cms.bool(True),
					pfcands		= cms.InputTag("packedPFCandidates"),
					chscands	= cms.InputTag("chs")
)
# Below module updated in future version, does not run
process.jets = cms.EDFilter('b2g_miniAodAnalyzer_jets',
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
	process.general
#	process.jets
)
print "---+++---+++---+++---"
