import FWCore.ParameterSet.Config as cms
###############################################
# Test V1
###############################################
#testFile =  'file:/eos/uscms/store/user/bchiari1/miniaod/ver1/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.00.root'
###############################################

process = cms.Process("B2G")

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
'/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/000A288B-8C70-E411-8830-20CF300E9EB6.root',
'/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3E123FA3-BC70-E411-8EF4-0025907B4ECA.root',
'/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/72721CA0-BC70-E411-9320-0025907B4F32.root',
'/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/78813892-BC70-E411-8128-00259073E34C.root',
'/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7E657093-BC70-E411-AE05-E0CB4E1A1163.root',
'/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9849B195-BC70-E411-AB45-E0CB4E4408DD.root',
'/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9E280490-BC70-E411-9911-E0CB4E19F9BC.root',
'/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C232BA9C-BC70-E411-A452-20CF3027A5AD.root',
'/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C820D7FD-8B70-E411-97A2-485B39800B62.root',
'/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E873348E-BC70-E411-BFA8-0025907B4FD6.root')
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.out = cms.OutputModule("PoolOutputModule",
		   fileName = cms.untracked.string("ntuple.root"),
		   SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p')),
		   outputCommands = cms.untracked.vstring('drop *','keep *_b2g_*_*'))
process.outpath = cms.EndPath(process.out)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.out.dropMetaData = cms.untracked.string("DROPPED")

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

# Can add a cut if desired
process.selectedca8PFJetsCHS = cms.EDFilter('PFJetSelector',
					src = cms.InputTag('ca8PFJetsCHS'),
					cut = cms.string('pt > 25 && abs(eta) < 2.4'))

process.ca8Njettiness = Njettiness.clone(src = 'ca8PFJetsCHS')
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

process.b2g = cms.EDFilter('b2g_miniAodAnalyzer_general',
					vertexToken	= cms.InputTag('offlineSlimmedPrimaryVertices'),
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

process.p = cms.Path(
	process.extraJetCols *
	process.b2g
)
