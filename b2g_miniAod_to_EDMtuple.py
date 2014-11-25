import FWCore.ParameterSet.Config as cms
###############################################
# Test V1
###############################################
testFile =  'file:./miniaod_plus.root'
###############################################

process = cms.Process("B2G") # <-- Start Process here

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(testFile),
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
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
					ca8TopTagInfoToken = cms.InputTag("caTopTagInfos"),
				    )

process.p = cms.Path(	
			process.b2g
		    )
process.out = cms.OutputModule("PoolOutputModule",
		   fileName = cms.untracked.string("diffmo_susy4body_allPF.root"),
		   SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p')),
		   outputCommands = cms.untracked.vstring('drop *','keep *_b2g*_*_*'))
process.outpath = cms.EndPath(process.out)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.out.dropMetaData = cms.untracked.string("DROPPED")
