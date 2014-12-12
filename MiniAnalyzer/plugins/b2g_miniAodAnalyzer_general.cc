#include <iostream>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/JetReco/interface/CATopJetTagInfo.h"

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include "fastjet/tools/Filter.hh"
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>

namespace B2gMiniToEdm
{

  	typedef std::vector<reco::Candidate::PolarLorentzVector> p4_vector;

}

using namespace B2gMiniToEdm;

//
// class declaration
//
class b2g_miniAodAnalyzer_general : public edm::EDFilter  
{
   public:
      	explicit b2g_miniAodAnalyzer_general(const edm::ParameterSet&);
      	~b2g_miniAodAnalyzer_general();
   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&); // this is essentially the Ntuplizer code
      virtual void endJob() ;
     // ----------member data ---------------------------
      	edm::EDGetTokenT<reco::VertexCollection> vertecies_;
      	edm::EDGetTokenT<pat::METCollection> met_;
      	edm::EDGetTokenT<pat::ElectronCollection> electrons_;
      	edm::EDGetTokenT<pat::MuonCollection> muons_;
      	edm::EDGetTokenT<pat::TauCollection> taus_;
      	edm::EDGetTokenT<pat::PhotonCollection> photons_;
      	edm::EDGetTokenT<pat::JetCollection> ak4slimmed_;
      	edm::EDGetTokenT<pat::JetCollection> ak8slimmed_;
	edm::EDGetTokenT<std::vector<pat::PackedCandidate> > pfcands_;
	bool ak8grommedMasses_;
};

b2g_miniAodAnalyzer_general::b2g_miniAodAnalyzer_general(const edm::ParameterSet& iConfig):
	// inputs from configuration:
    	vertecies_ 	(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertecies"))),
	met_		(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("met"))),
   	electrons_	(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
   	muons_		(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
   	taus_		(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
    	photons_	(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
    	ak4slimmed_	(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak4slimmed"))),
    	ak8slimmed_	(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak8slimmed"))),
	pfcands_	(consumes<std::vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("pfcands"))),
	ak8grommedMasses_	(iConfig.getParameter<bool>("ak8grommedMasses"))
{
	// Declare products
	produces<unsigned int>    ("numPV");
	produces<unsigned int>    ("numGoodPV");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("met");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("electrons");
	produces<std::vector<float> > ("electronRelIso");
	produces<std::vector<float> > ("electronRelIsoBetaCorr");
	produces<std::vector<float> > ("electronSuperClusterEta");
	produces<std::vector<float> > ("electronFull5by5");
	produces<std::vector<bool> > ("electronpassConversionVeto");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("muons");
	produces<std::vector<float> > ("muonRelIso");
	produces<std::vector<float> > ("muonRelIsoBetaCorr");
	produces<std::vector<bool> > ("muonIsLoose");
	produces<std::vector<bool> > ("muonIsSoft");
	produces<std::vector<bool> > ("muonIsTight");
	produces<std::vector<bool> > ("muonUserID");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("taus");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("photons");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("AK4slimmedjets");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("AK8slimmedjets");
	produces<std::vector<double> > ("AK4jetsCSV");
	produces<std::vector<double> > ("AK4jetsISV");
	if (ak8grommedMasses_) {
	    produces<std::vector<float> > ("AK8prunedMass");
	    produces<std::vector<float> > ("AK8trimmedMass");
	    produces<std::vector<float> > ("AK8filteredMass");
	    produces<std::vector<float> > ("AK8topMass");
        }
}
void b2g_miniAodAnalyzer_general::beginJob()
{

}

bool b2g_miniAodAnalyzer_general::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) //Make cuts/adjustments/corrections:
{
	// Initialize products
  	std::auto_ptr<unsigned int> npv( new unsigned int() );
  	std::auto_ptr<unsigned int> ngpv( new unsigned int() );
	std::auto_ptr<p4_vector> met( new p4_vector() );
	std::auto_ptr<p4_vector> electron( new p4_vector() );
	std::auto_ptr<std::vector<float>> electronRelIso( new std::vector<float> );
	std::auto_ptr<std::vector<float>> electronRelIsoBetaCorr( new std::vector<float> );
	std::auto_ptr<std::vector<float>> electronSuperClusterEta( new std::vector<float> );
	std::auto_ptr<std::vector<float>> electronFull5by5( new std::vector<float> );
	std::auto_ptr<std::vector<bool>> electronpassConversionVeto( new std::vector<bool> );
	std::auto_ptr<p4_vector> muon( new p4_vector() );
	std::auto_ptr<std::vector<float>> muonRelIso( new std::vector<float> );
	std::auto_ptr<std::vector<float>> muonRelIsoBetaCorr( new std::vector<float> );
	std::auto_ptr<std::vector<bool>> muonIsLoose( new std::vector<bool> );
	std::auto_ptr<std::vector<bool>> muonIsSoft( new std::vector<bool> );
	std::auto_ptr<std::vector<bool>> muonIsTight( new std::vector<bool> );
	std::auto_ptr<std::vector<bool>> muonUserID( new std::vector<bool> );
	std::auto_ptr<p4_vector> tau( new p4_vector() );
	std::auto_ptr<p4_vector> photon( new p4_vector() );
	std::auto_ptr<p4_vector> ak4slimmedjet( new p4_vector() );
	std::auto_ptr<p4_vector> ak8slimmedjet( new p4_vector() );
	std::auto_ptr<std::vector<double>> AK4jetsCSV( new std::vector<double> );
	std::auto_ptr<std::vector<double>> AK4jetsISV( new std::vector<double> );
	std::auto_ptr<std::vector<float>> AK8prunedMass( new std::vector<float> );
	std::auto_ptr<std::vector<float>> AK8trimmedMass( new std::vector<float> );
	std::auto_ptr<std::vector<float>> AK8filteredMass( new std::vector<float> );
	std::auto_ptr<std::vector<float>> AK8topMass( new std::vector<float> );
	// Read from Event:
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(vertecies_, vertices);

  	edm::Handle<pat::METCollection> mets;
  	iEvent.getByToken(met_, mets);

	edm::Handle<pat::PhotonCollection> photons;
	iEvent.getByToken(photons_, photons);

	edm::Handle<pat::ElectronCollection> electrons;
	iEvent.getByToken(electrons_, electrons);

  	edm::Handle<pat::MuonCollection> muons;
  	iEvent.getByToken(muons_, muons);

	edm::Handle<pat::TauCollection> taus;
	iEvent.getByToken(taus_, taus);

	edm::Handle<pat::JetCollection> ak4slimmedjets;
  	iEvent.getByToken(ak4slimmed_, ak4slimmedjets);

	edm::Handle<pat::JetCollection> ak8slimmedjets;
  	iEvent.getByToken(ak8slimmed_, ak8slimmedjets);

	edm::Handle<std::vector<pat::PackedCandidate> > pfcands;
  	iEvent.getByToken(pfcands_, pfcands);
	
	/////////////
	// Vertexes:
	if (vertices->empty()) return false; // skip the event if no PV found
  	unsigned int num_pv = 0;
	unsigned int num_gpv = 0;          
	for(reco::VertexCollection::const_iterator v=vertices->begin();v!=vertices->end(); ++v)
  	{
      		num_pv++;
      		if ( v->ndof() > 4 && fabs(v->z()) < 24 && fabs(v->position().rho()) < 2 ) num_gpv++;
  	}
	*npv = num_pv;
	*ngpv = num_gpv;
	const reco::Vertex &PV = vertices->front();

	/////////////
	// MET:
	const pat::MET &ievtmet = mets->front();
	reco::Candidate::PolarLorentzVector evtmet4 (ievtmet.pt(), (float) 0.0, ievtmet.phi(), (float) 0.0);
	met->push_back(evtmet4);

	/////////////
	// PHOTONS:
	for (const pat::Photon &ipho : *photons)
	{
		reco::Candidate::PolarLorentzVector pho4 (ipho.pt(), ipho.eta(), ipho.phi(), (float) 0.0);
		photon->push_back(pho4);

		if (ipho.pt() < 20 or ipho.chargedHadronIso()/ipho.pt() > 0.3) continue;
                ipho.pt();
		ipho.superCluster()->eta();
		ipho.sigmaIetaIeta();
	}

	/////////////
	// ELECTRONS:
	for (const pat::Electron &iel : *electrons)
	{
		reco::Candidate::PolarLorentzVector elec4 (iel.pt(), iel.eta(), iel.phi(), (float) 0.00051);
		electron->push_back(elec4);
		if (iel.pt() < 5) 
		{
			electronSuperClusterEta->push_back(0);
			electronFull5by5->push_back(0);
			electronpassConversionVeto->push_back(false);
		} else {
		// save electron variables
		electronSuperClusterEta->push_back(iel.superCluster()->eta());
		electronFull5by5->push_back(iel.full5x5_sigmaIetaIeta());
		electronpassConversionVeto->push_back(iel.passConversionVeto());
		// Compute electron isolation
	        double charged = 0, neutral = 0, pileup  = 0;
	        // get a list of the PF candidates used to build this lepton, so to exclude them
	        std::vector<reco::CandidatePtr> footprint;
	        for (unsigned int i = 0, n = iel.numberOfSourceCandidatePtrs(); i < n; ++i) {
	        	footprint.push_back(iel.sourceCandidatePtr(i));
	        }
	        // now loop on pf candidates
	        for (unsigned int i = 0, n = pfcands->size(); i < n; ++i) {
	        	const pat::PackedCandidate &pf = (*pfcands)[i];
	        	if (deltaR(pf,iel) < 0.3) {
	                // pfcandidate-based footprint removal
	                if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfcands,i)) != footprint.end()) {
	                	continue;
	                }
	                if (pf.charge() == 0) {
	                	if (pf.pt() > 0.5) neutral += pf.pt();
	                } else if (pf.fromPV() >= 2) {
	                	charged += pf.pt();
	                } else {
	                	if (pf.pt() > 0.5) pileup += pf.pt();
	                }
	            }
	        }
		// calculate desired iso
		float reliso = (charged+neutral) / iel.pt();
	        float reliso_corr = (charged + std::max(0.0, neutral-0.5*pileup)) / (iel.pt());
		electronRelIso->push_back(reliso);
		electronRelIsoBetaCorr->push_back(reliso_corr);	
		}
	}

	// MUONS:
	for (const pat::Muon &imu : *muons)
	{
		reco::Candidate::PolarLorentzVector muon4 (imu.pt(), imu.eta(), imu.phi(), (float) 0.1057);
		muon->push_back(muon4);
        	if (imu.pt() < 5 && !imu.isLooseMuon())
		{
			muonIsLoose->push_back( false );
			muonIsSoft->push_back( false );
			muonIsTight->push_back( false );
			muonUserID->push_back( false );
			muonRelIso->push_back(0);
			muonRelIsoBetaCorr->push_back(0);
		} else {	
		// POG muon IDs
		muonIsLoose->push_back( imu.isLooseMuon() );
		muonIsSoft->push_back( imu.isSoftMuon(PV) );
		muonIsTight->push_back( imu.isTightMuon(PV) );
		// User muon ID
		if ( 
		imu.isPFMuon() &&
		imu.isGlobalMuon() &&
		imu.isTrackerMuon() &&
		imu.normChi2() < 10 &&
		imu.track()->hitPattern().trackerLayersWithMeasurement() > 5 &&
		imu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
		fabs(imu.dB()) < 0.2 &&
		fabs(imu.vertex().z() - PV.z()) < 0.5 &&
		imu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
		imu.numberOfMatchedStations() > 1 )
			muonUserID->push_back( true );
		else
			muonUserID->push_back( false );
		// Compute muon isolation
	        double charged = 0, neutral = 0, pileup  = 0;
	        // get a list of the PF candidates used to build this lepton, so to exclude them
	        std::vector<reco::CandidatePtr> footprint;
	        for (unsigned int i = 0, n = imu.numberOfSourceCandidatePtrs(); i < n; ++i) {
	        	footprint.push_back(imu.sourceCandidatePtr(i));
	        }
	        // now loop on pf candidates
	        for (unsigned int i = 0, n = pfcands->size(); i < n; ++i) {
	        	const pat::PackedCandidate &pf = (*pfcands)[i];
	        	if (deltaR(pf,imu) < 0.3) {
	                // pfcandidate-based footprint removal
	                if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfcands,i)) != footprint.end()) {
	                	continue;
	                }
	                if (pf.charge() == 0) {
	                	if (pf.pt() > 0.5) neutral += pf.pt();
	                } else if (pf.fromPV() >= 2) {
	                	charged += pf.pt();
	                } else {
	                	if (pf.pt() > 0.5) pileup += pf.pt();
	                }
	            }
	        }
		// calculate desired iso
		float reliso = (charged+neutral) / imu.pt();
	        float reliso_corr = (charged + std::max(0.0, neutral-0.5*pileup)) / (imu.pt());
		muonRelIso->push_back(reliso);
		muonRelIsoBetaCorr->push_back(reliso_corr);
		}
	}

	/////////////
	// TAUS:
	for (const pat::Tau &itau : *taus)
	{
		reco::Candidate::PolarLorentzVector tau4 (itau.pt(), itau.eta(), itau.phi(), (float) 1.777);
		tau->push_back(tau4);
	}

	/////////////
	// JETS:
	/// slimmed jets (ak4)
	for (const pat::Jet &iak4jet : *ak4slimmedjets)
	{
		reco::Candidate::PolarLorentzVector ak4jet4 (iak4jet.pt(), iak4jet.eta(), iak4jet.phi(), iak4jet.mass());
		ak4slimmedjet->push_back(ak4jet4);
		AK4jetsCSV->push_back(std::max(0.f,iak4jet.bDiscriminator("combinedSecondaryVertexBJetTags")));
		AK4jetsISV->push_back(std::max(0.f,iak4jet.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags")));
	}

	/// slimmed fat jets (ak8)
	for (const pat::Jet &iak8jet : *ak8slimmedjets)
	{
		reco::Candidate::PolarLorentzVector ak8jet4 (iak8jet.pt(), iak8jet.eta(), iak8jet.phi(), iak8jet.mass());
		ak8slimmedjet->push_back(ak8jet4);
		if (ak8grommedMasses_) {
			AK8prunedMass->push_back(iak8jet.userFloat("ak8PFJetsCHSPrunedLinks"));
			AK8trimmedMass->push_back(iak8jet.userFloat("ak8PFJetsCHSTrimmedLinks"));
			AK8filteredMass->push_back(iak8jet.userFloat("ak8PFJetsCHSFilteredLinks"));
			AK8topMass->push_back(iak8jet.userFloat("cmsTopTagPFJetsCHSLinksAK8"));
                }
	}

	// Fill Everything:
  	iEvent.put( npv, "numPV");
  	iEvent.put( ngpv, "numGoodPV");
  	iEvent.put( met, "met");
  	iEvent.put( photon, "photons");
  	iEvent.put( electron, "electrons");
  	iEvent.put( electronRelIso, "electronRelIso");
  	iEvent.put( electronRelIsoBetaCorr, "electronRelIsoBetaCorr");
  	iEvent.put( electronSuperClusterEta, "electronSuperClusterEta");
  	iEvent.put( electronFull5by5, "electronFull5by5");
  	iEvent.put( electronpassConversionVeto, "electronpassConversionVeto");
  	iEvent.put( muon, "muons");
  	iEvent.put( muonRelIsoBetaCorr, "muonRelIsoBetaCorr");
  	iEvent.put( muonRelIso, "muonRelIso");
  	iEvent.put( muonIsLoose, "muonIsLoose");
  	iEvent.put( muonIsSoft, "muonIsSoft");
  	iEvent.put( muonIsTight, "muonIsTight");
  	iEvent.put( muonUserID, "muonUserID");
  	iEvent.put( tau, "taus");
  	iEvent.put( ak4slimmedjet, "AK4slimmedjets");
  	iEvent.put( ak8slimmedjet, "AK8slimmedjets");
  	iEvent.put( AK4jetsCSV, "AK4jetsCSV");
  	iEvent.put( AK4jetsISV, "AK4jetsISV");
        if (ak8grommedMasses_) {
  		iEvent.put( AK8prunedMass, "AK8prunedMass");
	  	iEvent.put( AK8trimmedMass, "AK8trimmedMass");
	  	iEvent.put( AK8filteredMass, "AK8filteredMass");
	  	iEvent.put( AK8topMass, "AK8topMass");
	}
	return true;
}


void b2g_miniAodAnalyzer_general::endJob() 
{

}


b2g_miniAodAnalyzer_general::~b2g_miniAodAnalyzer_general(){}


//define this as a plug-in
DEFINE_FWK_MODULE(b2g_miniAodAnalyzer_general);
