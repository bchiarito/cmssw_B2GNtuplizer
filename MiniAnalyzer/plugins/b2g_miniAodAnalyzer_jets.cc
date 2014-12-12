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
class b2g_miniAodAnalyzer_jets : public edm::EDFilter  
{
   public:
      	explicit b2g_miniAodAnalyzer_jets(const edm::ParameterSet&);
      	~b2g_miniAodAnalyzer_jets();
   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
     // ----------member data ---------------------------
      	edm::EDGetTokenT<reco::VertexCollection> vertecies_;
      	edm::EDGetTokenT<std::vector<pat::PackedCandidate> > pfcands_;
	std::string 	jetAlgo_;
	double 		jetR_;
	double 		jetPtmin_;
};

b2g_miniAodAnalyzer_jets::b2g_miniAodAnalyzer_jets(const edm::ParameterSet& iConfig):
	// inputs:
    	vertecies_ 	(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertecies"))),
    	pfcands_	(consumes<std::vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("pfcands"))),
	jetAlgo_	(iConfig.getParameter<std::string>("jetAlgo")),
	jetR_		(iConfig.getParameter<double>("jetR")),
	jetPtmin_	(iConfig.getParameter<double>("jetPtmin"))
{
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("fjjets");
}
void b2g_miniAodAnalyzer_jets::beginJob()
{

}

bool b2g_miniAodAnalyzer_jets::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) //Make cuts/adjustments/corrections:
{
	std::auto_ptr<p4_vector> fjjet( new p4_vector() );
	// Read from miniAOD:
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(vertecies_, vertices);

	edm::Handle<std::vector<pat::PackedCandidate> > pfcands;
  	iEvent.getByToken(pfcands_, pfcands);
	
	/// jets built from scratch with fastjet
	//make collection of pf candidates
	std::vector<fastjet::PseudoJet> FJparticles;
	for (const pat::PackedCandidate &pfcand : *pfcands) {
		if (pfcand.fromPV()>1)
        		FJparticles.push_back( fastjet::PseudoJet( pfcand.px(), pfcand.py(), pfcand.pz(), pfcand.energy() ) );	
	}
	//let fastjet cluster pf candidates
        std::vector<fastjet::PseudoJet> fjjets;	
	if (jetAlgo_.compare("CA") == 0) {
		fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, jetR_);
        	fastjet::ClusterSequence clust_seq(FJparticles, jet_def);
        	fjjets = sorted_by_pt(clust_seq.inclusive_jets(jetPtmin_));	
	}
	if (jetAlgo_.compare("AK") == 0) {
		fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jetR_);
        	fastjet::ClusterSequence clust_seq(FJparticles, jet_def);
        	fjjets = sorted_by_pt(clust_seq.inclusive_jets(jetPtmin_));	
	}
	if (jetAlgo_.compare("KT") == 0) {
		fastjet::JetDefinition jet_def(fastjet::kt_algorithm, jetR_);
        	fastjet::ClusterSequence clust_seq(FJparticles, jet_def);
        	fjjets = sorted_by_pt(clust_seq.inclusive_jets(jetPtmin_));	
	}
	//get collection of clustered jets from fastjet
	for (const fastjet::PseudoJet &ifjjet : fjjets)
	{
		reco::Candidate::PolarLorentzVector fjjet4 (ifjjet.pt(), ifjjet.eta(), ifjjet.phi(), ifjjet.m());
		fjjet->push_back(fjjet4);
	}

	// Fill:
  	iEvent.put( fjjet, "fjjets" );
	return true;
}


void b2g_miniAodAnalyzer_jets::endJob() 
{

}


b2g_miniAodAnalyzer_jets::~b2g_miniAodAnalyzer_jets(){}


//define this as a plug-in
DEFINE_FWK_MODULE(b2g_miniAodAnalyzer_jets);
