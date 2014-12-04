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
      virtual bool filter(edm::Event&, const edm::EventSetup&); // this is essentially the Ntuplizer code
      virtual void endJob() ;
     // ----------member data ---------------------------
      	edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
      	edm::EDGetTokenT<pat::METCollection> metToken_;
      	edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      	edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      	edm::EDGetTokenT<pat::TauCollection> tauToken_;
      	edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
      	edm::EDGetTokenT<pat::JetCollection> ak4slimmedToken_;
      	edm::EDGetTokenT<pat::JetCollection> ak8slimmedToken_;
      	edm::EDGetTokenT<std::vector<pat::PackedCandidate> > pfcandsToken_;
      	edm::EDGetTokenT<edm::PtrVector<reco::Candidate> > chscandsToken_;
      	edm::EDGetTokenT<std::vector<reco::PFJet> > ca8Token_;
      	edm::EDGetTokenT<edm::ValueMap<float> > ca8tau1Token_;
      	edm::EDGetTokenT<edm::ValueMap<float> > ca8tau2Token_;
      	edm::EDGetTokenT<edm::ValueMap<float> > ca8tau3Token_;
      	edm::EDGetTokenT<std::vector<reco::BasicJet> > ca8prunedToken_;
      	edm::EDGetTokenT<std::vector<reco::PFJet> > ca8subjetsToken_;
      	edm::EDGetTokenT<std::vector<reco::BasicJet> > ca8TopTagToken_;
      	edm::EDGetTokenT<std::vector<reco::CATopJetTagInfo> > ca8TopTagInfoToken_;
      	edm::EDGetTokenT<std::vector<reco::PFJet> > ak4Token_;
      	edm::EDGetTokenT<std::vector<reco::PFJet> > ak8Token_;
	std::string 	CMSTopTagString_;
};

b2g_miniAodAnalyzer_jets::b2g_miniAodAnalyzer_jets(const edm::ParameterSet& iConfig):
	// inputs:
    	vertexToken_ 	(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexToken"))),
	metToken_	(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metToken"))),
   	electronToken_	(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronToken"))),
   	muonToken_	(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonToken"))),
   	tauToken_	(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tauToken"))),
    	photonToken_	(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photonToken"))),
    	ak4slimmedToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak4slimmedToken"))),
    	ak8slimmedToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak8slimmedToken"))),
    	pfcandsToken_	(consumes<std::vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("pfcandsToken"))),
    	chscandsToken_	(consumes<edm::PtrVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("chscandsToken"))),
    	ca8Token_	(consumes<std::vector<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("ca8Token"))),
    	ca8tau1Token_	(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("ca8tau1Token"))),
    	ca8tau2Token_	(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("ca8tau2Token"))),
    	ca8tau3Token_	(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("ca8tau3Token"))),
    	ca8prunedToken_	(consumes<std::vector<reco::BasicJet> >(iConfig.getParameter<edm::InputTag>("ca8prunedToken"))),
    	ca8subjetsToken_(consumes<std::vector<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("ca8subjetsToken"))),
    	ca8TopTagToken_	(consumes<std::vector<reco::BasicJet> >(iConfig.getParameter<edm::InputTag>("ca8TopTagToken"))),
    	ca8TopTagInfoToken_	(consumes<std::vector<reco::CATopJetTagInfo> >(iConfig.getParameter<edm::InputTag>("ca8TopTagInfoToken")))
{
	produces<unsigned int>    ("numPV");
	produces<unsigned int>    ("numGoodPV");
	produces<unsigned int>    ("numElecs");
	produces<unsigned int>    ("numMuons");
	produces<unsigned int>    ("numTaus");
	produces<unsigned int>    ("numPhotons");
	produces<unsigned int>    ("numAK4jets");
	produces<unsigned int>    ("numAK8jets");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("met");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("electrons");
	produces<std::vector<float> > ("electronsRelIso");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("muons");
	produces<std::vector<float> > ("muonsRelIso");
	produces<std::vector<bool> > ("muonsIsTight");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("taus");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("photons");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("AK4slimmedjets");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("AK8slimmedjets");
	produces<std::vector<float> > ("AK8prunedMass");
	produces<std::vector<float> > ("AK8trimmedMass");
	produces<std::vector<float> > ("AK8filteredMass");
	produces<std::vector<float> > ("AK8topMass");
	produces<std::vector<double> > ("AK4jetsCSV");
	produces<std::vector<double> > ("AK4jetsISV");
	//
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("CA8fjjets");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("AK4fjjets");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("AK8fjjets");
	//
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("CA8jets");
	produces<std::vector<float> > ("CA8tau1");
	produces<std::vector<float> > ("CA8tau2");
	produces<std::vector<float> > ("CA8tau3");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("CA8prunedjets");
	produces<std::vector<reco::Candidate::PolarLorentzVector> > ("TopTagjets");
	produces<std::vector<unsigned int> > ("TopTagnsub");
	produces<std::vector<float> > ("TopTagminmass");
	produces<std::vector<float> > ("TopTagWmass");
	produces<std::vector<float> > ("TopTagtopmass");
}
void b2g_miniAodAnalyzer_jets::beginJob()
{

}

bool b2g_miniAodAnalyzer_jets::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) //Make cuts/adjustments/corrections:
{
  	std::auto_ptr<unsigned int> npv( new unsigned int() );
  	std::auto_ptr<unsigned int> ngpv( new unsigned int() );
  	std::auto_ptr<unsigned int> nel( new unsigned int() );
  	std::auto_ptr<unsigned int> nmu( new unsigned int() );
  	std::auto_ptr<unsigned int> ntau( new unsigned int() );
  	std::auto_ptr<unsigned int> np( new unsigned int() );
  	std::auto_ptr<unsigned int> nak4( new unsigned int() );
  	std::auto_ptr<unsigned int> nak8( new unsigned int() );
	std::auto_ptr<p4_vector> met( new p4_vector() );
	std::auto_ptr<p4_vector> electron( new p4_vector() );
	std::auto_ptr<std::vector<float>> electronRelIso( new std::vector<float> );
	std::auto_ptr<p4_vector> muon( new p4_vector() );
	std::auto_ptr<std::vector<float>> muonRelIso( new std::vector<float> );
	std::auto_ptr<std::vector<bool>> muonIsTight( new std::vector<bool> );
	std::auto_ptr<p4_vector> tau( new p4_vector() );
	std::auto_ptr<p4_vector> photon( new p4_vector() );
	std::auto_ptr<p4_vector> ak4slimmedjet( new p4_vector() );
	std::auto_ptr<p4_vector> ak8slimmedjet( new p4_vector() );
	std::auto_ptr<std::vector<float>> AK8prunedMass( new std::vector<float> );
	std::auto_ptr<std::vector<float>> AK8trimmedMass( new std::vector<float> );
	std::auto_ptr<std::vector<float>> AK8filteredMass( new std::vector<float> );
	std::auto_ptr<std::vector<float>> AK8topMass( new std::vector<float> );
	std::auto_ptr<std::vector<double>> AK4jetsCSV( new std::vector<double> );
	std::auto_ptr<std::vector<double>> AK4jetsISV( new std::vector<double> );
	//
	std::auto_ptr<p4_vector> ca8fjjet( new p4_vector() );
	std::auto_ptr<p4_vector> ak4fjjet( new p4_vector() );
	std::auto_ptr<p4_vector> ak8fjjet( new p4_vector() );
	//
	std::auto_ptr<p4_vector> ca8jet( new p4_vector() );
	std::auto_ptr<std::vector<float>> ca8tau1( new std::vector<float> );
	std::auto_ptr<std::vector<float>> ca8tau2( new std::vector<float> );
	std::auto_ptr<std::vector<float>> ca8tau3( new std::vector<float> );
	std::auto_ptr<p4_vector> ca8prunedjet( new p4_vector() );
	std::auto_ptr<p4_vector> toptagjet( new p4_vector() );
	std::auto_ptr<std::vector<unsigned int>> toptagnsub( new std::vector<unsigned int> );
	std::auto_ptr<std::vector<float>> toptagminmass( new std::vector<float> );
	std::auto_ptr<std::vector<float>> toptagWmass( new std::vector<float> );
	std::auto_ptr<std::vector<float>> toptagtopmass( new std::vector<float> );
	// Read from miniAOD:
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(vertexToken_, vertices);

  	edm::Handle<pat::METCollection> mets;
  	iEvent.getByToken(metToken_, mets);

	edm::Handle<pat::PhotonCollection> photons;
	iEvent.getByToken(photonToken_, photons);

	edm::Handle<pat::ElectronCollection> electrons;
	iEvent.getByToken(electronToken_, electrons);

  	edm::Handle<pat::MuonCollection> muons;
  	iEvent.getByToken(muonToken_, muons);

	edm::Handle<pat::TauCollection> taus;
	iEvent.getByToken(tauToken_, taus);

	edm::Handle<pat::JetCollection> ak4slimmedjets;
  	iEvent.getByToken(ak4slimmedToken_, ak4slimmedjets);

	edm::Handle<pat::JetCollection> ak8slimmedjets;
  	iEvent.getByToken(ak8slimmedToken_, ak8slimmedjets);

	edm::Handle<std::vector<pat::PackedCandidate> > pfcands;
  	iEvent.getByToken(pfcandsToken_, pfcands);
	
	edm::Handle<edm::PtrVector<reco::Candidate> > chscands;
  	iEvent.getByToken(chscandsToken_, chscands);
	
	edm::Handle<std::vector<reco::PFJet> > ca8jets;
  	iEvent.getByToken(ca8Token_, ca8jets);
	
	edm::Handle<edm::ValueMap<float> > ca8tau1s;
  	iEvent.getByToken(ca8tau1Token_, ca8tau1s);
	
	edm::Handle<edm::ValueMap<float> > ca8tau2s;
  	iEvent.getByToken(ca8tau2Token_, ca8tau2s);
	
	edm::Handle<edm::ValueMap<float> > ca8tau3s;
  	iEvent.getByToken(ca8tau3Token_, ca8tau3s);
	
	edm::Handle<std::vector<reco::BasicJet> > ca8prunedjets;
  	iEvent.getByToken(ca8prunedToken_, ca8prunedjets);

	edm::Handle<std::vector<reco::BasicJet> > toptagjets;
  	iEvent.getByToken(ca8TopTagToken_, toptagjets);
	
	edm::Handle<std::vector<reco::CATopJetTagInfo> > toptaginfos;
  	iEvent.getByToken(ca8TopTagInfoToken_, toptaginfos);
	
	// Here is the work:
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
	/////////////
	// MET:
	const pat::MET &ievtmet = mets->front();
	reco::Candidate::PolarLorentzVector evtmet4 (ievtmet.pt(), (float) 0.0, ievtmet.phi(), (float) 0.0);
	met->push_back(evtmet4);

	/////////////
	// PHOTONS:

	*np = photons->size();
	for (const pat::Photon &ipho : *photons)
	{
		reco::Candidate::PolarLorentzVector pho4 (ipho.pt(), ipho.eta(), ipho.phi(), (float) 0.0);
		photon->push_back(pho4);
	}

	/////////////
	// LEPTONS (except taus):
	// Electons:
	*nel = electrons->size();
	for (const pat::Electron &iel : *electrons)
	{
		reco::Candidate::PolarLorentzVector elec4 (iel.pt(), iel.eta(), iel.phi(), (float) 0.00051);
		electron->push_back(elec4);
		float caloiso = iel.caloIso(); // sum of ecalIs and hcalIso, summed Et of all recHits in a cone of deltaR<0.4
		float trackiso = iel.trackIso(); // summed track pt in a cone of deltaR<0.4
		float reliso = (caloiso + trackiso) / (iel.pt());
		electronRelIso->push_back(reliso);
	}

	// Muons:
	*nmu = muons->size();
	for (const pat::Muon &imu : *muons)
	{
		reco::Candidate::PolarLorentzVector muon4 (imu.pt(), imu.eta(), imu.phi(), (float) 0.1057);
		muon->push_back(muon4);
		float caloiso = imu.caloIso(); // sum of ecalIs and hcalIso, summed Et of all recHits in a cone of deltaR<0.4
		float trackiso = imu.trackIso(); // summed track pt in a cone of deltaR<0.4
		float reliso = (caloiso + trackiso) / (imu.pt());
		muonRelIso->push_back(reliso);
		muonIsTight->push_back(imu.isTightMuon(vertices->at(0)));
	}
	/////////////
	// TAUS:
	*ntau = taus->size();
	for (const pat::Tau &itau : *taus)
	{
		reco::Candidate::PolarLorentzVector tau4 (itau.pt(), itau.eta(), itau.phi(), (float) 1.777);
		tau->push_back(tau4);
	}

	/////////////
	// JETS:
	/// slimmed jets (ak4)
	*nak4 = ak4slimmedjets->size();
	for (const pat::Jet &iak4jet : *ak4slimmedjets)
	{
		reco::Candidate::PolarLorentzVector ak4jet4 (iak4jet.pt(), iak4jet.eta(), iak4jet.phi(), iak4jet.mass());
		ak4slimmedjet->push_back(ak4jet4);
		AK4jetsCSV->push_back(std::max(0.f,iak4jet.bDiscriminator("combinedSecondaryVertexBJetTags")));
		AK4jetsISV->push_back(std::max(0.f,iak4jet.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags")));
	}

	/// slimmed fat jets (ak8)
	*nak8 = ak8slimmedjets->size();
	for (const pat::Jet &iak8jet : *ak8slimmedjets)
	{
		reco::Candidate::PolarLorentzVector ak8jet4 (iak8jet.pt(), iak8jet.eta(), iak8jet.phi(), iak8jet.mass());
		ak8slimmedjet->push_back(ak8jet4);
		AK8prunedMass->push_back(iak8jet.userFloat("ak8PFJetsCHSPrunedLinks"));
		AK8trimmedMass->push_back(iak8jet.userFloat("ak8PFJetsCHSTrimmedLinks"));
		AK8filteredMass->push_back(iak8jet.userFloat("ak8PFJetsCHSFilteredLinks"));
		AK8topMass->push_back(iak8jet.userFloat("cmsTopTagPFJetsCHSLinksAK8"));
	}

	/// ca8 jets
	for (const reco::PFJet &ica8jet : *ca8jets)
	{
		reco::Candidate::PolarLorentzVector ca8jet4 (ica8jet.pt(), ica8jet.eta(), ica8jet.phi(), ica8jet.mass());
		ca8jet->push_back(ca8jet4);
	}

	/// ca8 taus
	edm::ValueMap<float>::const_iterator iter;
	iter = ca8tau1s->begin();
	for (unsigned int i = 0; i < iter.size(); i++ )
	{
		ca8tau1->push_back( iter[i] );
	}
	iter = ca8tau2s->begin();
	for (unsigned int i = 0; i < iter.size(); i++ )
	{
		ca8tau2->push_back( iter[i] );
	}
	iter = ca8tau3s->begin();
	for (unsigned int i = 0; i < iter.size(); i++ )
	{
		ca8tau3->push_back( iter[i] );
	}

	/// ca8 pruned jets
	for (const reco::BasicJet &ijet : *ca8prunedjets)
	{
		reco::Candidate::PolarLorentzVector jet4 (ijet.pt(), ijet.eta(), ijet.phi(), ijet.mass());
		//if (jet4.Pt() > 20 && jet4.Eta() < 2.4 && jet4.Eta() > -2.4) 
		//{
			ca8prunedjet->push_back(jet4);
		//} // kinematic cut on jets
	}
	
	/// toptag jets and info
	std::vector<reco::BasicJet>::const_iterator iter_jet;
	std::vector<reco::CATopJetTagInfo>::const_iterator iter_info;
	iter_jet = toptagjets->begin();
	iter_info = toptaginfos->begin();
	for (unsigned int i = 0; i < toptagjets->size(); i++ )
	{
		reco::Candidate::PolarLorentzVector jet4 (iter_jet[i].pt(), iter_jet[i].eta(), iter_jet[i].phi(), iter_jet[i].mass());
		const reco::CATopJetTagInfo info = iter_info[i];
		//if (jet4.Pt() > 20 && jet4.Eta() < 2.4 && jet4.Eta() > -2.4)
		//{
			toptagjet->push_back(jet4);
			toptagnsub->push_back(info.properties().nSubJets);
			toptagminmass->push_back(info.properties().minMass);
			toptagWmass->push_back(info.properties().wMass);
			toptagtopmass->push_back(info.properties().topMass);
		//} // kinematic cut on jets
	}

	/// ca8 jets built from scratch with fastjet
	//make collection of pf candidates
	std::vector<fastjet::PseudoJet> FJparticles;
	for (const reco::Candidate &pfcand : *pfcands) {
        	FJparticles.push_back( fastjet::PseudoJet( pfcand.px(), pfcand.py(), pfcand.pz(), pfcand.energy() ) );	
	}
	std::vector<fastjet::PseudoJet> FJCHSparticles;
	for (const edm::Ptr<reco::Candidate> &pfcand : *chscands) {
        	FJCHSparticles.push_back( fastjet::PseudoJet( pfcand->px(), pfcand->py(), pfcand->pz(), pfcand->energy() ) );	
	}
	//let fastjet cluster pf candidates
	fastjet::JetDefinition jet_def_ca8(fastjet::cambridge_algorithm, 0.8);
        fastjet::ClusterSequence clust_seq_ca8(FJCHSparticles, jet_def_ca8);
        double ptmin = 20.0;
	//get collection of clustered jets from fastjet
        std::vector<fastjet::PseudoJet> ca8fjjets = sorted_by_pt(clust_seq_ca8.inclusive_jets(ptmin));	
	//fill collections
	for (const fastjet::PseudoJet &ica8fjjet : ca8fjjets)
	{
		reco::Candidate::PolarLorentzVector ca8fjjet4 (ica8fjjet.pt(), ica8fjjet.eta(), ica8fjjet.phi(), ica8fjjet.m());
		ca8fjjet->push_back(ca8fjjet4);
	}

	/// ak4 jets built from scratch
	//let fastjet cluster pf candidates
	fastjet::JetDefinition jet_def_ak4(fastjet::antikt_algorithm, 0.4);
        fastjet::ClusterSequence clust_seq_ak4(FJCHSparticles, jet_def_ak4);
        ptmin = 10.0;
	//get collection of clustered jets from fastjet
        std::vector<fastjet::PseudoJet> ak4fjjets = sorted_by_pt(clust_seq_ak4.inclusive_jets(ptmin));	
	//fill collections
	for (const fastjet::PseudoJet &iak4fjjet : ak4fjjets)
	{
		reco::Candidate::PolarLorentzVector ak4fjjet4 (iak4fjjet.pt(), iak4fjjet.eta(), iak4fjjet.phi(), iak4fjjet.m());
		ak4fjjet->push_back(ak4fjjet4);
	}

	/// ak8 jets built from scratch
	//let fastjet cluster pf candidates
	fastjet::JetDefinition jet_def_ak8(fastjet::antikt_algorithm, 0.8);
        fastjet::ClusterSequence clust_seq_ak8(FJCHSparticles, jet_def_ak8);
        ptmin = 100.0;
	//get collection of clustered jets from fastjet
        std::vector<fastjet::PseudoJet> ak8fjjets = sorted_by_pt(clust_seq_ak8.inclusive_jets(ptmin));	
	//fill collections
	for (const fastjet::PseudoJet &iak8fjjet : ak8fjjets)
	{
		reco::Candidate::PolarLorentzVector ak8fjjet4 (iak8fjjet.pt(), iak8fjjet.eta(), iak8fjjet.phi(), iak8fjjet.m());
		ak8fjjet->push_back(ak8fjjet4);
	}

	// Fill Everything:
  	iEvent.put( npv, "numPV");
  	iEvent.put( ngpv, "numGoodPV");
  	iEvent.put( nel, "numElecs");
  	iEvent.put( nmu, "numMuons");
  	iEvent.put( ntau, "numTaus");
  	iEvent.put( nak4, "numAK4jets");
  	iEvent.put( nak8, "numAK8jets");
  	iEvent.put( met, "met");
  	iEvent.put( photon, "photons");
  	iEvent.put( electron, "electrons");
  	iEvent.put( electronRelIso, "electronsRelIso");
  	iEvent.put( muon, "muons");
  	iEvent.put( muonRelIso, "muonsRelIso");
  	iEvent.put( muonIsTight, "muonsIsTight");
  	iEvent.put( tau, "taus");
  	iEvent.put( ak4slimmedjet, "AK4slimmedjets");
  	iEvent.put( ak8slimmedjet, "AK8slimmedjets");
  	iEvent.put( AK4jetsCSV, "AK4jetsCSV");
  	iEvent.put( AK4jetsISV, "AK4jetsISV");
  	iEvent.put( AK8prunedMass, "AK8prunedMass");
  	iEvent.put( AK8trimmedMass, "AK8trimmedMass");
  	iEvent.put( AK8filteredMass, "AK8filteredMass");
  	iEvent.put( AK8topMass, "AK8topMass");
	//
  	iEvent.put( ca8fjjet, "CA8fjjets");
  	iEvent.put( ak4fjjet, "AK4fjjets");
  	iEvent.put( ak8fjjet, "AK8fjjets");
	//
  	iEvent.put( ca8jet, "CA8jets");
  	iEvent.put( ca8tau1, "CA8tau1");
  	iEvent.put( ca8tau2, "CA8tau2");
  	iEvent.put( ca8tau3, "CA8tau3");
  	iEvent.put( ca8prunedjet, "CA8prunedjets");
  	iEvent.put( toptagjet, "TopTagjets");
  	iEvent.put( toptagnsub, "TopTagnsub");
  	iEvent.put( toptagminmass, "TopTagminmass");
  	iEvent.put( toptagWmass, "TopTagWmass");
  	iEvent.put( toptagtopmass, "TopTagtopmass");
	return true;
}


void b2g_miniAodAnalyzer_jets::endJob() 
{

}


b2g_miniAodAnalyzer_jets::~b2g_miniAodAnalyzer_jets(){}


//define this as a plug-in
DEFINE_FWK_MODULE(b2g_miniAodAnalyzer_jets);
