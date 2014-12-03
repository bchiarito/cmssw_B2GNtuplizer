#include <memory>
#include <cmath>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

class b2g_miniAodAnalyzer_trigger : public edm::EDFilter {
   public:
      explicit b2g_miniAodAnalyzer_trigger(const edm::ParameterSet&);
      ~b2g_miniAodAnalyzer_trigger();

   private:
      virtual void beginJob();
      virtual void endJob();
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      // ----------member data ---------------------------
      bool printAll_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      bool isData_;
      bool useTriggerList_;
      std::vector<std::string> triggerList_;
};

b2g_miniAodAnalyzer_trigger::b2g_miniAodAnalyzer_trigger(const edm::ParameterSet& iConfig):
    printAll_(iConfig.getParameter<bool>("printAll")),
    triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
    triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
    isData_(iConfig.getParameter<bool>("isData")),
    useTriggerList_(iConfig.getParameter<bool>("useTriggerList")),
    triggerList_(iConfig.getParameter<std::vector<std::string>>("triggerList"))
{
    produces<std::vector<std::string> >  ("triggerNames");
    produces<std::vector<bool> >    	 ("triggerBits");
    produces<std::vector<unsigned int> > ("triggerPrescales");
    produces<std::vector<pat::TriggerObjectStandAlone> > ("triggerObjects");
}

void b2g_miniAodAnalyzer_trigger::beginJob()
{
}

void b2g_miniAodAnalyzer_trigger::endJob()
{
}

bool b2g_miniAodAnalyzer_trigger::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<std::vector<std::string>> triggerName( new std::vector<std::string> );
    std::auto_ptr<std::vector<bool>> triggerBit( new std::vector<bool> );
    std::auto_ptr<std::vector<unsigned int>> triggerPrescale( new std::vector<unsigned int> );
    std::auto_ptr<std::vector<pat::TriggerObjectStandAlone>> triggerObject( new std::vector<pat::TriggerObjectStandAlone> );

    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

    iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(triggerObjects_, triggerObjects);
    iEvent.getByToken(triggerPrescales_, triggerPrescales);

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    if (printAll_) std::cout << "\n === TRIGGER PATHS === " << std::endl;
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
        if (printAll_) {        
            std::cout << "Trigger " << names.triggerName(i) << 
            ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
            ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") 
            << std::endl;
        }
        if (useTriggerList_) {
	    if (std::find(triggerList_.begin(), triggerList_.end(), names.triggerName(i)) != triggerList_.end()) {
                triggerName->push_back(names.triggerName(i));
                triggerBit->push_back(triggerBits->accept(i));
                triggerPrescale->push_back(triggerPrescales->getPrescaleForIndex(i)); }
        } else {
            triggerName->push_back(names.triggerName(i));
            triggerBit->push_back(triggerBits->accept(i));
            triggerPrescale->push_back(triggerPrescales->getPrescaleForIndex(i)); }

    }
    if (isData_) {
    if (printAll_) std::cout << "\n === TRIGGER OBJECTS === " << std::endl;
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
        obj.unpackPathNames(names);
	std::vector<std::string> pathNamesAll  = obj.pathNames(false);
	std::vector<std::string> pathNamesLast = obj.pathNames(true);
        if (printAll_) {
            std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
            // Print trigger object collection and type
            std::cout << "\t   Collection: " << obj.collection() << std::endl;
	    std::cout << "\t   Type IDs:   ";
	    for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
	    std::cout << std::endl;
	    // Print associated trigger filters
	    std::cout << "\t   Filters:    ";
	    for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
	    std::cout << std::endl;
	    // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
	    // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
	    // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
	    std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
	    for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
		bool isBoth = obj.hasPathName( pathNamesAll[h], true, true ); 
		bool isL3   = obj.hasPathName( pathNamesAll[h], false, true ); 
		bool isLF   = obj.hasPathName( pathNamesAll[h], true, false ); 
		bool isNone = obj.hasPathName( pathNamesAll[h], false, false ); 
		std::cout << "   " << pathNamesAll[h];
		if (isBoth) std::cout << "(L,3)";
		if (isL3 && !isBoth) std::cout << "(*,3)";
		if (isLF && !isBoth) std::cout << "(L,*)";
		if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
	    }
	    std::cout << std::endl;
        }
        if (useTriggerList_) {
            for(unsigned int i = 0; i < pathNamesAll.size(); i++) {
	        if (std::find(triggerList_.begin(), triggerList_.end(), pathNamesAll[i]) != triggerList_.end()) {
                    triggerObject->push_back(obj); }
            }
        } else {
           triggerObject->push_back(obj);
        }
    }
    if (printAll_) std::cout << std::endl;
    }

    // Fill Collections
    iEvent.put( triggerName, "triggerNames" );
    iEvent.put( triggerBit, "triggerBits" );
    iEvent.put( triggerPrescale, "triggerPrescales" );
    iEvent.put(triggerObject, "triggerObjects");
    return true;
}

b2g_miniAodAnalyzer_trigger::~b2g_miniAodAnalyzer_trigger(){}

//define this as a plug-in
DEFINE_FWK_MODULE(b2g_miniAodAnalyzer_trigger);
