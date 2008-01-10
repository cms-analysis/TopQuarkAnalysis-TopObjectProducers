//
// $Id: TopJetProducer.cc,v 1.36.2.2 2007/11/25 19:03:40 lowette Exp $
//

#include "TopQuarkAnalysis/TopObjectProducers/interface/TopJetProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/TrackProbabilityTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackProbabilityTagInfoFwd.h"
#include "DataFormats/BTauReco/interface/TrackCountingTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackCountingTagInfoFwd.h"
#include "DataFormats/BTauReco/interface/SoftLeptonTagInfo.h"
#include "DataFormats/BTauReco/interface/SoftLeptonTagInfoFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "PhysicsTools/Utilities/interface/DeltaR.h"

#include "TopQuarkAnalysis/TopObjectResolutions/interface/TopObjectResolutionCalc.h"

#include <vector>
#include <memory>


//
// constructors and destructor
//

TopJetProducer::TopJetProducer(const edm::ParameterSet& iConfig) {
  // initialize the configurables
  jetsSrc_                 = iConfig.getParameter<edm::InputTag>            ( "jetSource" );
  // TEMP Jet cleaning from electrons
  doJetCleaning_           = iConfig.getParameter<bool>                     ( "doJetCleaning" );
  topElectronsLabel_       = iConfig.getParameter<edm::InputTag>            ( "topElectronsInput" );
  topMuonsLabel_           = iConfig.getParameter<edm::InputTag>            ( "topMuonsInput" );
  // TEMP End
  getJetMCFlavour_         = iConfig.getParameter<bool>                     ( "getJetMCFlavour" );
  jetPartonMapSource_      = iConfig.getParameter<edm::InputTag>            ( "JetPartonMapSource" );
  addGenPartonMatch_       = iConfig.getParameter<bool>                     ( "addGenPartonMatch" );
  genPartonSrc_            = iConfig.getParameter<edm::InputTag>            ( "genPartonSource" );
  addGenJetMatch_          = iConfig.getParameter<bool>                     ( "addGenJetMatch" );
  genJetSrc_               = iConfig.getParameter<edm::InputTag>            ( "genJetSource" );
  addPartonJetMatch_       = iConfig.getParameter<bool>                     ( "addPartonJetMatch" );
  partonJetSrc_            = iConfig.getParameter<edm::InputTag>            ( "partonJetSource" );
  addResolutions_          = iConfig.getParameter<bool>                     ( "addResolutions" );
  useNNReso_               = iConfig.getParameter<bool>                     ( "useNNResolutions" );
  caliJetResoFile_         = iConfig.getParameter<std::string>              ( "caliJetResoFile" );
  caliBJetResoFile_        = iConfig.getParameter<std::string>              ( "caliBJetResoFile" );
  addBTagInfo_             = iConfig.getParameter<bool>                     ( "addBTagInfo" );
  addDiscriminators_       = iConfig.getParameter<bool>                     ( "addDiscriminators" );
  addJetTagRefs_           = iConfig.getParameter<bool>                     ( "addJetTagRefs" );
  tagModuleLabelsToKeep_   = iConfig.getParameter<std::vector<std::string> >( "tagModuleLabelsToKeep" );
  addAssociatedTracks_     = iConfig.getParameter<bool>                     ( "addAssociatedTracks" ); 
  trackAssociationPSet_    = iConfig.getParameter<edm::ParameterSet>        ( "trackAssociation" );
  addJetCharge_            = iConfig.getParameter<bool>                     ( "addJetCharge" ); 
  jetChargePSet_           = iConfig.getParameter<edm::ParameterSet>        ( "jetCharge" );

  // TEMP Jet cleaning from electrons
  LEPJETDR_=0.3;//deltaR cut used to associate a jet to an electron for jet cleaning.  Make it configurable?
  ELEISOCUT_=2.0;//cut on electron isolation for jet cleaning, because Jim says so
  MUISOCUT_=2.0;//cut on muon isolation for jet cleaning
  // TEMP End
    
  // construct resolution calculator
  if (addResolutions_) {
    theResoCalc_ = new TopObjectResolutionCalc(edm::FileInPath(caliJetResoFile_).fullPath(), useNNReso_);
    theBResoCalc_ = new TopObjectResolutionCalc(edm::FileInPath(caliBJetResoFile_).fullPath(), useNNReso_);
  }

  // construct Jet Track Associator
  simpleJetTrackAssociator_ = helper::SimpleJetTrackAssociator(trackAssociationPSet_);
  // construct Jet Charge Computer
  if (addJetCharge_) jetCharge_ = new JetCharge(jetChargePSet_);
 
  // produces vector of jets
  produces<std::vector<TopJet> >();
}


TopJetProducer::~TopJetProducer() {
  if (addResolutions_) {
    delete theResoCalc_;
    delete theBResoCalc_;
  }
  if (addJetCharge_) delete jetCharge_;
}


//
// member functions
//

void TopJetProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // Get the vector of jets
  edm::Handle<std::vector<TopJetType> > jets;
  iEvent.getByLabel(jetsSrc_, jets);
  // TEMP Jet cleaning from electrons
  edm::Handle<std::vector<TopElectron> > electronsHandle;
  iEvent.getByLabel(topElectronsLabel_, electronsHandle);
  std::vector<TopElectron> electrons=*electronsHandle;
  edm::Handle<std::vector<TopMuon> > muonsHandle;
  iEvent.getByLabel(topMuonsLabel_, muonsHandle);
  std::vector<TopMuon> muons=*muonsHandle;
  // TEMP End

  if (doJetCleaning_) {
    // TEMP Jet cleaning from electrons
    //select isolated leptons to remove from jets collection
    electrons=selectIsolated(electrons,ELEISOCUT_);
    muons=selectIsolated(muons,MUISOCUT_);
    // TEMP End
  }

  // for jet flavour
  edm::Handle<reco::CandMatchMap> JetPartonMap;
  if (getJetMCFlavour_) iEvent.getByLabel (jetPartonMapSource_, JetPartonMap);

  // Get the vector of generated particles from the event if needed
  edm::Handle<reco::CandidateCollection> particles;
  if (addGenPartonMatch_) iEvent.getByLabel(genPartonSrc_, particles);
  // Get the vector of GenJets from the event if needed
  edm::Handle<reco::GenJetCollection> genJets;
  if (addGenJetMatch_) iEvent.getByLabel(genJetSrc_, genJets);
/* TO BE IMPLEMENTED FOR >= 1_5_X
  // Get the vector of PartonJets from the event if needed
  edm::Handle<reco::SomePartonJetType> particles;
  if (addPartonJetMatch_) iEvent.getByLabel(partonJetSrc_, partonJets);
*/

  // Get the vector of jet tags with b-tagging info
  std::vector<edm::Handle<std::vector<reco::JetTag> > > jetTags_testManyByType ;
  iEvent.getManyByType(jetTags_testManyByType);
  // Define the handles for the specific algorithms
  edm::Handle<reco::SoftLeptonTagInfoCollection> jetsInfoHandle_sl;
  edm::Handle<reco::TrackProbabilityTagInfoCollection> jetsInfoHandleTP;
  edm::Handle<reco::TrackCountingTagInfoCollection> jetsInfoHandleTC;

  // tracks Jet Track Association, by hand in CMSSW_1_3_X
  edm::Handle<reco::TrackCollection> hTracks;
  iEvent.getByLabel(trackAssociationPSet_.getParameter<edm::InputTag>("tracksSource"), hTracks);


  // loop over jets
  std::vector<TopJet> * topJets = new std::vector<TopJet>(); 
  for (size_t j = 0; j < jets->size(); j++) {

    if (doJetCleaning_) {
    // TEMP Jet cleaning from electrons
      //check that the jet doesn't match in deltaR with an isolated lepton
      //if it does, then it needs to be cleaned (ie don't put it in the TopJet collection)
      //FIXME: don't do muons until have a sensible cut value on their isolation
      float mindr=9999.;
      for (size_t ie=0; ie<electrons.size(); ie++) {
        float dr=DeltaR<reco::Candidate>()((*jets)[j],electrons[ie]);
        if (dr<mindr) {
          mindr=dr;
        }
      }
      //if the jet is closely matched in dR to electron, skip it
      if (mindr<LEPJETDR_) {
        continue;
      }
    // TEMP End
    }

    // define the jet correctors
    const JetCorrector * defaultJetCorr = JetCorrector::getJetCorrector("MCJetCorrectorIcone5", iSetup);
    const JetCorrector * udsJetCorr     = JetCorrector::getJetCorrector("L5FlavorJetCorrectorUds", iSetup);
    const JetCorrector * gluJetCorr     = JetCorrector::getJetCorrector("L5FlavorJetCorrectorGluon", iSetup);
    const JetCorrector * cJetCorr       = JetCorrector::getJetCorrector("L5FlavorJetCorrectorC", iSetup);
    const JetCorrector * bJetCorr       = JetCorrector::getJetCorrector("L5FlavorJetCorrectorB", iSetup);
    // calculate the energy correction factors
    double scaleDefault = defaultJetCorr->correction((*jets)[j]);
    double scaleUds     = scaleDefault * udsJetCorr->correction((*jets)[j]);
    double scaleGlu     = scaleDefault * gluJetCorr->correction((*jets)[j]);
    double scaleC       = scaleDefault * cJetCorr->correction((*jets)[j]);
    double scaleB       = scaleDefault * bJetCorr->correction((*jets)[j]);

    // construct the TopJet
    TopJet ajet((*jets)[j]);
    ajet.setP4(scaleDefault*(*jets)[j].p4());
    ajet.setScaleCalibFactors(1./scaleDefault, scaleUds, scaleGlu, scaleC, scaleB);

    // get the MC flavour information for this jet
    if (getJetMCFlavour_) {
      for (reco::CandMatchMap::const_iterator f = JetPartonMap->begin(); f != JetPartonMap->end(); f++) {
        const reco::Candidate * jetClone = f->key->masterClone().get();
        // if (jetClone == &((*jets)[j])) { // comparison by address doesn't work
        if (fabs(jetClone->eta() - (*jets)[j].eta()) < 0.001 &&
            fabs(jetClone->phi() - (*jets)[j].phi()) < 0.001) {
          ajet.setPartonFlavour(f->val->pdgId());
        }
      }
    }
    // do the parton matching
    if (addGenPartonMatch_) {
      // initialize best match as null
      reco::GenParticleCandidate bestParton(0, reco::Particle::LorentzVector(0, 0, 0, 0), reco::Particle::Point(0,0,0), 0, 0, true);
      float bestDR = 0;
      // find the closest parton
      for (reco::CandidateCollection::const_iterator itParton = particles->begin(); itParton != particles->end(); ++itParton) {
        reco::GenParticleCandidate aParton = *(dynamic_cast<reco::GenParticleCandidate *>(const_cast<reco::Candidate *>(&*itParton)));
        if (aParton.status()==3 &&
            (abs(aParton.pdgId())==1 || abs(aParton.pdgId())==2 ||
             abs(aParton.pdgId())==3 || abs(aParton.pdgId())==4 ||
             abs(aParton.pdgId())==5 || abs(aParton.pdgId())==21)) {
          float currDR = DeltaR<reco::Candidate>()(aParton, ajet);
          // matching with hard-cut at 0.4
          // can be improved a la muon-electron, such that each parton
          // maximally matches 1 jet, but this requires two loops
          if (bestDR == 0 || (currDR < bestDR || currDR < 0.4)) {
            bestParton = aParton;
            bestDR = currDR;
          }
        }
      }
      ajet.setGenParton(bestParton);
    }
    // do the GenJet matching
    if (addGenJetMatch_) {
      // initialize best match as null
//NEED TO INITIALIZE TO ZERO
      reco::GenJet bestGenJet;//0, reco::Particle::LorentzVector(0, 0, 0, 0), reco::Particle::Point(0,0,0), 0, 0);
      float bestDR = 0;
      // find the closest parton
      for (reco::GenJetCollection::const_iterator itGenJet = genJets->begin(); itGenJet != genJets->end(); ++itGenJet) {
// do we need some criteria?      if (itGenJet->status()==3) {
          float currDR = DeltaR<reco::Candidate>()(*itGenJet, ajet);
          // matching with hard-cut at 0.4
          // can be improved a la muon-electron, such that each genjet
          // maximally matches 1 jet, but this requires two loops
          if (bestDR == 0 || (currDR < bestDR || currDR < 0.4)) {
            bestGenJet = *itGenJet;
            bestDR = currDR;
          }
//          }
      }
      ajet.setGenJet(bestGenJet);
    }
    // TO BE IMPLEMENTED FOR >=1_5_X: do the PartonJet matching
    if (addPartonJetMatch_) {
    }

    // add resolution info if demanded
    if (addResolutions_) {
      (*theResoCalc_)(ajet);
      TopJet abjet(ajet.getBCorrJet());
      (*theBResoCalc_)(abjet);
      ajet.setBResolutions(abjet.getResET(), abjet.getResEta(), abjet.getResPhi(), abjet.getResA(), abjet.getResB(), abjet.getResC(), abjet.getResD(), abjet.getResTheta());
    }

    // add b-tag info if available & required
    if (addBTagInfo_) {
      for (size_t k=0; k<jetTags_testManyByType.size(); k++) {
        edm::Handle<std::vector<reco::JetTag> > jetTags = jetTags_testManyByType[k];

        //get label and module names
        std::string moduleLabel = (jetTags).provenance()->moduleLabel();
	
	
	//look only at the tagger present in tagModuleLabelsToKeep_
	for (unsigned int i = 0; i < tagModuleLabelsToKeep_.size(); ++i) {
	  if (moduleLabel == tagModuleLabelsToKeep_[i]) {
	    for (size_t t = 0; t < jetTags->size(); t++) {
	      edm::RefToBase<reco::Jet> jet_p = (*jetTags)[t].jet();
	      if (jet_p.isNull()) {
		/*std::cout << "-----------> JetTag::jet() returned null reference" << std::endl; */
		continue;
	      }
	      if (DeltaR<reco::Candidate>()( (*jets)[j], *jet_p ) < 0.00001) {
		//********store discriminators*********
		if (addDiscriminators_) {
		  //look only at the tagger present in tagModuleLabelsToKeep_
                  std::pair<std::string, double> pairDiscri;
                  pairDiscri.first = moduleLabel;
                  pairDiscri.second = (*jetTags)[t].discriminator();
                  ajet.addBDiscriminatorPair(pairDiscri);
                  continue;
                }
              }
	      
	      //********store jetTagRef*********
	      if (addJetTagRefs_) {
		std::pair<std::string, reco::JetTagRef> pairjettagref;
		pairjettagref.first = moduleLabel;
		pairjettagref.second = reco::JetTagRef(jetTags, t);
		ajet.addBJetTagRefPair(pairjettagref);
	      }
            }
          }
        }
      }
    }
    
    // Associate tracks with jet (at least temporary)
    simpleJetTrackAssociator_.associate(ajet.momentum(), hTracks, ajet.associatedTracks_);
    
    // PUT HERE EVERYTHING WHICH NEEDS TRACKS
    if (addJetCharge_) {
      ajet.jetCharge_ = static_cast<float>(jetCharge_->charge(ajet.p4(), ajet.associatedTracks_));
    }

    // drop jet track association if the user does not want it
    if (!addAssociatedTracks_) ajet.associatedTracks_.clear();

    // end of TopObjectProducer loop
    topJets->push_back(ajet);
  }

  // sort jets in ET
  std::sort(topJets->begin(), topJets->end(), eTComparator_);

  // put genEvt  in Event
  std::auto_ptr<std::vector<TopJet> > myTopJetProducer(topJets);
  iEvent.put(myTopJetProducer);

}


// TEMP Jet cleaning from electrons
//takes a vector of electrons and returns a vector that only contains the ones that are isolated
//The second argument is the isolation cut to use
std::vector<TopElectron> TopJetProducer::selectIsolated(const std::vector<TopElectron> &electrons, float isoCut) {
  std::vector<TopElectron> output;
  for (size_t ie=0; ie<electrons.size(); ie++) {
    if (electrons[ie].getTrackIso() < isoCut) {
      output.push_back(electrons[ie]);
    }
  }
  return output;
}
// TEMP End


// TEMP Jet cleaning from electrons
//takes a vector of muons and returns a vector that only contains the ones that are isolated
//The second argument is the isolation cut to use
//FIXME I could combine this with the one for electrons using templates?
std::vector<TopMuon> TopJetProducer::selectIsolated(const std::vector<TopMuon> &muons, float isoCut) {
  std::vector<TopMuon> output;
  for (size_t iu=0; iu<muons.size(); iu++) {
    if (muons[iu].getTrackIso() < isoCut) {
      output.push_back(muons[iu]);
    }
  }
  return output;
}
// TEMP End


