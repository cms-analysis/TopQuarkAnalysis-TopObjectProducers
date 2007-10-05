//
// $Id: TopMuonProducer.cc,v 1.7.2.5 2007/10/02 16:55:23 lowette Exp $
//

#include "TopQuarkAnalysis/TopObjectProducers/interface/TopMuonProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "PhysicsTools/Utilities/interface/DeltaR.h"

#include "TopQuarkAnalysis/TopObjectResolutions/interface/TopObjectResolutionCalc.h"
#include "TopQuarkAnalysis/TopLeptonSelection/interface/TrackerIsolationPt.h"
#include "TopQuarkAnalysis/TopLeptonSelection/interface/CaloIsolationEnergy.h"
#include "TopQuarkAnalysis/TopLeptonSelection/interface/TopLeptonLRCalc.h"

#include "TMath.h"

#include <vector>
#include <memory>


TopMuonProducer::TopMuonProducer(const edm::ParameterSet & iConfig) {

  // general configurables
  muonSrc_       = iConfig.getParameter<edm::InputTag>( "muonSource" );
  // MC matching configurables
  addGenMatch_    = iConfig.getParameter<bool>        ( "addGenMatch" );
  genPartSrc_    = iConfig.getParameter<edm::InputTag>( "genParticleSource" );
  maxDeltaR_     = iConfig.getParameter<double>       ( "maxDeltaR" );
  minRecoOnGenEt_= iConfig.getParameter<double>       ( "minRecoOnGenEt" );
  maxRecoOnGenEt_= iConfig.getParameter<double>       ( "maxRecoOnGenEt" );
  // resolution configurables
  addResolutions_= iConfig.getParameter<bool>         ( "addResolutions" );
  useNNReso_     = iConfig.getParameter<bool>         ( "useNNResolutions" );
  muonResoFile_  = iConfig.getParameter<std::string>  ( "muonResoFile" );
  // isolation configurables
  doTrkIso_      = iConfig.getParameter<bool>         ( "doTrkIsolation" );
  tracksSrc_     = iConfig.getParameter<edm::InputTag>( "tracksSource" );
  doCalIso_      = iConfig.getParameter<bool>         ( "doCalIsolation" );
  towerSrc_      = iConfig.getParameter<edm::InputTag>( "towerSource" );
  // muon ID configurables
  addMuonID_     = iConfig.getParameter<bool>         ( "addMuonID" );
//  muonIDSrc_     = iConfig.getParameter<edm::InputTag>( "muonIDSource" );
  // likelihood ratio configurables
  addLRValues_   = iConfig.getParameter<bool>         ( "addLRValues" );
  tracksSrc_     = iConfig.getParameter<edm::InputTag>( "tracksSource" );
  muonLRFile_    = iConfig.getParameter<std::string>  ( "muonLRFile" );

  // construct resolution calculator
  if (addResolutions_) {
    theResoCalc_ = new TopObjectResolutionCalc(edm::FileInPath(muonResoFile_).fullPath(), useNNReso_);
  }

  // produces vector of muons
  produces<std::vector<TopMuon > >();

}


TopMuonProducer::~TopMuonProducer() {
  if (addResolutions_) delete theResoCalc_;
}


void TopMuonProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // Get the collection of muons from the event
  edm::Handle<std::vector<TopMuonType> > muonsHandle;
  iEvent.getByLabel(muonSrc_, muonsHandle);
  std::vector<TopMuonType> muons = *muonsHandle;

  // prepare the MC matching
  edm::Handle<reco::CandidateCollection> particles;
  if (addGenMatch_) {
    iEvent.getByLabel(genPartSrc_, particles);
    matchTruth(*particles, muons);
  }

  // prepare isolation calculation
  edm::Handle<reco::TrackCollection> trackHandle;
  if (doTrkIso_) {
    trkIsolation_= new TrackerIsolationPt();
    iEvent.getByLabel(tracksSrc_, trackHandle);
  }
  std::vector<CaloTower> towers;
  if (doCalIso_) {
    calIsolation_= new CaloIsolationEnergy();
    edm::Handle<CaloTowerCollection> towerHandle;
    iEvent.getByLabel(towerSrc_, towerHandle);
    CaloTowerCollection towerColl = *(towerHandle.product());
    for (CaloTowerCollection::const_iterator itTower = towerColl.begin(); itTower != towerColl.end(); itTower++) {
      towers.push_back(*itTower);
    }
  }

  // prepare LR calculation
  if (addLRValues_) {
    theLeptonLRCalc_ = new TopLeptonLRCalc(iSetup, "", edm::FileInPath(muonLRFile_).fullPath(), "");
  }

  // loop over muons
  std::vector<TopMuon> * topMuons = new std::vector<TopMuon>();
  for (size_t m = 0; m < muons.size(); ++m) {
    // construct the TopMuon
    TopMuon aMuon(muons[m]);
    // match to generated final state muons
    if (addGenMatch_) {
      aMuon.setGenLepton(findTruth(*particles, muons[m]));
    }
    // add resolution info
    if (addResolutions_) {
      (*theResoCalc_)(aMuon);
    }
    // do tracker isolation
    if (doTrkIso_) {
      aMuon.setTrackIso(trkIsolation_->calculate(aMuon, *trackHandle));
    }
    // do calorimeter isolation
    if (doCalIso_) {
      aMuon.setCaloIso(calIsolation_->calculate(aMuon, towers));
    }
    // add muon ID info
    if (addMuonID_) {
      aMuon.setLeptonID((double) TMath::Prob((Double_t) muons[m].combinedMuon()->chi2(), (Int_t) muons[m].combinedMuon()->ndof()));
    }
    // add lepton LR info
    if (addLRValues_) {
      theLeptonLRCalc_->calcLikelihood(aMuon, trackHandle, iEvent);
    }
    // add sel to selected
    topMuons->push_back(TopMuon(aMuon));
  }

  // sort muons in pt
  std::sort(topMuons->begin(), topMuons->end(), pTComparator_);

  // put genEvt object in Event
  std::auto_ptr<std::vector<TopMuon> > ptr(topMuons);
  iEvent.put(ptr);

  if (doTrkIso_ ) delete trkIsolation_;
  if (doCalIso_ ) delete calIsolation_;
  if (addLRValues_) delete theLeptonLRCalc_;

}


reco::GenParticleCandidate TopMuonProducer::findTruth(const reco::CandidateCollection & parts, const TopMuonType & muon) {
  reco::GenParticleCandidate gen;
  for(unsigned int idx=0; idx!=pairGenRecoMuonsVector_.size(); ++idx){
    std::pair<const reco::Candidate*, TopMuonType*> pairGenRecoMuons;
    pairGenRecoMuons = pairGenRecoMuonsVector_[idx];
    double dR = DeltaR<reco::Candidate>()( muon, *(pairGenRecoMuons.second));
    if( !(dR > 0) ){
      gen = *(dynamic_cast<const reco::GenParticleCandidate*>( pairGenRecoMuons.first ) );
    }
  }
  return gen;
}


void TopMuonProducer::matchTruth(const reco::CandidateCollection & parts, std::vector<TopMuonType> & muons) {
  pairGenRecoMuonsVector_.clear();
  reco::CandidateCollection::const_iterator part = parts.begin();
  for( ; part != parts.end(); ++part ){
    reco::GenParticleCandidate gen = *(dynamic_cast<const reco::GenParticleCandidate*>( &(*part)) );
    if( abs(gen.pdgId())==13 && gen.status()==1 ){
      bool   found = false;
      double minDR = 99999;
      TopMuonType* rec = 0;
      TopMuonTypeCollection::iterator muon = muons.begin();
      for ( ; muon !=muons.end(); ++muon){
	double dR = DeltaR<reco::Candidate>()( gen, *muon);
	double ptRecOverGen = muon->pt()/gen.pt();
	if ( ( ptRecOverGen > minRecoOnGenEt_ ) && 
	     ( ptRecOverGen < maxRecoOnGenEt_ ) && 
	     ( dR < maxDeltaR_) ){
	  if ( dR < minDR ){
	    rec = &(*muon);
	    minDR = dR;
	    found = true;
	  }
	}
      }
      if( found == true ){
	pairGenRecoMuonsVector_.push_back( std::pair<const reco::Candidate*,TopMuonType*>(&(*part), rec ) );
      }
    }
  }
}
