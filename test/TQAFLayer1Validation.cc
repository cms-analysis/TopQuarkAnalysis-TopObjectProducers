// -*- C++ -*-
//
// Package:    TQAFLayer1Validation
// Class:      TQAFLayer1Validation
// 
/**\class TQAFLayer1Validation TQAFLayer1Validation.cc TopQuarkAnalysis/TQAFLayer1Validation/src/TQAFLayer1Validation.cc

 Description: validatation for TQAF Layer1 code

 Implementation:
 List of checks:
 1) verify no duplicate electrons in TopElectron collection
 2) calculate reco efficiency from TopElectrons
 3) check that all TopElectrons are in the original PMGsfElectron collection, and that they have the same electronId
 4) the fraction of electrons that pass the electron id
 5) the fraction of electrons that, having passed electron id, pass some specified cuts on isolation
 NB obviously the results of some of these checks are sample-dependent
*/
//
// Original Author:  James LAMB
//         Created:  Thu Oct 25 17:44:42 CEST 2007
// $Id$
//
//


// system include files
#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>


// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"



#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/PMGsfElectronIsoCollection.h"
#include "DataFormats/EgammaCandidates/interface/PMGsfElectronIsoNumCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidateFwd.h"

#include "AnalysisDataFormats/Egamma/interface/ElectronID.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronIDFwd.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronIDAssociation.h"
#include "AnalysisDataFormats/TopObjects/interface/TopElectron.h"

#include "TH1F.h"
#include "TFile.h"


class TQAFLayer1Validation : public edm::EDAnalyzer {
public:
  explicit TQAFLayer1Validation(const edm::ParameterSet&);
  ~TQAFLayer1Validation();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool findDuplicates(edm::Handle<std::vector<TopElectron> > topElectrons);
  std::vector<reco::GenParticleCandidate> getGenPrimEles(edm::Handle<reco::CandidateCollection> genParHandle);
  std::vector<uint32_t> findNonDupeIndices(edm::Handle<std::vector<reco::PixelMatchGsfElectron> > PMGsfElectronsH);

  //std::vector<reco::Candidate> getGenPrimEles(edm::Handle<reco::CandidateCollection> genParHandle);
  std::pair<float,size_t> calcMinDR(reco::Particle par, std::vector<TopElectron> eles);  
  void initHistos();
  void saveHistos();
  
  
  // ----------member data ---------------------------

  //for electrons
  uint32_t nEventsWithDuplicates_;  
  uint32_t nEleMisMatchId_;
  uint32_t nTopElectrons_;
  uint32_t nTopElectronsPassEleId_;
  uint32_t nTopElectronsPassEleIdPassTqaf_tkIso_;
  uint32_t nTopElectronsPassEleIdPassTqaf_caloIso_;
  uint32_t nPMGsfElectrons_;
  uint32_t nPMGsfElectronsPassEleId_;

  TH1F *allGenEleVsEt_;
  TH1F *allGenEleVsEta_;
  TH1F *matchedGenEleVsEt_;
  TH1F *matchedGenEleVsEta_;


  
  //utilities
  uint32_t nEvt_;
  TH1F *theNumbers_;
  std::string outputFileName_;
  TFile *outputRootFile_;

  

};


TQAFLayer1Validation::TQAFLayer1Validation(const edm::ParameterSet& iConfig) {

  nEventsWithDuplicates_=0;
  nEleMisMatchId_=0;
  nTopElectrons_=0;
  nTopElectronsPassEleId_=0;
  nTopElectronsPassEleIdPassTqaf_tkIso_=0;
  nTopElectronsPassEleIdPassTqaf_caloIso_=0;
  nPMGsfElectrons_=0;
  nPMGsfElectronsPassEleId_=0;

  outputFileName_=iConfig.getUntrackedParameter<std::string>("outputFileName");
  nEvt_=0;
  initHistos();
  
}


TQAFLayer1Validation::~TQAFLayer1Validation() {

}


void
TQAFLayer1Validation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  using namespace edm;
  if (!(nEvt_%1000)) {
    std::cout <<"doing event: "<<nEvt_<<std::endl;
  }
  nEvt_++;
  /* get the necessary event data, some minor re-containerizing */
  Handle<std::vector<TopElectron> > topElectronsH;
  iEvent.getByLabel("allLayer1TopElectrons","",topElectronsH);
  Handle<reco::CandidateCollection> genParHandle;
  iEvent.getByLabel("genParticleCandidates",genParHandle);
  Handle<std::vector<reco::PixelMatchGsfElectron> > PMGsfElectronsH;
  iEvent.getByLabel("pixelMatchGsfElectrons",PMGsfElectronsH);
  Handle<reco::ElectronIDCollection> eidHandle;
  Handle<reco::ElectronIDAssociationCollection> eidAHandle;
  iEvent.getByLabel("electronId","",eidHandle);
  iEvent.getByLabel("electronId","",eidAHandle);
  /* end get the necessary event data*/

  //task 1: check for any duplicates in the TopElectron collection
  if (findDuplicates(topElectronsH)) {
    nEventsWithDuplicates_++;
  }
  
  //task 2: efficiency for electrons
  std::vector<reco::GenParticleCandidate> genPrimEles=getGenPrimEles(genParHandle);
  //std::vector<reco::Candidate> genPrimEles=getGenPrimEles(genParHandle);
  for (uint32_t ip=0;ip<genPrimEles.size();ip++) {
    
    allGenEleVsEt_->Fill(genPrimEles[ip].et());
    allGenEleVsEta_->Fill(genPrimEles[ip].eta());
    
    //find the nearest reco ele
    std::pair<float, size_t> match=calcMinDR(genPrimEles[ip], *topElectronsH);  
    if (match.first<0.2) {//criterion for matching is deltaR<0.2 - pretty loose
      matchedGenEleVsEt_->Fill(genPrimEles[ip].et());
      matchedGenEleVsEta_->Fill(genPrimEles[ip].eta());
    }
  }
  
  //task 3: check that all TopElectrons are found in the original collection, and that their 
  //electronID results are the same
  for (uint32_t ie=0;ie<topElectronsH->size();ie++) {
    reco::GsfTrackRef refTrack=topElectronsH->at(ie).gsfTrack();
    reco::SuperClusterRef refSC=topElectronsH->at(ie).superCluster();
    double topEleID=topElectronsH->at(ie).getLeptonID();
    bool isMatched=false;
    for (uint32_t je=0;je<PMGsfElectronsH->size();je++) {
      reco::GsfTrackRef cmpTrack=PMGsfElectronsH->at(je).gsfTrack();
      reco::SuperClusterRef cmpSC=PMGsfElectronsH->at(je).superCluster();
      if (refTrack==cmpTrack && refSC==cmpSC) {
	isMatched=true;
	edm::Ref<std::vector<reco::PixelMatchGsfElectron> > elecsRef( PMGsfElectronsH, je );
	reco::ElectronIDAssociationCollection::const_iterator elecID = eidAHandle->find(elecsRef);
	const reco::ElectronIDRef& id = elecID->val;
	double gsfEleID=id->cutBasedDecision();
	if (gsfEleID!=topEleID) {
	  nEleMisMatchId_++;
	}
      }
    }
  }
  
  //task 4: check fractions of events passing electron id
  nTopElectrons_+=topElectronsH->size();
  for (uint32_t ie=0;ie<topElectronsH->size();ie++) {
    double topEleID=topElectronsH->at(ie).getLeptonID();
    if (topEleID) {
      nTopElectronsPassEleId_++;
    }
  }
  
  //task 5: check fractions of electrons that, having passed electron id, 
  // pass specified electron cuts
  float tqaf_tkIsoCut=3;
  float tqaf_caloIsoCut=6;
  //stick with just these two for now, add egamma-produced isolation soon
  for (uint32_t ie=0;ie<topElectronsH->size();ie++) {
    double topEleID=topElectronsH->at(ie).getLeptonID();
    if (topEleID) {
      float tqaf_tkIso=topElectronsH->at(ie).getTrackIso();
      float tqaf_caloIso=topElectronsH->at(ie).getCaloIso();
      if (tqaf_tkIso<tqaf_tkIsoCut) nTopElectronsPassEleIdPassTqaf_tkIso_++;
      if (tqaf_caloIso<tqaf_caloIsoCut) nTopElectronsPassEleIdPassTqaf_caloIso_++;
    }
    
  }
  
  //task 6: compare numbers of PMGsfElectrons and TopElectrons, not counting the duplicates
  // in the PMGsfElectron collection
   std::vector<uint32_t> goodIndices=findNonDupeIndices(PMGsfElectronsH);
   for (uint32_t ie=0;ie<PMGsfElectronsH->size();ie++) {
     if (find(goodIndices.begin(),goodIndices.end(),ie)==goodIndices.end()) continue;
     nPMGsfElectrons_++;
     edm::Ref<std::vector<reco::PixelMatchGsfElectron> > elecsRef( PMGsfElectronsH, ie );
     reco::ElectronIDAssociationCollection::const_iterator elecID = eidAHandle->find(elecsRef);
     const reco::ElectronIDRef& id = elecID->val;
     double gsfEleID=id->cutBasedDecision();
     if (gsfEleID) {
       nPMGsfElectronsPassEleId_++;
     }
   }
   

   
}


void TQAFLayer1Validation::beginJob(const edm::EventSetup&) {
}

void TQAFLayer1Validation::endJob() {
  
  theNumbers_->SetBinContent(1,nEventsWithDuplicates_);
  theNumbers_->SetBinContent(2,nEleMisMatchId_);
  theNumbers_->SetBinContent(3,nTopElectrons_);
  theNumbers_->SetBinContent(4,nTopElectronsPassEleId_);
  theNumbers_->SetBinContent(5,nTopElectronsPassEleIdPassTqaf_tkIso_);
  theNumbers_->SetBinContent(6,nTopElectronsPassEleIdPassTqaf_caloIso_);
  theNumbers_->SetBinContent(7,nPMGsfElectrons_);
  theNumbers_->SetBinContent(8,nPMGsfElectronsPassEleId_);


  saveHistos();

}

//returns true if the vector contains duplicate electrons
//duplicate electrons, more precisely said, are two (or more) electrons with the same 
//supercluster or track
bool TQAFLayer1Validation::findDuplicates(edm::Handle<std::vector<TopElectron> > topElectrons) {
  
  bool return_value=false;

  for (uint32_t ie=0;ie<topElectrons->size();ie++) {
    reco::GsfTrackRef refTrack=topElectrons->at(ie).gsfTrack();
    reco::SuperClusterRef refSC=topElectrons->at(ie).superCluster();
    for (uint32_t je=ie+1;je<topElectrons->size();je++) {
      reco::GsfTrackRef cmpTrack=topElectrons->at(je).gsfTrack();
      reco::SuperClusterRef cmpSC=topElectrons->at(je).superCluster();
      if (refTrack==cmpTrack || refSC==cmpSC) return_value=true;
    }
  }
  
  return return_value;
  
}

//get a vector of the primary electrons out of the event
//that means the electrons that come from the W
//in practice that means the electrons that come from the electrons that come from the W
// (so after any FSR)
std::vector<reco::GenParticleCandidate> TQAFLayer1Validation::getGenPrimEles(edm::Handle<reco::CandidateCollection> genParHandle) {
  
  std::vector<reco::GenParticleCandidate> output;
  for (uint32_t ip=0;ip<genParHandle->size();ip++) {
    if (abs((*genParHandle)[ip].pdgId())!=11) continue;
    if ((*genParHandle)[ip].status()!=1) continue;
    if ((*genParHandle)[ip].numberOfMothers()<1) continue;
    if (abs((*genParHandle)[ip].mother()->pdgId())!=11) continue;
    if ((*genParHandle)[ip].mother()->numberOfMothers()<1) continue;
    if (abs((*genParHandle)[ip].mother()->mother()->pdgId())!=24) continue;
    //yikes, now we finally have electrons whose mothers are electrons whose mothers are Ws
    //reco::GenParticleCandidate tmp=static_cast<reco::GenParticleCandidate> ((*genParHandle)[ip]);
    //output.push_back((*genParHandle)[ip]);
    //reco::GenParticleCandidate tmp=(reco::GenParticleCandidate) ((*genParHandle)[ip]);
    reco::GenParticleCandidate tmp = *(dynamic_cast<reco::GenParticleCandidate *>(const_cast<reco::Candidate *>(&(*genParHandle)[ip])));
    output.push_back(tmp);
  }
  return output;

}

//return indices of PMGsfElectrons which are _not_ duplicates
//in the case where duplicates exist in the original collection, the ones with best E/P are selected
std::vector<uint32_t> TQAFLayer1Validation::findNonDupeIndices(edm::Handle<std::vector< reco::PixelMatchGsfElectron> > PMGsfElectronsH) {
  
  std::vector<uint32_t> goodIndices;
  for (uint32_t ie=0;ie<PMGsfElectronsH->size();ie++) {
    goodIndices.push_back(ie);
  }
  for (uint32_t ie=0;ie<PMGsfElectronsH->size();ie++) {
    //skip this if it is already removed from the vector
    if (std::find(goodIndices.begin(),goodIndices.end(),ie)==goodIndices.end()) continue;
    reco::GsfTrackRef refTrack=PMGsfElectronsH->at(ie).gsfTrack();
    reco::SuperClusterRef refSC=PMGsfElectronsH->at(ie).superCluster();
    for (uint32_t je=0;je<PMGsfElectronsH->size();je++) {
      if (je==ie) continue;
      if (std::find(goodIndices.begin(),goodIndices.end(),je)==goodIndices.end()) continue;
      reco::GsfTrackRef cmpTrack=PMGsfElectronsH->at(je).gsfTrack();
      reco::SuperClusterRef cmpSC=PMGsfElectronsH->at(je).superCluster();
      if (refTrack==cmpTrack||refSC==cmpSC) {
	//have match, now remove one from the vector of good indices
	float diff1=fabs(PMGsfElectronsH->at(ie).eSuperClusterOverP()-1);
	float diff2=fabs(PMGsfElectronsH->at(je).eSuperClusterOverP()-1);
	if (diff1>diff2) {
	  if (std::find(goodIndices.begin(),goodIndices.end(),ie)!=goodIndices.end()) {
	    goodIndices.erase(std::find(goodIndices.begin(),goodIndices.end(),ie));
	  }
	} else {
	  if (std::find(goodIndices.begin(),goodIndices.end(),je)!=goodIndices.end()) {
	    goodIndices.erase(std::find(goodIndices.begin(),goodIndices.end(),je));
	  }
	}
      }
    }
  }

  return goodIndices;
}


std::pair<float,size_t> TQAFLayer1Validation::calcMinDR(reco::Particle par, std::vector<TopElectron> eles) {

  float mindr=9999.;
  size_t indexMin=0;
  
  for (size_t ie=0;ie<eles.size();ie++) {
    
    float dr=deltaR(par,eles[ie]);
    if (dr<mindr) {
      mindr=dr;
      indexMin=ie;
    }
    
  }
  return std::make_pair(mindr,indexMin);
}


void TQAFLayer1Validation::initHistos() {

  theNumbers_=new TH1F("theNumbers_","Misc. Quantities",10,0,10);
  allGenEleVsEt_=new TH1F("allGenEleVsEt_","Et Spectrum of Gen-level electrons",500,0,500);
  allGenEleVsEta_=new TH1F("allGenEleVsEta_","Eta Spectrum of Gen-level electrons",100,-5,5);
  matchedGenEleVsEt_=new TH1F("matchedGenEleVsEt_","Et Spectrum of Gen-level electrons found in Reco",500,0,500);
  matchedGenEleVsEta_=new TH1F("matchedGenEleVsEta_","Eta Spectrum of Gen-level electrons found in Reco",100,-5,5);

  theNumbers_->SetDirectory(0);
  allGenEleVsEt_->SetDirectory(0);
  allGenEleVsEta_->SetDirectory(0);
  matchedGenEleVsEt_->SetDirectory(0);
  matchedGenEleVsEta_->SetDirectory(0);

  

}

void TQAFLayer1Validation::saveHistos() {

  outputRootFile_=new TFile(outputFileName_.c_str(),"RECREATE");

  theNumbers_->SetDirectory(outputRootFile_);
  allGenEleVsEt_->SetDirectory(outputRootFile_);
  allGenEleVsEta_->SetDirectory(outputRootFile_);
  matchedGenEleVsEt_->SetDirectory(outputRootFile_);
  matchedGenEleVsEta_->SetDirectory(outputRootFile_);

  outputRootFile_->Write();
  outputRootFile_->Close();
  delete outputRootFile_;

}


//define this as a plug-in
DEFINE_FWK_MODULE(TQAFLayer1Validation);
