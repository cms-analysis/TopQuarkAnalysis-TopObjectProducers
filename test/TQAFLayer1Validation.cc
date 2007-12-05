// -*- C++ -*-
//
// Package:    TQAFLayer1Validation
// Class:      TQAFLayer1Validation
// 
/**\class TQAFLayer1Validation TQAFLayer1Validation.cc TopQuarkAnalysis/TQAFLayer1Validation/src/TQAFLayer1Validation.cc

 Description: validatation for TQAF Layer1 code. 

 Implementation:
 List of checks:
 1) verify no duplicate electrons in TopElectron collection
 2) calculate reco, selection, and isolation efficiency from TopElectrons
 3) check that all TopElectrons are in the original PMGsfElectron collection, and that they have the same electronId
 4) the fraction of electrons that pass the electron id
 5) the fraction of electrons that, having passed electron id, pass some specified cuts on isolation
 6) compare numbers of PMGsfElectrons and TopElectrons, not counting the duplicates in the PMGsfElectron collection
 7) some plots for MET
 8) some plots for jets
 NB obviously the results of some of these checks are sample-dependent
*/
//
// Original Author:  James LAMB
//         Created:  Thu Oct 25 17:44:42 CEST 2007
// $Id: TQAFLayer1Validation.cc,v 1.3 2007/11/26 09:01:27 jlamb Exp $
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
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"

#include "AnalysisDataFormats/Egamma/interface/ElectronID.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronIDFwd.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronIDAssociation.h"
#include "AnalysisDataFormats/TopObjects/interface/TopElectron.h"
#include "AnalysisDataFormats/TopObjects/interface/TopMET.h"
#include "AnalysisDataFormats/TopObjects/interface/TopJet.h"
#include "AnalysisDataFormats/TopObjects/interface/TopMuon.h"

#include "TH1F.h"
#include "TH2F.h"
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
  std::vector<reco::GenParticleCandidate> getGenWDecayProds(edm::Handle<reco::CandidateCollection> genParHandle);
  std::vector<uint32_t> findNonDupeIndices(edm::Handle<std::vector<reco::PixelMatchGsfElectron> > PMGsfElectronsH);
  std::pair<float,size_t> calcMinDR(reco::Particle par, std::vector<TopElectron> eles);  
  std::pair<float,size_t> calcMinDR(reco::Particle par, std::vector<TopJet> jets);
  uint32_t classifyEvent(const reco::CandidateCollection &genPars);
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
  TH1F *recoMatchedGenEleVsEt_;
  TH1F *recoMatchedGenEleVsEta_;
  TH1F *selMatchedGenEleVsEt_;
  TH1F *selMatchedGenEleVsEta_;
  TH1F *isoMatchedGenEleVsEt_;
  TH1F *isoMatchedGenEleVsEta_;
  TH1F *deltaRReco_;
  TH1F *deltaRSel_;
  TH1F *deltaRIso_;
  
  TH1F *corrCaloMET_;
  TH1F *topMET_;
  TH1F *genMET_;
  TH1F *tqafGenMET_;

  TH1F *jetsInvMass_;
  TH1F *jetsEt_;
  TH1F *jetsEta_;
  TH2F* jetsEtVsEta_;
  
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
  Handle<std::vector<TopMET> > topMETHandle;
  iEvent.getByLabel("allLayer1TopMETs",topMETHandle);
  Handle<reco::GenMETCollection> genMETHandle;
  iEvent.getByLabel("genMet",genMETHandle);
  TopMET topMET=(*topMETHandle).at(0);
  reco::GenMET genMET=(*genMETHandle).at(0);
  Handle<std::vector<TopMuon> > muonsHandle;
  iEvent.getByLabel("allLayer1TopMuons",muonsHandle);
  Handle<reco::CaloMETCollection> corrCaloMETHandle;
  iEvent.getByLabel("corMetType1Icone5",corrCaloMETHandle);
  reco::CaloMET corrCaloMET=(*corrCaloMETHandle).at(0);
  Handle<std::vector<TopJet> > jetsHandle;
  iEvent.getByLabel("allLayer1TopJets",jetsHandle);
  /* end get the necessary event data*/
  

  //task 1: check for any duplicates in the TopElectron collection
  if (findDuplicates(topElectronsH)) {
    nEventsWithDuplicates_++;
  }
  
  //task 2: efficiency for electrons
  //select subsets of TopElectrons passing ID cuts, and passing ID cuts and isolation
  std::vector<TopElectron> selEle;
  std::vector<TopElectron> isoEle;
  for (uint32_t ie=0;ie<topElectronsH->size();ie++) {
    //for consistency with fake rate (defined elsewhere), require to be within eta<2.4
    if (fabs(topElectronsH->at(ie).eta()) > 2.4) continue;
    double topEleID=topElectronsH->at(ie).getLeptonID();
    double tkIso=topElectronsH->at(ie).getTrackIso();
    double caloIso=topElectronsH->at(ie).getCaloIso();
    if (topEleID!=0) selEle.push_back(topElectronsH->at(ie));
    if (topEleID!=0 && tkIso<1 && caloIso < 3) isoEle.push_back(topElectronsH->at(ie));
  }
  std::vector<reco::GenParticleCandidate> genPrimEles=getGenPrimEles(genParHandle);
  //std::vector<reco::Candidate> genPrimEles=getGenPrimEles(genParHandle);
  for (uint32_t ip=0;ip<genPrimEles.size();ip++) {
    
    //require the gen ele to be within eta 2.4
    if (fabs(genPrimEles[ip].eta()) > 2.4) continue;

    allGenEleVsEt_->Fill(genPrimEles[ip].et());
    allGenEleVsEta_->Fill(genPrimEles[ip].eta());
    
    //find the nearest reco ele
    std::pair<float, size_t> match=calcMinDR(genPrimEles[ip], *topElectronsH);  
    deltaRReco_->Fill(match.first);
    if (match.first<0.2) {//criterion for matching is deltaR<0.2 - pretty loose
      recoMatchedGenEleVsEt_->Fill(genPrimEles[ip].et());
      recoMatchedGenEleVsEta_->Fill(genPrimEles[ip].eta());
    }
    //find nearest selected ele
    std::pair<float, size_t> matchSel=calcMinDR(genPrimEles[ip], selEle);  
    deltaRSel_->Fill(matchSel.first);
    if (matchSel.first<0.2) {//criterion for matching is deltaR<0.2 - pretty loose
      selMatchedGenEleVsEt_->Fill(genPrimEles[ip].et());
      selMatchedGenEleVsEta_->Fill(genPrimEles[ip].eta());
    }
    //find nearest isolated ele
    std::pair<float, size_t> matchIso=calcMinDR(genPrimEles[ip], isoEle);  
    deltaRIso_->Fill(matchIso.first);
    if (matchIso.first<0.2) {//criterion for matching is deltaR<0.2 - pretty loose
      isoMatchedGenEleVsEt_->Fill(genPrimEles[ip].et());
      isoMatchedGenEleVsEta_->Fill(genPrimEles[ip].eta());
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
   
   //task 7: compare TopMET and GenMET
   topMET_->Fill(topMET.et());
   genMET_->Fill(genMET.et());
   tqafGenMET_->Fill(topMET.getGenMET().et());
   corrCaloMET_->Fill(corrCaloMET.et());

   //task 8: some plots for jets
   //invariant mass of all jets wrt each other
   for (uint32_t ij=0;ij<jetsHandle->size();ij++) {
     reco::Particle::LorentzVector vi=(jetsHandle->at(ij)).p4();
     for (uint32_t jj=ij+1;jj<jetsHandle->size();jj++) {
       reco::Particle::LorentzVector vj=(jetsHandle->at(jj)).p4();
       float mass=(vi+vj).M();
       jetsInvMass_->Fill(mass);
     }
   }
   //jets Et, eta, et vs. eta
   for (uint32_t ij=0;ij<jetsHandle->size();ij++) {
     TopJet jet=jetsHandle->at(ij);
     jetsEt_->Fill(jet.et());
     jetsEta_->Fill(jet.eta());
     jetsEtVsEta_->Fill(jet.eta(),jet.et());
   }
   //W invariant mass
   std::vector<reco::GenParticleCandidate> genWDP=getGenWDecayProds(genParHandle);
   if (genWDP.size()==4) {//in case this code is run on non-top, this block won't run correctly on non-top
     if (abs(genWDP.at(0).pdgId())>0 && abs(genWDP.at(0).pdgId())<6) {//it's a quark
       //therefor W+ decayed as quark
       //find closest jet 
       
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


//get the vector of gen-level particles who's mother's mother's are W's.  So this should be the W decay products after FSR.
//the vector is guaranteed to be ordered such that the decay products of the W+ come before that of the W-
std::vector<reco::GenParticleCandidate> TQAFLayer1Validation::getGenWDecayProds(edm::Handle<reco::CandidateCollection> genParHandle) {
  
  std::vector<reco::GenParticleCandidate> output;
  //first the W+
  for (uint32_t ip=0;ip<genParHandle->size();ip++) {
    if ((*genParHandle)[ip].status()!=1) continue;
    if ((*genParHandle)[ip].numberOfMothers()<1) continue;
    if ((*genParHandle)[ip].mother()->numberOfMothers()<1) continue;
    if ((*genParHandle)[ip].mother()->mother()->pdgId()!=24) continue;

    reco::GenParticleCandidate tmp = *(dynamic_cast<reco::GenParticleCandidate *>(const_cast<reco::Candidate *>(&(*genParHandle)[ip])));
    output.push_back(tmp);
    //break the first loop if got the two W+ decay products
    if (output.size()==2) break;
  }
  //now the W-
  for (uint32_t ip=0;ip<genParHandle->size();ip++) {
    if ((*genParHandle)[ip].status()!=1) continue;
    if ((*genParHandle)[ip].numberOfMothers()<1) continue;
    if ((*genParHandle)[ip].mother()->numberOfMothers()<1) continue;
    if ((*genParHandle)[ip].mother()->mother()->pdgId()!=-24) continue;

    reco::GenParticleCandidate tmp = *(dynamic_cast<reco::GenParticleCandidate *>(const_cast<reco::Candidate *>(&(*genParHandle)[ip])));
    output.push_back(tmp);
    //break the first loop if got the two W- products
    if (output.size()==4) break;
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


std::pair<float,size_t> TQAFLayer1Validation::calcMinDR(reco::Particle par, std::vector<TopJet> jets) {

  float mindr=9999.;
  size_t indexMin=0;
  
  for (size_t ij=0;ij<jets.size();ij++) {
    
    float dr=deltaR(par,jets[ij]);
    if (dr<mindr) {
      mindr=dr;
      indexMin=ij;
    }
    
  }
  return std::make_pair(mindr,indexMin);
}




//this function "classifies" the event at generator level, returning an integer code according to the production
//process (QCD, W, Z, ttbar) and decay mode.  returns 0 for anything it can't figure out.
uint32_t TQAFLayer1Validation::classifyEvent(const reco::CandidateCollection &genPars) {
  
  //first, check for tops and W in the doc lines
  int nTop=0;
  int nW=0;
  int nWEle=0;
  int nWMu=0;
  int nWTau=0;
  int nZ=0;
  int nZEle=0;
  int nZMu=0;
  int nZTau=0;

  
  //loop over the doc lines (assumed not to extend past index 99
  for (uint32_t ip=0;ip<100;ip++) {
    
    //skip anything not a doc line:
    if (genPars[ip].status()!=3) continue;
    
    //count tops and Ws
    if (abs(genPars[ip].pdgId())==6) {
      nTop++;
    }
    if (abs(genPars[ip].pdgId())==24) {
      nW++;
    }
    if (abs(genPars[ip].pdgId())==23) {
      nZ++;
    }
    
    //look at W/Z decay products
    const reco::Candidate *firstMother=genPars[ip].mother(0);
    int firstMotherID=(firstMother!=0)?firstMother->pdgId():0;
    if (abs(firstMotherID)==24 && abs(genPars[ip].pdgId())==11) {
      nWEle++;
    } else if (abs(firstMotherID)==24 && abs(genPars[ip].pdgId())==13) {
      nWMu++;
    } else if (abs(firstMotherID)==24 && abs(genPars[ip].pdgId())==15) {
      nWTau++;
    } else if (abs(firstMotherID)==23 && abs(genPars[ip].pdgId())==11) {
      nZEle++;
    } else if (abs(firstMotherID)==23 && abs(genPars[ip].pdgId())==13) {
      nZMu++;
    } else if (abs(firstMotherID)==23 && abs(genPars[ip].pdgId())==15) {
      nZTau++;
    }
  }
  if (nTop==0 && nW==0 && nZ==0) {//it's qcd
    return 4;
  } else if (nTop==0 && nZ==0 && nW>0) {//it's W
    if (nWEle==1 && nWMu==0 && nWTau==0) {
      return 1;
    } else if (nWEle==0 && nWMu==1 && nWTau==0) {
      return 2;
    } else if (nWEle==0 && nWMu==0 && nWTau==1) {
      return 3;
    } else {
      return 0;//to signal an error
    }
  } else if (nTop==0 && nZ>0 && nW==0) {//it's Z
    if (nZEle==1 && nZMu==0 && nZTau==0) {
      return 15;
    } else if (nZEle==0 && nZMu==1 && nZTau==0) {
      return 16;
    } else if (nZEle==0 && nZMu==0 && nZTau==1) {
      return 17;
    } else {
      return 0;//to signal an error
    }
  } else if (nTop==2) {
    int nLep=nWEle+nWMu+nWTau;
    if (nLep==0) {
      return 14;
    } else if (nLep==1 && nWEle==1) {
      return 5;
    } else if (nLep==1 && nWMu==1) {
      return 6;
    } else if (nLep==1 && nWTau==1) {
      return 7;
    } else if (nLep==2 && nWEle==2) {
      return 8;
    } else if (nLep==2 && nWEle==1 && nWMu==1) {
      return 9;
    } else if (nLep==2 && nWEle==1 && nWTau==1) {
      return 10;
    } else if (nLep==2 && nWMu==2) {
      return 11;
    } else if (nLep==2 && nWMu==1 && nWTau==1) {
      return 12;
    } else if (nLep==2 && nWTau==2) {
      return 13;
    } else {
      return 0;//to signal an error
    }
  } 
  
  //for the compiler warning... 
  //if control reaches here I am surprised
  return 0;
}





void TQAFLayer1Validation::initHistos() {

  theNumbers_=new TH1F("theNumbers_","Misc. Quantities",10,0,10);
  allGenEleVsEt_=new TH1F("allGenEleVsEt_","Et Spectrum of Gen-level electrons",500,0,500);
  allGenEleVsEta_=new TH1F("allGenEleVsEta_","Eta Spectrum of Gen-level electrons",100,-5,5);
  recoMatchedGenEleVsEt_=new TH1F("recoMatchedGenEleVsEt_","Et Spectrum of Gen-level electrons found in Reco",500,0,500);
  recoMatchedGenEleVsEta_=new TH1F("recoMatchedGenEleVsEta_","Eta Spectrum of Gen-level electrons found in Reco",100,-5,5);
  selMatchedGenEleVsEt_=new TH1F("selMatchedGenEleVsEt_","Et Spectrum of Gen-level electrons found in Reco, selected \"loose\"",500,0,500);
  selMatchedGenEleVsEta_=new TH1F("selMatchedGenEleVsEta_","Eta Spectrum of Gen-level electrons found in Reco, selected \"loose\"",100,-5,5);
  isoMatchedGenEleVsEt_=new TH1F("isoMatchedGenEleVsEt_","Et Spectrum of Gen-level electrons found in Reco, selected \"loose\", and isolated",500,0,500);
  isoMatchedGenEleVsEta_=new TH1F("isoMatchedGenEleVsEta_","Eta Spectrum of Gen-level electrons found in Reco, selected \"loose\", and isolated",100,-5,5);
  deltaRReco_=new TH1F("deltaRReco_","deltaR between gen-level electron and nearest reco electron",1000,0,1);
  deltaRSel_=new TH1F("deltaRSel_","deltaR between gen-level electron and nearest \"loose\" selected electron",1000,0,1);
  deltaRIso_=new TH1F("deltaRIso_","deltaR between gen-level electron and nearest \"loose\" selected and isolated electron",1000,0,1);

  topMET_=new TH1F("topMET_","TopMET",200,0,200);
  genMET_=new TH1F("genMET_","GenMET",200,0,200);
  tqafGenMET_=new TH1F("tqafGenMET_","Gen MET from TQAF",200,0,200);
  corrCaloMET_=new TH1F("corrCaloMET_","Corrected Calo MET",200,0,200);

  jetsInvMass_=new TH1F("jetsInvMass_","Invariant Mass of All Pairs of Jets",200,0,200);
  jetsEt_=new TH1F("jetsEt_","TopJets Et",200,0,200);
  jetsEta_=new TH1F("jetsEta_","TopJets Eta",100,-5,5);
  jetsEtVsEta_=new TH2F("jetsEtVsEta_","TopJets Et vs. Eta",100,-5,5,200,0,200);  

  theNumbers_->SetDirectory(0);
  allGenEleVsEt_->SetDirectory(0);
  allGenEleVsEta_->SetDirectory(0);
  recoMatchedGenEleVsEt_->SetDirectory(0);
  recoMatchedGenEleVsEta_->SetDirectory(0);
  selMatchedGenEleVsEt_->SetDirectory(0);
  selMatchedGenEleVsEta_->SetDirectory(0);
  isoMatchedGenEleVsEt_->SetDirectory(0);
  isoMatchedGenEleVsEta_->SetDirectory(0);
  deltaRReco_->SetDirectory(0);
  deltaRSel_->SetDirectory(0);
  deltaRIso_->SetDirectory(0);
  topMET_->SetDirectory(0);
  genMET_->SetDirectory(0);
  tqafGenMET_->SetDirectory(0);
  corrCaloMET_->SetDirectory(0);
  jetsInvMass_->SetDirectory(0);
  jetsEt_->SetDirectory(0);
  jetsEta_->SetDirectory(0);
  jetsEtVsEta_->SetDirectory(0);
}

void TQAFLayer1Validation::saveHistos() {

  outputRootFile_=new TFile(outputFileName_.c_str(),"RECREATE");

  theNumbers_->SetDirectory(outputRootFile_);
  allGenEleVsEt_->SetDirectory(outputRootFile_);
  allGenEleVsEta_->SetDirectory(outputRootFile_);
  recoMatchedGenEleVsEt_->SetDirectory(outputRootFile_);
  recoMatchedGenEleVsEta_->SetDirectory(outputRootFile_);
  selMatchedGenEleVsEt_->SetDirectory(outputRootFile_);
  selMatchedGenEleVsEta_->SetDirectory(outputRootFile_);
  isoMatchedGenEleVsEt_->SetDirectory(outputRootFile_);
  isoMatchedGenEleVsEta_->SetDirectory(outputRootFile_);
  deltaRReco_->SetDirectory(outputRootFile_);
  deltaRSel_->SetDirectory(outputRootFile_);
  deltaRIso_->SetDirectory(outputRootFile_);
  topMET_->SetDirectory(outputRootFile_);
  genMET_->SetDirectory(outputRootFile_);
  tqafGenMET_->SetDirectory(outputRootFile_);
  corrCaloMET_->SetDirectory(outputRootFile_);
  jetsInvMass_->SetDirectory(outputRootFile_);
  jetsEt_->SetDirectory(outputRootFile_);
  jetsEta_->SetDirectory(outputRootFile_);
  jetsEtVsEta_->SetDirectory(outputRootFile_);

  outputRootFile_->Write();
  outputRootFile_->Close();
  delete outputRootFile_;
}


//define this as a plug-in
DEFINE_FWK_MODULE(TQAFLayer1Validation);
