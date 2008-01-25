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
 9) calculate reco and isolation efficiency from TopMuons
 NB obviously the results of some of these checks are sample-dependent
*/
//
// Original Author:  James LAMB
//         Created:  Thu Oct 25 17:44:42 CEST 2007
// $Id: TQAFLayer1Validation.cc,v 1.7 2008/01/22 15:36:38 jlamb Exp $
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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
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

const float MUON_TKISOCUT=3;
const float MUON_CALOISOCUT=5;
const float MUON_ETA_MAX=2.4;
const float MINDR_MUONGENMATCH=0.2;

class TQAFLayer1Validation : public edm::EDAnalyzer {
public:
  explicit TQAFLayer1Validation(const edm::ParameterSet&);
  ~TQAFLayer1Validation();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool findDuplicates(edm::Handle<std::vector<TopElectron> > topElectrons);
  std::vector<reco::GenParticle> getGenPrimEles(edm::Handle<reco::GenParticleCollection> genParHandle);
  std::vector<reco::GenParticle> getGenPrimMuons(edm::Handle<reco::GenParticleCollection> genParHandle);
  std::vector<reco::GenParticle> getGenWDecayProds(edm::Handle<reco::GenParticleCollection> genParHandle);
  std::vector<uint32_t> findNonDupeIndices(edm::Handle<std::vector<reco::PixelMatchGsfElectron> > PMGsfElectronsH);
  std::pair<float,size_t> calcMinDR(reco::Particle par, std::vector<TopElectron> eles);  
  std::pair<float,size_t> calcMinDR(reco::Particle par, std::vector<TopMuon> muons);  
  std::pair<float,size_t> calcMinDR(reco::Particle par, std::vector<TopJet> topJets);
  uint32_t classifyEvent(const reco::GenParticleCollection &genPars);
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

  TH1F *allGenMuonVsEt_;
  TH1F *allGenMuonVsEta_;
  TH1F *recoMatchedGenMuonVsEt_;
  TH1F *recoMatchedGenMuonVsEta_;
  TH1F *selMatchedGenMuonVsEt_;
  TH1F *selMatchedGenMuonVsEta_;
  TH1F *isoMatchedGenMuonVsEt_;
  TH1F *isoMatchedGenMuonVsEta_;
  TH1F *deltaRRecoMuon_;
  TH1F *deltaRIsoMuon_;

  
  TH1F *corrCaloMET_;
  TH1F *topMET_;
  TH1F *genMET_;
  TH1F *tqafGenMET_;

  TH1F *topJetsEt_;
  TH1F *topJetsEta_;
  TH2F *topJetsEtVsEta_;
  TH1F *caloJetsEt_;
  TH1F *caloJetsEta_;
  TH2F *caloJetsEtVsEta_;
  TH1F *corrCaloJetsEt_;
  TH1F *corrCaloJetsEta_;
  TH2F *corrCaloJetsEtVsEta_;
  TH1F *genJetsEt_;
  TH1F *genJetsEta_;
  TH2F *genJetsEtVsEta_;
  TH1F *nTopJetsPerEvt_;
  TH1F *nCaloJetsPerEvt_;
  TH1F *nCaloJetsNoAcceptCutsPerEvt_;
  TH1F *nCorrCaloJetsPerEvt_;
  TH1F *nGenJetsPerEvt_;
  
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
  Handle<std::vector<TopMuon> > topMuonsH;
  iEvent.getByLabel("allLayer1TopMuons","",topMuonsH);
  Handle<reco::GenParticleCollection> genParHandle;
  iEvent.getByLabel("genParticles",genParHandle);
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
//   Handle<std::vector<reco::CaloJet> > corrCaloJetsHandle;
//   iEvent.getByLabel("MCJetCorJetIcone5",corrCaloJetsHandle);
  Handle<std::vector<reco::CaloJet> > caloJetsHandle;
  iEvent.getByLabel("iterativeCone5CaloJets",caloJetsHandle);
  Handle<std::vector<reco::GenJet> > genJetsHandle;
  iEvent.getByLabel("iterativeCone5GenJetsPt10",genJetsHandle);
  Handle<std::vector<TopJet> > topJetsHandle;
  iEvent.getByLabel("allLayer1TopJets",topJetsHandle);
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
  std::vector<reco::GenParticle> genPrimEles=getGenPrimEles(genParHandle);
  //std::vector<reco::GenParticle> genPrimEles=getGenPrimEles(genParHandle);
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
   //topJets Et, eta, et vs. eta
   int nTopJets=0;
   int nCaloJets=0;
   int nCorrCaloJets=0;
   int nGenJets=0;
   for (uint32_t ij=0;ij<topJetsHandle->size();ij++) {
     TopJet topJet=topJetsHandle->at(ij);
     topJetsEt_->Fill(topJet.et());
     topJetsEta_->Fill(topJet.eta());
     topJetsEtVsEta_->Fill(topJet.eta(),topJet.et());
     if (topJet.et()>5 && fabs(topJet.eta()<5)) nTopJets++;
   }
   //uncorrected reco::CaloJets Et, eta, et vs. eta
   for (uint32_t ij=0;ij<caloJetsHandle->size();ij++) {
     reco::CaloJet caloJet=caloJetsHandle->at(ij);
     caloJetsEt_->Fill(caloJet.et());
     caloJetsEta_->Fill(caloJet.eta());
     caloJetsEtVsEta_->Fill(caloJet.eta(),caloJet.et());
     if (caloJet.et()>5 && fabs(caloJet.eta()<5)) nCaloJets++;
   }
   //corrected reco::CaloJets Et, eta, et vs. eta
//    for (uint32_t ij=0;ij<corrCaloJetsHandle->size();ij++) {
//      reco::CaloJet corrCaloJet=corrCaloJetsHandle->at(ij);
//      corrCaloJetsEt_->Fill(corrCaloJet.et());
//      corrCaloJetsEta_->Fill(corrCaloJet.eta());
//      corrCaloJetsEtVsEta_->Fill(corrCaloJet.eta(),corrCaloJet.et());
//      if (corrCaloJet.et()>5 && fabs(corrCaloJet.eta()<5)) nCorrCaloJets++;
//    }
   for (uint32_t ij=0;ij<genJetsHandle->size();ij++) {
     reco::GenJet genJet=genJetsHandle->at(ij);
     genJetsEt_->Fill(genJet.et());
     genJetsEta_->Fill(genJet.eta());
     genJetsEtVsEta_->Fill(genJet.eta(),genJet.et());
     // switch to 'genParticle here, as soon as this also happens in the 'GenJet' class:
//      std::vector<const reco::GenParticle *> genJetConst=genJet.getConstituents();
     std::vector<const reco::GenParticleCandidate *> genJetConst=genJet.getConstituents();
     if (genJet.et()>5 && fabs(genJet.eta()<5)) nGenJets++;
   }
   nTopJetsPerEvt_->Fill(nTopJets);
   nCaloJetsPerEvt_->Fill(nCaloJets);
   nCorrCaloJetsPerEvt_->Fill(nCorrCaloJets);
   nGenJetsPerEvt_->Fill(nGenJets);
   
   //task 9: efficiency for muons
   std::vector<TopMuon> isoMuons;
   for (uint32_t im=0;im<topMuonsH->size();im++) {
     //for consistency with fake rate (defined elsewhere), require to be within eta<2.4
     if (fabs(topMuonsH->at(im).eta()) > MUON_ETA_MAX) continue;
     double tkIso=topMuonsH->at(im).getTrackIso();
     double caloIso=topMuonsH->at(im).getCaloIso();
     if (tkIso<MUON_TKISOCUT && caloIso < MUON_CALOISOCUT) isoMuons.push_back(topMuonsH->at(im));
   }
   std::vector<reco::GenParticle> genPrimMuons=getGenPrimMuons(genParHandle);
   for (uint32_t ip=0;ip<genPrimMuons.size();ip++) {
     
     //require the gen muon to be within eta acceptance
     if (fabs(genPrimMuons[ip].eta()) > MUON_ETA_MAX) continue;
     
     allGenMuonVsEt_->Fill(genPrimMuons[ip].et());
     allGenMuonVsEta_->Fill(genPrimMuons[ip].eta());
     
     //find the nearest reco muon
     std::pair<float, size_t> match=calcMinDR(genPrimMuons[ip], *topMuonsH);  
     deltaRRecoMuon_->Fill(match.first);
     if (match.first<MINDR_MUONGENMATCH) {
       recoMatchedGenMuonVsEt_->Fill(genPrimMuons[ip].et());
       recoMatchedGenMuonVsEta_->Fill(genPrimMuons[ip].eta());
     }
     //find nearest isolated muon
     std::pair<float, size_t> matchIso=calcMinDR(genPrimMuons[ip], isoMuons);
     deltaRIsoMuon_->Fill(matchIso.first);
     if (matchIso.first<MINDR_MUONGENMATCH) {
       isoMatchedGenMuonVsEt_->Fill(genPrimMuons[ip].et());
       isoMatchedGenMuonVsEta_->Fill(genPrimMuons[ip].eta());
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
  theNumbers_->SetBinContent(9,nEvt_);

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
//that means the electrons that come from the W or Z
//in practice that means the electrons that come from the electrons that come from the W of Z
// (so after any FSR)
std::vector<reco::GenParticle> TQAFLayer1Validation::getGenPrimEles(edm::Handle<reco::GenParticleCollection> genParHandle) {
  
  std::vector<reco::GenParticle> output;
  for (uint32_t ip=0;ip<genParHandle->size();ip++) {
    if (abs((*genParHandle)[ip].pdgId())!=11) continue;
    if ((*genParHandle)[ip].status()!=1) continue;
    if ((*genParHandle)[ip].numberOfMothers()<1) continue;
    if (abs((*genParHandle)[ip].mother(0)->pdgId())!=11) continue;
    if ((*genParHandle)[ip].mother(0)->numberOfMothers()<1) continue;
    if (abs((*genParHandle)[ip].mother(0)->mother(0)->pdgId())!=23 &&
	abs((*genParHandle)[ip].mother(0)->mother(0)->pdgId())!=24) continue;
    //yikes, now we finally have electrons whose mothers are electrons whose mothers are Ws or Zs
    //reco::GenParticle tmp=static_cast<reco::GenParticle> ((*genParHandle)[ip]);
    //output.push_back((*genParHandle)[ip]);
    //reco::GenParticle tmp=(reco::GenParticle) ((*genParHandle)[ip]);
    reco::GenParticle tmp = *(dynamic_cast<reco::GenParticle *>(const_cast<reco::GenParticle *>(&(*genParHandle)[ip])));
    output.push_back(tmp);
  }
  return output;

}

//get a vector of the primary muons out of the event
//that means the muons that come from the W or Z
//in practice that means the muons that come from the muons that come from the W or the Z
// (so after any FSR)
std::vector<reco::GenParticle> TQAFLayer1Validation::getGenPrimMuons(edm::Handle<reco::GenParticleCollection> genParHandle) {
  
  std::vector<reco::GenParticle> output;
  for (uint32_t ip=0;ip<genParHandle->size();ip++) {
    if (abs((*genParHandle)[ip].pdgId())!=13) continue;
    if ((*genParHandle)[ip].status()!=1) continue;
    if ((*genParHandle)[ip].numberOfMothers()<1) continue;
    if (abs((*genParHandle)[ip].mother(0)->pdgId())!=13) continue;
    if ((*genParHandle)[ip].mother(0)->numberOfMothers()<1) continue;
    if (abs((*genParHandle)[ip].mother(0)->mother(0)->pdgId())!=23 &&
	abs((*genParHandle)[ip].mother(0)->mother(0)->pdgId())!=24
	) continue;
    //yikes, now we finally have muons whose mothers are muons whose mothers are Ws
    reco::GenParticle tmp = *(dynamic_cast<reco::GenParticle *>(const_cast<reco::GenParticle *>(&(*genParHandle)[ip])));
    output.push_back(tmp);
  }
  return output;

}


//get the vector of gen-level particles who's mother's mother's are W's.  So this should be the W decay products after FSR.
//the vector is guaranteed to be ordered such that the decay products of the W+ come before that of the W-
std::vector<reco::GenParticle> TQAFLayer1Validation::getGenWDecayProds(edm::Handle<reco::GenParticleCollection> genParHandle) {
  
  std::vector<reco::GenParticle> output;
  //first the W+
  for (uint32_t ip=0;ip<genParHandle->size();ip++) {
    if ((*genParHandle)[ip].status()!=1) continue;
    if ((*genParHandle)[ip].numberOfMothers()<1) continue;
    if ((*genParHandle)[ip].mother(0)->numberOfMothers()<1) continue;
    if ((*genParHandle)[ip].mother(0)->mother(0)->pdgId()!=24) continue;

    reco::GenParticle tmp = *(dynamic_cast<reco::GenParticle *>(const_cast<reco::GenParticle *>(&(*genParHandle)[ip])));
    output.push_back(tmp);
    //break the first loop if got the two W+ decay products
    if (output.size()==2) break;
  }
  //now the W-
  for (uint32_t ip=0;ip<genParHandle->size();ip++) {
    if ((*genParHandle)[ip].status()!=1) continue;
    if ((*genParHandle)[ip].numberOfMothers()<1) continue;
    if ((*genParHandle)[ip].mother(0)->numberOfMothers()<1) continue;
    if ((*genParHandle)[ip].mother(0)->mother(0)->pdgId()!=-24) continue;

    reco::GenParticle tmp = *(dynamic_cast<reco::GenParticle *>(const_cast<reco::GenParticle *>(&(*genParHandle)[ip])));
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
  
std::pair<float,size_t> TQAFLayer1Validation::calcMinDR(reco::Particle par, std::vector<TopMuon> muons) {

  float mindr=9999.;
  size_t indexMin=0;
  
  for (size_t im=0;im<muons.size();im++) {
    
    float dr=deltaR(par,muons[im]);
    if (dr<mindr) {
      mindr=dr;
      indexMin=im;
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
uint32_t TQAFLayer1Validation::classifyEvent(const reco::GenParticleCollection &genPars) {
  
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
  for (uint32_t ip=0;ip<100 && ip<genPars.size();ip++) {
    
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

  allGenMuonVsEt_=new TH1F("allGenMuonVsEt_","Et Spectrum of Gen-level muons",500,0,500);
  allGenMuonVsEta_=new TH1F("allGenMuonVsEta_","Eta Spectrum of Gen-level muons",100,-5,5);
  recoMatchedGenMuonVsEt_=new TH1F("recoMatchedGenMuonVsEt_","Et Spectrum of Gen-level muons found in Reco",500,0,500);
  recoMatchedGenMuonVsEta_=new TH1F("recoMatchedGenMuonVsEta_","Eta Spectrum of Gen-level muons found in Reco",100,-5,5);
  isoMatchedGenMuonVsEt_=new TH1F("isoMatchedGenMuonVsEt_","Et Spectrum of Gen-level muons found in Reco and isolated",500,0,500);
  isoMatchedGenMuonVsEta_=new TH1F("isoMatchedGenMuonVsEta_","Eta Spectrum of Gen-level muons found in Reco and isolated",100,-5,5);
  deltaRRecoMuon_=new TH1F("deltaRRecoMuon_","deltaR between gen-level electron and nearest reco muon",1000,0,1);
  deltaRIsoMuon_=new TH1F("deltaRIsoMuon_","deltaR between gen-level muon and nearest \"loose\" selected and isolated muon",1000,0,1);

  topMET_=new TH1F("topMET_","TopMET",200,0,200);
  genMET_=new TH1F("genMET_","GenMET",200,0,200);
  tqafGenMET_=new TH1F("tqafGenMET_","Gen MET from TQAF",200,0,200);
  corrCaloMET_=new TH1F("corrCaloMET_","Corrected Calo MET",200,0,200);
  topJetsEt_=new TH1F("topJetsEt_","TopJets Et",200,0,200);
  topJetsEta_=new TH1F("topJetsEta_","TopJets Eta",100,-5,5);
  topJetsEtVsEta_=new TH2F("topJetsEtVsEta_","TopJets Et vs. Eta",100,-5,5,200,0,200);  
  corrCaloJetsEt_=new TH1F("corrCaloJetsEt_","Corrected CaloJets Et",200,0,200);
  corrCaloJetsEta_=new TH1F("corrCaloJetsEta_","Corrected CaloJets Eta",100,-5,5);
  corrCaloJetsEtVsEta_=new TH2F("corrCaloJetsEtVsEta_","Corrected CaloJets Et vs. Eta",100,-5,5,200,0,200);  
  caloJetsEt_=new TH1F("caloJetsEt_","Uncorrected CaloJets Et",200,0,200);
  caloJetsEta_=new TH1F("caloJetsEta_","Uncorrected CaloJets Eta",100,-5,5);
  caloJetsEtVsEta_=new TH2F("caloJetsEtVsEta_","Uncorrected CaloJets Et vs. Eta",100,-5,5,200,0,200);  
  genJetsEt_=new TH1F("genJetsEt_","GenJets Et",200,0,200);
  genJetsEta_=new TH1F("genJetsEta_","GenJets Eta",100,-5,5);
  genJetsEtVsEta_=new TH2F("genJetsEtVsEta_","GenJets Et vs. Eta",100,-5,5,200,0,200);  
  nTopJetsPerEvt_=new TH1F("nTopJetsPerEvt_","Number of Top Jets per Events, et>5 abs(eta)<5",100,0,100);
  nCaloJetsPerEvt_=new TH1F("nCaloJetsPerEvt_","Number of Uncorrected Calo Jets per Events, et>5 abs(eta)<5",100,0,100);
  nCorrCaloJetsPerEvt_=new TH1F("nCorrCaloJetsPerEvt_","Number of Corrected Calo Jets per Events, et>5 abs(eta)<5",100,0,100);
  nGenJetsPerEvt_=new TH1F("nGenJetsPerEvt_","Number of Gen Jets per Events, et>5 abs(eta)<5",100,0,100);

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

  allGenMuonVsEt_->SetDirectory(0);
  allGenMuonVsEta_->SetDirectory(0);
  recoMatchedGenMuonVsEt_->SetDirectory(0);
  recoMatchedGenMuonVsEta_->SetDirectory(0);
  isoMatchedGenMuonVsEt_->SetDirectory(0);
  isoMatchedGenMuonVsEta_->SetDirectory(0);
  deltaRRecoMuon_->SetDirectory(0);
  deltaRIsoMuon_->SetDirectory(0);

  topMET_->SetDirectory(0);
  genMET_->SetDirectory(0);
  tqafGenMET_->SetDirectory(0);
  corrCaloMET_->SetDirectory(0);
  topJetsEt_->SetDirectory(0);
  topJetsEta_->SetDirectory(0);
  topJetsEtVsEta_->SetDirectory(0);
  corrCaloJetsEt_->SetDirectory(0);
  corrCaloJetsEta_->SetDirectory(0);
  corrCaloJetsEtVsEta_->SetDirectory(0);
  caloJetsEt_->SetDirectory(0);
  caloJetsEta_->SetDirectory(0);
  caloJetsEtVsEta_->SetDirectory(0);
  genJetsEt_->SetDirectory(0);
  genJetsEta_->SetDirectory(0);
  genJetsEtVsEta_->SetDirectory(0);
  nTopJetsPerEvt_->SetDirectory(0);
  nCaloJetsPerEvt_->SetDirectory(0);
  nCorrCaloJetsPerEvt_->SetDirectory(0);
  nGenJetsPerEvt_->SetDirectory(0);

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

  allGenMuonVsEt_->SetDirectory(outputRootFile_);
  allGenMuonVsEta_->SetDirectory(outputRootFile_);
  recoMatchedGenMuonVsEt_->SetDirectory(outputRootFile_);
  recoMatchedGenMuonVsEta_->SetDirectory(outputRootFile_);
  isoMatchedGenMuonVsEt_->SetDirectory(outputRootFile_);
  isoMatchedGenMuonVsEta_->SetDirectory(outputRootFile_);
  deltaRRecoMuon_->SetDirectory(outputRootFile_);
  deltaRIsoMuon_->SetDirectory(outputRootFile_);


  topMET_->SetDirectory(outputRootFile_);
  genMET_->SetDirectory(outputRootFile_);
  tqafGenMET_->SetDirectory(outputRootFile_);
  corrCaloMET_->SetDirectory(outputRootFile_);
  topJetsEt_->SetDirectory(outputRootFile_);
  topJetsEta_->SetDirectory(outputRootFile_);
  topJetsEtVsEta_->SetDirectory(outputRootFile_);
  caloJetsEt_->SetDirectory(outputRootFile_);
  caloJetsEta_->SetDirectory(outputRootFile_);
  caloJetsEtVsEta_->SetDirectory(outputRootFile_);
  corrCaloJetsEt_->SetDirectory(outputRootFile_);
  corrCaloJetsEta_->SetDirectory(outputRootFile_);
  corrCaloJetsEtVsEta_->SetDirectory(outputRootFile_);
  genJetsEt_->SetDirectory(outputRootFile_);
  genJetsEta_->SetDirectory(outputRootFile_);
  genJetsEtVsEta_->SetDirectory(outputRootFile_);
  nTopJetsPerEvt_->SetDirectory(outputRootFile_);
  nCaloJetsPerEvt_->SetDirectory(outputRootFile_);
  nCorrCaloJetsPerEvt_->SetDirectory(outputRootFile_);
  nGenJetsPerEvt_->SetDirectory(outputRootFile_);



  outputRootFile_->Write();
  outputRootFile_->Close();
  delete outputRootFile_;
}


//define this as a plug-in
DEFINE_FWK_MODULE(TQAFLayer1Validation);
