//
// $Id$
//

#ifndef TopObjectProducers_TopJetProducer_h
#define TopObjectProducers_TopJetProducer_h

/**
  \class    TopJetProducer TopJetProducer.h "TopQuarkAnalysis/TopObjectProducers/interface/TopJetProducer.h"
  \brief    Produces TopJet's

   TopJetProducer produces TopJet's starting from a JetType collection,
   with possible adding of resolutions and more things to come

  \author   Jan Heyninck
  \version  $Id$
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "PhysicsTools/Utilities/interface/EtComparator.h"

#include "AnalysisDataFormats/TopObjects/interface/TopJet.h"
#include "AnalysisDataFormats/TopObjects/interface/TopElectron.h"
#include "AnalysisDataFormats/TopObjects/interface/TopMuon.h"

#include "PhysicsTools/JetCharge/interface/JetCharge.h"
#include "TopQuarkAnalysis/TopObjectProducers/interface/SimpleJetTrackAssociator.h"


class JetFlavourIdentifier;
class TopObjectResolutionCalc;


class TopJetProducer : public edm::EDProducer {

  public:

    explicit TopJetProducer(const edm::ParameterSet & iConfig);
    ~TopJetProducer();

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  private:

    // TEMP Jet cleaning from electrons
    std::vector<TopElectron> selectIsolated(const std::vector<TopElectron> & electrons, float isoCut);
    std::vector<TopMuon>     selectIsolated(const std::vector<TopMuon> & muons,         float isoCut);
    // TEMP End

    // configurables
    edm::InputTag     caliJetsSrc_;
    edm::InputTag     recJetsSrc_;
    // TEMP Jet cleaning from electrons
    bool              doJetCleaning_;
    edm::InputTag     topElectronsLabel_;
    edm::InputTag     topMuonsLabel_;
    float             LEPJETDR_;
    float             ELEISOCUT_;
    float             MUISOCUT_;
    // TEMP End
    bool              getJetMCFlavour_;
    bool              addGenPartonMatch_;
    edm::InputTag     genPartonSrc_;
    bool              addGenJetMatch_;
    edm::InputTag     genJetSrc_;
    bool              addPartonJetMatch_;
    edm::InputTag     partonJetSrc_;
    bool              addResolutions_;
    bool              useNNReso_;
    std::string       caliJetResoFile_;
    bool              addBTagInfo_;
    bool              ignoreTrackCountingFromAOD_;
    bool              ignoreTrackProbabilityFromAOD_;
    bool              ignoreSoftMuonFromAOD_;
    bool              ignoreSoftElectronFromAOD_;
    bool              addDiscriminators_; 
    bool              addJetTagRefs_;
    bool              addAssociatedTracks_;
    edm::ParameterSet trackAssociationPSet_;
    bool              addJetCharge_;
    edm::ParameterSet jetChargePSet_;
    // tools
    JetFlavourIdentifier                   * jetFlavId_;
    TopObjectResolutionCalc                * theResoCalc_;
    reco::helper::SimpleJetTrackAssociator   simpleJetTrackAssociator_;
    JetCharge                              * jetCharge_;
    EtInverseComparator<TopJet>                      eTComparator_;

};


#endif
