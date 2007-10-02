//
// Author:  Jan Heyninck, Steven Lowette
// Created: Tue Apr  10 12:01:49 CEST 2007
//
// $Id$
//

#ifndef TopObjectProducers_TopMETProducer_h
#define TopObjectProducers_TopMETProducer_h

/**
  \class    TopMETProducer TopMETProducer.h "TopQuarkAnalysis/TopObjectProducers/interface/TopMETProducer.h"
  \brief    Produces TopMET's

   TopMETProducer produces TopMET's starting from a METType collection,
   with possible matching to generator level, addition of resolutions and
   correction for a given collection of muons.

  \author   Jan Heyninck, Steven Lowette
  \version  $Id$
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/Utilities/interface/EtComparator.h"

#include "AnalysisDataFormats/TopObjects/interface/TopMET.h"


class TopObjectResolutionCalc;


class TopMETProducer : public edm::EDProducer {

  public:

    explicit TopMETProducer(const edm::ParameterSet & iConfig);
    ~TopMETProducer();

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  private:

    // configurables
    edm::InputTag metSrc_;
    bool          addGenMET_;
    edm::InputTag genPartSrc_;
    bool          addResolutions_;
    bool          useNNReso_;
    std::string   metResoFile_;
    bool          addMuonCorr_;
    edm::InputTag muonSrc_;
    // tools
    TopObjectResolutionCalc * metResoCalc_;
    GreaterByEt<TopMET> eTComparator_;

};


#endif
