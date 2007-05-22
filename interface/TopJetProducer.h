// -*- C++ -*-
//
// Package:    TopJetProducer
// Class:      TopJetProducer
// 
/**\class TopJetProducer TopJetProducer.cc Top/TopEventProducers/src/TopJetProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jan Heyninck
//         Created:  Tue Apr  10 12:01:49 CEST 2007
// $Id: TopJetProducer.h,v 1.2 2007/05/08 14:01:21 heyninck Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "AnalysisDataFormats/TopObjects/interface/TopJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "TopQuarkAnalysis/TopObjectResolutions/interface/TopObjectResolutionCalc.h"
#include "PhysicsTools/Utilities/interface/EtComparator.h"


#include <vector>
#include <Math/VectorUtil.h>

using namespace std;

//
// class decleration
//

class TopJetProducer : public edm::EDProducer {
   public:
      explicit TopJetProducer(const edm::ParameterSet&);
      ~TopJetProducer();

      virtual void produce(edm::Event&, const edm::EventSetup&);
   private:
     string jetTagsLabel_;
     string recJetsLabel_;
     string caliJetsLabel_;
     string caliJetResoFile_;
     double recJetETcut_;
     double jetEtaCut_;
     int minNrConstis_;
     bool addResolutions_;  
     EtInverseComparator<TopJet> eTComparator;
     TopObjectResolutionCalc *jetsResCalc;

};