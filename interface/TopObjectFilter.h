//
// Author:  Steven Lowette
// Created: Tue Oct  2 19:55:00 CEST 2007
//
// $Id: TopObjectFilter.h,v 1.1.2.3 2007/11/01 01:11:32 lowette Exp $
//

#ifndef TopObjectProducer_TopObjectFilter_h
#define TopObjectProducer_TopObjectFilter_h


#include "PhysicsTools/UtilAlgos/interface/AnySelector.h"
#include "PhysicsTools/UtilAlgos/interface/MinNumberSelector.h"
#include "TopQuarkAnalysis/TopObjectProducers/interface/MaxNumberSelector.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectCountFilter.h"

#include "AnalysisDataFormats/TopObjects/interface/TopElectron.h"
#include "AnalysisDataFormats/TopObjects/interface/TopMuon.h"
#include "AnalysisDataFormats/TopObjects/interface/TopTau.h"
#include "AnalysisDataFormats/TopObjects/interface/TopJet.h"
#include "AnalysisDataFormats/TopObjects/interface/TopMET.h"
#include "AnalysisDataFormats/TopObjects/interface/TopParticle.h"

#include <vector>


typedef ObjectCountFilter<std::vector<TopElectron>, AnySelector<std::vector<TopElectron>::value_type>, MinNumberSelector> TopElectronMinFilter;
typedef ObjectCountFilter<std::vector<TopMuon>,     AnySelector<std::vector<TopMuon>::value_type>,     MinNumberSelector> TopMuonMinFilter;
typedef ObjectCountFilter<std::vector<TopTau>,      AnySelector<std::vector<TopTau>::value_type>,      MinNumberSelector> TopTauMinFilter;
typedef ObjectCountFilter<std::vector<TopJet>,      AnySelector<std::vector<TopJet>::value_type>,      MinNumberSelector> TopJetMinFilter;
typedef ObjectCountFilter<std::vector<TopMET>,      AnySelector<std::vector<TopMET>::value_type>,      MinNumberSelector> TopMETMinFilter;
typedef ObjectCountFilter<std::vector<TopParticle>, AnySelector<std::vector<TopParticle>::value_type>, MinNumberSelector> TopParticleMinFilter;

typedef ObjectCountFilter<std::vector<TopElectron>, AnySelector<std::vector<TopElectron>::value_type>, MaxNumberSelector> TopElectronMaxFilter;
typedef ObjectCountFilter<std::vector<TopMuon>,     AnySelector<std::vector<TopMuon>::value_type>    , MaxNumberSelector> TopMuonMaxFilter;
typedef ObjectCountFilter<std::vector<TopTau>,      AnySelector<std::vector<TopTau>::value_type>     , MaxNumberSelector> TopTauMaxFilter;
typedef ObjectCountFilter<std::vector<TopJet>,      AnySelector<std::vector<TopJet>::value_type>     , MaxNumberSelector> TopJetMaxFilter;
typedef ObjectCountFilter<std::vector<TopMET>,      AnySelector<std::vector<TopMET>::value_type>     , MaxNumberSelector> TopMETMaxFilter;
typedef ObjectCountFilter<std::vector<TopParticle>, AnySelector<std::vector<TopParticle>::value_type>, MaxNumberSelector> TopParticleMaxFilter;


#endif
