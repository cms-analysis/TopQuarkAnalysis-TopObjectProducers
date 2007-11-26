#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_SEAL_MODULE();


#include "TopQuarkAnalysis/TopObjectProducers/interface/TopElectronProducer.h"
#include "TopQuarkAnalysis/TopObjectProducers/interface/TopMuonProducer.h"
#include "TopQuarkAnalysis/TopObjectProducers/interface/TopTauProducer.h"
#include "TopQuarkAnalysis/TopObjectProducers/interface/TopJetProducer.h"
#include "TopQuarkAnalysis/TopObjectProducers/interface/TopMETProducer.h"

DEFINE_ANOTHER_FWK_MODULE(TopElectronProducer);
DEFINE_ANOTHER_FWK_MODULE(TopMuonProducer);
DEFINE_ANOTHER_FWK_MODULE(TopTauProducer);
DEFINE_ANOTHER_FWK_MODULE(TopJetProducer);
DEFINE_ANOTHER_FWK_MODULE(TopMETProducer);


#include "TopQuarkAnalysis/TopObjectProducers/interface/TopObjectSelector.h"

DEFINE_ANOTHER_FWK_MODULE(TopElectronSelector);
DEFINE_ANOTHER_FWK_MODULE(TopMuonSelector);
DEFINE_ANOTHER_FWK_MODULE(TopTauSelector);
DEFINE_ANOTHER_FWK_MODULE(TopJetSelector);
DEFINE_ANOTHER_FWK_MODULE(TopMETSelector);
DEFINE_ANOTHER_FWK_MODULE(TopParticleSelector);


#include "TopQuarkAnalysis/TopObjectProducers/interface/TopLeptonCountFilter.h"

DEFINE_ANOTHER_FWK_MODULE(TopLeptonCountFilter);

#include "TopQuarkAnalysis/TopObjectProducers/interface/TopObjectFilter.h"

DEFINE_ANOTHER_FWK_MODULE(TopElectronMinFilter);
DEFINE_ANOTHER_FWK_MODULE(TopMuonMinFilter);
DEFINE_ANOTHER_FWK_MODULE(TopTauMinFilter);
DEFINE_ANOTHER_FWK_MODULE(TopJetMinFilter);
DEFINE_ANOTHER_FWK_MODULE(TopMETMinFilter);
DEFINE_ANOTHER_FWK_MODULE(TopParticleMinFilter);

DEFINE_ANOTHER_FWK_MODULE(TopElectronMaxFilter);
DEFINE_ANOTHER_FWK_MODULE(TopMuonMaxFilter);
DEFINE_ANOTHER_FWK_MODULE(TopTauMaxFilter);
DEFINE_ANOTHER_FWK_MODULE(TopJetMaxFilter);
DEFINE_ANOTHER_FWK_MODULE(TopMETMaxFilter);
DEFINE_ANOTHER_FWK_MODULE(TopParticleMaxFilter);


#include "TopQuarkAnalysis/TopObjectProducers/interface/TopObjectEnergyScale.h"

typedef TopObjectEnergyScale<TopElectron> TopElectronEnergyScale;
typedef TopObjectEnergyScale<TopMuon>     TopMuonEnergyScale;
typedef TopObjectEnergyScale<TopJet>      TopJetEnergyScale;
typedef TopObjectEnergyScale<TopMET>      TopMETEnergyScale;

DEFINE_ANOTHER_FWK_MODULE(TopElectronEnergyScale);
DEFINE_ANOTHER_FWK_MODULE(TopMuonEnergyScale);
DEFINE_ANOTHER_FWK_MODULE(TopJetEnergyScale);
DEFINE_ANOTHER_FWK_MODULE(TopMETEnergyScale);
