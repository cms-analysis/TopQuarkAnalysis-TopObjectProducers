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

DEFINE_ANOTHER_FWK_MODULE(CaloJetSelector);
DEFINE_ANOTHER_FWK_MODULE(TopElectronSelector);
DEFINE_ANOTHER_FWK_MODULE(TopTauSelector);
DEFINE_ANOTHER_FWK_MODULE(TopMuonSelector);
DEFINE_ANOTHER_FWK_MODULE(TopJetSelector);
DEFINE_ANOTHER_FWK_MODULE(TopMETSelector);
DEFINE_ANOTHER_FWK_MODULE(TopParticleSelector);

#include "TopQuarkAnalysis/TopObjectProducers/interface/TopObjectFilter.h"

DEFINE_ANOTHER_FWK_MODULE(TopElectronCountFilter);
DEFINE_ANOTHER_FWK_MODULE(TopMuonCountFilter);
DEFINE_ANOTHER_FWK_MODULE(TopTauCountFilter);
DEFINE_ANOTHER_FWK_MODULE(TopJetCountFilter);
DEFINE_ANOTHER_FWK_MODULE(TopMETCountFilter);
DEFINE_ANOTHER_FWK_MODULE(TopParticleCountFilter);
