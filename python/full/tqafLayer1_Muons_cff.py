import FWCore.ParameterSet.Config as cms

allLayer1Muons.muonSource = 'allLayer0Muons'
allLayer1Muons.addGenMatch = True
allLayer1Muons.genParticleMatch = 'muonMatch'
allLayer1Muons.addTrigMatch = True
allLayer1Muons.trigPrimMatch = ['muonTrigMatchHLT1MuonNonIso', 'muonTrigMatchHLT1MET65']
allLayer1Muons.addResolutions = True
allLayer1Muons.useNNResolutions = False
allLayer1Muons.muonResoFile = 'PhysicsTools/PatUtils/data/Resolutions_muon.root'
allLayer1Muons.isolation.tracker.src = 'layer0MuonIsolations:muIsoDepositTk'
allLayer1Muons.isolation.tracker.deltaR = 0.3
allLayer1Muons.isolation.ecal.src = 'layer0MuonIsolations:muIsoDepositCalByAssociatorTowersecal'
allLayer1Muons.isolation.hcal.deltaR = 0.3
allLayer1Muons.isolation.hcal.src = 'layer0MuonIsolations:muIsoDepositCalByAssociatorTowershcal'
allLayer1Muons.isolation.ecal.deltaR = 0.3
allLayer1Muons.isoDeposits.tracker = 'layer0MuonIsolations:muIsoDepositTk'
allLayer1Muons.isoDeposits.ecal = 'layer0MuonIsolations:muIsoDepositCalByAssociatorTowersecal'
allLayer1Muons.isoDeposits.hcal = 'layer0MuonIsolations:muIsoDepositCalByAssociatorTowershcal'
allLayer1Muons.isoDeposits.user = ['layer0MuonIsolations:muIsoDepositCalByAssociatorTowersho', 'layer0MuonIsolations:muIsoDepositJets']
allLayer1Muons.addMuonID = True
selectedLayer1Muons.src = 'allLayer1Muons'
selectedLayer1Muons.cut = 'pt > 10. & abs(eta) < 2.4 & trackIso < 3 & caloIso < 5'

