import FWCore.ParameterSet.Config as cms

#
# L0 input
#

## import module
from PhysicsTools.PatAlgos.cleaningLayer0.electronCleaner_cfi import allLayer0Electrons

## configure for tqaf
allLayer0Electrons.electronSource        = 'pixelMatchGsfElectrons'
allLayer0Electrons.removeDuplicates      = True
allLayer0Electrons.selection             = cms.PSet( type = cms.string('none') )
allLayer0Electrons.isolation.tracker.src = 'patAODElectronIsolations:eleIsoDepositTk'
allLayer0Electrons.isolation.tracker.cut = 3.0
allLayer0Electrons.isolation.ecal.src    = 'patAODElectronIsolations:eleIsoDepositEcalFromClusts'
allLayer0Electrons.isolation.ecal.cut    = 9999.0
allLayer0Electrons.isolation.hcal.src    = 'patAODElectronIsolations:eleIsoDepositHcalFromTowers'
allLayer0Electrons.isolation.hcal.cut    = 9999.0

#
# genMatch
#

## import module
from PhysicsTools.PatAlgos.mcMatchLayer0.electronMatch_cfi import electronMatch

## configure for tqaf
electronMatch.src                        = 'allLayer0Electrons'
electronMatch.matched                    = 'genParticles'
electronMatch.maxDeltaR                  = 0.5
electronMatch.maxDPtRel                  = 0.5
electronMatch.resolveAmbiguities         = True
electronMatch.resolveByMatchQuality      = False
electronMatch.checkCharge                = True
electronMatch.mcPdgId                    = [11]
electronMatch.mcStatus                   =  [1]



