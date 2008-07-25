import FWCore.ParameterSet.Config as cms

#-------------------------------------------------
# test cfg file for tqaflayer1 production from
# fullsim
#-------------------------------------------------
process = cms.Process("TQAF")

## add message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

#-------------------------------------------------
# process configuration
#-------------------------------------------------

## define input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_1_X_2008-07-08_STARTUP_V4-AODSIM.100.root')
)

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)

## configure geometry
process.load("Configuration.StandardSequences.Geometry_cff")

## configure conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V4::All')

## load magnetic field
process.load("Configuration.StandardSequences.MagneticField_cff")


#-------------------------------------------------
# tqaf configuration
#-------------------------------------------------

## std sequence for tqaf layer1
process.load("TopQuarkAnalysis.TopObjectProducers.tqafLayer1_full_cff")


## add more jet collections to the tqaf layer1
from PhysicsTools.PatAlgos.tools.jetTools import *

## add sisCone5CaloJets
addJetCollection(process,
                 'sisCone5CaloJets',
                 'SisCone5Calo',
                 runCleaner="CaloJet",
                 doJTA=True,doBTagging=True,
                 jetCorrLabel='Scone5',
                 doType1MET=True,
                 doL1Counters=False
                 )

## add kt6CaloJets
addJetCollection(process,
                 'kt6CaloJets',
                 'Kt6Calo',
                 runCleaner=None,
                 doJTA=True,
                 doBTagging=False,
                 jetCorrLabel=None,
                 doType1MET=True,
                 doL1Counters=False
                 )

## add icone5PFJets
addJetCollection(process,
                 'iterativeCone5PFJets',
                 'Icone5PFlow',
                 runCleaner=None,
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=None,
                 doType1MET=True,
                 doL1Counters=False
                 )

#-------------------------------------------------
# process paths;
#-------------------------------------------------

## process path
process.p = cms.Path(process.tqafLayer1)


#-------------------------------------------------
# tqaf event content; first ALL objects are
# dropped in this process; then tqafLayer1
# concent is added
#-------------------------------------------------

## define event content
process.tqafEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *')
)

## define tqaf layer1 event content
process.load("TopQuarkAnalysis.TopObjectProducers.tqafLayer1_EventContent_cff")
process.tqafEventContent.outputCommands.extend(process.patLayer1EventContent.outputCommands)
## keep the additionally produced jet
## collections in the event record
process.tqafEventContent.outputCommands.extend(["keep *_selectedLayer1Jets*_*_*"])

#-------------------------------------------------
# process output; first the event selection is
# defined: only those events that have passed the
# full production path are selected and written
# to file; the event content has been defined
# above
#-------------------------------------------------

## define event selection
process.EventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )
)

## configure output module
process.out = cms.OutputModule("PoolOutputModule",
    process.EventSelection,
    process.tqafEventContent,
    verbose = cms.untracked.bool(True),
    fileName = cms.untracked.string('TQAFLayer1_Output.fromAOD_full.root')
)


#-------------------------------------------------
# output paths; in order not to write the
# persistent output to file comment the output
# path
#-------------------------------------------------

## output
process.outpath = cms.EndPath(process.out)
