import FWCore.ParameterSet.Config as cms

#-------------------------------------------------
# test cfg file for tqaflayer1 production from
# fullsim
#-------------------------------------------------
process = cms.Process("TEST")

## add message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

#-------------------------------------------------
# process configuration
#-------------------------------------------------

## define input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
   #PAT test sample
   #'file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_1_X_2008-07-08_STARTUP_V4-AODSIM.100.root'
   #210 RelVal sample
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_2_1_0/RelValTTbar/GEN-SIM-RECO/STARTUP_V4_v3/0001/061DC5C9-8962-DD11-AB87-001617C3B5F4.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_2_1_0/RelValTTbar/GEN-SIM-RECO/STARTUP_V4_v3/0001/1846FB92-8B62-DD11-BF46-001617C3B5D8.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_2_1_0/RelValTTbar/GEN-SIM-RECO/STARTUP_V4_v3/0001/28BA9967-8A62-DD11-8CBC-001617C3B6CC.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_2_1_0/RelValTTbar/GEN-SIM-RECO/STARTUP_V4_v3/0001/3CE74890-8A62-DD11-A309-001617C3B79A.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_2_1_0/RelValTTbar/GEN-SIM-RECO/STARTUP_V4_v3/0001/6CE93E47-E262-DD11-99D9-000423D6BA18.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_2_1_0/RelValTTbar/GEN-SIM-RECO/STARTUP_V4_v3/0001/8404EE20-8B62-DD11-A6AD-001617C3B6C6.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_2_1_0/RelValTTbar/GEN-SIM-RECO/STARTUP_V4_v3/0001/A2111BED-8E62-DD11-9AB8-000423D98804.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_2_1_0/RelValTTbar/GEN-SIM-RECO/STARTUP_V4_v3/0001/DEFB6B46-8C62-DD11-8643-001617C3B79A.root'
    )
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

from TopQuarkAnalysis.TopObjectProducers.tqafLayer1_genParticles_cff import *   ## pruned genParticles which contain only information
tqafLayer1GenParticles(process)                                                 ## relevant for the TopGenEvnet and stable particles

from TopQuarkAnalysis.TopObjectProducers.tqafLayer1_jetCollections_cff import * ## jet collections of interest for the Top PAG and top
tqafLayer1JetCollections(process)                                               ## analyses


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
