import FWCore.ParameterSet.Config as cms

#-------------------------------------------------
# test cfg file for tqaflayer1 production from
# fullsim
#-------------------------------------------------
process = cms.Process("USER")

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
   'file:/afs/cern.ch/user/r/rwolf/pccmsuhh06/testPatTuple_prov.root'
    )
)

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)

## configure geometry
process.load("Configuration.StandardSequences.Geometry_cff")

## configure conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V7::All')

## load magnetic field
process.load("Configuration.StandardSequences.MagneticField_cff")


#-------------------------------------------------
# do the jet re-reconstruction from caloTowers
# on the patTuple
#-------------------------------------------------

## add geometry needed for pileup corrections for jet reco
process.load("TopQuarkAnalysis.TopObjectProducers.patTuple_JetReco_cff")

#-------------------------------------------------
# process paths;
#-------------------------------------------------

## process path
process.p = cms.Path(process.patTupleJetReco)

#-------------------------------------------------
# pat tuple event content; first ALL objects
# are dropped in this process; then patTuple
# content is added
#-------------------------------------------------

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
    verbose = cms.untracked.bool(True),
    dropMetaDataForDroppedData = cms.untracked.bool(False),                           
    fileName = cms.untracked.string('testPatTuple_prov.root')
)

#-------------------------------------------------
# output paths; in order not to write the
# persistent output to file comment the output
# path
#-------------------------------------------------

## output
process.outpath = cms.EndPath(process.out)
