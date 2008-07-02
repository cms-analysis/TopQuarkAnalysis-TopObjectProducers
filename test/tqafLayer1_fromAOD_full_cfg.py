import FWCore.ParameterSet.Config as cms

process = cms.Process("TQAF")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("TopQuarkAnalysis.TopObjectProducers.tqafLayer1_full_cff")

process.load("TopQuarkAnalysis.TopObjectProducers.tqafLayer1_EventContent_cff")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTbar-210p5.1-AODSIM.100.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.tqafEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *')
)
process.EventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )
)
process.out = cms.OutputModule("PoolOutputModule",
    process.EventSelection,
    process.tqafEventContent,
    verbose = cms.untracked.bool(True),
    fileName = cms.untracked.string('TQAFLayer1_Output.fromAOD_full.root')
)

process.p = cms.Path(process.tqafLayer1)
process.outpath = cms.EndPath(process.out)
process.MessageLogger.cerr.threshold = 'INFO'
process.tqafEventContent.outputCommands.extend(process.patLayer1EventContent.outputCommands)

