import FWCore.ParameterSet.Config as cms

process = cms.Process("TQAF")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("PhysicsTools.PatAlgos.famos.boostrapWithFamos_cff")

process.load("TopQuarkAnalysis.TopObjectProducers.tqafLayer1_fast_cff")

process.load("TopQuarkAnalysis.TopObjectProducers.tqafLayer1_EventContent_cff")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
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
    verbose = cms.untracked.bool(False),
    fileName = cms.untracked.string('TQAFLayer1_Output.fromScratch_fast.root')
)

process.p = cms.Path(process.famosWithEverythingPAT*process.tqafLayer1_withoutTrigMatch)
process.outpath = cms.EndPath(process.out)
process.MessageLogger.cerr.threshold = 'INFO'
process.tqafEventContent.outputCommands.extend(process.patLayer1EventContent.outputCommands)

