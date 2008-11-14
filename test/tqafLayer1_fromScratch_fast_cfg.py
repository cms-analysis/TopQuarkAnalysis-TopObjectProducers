import FWCore.ParameterSet.Config as cms

#-------------------------------------------------
# test cfg file for tqaflayer1 production from
# fastsim
#-------------------------------------------------
process = cms.Process("TEST")

## add message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

#-------------------------------------------------
# process configuration; be aware that here there
# is no source; this is called within the fastsim
# bootstrap
#-------------------------------------------------

## invoke fastsim bootstrap
process.load("PhysicsTools.PatAlgos.famos.boostrapWithFamos_cff")

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

## load magnetic field
process.load("Configuration.StandardSequences.MagneticField_cff")

#-------------------------------------------------
# tqaf configuration
#-------------------------------------------------

# FAMOS source
process.load("FastSimulation.Configuration.ttbar_cfi")

## std sequence for tqaf layer1
process.load("TopQuarkAnalysis.TopObjectProducers.tqafLayer1_fast_cff")

#-------------------------------------------------
# process paths;
#-------------------------------------------------

## process path
process.p = cms.Path(process.famosWithEverythingPAT *
                     process.tqafLayer1_withoutTrigMatch)

#-------------------------------------------------
# tqaf event content; first ALL objects are
# dropped in this process; then tqafLayer1
# concent is added
#-------------------------------------------------

## define tqaf layer1 event content
from TopQuarkAnalysis.TopObjectProducers.tqafLayer1_EventContent_cff import *
tqafLayer1EventContent(process)

#
# more or changed jet collections
#
from PhysicsTools.PatAlgos.tools.jetTools import *

switchJetCollection(process, 
                    'sisCone5CaloJets',      # jet collection; must be already in the event when patLayer0 sequence is executed
                    layers       = [0,1],    # if you're not running patLayer1, set 'layers=[0]' 
                    runCleaner   = "CaloJet",# =None if not to clean
                    doJTA        = True,     # run jet-track association & JetCharge
                    doBTagging   = True,     # run b-tagging
                    jetCorrLabel = 'Scone5', # example jet correction name; set to None for no JEC
                    doType1MET   = True      # recompute Type1 MET using these jets
                    )

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
    verbose = cms.untracked.bool(False),
    fileName = cms.untracked.string('TQAFLayer1_Output.fromScratch_fast.root')
)

#-------------------------------------------------
# output paths; in order not to write the
# persistent output to file comment the the output
# path
#-------------------------------------------------

## output
process.outpath = cms.EndPath(process.out)
