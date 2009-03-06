import FWCore.ParameterSet.Config as cms

#-------------------------------------------------
# test cfg file for tqaflayer1 production from
# fastsim
#-------------------------------------------------
process = cms.Process("TQAF")

## add message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

#-------------------------------------------------
# process configuration; be aware that here there
# is no source; this is called within the fastsim
# bootstrap
#-------------------------------------------------

# FAMOS source
process.load("FastSimulation.Configuration.ttbar_cfi")

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

## std sequence for tqaf layer1
process.load("TopQuarkAnalysis.TopObjectProducers.tqafLayer1_cff")

## switch off trigger matching in pat
from PhysicsTools.PatAlgos.tools.trigTools import switchLayer1Off
switchLayer1Off(process)

#-------------------------------------------------
# process paths;
#-------------------------------------------------

## process path
process.p = cms.Path(process.famosWithEverything *
                     process.tqafLayer1_withoutTrigMatch
                     )

#-------------------------------------------------
# tqaf event content; first ALL objects are
# dropped in this process; then tqafLayer1
# concent is added
#-------------------------------------------------

## define tqaf layer1 event content
from TopQuarkAnalysis.TopObjectProducers.tqafLayer1_EventContent_cff import *
makeTqafLayer1EventContent(process)

## uncomment the following two lines and edit the
## corresponding file to replacegenbParticles by
## pruned genParticles

#from TopQuarkAnalysis.TopObjectProducers.tools.pruneGenParticles_cff import *
#tqafPruneGenParticles(process)

#-------------------------------------------------
# add more and/or change jet collection(s)
#-------------------------------------------------
from PhysicsTools.PatAlgos.tools.jetTools import *

## uncomment the following two lines and edit the
## corresponding file to add more jet collections

#from TopQuarkAnalysis.TopObjectProducers.tools.addJetCollections_cff import *
#tqafAddJetCollections(process)

## switch pat default jet collection from IC5 to
## SC5
switchJetCollection(process, 
                    'sisCone5CaloJets',             # jet collection; must be already in the event when patLayer0 sequence is executed
                    layers       = [0,1],           # if you're not running patLayer1, set 'layers=[0]' 
                    runCleaner   = "CaloJet",       # =None if not to clean
                    doJTA        = True,            # run jet-track association & JetCharge
                    doBTagging   = True,            # run b-tagging
                    jetCorrLabel = None,            # example jet correction name; set to None for no JEC
                    doType1MET   = True,            # recompute Type1 MET using these jets
                    genJetCollection = "sisCone5GenJets"                       
                    )

## now set JEC by hand
process.jetCorrFactors.jetSource = cms.InputTag("sisCone5CaloJets")
process.jetCorrFactors.L1Offset  = cms.string('none')
process.jetCorrFactors.L2Relative= cms.string('Summer08_L2Relative_SC5Calo')
process.jetCorrFactors.L3Absolute= cms.string('Summer08_L3Absolute_SC5Calo')
process.jetCorrFactors.L4EMF     = cms.string('none')
process.jetCorrFactors.L5Flavor  = cms.string('none')
process.jetCorrFactors.L6UE      = cms.string('none')
process.jetCorrFactors.L7Parton  = cms.string('none')

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
    dropMetaDataForDroppedData = cms.untracked.bool(True),                                
    fileName = cms.untracked.string('TQAFLayer1_Output.fromAOD_fast.root')
)

#-------------------------------------------------
# output paths; in order not to write the
# persistent output to file comment the output
# path
#-------------------------------------------------

## output
process.outpath = cms.EndPath(process.out)
