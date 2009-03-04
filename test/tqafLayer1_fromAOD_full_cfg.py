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
    fileNames = cms.untracked.vstring(
   #PAT test sample for 2.2.X
   #'file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_2_X_2008-11-03-STARTUP_V7-AODSIM.100.root'
   #PAT test sample for 2.1.X
    'file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_1_X_2008-07-08_STARTUP_V4-AODSIM.100.root'    
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
process.GlobalTag.globaltag = cms.string('STARTUP_V7::All')

## load magnetic field
process.load("Configuration.StandardSequences.MagneticField_cff")


#-------------------------------------------------
# tqaf configuration
#-------------------------------------------------

## std sequence for tqaf layer1
process.load("TopQuarkAnalysis.TopObjectProducers.tqafLayer1_cff")

## necessary fixes to run 2.2.X on 2.1.X data
## comment this when running on samples produced
## with 22X
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)

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
                    'sisCone5CaloJets',                 # jet collection; must be already in the event when patLayer0 sequence is executed
                    layers           = [0,1],           # if you're not running patLayer1, set 'layers=[0]' 
                    runCleaner       = "CaloJet",       # =None if not to clean
                    doJTA            = True,            # run jet-track association & JetCharge
                    doBTagging       = True,            # run b-tagging
                    jetCorrLabel     = None,            # example jet correction name; set to None for no JEC
                    doType1MET       = True,            # recompute Type1 MET using these jets
                    genJetCollection = "sisCone5GenJets"                    
                    )

# now set JEC by hand
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
    fileName = cms.untracked.string('TQAFLayer1_Output.fromAOD_full.root')
)


#-------------------------------------------------
# output paths; in order not to write the
# persistent output to file comment the output
# path
#-------------------------------------------------

## output
process.outpath = cms.EndPath(process.out)
