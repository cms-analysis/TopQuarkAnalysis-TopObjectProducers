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
  # PAT test sample for 2.2.X
   #'file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_2_X_2008-11-03-STARTUP_V7-AODSIM.100.root'
  # PAT test sample for 2.1.X
   #'file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_1_X_2008-07-08_STARTUP_V4-AODSIM.100.root'
  # 2110 RelVal sample
  #'/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/0609A88C-F69A-DD11-AE42-001731AF6A8D.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/08DB9384-FD9A-DD11-ACFD-003048678DD6.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/12E0A583-F69A-DD11-9B03-003048678B38.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/1AFBB319-F89A-DD11-B491-001A928116C8.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/1CC22E59-069B-DD11-80F0-0018F3D09634.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/1EDF3F84-F69A-DD11-BEDA-001A928116CC.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/2864381E-F89A-DD11-9E2D-001A92810AAA.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/28B5240A-FC9A-DD11-A178-003048678C3A.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/2EB12E11-F89A-DD11-8D9A-003048678A7E.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/362701E1-4A9B-DD11-A9A4-001A92971BDA.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/446FE812-099B-DD11-B3CC-001A92810AEE.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/447DCD63-069B-DD11-8E7A-001731AF6A89.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/52349E81-FD9A-DD11-AF48-00304867918E.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/62085221-F89A-DD11-8F60-001731AF6BD3.root',
  #'/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/6419021D-F89A-DD11-8F79-0018F3D096AE.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/669F3A16-F89A-DD11-805E-003048679000.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/6A476935-F99A-DD11-8D7C-003048678AF4.root',
  #'/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/76815C20-F89A-DD11-A082-001731230A77.root',
  #'/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/7AF9CE67-359B-DD11-8EC2-001731A283E1.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/7E9AE61A-F89A-DD11-BEC2-0018F3D096C2.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/8E423E15-F89A-DD11-BEA4-0030486790C0.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/9097AC1D-F89A-DD11-937E-001A928116B2.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/90B8C28C-019B-DD11-AD6C-003048678B1C.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/9421B73D-F99A-DD11-97F1-0018F3D096C2.root',
  #'/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/984A721A-F89A-DD11-A351-0018F3D09708.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/9CB1827D-F69A-DD11-B362-003048678B0C.root',
  #'/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/9CDF9A1C-F89A-DD11-A0D8-001A92810AD8.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/A8D04310-099B-DD11-9A5D-001A928116E2.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/ACFF381C-F89A-DD11-9ABB-001A92971B08.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/B07E4F36-F99A-DD11-9C17-001A92971B9C.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/B2D4E439-F99A-DD11-B2C8-0018F3D09650.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/B42DD017-F89A-DD11-911E-001A92971B3A.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/B88DDA1C-F89A-DD11-A438-0018F3D096D4.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/C4CE1E88-F69A-DD11-8951-001A92811726.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/CE04991C-F89A-DD11-BB9E-001A92971B28.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/CEF0CA1C-F89A-DD11-ACDE-001A92971B8C.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/D253EE10-FC9A-DD11-A409-001A92971ACE.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/DA914143-F99A-DD11-8859-0017312B5A75.root',
  #'/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/DADC371C-F89A-DD11-8F20-001A92971B26.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/E8A52D35-F99A-DD11-8A52-0018F3D09642.root',
  #'/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/FCFD2C27-009B-DD11-BAC9-001A92971B9C.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/FECC8DE7-039B-DD11-89CD-003048769E65.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0002/04FB6353-0B9B-DD11-9FDF-001731AF66A5.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0002/1A7C8912-0D9B-DD11-92A4-00304876A13B.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0002/2C2D0EF5-5F9B-DD11-BF60-0018F3D095FA.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0002/2EA0E20C-099B-DD11-B19A-001A92971ACE.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0002/381C3CF2-579B-DD11-8620-0018F3D09702.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0002/427F0A35-489B-DD11-ACEA-003048767DF9.root'
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0002/747C2956-0B9B-DD11-B226-001731AF67EF.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0002/7C492256-0B9B-DD11-8C48-001731AF66C1.root',
   '/store/relval/CMSSW_2_1_10/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0002/F87A0845-0B9B-DD11-BB8A-003048678BAC.root'
    )
)

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
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
# patTuple configuration
#-------------------------------------------------

## std sequence for tqaf layer1
process.load("TopQuarkAnalysis.TopObjectProducers.patTuple_cff")

## necessary fixes to run 2.2.X on 2.1.X data
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)

## switch from clusters to rec hits in electron isolation
from PhysicsTools.PatAlgos.recoLayer0.electronIsolation_cff import useElectronRecHitIsolation
useElectronRecHitIsolation(process)

#-------------------------------------------------
# process paths;
#-------------------------------------------------

## process path
process.p = cms.Path(process.patTuple)


#-------------------------------------------------
# pat tuple event content; first ALL objects
# are dropped in this process; then patTuple
# content is added
#-------------------------------------------------

## define pat tuple event content
from TopQuarkAnalysis.TopObjectProducers.patTuple_EventContent_cff import *
makePatTupleEventContent(process)

## change jet collection
from PhysicsTools.PatAlgos.tools.jetTools import *

switchJetCollection(process, 
                    'sisCone5CaloJets',             # jet collection; must be already in the event when patLayer0 sequence is executed
                    layers       = [0,1],           # if you're not running patLayer1, set 'layers=[0]' 
                    runCleaner   = "CaloJet",       # =None if not to clean
                    doJTA        = True,            # run jet-track association & JetCharge
                    doBTagging   = True,            # run b-tagging
                    jetCorrLabel = ('SC5', 'Calo'), # example jet correction name; set to None for no JEC
                    doType1MET   = True             # recompute Type1 MET using these jets
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
    process.patTupleEventContent,
    verbose = cms.untracked.bool(True),
    dropMetaDataForDroppedData = cms.untracked.bool(True),                           
    fileName = cms.untracked.string('/afs/cern.ch/user/r/rwolf/pccmsuhh06/testPatTuple.root')
##  fileName = cms.untracked.string('testPatTuple.root')
)


#-------------------------------------------------
# output paths; in order not to write the
# persistent output to file comment the output
# path
#-------------------------------------------------

## output
process.outpath = cms.EndPath(process.out)
