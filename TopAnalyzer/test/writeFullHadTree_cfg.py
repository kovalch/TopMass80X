import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys,os

# setup 'standard' options
options = VarParsing.VarParsing ('standard')
## decide whether to run on:  * data *, * sig *, * bkg *, * qcd *, * allRedigi * or * all *
options.register('eventFilter', 'data', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "kind of data to be processed")
## choose whether to use PF or not
options.register('usePF'      , True , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool , "use PF for processing")
## set the jet energy scale variation type (abs, jes)
options.register('jesType'  ,  'jes' , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "jet energy scale type")
## set the jet energy scale factor for MC
options.register('jesFactor'  ,  1.0 , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.float, "jet energy scale factor for MC")
## set the jet energy resolution smear factors for MC (-1,0,+1,9,10,11) (down,default,up,now down,new default, new up)
options.register('jetEResol', 10 , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "jet energy resol smear for MC")
## set the resolution fractor for the kinematic fit
options.register('fitResol'   ,  10 , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "kinFit resolution factor")
## choose whether to write output to disk or not
options.register('writeOutput', False , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "write events surviving all cuts to disk")
## setup the ptHatFilter in case 'eventFilter' is chosen to be qcd
options.register('maxPtHat', 999999., VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.float, "maxPtHat to be processed")
options.register('minPtHat', 0.     , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.float, "minPtHat to be processed")
## include PDF uncertainty histograms / tree
options.register('pdfUn'   , 1 , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "include hists/tree for pdf uncertainty")
## choose which background estimate should be done (0: no estimate, 1: b-tag weighting, 2: event mixing)
options.register('backgroundEstimation', 0 , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "decide which background estimate should be done")
## run directly on AOD, do everything on the fly
options.register('runOnAOD', True , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "run on AOD, do everything on the fly")
## which b-tag algo and WP should be used
options.register('bTagAlgoWP', 'CSVT', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "which b-tag algo and WP should be used")
## do mva selection instead of old style selection
options.register('mvaSelection', True , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "mva selection")
## weight of events (should be used only for MC samples)
options.register('mcWeight', -1.0 , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.float, "MC sample event weight")
## which PU scenario for PU reweighting
options.register('PUscenario', '2012', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "PU distribution used for MC PUweight calculation")
## trigger results tag
options.register('triggerTag', 'HLT', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "tag of trigger Results")
## use the *totalKinematicsFilter*
options.register('useTotalKinFilter', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "use totalKinematicsFilter")
## do only skimming, no analysis to be done
options.register('onlySkimming', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "do only skimming, no analysis to be done")
## use the skimmed event content (usually only applies for event mixing)
options.register('useSkimmedEventContent', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "use skimmed event content")
## derive b-tagging efficiencies
options.register('bTagEfficiencyDetermination', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "derive b-tagging efficiencies")


# get and parse the command line arguments
if hasattr(sys, "argv") :
    for args in sys.argv :
        arg = args.split(',')
        for val in arg:
            val = val.split('=')
            if len(val)==2 :
                setattr(options, val[0], val[1])

print "eventFilter  . . . . . . . :", options.eventFilter
print "usePF  . . . . . . . . . . :", options.usePF
print "jesType  . . . . . . . . . :", options.jesType
print "jesFactor  . . . . . . . . :", options.jesFactor
print "jetEResol  . . . . . . . . :", options.jetEResol
print "fitResol . . . . . . . . . :", options.fitResol
print "writeOutput  . . . . . . . :", options.writeOutput
print "maxPtHat . . . . . . . . . :", options.maxPtHat
print "minPtHat . . . . . . . . . :", options.minPtHat
print "pdfUncertainty . . . . . . :", options.pdfUn
print "backgroundEstimation . . . :", options.backgroundEstimation
print "runOnAOD . . . . . . . . . :", options.runOnAOD
print "bTagAlgoWP . . . . . . . . :", options.bTagAlgoWP
print "mvaSelection . . . . . . . :", options.mvaSelection
print "mcWeight . . . . . . . . . :", options.mcWeight
print "PUscenario . . . . . . . . :", options.PUscenario
print "triggerTag . . . . . . . . :", options.triggerTag
print "useTotalKinFilter  . . . . :", options.useTotalKinFilter
print "onlySkimming . . . . . . . :", options.onlySkimming
print "useSkimmedEventContent . . :", options.useSkimmedEventContent
print "bTagEfficiencyDetermination:", options.bTagEfficiencyDetermination

## use the FullHadTreeWriter to produce a TTree with the desired information
if options.useSkimmedEventContent :
    process = cms.Process("FullHadTreeWriterSkimmed")
else :
    process = cms.Process("FullHadTreeWriter")

## configure message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# use only when PATTrigger is needed on data to suppress useless error messages
#process.MessageLogger.cerr.default.limit = 0

## define input
process.source = cms.Source("PoolSource",
                            ## add your favourite file here
                            fileNames = cms.untracked.vstring(
        #'file:patTuple.root',
        #'file:patTuple_selected.root',
        #'file:/tmp/eschliec/tmp.root',
        '/store/user/eschliec/Skim_03_Data01/Skim_MJP12D1_v1_data/patTuple_selected_2066.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/209D26E9-AEE1-E111-BAA6-0030487D5D8D.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/54145FEE-1AE2-E111-8B8E-003048C69408.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/8CA6B320-CFE1-E111-A04C-003048D2BB22.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/B6D6C237-DBE1-E111-B65F-002481E0DC4C.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/E88CC8A4-85E1-E111-9EFF-0030487F1A55.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/2289830F-84E1-E111-9B1A-0030487F1717.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/56A8FFAD-29E2-E111-BB50-003048F0E55A.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/8EFC5ECD-D3E1-E111-BC9F-0030487F1737.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/BA1CA24C-CFE1-E111-8991-003048D436D4.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/EA0F572E-93E1-E111-9521-0025901D4936.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/28609559-8CE1-E111-9836-003048F0E18E.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/5898B1DB-6CE1-E111-86F5-00266CFB8D74.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/92900233-DDE1-E111-A7F8-0030487D710F.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/BCD8A9A3-EEE1-E111-B148-003048C69402.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/ECCFBCDD-89E1-E111-96C9-0025901D4936.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/2ADEA5DB-00E2-E111-AF20-0030487E4EBF.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/5CEE2644-2AE2-E111-AAE8-0030487F9151.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/9636845B-06E2-E111-8F7B-0030487DA061.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/C0E610CD-CEE1-E111-AEF3-00266CF1074C.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/F097145C-E7E1-E111-A64E-0030487E0A2D.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/2E5C881C-75E1-E111-AC25-0025901D4AF0.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/624FA788-7CE1-E111-B1EC-0030487FA607.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/9899EDC7-8FE1-E111-ADBD-0025904B1446.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/C4169660-82E1-E111-BFA3-0030487F9297.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/F60F0652-C8E1-E111-BE82-002481E94B26.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/324A56C4-94E1-E111-87C2-003048C68A90.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/66279F97-99E1-E111-83BC-003048C68A8A.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/9A494832-27E2-E111-983A-0030487F1BCF.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/C6DA98EE-80E1-E111-8A16-002481E0D480.root',
        #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/FA1056FA-06E2-E111-BE82-0030487EBB25.root',
        #'/store/mc/Summer12_DR53X/TT_8TeV-mcatnlo/AODSIM/PU_S10_START53_V7A-v1/00003/7687BF3C-5E19-E211-B850-003048679180.root'
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/001A6A4D-E95B-E211-9C65-002590200B78.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/425D1632-885B-E211-914A-001E67396C9D.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/7E9D8C16-1B5B-E211-AB4C-002590200A1C.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/00CDCBFA-545B-E211-8A2A-001E672CC1E7.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/42A54E74-255B-E211-B67C-001E673975EE.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/803354FE-1A5A-E211-8019-003048673F54.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/023141E4-7B5B-E211-94CE-003048673FEA.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/442C6826-595A-E211-AFDD-001E673969FA.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/8076E800-125C-E211-B770-001E67397D55.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/02A75FFB-255B-E211-BA26-001E67396AC2.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/44B95CED-7D5B-E211-A945-001E67396D51.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/8204B281-9C5B-E211-BDE6-001E6739801B.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/0420CCF5-7F5B-E211-9F50-001E67396A22.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/44EEF52A-3B5B-E211-B538-003048D47A6A.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/8230EB0E-285B-E211-A269-002590200AB8.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/04911EE3-515B-E211-9584-002590200878.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/465ED53D-955B-E211-81FA-003048D460FA.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/82DAA0CE-3F5B-E211-950F-0025B3E05DBE.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/06029B85-235C-E211-BB9A-003048D462D2.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/46A0A04F-165C-E211-B098-002590200850.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/842161B1-855B-E211-A36E-003048D47A6A.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/06532BFD-245C-E211-AAC6-002590200A94.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/481ED799-3E5B-E211-9D98-001E673989DF.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/8454E820-105C-E211-A060-003048D45F5C.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/0815D818-995B-E211-B6BC-001E67396A22.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/48D4F772-1B5B-E211-84B8-003048D45FB8.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/84BAAD77-7A5B-E211-9AAB-0025B3E05CF2.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/086421D1-FC5B-E211-883F-002590200828.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/4A2FB095-1A5B-E211-9DB1-001E67396C9D.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/861550D4-1E5B-E211-BB1C-001E67397094.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/08C77995-0A5B-E211-90C2-0025B3E0661E.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/4A580CD8-335C-E211-A094-001E67397BC5.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/86748E6C-F25B-E211-9A70-002590200AC4.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/0A22AB6F-3D5B-E211-9F69-0025902008B8.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/4A952395-2C5C-E211-AF85-001E67397314.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/868DB624-735A-E211-8947-001E67397D00.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/0AF644BA-265A-E211-8C20-003048D46050.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/4AF904C2-0B5B-E211-81F2-0025902009B4.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/86E7A3E1-405B-E211-98D3-003048D47A14.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/0C5329EB-6F5B-E211-B8E7-001E67396D10.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/4C6C96AD-BE5B-E211-93EA-0030486740BA.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/888DCD34-F35B-E211-9B1D-002590200878.root',
        #'/store/data/Run2012C/MultiJet1Parked/AOD/part2_05Nov2012-v2/10001/0E359FB4-8D5B-E211-93DC-001E6739801B.root',
        #"/store/mc/Summer12_DR53X/TTJets_HadronicMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/00001/F67C1297-9618-E211-AA98-001EC9D8A8D4.root",
        #"/store/mc/Summer12_DR53X/TTJets_HadronicMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/00001/FA7398AB-A618-E211-A744-90E6BA19A227.root",
        #"/store/mc/Summer12_DR53X/TTJets_HadronicMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/00001/36680301-9D18-E211-BDCA-00259073E496.root",
        #"/store/mc/Summer12_DR53X/TTJets_HadronicMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/00001/6A8BEA61-B218-E211-A181-E0CB4E29C4D9.root",
        #"/store/mc/Summer12_DR53X/TTJets_HadronicMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/00001/665A64E6-E818-E211-B986-E0CB4E19F9BD.root"
        #"/store/user/eschliec/Skim_02_Data01/Skim_MJP12D2_v1_data/patTuple_selected_3266.root",
        #"/store/user/eschliec/Skim_02_Data01/Skim_MJP12D2_v1_data/patTuple_selected_3259.root",
        #"/store/user/eschliec/Skim_02_Data01/Skim_MJP12D2_v1_data/patTuple_selected_3260.root",
        #"/store/user/eschliec/Skim_02_Data01/Skim_MJP12D2_v1_data/patTuple_selected_3261.root",
        #"/store/user/eschliec/Skim_02_Data01/Skim_MJP12D2_v1_data/patTuple_selected_3262.root",
        #"/store/user/eschliec/Skim_02_Data01/Skim_MJP12D2_v1_data/patTuple_selected_3263.root",
        #"/store/user/eschliec/Skim_02_Data01/Skim_MJP12D2_v1_data/patTuple_selected_3264.root",
        #"/store/user/eschliec/Skim_02_Data01/Skim_MJP12D2_v1_data/patTuple_selected_3265.root",
        #"/store/user/eschliec/Skim_02_Data01/Skim_MJP12D2_v1_data/patTuple_selected_3266.root",
        #"/store/user/eschliec/Skim_02_Data01/Skim_MJP12D2_v1_data/patTuple_selected_3267.root",
        #"/store/user/eschliec/Skim_02_Data01/Skim_MJP12D2_v1_data/patTuple_selected_3268.root",
        #"/store/user/eschliec/Skim_02_Data01/Skim_MJP12D2_v1_data/patTuple_selected_3269.root",
        #"/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00002/009B74EF-A84A-E311-AEBD-848F69FD2997.root",
        #"/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00002/D0680DCA-B04A-E311-9C3C-848F69FD45A7.root",
        #"/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00002/40583575-B04A-E311-BA57-848F69FD298E.root",
        #"/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00002/5C77C9EF-B14A-E311-9C20-00266CF2506C.root",
        #"/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00002/E0D3C157-B24A-E311-A040-00266CF9B970.root",
        #"/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00002/4279F970-B44A-E311-A11E-001D09FDD91E.root",
        #"/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00002/A87145D0-B44A-E311-B10C-7845C4FC3C11.root",
        #"/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00002/74C5B553-B54A-E311-8817-7845C4FC3767.root"
),
                            skipEvents = cms.untracked.uint32(0)
                            )

## process / skip specific events
#process.source.eventsToProcess = cms.untracked.VEventRange('1:1-1:6','2:100-3:max')
#process.source.eventsToSkip = cms.untracked.VEventRange('1:1-1:6','2:100-3:max')
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('1:1-1:6','2:100-3:max')
#process.source.lumisToSkip = cms.untracked.VLuminosityBlockRange('1:1-1:6','2:100-3:max')
#process.source.eventsToProcess = cms.untracked.VEventRange('1:26916229','1:41631988','1:41647546')

## print memory infos to check for modules with memory leaks
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck", ignoreTotal = cms.untracked.int32(0))

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('UnimplementedFeature')
)

## register TFileService
process.TFileService = cms.Service("TFileService",
            fileName = cms.string('writeFullHadTree_test.root')
)

## configure geometry & conditions
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if(options.eventFilter=='data'):
    if os.getenv('CMSSW_VERSION').startswith('CMSSW_5_3_'):
        #process.GlobalTag.globaltag = cms.string('GR_R_53_V19::All')
        process.GlobalTag.globaltag = cms.string('FT53_V21A_AN6::All')
    elif os.getenv('CMSSW_VERSION').startswith('CMSSW_4_2_'):
        #process.GlobalTag.globaltag = cms.string('GR_R_38X_V15::All')
        #process.GlobalTag.globaltag = cms.string('GR_R_42_V14::All')
        #process.GlobalTag.globaltag = cms.string('FT_R_42_V18A::All')
        #process.GlobalTag.globaltag = cms.string('GR_R_42_V19::All')
        #process.GlobalTag.globaltag = cms.string('GR_R_42_V21A::All')
        process.GlobalTag.globaltag = cms.string('GR_R_42_V23::All')
        
else:
    if os.getenv('CMSSW_VERSION').startswith('CMSSW_5_3_'):
        process.GlobalTag.globaltag = cms.string('START53_V26::All')
    if os.getenv('CMSSW_VERSION').startswith('CMSSW_5_3_11'):
        process.GlobalTag.globaltag = cms.string('START53_V27::All')
    elif os.getenv('CMSSW_VERSION').startswith('CMSSW_4_2_'):
        #process.GlobalTag.globaltag = cms.string('START42_V12::All')
        #process.GlobalTag.globaltag = cms.string('START42_V13::All')
        #process.GlobalTag.globaltag = cms.string('START42_V15B::All')
        process.GlobalTag.globaltag = cms.string('START42_V17::All')
    elif os.getenv('CMSSW_VERSION').startswith('CMSSW_4_1_'):
        #process.GlobalTag.globaltag = cms.string('START38_V14::All')
        process.GlobalTag.globaltag = cms.string('START41_V0::All')
        

### save JEC and uncertainties to files
#process.readAK5PF = cms.EDAnalyzer( 'JetCorrectorDBReader'
#                                    # below is the communication to the database 
#                                  , payloadName    = cms.untracked.string('AK5PFchs')
#                                  #, payloadName    = cms.untracked.string('AK5Calo')
#                                    # this is used ONLY for the name of the printed txt files. You can use any name that you like, 
#                                    # but it is recommended to use the GT name that you retrieved the files from.
#                                  , globalTag      = cms.untracked.string('GR_R_42_V23')
#                                  , printScreen    = cms.untracked.bool(False)
#                                  , createTextFile = cms.untracked.bool(True)
#                                  )

## std sequence to produce the ttGenEvt
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

## filter for full-hadronic 
process.load("TopQuarkAnalysis.TopEventProducers.producers.TtDecaySelection_cfi")
process.ttFullHadronicFilter = process.ttDecaySelection.clone()

## ptHat filter
process.load("TopAnalysis.TopFilter.filters.PtHatFilter_cfi")

## additional jet energy smearing for MC
if options.jetEResol == 0:
    resolutions = [1.1, 1.1, 1.1]
    etaRanges = [0.0,1.5,1.5,2.0,2.0,-1.]
elif options.jetEResol == 1:
    resolutions = [1.2, 1.25, 1.3]
    etaRanges = [0.0,1.5,1.5,2.0,2.0,-1.]
elif options.jetEResol == -1:
    resolutions = [1.0, 0.95, 0.9]
    etaRanges = [0.0,1.5,1.5,2.0,2.0,-1.]
elif options.jetEResol == 9:
    resolutions = [0.990, 1.001, 1.032, 1.042, 1.089]
    etaRanges = [0.0,0.5,0.5,1.1,1.1,1.7,1.7,2.3,2.3,-1.]
elif options.jetEResol == 10:
    resolutions = [1.052, 1.057, 1.096, 1.134, 1.288]
    etaRanges = [0.0,0.5,0.5,1.1,1.1,1.7,1.7,2.3,2.3,-1.]
elif options.jetEResol == 11:
    resolutions = [1.115, 1.114, 1.161, 1.228, 1.488]
    etaRanges = [0.0,0.5,0.5,1.1,1.1,1.7,1.7,2.3,2.3,-1.]
else:
    exit('jetEResol: ' + str(options.jetEResol) + ' not supported, choose -1/0/+1/9/10/11')

from TopAnalysis.TopUtils.JetEnergyScale_cfi import scaledJetEnergy
process.scaledJetEnergy = scaledJetEnergy.clone( inputJets            = "selectedPatJetsAK5PF",
                                                 inputMETs            = "patMETsPF",
                                                 payload              = "AK5PFchs",
                                                 scaleFactor          = options.jesFactor,
                                                 scaleType            = "abs", #abs, rel, top:up, top:down, jes:up, jes:down
                                                 jetPTThresholdForMET = 20.,
                                                 jetEMLimitForMET     = 0.9,
                                                 resolutionFactors    = resolutions,
                                                 resolutionEtaRanges  = etaRanges,
                                                 JECUncSrcFile        = "TopAnalysis/TopUtils/data/Summer13_V5_DATA_UncertaintySources_AK5PFchs.txt",
                                                 )

## don't scale MET if skimmed event content is used
if options.useSkimmedEventContent :
    process.scaledJetEnergy.inputMETs = ""

## take proper JEC uncertainties for 2011 data
if os.getenv('CMSSW_VERSION').startswith('CMSSW_4_'):
    process.scaledJetEnergy.JECUncSrcFile = "TopAnalysis/TopUtils/data/JEC11_V10_AK5PF_UncertaintySources.txt"

## set set energy scaling factors
if options.jesType == 'jes' :
    if options.jesFactor > 1.0 :
        process.scaledJetEnergy.scaleType = "jes:up"
    elif options.jesFactor < 1.0 :
        process.scaledJetEnergy.scaleType = "jes:down"
## set set energy scaling factors for flavor ONLY
elif options.jesType == 'flavor' :
    process.scaledJetEnergy.sourceName = "FlavorPureBottom"
    if options.jesFactor > 1.0 :
        process.scaledJetEnergy.scaleType = "flavor:up"
    elif options.jesFactor < 1.0 :
        process.scaledJetEnergy.scaleType = "flavor:down"
## set energy scaling factors for a single source
elif options.jesType == 'CorrelationGroupMPFInSitu' or options.jesType == 'CorrelationGroupFlavor' or options.jesType == 'CorrelationGroupIntercalibration' or options.jesType == 'CorrelationGroupUncorrelated' or options.jesType == 'CorrelationGroupbJES' or options.jesType == 'TotalNoFlavor' or options.jesType == 'FlavorPureQuark' or options.jesType == 'PileUpPtBB' or options.jesType == 'PileUpPtEC':
    process.scaledJetEnergy.sourceName = options.jesType
    if options.jesFactor > 1.0 :
        process.scaledJetEnergy.scaleType = "source:up"
    elif options.jesFactor < 1.0 :
        process.scaledJetEnergy.scaleType = "source:down"

## residual jet corrector for data
from TopAnalysis.TopUtils.ResidualJetCorrector_cfi import residualCorrectedJets
process.residualCorrectedJets = residualCorrectedJets.clone()

if options.eventFilter=='data' :
    ## dummy sequence, will be removed later on in the process
    process.filterSequence = cms.Sequence(process.residualCorrectedJets
                                          )
    
elif options.eventFilter=='sig' :
    ## sequence with fullHad ttbar filter
    process.filterSequence = cms.Sequence(process.makeGenEvt *
                                          process.ttFullHadronicFilter *
                                          process.scaledJetEnergy
                                          )

elif options.eventFilter=='bkg' :
    ## invert fullHad filter
    process.ttFullHadronicFilter.invert = True
    ## sequence with non-fullHad ttbar filter
    process.filterSequence = cms.Sequence(process.makeGenEvt *
                                          process.ttFullHadronicFilter *
                                          process.scaledJetEnergy
                                          )

elif options.eventFilter=='qcd' :
    if options.maxPtHat<999999. :
        ## ptHat filter
        process.filterPtHat.maxPtHat = options.maxPtHat
        print "included ptHatFilter with 'maxPtHat' =",
        print process.filterPtHat.maxPtHat

    if options.minPtHat>0. :
        process.filterPtHat.minPtHat = options.minPtHat
        print "included ptHatFilter with 'minPtHat' =",
        print process.filterPtHat.minPtHat

    ## sequence with ptHat filter
    process.filterSequence = cms.Sequence(process.filterPtHat *
                                          process.scaledJetEnergy)

elif options.eventFilter=='all' :
    process.filterSequence = cms.Sequence(process.scaledJetEnergy)
    
elif options.eventFilter=='toyMC' :
    ## convert GenJets to PatJets
    process.load("TopAnalysis.TopUtils.convertGenToPatJets_cff")
    process.filterSequence = cms.Sequence(process.fullHadGenJetSelection*process.convertGenToPatJets)
    process.combinedSecondaryVertexBJetTags.jetsToTag = [2,3]
    
else:
    raise NameError, "'"+options.eventFilter+"' is not a proper eventFilter!"

# adapt output filename
process.TFileService.fileName = 'writeFullHadTree_test_'+options.eventFilter+'.root'

## fully hadronic selection
if options.mvaSelection:
    process.load("TopAnalysis.TopFilter.sequences.fullHadronicSelectionMVA_cff")
    from TopAnalysis.TopFilter.sequences.fullHadronicSelectionMVA_cff import *
else:
    process.load("TopAnalysis.TopFilter.sequences.fullHadronicSelection_cff")
    from TopAnalysis.TopFilter.sequences.fullHadronicSelection_cff import *

## do kinematic fit needed for fully hadronic selection
from TopQuarkAnalysis.TopEventProducers.sequences.ttFullHadEvtBuilder_cff import *

## add a second genMatch which does not care about ambiguities
cloneTtFullHadEvent(process)
process.ttFullHadJetPartonMatch2.algorithm = 'totalMinDist'

## add the kinematic fit to the event
addTtFullHadHypotheses(process,
                       ["kKinFit"]
                       )

## remove the genMatch if not running on signal events
if options.useSkimmedEventContent or not options.eventFilter=='sig' :
    removeTtFullHadHypGenMatch(process)

## load the background estimation modules
runAsBackgroundEstimation(process, options.backgroundEstimation)

## changing bTagger, possible are: TCHEL, TCHEM, TCHPM, TCHPT, SSVHEM, SSVHPT, CSVM, CSVT, CSVMVA
## CSVMVA has NO officialy blessed WP
if options.bTagAlgoWP == "TCHEM" :
    switchToTCHEM(process)
elif options.bTagAlgoWP == "TCHEL" :
    switchToTCHEL(process)
elif options.bTagAlgoWP == "TCHPM" :
    switchToTCHPM(process)
elif options.bTagAlgoWP == "TCHPT" :
    switchToTCHPT(process)
elif options.bTagAlgoWP == "SSVHEM" :
    switchToSSVHEM(process)
elif options.bTagAlgoWP == "SSVHPT" :
    switchToSSVHPT(process)
elif options.bTagAlgoWP == "CSVM" :
    switchToCSVM(process)
elif options.bTagAlgoWP == "CSVT" :
    switchToCSVT(process)
elif options.bTagAlgoWP == "CSVTSysUp" :
    switchToCSVTSysUp(process)
elif options.bTagAlgoWP == "CSVTSysDown" :
    switchToCSVTSysDown(process)
elif options.bTagAlgoWP == "CSVMVA" :
    switchToCSVMVA(process)
else :
    exit('bTagAlgoWP * ' + options.bTagAlgoWP + " * NOT SUPPORTED, STOP PROCESSING")
    
## selection should be run on PFJets instead of caloJets
if not options.usePF: 
    runOnCalo(process)

## if running on real data, do everything needed for this
if options.eventFilter=='data' :
    runOnData(process)
    ## needed because the patTriggerEvent is to chatty to be used and
    ## data is already skimmed with trigger requirement
    #removeTrigger(process)

## increase the resolution of the kinematic fit
if options.fitResol == 1 :
    fitEtaRages = [0.0,-1.0]
    fitResolutions = [1.1]
elif options.fitResol == 10 :
    fitEtaRages = [0.0,0.5,1.1,1.7,2.3,-1.]
    fitResolutions = [1.052, 1.057, 1.096, 1.134, 1.288]
else:
    exit('fitResol: ' + str(options.fitResol) + ' not supported, choose 1/10')
increaseKinFitResolution(process, fitEtaRages, fitResolutions)

## remove cutflow monitoring
removeMonitoringOfCutflow(process)

## THE treeWriter itself
process.load("TopAnalysis.TopAnalyzer.FullHadTreeWriter_cfi")

## register TreeRegistryService
process.load("TopMass.TopEventTree.TreeRegistryService_cfi")
process.TreeRegistryService.treeName  = "eventTree"
process.TreeRegistryService.treeTitle = "Tree for UHH top-quark analysis.\n//Particles are in order {TTBar0, Top1, Top2, W1, W2, B1, LightQ1, LightQBar1, B2, LightQ2, LightQBar2}"

## import QuarkGluonTagger
#from QuarkGluonTagger.EightTeV.QGTagger_RecoJets_cff import *
   
## load HypothesisAnalyzer
from TopMass.TopEventTree.EventHypothesisAnalyzer_cfi import analyzeHypothesis
process.analyzeKinFit = analyzeHypothesis.clone(hypoClassKey = "ttFullHadHypKinFit:Key", ttEvent = "ttFullHadEvent", ttEventGen2 = "ttFullHadEvent2", jets = "tightLeadingJets", gluonTagSrc = cms.InputTag('QGTagger', 'qgLikelihood'), maxNJets = cms.int32(20))
from TopMass.TopEventTree.JetEventAnalyzer_cfi import analyzeJets
process.analyzeJets = analyzeJets.clone(jets = "tightLeadingJets",
                                        alternativeJets = "selectedPatJetsAK5Calo",
                                        gluonTagSrc = cms.InputTag('QGTagger', 'qgLikelihood'))
from TopMass.TopEventTree.BRegJetEventAnalyzer_cfi import analyzeBRegJets
process.analyzeBRegJets = analyzeBRegJets.clone(jets = "tightLeadingJets")
from TopMass.TopEventTree.WeightEventAnalyzer_cfi import analyzeWeights
process.analyzeWeights = analyzeWeights.clone(jets = "tightLeadingJets",
                                              mcWeight = options.mcWeight,
                                              puWeightSrc     = cms.InputTag("eventWeightPU"       , "eventWeightPU"),
                                              puWeightUpSrc   = cms.InputTag("eventWeightPUSysUp"  , "eventWeightPU"),
                                              puWeightDownSrc = cms.InputTag("eventWeightPUSysDown", "eventWeightPU"),
                                              savePDFWeights = True)

## switch tree writer to appropriate jet source
process.FullHadTreeWriter = process.writeFullHadTree.clone(JetSrc = "tightLeadingJets")

## set correct MC weight
process.FullHadTreeWriter.MCweight = options.mcWeight

## calculate PU weight factor
process.load("TopAnalysis.TopUtils.EventWeightPU_cfi")
if os.getenv('CMSSW_VERSION').startswith('CMSSW_5_3_'):
    process.eventWeightPU = process.eventWeightPU.clone(MCSampleFile = "TopAnalysis/TopUtils/data/MC_PUDist_Summer12_S10.root", MCSampleHistoName = "puhisto")
elif os.getenv('CMSSW_VERSION').startswith('CMSSW_4_2_'):
    #process.eventWeightPU = process.eventWeightPU.clone(MCSampleFile = "TopAnalysis/TopUtils/data/MC_PUDist_Summer11_TTJets_TuneZ2_7TeV_madgraph_tauola.root")
    #process.eventWeightPU = process.eventWeightPU.clone(MCSampleFile = "TopAnalysis/TopUtils/data/Fall11/MC_PUDist_Fall11_TTJets_TuneZ2_7TeV_MadGraph.root")
    #rocess.eventWeightPU = process.eventWeightPU.clone(MCSampleFile = "TopAnalysis/TopUtils/data/MC_PUDist_Default2011.root", MCSampleHistoName = "histoSummer11_flat_observed")
     process.eventWeightPU = process.eventWeightPU.clone(MCSampleFile = "TopAnalysis/TopUtils/data/MC_PUDist_Default2011.root", MCSampleHistoName = "histo_Fall11_true")
    #process.eventWeightPU = process.eventWeightPU.clone(MCSampleFile = "TopAnalysis/TopUtils/data/MC_PUDist_WJets_Spring11.root", MCSampleHistoName = "pileup")
elif os.getenv('CMSSW_VERSION').startswith('CMSSW_4_1_'):
    process.eventWeightPU = process.eventWeightPU.clone(MCSampleFile = "TopAnalysis/TopUtils/data/MC_PUDist_WJets_Spring11.root", MCSampleHistoName = "pileup")

if options.PUscenario == '2012':
    process.eventWeightPU.DataFile = "TopAnalysis/TopUtils/data/MyPileup_HLT_QuadJet50_69400_SysNo.root"
    process.eventWeightPU.DataHistoName = "pileup"
    process.eventWeightPUSysUp   = process.eventWeightPU.clone(DataFile = "TopAnalysis/TopUtils/data/MyPileup_HLT_QuadJet50_72870_SysUp.root")
    process.eventWeightPUSysDown = process.eventWeightPU.clone(DataFile = "TopAnalysis/TopUtils/data/MyPileup_HLT_QuadJet50_65930_SysDown.root")
elif options.PUscenario == '11_178078':
    #process.eventWeightPU.DataFile = "TopAnalysis/TopUtils/data/MyDataPileupHistogram_precaleWeighted_new.root"
    #process.eventWeightPU.DataFile = "TopAnalysis/TopUtils/data/MyDataPileupHistogram_precaleWeighted_3545fb-1.root"
    #process.eventWeightPU.DataFile = "TopAnalysis/TopUtils/data/MyDataPileupHistogram_precaleWeighted_3545fb-1_-8.root"
    #process.eventWeightPU.DataFile = "TopAnalysis/TopUtils/data/MyDataPileupHistogram_precaleWeighted_3545fb-1_+8.root"
    #process.eventWeightPU.DataFile = "TopAnalysis/TopUtils/data/MyDataPileupHistogram_observed.root"
    process.eventWeightPU.DataFile = "TopAnalysis/TopUtils/data/MyDataPileupHistogram_true.root"
    process.eventWeightPU.DataHistoName = "pileup"
elif options.PUscenario == '11_167913':
#if options.PUscenario == '11_167913':
    process.eventWeightPU.DataFile = "TopAnalysis/TopUtils/data/Data_PUDist_160404-163869_7TeV_May10ReReco_Collisions11_v2_and_165088-167913_7TeV_PromptReco_Collisions11.root"
elif options.PUscenario == '11_166861':
    process.eventWeightPU.DataFile = "TopAnalysis/TopUtils/data/Data_PUDist_160404-166861_7TeV_PromptReco_Collisions11.root"
elif options.PUscenario == '11_166502':
    process.eventWeightPU.DataFile = "TopAnalysis/TopUtils/data/Data_PUDist_160404-166502_7TeV_PromptReco_Collisions11.root"
elif options.PUscenario == '11_May10':
    process.eventWeightPU.DataFile = "TopAnalysis/TopUtils/data/Data_PUDist_160404-163869_7TeV_May10ReReco_Collisions11.root"
elif options.PUscenario == '10_Apr21':
    process.eventWeightPU.DataFile = "TopAnalysis/TopUtils/data/Data_PUDist_136033-149442_7TeV_Apr21ReReco_Collisions10.root"
elif options.PUscenario == 'EPS':
    process.eventWeightPU.DataFile = "TopAnalysis/TopUtils/data/Data_PUDist_160404-163869_7TeV_May10ReReco_Collisions11_v2_and_165088-167913_7TeV_PromptReco_Collisions11.root"
elif not options.eventFilter == 'data':
    exit('PU SCENARIO * ' + options.PUscenario + " * NOT SUPPORTED, STOP PROCESSING")

## calculate btag weight factor
from TopAnalysis.TopUtils.BTagSFEventWeight_cfi import bTagSFEventWeight
process.bTagSFEventWeight = bTagSFEventWeight.clone(jets = "tightLeadingJets", bTagAlgo = "CSVT", version  = "2012", filename = "TopAnalysis/Configuration/data/analyzeBTagEfficiency"+options.bTagAlgoWP+".root", noHistograms = True)
process.bTagSFEventWeightBTagSFUp     = process.bTagSFEventWeight.clone(sysVar = "bTagSFUp")
process.bTagSFEventWeightBTagSFDown   = process.bTagSFEventWeight.clone(sysVar = "bTagSFDown")
process.bTagSFEventWeightMisTagSFUp   = process.bTagSFEventWeight.clone(sysVar = "misTagSFUp")
process.bTagSFEventWeightMisTagSFDown = process.bTagSFEventWeight.clone(sysVar = "misTagSFDown")

## calculate b tag efficiencies
process.load("TopAnalysis.TopAnalyzer.BTagEfficiencyAnalyzer_cfi")
## NOTE: process needs to be named bTagEff, so that BTagSFEventWeight.cc can find the histo
process.bTagEff = process.analyzeBTagEfficiency.clone(jets         = "tightLeadingJets",
                                                      bTagAlgo     = "combinedSecondaryVertexBJetTags",
                                                      bTagDiscrCut = 0.898, ## CVST
                                                      weight       = cms.InputTag("eventWeightPU", "eventWeightPU")
                                                     )
## remove PDF uncertainties
removePDFUncertainties(process)
if not options.pdfUn==0 :
    process.FullHadTreeWriter.DoPDFUncertainty = True

## get b-tag efficiency infos
process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")
process.load ("RecoBTag.PerformanceDB.BTagPerformanceDB1107")

## ---
##    run the final sequence
## ---
process.eventContentAnalyzer = cms.EDAnalyzer("EventContentAnalyzer")
process.p1 = cms.Path(## do the genEvent selection
                      process.filterSequence
                      ## do the filtering
                      *process.createJetCollections
                      *process.leadingJetSelection
                      *process.eventWeightPU
                      *process.eventWeightPUSysUp
                      *process.eventWeightPUSysDown
                      *process.bTagSFEventWeight
                      *process.bTagSFEventWeightBTagSFUp
                      *process.bTagSFEventWeightBTagSFDown
                      *process.bTagSFEventWeightMisTagSFUp
                      *process.bTagSFEventWeightMisTagSFDown
                      *process.makeTtFullHadEvent
                      ## write the tree
                      #*process.FullHadTreeWriter
                      *process.analyzeKinFit
                      *process.analyzeJets
                      #*process.analyzeBRegJets
                      *process.analyzeWeights
                      #*process.eventContentAnalyzer
                      )

if options.mvaSelection:
    process.p1.replace(process.leadingJetSelection,process.leadingJetSelection*process.mvaDisc)
    process.p1.remove(process.trackCountingHighEffBJets)
    process.p1.remove(process.tightBottomJets)

if options.backgroundEstimation != 0:
    if options.backgroundEstimation == 1:
        process.p1.remove(process.makeTtFullHadEvent)
        process.p1.remove(process.FullHadTreeWriter)
        if not options.eventFilter=='data' :
            process.analyzeFullHadQCDEstimation.bTagAlgoWP = "TCHEM40MC"
            process.analyzeFullHadQCDEstimation.MCweight = options.mcWeight
        
        if options.usePF:
            getattr(process, process.tightBottomJets.src.getModuleLabel()).src = cms.InputTag("goodJetsPF")
                     
        process.p1 += cms.Sequence( getattr(process, process.tightBottomJets.src.getModuleLabel())
                                  * process.tightBottomJets
                                  * process.tightBottomJetSelection
                                  * process.analyzeFullHadQCDEstimation
                                  )

        process.TFileService.fileName = 'QCDEstimation_1_'+options.eventFilter+'.root'

    if options.backgroundEstimation == 2:
        process.p1.remove(process.FullHadTreeWriter)
        #if options.usePF:
        #    getattr(process, process.tightBottomJets.src.getModuleLabel()).src = cms.InputTag("goodJetsPF")
        #             
        #if(options.eventFilter=='toyMC'):
        #    process.analyzeFullHadEventMixer.speedUp = 10.0
        #process.p1 += cms.Sequence( getattr(process, process.tightBottomJets.src.getModuleLabel())
        #                          * process.tightBottomJets
        #                          * process.tightBottomJetSelection
        #                          * process.analyzeFullHadEventMixer
        #                          )

        if hasattr(process,'scaledJetEnergy'):
            process.p1.remove(process.scaledJetEnergy)
        if hasattr(process,'goodJetsPF'):
            process.p1.remove(process.goodJetsPF)
        if hasattr(process,'tight4LeadingJets'):
            process.p1.remove(process.tight4LeadingJets)
        if hasattr(process,'tight5LeadingJets'):
            process.p1.remove(process.tight5LeadingJets)
        if hasattr(process,'tight6LeadingJets'):
            process.p1.remove(process.tight6LeadingJets)
        if hasattr(process,'tightLeadingJets'):
            process.p1.remove(process.tightLeadingJets)
        if hasattr(process,'leading4JetSelection'):
            process.p1.remove(process.leading4JetSelection)
        if hasattr(process,'leading5JetSelection'):
            process.p1.remove(process.leading5JetSelection)
        if hasattr(process,'leading6JetSelection'):
            process.p1.remove(process.leading6JetSelection)

        from TopMass.TopEventTree.JetEventMixer_cfi import mixJets
        process.tightLeadingJets = mixJets.clone()
        process.tightLeadingJets.input.fileNames = process.source.fileNames
        process.tightLeadingJets.nMix = 8
        process.tightLeadingJets.speedUp = 1
        #process.tightLeadingJets.input.skipEvents = cms.untracked.uint32(0) ## DOES NOT WORK
        process.source = cms.Source("EmptySource")
        #process.p1.replace(process.leading6JetSelection, process.leading6JetSelection*process.mixJets)
        #process.p1.replace(process.gen4JetSelector, process.mixJets*process.gen4JetSelector)
        process.eventContentAnalyzer = cms.EDAnalyzer("EventContentAnalyzer")
        process.p1.insert(0, process.tightLeadingJets)#*process.eventContentAnalyzer)
        #process.p1 = cms.Path(process.tightLeadingJets*process.eventContentAnalyzer)
        PUSource = cms.InputTag("tightLeadingJets", "addPileupInfo")
        process.eventWeightPU.PUSource = PUSource
        process.analyzeWeights.puSrc = PUSource
        
        process.analyzeJets.alternativeJets = cms.InputTag("tightLeadingJets","calo")

        process.TFileService.fileName = 'QCDEstimation_2_'+options.eventFilter+'.root'

if(options.pdfUn==2):
    print "Writing TTree for all", options.eventFilter, "events, no further selection done!"
    #process.p1.remove(process.analyseFullHadronicSelection)
    process.p1 = cms.Path(process.FullHadTreeWriter)
    process.FullHadTreeWriter.MultiJetMVADiscSrc = ""
    process.FullHadTreeWriter.PUweightSrc        = ""
    process.FullHadTreeWriter.ElectronSrc        = ""
    process.FullHadTreeWriter.FitSrc             = ""
    process.FullHadTreeWriter.METSrc             = ""
    process.FullHadTreeWriter.PUSrc              = ""
    process.FullHadTreeWriter.VertexSrc          = ""
    process.FullHadTreeWriter.MuonSrc            = ""
    process.FullHadTreeWriter.JetSrc             = ""
    process.FullHadTreeWriter.bTagVal            = cms.vstring()
    process.FullHadTreeWriter.bTagName           = cms.vstring()
    

if options.eventFilter=='data':
    process.p1.remove(process.filterSequence)
    process.p1.remove(process.eventWeightPU)
    process.p1.remove(process.eventWeightPUSysUp)
    process.p1.remove(process.eventWeightPUSysDown)
    process.p1.remove(process.bTagSFEventWeight)
    process.p1.remove(process.bTagSFEventWeightBTagSFUp)
    process.p1.remove(process.bTagSFEventWeightBTagSFDown)
    process.p1.remove(process.bTagSFEventWeightMisTagSFUp)
    process.p1.remove(process.bTagSFEventWeightMisTagSFDown)
else:
    # different energy resolution of jets in simulation (has to be propagated to MET)
    if(hasattr(process, 'goodJets') & hasattr(process, 'scaledJetEnergy')):
        process.goodJets.src   = cms.InputTag('scaledJetEnergy', 'selectedPatJets', process.name_())
        process.FullHadTreeWriter.METSrc = cms.InputTag('scaledJetEnergy', 'patMETs', process.name_())
    if(hasattr(process, 'goodJetsPF') & hasattr(process, 'scaledJetEnergy')):
        process.goodJetsPF.src = cms.InputTag('scaledJetEnergy', 'selectedPatJetsAK5PF', process.name_())
        process.FullHadTreeWriter.METSrc = cms.InputTag('scaledJetEnergy', 'patMETsPF', process.name_())
    # change input for toyMC
    if(options.eventFilter=='toyMC'):
        from TopAnalysis.TopUtils.convertGenToPatJets_cff import selectedPatJets as mySelectedPatJets
        process.selectedPatJets = mySelectedPatJets
        process.goodJetsPF.src = cms.InputTag('selectedPatJets', '', process.name_())
        process.goodJetsPF.cut = "abs(eta) < 2.4 & pt > 20."
        process.p1.remove(process.eventWeightPU)
        process.p1.remove(process.eventWeightPUSysUp)
        process.p1.remove(process.eventWeightPUSysDown)
        process.p1.remove(process.bTagSFEventWeight)
        process.p1.remove(process.bTagSFEventWeightBTagSFUp)
        process.p1.remove(process.bTagSFEventWeightBTagSFDown)
        process.p1.remove(process.bTagSFEventWeightMisTagSFUp)
        process.p1.remove(process.bTagSFEventWeightMisTagSFDown)
        #process.p1.remove(process.FullHadTreeWriter)
        process.FullHadTreeWriter.MultiJetMVADiscSrc = ""
        process.FullHadTreeWriter.PUweightSrc        = ""
        process.FullHadTreeWriter.ElectronSrc        = ""
        #process.FullHadTreeWriter.FitSrc             = ""
        process.FullHadTreeWriter.METSrc             = ""
        process.FullHadTreeWriter.PUSrc              = ""
        process.FullHadTreeWriter.VertexSrc          = ""
        process.FullHadTreeWriter.MuonSrc            = ""
        #process.FullHadTreeWriter.bTagVal            = cms.vstring()
        #process.FullHadTreeWriter.bTagName           = cms.vstring()

if options.useSkimmedEventContent:
    process.FullHadTreeWriter.ElectronSrc = ""
    process.FullHadTreeWriter.METSrc      = ""
    process.FullHadTreeWriter.VertexSrc   = ""
    process.FullHadTreeWriter.MuonSrc     = ""
    process.FullHadTreeWriter.GenSrc      = ""

    process.analyzeWeights.genEventSrc = ""

if(not options.pdfUn==2 and not options.eventFilter=='toyMC'):
    
     ## prepend the sequence needed to run on the NEW PAT tuples
    from TopAnalysis.TopUtils.usePatTupleWithParticleFlow_cff import prependPF2PATSequence
    pf2patOptions = {'runOnOLDcfg': True}
    pf2patOptions['pfIsoConeMuon'] = 0.3
    pf2patOptions['pfIsoConeElec'] = 0.3
    pf2patOptions['pfIsoValMuon'] = 0.2
    pf2patOptions['pfIsoValElec'] = 0.2
    pf2patOptions['electronIDs'] = ''
    pf2patOptions['excludeElectronsFromWsFromGenJets'] = True
    #pf2patOptions['noMuonTopProjection'] = True
    #pf2patOptions['noElecTopProjection'] = True
    pf2patOptions['METCorrectionLevel'] = 1 
    pf2patOptions['addResolutions'] = True
    pf2patOptions['resolutionsVersion'] = 'fall11'
    if not os.getenv('CMSSW_VERSION').startswith('CMSSW_5_3_'):
        pf2patOptions['JECEra'] = 'Jec11_V10'
        pf2patOptions['JECFile'] = 'TopAnalysis/TopUtils/data/Jec11_V10_L7PartonPrivate.db'
        if os.getenv('GC_CONF'):
            print "Running with GC, resetting address of JECFile!"
            pf2patOptions['JECFile'] = '../src/TopAnalysis/TopUtils/data/Jec11_V10_L7PartonPrivate.db'
        pf2patOptions['additionalJECLevels'] = ['L7Parton']
    if options.eventFilter=='data':
        pf2patOptions['runOnMC'] = False
    if options.runOnAOD:
        pf2patOptions['runOnAOD'] = True
    if options.writeOutput:
        pf2patOptions['switchOffEmbedding'] = False
    pf2patOptions['addGluonTags'] = True

    if not options.useSkimmedEventContent:
        prependPF2PATSequence(process, options = pf2patOptions)
        
        ## choose correct set of jec levels for MC and data
        if options.eventFilter=='data':
            jecLevels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']
            myGenJetCollection = None
        else:
            jecLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
            myGenJetCollection = cms.InputTag("ak5GenJets")

        if not options.bTagEfficiencyDetermination:
            ## add AK5Calo Jets
            preSeq = process.patDefaultSequence.copy()
            from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
            addJetCollection(process,cms.InputTag('ak5CaloJets'),
                 'AK5', 'Calo',
                 doJTA        = False,
                 doBTagging   = False,
                 jetCorrLabel = ('AK5Calo', cms.vstring(jecLevels)),
                 doType1MET   = False,
                 doL1Cleaning = False,
                 doL1Counters = False,
                 genJetCollection=myGenJetCollection,
                 doJetID      = True,
                 jetIdLabel   = "ak5",
                 outputModules = []
                 )
            for mod in preSeq.moduleNames() :
                process.patDefaultSequence.remove(getattr(process,mod))
        
            process.patJetsAK5Calo.addResolutions = False
            process.patJetsAK5Calo.resolutions = cms.PSet()
            process.patJetsAK5Calo.addTagInfos = False
            process.patJetsAK5Calo.userData.userFloats.src = cms.VInputTag("")
            process.patJetsAK5Calo.addDiscriminators = False
        
            process.selectedPatJetsAK5Calo.cut = 'pt > 10. & abs(eta) < 2.4'
            process.p1.replace(process.leadingJetSelection, process.leadingJetSelection*process.patDefaultSequence)
        
    ## load bugfix for wrongly generated MC files (MadGraph / Powheg + Pythia6)
    process.load("GeneratorInterface.GenFilters.TotalKinematicsFilter_cfi")
    ## test higher threshold to avoid bias
    process.totalKinematicsFilter.tolerance = 5.0

    ## adaptions (re-aranging of modules) to speed up processing
    pathnames = process.paths_().keys()
    for pathname in pathnames:
        ## move the ttGenEvent filter to the beginning of the sequence
        if options.eventFilter=='sig' or options.eventFilter=='bkg' :
            getattr(process, pathname).remove(process.ttFullHadronicFilter)
            getattr(process, pathname).remove(process.genEvt)
            getattr(process, pathname).remove(process.decaySubset)
            getattr(process, pathname).remove(process.initSubset)
            ## if skimmed event content is used TtGenEvent is not needed
            if not options.useSkimmedEventContent :
                #getattr(process, pathname).insert(0,process.ttFullHadronicFilter)
                getattr(process, pathname).insert(0,process.genEvt)
                getattr(process, pathname).insert(0,process.decaySubset)
                getattr(process, pathname).insert(0,process.initSubset)
        ## move the trigger to the beginning of the sequence
        getattr(process, pathname).remove(process.trigger)
        ## if skimmed event content is used triggering was done already
        if not options.useSkimmedEventContent:
            getattr(process, pathname).insert(0,process.trigger)
        if hasattr(process,'readAK5PF') :
            getattr(process, pathname).insert(0,process.readAK5PF)
        
        ## add bugfix for wrongly generated MC files
        if not options.eventFilter == 'data' and options.useTotalKinFilter == True:
            getattr(process, pathname).insert(0,process.totalKinematicsFilter)

        ## add filter for hcal laser events
        if options.eventFilter == 'data' and os.getenv('CMSSW_VERSION').startswith('CMSSW_5_') and not options.useSkimmedEventContent == True:
            process.hcalLaserFilter = hltHighLevel.clone( TriggerResultsTag = cms.InputTag('TriggerResults')
                                                        , HLTPaths = cms.vstring('user_step')
                                                        )
            getattr(process, pathname).insert(0,process.hcalLaserFilter)
 
        from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
        ## change tag of triggerResults
        if not options.triggerTag=='HLT' :
            massSearchReplaceAnyInputTag(getattr(process, pathname), cms.InputTag("TriggerResults","","HLT"), cms.InputTag("TriggerResults","",options.triggerTag))
            process.patTrigger.processName      = options.triggerTag
            process.patTriggerEvent.processName = options.triggerTag

        ## adaption for 38X MC files, they don't have PUSummaryInfo
        if os.getenv('CMSSW_VERSION').startswith('CMSSW_4_1_') and (options.triggerTag=='HLT' or options.triggerTag=='REDIGI38X') :
            getattr(process, pathname).remove(process.eventWeightPU)
         
        ## do the replaceings, usually done in PF2PAT, here (again) in case PF2PAT was not run
        massSearchReplaceAnyInputTag(getattr(process,pathname), 'selectedPatJetsAK5PF', 'selectedPatJets')
        massSearchReplaceAnyInputTag(getattr(process,pathname), 'patMETsPF', 'patMETs')
        massSearchReplaceAnyInputTag(getattr(process,pathname), cms.InputTag('scaledJetEnergy', 'selectedPatJetsAK5PF', process.name_()), cms.InputTag('scaledJetEnergy', 'selectedPatJets', process.name_()))

    ## if skimmed event content is used or
    ## bTag efficiency determination is done,
    ## many modules are not needed
    if not options.useSkimmedEventContent :
        if options.onlySkimming or options.bTagEfficiencyDetermination:
            if hasattr(process,'mvaDisc')                        :
                process.p1.remove(process.mvaDisc)                          
            if hasattr(process,'ttFullHadJetPartonMatch')        :
                process.p1.remove(process.ttFullHadJetPartonMatch)          
            if hasattr(process,'ttFullHadHypGenMatch')           :
                process.p1.remove(process.ttFullHadHypGenMatch)             
            if hasattr(process,'kinFitTtFullHadEventHypothesis') :
                process.p1.remove(process.kinFitTtFullHadEventHypothesis)   
            if hasattr(process,'ttFullHadHypKinFit')             :
                process.p1.remove(process.ttFullHadHypKinFit)               
            if hasattr(process,'ttFullHadEvent')                 :
                process.p1.remove(process.ttFullHadEvent)                   
            if hasattr(process,'ttFullHadJetPartonMatch2')       :
                process.p1.remove(process.ttFullHadJetPartonMatch2)         
            if hasattr(process,'ttFullHadHypGenMatch2')          :
                process.p1.remove(process.ttFullHadHypGenMatch2)            
            if hasattr(process,'ttFullHadEvent2')                :
                process.p1.remove(process.ttFullHadEvent2)                  
            if hasattr(process,'bTagSFEventWeight')              :
                process.p1.remove(process.bTagSFEventWeight)
            if hasattr(process,'bTagSFEventWeightBTagSFUp')      :
                process.p1.remove(process.bTagSFEventWeightBTagSFUp)
            if hasattr(process,'bTagSFEventWeightBTagSFDown')    :
                process.p1.remove(process.bTagSFEventWeightBTagSFDown)
            if hasattr(process,'bTagSFEventWeightMisTagSFUp')    :
                process.p1.remove(process.bTagSFEventWeightMisTagSFUp)
            if hasattr(process,'bTagSFEventWeightMisTagSFDown')  :
                process.p1.remove(process.bTagSFEventWeightMisTagSFDown)
            if hasattr(process,'FullHadTreeWriter')              :
                process.p1.remove(process.FullHadTreeWriter)                
            if hasattr(process,'analyzeKinFit')                  :
                process.p1.remove(process.analyzeKinFit)                    
            if hasattr(process,'analyzeJets')                    :
                process.p1.remove(process.analyzeJets)                      
            if hasattr(process,'analyzeWeights')                 :
                process.p1.remove(process.analyzeWeights)
            ## when skimming is done, remove PUweight, too
            if options.onlySkimming:
                if hasattr(process,'eventWeightPU'):
                    process.p1.remove(process.eventWeightPU)
                if hasattr(process,'eventWeightPUSysUp'):
                    process.p1.remove(process.eventWeightPUSysUp)
                if hasattr(process,'eventWeightPUSysDown'):
                    process.p1.remove(process.eventWeightPUSysDown)
                if hasattr(process,'TFileService'):
                    delattr(process, 'TFileService')
             ## when bTag efficiency is determined PUweight and bTagEff are needed
            if options.bTagEfficiencyDetermination:
                if hasattr(process,'initSubset'):
                    process.p1.remove(process.initSubset)
                if hasattr(process,'decaySubset'):
                    process.p1.remove(process.decaySubset)
                if hasattr(process,'genEvt'):
                    process.p1.remove(process.genEvt)
                if hasattr(process,'eventWeightPUSysUp'):
                    process.p1.remove(process.eventWeightPUSysUp)
                if hasattr(process,'eventWeightPUSysDown'):
                    process.p1.remove(process.eventWeightPUSysDown)
                process.p1 += process.bTagEff
                process.TFileService.fileName = 'bTagEff_test.root'


#process.kinFitTtFullHadEventHypothesis.constraints = [1,2,3,4]
#process.kinFitTtFullHadEventHypothesis.mTop = 190

## Output Module Configuration
if options.writeOutput:
    process.outModule = cms.OutputModule("PoolOutputModule",
                                         fileName = cms.untracked.string('patTuple_selected.root'),
                                         # save only events passing the full path
                                         SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p1') ),
                                         # drop meta data for dropped events
                                         dropMetaData = cms.untracked.string('DROPPED'),
                                         # save output (comment to keep everything...)
                                         outputCommands = cms.untracked.vstring('drop *',
                                                                                #'keep *_generator_*_*',
                                                                                #'keep *_TriggerResults_*_HLT',
                                                                                'keep *_addPileupInfo_*_*',
                                                                                #'keep *_genEvt_*_FullHadTreeWriter',
                                                                                #'keep *_decaySubset_*_FullHadTreeWriter',
                                                                                #'keep *_initSubset_*_FullHadTreeWriter',
                                                                                #'keep *_ttFullHadEvent_*_FullHadTreeWriter',
                                                                                #'keep *_mvaDisc_*_*',
                                                                                #'keep *_eventWeightPU_*_*',
                                                                                #'keep *_selectedPatJets_*_*',
                                                                                'keep *_tightLeadingJets_*_*',
                                                                                'keep *_selectedPatJetsAK5Calo__*',
                                                                                #'keep *_selectedPatElectrons_*_*',
                                                                                #'keep *_selectedPatMuons_*_*',
                                                                                #'keep *_patMETs_*_*',
                                                                                'keep *_goodOfflinePrimaryVertices_*_*',
                                                                                )
                                         )
    process.outpath = cms.EndPath(process.outModule)

