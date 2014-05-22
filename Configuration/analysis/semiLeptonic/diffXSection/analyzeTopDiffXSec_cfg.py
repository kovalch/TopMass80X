import FWCore.ParameterSet.Config as cms

## ---
##   use this file to study different distributions for measurement of differential Cross Section
##   does also a monitoring of the cuts before they are applied
## ---

## ---
##    example to run this cfg file
## ---
## cmsRun analyzeTopDiffXSec_cfg.py triggerTag=HLT,sample=ttbar,lepton=muon,mctag=Summer12,eventsToProcess=5000

## ---
##    options
## ---

## introduce interface with crab file
## allows to choose registered parameter for every job individually
import FWCore.ParameterSet.VarParsing as VarParsing
import sys
import os
options = VarParsing.VarParsing ('standard')

# create object triggerTag with default value HLT of type singleton and string
options.register('triggerTag', 'HLT',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "chosen trigger tag")
# create sample label with default value data
# for Summer12 MC one can choose: ttbar, wjets, zjets, singleAntiTopS, singleTopT, singleAntiTopT, singleTopTw, singleAntiTopTw, singleTopS WW, WZ, ZZ, zprime_m1000gev_w100gev, qcd (for muon channel); qcdEM1, qcdEM2, qcdEM3, qcdBCE1, qcdBCE2, qcdBCE3 (for electron channel)
# for systematic samples see list for each MC sample
options.register('sample', 'none',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "chosen sample")
# create lepton channel label 
options.register('lepton', 'unset',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "chosen decay channel")
# create label to select version of MC samples (Summer11, Fall11)
# a) important for trigger
# b) also used to automatically select the mode for PU event reweighting
options.register('mctag', 'Summer12',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "label for version of MC samples (Summer12, Fall11, Summer11")
# create variable to indicate number of processed events
# -42 means nothing is changed wrt number typed below
options.register('eventsToProcess', -42,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "events to process")
# create label to apply/deactivate MC event weights (needed for MC@NLO)
options.register('MCweighting', 'unset',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "processing MC weights?")
# create tag for jet energy variations (JERUp, JERDown, JESUp, JESDown)
options.register('JEtag', 'none',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "chosen jet energy tag")
# closure test shape distortions
options.register('sysDistort', 'unset',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "ttbar shape distortion")
# value for top mass constraint in KinFit
options.register('massfix', 172.5 ,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "top mass constraint used in double kinFit approach")

# define the syntax for parsing
# you need to enter in the cfg file:
# search for arguments entered after cmsRun
if( hasattr(sys, "argv") ):
    # split arguments by comma - seperating different variables
    for args in sys.argv :
        arg = args.split(',')
        # split further by = to separate variable name and value
        for val in arg:
            val = val.split('=')
            # set variable var to value val (expected crab syntax: var=val)
            if(len(val)==2):
                setattr(options,val[0], val[1])

## choose jet collection and corresponding MET
if(not globals().has_key('jetType')):
    jetType =  'particleFlow' # 'Calo'

## run PF2PAT?
## only possible for special pat tuples!!!
if(not globals().has_key('pfToPAT')):
    pfToPAT = True #False

## choose the semileptonic decay channel (electron or muon)
#decayChannel=options.lepton
if(options.lepton=='unset'): 
    if(not globals().has_key('decayChannel')):
        decayChannel = 'muon' # 'electron'
else:
    decayChannel=options.lepton

# switch to run on data and remove all gen plots (type 'MC' or 'data')
if(not globals().has_key('runningOnData')): 
    runningOnData = "MC"
if("data" in options.sample):
    runningOnData = "data"
    
# print chosen sample (e.g. ttbar)
# value is known from external parsing
# if set, switches runOnAOD in PF2PAT to true
#print "Chosen sample to run over: ", options.sample

## choose JSON file for data
if(not globals().has_key('jsonFile')):
    jsonFile =  ''

## choose whether synchronisation exercise should be done
if(not globals().has_key('cutflowSynch')): 
    cutflowSynch  = False # True
    if(options.sample=="synch"):
        cutflowSynch  = True

## enable/ disable PU event reweighting
if(not globals().has_key('PUreweigthing')):
    PUreweigthing = True # False
    # take care of data
    if (cutflowSynch or not runningOnData == "MC"):
        PUreweigthing = False

## enable/ disable btag SF event reweighting
if(not globals().has_key('BtagReweigthing')):
    BtagReweigthing = True # False
    # take care of data
    if (cutflowSynch or not runningOnData == "MC"):
        BtagReweigthing = False

## choose b tag algo
if(not globals().has_key('bTagAlgo')):
    bTagAlgo =  'combinedSecondaryVertexBJetTags' # 'simpleSecondaryVertexHighEffBJetTags'

## choose b tag working point (discriminator cut value)
if(not globals().has_key('bTagDiscrCut')):
    bTagDiscrCut =  0.679  # 1.74 for SSVHEM

## enable/ disable efficiency SF event reweighting
if(not globals().has_key('effSFReweigthing')):
    effSFReweigthing = True # False
    # take care of data
    if (cutflowSynch or not runningOnData == "MC"):
        effSFReweigthing = False

# choose jet correction level shown in plots
# L3Absolute for MC, L2L3Residual for data
if(not globals().has_key('corrLevel')):
    corrLevel='L3Absolute'
if("data" in options.sample):
    corrLevel="L2L3Residual"

## run kinematic fit?
## ATTENTION: until the new parameter jetResolutionSmearFactor
## is implemented in a higher version of the TKinFitter package you
## need to check out the head:
## cvs co TopQuarkAnalysis/TopKinFitter

if(not globals().has_key('applyKinFit')):
    applyKinFit = True # False
    if(cutflowSynch):
        applyKinFit = False
        
## choose whether you want to use the double KinFit approach
## fit1: mtop=172.5GeV fix -> jet permutation
## -> fit2: mtop free, jet permutation from fit1 -> result                
doubleKinFit = True
if(applyKinFit==False):
    doubleKinFit=False

## choose whether you want to do the probability cut
## already before filling any result plot
ProbCutByDefault = False
if(applyKinFit==False):
    ProbCutByDefault=False
    
## choose whether you want a pat tuple as output
if(not globals().has_key('writeOutput')): 
    writeOutput  = False # True
    
## remove all ttbar specific gen level filter - used by other _cfg based on this file
if(not globals().has_key('removeGenTtbar')):
    removeGenTtbar = False 

## eventfilter is to select a special ttbar decay channel from ttbarSample by genmatching (ttbar MC only, other MC: choose 'all')
if(not globals().has_key('eventFilter')):
    eventFilter  = 'signal only' # 'background only' # 'all' # 'signal only' # 'semileptonic electron only' # 'dileptonic electron only' # 'dileptonic muon only' # 'fullhadronic' # 'dileptonic muon + electron only' # 'via single tau only' # 'dileptonic via tau only'
if("BG" in options.sample):
    eventFilter='background only'
    print "ttbar decay subset filter is inverted semileptonic muon or electron decay"
if(not "ttbar" in options.sample and not "powheg" in options.sample and not "perugia" in options.sample and not "mcatnlo" in options.sample and not "zprime" in options.sample):
    removeGenTtbar = True
    eventFilter='all'
if (cutflowSynch):
    eventFilter  = 'all'
    
## choose if you want to include tau (->e/mu) events
if(not globals().has_key('tau')):
    tau = False # True
# use these tau events only in SG/SB configuration not for single decay modes or all
if(not eventFilter=='signal only' and not eventFilter=='background only'):
    tau= False

# choose if you want to have only the gen paths
# and use the full statistics sample
if(not globals().has_key('genFull')):
    genFull=False
# makes no sense for non ttbar signal
if(not eventFilter=='signal only' or not options.sample=="ttbar"):
    genFull=False

# choose if you want to reduce the content
# (removes some control plots)
if(not globals().has_key('reduced')):
    reduced=True
    
## choose whether additional event weights should be applied
if(options.JEtag=="none"):
    additionalEventWeights  = True
else:
    additionalEventWeights  = False

if(not globals().has_key('additionalEventWeights')): 
    additionalEventWeights  = True

## enable/ disable systematic ttbar shape distortion event reweighting
if(not options.sysDistort=='unset'):
    sysDistort=options.sysDistort
if(not globals().has_key('sysDistort')):
    sysDistort =  ''
    #sysDistort =  'data' // -> topPt wrt measurement
    #sysDistort =  'ttbarMassUp'   
    #sysDistort =  'ttbarMassDown'
    #sysDistort =  'topPtUp'   
    #sysDistort =  'topPtDown'
# only done for ttbar
if(not "ttbar" in options.sample and not "mcatnlo" in options.sample and not "powheg" in options.sample and not "perugia" in options.sample ):
    sysDistort=''
# coupled to PU weight, therefore not applicable without
if(not PUreweigthing):
    sysDistort=''
## enable/ disable MC event weights (needed only for MC@NLO)
MCweighting = False # True
if("mcatnlo" in options.sample):
    MCweighting = True
if(not options.MCweighting=='unset'):
    MCweighting=options.MCweighting
# take care of data
if (not runningOnData == "MC"):
    MCweighting = False

## print gen particle list?
printGenParticleList=False
if(not "ttbar" in options.sample):
    printGenParticleList=False

# differetial xSec Analysis
process = cms.Process("topDifferentialXSec")

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.MessageLogger.categories.append('TtSemiLepKinFitter')
process.MessageLogger.categories.append('KinFitter')
process.MessageLogger.categories.append('GenCandSelector')
process.MessageLogger.cerr.TtSemiLepKinFitter = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)
process.MessageLogger.cerr.KinFitter = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)
process.MessageLogger.cerr.GenCandSelector = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)
## print memory infos to check for modules with memory leaks
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck", ignoreTotal = cms.untracked.int32(0)) 

## define input
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(    
    ## add your favourite file here
    #'/store/data/Run2011A/ElectronHad/AOD/PromptReco-v4/000/165/093/2C186C6C-C27F-E011-A1C7-001617E30F58.root'
    #'/store/data/Run2011A/MuHad/AOD/PromptReco-v4/000/165/205/0C569F2A-D382-E011-B122-00304879FBB2.root' 
    #'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/FEEE3638-F297-E011-AAF8-00304867BEC0.root'
    #'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/1204EE92-F397-E011-99E8-003048679012.root'
    # '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/FCB7FB42-ACE1-E111-B8AB-0025901D4936.root'
    #'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/DA4A7A9F-8DE1-E111-82D9-002481E0D480.root'
    #'/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/00001/FC5CECAE-8B14-E211-8578-0025B3E0652A.root'
    )
)

## create output file name
outputFileName = ""
if(decayChannel=="muon"):
    outputFileName = "muon"
elif(decayChannel=="electron"):
    outputFileName = "elec"
## will be continued later

## Flag required to compensate bug in energy-momentum conservation for non Pythia samples
## Modified in the following depending on the input sample, but only for the Pythia samples!
## Default value: False

PythiaSample="False"
usedSample="none"

## automatically load the correct (AOD) .root file list for each MC sample
if(not options.sample=="none"):
    outputFileName+="DiffXSec"
    if(options.sample=="dataABCD"):
        if(decayChannel=="muon"):
            usedSample="TopAnalysis/Configuration/samples/SingleMu_Run2012ABCDJan22ReReco_cff"
        elif(decayChannel=="electron"):
            usedSample="TopAnalysis/Configuration/samples/SingleElectron_Run2012ABCDJan22ReReco_cff"
        outputFileName+="Data2012ABCDJanReReco"
    elif(options.sample=="dataA"):
        if(decayChannel=="muon"):
            usedSample="TopAnalysis/Configuration/samples/SingleMu_Run2012AJan22ReReco_cff"
        elif(decayChannel=="electron"):
            usedSample="TopAnalysis/Configuration/samples/SingleElectron_Run2012AJan22ReReco_cff"
        outputFileName+="Data2012AJanReReco"
    elif(options.sample=="dataB"):
        if(decayChannel=="muon"):
            usedSample="TopAnalysis/Configuration/samples/SingleMu_Run2012BJan22ReReco_cff"
        elif(decayChannel=="electron"):
            usedSample="TopAnalysis/Configuration/samples/SingleElectron_Run2012BJan22ReReco_cff"
        outputFileName+="Data2012BJanReReco"
    elif(options.sample=="dataC"):
        if(decayChannel=="muon"):
            usedSample="TopAnalysis/Configuration/samples/SingleMu_Run2012CJan22ReReco_cff"
        elif(decayChannel=="electron"):
            usedSample="TopAnalysis/Configuration/samples/SingleElectron_Run2012CJan22ReReco_cff"
        outputFileName+="Data2012CJanReReco"
    elif(options.sample=="dataD"):
        if(decayChannel=="muon"):
            usedSample="TopAnalysis/Configuration/samples/SingleMu_Run2012DJan22ReReco_cff"
        elif(decayChannel=="electron"):
            usedSample="TopAnalysis/Configuration/samples/SingleElectron_Run2012DJan22ReReco_cff"
        outputFileName+="Data2012DJanReReco"
    elif("ttbar" in options.sample):
        # limited statistics (wo spin correlation inclusive sample)
        #usedSample="TopAnalysis/Configuration/Summer12/TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff" #////
        usedSample="TopAnalysis/Configuration/Summer12/TTJets_MSDecays_central_TuneZ2star_8TeV_madgraph_tauola__Summer12_DR53X_PU_S10_START53_V19_v1_cff"
        # full statistics: FIXME not available yet
        if(genFull):
            usedSample="TopAnalysis/Configuration/Summer12/TTJets_SemiLeptMGDecays_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
	if(eventFilter=='signal only'):
	    outputFileName+="Sig"
	elif(eventFilter=='background only'):
	    outputFileName+="Bkg"
	if(sysDistort!=""):
            if(sysDistort!="data" and sysDistort!="topptTTsys"):
                additionalEventWeights=False # if set to false no variations (SF+-, ...) are done
	    outputFileName+="SysDistort"+sysDistort
        if("ttbarMatchingDown" in options.sample):        
            usedSample="TopAnalysis/Configuration/Summer12/TTJets_MSDecays_matchingdown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V19_v2_cff"
            additionalEventWeights=False
            outputFileName+="MatchDown"
        elif("ttbarMatchingUp" in options.sample):        
            usedSample="TopAnalysis/Configuration/Summer12/TTJets_MSDecays_matchingup_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V19_v2_cff"
            additionalEventWeights=False
            outputFileName+="MatchUp"
        elif("ttbarScaleDown" in options.sample):        
            usedSample="TopAnalysis/Configuration/Summer12/TTJets_MSDecays_scaledown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V19_v1_cff"
            additionalEventWeights=False
            outputFileName+="ScaleDown"
        elif("ttbarScaleUp" in options.sample):        
            usedSample="TopAnalysis/Configuration/Summer12/TTJets_MSDecays_scaleup_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V19_v1_cff"
            additionalEventWeights=False
            outputFileName+="ScaleUp"
        elif("ttbarMassDown3" in options.sample):        
            usedSample="TopAnalysis/Configuration/Summer12/TTJets_MSDecays_mass166_5_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V19_v1_cff" 
            additionalEventWeights=False
            outputFileName+="TopMassDown3"
        elif("ttbarMassDown2" in options.sample):        
            usedSample="TopAnalysis/Configuration/Summer12/TTJets_MSDecays_mass169_5_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V19_v1_cff"
            additionalEventWeights=False
            outputFileName+="TopMassDown2"
        elif("ttbarMassDown" in options.sample):        
            usedSample="TopAnalysis/Configuration/Summer12/TTJets_MSDecays_mass171_5_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V19_v1_cff"
            additionalEventWeights=False
            outputFileName+="TopMassDown"
        elif("ttbarMassUp3" in options.sample):        
            usedSample="TopAnalysis/Configuration/Summer12/TTJets_MSDecays_mass178_5_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V19_v1_cff"
            additionalEventWeights=False
            outputFileName+="TopMassUp3"
        elif("ttbarMassUp2" in options.sample):        
            usedSample="TopAnalysis/Configuration/Summer12/TTJets_MSDecays_mass175_5_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V19_v1_cff"
            additionalEventWeights=False
            outputFileName+="TopMassUp2"
        elif("ttbarMassUp" in options.sample):        
            usedSample="TopAnalysis/Configuration/Summer12/TTJets_MSDecays_mass173_5_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V19_v1_cff"
            additionalEventWeights=False
            outputFileName+="TopMassUp"
        elif("MadGraph" in options.sample):
            usedSample="TopAnalysis/Configuration/Summer12/TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
            additionalEventWeights=False
            outputFileName+="MadGraph"
        elif("ttbar"==options.sample or "ttbarBG"==options.sample):        
            usedSample="TopAnalysis/Configuration/Summer12/TTJets_MSDecays_central_TuneZ2star_8TeV_madgraph_tauola__Summer12_DR53X_PU_S10_START53_V19_v1_cff" 
            additionalEventWeights=False
            outputFileName+="MadSpin"
        if(not options.massfix==172.5):        
            outputFileName+="TopMassConstraint"
            outputFileName+=((str(options.massfix)).replace(".", "p"))
    elif(options.sample=="synch"):
        usedSample="TopAnalysis/Configuration/Summer12/TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_synch2_cff"
        outputFileName+="Synch"
    elif("perugianoCR" in options.sample):       
        usedSample="TopAnalysis/Configuration/Summer12/TTJets_SemiLeptMGDecays_TuneP11noCR_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V19_v1_cff"
        if(eventFilter=='signal only'):
            outputFileName+="SigPerugianoCR"
        elif(eventFilter=='background only'):
            outputFileName+="BkgPerugianoCR"
        if(sysDistort!=""):
            if(sysDistort!="data"):
                additionalEventWeights=False # if set to false no variations (SF+-, ...) are done
            outputFileName+="SysDistort"+sysDistort
    elif("perugiampiHi" in options.sample):       
        usedSample="TopAnalysis/Configuration/Summer12/TTJets_SemiLeptMGDecays_TuneP11mpiHi_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V19_v1_cff"
        if(eventFilter=='signal only'):
            outputFileName+="SigPerugiampiHi"
        elif(eventFilter=='background only'):
            outputFileName+="BkgPerugiampiHi"
        if(sysDistort!=""):
            if(sysDistort!="data"):
                additionalEventWeights=False # if set to false no variations (SF+-, ...) are done
            outputFileName+="SysDistort"+sysDistort
    elif("perugiaTeV" in options.sample):
        usedSample="TopAnalysis/Configuration/Summer12/TTJets_SemiLeptMGDecays_TuneP11TeV_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V19_v1_cff"
        if(eventFilter=='signal only'):
            outputFileName+="SigPerugiaTeV"
        elif(eventFilter=='background only'):
            outputFileName+="BkgPerugiaTeV"
        if(sysDistort!=""):
            if(sysDistort!="data"):
                additionalEventWeights=False # if set to false no variations (SF+-, ...) are done
            outputFileName+="SysDistort"+sysDistort
    elif("perugia" in options.sample):
        usedSample="TopAnalysis/Configuration/Summer12/TTJets_SemiLeptMGDecays_TuneP11_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V19_v1_cff" 
        if(eventFilter=='signal only'):
            outputFileName+="SigPerugia"
        elif(eventFilter=='background only'):
            outputFileName+="BkgPerugia"
        if(sysDistort!=""):
            if(sysDistort!="data"):
                additionalEventWeights=False # if set to false no variations (SF+-, ...) are done
	    outputFileName+="SysDistort"+sysDistort
    elif("powhegherwig" in options.sample):
        usedSample="TopAnalysis/Configuration/Summer12/TT_CT10_AUET2_8TeV_powheg_herwig_Summer12_DR53X-PU_S10_START53_V19_v1_cff" 
        if(eventFilter=='signal only'):
            outputFileName+="SigPowhegHerwig"
        elif(eventFilter=='background only'):
            outputFileName+="BkgPowhegHerwig"
        if(sysDistort!=""):
            if(sysDistort!="data"):
                additionalEventWeights=False # if set to false no variations (SF+-, ...) are done
	    outputFileName+="SysDistort"+sysDistort
    elif("powheg" in options.sample):
        usedSample="TopAnalysis/Configuration/Summer12/TT_CT10_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v2_cff" 
        if(eventFilter=='signal only'):
            outputFileName+="SigPowheg"
        elif(eventFilter=='background only'):
            outputFileName+="BkgPowheg"
        if(sysDistort!=""):
            if(sysDistort!="data"):
                additionalEventWeights=False # if set to false no variations (SF+-, ...) are done
	    outputFileName+="SysDistort"+sysDistort
    elif("mcatnlo" in options.sample):
        usedSample="TopAnalysis/Configuration/Summer12/TT_8TeV_mcatnlo_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        PythiaSample="True"
        if(eventFilter=='signal only'):
            outputFileName+="SigMcatnlo"
        elif(eventFilter=='background only'):
            outputFileName+="BkgMcatnlo"
        if(sysDistort!=""):
            if(sysDistort!="data"):
                additionalEventWeights=False # if set to false no variations (SF+-, ...) are done
	    outputFileName+="SysDistort"+sysDistort
    elif(options.sample=="wjets"):        
        usedSample="TopAnalysis/Configuration/Summer12/WJetsToLNu_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v2_cff"
	outputFileName+="Wjets"
    elif(options.sample=="wjetsMatchingUp"  ):        
        usedSample="TopAnalysis/Configuration/Summer12/WJetsToLNu_matchingup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
 	additionalEventWeights=False
 	outputFileName+="WjetsMatchUp"
    elif(options.sample=="wjetsMatchingDown"):        
        usedSample="TopAnalysis/Configuration/Summer12/WJetsToLNu_matchingdown_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
 	additionalEventWeights=False
 	outputFileName+="WjetsMatchDown"
    elif(options.sample=="wjetsScaleUp"     ):        
        usedSample="TopAnalysis/Configuration/Summer12/WJetsToLNu_scaleup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v2_cff"
 	additionalEventWeights=False
 	outputFileName+="WjetsScaleUp"
    elif(options.sample=="wjetsScaleDown"   ):
        usedSample="TopAnalysis/Configuration/Summer12/WJetsToLNu_scaledown_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
 	additionalEventWeights=False        
 	outputFileName+="WjetsScaleDown"
    elif(options.sample=="zjets"):        
        usedSample="TopAnalysis/Configuration/Summer12/DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        outputFileName+="Zjets"
    elif(options.sample=="zjetsMatchingUp"):        
        usedSample="TopAnalysis/Configuration/Summer12/DYJetsToLL_M_50_matchingup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
 	additionalEventWeights=False
 	outputFileName+="ZjetsMatchUp"
    elif(options.sample=="zjetsMatchingDown"):        
        usedSample="TopAnalysis/Configuration/Summer12/DYJetsToLL_M_50_matchingdown_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
	additionalEventWeights=False
	outputFileName+="ZjetsMatchDown"
    elif(options.sample=="zjetsScaleUp"):        
        usedSample="TopAnalysis/Configuration/Summer12/DYJetsToLL_M_50_scaleup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
	additionalEventWeights=False
	outputFileName+="ZjetsScaleUp"
    elif(options.sample=="zjetsScaleDown"):        
        usedSample="TopAnalysis/Configuration/Summer12/DYJetsToLL_M_50_scaledown_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
	additionalEventWeights=False      
	outputFileName+="ZjetsScaleDown"
    # INFO: single top scale systematics are exclusive samples
    elif(options.sample=="singleTopS"):        
        usedSample="TopAnalysis/Configuration/Summer12/T_s_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
	outputFileName+="SingleTopS"
    elif(options.sample=="singleTopSScaleDown"):
        usedSample="TopAnalysis/Configuration/Summer12/TToLeptons_s_channel_scaledown_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        additionalEventWeights=False   
        outputFileName+="SingleTopSScaleDown"
    elif(options.sample=="singleTopSScaleUp"):        
        usedSample="TopAnalysis/Configuration/Summer12/TToLeptons_s_channel_scaleup_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        additionalEventWeights=False
        outputFileName+="SingleTopSScaleUp"        
    elif(options.sample=="singleAntiTopS"):        
        usedSample="TopAnalysis/Configuration/Summer12/Tbar_s_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
  	outputFileName+="SingleAntiTopS"        
    elif(options.sample=="singleAntiTopSScaleDown"):        
        usedSample="TopAnalysis/Configuration/Summer12/TBarToLeptons_s_channel_scaledown_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        additionalEventWeights=False
        outputFileName+="SingleAntiTopSScaleDown"
    elif(options.sample=="singleAntiTopSScaleUp"):        
        usedSample="TopAnalysis/Configuration/Summer12/TBarToLeptons_s_channel_scaleup_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        additionalEventWeights=False
	outputFileName+="SingleAntiTopSScaleUp"
    elif(options.sample=="singleTopT"):        
        usedSample="TopAnalysis/Configuration/Summer12/T_t_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
 	outputFileName+="SingleTopT"
    elif(options.sample=="singleTopTScaleDown"):        
        usedSample="TopAnalysis/Configuration/Summer12/TToLeptons_t_channel_scaledown_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
   	additionalEventWeights=False
   	outputFileName+="SingleTopTScaleDown"
    elif(options.sample=="singleTopTScaleUp"):        
        usedSample="TopAnalysis/Configuration/Summer12/TToLeptons_t_channel_scaleup_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"        
   	additionalEventWeights=False
   	outputFileName+="SingleTopTScaleUp"
    elif(options.sample=="singleAntiTopT"):       
        usedSample="TopAnalysis/Configuration/Summer12/Tbar_t_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        outputFileName+="SingleAntiTopT"
    elif(options.sample=="singleAntiTopTScaleDown"):        
        usedSample="TopAnalysis/Configuration/Summer12/TBarToLeptons_t_channel_scaledown_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
   	additionalEventWeights=False
   	outputFileName+="SingleAntiTopTScaleDown"
    elif(options.sample=="singleAntiTopTScaleUp"):        
        usedSample="TopAnalysis/Configuration/Summer12/TBarToLeptons_t_channel_scaleup_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
   	additionalEventWeights=False
   	outputFileName+="SingleAntiTopTScaleUp"        
    elif(options.sample=="singleTopTw"):        
        usedSample="TopAnalysis/Configuration/Summer12/T_tW_channel_DR_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
	outputFileName+="SingleTopTW"
    elif(options.sample=="singleTopTwScaleDown1"):        
        usedSample="TopAnalysis/Configuration/Summer12/TToDilepton_tW_channel_DR_scaledown_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
   	additionalEventWeights=False
   	outputFileName+="SingleTopTWScaleDown1"   
    elif(options.sample=="singleTopTwScaleUp1"):
        usedSample="TopAnalysis/Configuration/Summer12/TToDilepton_tW_channel_DR_scaleup_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
   	additionalEventWeights=False
   	outputFileName+="SingleTopTWScaleUp1"
    elif(options.sample=="singleTopTwScaleDown2"):        
        usedSample="TopAnalysis/Configuration/Summer12/TToThadWlep_tW_channel_DR_scaledown_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
   	additionalEventWeights=False             
   	outputFileName+="SingleTopTWScaleDown2"   
    elif(options.sample=="singleTopTwScaleUp2"):
        usedSample="TopAnalysis/Configuration/Summer12/TToThadWlep_tW_channel_DR_scaleup_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
   	additionalEventWeights=False
   	outputFileName+="SingleTopTWScaleUp2"
    elif(options.sample=="singleTopTwScaleDown3"):        
        usedSample="TopAnalysis/Configuration/Summer12/TToTlepWhad_tW_channel_DR_scaledown_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
   	additionalEventWeights=False             
   	outputFileName+="SingleTopTWScaleDown3"   
    elif(options.sample=="singleTopTwScaleUp3"):
        usedSample="TopAnalysis/Configuration/Summer12/TToTlepWhad_tW_channel_DR_scaleup_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
   	additionalEventWeights=False
   	outputFileName+="SingleTopTWScaleUp3"
    elif(options.sample=="singleAntiTopTw"):        
        usedSample="TopAnalysis/Configuration/Summer12/Tbar_tW_channel_DR_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
	outputFileName+="SingleAntiTopTW"
    elif(options.sample=="singleAntiTopTwScaleDown1"):        
        usedSample="TopAnalysis/Configuration/Summer12/TBarToDilepton_tW_channel_DR_scaledown_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
	additionalEventWeights=False
	outputFileName+="SingleAntiTopTWScaleDown1"
    elif(options.sample=="singleAntiTopTwScaleUp1"):        
        usedSample="TopAnalysis/Configuration/Summer12/TBarToDilepton_tW_channel_DR_scaleup_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
	additionalEventWeights=False
	outputFileName+="SingleAntiTopTWScaleUp1"
    elif(options.sample=="singleAntiTopTwScaleDown2"):        
        usedSample="TopAnalysis/Configuration/Summer12/TBarToThadWlep_tW_channel_DR_scaledown_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
	additionalEventWeights=False
	outputFileName+="SingleAntiTopTWScaleDown2"
    elif(options.sample=="singleAntiTopTwScaleUp2"):        
        usedSample="TopAnalysis/Configuration/Summer12/TBarToThadWlep_tW_channel_DR_scaleup_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
	additionalEventWeights=False
	outputFileName+="SingleAntiTopTWScaleUp2"
    elif(options.sample=="singleAntiTopTwScaleDown3"):        
        usedSample="TopAnalysis/Configuration/Summer12/TBarToTlepWhad_tW_channel_DR_scaledown_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
	additionalEventWeights=False
	outputFileName+="SingleAntiTopTWScaleDown3"
    elif(options.sample=="singleAntiTopTwScaleUp3"):        
        usedSample="TopAnalysis/Configuration/Summer12/TBarToTlepWhad_tW_channel_DR_scaleup_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
	additionalEventWeights=False
	outputFileName+="SingleAntiTopTWScaleUp3"
    elif("zprime_m1000gev_w100gev" in options.sample):
        usedSample="TopAnalysis/Configuration/Summer12/ZPrimeToTTJets_M1000GeV_W100GeV_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_cff"
        additionalEventWeights=False
        if(eventFilter=='signal only'):
            outputFileName+="Sig"
        elif(eventFilter=='background only'):
            outputFileName+="Bkg"
        outputFileName+="ZprimeM1000W100"
    elif(options.sample=="qcd" and decayChannel=='muon'):
        usedSample="TopAnalysis/Configuration/Summer12/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v3_cff"
        PythiaSample="True"
        outputFileName+="QCD"
    elif(options.sample=="qcdEM1" and decayChannel=='electron'):
        usedSample="TopAnalysis/Configuration/Summer12/QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        PythiaSample="True"
        outputFileName+="QCDEM1"
    elif(options.sample=="qcdEM2" and decayChannel=='electron'):
        usedSample="TopAnalysis/Configuration/Summer12/QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        PythiaSample="True"
        outputFileName+="QCDEM2"
    elif(options.sample=="qcdEM3" and decayChannel=='electron'):
        usedSample="TopAnalysis/Configuration/Summer12/QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        PythiaSample="True"
        outputFileName+="QCDEM3"
    elif(options.sample=="qcdBCE1" and decayChannel=='electron'):
        usedSample="TopAnalysis/Configuration/Summer12/QCD_Pt_20_30_BCtoE_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        PythiaSample="True"
        outputFileName+="QCDBCE1"
    elif(options.sample=="qcdBCE2" and decayChannel=='electron'):
        usedSample="TopAnalysis/Configuration/Summer12/QCD_Pt_30_80_BCtoE_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        PythiaSample="True"
        outputFileName+="QCDBCE2"
    elif(options.sample=="qcdBCE3" and decayChannel=='electron'):
        usedSample="TopAnalysis/Configuration/Summer12/QCD_Pt_80_170_BCtoE_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        PythiaSample="True"
        outputFileName+="QCDBCE3"
    elif(options.sample=="WW"):        
        usedSample="TopAnalysis/Configuration/Summer12/WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        PythiaSample="True"
	outputFileName+="WW"
    elif(options.sample=="WZ"):        
        usedSample="TopAnalysis/Configuration/Summer12/WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        PythiaSample="True"
	outputFileName+="WZ"
    elif(options.sample=="ZZ"):        
        usedSample="TopAnalysis/Configuration/Summer12/ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_cff"
        PythiaSample="True"
	outputFileName+="ZZ"
    else:
         usedSample="none"
        
    if(not usedSample=='none'):
        process.load(usedSample)
    else:
        raise NameError, "The sample '"+options.sample+"' is not known, no poolsource file loaded!"


## TEST a specific file, careful- needs to be commented out again afterwards!!!
#process.source.fileNames = cms.untracked.vstring(    
#    ## add your test file here
#    '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/DA4A7A9F-8DE1-E111-82D9-002481E0D480.root'
#    )

#process.source.eventsToProcess = cms.untracked.VEventRange('1:1258499') 
outputFileNamePart=outputFileName
if(options.JEtag!="none"):
    outputFileNamePart+=options.JEtag
if(not "data" in options.sample):
    outputFileName=outputFileNamePart+options.mctag+"PF"
outputFileName+=".root"

#### =================================================
####  Print out summary of configuration parameters
#### =================================================

print "                                     "
print " =================================== "
print "   Basic Configuration Parameters    "
print " =================================== "
print "                                     "
print " Trigger path:                       ","TriggerResults::"+options.triggerTag
print " Jet Type:                           ",jetType 
print " Run PF2PAT:                         ",pfToPAT," ---- Not working if the file does not contain the necessary information"
print " RunningOnData:                      ",runningOnData
print " JsonFile:                           ",jsonFile
print " Label for MC samples:               ",options.mctag
print " Chosen sample                       ",options.sample
print " Lepton decay channel:               ",decayChannel
print " ttbar filter:                       ",eventFilter," ---- Obsolete if RunningOnData != 'MC' "
print " Include tau->l events in SG         ",tau
print " B-tag algo:                         ",bTagAlgo
print " B-tag discriminator cut value:      ",bTagDiscrCut
print " Apply PU reweighting:               ",PUreweigthing
print " Apply Btag reweighting:             ",BtagReweigthing
print " Apply effSF reweighting:            ",effSFReweigthing
print " JEC level in jetKinematics:         ",corrLevel
print " Apply kinematic Fit:                ",applyKinFit," ---- Programme execution slowed down if 'True'"
print " Use double KinFit approach:         ",doubleKinFit
print " mass constraintin (double) KinFit:  ",options.massfix
print " Prob selection for all result plots: ",ProbCutByDefault
print " Write PAT tuples:                   ",writeOutput
print " Analyzed sample:                    ",usedSample+".py"
print " Output file name:                   ",outputFileName
print " Synchronization exercise:           ",cutflowSynch
print " Additional event weights:           ",additionalEventWeights," ---- If 'True' weights are applied to the KinFit analyzers for monitoring, PU, b-tag and lepton eff. variations"
print " Apply MC weights(MC@NLO):           ",MCweighting 
print " Apply shape reweighting, variation: ",sysDistort
print " reduced module content?:            ",reduced
print " rm rec path& use full ttbar sample: ",genFull
print " "
if(removeGenTtbar==True):
    print " All gen level filters using ttbar decay subset are removed" 
    if(runningOnData=="MC"):
        print " Selection for gen plots only via TopAnalysis.TopFilter.sequences.genSelection_cff"
print " "

#### =================================================

## register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(outputFileName))

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source.skipEvents = cms.untracked.uint32(0)
if(not options.eventsToProcess==-42):
    process.maxEvents.input=options.eventsToProcess

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#global tag:
#(see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions)
# 5_3:
#from Configuration.PyReleaseValidation.autoCond import autoCond
if(runningOnData=="MC"):
    process.GlobalTag.globaltag = cms.string('START53_V27::All') # Fall 2012 JEC for full 8 TeV data
else: 
    process.GlobalTag.globaltag = cms.string('FT_53_V21_AN5::All') # Fall 2012 JEC for full 8 TeV data
    
## Needed for redoing the ak5GenJets
#if(runningOnData=="MC" and pfToPAT==False):
#    process.load("TopAnalysis.TopUtils.GenJetParticles_cff")
#    process.load("RecoJets.Configuration.RecoGenJets_cff")  
#    process.p0 = cms.Path(## redo genjets without mu/nu from tau
#                          process.genJetParticles *
#                          process.ak5GenJets
#                          )

## load JSON file for data
if("data" in options.sample):
    jsonFile = 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt'
if(runningOnData=="data"):
    if(jsonFile!=""):
        import FWCore.PythonUtilities.LumiList as LumiList
        import FWCore.ParameterSet.Types as CfgTypes
	print "The following JSON file is used in the cfg file: ", jsonFile
	myLumis = LumiList.LumiList(filename = jsonFile).getCMSSWString().split(',')
	process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
	process.source.lumisToProcess.extend(myLumis)
	## ATTENTION!!! At the moment myLumis are filled in the separate data_cfg files again
	## as otherwise overwritten by load("data_cff")
    else:
        print "No JSON file specified in cfg file (but possibly via CRAB)."
    
    
## ---
##    configure the cutflow scenario
## ---
## std sequence to produce the ttGenEvt
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

## high level trigger filter
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.hltFilter = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::"+options.triggerTag, HLTPaths = ["HLT_IsoMu24_eta2p1_v*"], throw=False)
if(options.mctag=="Summer11"):
    process.hltFilter.HLTPaths=["HLT_IsoMu24_v*"]
if(options.mctag=="Summer12"):
    process.hltFilter.HLTPaths=["HLT_IsoMu24_eta2p1_v*"]
#process.hltFilter = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = ["HLT_Mu15_v*"], throw=False)
#process.hltFilter.HLTPaths = ["HLT_Mu17_TriCentralJet30_v*"]
# check availability of HLT path via "edmProvDump file.root"

## semileptonic selection
process.load("TopAnalysis.TopFilter.sequences.semiLeptonicSelection_cff")
## generator matching
process.load("TopAnalysis.TopFilter.sequences.generatorMatching_cff")
## muon selection
process.load("TopAnalysis.TopFilter.sequences.muonSelection_cff")
## jet selection
process.load("TopAnalysis.TopFilter.sequences.jetSelection_cff")
## redefine veto jets to be sure it is also replaced when running on PF
from TopAnalysis.TopFilter.sequences.jetSelection_cff import goodJets
process.vetoJets.src="goodJets"
process.vetoJets.cut=''

## tool to select muons from gen Particles and save them as new collection
process.load("TopAnalysis.TopUtils.GenCandSelector_cfi")
## use status 3 leptons
process.isolatedGenMuons.target.status = cms.int32( 3)
process.isolatedGenElectrons.target.status = cms.int32( 3)
## generator level based collections and semileptonic selection (muon and jets)
process.load("TopAnalysis.TopFilter.sequences.genSelection_cff")
## generator filter based on direct cuts on lepton and parton pt and eta
## at ttGenEventLevel
from TopAnalysis.TopFilter.filters.SemiLeptonicGenPhaseSpaceFilter_cfi import filterSemiLeptonicGenPhaseSpace
process.filterGenPhaseSpace = filterSemiLeptonicGenPhaseSpace.clone(src = "genEvt")
process.filterLeptonPhaseSpace = filterSemiLeptonicGenPhaseSpace.clone(src = "genEvt")
process.filterLeptonPhaseSpace.partonMaxEta = cms.double(999.0)
process.filterLeptonPhaseSpace.partonMinPt  = cms.double(0.0)

## Generator kinematics selection (https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/1489.html)
## Currently set to 5.0 to prevent bias in throwing out heavy flavours, see:
## https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/1489/4/1/1/1/1/1.html
process.load("GeneratorInterface.GenFilters.TotalKinematicsFilter_cfi")
#process.totalKinematicsFilterDefault = process.totalKinematicsFilter.clone(tolerance = 0.5)
process.totalKinematicsFilterDefault = process.totalKinematicsFilter.clone(tolerance = 5.0)

## ---
## including analysis tools
## ---
## MC weights (needed for MC@NLO)
process.load("TopAnalysis.TopUtils.EventWeightMC_cfi")
## cross section module
process.load("TopAnalysis.TopAnalyzer.MuonCrossSection_cfi")
## jet kinematics analyzer
process.load("TopAnalysis.TopAnalyzer.JetKinematics_cfi")
## muon kinematics analyzer
process.load("TopAnalysis.TopAnalyzer.MuonKinematics_cfi")
## jet quality analyzer
process.load("TopAnalysis.TopAnalyzer.JetQuality_cfi")
## muon quality analyzer
process.load("TopAnalysis.TopAnalyzer.MuonQuality_cfi")
## electron quality analyzer
process.load("TopAnalysis.TopAnalyzer.ElectronQuality_cfi")
## MET analyzer
process.load("TopAnalysis.TopAnalyzer.METKinematics_cfi")
## electron kinematics analyzer
process.load("TopAnalysis.TopAnalyzer.ElectronKinematics_cfi")
## muon jet kinematics analyzer
process.load("TopAnalysis.TopAnalyzer.MuonJetKinematics_cfi")
## muon vertex analyzer
process.load("TopAnalysis.TopAnalyzer.MuonVertexKinematics_cfi")
## mixed object analyzer
process.load("TopAnalysis.TopAnalyzer.MixedObjectsAnalyzer_cfi")
## event listing
process.load("TopAnalysis.TopUtils.EventListing_cfi")
## multi event filter
from RecoMET.METFilters.multiEventFilter_cfi import multiEventFilter

## investigate specific events
#process.multiEventFilterElec = process.multiEventFilter.clone(file = 'TopAnalysis/Configuration/analysis/semiLeptonic/diffXSection/EventsEleSync.txt')
#process.multiEventFilterMuon = process.multiEventFilter.clone(file = 'TopAnalysis/Configuration/analysis/semiLeptonic/diffXSection/EventsMuonSync.txt')
#delattr(process,"multiEventFilter")
#delattr(process,"vetoIncMuons")

## ---
##    set up vertex filter
## ---
#process.PVSelection = cms.EDFilter("PrimaryVertexFilter",
#                                   pvSrc   = cms.InputTag('offlinePrimaryVertices'),
#                                   minNdof = cms.double(4.0),
#                                   maxZ    = cms.double(24.0),
#                                   maxRho  = cms.double(2.0)
#                                   )
## ---
##    set up filter for different ttbar decay channels
## ---
from TopQuarkAnalysis.TopSkimming.TtDecayChannelFilter_cfi import ttDecayChannelFilter
process.ttSemiLeptonicFilter = ttDecayChannelFilter.clone()
process.ttSemiLeptonicFilter.src = "genEvt"
process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon = True
if(tau==True):
    process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.tau  = True
    process.ttSemiLeptonicFilter.restrictTauDecays = cms.PSet(
        leptonic   = cms.bool(True ),
        oneProng   = cms.bool(False),
        threeProng = cms.bool(False)
        )
# take care of electron channel
if(decayChannel=='electron'):
    process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = False
    process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.electron = True
if(not eventFilter=='all'):
    ## adapt output filename
    if(eventFilter=='signal only'):
        print " " # Dummy output to avoid programme crash
    elif(eventFilter=='background only'):
        process.ttSemiLeptonicFilter.invert = True
    elif(eventFilter=='semileptonic electron only'):
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = False
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.electron = True
    elif(eventFilter=='dileptonic electron only'):
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = False
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.electron = True
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchB.electron = True
    elif(eventFilter=='dileptonic muon only'):
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = True
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchB.muon     = True
    elif(eventFilter=='fullhadronic'):
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = False
    elif(eventFilter=='dileptonic muon + electron only'):
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = True
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchB.electron = True
    elif(eventFilter=='dileptonic via tau only'):
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = False
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.tau      = True
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchB.electron = True
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchB.muon     = True
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchB.tau      = True
    elif(eventFilter=='via single tau only'):
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = False
        process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.tau      = True
    else:
        raise NameError, "'"+eventFilter+"' is not a proper eventFilter name choose: 'all', 'signal only', 'background only', 'semileptonic electron only', 'dileptonic electron only', 'dileptonic muon only', 'fullhadronic', 'via single tau only', 'dileptonic via tau only' or 'dileptonic muon + electron only'"
    ## sequence with filter for decay channel and trigger selection hltFilter
    process.filterSequence = cms.Sequence(process.makeGenEvt *
                                          process.ttSemiLeptonicFilter *
                                          process.hltFilter
                                          )
## sequence without filter (only trigger selection hltFilter) - done when 'all' is chosen, removeGenTtbar=True or runningOnData = "data"
if(eventFilter=='all') or (removeGenTtbar==True) or (runningOnData=="data"):
    process.filterSequence = cms.Sequence(process.hltFilter)
     
process.genFilterSequence = cms.Sequence(process.makeGenEvt                      *
                                         process.ttSemiLeptonicFilter             )

## define ordered jets
uds0    = cms.PSet(index = cms.int32(0), correctionLevel = cms.string(corrLevel), useTree = cms.bool(False) )
uds1    = cms.PSet(index = cms.int32(1), correctionLevel = cms.string(corrLevel), useTree = cms.bool(False) )
uds2    = cms.PSet(index = cms.int32(2), correctionLevel = cms.string(corrLevel), useTree = cms.bool(False) )
uds3    = cms.PSet(index = cms.int32(3), correctionLevel = cms.string(corrLevel), useTree = cms.bool(False) )
udsAll  = cms.PSet(index = cms.int32(-1),correctionLevel = cms.string(corrLevel), useTree = cms.bool(False) )

## ---
##    Set up selection steps for muon selection
## ---
process.combinedMuonsSelection        = process.muonSelection.clone (src = 'combinedMuons'       , minNumber = 1, maxNumber = 99999999)
process.kinematicMuonsSelection       = process.muonSelection.clone (src = 'kinematicMuons'      , minNumber = 1, maxNumber = 99999999)
process.trackMuonsSelection           = process.muonSelection.clone (src = 'trackMuons'          , minNumber = 1, maxNumber = 99999999)
process.highPtMuonsSelection          = process.muonSelection.clone (src = 'highPtMuons'         , minNumber = 1, maxNumber = 99999999)
process.goldenMuonsSelection          = process.muonSelection.clone (src = 'goldenMuons'         , minNumber = 1, maxNumber = 99999999)
process.tightMuonsSelection           = process.muonSelection.clone (src = 'tightMuons'          , minNumber = 1, maxNumber = 99999999)

process.muonCuts = cms.Sequence(process.combinedMuonsSelection        +
                                process.highPtMuonsSelection          +
                                process.kinematicMuonsSelection       +
                                process.trackMuonsSelection           +
                                process.goldenMuonsSelection          +
                                process.tightMuonsSelection           +
                                process.muonSelection
                                )

## ---
##    Set up selection steps for different jet multiplicities
## ---
process.leadingJetSelectionNjets1 = process.leadingJetSelection.clone (src = 'tightLeadingPFJets', minNumber = 1)
process.leadingJetSelectionNjets2 = process.leadingJetSelection.clone (src = 'tightLeadingPFJets', minNumber = 2)
process.leadingJetSelectionNjets3 = process.leadingJetSelection.clone (src = 'tightLeadingPFJets', minNumber = 3)
process.leadingJetSelectionNjets4 = process.leadingJetSelection.clone (src = 'tightLeadingPFJets', minNumber = 4)

## ---
##    Set up low level selection steps (sligthly above trigger) for monitoring
## ---
from TopAnalysis.TopFilter.sequences.jetSelection_cff import selectedPatJets
process.looseCentralJets = selectedPatJets.clone(src = 'selectedPatJetsAK5PF', cut = 'abs(eta) < 2.5 & pt > 30.')
process.looseJetSelectionNjets3 = process.leadingJetSelection.clone (src = 'looseCentralJets', minNumber = 3)
process.looseCuts = cms.Sequence(process.looseCentralJets*process.kinematicMuonsSelection+process.looseJetSelectionNjets3)
                               
## ---
##    Set up selection steps for different (gen)-jet multiplicities
## ---
process.leadingGenJetSelectionNjets1 = process.leadingGenJetSelection.clone (src = 'selectedGenJetCollection', minNumber = 1)
process.leadingGenJetSelectionNjets2 = process.leadingGenJetSelection.clone (src = 'selectedGenJetCollection', minNumber = 2)
process.leadingGenJetSelectionNjets3 = process.leadingGenJetSelection.clone (src = 'selectedGenJetCollection', minNumber = 3)
process.leadingGenJetSelectionNjets4 = process.leadingGenJetSelection.clone (src = 'selectedGenJetCollection', minNumber = 4)

#from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
## remove gen jets overlapping with leptons and also potentially radiated photons in the electron channel
process.cleanedGenJetCollection = cms.EDProducer("PATGenJetCleaner",
    src = cms.InputTag("ak5GenJets"),
    ## preselection (any string-based cut on pat::Jet)
    preselection = cms.string(''),
    ## overlap checking configurables
    checkOverlaps = cms.PSet(
       muons = cms.PSet(
       src       = cms.InputTag("selectedGenMuonCollection"),
       algorithm = cms.string("byDeltaR"),
       preselection        = cms.string(''),
       deltaR              = cms.double(0.4),
       checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
       pairCut             = cms.string(""),
       requireNoOverlaps   = cms.bool(True), # overlaps don't cause the jet to be discared
       )
       ),
    ## finalCut (any string-based cut on pat::Jet)
    finalCut = cms.string(''),
    )

process.noOverlapGenJetCollection=cms.EDFilter("CommonGenJetSelector",
                                               src = cms.InputTag("cleanedGenJetCollection"),
                                               cut = cms.string('abs(eta) < 2.4 & pt > 30.')
                                               )

process.leadingCleanedGenJetSelectionNjets4 = process.leadingGenJetSelection.clone (src = 'noOverlapGenJetCollection', minNumber = 4)


process.selectedGenMuonCollection.cut    =cms.string('abs(eta) < 2.1 & pt > 30.')
process.selectedGenElectronCollection.cut=cms.string('abs(eta) < 2.1 & pt > 30.')
process.genAllMuonKinematics = process.analyzeMuonKinematics.clone    (src = 'isolatedGenMuons')
process.genAllElectronKinematics = process.analyzeElectronKinematics.clone(src = 'isolatedGenElectrons')
process.genAllJetKinematics  = process.analyzeJetKinematics.clone(src = 'ak5GenJets', analyze = udsAll)
process.genSelJetKinematics  = process.analyzeJetKinematics.clone(src = 'selectedGenJetCollection', analyze = udsAll)
process.hadLvObjectMonitoring = cms.Sequence(process.genAllElectronKinematics *
                                             process.genAllMuonKinematics     *
                                             process.genAllJetKinematics      *
                                             process.genSelJetKinematics      
                                             )
## ---
##    Set up selection for b-jet multiplicity
## ---

## switch to desired btagging value
process.tightBottomPFJets.src = 'goodJetsPF30'
process.tightBottomPFJets.cut=cms.string('bDiscriminator("'+str(bTagAlgo)+'")>'+str(bTagDiscrCut))
## select number of b tags

process.btagSelection1 = process.bottomJetSelection.clone(src = 'tightBottomPFJets', minNumber = 1, maxNumber = 99999)
process.btagSelection = process.bottomJetSelection.clone(src = 'tightBottomPFJets', minNumber = 2, maxNumber = 99999)
process.btagSelectionSSV=process.btagSelection.clone(src = 'simpleSecondaryVertexHighEffBJets')

## ---
##    Set up selection for light jets (because kinfir needs 2 light jets)
## ---
process.tightLightPFJets  = process.selectedPatJets.clone(src = 'goodJetsPF30', cut=cms.string('bDiscriminator("'+str(bTagAlgo)+'")<'+str(bTagDiscrCut)))
process.lightJetSelection = process.bottomJetSelection.clone(src = 'tightLightPFJets', minNumber = 2, maxNumber = 99999)

## kinematic contributions
## muon
process.tightMuonKinematics        = process.analyzeMuonKinematics.clone (src = 'tightMuons')
process.tightMuonQuality           = process.analyzeMuonQuality.clone    (src = 'tightMuons')

process.tightMuonKinematicsTagged  = process.tightMuonKinematics.clone();
process.tightMuonQualityTagged     = process.tightMuonQuality.clone();
process.tightMuonKinematicsSSV  = process.analyzeMuonKinematics.clone (src = 'tightMuons'    )
process.tightMuonQualitySSV     = process.analyzeMuonQuality.clone    (src = 'tightMuons'    )

process.tightMuonKinematicsNjets1 = process.tightMuonKinematics.clone() 
process.tightMuonQualityNjets1    = process.tightMuonQuality.clone   ()    
process.tightMuonKinematicsNjets2 = process.tightMuonKinematics.clone() 
process.tightMuonQualityNjets2    = process.tightMuonQuality.clone   ()
process.tightMuonKinematicsNjets3 = process.tightMuonKinematics.clone() 
process.tightMuonQualityNjets3    = process.tightMuonQuality.clone   ()

## jets
process.tightLead_0_JetKinematics = process.analyzeJetKinematics.clone(src = 'tightLeadingPFJets', analyze = uds0  )
process.tightLead_1_JetKinematics = process.analyzeJetKinematics.clone(src = 'tightLeadingPFJets', analyze = uds1  )
process.tightLead_2_JetKinematics = process.analyzeJetKinematics.clone(src = 'tightLeadingPFJets', analyze = uds2  )
process.tightLead_3_JetKinematics = process.analyzeJetKinematics.clone(src = 'tightLeadingPFJets', analyze = uds3  )
process.tightJetKinematics        = process.analyzeJetKinematics.clone(src = 'tightLeadingPFJets', analyze = udsAll)
process.tightJetQuality           = process.analyzeJetQuality.clone   (src = 'tightLeadingPFJets')
process.bottomJetKinematics       = process.analyzeJetKinematics.clone(src = 'tightBottomPFJets' , analyze = udsAll)
process.bottomJetQuality          = process.analyzeJetQuality.clone   (src = 'tightBottomPFJets' )
process.tightJetKinematicsSSV     = process.analyzeJetKinematics.clone(src = 'tightLeadingPFJets', analyze = udsAll)
process.tightJetQualitySSV        = process.analyzeJetQuality.clone   (src = 'tightLeadingPFJets')
process.tightlightJetKinematics   = process.analyzeJetKinematics.clone(src = 'tightLightPFJets'  )

process.tightJetKinematicsNjets1=process.tightJetKinematics.clone()
process.tightJetQualityNjets1   =process.tightJetQuality   .clone()
process.tightJetKinematicsNjets2=process.tightJetKinematics.clone()
process.tightJetQualityNjets2   =process.tightJetQuality   .clone()
process.tightJetKinematicsNjets3=process.tightJetKinematics.clone()
process.tightJetQualityNjets3   =process.tightJetQuality   .clone()

process.tightLead_0_JetKinematicsTagged  = process.tightLead_0_JetKinematics.clone()
process.tightLead_1_JetKinematicsTagged  = process.tightLead_1_JetKinematics.clone()
process.tightLead_2_JetKinematicsTagged  = process.tightLead_2_JetKinematics.clone()
process.tightLead_3_JetKinematicsTagged  = process.tightLead_3_JetKinematics.clone()
process.tightJetKinematicsTagged         = process.tightJetKinematics.clone()
process.tightJetKinematicsKinFit         = process.tightJetKinematics.clone()
process.tightJetQualityTagged            = process.tightJetQuality.clone()
process.tightJetQualityKinFit            = process.tightJetQuality.clone()
process.bottomJetKinematicsTagged        = process.bottomJetKinematics.clone()
process.bottomJetKinematicsKinFit        = process.bottomJetKinematics.clone()
process.bottomJetQualityTagged           = process.bottomJetQuality.clone()
process.bottomJetQualityKinFit           = process.bottomJetQuality.clone()
process.bottomLead_0_JetKinematicsTagged = process.analyzeJetKinematics.clone (src = 'tightBottomPFJets', analyze = uds0 )
process.bottomLead_1_JetKinematicsTagged = process.analyzeJetKinematics.clone (src = 'tightBottomPFJets', analyze = uds1 )

## muon&jets
process.tightMuontightJetsKinematics       = process.analyzeMuonJetKinematics.clone(srcA = 'tightMuons', srcB = 'tightLeadingPFJets')
process.tightMuontightJetsKinematicsTagged = process.tightMuontightJetsKinematics.clone()
process.tightMuontightJetsKinematicsSSV = process.analyzeMuonJetKinematics.clone(srcA = 'tightMuons', srcB = 'vetoJets')

## multiple objects
process.compositedKinematics  = process.analyzeCompositedObjects.clone(
                                  JetSrc = 'tightLeadingPFJets',
                                  METSrc = 'patMETs',
                                  MuonSrc = 'tightMuons',
                                  ElectronSrc = 'goodElectronsEJ',
                                  GenJetSrc = cms.InputTag('ak5GenJets'),
                                  addGenJetSrc = cms.InputTag('additionalGenJet'),
                                  GenMETSrc = 'genMetTrue',
                                  GenLepSrc = 'selectedGenMuonCollection',
                                  ingenPS = cms.InputTag("visibleHadronPS", "inVisiblePS"),
                                  weight = "",
                                  VertexSrc = "goodOfflinePrimaryVertices",
                                  semiLepEvent = cms.InputTag(""),
                                  hypoKey = cms.string("kKinFit"),
                                  btagAlgo=cms.string("combinedSecondaryVertexBJetTags"),
                                  btagDiscr=cms.double(0.679),
                                  addJetPt=cms.double(50.)
                                  )
if(decayChannel=='electron'):
    process.compositedKinematics.GenLepSrc = 'selectedGenElectronCollection'
process.compositedKinematicsKinFit = process.compositedKinematics.clone()
process.compositedKinematicsKinFit.semiLepEvent = cms.InputTag("ttSemiLepEvent")
process.compositedKinematicsKinFit.GenJetSrc= cms.InputTag('noOverlapGenJetCollection')
process.compositedKinematicsProbSel=process.compositedKinematicsKinFit.clone()
process.compositedKinematics.btagDiscr=cms.double(0.244) # loose WP for untagged selection
# gen
process.compositedPartonGen          =process.compositedKinematicsKinFit.clone(addGenJetSrc = cms.InputTag('additionalGenJet'))
process.compositedPartonGenPhaseSpace=process.compositedPartonGen.clone()
process.compositedHadronGenPhaseSpace=process.compositedPartonGen.clone()

if(removeGenTtbar==True or eventFilter!="signal only"):
    process.compositedKinematics.GenJetSrc = cms.InputTag('')
    process.compositedKinematics.addGenJetSrc = cms.InputTag('')
    process.compositedKinematics.GenMETSrc = ''
    process.compositedKinematics.GenLepSrc = ''
    process.compositedKinematicsKinFit.GenJetSrc = cms.InputTag('')
    process.compositedKinematicsKinFit.addGenJetSrc = cms.InputTag('')
    process.compositedKinematicsKinFit.GenMETSrc = ''
    process.compositedKinematicsKinFit.GenLepSrc = ''
    process.compositedKinematicsProbSel.GenJetSrc = cms.InputTag('')
    process.compositedKinematicsProbSel.addGenJetSrc = cms.InputTag('')
    process.compositedKinematicsProbSel.GenMETSrc = ''
    process.compositedKinematicsProbSel.GenLepSrc = ''
    process.compositedKinematics.semiLepEvent = cms.InputTag('')
    
## electrons
process.tightElectronKinematics        = process.analyzeElectronKinematics.clone( src = 'goodElectronsEJ'  )
process.tightElectronQuality           = process.analyzeElectronQuality.clone   ( src = 'goodElectronsEJ'  )
process.tightElectronKinematicsTagged  = process.analyzeElectronKinematics.clone( src = 'goodElectronsEJ'  )
process.tightElectronQualityTagged     = process.analyzeElectronQuality.clone   ( src = 'goodElectronsEJ'  )
process.tightElectronKinematicsSSV     = process.analyzeElectronKinematics.clone( src = 'goodElectronsEJ'  )
process.tightElectronQualitySSV        = process.analyzeElectronQuality.clone   ( src = 'goodElectronsEJ'  )

process.tightElectronKinematicsNjets1  = process.tightElectronKinematics.clone()
process.tightElectronQualityNjets1     = process.tightElectronQuality   .clone()
process.tightElectronKinematicsNjets2  = process.tightElectronKinematics.clone()
process.tightElectronQualityNjets2     = process.tightElectronQuality   .clone()
process.tightElectronKinematicsNjets3  = process.tightElectronKinematics.clone()
process.tightElectronQualityNjets3     = process.tightElectronQuality   .clone()

## MET
process.analyzeMETMuon       = process.analyzeMETCorrelations.clone(srcA = 'patMETs', srcB='tightMuons')
process.analyzeMETMuonTagged = process.analyzeMETMuon.clone()

## ---
##    MC Event Weight monitoring
## ---
process.load("TopAnalysis.TopAnalyzer.EventWeightDistributions_cfi")
process.EventWeightDistributions.EventWeightSrc           = cms.InputTag("eventWeightMC")

## ---
##    PU reweighting monitoring
## ---
process.load("TopAnalysis.TopAnalyzer.PUControlDistributions_cfi")
process.PUControlDistributions.PUSource                   = cms.InputTag("addPileupInfo"                       )
process.PUControlDistributions.PVertexSource              = cms.InputTag("goodOfflinePrimaryVertices"          )
process.PUControlDistributions.PUEventWeightSource        = cms.InputTag("eventWeightPUsysNo",  "eventWeightPU"       )
process.PUControlDistributions.PUEventWeightUpSource      = cms.InputTag("eventWeightPUsysUp",  "eventWeightPUUp"     )
process.PUControlDistributions.PUEventWeightDownSource    = cms.InputTag("eventWeightPUsysDown","eventWeightPUDown"   )
process.PUControlDistributions.MCSampleTag                = options.mctag

process.PUControlDistributionsDefault        = process.PUControlDistributions.clone()
process.PUControlDistributionsBeforeBtagging = process.PUControlDistributions.clone()
process.PUControlDistributionsAfterBtagging  = process.PUControlDistributions.clone()

## collect kinematics
process.monitorKinematicsNjets1a = cms.Sequence(process.tightJetKinematicsNjets1  +
                                                process.tightJetQualityNjets1
                                                )

process.monitorKinematicsNjets2a = cms.Sequence(process.tightJetKinematicsNjets2  +
                                                process.tightJetQualityNjets2
                                                )

process.monitorKinematicsNjets3a = cms.Sequence(process.tightJetKinematicsNjets3  +
                                                process.tightJetQualityNjets3
                                                )


if(decayChannel =='muon'):
    process.monitorKinematicsNjets1 = cms.Sequence(process.monitorKinematicsNjets1a  +
                                                   process.tightMuonKinematicsNjets1 +
                                                   process.tightMuonQualityNjets1    
                                                   )
    process.monitorKinematicsNjets2 = cms.Sequence(process.monitorKinematicsNjets2a  +
                                                   process.tightMuonKinematicsNjets2 +
                                                   process.tightMuonQualityNjets2    
                                                   )
    process.monitorKinematicsNjets3 = cms.Sequence(process.monitorKinematicsNjets3a  +
                                                   process.tightMuonKinematicsNjets3 +
                                                   process.tightMuonQualityNjets3    
                                                   )
elif(decayChannel =='electron'):
    process.monitorKinematicsNjets1 = cms.Sequence(process.monitorKinematicsNjets1a      +
                                                   process.tightElectronKinematicsNjets1 +
                                                   process.tightElectronQualityNjets1    
                                                   )
    process.monitorKinematicsNjets2 = cms.Sequence(process.monitorKinematicsNjets2a      +
                                                   process.tightElectronKinematicsNjets2 +
                                                   process.tightElectronQualityNjets2    
                                                   )
    process.monitorKinematicsNjets3 = cms.Sequence(process.monitorKinematicsNjets3a      +
                                                   process.tightElectronKinematicsNjets3 +
                                                   process.tightElectronQualityNjets3    
                                                   )
     
    
process.monitorKinematicsBeforeBtagging = cms.Sequence(process.tightMuonKinematics          +
                                                       process.tightMuonQuality             +
                                                       process.tightLead_0_JetKinematics    +
                                                       process.tightLead_1_JetKinematics    +
                                                       process.tightLead_2_JetKinematics    +
                                                       process.tightLead_3_JetKinematics    +
                                                       process.tightJetKinematics           +
                                                       process.tightJetQuality              +
                                                       process.bottomJetKinematics          +
						       process.bottomJetQuality             +
                                                       process.analyzeMETMuon               +
                                                       process.tightMuontightJetsKinematics 
                                                       #process.compositedKinematics
                                                       )


process.monitorKinematicsAfterBtagging = cms.Sequence(process.tightMuonKinematicsTagged          +
                                                      process.tightMuonQualityTagged             +
                                                      process.tightLead_0_JetKinematicsTagged    +
                                                      process.tightLead_1_JetKinematicsTagged    +
                                                      process.tightLead_2_JetKinematicsTagged    +
                                                      process.tightLead_3_JetKinematicsTagged    +
                                                      process.tightJetKinematicsTagged           +
                                                      process.tightJetQualityTagged              +
                                                      process.bottomJetKinematicsTagged          +
                                                      process.bottomJetQualityTagged             +
                                                      process.analyzeMETMuonTagged               +
                                                      process.tightMuontightJetsKinematicsTagged +
                                                      process.bottomLead_0_JetKinematicsTagged   +
                                                      process.bottomLead_1_JetKinematicsTagged
                                                      )

process.monitorJetKinematicsAfterKinFit = cms.Sequence(process.tightJetKinematicsKinFit           +
                                                       process.tightJetQualityKinFit              +
                                                       process.bottomJetKinematicsKinFit          +
                                                       process.bottomJetQualityKinFit             
                                                       )

process.SSVMonitoring = cms.Sequence(process.tightMuontightJetsKinematicsSSV +
                                     process.tightMuonKinematicsSSV          +
                                     process.tightMuonQualitySSV             +
                                     process.tightJetKinematicsSSV           +
                                     process.tightJetQualitySSV              
                                     )

process.monitorElectronKinematicsBeforeBtagging = cms.Sequence(process.tightElectronKinematics+
                                                               process.tightElectronQuality   )
process.monitorElectronKinematicsAfterBtagging  = cms.Sequence(process.tightElectronKinematicsTagged+
                                                               process.tightElectronQualityTagged   )

# combined jet selection+monitoring
process.jetSelection = cms.Sequence(process.leadingJetSelectionNjets1 +
                                    #process.monitorKinematicsNjets1   +
                                    process.leadingJetSelectionNjets2 +
                                    #process.monitorKinematicsNjets2   +
                                    process.leadingJetSelectionNjets3 +
                                    #process.monitorKinematicsNjets3   +
                                    process.leadingJetSelectionNjets4 
                                    )

## ---
##    configure Kinematic fit
## ---

## produce top reconstructed event
process.load('TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff')
## process.ttSemiLepJetPartonMatch.verbosity = 1

# add hypothesis
from TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff import *
addTtSemiLepHypotheses(process,['kKinFit'])
if(decayChannel=='electron'):
    useElectronsForAllTtSemiLepHypotheses(process,'goodElectronsEJ')
if(eventFilter=='signal only') and (runningOnData=="MC"):
    if(applyKinFit==True):
        print ' kinFit: processing ttbar SG MC - build genmatch'
else:
    removeTtSemiLepHypGenMatch(process)
    if(applyKinFit==True):
        print ' kinFit: processing bkg or data - genmatch removed'

## choose collections
## in fitting procedure
process.kinFitTtSemiLepEventHypothesis.leps = 'tightMuons'
process.kinFitTtSemiLepEventHypothesis.jets = 'tightLeadingPFJets'
process.kinFitTtSemiLepEventHypothesis.mets = 'patMETsPF'
## in genmatch
process.ttSemiLepHypGenMatch.jets = 'tightLeadingPFJets'
process.ttSemiLepHypGenMatch.leps = 'tightMuons'
process.ttSemiLepHypGenMatch.mets = 'patMETsPF'
process.ttSemiLepHypGenMatch.jetCorrectionLevel=corrLevel
## in event hypothesis used for 
process.ttSemiLepHypKinFit.jets = 'tightLeadingPFJets'
process.ttSemiLepHypKinFit.leps = 'tightMuons'
process.ttSemiLepHypKinFit.mets = 'patMETsPF'
process.ttSemiLepEvent.verbosity=0

# maximum number of jets to be considered in the jet combinatorics
# (has to be >= 4, can be set to -1 if you want to take all)
process.kinFitTtSemiLepEventHypothesis.maxNJets = 5

# maximum number of jet combinations finally written into the event, starting from the "best"
# (has to be >= 1, can be set to -1 if you want to take all)
process.kinFitTtSemiLepEventHypothesis.maxNComb = 3

# set constraints:: 1: Whad-mass, 2: Wlep-mass, 3: thad-mass, 4: tlep-mass, 5: nu-mass, 6: equal t-masses, 7: Pt balance
process.kinFitTtSemiLepEventHypothesis.constraints = [1, 2, 6]
#process.kinFitTtSemiLepEventHypothesis.constraints = [1, 2, 3, 4]

# the number of degrees of freedom is the number of constraints reduced by one due to the unmeasured eta of the neutrino
# the cut on the fit probability is not working properly as the inherent ndof is wrong. Therefore, the cut (fitProb>=0.02) is converted to a cut on chi^2.
if(process.kinFitTtSemiLepEventHypothesis.constraints == [1, 2, 3, 4]):
    #process.kinFitTtSemiLepEventHypothesis.mTop = 172.5
    process.kinFitTtSemiLepEventHypothesis.mTop=options.massfix
    degreesOfFreedom = cms.int32(3)
    chi2cut = "fitChi2<=9.8374093111925935418"
    print "ndof set to 3"
elif(process.kinFitTtSemiLepEventHypothesis.constraints == [1, 2, 6]):
    degreesOfFreedom = cms.int32(2)
    chi2cut = "fitChi2<=7.8240460108562923658"
    print "ndof set to 2"
else:
    raise ValueError, "please configure the number of degrees of freedom corresponding to your constraints"

# consider b-tagging in event reconstruction
process.kinFitTtSemiLepEventHypothesis.bTagAlgo = bTagAlgo # "simpleSecondaryVertexHighEffBJetTags"

# TCHE  discr.values 7TeV: 1.70,  3.30  , 10.20
# SSVHE discr.values 7TeV: x.xx,  1.74  ,  x.xx
# CSV   discr.values 7TeV: 0.244, 0.679 , 0.898
process.kinFitTtSemiLepEventHypothesis.minBDiscBJets     = bTagDiscrCut
process.kinFitTtSemiLepEventHypothesis.maxBDiscLightJets = bTagDiscrCut
process.kinFitTtSemiLepEventHypothesis.useBTagging       = True

# eta-dependent scaling of JER (input to kinFit)
if(applyKinFit==True):
    process.kinFitTtSemiLepEventHypothesis.jetEnergyResolutionScaleFactors = cms.vdouble( 1.052 , 1.057 , 1.096 , 1.134 , 1.288   )
    process.kinFitTtSemiLepEventHypothesis.jetEnergyResolutionEtaBinning = cms.vdouble(0.0  ,  0.5  ,  1.1  ,  1.7  ,  2.3  ,  -1.)

## keep only events with unambigues parton matches
## (no other partons exist in dR=0.3 cone) 
## attention: improves purity but reduces efficiency
if(eventFilter=='signal only') and (runningOnData=="MC"):
    process.ttSemiLepJetPartonMatch.algorithm = "unambiguousOnly"
    #process.ttSemiLepJetPartonMatch.algorithm = "totalMinDist"
    process.ttSemiLepJetPartonMatch.useMaxDist = True
    ## set number of jets considered in jet-parton matching
    process.ttSemiLepJetPartonMatch.maxNJets=-1
    ## choose jet collection considered in jet-parton matching
    process.ttSemiLepJetPartonMatch.jets='tightLeadingPFJets'

## ---
##    configure KinFit Analyzers
## ---

## kinfit succeeded?
process.load("TopQuarkAnalysis.TopEventProducers.producers.TtSemiLepEvtFilter_cfi")
process.filterRecoKinFit  = process.ttSemiLepEventFilter.clone( cut = cms.string("isHypoValid('kKinFit')"  ) )
process.filterProbKinFit  = process.ttSemiLepEventFilter.clone( cut = cms.string("isHypoValid('kKinFit') && "+chi2cut ) )
process.filterMatchKinFit = process.ttSemiLepEventFilter.clone( cut = cms.string("isHypoValid('kGenMatch')") )

## add a second kinematic fit which calculates the jet permutation as input for the first kin fit
if(doubleKinFit==True):
    cloneTtSemiLepEvent(process)
    process.kinFitTtSemiLepEventHypothesis2.constraints = [1, 2, 3, 4]
    process.kinFitTtSemiLepEventHypothesis2.mTop = options.massfix
    process.kinFitTtSemiLepEventHypothesis.match = cms.InputTag("kinFitTtSemiLepEventHypothesis2")
    process.kinFitTtSemiLepEventHypothesis.useOnlyMatch = cms.bool(True)
    
    process.filterRecoKinFit2 = process.ttSemiLepEventFilter.clone( src = cms.InputTag("ttSemiLepEvent2"),
                                                                    cut = cms.string("isHypoValid('kKinFit')"  ) )
    
    process.makeTtSemiLepEvent += process.filterRecoKinFit2
    
    process.makeTtSemiLepEvent.replace(process.filterRecoKinFit2,
                                       process.filterRecoKinFit2*process.makeTtSemiLepEventBase)
    
    process.makeTtSemiLepEvent.remove(process.makeTtSemiLepEventBase)

## configure top reconstruction analyzers & define PSets
## A) for top reconstruction analyzer
process.load("TopAnalysis.TopAnalyzer.TopKinematics_cfi")
## 1)  plots built from event hypothesis kinFit after reco selection
recoKinFit        = cms.PSet(hypoKey=cms.string('kKinFit'  ), lepton=cms.string(decayChannel), useTree=cms.bool(True),
                             matchForStabilityAndPurity=cms.bool(False), ttbarInsteadOfLepHadTop = cms.bool(False),
                             maxNJets = process.kinFitTtSemiLepEventHypothesis.maxNJets, ndof = degreesOfFreedom)
process.analyzeTopRecoKinematicsKinFit = process.analyzeTopRecKinematics.clone(analyze=recoKinFit)
## 2)  same as 1) but for top/antitop instead of leptonic/hadronic top
recoKinFitTopAntitop = cms.PSet(hypoKey=cms.string('kKinFit'  ), lepton=cms.string(decayChannel), useTree=cms.bool(True),
                                matchForStabilityAndPurity=cms.bool(False), ttbarInsteadOfLepHadTop = cms.bool(True),
                                maxNJets = process.kinFitTtSemiLepEventHypothesis.maxNJets, ndof = degreesOfFreedom)
process.analyzeTopRecoKinematicsKinFitTopAntitop = process.analyzeTopRecKinematics.clone(analyze=recoKinFitTopAntitop)
## 3)  plots built from event hypothesis of objects from genmatch to partons (ttSemiLepJetPartonMatch) after reco selection
recoGenMatch      = cms.PSet(hypoKey=cms.string('kGenMatch'), lepton=cms.string(decayChannel), useTree=cms.bool(True),
                             matchForStabilityAndPurity=cms.bool(False), ttbarInsteadOfLepHadTop = cms.bool(False),
                             maxNJets = process.kinFitTtSemiLepEventHypothesis.maxNJets, ndof = degreesOfFreedom)
process.analyzeTopRecoKinematicsGenMatch      = process.analyzeTopRecKinematics.clone(analyze=recoGenMatch)
## 4) plots built from parton level objects
## a) after parton level phase space selection
genTtbarSemiMu    = cms.PSet(hypoKey=cms.string("None"     ), lepton=cms.string(decayChannel), useTree=cms.bool(True),
                             matchForStabilityAndPurity=cms.bool(False), ttbarInsteadOfLepHadTop = cms.bool(False),
                             maxNJets = process.kinFitTtSemiLepEventHypothesis.maxNJets, ndof = degreesOfFreedom)
process.analyzeTopPartonLevelKinematics = process.analyzeTopGenKinematics.clone(analyze=genTtbarSemiMu)
## b) without phase space selection
process.analyzeTopPartonLevelKinematicsPhaseSpace = process.analyzeTopGenKinematics.clone(analyze=genTtbarSemiMu)
## c) after hadron level phase space selection
process.analyzeTopHadronLevelKinematicsPhaseSpace = process.analyzeTopGenKinematics.clone(analyze=genTtbarSemiMu)
                                             
## configure Kin Fit performance analyzers
process.load("TopAnalysis.TopAnalyzer.HypothesisKinFit_cfi"    )
hypoKinFit = cms.PSet(hypoKey = cms.string("kKinFit"),
                      lepton  = cms.string(decayChannel),
                      wantTree = cms.bool(True),
                      maxNJets = process.kinFitTtSemiLepEventHypothesis.maxNJets,
                      ndof = degreesOfFreedom)
process.analyzeHypoKinFit = process.analyzeHypothesisKinFit.clone(analyze=hypoKinFit)

## B object reconstruction quality
# 1) Psets:
# a) all
analyzeAll = cms.PSet( corrPerm=cms.bool(False), maxChi2=cms.double(10000) )
# b) correct permutation
analyzeCorrect = cms.PSet( corrPerm=cms.bool(True), maxChi2=cms.double(10000) )
# 2) lepton
process.load("TopAnalysis.TopAnalyzer.HypothesisKinFitLepton_cfi")
process.analyzeHypoKinFitLepton      = process.analyzeHypothesisKinFitLepton.clone(srcA = "ttSemiLepEvent", srcB = "tightMuons", analyze=analyzeAll    )
process.analyzeHypoKinFitLeptonCorr  = process.analyzeHypothesisKinFitLepton.clone(srcA = "ttSemiLepEvent", srcB = "tightMuons", analyze=analyzeCorrect)
# 3) jets
process.load("TopAnalysis.TopAnalyzer.HypothesisKinFitJets_cfi")
process.analyzeHypoKinFitJets      = process.analyzeHypothesisKinFitJets.clone(srcA = "ttSemiLepEvent", srcB = "tightLeadingPFJets", analyze=analyzeAll    )
process.analyzeHypoKinFitJetsCorr  = process.analyzeHypothesisKinFitJets.clone(srcA = "ttSemiLepEvent", srcB = "tightLeadingPFJets", analyze=analyzeCorrect)
# 4) neutrino/MET
process.load("TopAnalysis.TopAnalyzer.HypothesisKinFitMET_cfi" )
process.analyzeHypoKinFitMET      = process.analyzeHypothesisKinFitMET.clone(srcA = "ttSemiLepEvent", srcB = "patMETs", analyze=analyzeAll    )
process.analyzeHypoKinFitMETCorr  = process.analyzeHypothesisKinFitMET.clone(srcA = "ttSemiLepEvent", srcB = "patMETs", analyze=analyzeCorrect)

## ---
##    b hadron level distributions
## ---
## tool to identify b-jets from genjet collection
process.load("TopAnalysis.TopUtils.GenLevelBJetProducer_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi") # supplies PDG ID to real name resolution of MC particles, necessary for GenLevelBJetProducer
process.makeGenLevelBJets=process.produceGenLevelBJets.clone(
    ttGenEvent = cms.InputTag('genEvt'),
    #genJets = cms.InputTag('ak5GenJets','','HLT'),
    genJets = cms.InputTag("cleanedGenJetCollection"),
    deltaR = cms.double(0.5),
    resolveParticleName = cms.bool(False),
    requireTopBquark = cms.bool(True),
    noBBbarResonances = cms.bool(True),
    )
process.altermakeGenLevelBJets=process.produceGenLevelBJets.clone(
    ttGenEvent = cms.InputTag('genEvt'),
    #genJets = cms.InputTag('ak5GenJets','','HLT'),
    genJets = cms.InputTag("noOverlapGenJetCollection"),
    deltaR = cms.double(5.0),
    resolveParticleName = cms.bool(False),
    requireTopBquark = cms.bool(True),
    noBBbarResonances = cms.bool(True),
    )

## tool to select identified bjets from genJet collection
process.load("TopAnalysis.TopUtils.UhhGenJetSelector_cfi")

process.bjetGenJets=process.uhhSelectedGenJets.clone(
    genJet = cms.InputTag("noOverlapGenJetCollection"),
    BHadJetIndex     = cms.InputTag("altermakeGenLevelBJets", "BHadJetIndex"    ),
    AntiBHadJetIndex = cms.InputTag("altermakeGenLevelBJets", "AntiBHadJetIndex"),
    pt =cms.double(30.),
    eta=cms.double(2.4)                          
    )

# make b jet gen and rec plots using the identification from above
process.load("TopAnalysis.TopAnalyzer.SemiLepBjetAnalyzer_cfi")
process.analyzeTopRecoKinematicsBjets=process.analyzeSemiLepBJets.clone(
    semiLepEvent = cms.InputTag("ttSemiLepEvent"),
    hypoKey = cms.string("kKinFit"),
    #genJets = cms.InputTag('ak5GenJets','','HLT'),
    #genJets = cms.InputTag("noOverlapGenJetCollection"),
    genJets = cms.InputTag("bjetGenJets"),
    genLeptons = cms.InputTag('selectedGenMuonCollection'),
    bJetCollection = cms.bool(True),
    recoJets= cms.InputTag('tightLeadingPFJets'),
    useRecBjetsKinematicsBeforeFit= cms.bool(False),
    output = cms.int32(0),
    weight = cms.InputTag(""),
    genPlots = cms.bool(True),
    ingenPS = cms.InputTag("visibleHadronPS", "inVisiblePS"),
    recPlots = cms.bool(True),
    BHadJetIndex     = cms.InputTag("", ""),
    AntiBHadJetIndex = cms.InputTag("", ""),
    #BHadJetIndex     = cms.InputTag("makeGenLevelBJets", "BHadJetIndex"    ),
    #AntiBHadJetIndex = cms.InputTag("makeGenLevelBJets", "AntiBHadJetIndex"),
    #useClosestDrBs= cms.bool(False),
    useClosestDrBs= cms.bool(True),
    useTree = cms.bool(True)
    )

if(decayChannel=="electron"):
    process.analyzeTopRecoKinematicsBjets.genLeptons = cms.InputTag('selectedGenElectronCollection')

process.analyzeTopPartonLevelKinematicsBjets=process.analyzeTopRecoKinematicsBjets.clone(recPlots = cms.bool(False),useRecBjetsKinematicsBeforeFit= cms.bool(False))
process.analyzeTopPartonLevelKinematicsBjetsPhaseSpace=process.analyzeTopRecoKinematicsBjets.clone(recPlots = cms.bool(False),useRecBjetsKinematicsBeforeFit= cms.bool(False))
process.analyzeTopHadronLevelKinematicsBjetsPhaseSpace=process.analyzeTopRecoKinematicsBjets.clone(recPlots = cms.bool(False),useRecBjetsKinematicsBeforeFit= cms.bool(False))
process.analyzeTopHadronLevelKinematicsBjetsPhaseSpace.genJets = cms.InputTag("bjetGenJets")

#process.analyzeTopRecoKinematicsBjets.genJets = cms.InputTag('ak5GenJets','','HLT')
#process.analyzeTopRecoKinematicsBjets.bJetCollection = cms.bool(False)
#process.analyzeTopRecoKinematicsBjets.BHadJetIndex     = cms.InputTag("makeGenLevelBJets", "BHadJetIndex"    )
#process.analyzeTopRecoKinematicsBjets.AntiBHadJetIndex = cms.InputTag("makeGenLevelBJets", "AntiBHadJetIndex")

# b gen jet selection
process.bGenJetSelection    = process.leadingGenJetSelection.clone (src = 'bjetGenJets'   , minNumber = 2)
process.lightGenJetSelection = process.leadingGenJetSelection.clone (src = cms.InputTag("bjetGenJets", "lightJetsFromTop")  , minNumber = 2)

process.additionalGenJet=cms.EDFilter("CommonGenJetSelector",
                                      src = cms.InputTag("bjetGenJets", "lightJetsFromTop"),
                                      cut = cms.string('abs(eta) < 2.4 & pt > 30.')
                                      )

#cms.InputTag(b"jetGenJets", "lightJetsFromTop")


## tool to check whether the event is within the visible hadron level generator PS 
process.load("TopAnalysis.TopUtils.visiblePhaseSpaceFilter_cfi")

process.visibleHadronPS=process.visiblePhaseSpace.clone(
    ## input gen lepton collection     (should be of type reco::GenParticle)
    genLeptons   = cms.InputTag("selectedGenMuonCollection"),
    ## required lepton multiplicity
    nLeptonsMin   = cms.int32(1),
    nLeptonsMax   = cms.int32(1),
    ## input gen jet collection        (should be of type reco::GenJet)
    genJets      = cms.InputTag("noOverlapGenJetCollection"),  
    ## required jet multiplicity
    nJetsMin      = cms.int32(4),
    nJetsMax      = cms.int32(99999),            
    ## input gen b jet collection      (should be of type reco::GenJet)
    genbJets     = cms.InputTag("bjetGenJets"),
    ## required b jet multiplicity
    nbJetsMin     = cms.int32(2),
    nbJetsMax     = cms.int32(2),
    ## input gen light jet collection  (should be of type reco::GenJet)
    genlightJets = cms.InputTag("bjetGenJets", "lightJetsFromTop"),
    ## required light jet multiplicity
    nlightJetsMin = cms.int32(2),
    nlightJetsMax = cms.int32(2)
    )

if(decayChannel=="electron"):
    process.visibleHadronPS.genLeptons = cms.InputTag('selectedGenElectronCollection')
    

## ---
##    lepton hadron level distributions
## ---
process.load("TopAnalysis.TopAnalyzer.SemiLepLeptonAnalyzer_cfi")
process.analyzeTopRecoKinematicsLepton=process.analyzeSemiLepLepton.clone(
                                     semiLepEvent = cms.InputTag("ttSemiLepEvent"),
                                     hypoKey = cms.string("kKinFit"),
                                     recLeptons = cms.InputTag('tightMuons'),
                                     genLeptons = cms.InputTag('selectedGenMuonCollection'),
                                     output = cms.int32(0),
                                     weight = cms.InputTag(""),
                                     genPlots = cms.bool(True),
                                     ingenPS = cms.InputTag("visibleHadronPS", "inVisiblePS"),
                                     recPlots = cms.bool(True),
                                     useRecLeptonKinematicsBeforeFit= cms.bool(False),
                                     useTree = cms.bool(True)
                                     )
if(decayChannel=="electron"):
    process.analyzeTopRecoKinematicsLepton.recLeptons = cms.InputTag('tightElectronsEJ')
    process.analyzeTopRecoKinematicsLepton.genLeptons = cms.InputTag('selectedGenElectronCollection')

process.analyzeTopPartonLevelKinematicsLepton=process.analyzeTopRecoKinematicsLepton.clone(
    genPlots = cms.bool(True), 
    recPlots = cms.bool(False)
    )
process.analyzeTopPartonLevelKinematicsLeptonPhaseSpace=process.analyzeTopRecoKinematicsLepton.clone(
    genPlots = cms.bool(True), 
    recPlots = cms.bool(False)
    )
process.analyzeTopHadronLevelKinematicsLeptonPhaseSpace=process.analyzeTopRecoKinematicsLepton.clone(
    genPlots = cms.bool(True), 
    recPlots = cms.bool(False)
    )

## ---
##    collect hadron level jet selection
## ---   
process.genJetCuts = cms.Sequence(process.leadingGenJetSelectionNjets1 +
                                  process.leadingGenJetSelectionNjets2 +
                                  process.leadingGenJetSelectionNjets3 +
                                  process.leadingGenJetSelectionNjets4 +
                                  process.leadingCleanedGenJetSelectionNjets4       +
                                  # add bjet indices
                                  #process.makeGenLevelBJets                         +
                                  # add bjet selection
                                  process.bGenJetSelection                          +
                                  process.lightGenJetSelection
                                  )

## collect relevant analyzers after probability selection
process.analyzeTopRecoKinematicsKinFitProbSel = process.analyzeTopRecoKinematicsKinFit.clone()
process.analyzeTopRecoKinematicsBjetsProbSel  = process.analyzeTopRecoKinematicsBjets.clone()
process.analyzeTopRecoKinematicsLeptonProbSel = process.analyzeTopRecoKinematicsLepton.clone()

## ---
##    collect KinFit Analyzers depending on sample processed
## ---
## dummy to avoid empty sequences
process.dummy = process.hltFilter.clone()
## check if kinFit should be included in main path
if(applyKinFit==True):
    ## case 1: MC sample
    if(runningOnData=="MC"):
        ## case 1a): ttbar semileptonic mu-signal
        if(eventFilter=='signal only'):
            process.kinFit    = cms.Sequence(process.tightLightPFJets                        +
                                             process.lightJetSelection                       +
                                             process.tightlightJetKinematics                 +
                                             process.makeTtSemiLepEvent                      +
                                             process.filterRecoKinFit                        +
                                             process.analyzeTopRecoKinematicsKinFit          +
                                             process.analyzeTopRecoKinematicsKinFitTopAntitop+
                                             process.analyzeTopRecoKinematicsGenMatch        +
                                             process.analyzeHypoKinFit                       +
                                             process.compositedKinematicsKinFit              +
                                             process.analyzeHypoKinFitLepton                 +
                                             process.analyzeHypoKinFitLeptonCorr             +
                                             process.analyzeHypoKinFitMET                    +
                                             process.analyzeHypoKinFitMETCorr                +
                                             process.analyzeHypoKinFitJets                   +
                                             process.analyzeHypoKinFitJetsCorr               +
                                             process.analyzeTopRecoKinematicsBjets           +
                                             process.analyzeTopRecoKinematicsLepton          +
                                             process.filterProbKinFit                        +
                                             process.analyzeTopRecoKinematicsKinFitProbSel   +
                                             process.analyzeTopRecoKinematicsBjetsProbSel    +
                                             process.analyzeTopRecoKinematicsLeptonProbSel   +
                                             process.compositedKinematicsProbSel             +
                                             process.filterMatchKinFit
                                             )
            if(ProbCutByDefault):
                process.kinFit.replace(process.filterRecoKinFit,
                                       process.filterRecoKinFit+process.filterProbKinFit)
                
            process.kinFitGen           = cms.Sequence(process.analyzeTopPartonLevelKinematics     
                                                       # add bjet indices
                                                       #process.makeGenLevelBJets                   +
                                                       # add bjet analyzer
                                                       #process.analyzeTopPartonLevelKinematicsBjets +
                                                       # add lepton analyzer
                                                       #process.analyzeTopPartonLevelKinematicsLepton +
                                                       # add mixed object analyzer
                                                       #process.compositedPartonGen
                                                       )
            process.kinFitGenPhaseSpace = cms.Sequence(process.analyzeTopPartonLevelKinematicsPhaseSpace
                                                       # add bjet indices
                                                       #process.makeGenLevelBJets                        +
                                                       # add bjet analyzer
                                                       #process.analyzeTopPartonLevelKinematicsBjetsPhaseSpace +
                                                       # add lepton analyzer
                                                       #process.analyzeTopPartonLevelKinematicsLeptonPhaseSpace +
                                                       # add mixed object analyzer
                                                       #process.compositedPartonGenPhaseSpace
                                                       )
            process.kinFitGenPhaseSpaceHad = cms.Sequence(
                                                          # default analyzer module
                                                          #process.analyzeTopHadronLevelKinematicsPhaseSpace +
                                                          # add bjet analyzer
                                                          process.analyzeTopHadronLevelKinematicsBjetsPhaseSpace +
                                                          # add lepton analyzer
                                                          process.analyzeTopHadronLevelKinematicsLeptonPhaseSpace +
                                                          # add mixed object analyzer
                                                          process.compositedHadronGenPhaseSpace
                                                          )

        ## case 1b): other MC
        else:
            process.kinFit    = cms.Sequence(process.tightLightPFJets                        +
                                             process.lightJetSelection                       +
                                             process.tightlightJetKinematics                 +
                                             process.makeTtSemiLepEvent                      +
                                             process.filterRecoKinFit                        +
                                             process.analyzeTopRecoKinematicsKinFit          +
                                             process.analyzeTopRecoKinematicsKinFitTopAntitop+
                                             process.compositedKinematicsKinFit              +
                                             process.filterProbKinFit                        +
                                             process.analyzeTopRecoKinematicsKinFitProbSel   +
                                             process.compositedKinematicsProbSel             
                                             )
            if(ProbCutByDefault):
                process.kinFit.replace(process.filterRecoKinFit,
                                       process.filterRecoKinFit+process.filterProbKinFit)
            process.kinFitGen           = cms.Sequence(process.dummy)
            process.kinFitGenPhaseSpace = cms.Sequence(process.dummy)
            process.kinFitGenPhaseSpaceHad = cms.Sequence(process.dummy)
    ## case 2: data sample
    elif(runningOnData=="data"):
        process.kinFit    = cms.Sequence(process.tightLightPFJets                        +
                                         process.lightJetSelection                       +
                                         process.tightlightJetKinematics                 +
                                         process.makeTtSemiLepEvent                      +
                                         process.filterRecoKinFit                        +
                                         process.analyzeTopRecoKinematicsKinFit          +
                                         process.analyzeTopRecoKinematicsKinFitTopAntitop+
                                         process.compositedKinematicsKinFit              +
                                         process.filterProbKinFit                        +
                                         process.analyzeTopRecoKinematicsKinFitProbSel   +
                                         process.compositedKinematicsProbSel             
                                         )
        if(ProbCutByDefault):
            process.kinFit.replace(process.filterRecoKinFit,
                                   process.filterRecoKinFit+process.filterProbKinFit)
        process.kinFitGen           = cms.Sequence(process.dummy)
        process.kinFitGenPhaseSpace = cms.Sequence(process.dummy)
        process.kinFitGenPhaseSpaceHad = cms.Sequence(process.dummy)
    else:
         print "choose runningOnData= data or MC"
else:
    process.kinFit              = cms.Sequence(process.dummy)
    process.kinFitGen           = cms.Sequence(process.dummy)
    process.kinFitGenPhaseSpace = cms.Sequence(process.dummy)
    process.kinFitGenPhaseSpaceHad = cms.Sequence(process.dummy)

## ============================
##  MC PU reweighting
## ============================

process.load("TopAnalysis.TopUtils.EventWeightPU_cfi")

## Apply common setting before module is cloned for systematic studies

process.eventWeightPU.MCSampleTag = options.mctag

if (options.mctag == "Fall11"):
    process.eventWeightPU.MCSampleHistoName        = "histo_Fall11_true"
    process.eventWeightPU.DataHistoName            = "histoData_true"
elif (options.mctag == "Summer11"):    
    process.eventWeightPU.MCSampleHistoName        = "histoSummer11_flat_true"
    process.eventWeightPU.DataHistoName            = "histoData_true_fineBinning"
elif (options.mctag == "Summer12"):
    process.eventWeightPU.MCSampleHistoName        = "puhisto"
    process.eventWeightPU.DataHistoName            = "pileup"

    process.eventWeightPU.MCSampleFile             = "TopAnalysis/TopUtils/data/MC_PUDist_Summer12_S10.root"
    process.eventWeightPU.DataFile                 = "TopAnalysis/TopUtils/data/Data_PUDist_sysNo_69400_2012ABCD22JanReReco_190456-208686_8TeV.root"

process.eventWeightPUsysNo   = process.eventWeightPU.clone()
process.eventWeightPUsysUp   = process.eventWeightPU.clone()
process.eventWeightPUsysDown = process.eventWeightPU.clone()

#### Parameters 'CreateWeight3DHisto' and 'Weight3DHistoFile' required for cff-file, but actually not used for Fall11 samples
    
#### Configuration for Nominal PU Weights

process.eventWeightPUsysNo.WeightName          = "eventWeightPU"
if ("11" in options.mctag):
    process.eventWeightPUsysNo.DataFile        = "TopAnalysis/TopUtils/data/Data_PUDist_2011Full.root"
elif ("12" in options.mctag):
    process.eventWeightPUsysNo.DataFile        = "TopAnalysis/TopUtils/data/Data_PUDist_sysNo_69400_2012ABCD22JanReReco_190456-208686_8TeV.root"
else:
    print "only configured for 2011 and 2012 so far!"
process.eventWeightPUsysNo.CreateWeight3DHisto = False 
process.eventWeightPUsysNo.Weight3DHistoFile   = "TopAnalysis/TopUtils/data/DefaultWeight3D.root"

#### Configuration for PU Up Variations

process.eventWeightPUsysUp.WeightName          = "eventWeightPUUp"
if ("11" in options.mctag):
    process.eventWeightPUsysUp.DataFile        = "TopAnalysis/TopUtils/data/Data_PUDist_sysUp_2011Full.root"
elif ("12" in options.mctag):
    process.eventWeightPUsysUp.DataFile        = "TopAnalysis/TopUtils/data/Data_PUDist_sysUp_72870_2012ABCD22JanReReco_190456-208686_8TeV.root"
else:
    print "only configured for 2011 and 2012 so far!"
process.eventWeightPUsysUp.CreateWeight3DHisto = False
process.eventWeightPUsysUp.Weight3DHistoFile   = "TopAnalysis/TopUtils/data/DefaultWeight3DUp.root"

#### Configuration for PU Down Variations

process.eventWeightPUsysDown.WeightName          = "eventWeightPUDown"
if ("11" in options.mctag):
    process.eventWeightPUsysDown.DataFile        = "TopAnalysis/TopUtils/data/Data_PUDist_sysDown_2011Full.root"
elif ("12" in options.mctag):
    process.eventWeightPUsysDown.DataFile        = "TopAnalysis/TopUtils/data/Data_PUDist_sysDn_65930_2012ABCD22JanReReco_190456-208686_8TeV.root"
else:
    print "only configured for 2011 and 2012 so far!"
process.eventWeightPUsysDown.CreateWeight3DHisto = False
process.eventWeightPUsysDown.Weight3DHistoFile   = "TopAnalysis/TopUtils/data/DefaultWeight3DDown.root"

#### event weight sequence

process.makeEventWeightsPU = cms.Sequence(process.eventWeightPUsysNo   *
                                          process.eventWeightPUsysUp   *
                                          process.eventWeightPUsysDown  )

#### relevant PU event weights (potentially merged with shape distortion weights)

PUweightraw     = cms.InputTag("eventWeightPUsysNo",  "eventWeightPU")
PUweightrawUp   = cms.InputTag("eventWeightPUsysUp",  "eventWeightPUUp")
PUweightrawDown = cms.InputTag("eventWeightPUsysDown","eventWeightPUDown")

## ---
##    combine MC event weights with PU event weights
## ---
## tool to multiply event weights
process.load("TopAnalysis.TopUtils.EventWeightMultiplier_cfi")

weightlistMCandPU     = cms.VInputTag()
weightlistMCandPUUp   = cms.VInputTag()
weightlistMCandPUDown = cms.VInputTag()

weightlistMCandPU     .append("eventWeightMC")
weightlistMCandPUUp   .append("eventWeightMC")
weightlistMCandPUDown .append("eventWeightMC")

weightlistMCandPU     .append(PUweightraw)
weightlistMCandPUUp   .append(PUweightrawUp)
weightlistMCandPUDown .append(PUweightrawDown)

process.eventWeightMCandPUweight     = process.eventWeightMultiplier.clone(eventWeightTags = weightlistMCandPU    )
process.eventWeightMCandPUweightUp   = process.eventWeightMultiplier.clone(eventWeightTags = weightlistMCandPUUp  )
process.eventWeightMCandPUweightDown = process.eventWeightMultiplier.clone(eventWeightTags = weightlistMCandPUDown)

process.combineEventWeightsMCandPU = cms.Sequence( process.eventWeightMCandPUweight    *
                                                   process.eventWeightMCandPUweightUp  *
                                                   process.eventWeightMCandPUweightDown )

## =================================================
##    MC ttbar systematic variation reweighting
## =================================================

## load weight producer
process.load("TopAnalysis.TopUtils.EventWeightDileptonModelVariation_cfi")
# specify parameters
#process.eventWeightDileptonModelVariation.variation=cms.string(sysDistort)
#process.eventWeightDileptonModelVariation.ttGenEvent=cms.InputTag('genEvt')
process.eventWeightDileptonModelVariation.ttGenEvent = cms.InputTag('genEvt')
process.eventWeightDileptonModelVariation.weightVariable = cms.string('ttbarmass') #valid values: toppt, topeta, ttbarmass, ttbarmasslandau, data
process.eventWeightDileptonModelVariation.minWeight = cms.double(0.1) #low cut-off, at least 0.1 event weight
process.eventWeightDileptonModelVariation.maxWeight = cms.double(2.)  #high cut-off, at most 2 event weight
process.eventWeightDileptonModelVariation.landauMPV = cms.double(420)
process.eventWeightDileptonModelVariation.landauSigma = cms.double(34)

if(sysDistort=='data'):
    process.eventWeightDileptonModelVariation.weightVariable = cms.string('data')
    process.eventWeightDileptonModelVariation.weight1x = cms.double(112.766) 
    process.eventWeightDileptonModelVariation.slope    = cms.double(0.00141)
elif(sysDistort=='topptTTsys'):
    process.eventWeightDileptonModelVariation.weightVariable = cms.string('topptTTsys')
    process.eventWeightDileptonModelVariation.weight1x = cms.double(102.4) 
    process.eventWeightDileptonModelVariation.slope    = cms.double(0.00125)
elif( sysDistort.find('ttbarMass')>-1):
    process.eventWeightDileptonModelVariation.weightVariable = cms.string('ttbarmass')
    process.eventWeightDileptonModelVariation.weight1x = cms.double(350.) 
    if( sysDistort.find('Up'  )>-1):
        process.eventWeightDileptonModelVariation.slope = cms.double( 0.0008)
    elif( sysDistort.find('Down')>-1):
        process.eventWeightDileptonModelVariation.slope = cms.double(-0.0008)
elif( sysDistort.find('topPt')>-1):
    process.eventWeightDileptonModelVariation.weight1x = cms.double(0.318)    # equals 2*a from exp(a+bx) fit to 8TeV l+jets 12/fb data/MadGraph+Pythia
    process.eventWeightDileptonModelVariation.slope = cms.double(-0.00282)    # equals 2*b from exp(a+bx) fit to 8TeV l+jets 12/fb data/MadGraph+Pythia
    if( sysDistort.find('Up')>-1):
        process.eventWeightDileptonModelVariation.weightVariable = cms.string('topptHard')
    elif( sysDistort.find('Down')>-1):
        process.eventWeightDileptonModelVariation.weightVariable = cms.string('topptSoft')
        
# multiply with PU weight
eventWeightDileptonModelVariation=cms.InputTag("eventWeightDileptonModelVariation")
weightlistDistortPU=cms.VInputTag()
weightlistDistortPU.append(PUweightraw)
weightlistDistortPU.append(eventWeightDileptonModelVariation)
weightlistDistortPUup=cms.VInputTag()
weightlistDistortPUup.append(PUweightrawUp)
weightlistDistortPUup.append(eventWeightDileptonModelVariation)
weightlistDistortPUdown=cms.VInputTag()
weightlistDistortPUdown.append(PUweightrawDown)
weightlistDistortPUdown.append(eventWeightDileptonModelVariation)
process.eventWeightPUDistort     = process.eventWeightMultiplier.clone(eventWeightTags = weightlistDistortPU)
process.eventWeightPUupDistort   = process.eventWeightMultiplier.clone(eventWeightTags = weightlistDistortPUup)
process.eventWeightPUdownDistort = process.eventWeightMultiplier.clone(eventWeightTags = weightlistDistortPUdown)
if(sysDistort==''):
    # final PU weights without shape distortion weights
    PUweight    =cms.InputTag("eventWeightPUsysNo",  "eventWeightPU")
    PUweightUp  =cms.InputTag("eventWeightPUsysUp",  "eventWeightPUUp")
    PUweightDown=cms.InputTag("eventWeightPUsysDown","eventWeightPUDown")
else:
    # final PU x shape distortion weights
    PUweight    =cms.InputTag("eventWeightPUDistort")
    PUweightUp  =cms.InputTag("eventWeightPUupDistort")
    PUweightDown=cms.InputTag("eventWeightPUdownDistort")

## ---
##    MC B-tag reweighting
## ---
## load BTV database
process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")
process.load ("RecoBTag.PerformanceDB.BTagPerformanceDB1107")

process.load("TopAnalysis.TopUtils.BTagSFEventWeight_cfi")
## DB only accepts short tagger algo names
bTagAlgoShort = "CSVM"
if(bTagAlgo =='simpleSecondaryVertexHighEffBJetTags'):
    bTagAlgoShort = "SSVHEM"
process.bTagSFEventWeight.jets=cms.InputTag("tightLeadingPFJets")
process.bTagSFEventWeight.bTagAlgo=bTagAlgoShort
process.bTagSFEventWeight.version="2012"
process.bTagSFEventWeight.sysVar   = cms.string("") # bTagSFUp, bTagSFDown, misTagSFUp, misTagSFDown, 
                                                    # bTagSFShapeUpPt, bTagSFShapeDownPt, bTagSFShapeUpEta, bTagSFShapeDownEta possible;
#process.bTagSFEventWeight.filename= "TopAnalysis/Configuration/data/analyzeBTagEfficiency"+bTagAlgoShort+".root"
#FIXME: use the following file for the next round of resubmitting:
process.bTagSFEventWeight.filename= "TopAnalysis/Configuration/data/analyzeBTagEfficiency2012.root" 
process.bTagSFEventWeight.verbose=cms.int32(0)

## for SSV:
process.bTagSFEventWeightSSV                      =  process.bTagSFEventWeight.clone(bTagAlgo = "SSVHEM", 
                                                                       filename = "TopAnalysis/Configuration/data/analyzeBTagEfficiencySSVHEM.root")

process.bTagSFEventWeightBTagSFUp                  = process.bTagSFEventWeight.clone(sysVar = "bTagSFUp")
process.bTagSFEventWeightBTagSFDown                = process.bTagSFEventWeight.clone(sysVar = "bTagSFDown")
process.bTagSFEventWeightMisTagSFUp                = process.bTagSFEventWeight.clone(sysVar = "misTagSFUp")
process.bTagSFEventWeightMisTagSFDown              = process.bTagSFEventWeight.clone(sysVar = "misTagSFDown")
## per default SF is only varied by half of the BTV uncertainty for shape variations in BTagSFEventWeight_cfi 
## -> "HalfShape" not needed anymore (commented out later on)
process.bTagSFEventWeightBTagSFShapeUpPt65         = process.bTagSFEventWeight.clone(sysVar = "bTagSFShapeUpPt")
process.bTagSFEventWeightBTagSFShapeDownPt65       = process.bTagSFEventWeight.clone(sysVar = "bTagSFShapeDownPt")
process.bTagSFEventWeightBTagSFShapeUpEta0p7       = process.bTagSFEventWeight.clone(sysVar = "bTagSFShapeUpEta")
process.bTagSFEventWeightBTagSFShapeDownEta0p7     = process.bTagSFEventWeight.clone(sysVar = "bTagSFShapeDownEta")
process.bTagSFEventWeightBTagSFShapeUpPt100        = process.bTagSFEventWeight.clone(sysVar = "bTagSFShapeUpPt",    shapeVarPtThreshold   = 100.0)
process.bTagSFEventWeightBTagSFShapeDownPt100      = process.bTagSFEventWeight.clone(sysVar = "bTagSFShapeDownPt",  shapeVarPtThreshold   = 100.0)
process.bTagSFEventWeightBTagSFShapeUpEta1p2       = process.bTagSFEventWeight.clone(sysVar = "bTagSFShapeUpEta",   shapeVarEtaThreshold  =   1.2)
process.bTagSFEventWeightBTagSFShapeDownEta1p2     = process.bTagSFEventWeight.clone(sysVar = "bTagSFShapeDownEta", shapeVarEtaThreshold  =   1.2)
process.bTagSFEventWeightBTagSFFullShapeUpPt65     = process.bTagSFEventWeight.clone(sysVar = "bTagSFShapeUpPt",    shapeDistortionFactor =   1.)
process.bTagSFEventWeightBTagSFFullShapeDownPt65   = process.bTagSFEventWeight.clone(sysVar = "bTagSFShapeDownPt",  shapeDistortionFactor =   1.)
process.bTagSFEventWeightBTagSFFullShapeUpEta0p7   = process.bTagSFEventWeight.clone(sysVar = "bTagSFShapeUpEta",   shapeDistortionFactor =   1.)
process.bTagSFEventWeightBTagSFFullShapeDownEta0p7 = process.bTagSFEventWeight.clone(sysVar = "bTagSFShapeDownEta", shapeDistortionFactor =   1.)

## ---
##    MC eff SF reweighting
## ---
## scale factor for trigger and lepton selection efficiency
process.load("TopAnalysis.TopUtils.EffSFLepton2DEventWeight_cfi")
process.effSFMuonEventWeight=process.effSFLepton2DEventWeight.clone()
process.effSFMuonEventWeight.particles=cms.InputTag("tightMuons")
process.effSFMuonEventWeight.sysVar   = cms.string("")
process.effSFMuonEventWeight.filename= "TopAnalysis/Configuration/data/MuonEffSF2D2012.root"
process.effSFMuonEventWeight.verbose=cms.int32(0)
process.effSFMuonEventWeight.additionalFactor=1.0 ## lepton selection and trigger eff. SF both included in loaded histo
process.effSFMuonEventWeight.additionalFactorErr=0.01 ## 1.0% sys error to account for selection difference Z - ttbar
process.effSFMuonEventWeight.shapeDistortionFactor=-1

process.effSFMuonEventWeightPUup              = process.effSFMuonEventWeight.clone(sysVar = "PUup")
process.effSFMuonEventWeightPUdown            = process.effSFMuonEventWeight.clone(sysVar = "PUdown")
process.effSFMuonEventWeightFlatEffSF         = process.effSFMuonEventWeight.clone(sysVar = "FlatEffSF")
process.effSFMuonEventWeightEffSFNormUpStat   = process.effSFMuonEventWeight.clone(sysVar = "combinedEffSFNormUpStat")
process.effSFMuonEventWeightEffSFNormDownStat = process.effSFMuonEventWeight.clone(sysVar = "combinedEffSFNormDownStat")
process.effSFMuonEventWeightEffSFShapeUpEta   = process.effSFMuonEventWeight.clone(sysVar = "combinedEffSFShapeUpEta")
process.effSFMuonEventWeightEffSFShapeDownEta = process.effSFMuonEventWeight.clone(sysVar = "combinedEffSFShapeDownEta")
process.effSFMuonEventWeightEffSFShapeUpPt    = process.effSFMuonEventWeight.clone(sysVar = "combinedEffSFShapeUpPt")
process.effSFMuonEventWeightEffSFShapeDownPt  = process.effSFMuonEventWeight.clone(sysVar = "combinedEffSFShapeDownPt")
process.effSFMuonEventWeightEffSFShapeUpPt40  = process.effSFMuonEventWeight.clone(sysVar = "combinedEffSFShapeUpPt",   shapeVarPtThreshold=40.0)
process.effSFMuonEventWeightEffSFShapeDownPt40= process.effSFMuonEventWeight.clone(sysVar = "combinedEffSFShapeDownPt", shapeVarPtThreshold=40.0)
process.effSFMuonEventWeightEffSFNormUpSys    = process.effSFMuonEventWeight.clone(sysVar = "combinedEffSFNormUpSys")
process.effSFMuonEventWeightEffSFNormDownSys  = process.effSFMuonEventWeight.clone(sysVar = "combinedEffSFNormDownSys")

process.effSFElectronEventWeight=process.effSFLepton2DEventWeight.clone()
process.effSFElectronEventWeight.particles=cms.InputTag("goodElectronsEJ")
process.effSFElectronEventWeight.jets=cms.InputTag("tightLeadingPFJets")
process.effSFElectronEventWeight.sysVar   = cms.string("")
process.effSFElectronEventWeight.filename= "TopAnalysis/Configuration/data/EleEffSF2D2012.root"
process.effSFElectronEventWeight.verbose=cms.int32(0)
process.effSFElectronEventWeight.additionalFactor=1.0 ## lepton selection and trigger eff. SF both included in loaded histo
process.effSFElectronEventWeight.additionalFactorErr=0.01 ## 1.0% sys error to account for selection difference Z - ttbar
process.effSFElectronEventWeight.shapeDistortionFactor=-1

process.effSFElectronEventWeightPUup              = process.effSFElectronEventWeight.clone(sysVar = "PUup"  ) #load PU variation up results as SF histo from rootfile
process.effSFElectronEventWeightPUdown            = process.effSFElectronEventWeight.clone(sysVar = "PUdown") #load PU variation dn results as SF histo from rootfile
process.effSFElectronEventWeightFlatEffSF         = process.effSFElectronEventWeight.clone(sysVar = "FlatEffSF") # use flat SF as specified in additionalFactor
process.effSFElectronEventWeightEffSFNormUpStat   = process.effSFElectronEventWeight.clone(sysVar = "combinedEffSFNormUpStat"  ) # use SF+stat error up
process.effSFElectronEventWeightEffSFNormDownStat = process.effSFElectronEventWeight.clone(sysVar = "combinedEffSFNormDownStat") # use SF-stat error dn 
process.effSFElectronEventWeightEffSFShapeUpEta   = process.effSFElectronEventWeight.clone(sysVar = "combinedEffSFShapeUpEta"  ) # use SF+/-stat error up/dn for |eta|>/< shapeVarEtaThreshold
process.effSFElectronEventWeightEffSFShapeDownEta = process.effSFElectronEventWeight.clone(sysVar = "combinedEffSFShapeDownEta") # use SF-/+stat error dn/up for |eta|>/< shapeVarEtaThreshold
process.effSFElectronEventWeightEffSFShapeUpPt    = process.effSFElectronEventWeight.clone(sysVar = "combinedEffSFShapeUpPt"   ) # use SF+/-stat error up/dn for pT>/< shapeVarEtaThreshold
process.effSFElectronEventWeightEffSFShapeDownPt  = process.effSFElectronEventWeight.clone(sysVar = "combinedEffSFShapeDownPt" ) # use SF-/+stat error dn/up for pT>/< shapeVarEtaThreshold
process.effSFElectronEventWeightEffSFShapeUpPt40  = process.effSFElectronEventWeight.clone(sysVar = "combinedEffSFShapeUpPt"  , shapeVarPtThreshold=40.) # see above
process.effSFElectronEventWeightEffSFShapeDownPt40= process.effSFElectronEventWeight.clone(sysVar = "combinedEffSFShapeDownPt", shapeVarPtThreshold=40.) # see above
process.effSFElectronEventWeightEffSFNormUpSys    = process.effSFElectronEventWeight.clone(sysVar = "combinedEffSFNormUpSys"   ) # use SF+additionalFactorErr
process.effSFElectronEventWeightEffSFNormDownSys  = process.effSFElectronEventWeight.clone(sysVar = "combinedEffSFNormDownSys" ) # use SF-additionalFactorErr

## ---
##    collect all eventweights
## ---

## create weightlists
weightlistFinal                    =cms.VInputTag()
weightlistNoEffSFWeight            =cms.VInputTag()
weightlistNoBtagSFWeight           =cms.VInputTag()
weightlistNoPUWeight               =cms.VInputTag()
weightlistPUup                     =cms.VInputTag()
weightlistPUdown                   =cms.VInputTag()
weightlistFlatEffSF                =cms.VInputTag()
weightlistEffSFNormUpStat          =cms.VInputTag()
weightlistEffSFNormDownStat        =cms.VInputTag()
weightlistEffSFShapeUpEta          =cms.VInputTag()
weightlistEffSFShapeDownEta        =cms.VInputTag()
weightlistEffSFShapeUpPt           =cms.VInputTag()
weightlistEffSFShapeDownPt         =cms.VInputTag()
weightlistEffSFShapeUpPt40         =cms.VInputTag()
weightlistEffSFShapeDownPt40       =cms.VInputTag()
weightlistEffSFNormUpSys           =cms.VInputTag()
weightlistEffSFNormDownSys         =cms.VInputTag()
weightlistFinalSSV                 =cms.VInputTag()
weightlistBtagSFup                 =cms.VInputTag()
weightlistBtagSFdown               =cms.VInputTag()
weightlistMisTagSFup               =cms.VInputTag()
weightlistMisTagSFdown             =cms.VInputTag()
weightlistBTagSFShapeUpPt65        =cms.VInputTag()
weightlistBTagSFShapeDownPt65      =cms.VInputTag()
weightlistBTagSFShapeUpEta1p2      =cms.VInputTag()
weightlistBTagSFShapeDownEta1p2    =cms.VInputTag()
weightlistBTagSFFullShapeUpPt65    =cms.VInputTag()
weightlistBTagSFFullShapeDownPt65  =cms.VInputTag()
weightlistBTagSFFullShapeUpEta0p7  =cms.VInputTag()
weightlistBTagSFFullShapeDownEta0p7=cms.VInputTag()
weightlistBTagSFShapeUpPt100       =cms.VInputTag()
weightlistBTagSFShapeDownPt100     =cms.VInputTag()
weightlistBTagSFShapeUpEta0p7      =cms.VInputTag()
weightlistBTagSFShapeDownEta0p7    =cms.VInputTag()

if (MCweighting):
    weightlistFinal                    .append("eventWeightMC")
    weightlistNoEffSFWeight            .append("eventWeightMC")
    weightlistNoBtagSFWeight           .append("eventWeightMC")
    weightlistNoPUWeight               .append("eventWeightMC")
    weightlistPUup                     .append("eventWeightMC")
    weightlistPUdown                   .append("eventWeightMC")
    weightlistFlatEffSF                .append("eventWeightMC")
    weightlistEffSFNormUpStat          .append("eventWeightMC")
    weightlistEffSFNormDownStat        .append("eventWeightMC")
    weightlistEffSFShapeUpEta          .append("eventWeightMC")
    weightlistEffSFShapeDownEta        .append("eventWeightMC")
    weightlistEffSFShapeUpPt           .append("eventWeightMC")
    weightlistEffSFShapeDownPt         .append("eventWeightMC")
    weightlistEffSFShapeUpPt40         .append("eventWeightMC")
    weightlistEffSFShapeDownPt40       .append("eventWeightMC")
    weightlistEffSFNormUpSys           .append("eventWeightMC")
    weightlistEffSFNormDownSys         .append("eventWeightMC")
    weightlistFinalSSV                 .append("eventWeightMC")
    weightlistBtagSFup                 .append("eventWeightMC")
    weightlistBtagSFdown               .append("eventWeightMC")
    weightlistMisTagSFup               .append("eventWeightMC")
    weightlistMisTagSFdown             .append("eventWeightMC")
    weightlistBTagSFShapeUpPt65        .append("eventWeightMC") 
    weightlistBTagSFShapeDownPt65      .append("eventWeightMC") 
    weightlistBTagSFShapeUpEta1p2      .append("eventWeightMC") 
    weightlistBTagSFShapeDownEta1p2    .append("eventWeightMC") 
    weightlistBTagSFFullShapeUpPt65    .append("eventWeightMC") 
    weightlistBTagSFFullShapeDownPt65  .append("eventWeightMC") 
    weightlistBTagSFFullShapeUpEta0p7  .append("eventWeightMC") 
    weightlistBTagSFFullShapeDownEta0p7.append("eventWeightMC") 
    weightlistBTagSFShapeUpPt100       .append("eventWeightMC") 
    weightlistBTagSFShapeDownPt100     .append("eventWeightMC") 
    weightlistBTagSFShapeUpEta0p7      .append("eventWeightMC") 
    weightlistBTagSFShapeDownEta0p7    .append("eventWeightMC") 

if(PUreweigthing):
    weightlistFinal                    .append(PUweight)
    weightlistNoBtagSFWeight           .append(PUweight)
    weightlistPUup                     .append(PUweightUp)
    weightlistPUdown                   .append(PUweightDown)
    weightlistFlatEffSF                .append(PUweight)
    weightlistEffSFNormUpStat          .append(PUweight)
    weightlistEffSFNormDownStat        .append(PUweight)
    weightlistEffSFShapeUpEta          .append(PUweight)
    weightlistEffSFShapeDownEta        .append(PUweight)
    weightlistEffSFShapeUpPt           .append(PUweight)
    weightlistEffSFShapeDownPt         .append(PUweight)
    weightlistEffSFShapeUpPt40         .append(PUweight)
    weightlistEffSFShapeDownPt40       .append(PUweight)
    weightlistEffSFNormUpSys           .append(PUweight)
    weightlistEffSFNormDownSys         .append(PUweight)
    weightlistFinalSSV                 .append(PUweight)
    weightlistBtagSFup                 .append(PUweight)
    weightlistBtagSFdown               .append(PUweight)
    weightlistMisTagSFup               .append(PUweight)
    weightlistMisTagSFdown             .append(PUweight)
    weightlistBTagSFShapeUpPt65        .append(PUweight) 
    weightlistBTagSFShapeDownPt65      .append(PUweight) 
    weightlistBTagSFShapeUpEta1p2      .append(PUweight) 
    weightlistBTagSFShapeDownEta1p2    .append(PUweight) 
    weightlistBTagSFFullShapeUpPt65    .append(PUweight) 
    weightlistBTagSFFullShapeDownPt65  .append(PUweight) 
    weightlistBTagSFFullShapeUpEta0p7  .append(PUweight) 
    weightlistBTagSFFullShapeDownEta0p7.append(PUweight) 
    weightlistBTagSFShapeUpPt100       .append(PUweight) 
    weightlistBTagSFShapeDownPt100     .append(PUweight) 
    weightlistBTagSFShapeUpEta0p7      .append(PUweight) 
    weightlistBTagSFShapeDownEta0p7    .append(PUweight) 

if(effSFReweigthing and decayChannel=="muon"):
    weightlistFinal                    .append("effSFMuonEventWeight")
    weightlistNoBtagSFWeight           .append("effSFMuonEventWeight")
    weightlistNoPUWeight               .append("effSFMuonEventWeight")
    weightlistPUup                     .append("effSFMuonEventWeightPUup")
    weightlistPUdown                   .append("effSFMuonEventWeightPUdown")
    weightlistFlatEffSF                .append("effSFMuonEventWeightFlatEffSF")
    weightlistEffSFNormUpStat          .append("effSFMuonEventWeightEffSFNormUpStat")
    weightlistEffSFNormDownStat        .append("effSFMuonEventWeightEffSFNormDownStat")
    weightlistEffSFShapeUpEta          .append("effSFMuonEventWeightEffSFShapeUpEta")
    weightlistEffSFShapeDownEta        .append("effSFMuonEventWeightEffSFShapeDownEta")
    weightlistEffSFShapeUpPt           .append("effSFMuonEventWeightEffSFShapeUpPt")
    weightlistEffSFShapeDownPt         .append("effSFMuonEventWeightEffSFShapeDownPt")
    weightlistEffSFShapeUpPt40         .append("effSFMuonEventWeightEffSFShapeUpPt40")
    weightlistEffSFShapeDownPt40       .append("effSFMuonEventWeightEffSFShapeDownPt40")
    weightlistEffSFNormUpSys           .append("effSFMuonEventWeightEffSFNormUpSys")
    weightlistEffSFNormDownSys         .append("effSFMuonEventWeightEffSFNormDownSys")
    weightlistFinalSSV                 .append("effSFMuonEventWeight")
    weightlistBtagSFup                 .append("effSFMuonEventWeight")
    weightlistBtagSFdown               .append("effSFMuonEventWeight")
    weightlistMisTagSFup               .append("effSFMuonEventWeight")
    weightlistMisTagSFdown             .append("effSFMuonEventWeight")
    weightlistBTagSFShapeUpPt65        .append("effSFMuonEventWeight")
    weightlistBTagSFShapeDownPt65      .append("effSFMuonEventWeight")
    weightlistBTagSFShapeUpEta1p2      .append("effSFMuonEventWeight")
    weightlistBTagSFShapeDownEta1p2    .append("effSFMuonEventWeight")
    weightlistBTagSFFullShapeUpPt65    .append("effSFMuonEventWeight")
    weightlistBTagSFFullShapeDownPt65  .append("effSFMuonEventWeight")
    weightlistBTagSFFullShapeUpEta0p7  .append("effSFMuonEventWeight")
    weightlistBTagSFFullShapeDownEta0p7.append("effSFMuonEventWeight")
    weightlistBTagSFShapeUpPt100       .append("effSFMuonEventWeight")
    weightlistBTagSFShapeDownPt100     .append("effSFMuonEventWeight")
    weightlistBTagSFShapeUpEta0p7      .append("effSFMuonEventWeight")
    weightlistBTagSFShapeDownEta0p7    .append("effSFMuonEventWeight")
if(effSFReweigthing and decayChannel=="electron"):
    weightlistFinal                    .append("effSFElectronEventWeight")
    weightlistNoBtagSFWeight           .append("effSFElectronEventWeight")
    weightlistNoPUWeight               .append("effSFElectronEventWeight")
    weightlistPUup                     .append("effSFElectronEventWeightPUup")
    weightlistPUdown                   .append("effSFElectronEventWeightPUdown")
    weightlistFlatEffSF                .append("effSFElectronEventWeightFlatEffSF")
    weightlistEffSFNormUpStat          .append("effSFElectronEventWeightEffSFNormUpStat")
    weightlistEffSFNormDownStat        .append("effSFElectronEventWeightEffSFNormDownStat")
    weightlistEffSFShapeUpEta          .append("effSFElectronEventWeightEffSFShapeUpEta")
    weightlistEffSFShapeDownEta        .append("effSFElectronEventWeightEffSFShapeDownEta")
    weightlistEffSFShapeUpPt           .append("effSFElectronEventWeightEffSFShapeUpPt")
    weightlistEffSFShapeDownPt         .append("effSFElectronEventWeightEffSFShapeDownPt")
    weightlistEffSFShapeUpPt40         .append("effSFElectronEventWeightEffSFShapeUpPt40")
    weightlistEffSFShapeDownPt40       .append("effSFElectronEventWeightEffSFShapeDownPt40")
    weightlistEffSFNormUpSys           .append("effSFElectronEventWeightEffSFNormUpSys")
    weightlistEffSFNormDownSys         .append("effSFElectronEventWeightEffSFNormDownSys")
    weightlistFinalSSV                 .append("effSFElectronEventWeight")
    weightlistBtagSFup                 .append("effSFElectronEventWeight")
    weightlistBtagSFdown               .append("effSFElectronEventWeight")
    weightlistMisTagSFup               .append("effSFElectronEventWeight")
    weightlistMisTagSFdown             .append("effSFElectronEventWeight")
    weightlistBTagSFShapeUpPt65        .append("effSFElectronEventWeight")
    weightlistBTagSFShapeDownPt65      .append("effSFElectronEventWeight")
    weightlistBTagSFShapeUpEta1p2      .append("effSFElectronEventWeight")
    weightlistBTagSFShapeDownEta1p2    .append("effSFElectronEventWeight")
    weightlistBTagSFFullShapeUpPt65    .append("effSFElectronEventWeight")
    weightlistBTagSFFullShapeDownPt65  .append("effSFElectronEventWeight")
    weightlistBTagSFFullShapeUpEta0p7  .append("effSFElectronEventWeight")
    weightlistBTagSFFullShapeDownEta0p7.append("effSFElectronEventWeight")
    weightlistBTagSFShapeUpPt100       .append("effSFElectronEventWeight")
    weightlistBTagSFShapeDownPt100     .append("effSFElectronEventWeight")
    weightlistBTagSFShapeUpEta0p7      .append("effSFElectronEventWeight")
    weightlistBTagSFShapeDownEta0p7    .append("effSFElectronEventWeight")
if(BtagReweigthing):
    weightlistFinal                    .append("bTagSFEventWeight")
    weightlistNoPUWeight               .append("bTagSFEventWeight")
    weightlistPUup                     .append("bTagSFEventWeight")
    weightlistPUdown                   .append("bTagSFEventWeight")
    weightlistFlatEffSF                .append("bTagSFEventWeight")
    weightlistEffSFNormUpStat          .append("bTagSFEventWeight")
    weightlistEffSFNormDownStat        .append("bTagSFEventWeight")
    weightlistEffSFShapeUpEta          .append("bTagSFEventWeight")
    weightlistEffSFShapeDownEta        .append("bTagSFEventWeight")
    weightlistEffSFShapeUpPt           .append("bTagSFEventWeight")
    weightlistEffSFShapeDownPt         .append("bTagSFEventWeight")
    weightlistEffSFShapeUpPt40         .append("bTagSFEventWeight")
    weightlistEffSFShapeDownPt40       .append("bTagSFEventWeight")
    weightlistEffSFNormUpSys           .append("bTagSFEventWeight")
    weightlistEffSFNormDownSys         .append("bTagSFEventWeight")
    weightlistFinalSSV                 .append("bTagSFEventWeightSSV")
    weightlistBtagSFup                 .append("bTagSFEventWeightBTagSFUp")
    weightlistBtagSFdown               .append("bTagSFEventWeightBTagSFDown")
    weightlistMisTagSFup               .append("bTagSFEventWeightMisTagSFUp")
    weightlistMisTagSFdown             .append("bTagSFEventWeightMisTagSFDown")
    weightlistBTagSFShapeUpPt65        .append("bTagSFEventWeightBTagSFShapeUpPt65")
    weightlistBTagSFShapeDownPt65      .append("bTagSFEventWeightBTagSFShapeDownPt65")
    weightlistBTagSFShapeUpEta1p2      .append("bTagSFEventWeightBTagSFShapeUpEta1p2")
    weightlistBTagSFShapeDownEta1p2    .append("bTagSFEventWeightBTagSFShapeDownEta1p2")
    weightlistBTagSFFullShapeUpPt65    .append("bTagSFEventWeightBTagSFFullShapeUpPt65")
    weightlistBTagSFFullShapeDownPt65  .append("bTagSFEventWeightBTagSFFullShapeDownPt65")
    weightlistBTagSFFullShapeUpEta0p7  .append("bTagSFEventWeightBTagSFFullShapeUpEta0p7")
    weightlistBTagSFFullShapeDownEta0p7.append("bTagSFEventWeightBTagSFFullShapeDownEta0p7")
    weightlistBTagSFShapeUpPt100       .append("bTagSFEventWeightBTagSFShapeUpPt100")
    weightlistBTagSFShapeDownPt100     .append("bTagSFEventWeightBTagSFShapeDownPt100")
    weightlistBTagSFShapeUpEta0p7      .append("bTagSFEventWeightBTagSFShapeUpEta0p7")
    weightlistBTagSFShapeDownEta0p7    .append("bTagSFEventWeightBTagSFShapeDownEta0p7")

## multiply all event weights
## a) default
process.eventWeightNoBtagSFWeight           = process.eventWeightMultiplier.clone(eventWeightTags = weightlistNoBtagSFWeight)
process.eventWeightNoPUWeight               = process.eventWeightMultiplier.clone(eventWeightTags = weightlistNoPUWeight)
process.eventWeightFinal                    = process.eventWeightMultiplier.clone(eventWeightTags = weightlistFinal)

process.eventWeightFinalSSV                 = process.eventWeightMultiplier.clone(eventWeightTags = weightlistFinalSSV)

## b) for systematics
process.eventWeightPUup                     = process.eventWeightMultiplier.clone(eventWeightTags = weightlistPUup)
process.eventWeightPUdown                   = process.eventWeightMultiplier.clone(eventWeightTags = weightlistPUdown)
process.eventWeightFlatEffSF                = process.eventWeightMultiplier.clone(eventWeightTags = weightlistFlatEffSF)
process.eventWeightEffSFNormUpStat          = process.eventWeightMultiplier.clone(eventWeightTags = weightlistEffSFNormUpStat)
process.eventWeightEffSFNormDownStat        = process.eventWeightMultiplier.clone(eventWeightTags = weightlistEffSFNormDownStat)
process.eventWeightEffSFShapeUpEta          = process.eventWeightMultiplier.clone(eventWeightTags = weightlistEffSFShapeUpEta)
process.eventWeightEffSFShapeDownEta        = process.eventWeightMultiplier.clone(eventWeightTags = weightlistEffSFShapeDownEta)
process.eventWeightEffSFShapeUpPt           = process.eventWeightMultiplier.clone(eventWeightTags = weightlistEffSFShapeUpPt)
process.eventWeightEffSFShapeDownPt         = process.eventWeightMultiplier.clone(eventWeightTags = weightlistEffSFShapeDownPt)
process.eventWeightEffSFShapeUpPt40         = process.eventWeightMultiplier.clone(eventWeightTags = weightlistEffSFShapeUpPt40)
process.eventWeightEffSFShapeDownPt40       = process.eventWeightMultiplier.clone(eventWeightTags = weightlistEffSFShapeDownPt40)
process.eventWeightEffSFNormUpSys           = process.eventWeightMultiplier.clone(eventWeightTags = weightlistEffSFNormUpSys)
process.eventWeightEffSFNormDownSys         = process.eventWeightMultiplier.clone(eventWeightTags = weightlistEffSFNormDownSys)
process.eventWeightBtagSFup                 = process.eventWeightMultiplier.clone(eventWeightTags = weightlistBtagSFup)
process.eventWeightBtagSFdown               = process.eventWeightMultiplier.clone(eventWeightTags = weightlistBtagSFdown)
process.eventWeightMisTagSFup               = process.eventWeightMultiplier.clone(eventWeightTags = weightlistMisTagSFup)
process.eventWeightMisTagSFdown             = process.eventWeightMultiplier.clone(eventWeightTags = weightlistMisTagSFdown)
process.eventWeightBTagSFShapeUpPt65        = process.eventWeightMultiplier.clone(eventWeightTags = weightlistBTagSFShapeUpPt65)
process.eventWeightBTagSFShapeDownPt65      = process.eventWeightMultiplier.clone(eventWeightTags = weightlistBTagSFShapeDownPt65)
process.eventWeightBTagSFShapeUpEta1p2      = process.eventWeightMultiplier.clone(eventWeightTags = weightlistBTagSFShapeUpEta1p2)
process.eventWeightBTagSFShapeDownEta1p2    = process.eventWeightMultiplier.clone(eventWeightTags = weightlistBTagSFShapeDownEta1p2)
process.eventWeightBTagSFFullShapeUpPt65    = process.eventWeightMultiplier.clone(eventWeightTags = weightlistBTagSFFullShapeUpPt65)
process.eventWeightBTagSFFullShapeDownPt65  = process.eventWeightMultiplier.clone(eventWeightTags = weightlistBTagSFFullShapeDownPt65)
process.eventWeightBTagSFFullShapeUpEta0p7  = process.eventWeightMultiplier.clone(eventWeightTags = weightlistBTagSFFullShapeUpEta0p7)
process.eventWeightBTagSFFullShapeDownEta0p7= process.eventWeightMultiplier.clone(eventWeightTags = weightlistBTagSFFullShapeDownEta0p7)
process.eventWeightBTagSFShapeUpPt100       = process.eventWeightMultiplier.clone(eventWeightTags = weightlistBTagSFShapeUpPt100)
process.eventWeightBTagSFShapeDownPt100     = process.eventWeightMultiplier.clone(eventWeightTags = weightlistBTagSFShapeDownPt100)
process.eventWeightBTagSFShapeUpEta0p7      = process.eventWeightMultiplier.clone(eventWeightTags = weightlistBTagSFShapeUpEta0p7)
process.eventWeightBTagSFShapeDownEta0p7    = process.eventWeightMultiplier.clone(eventWeightTags = weightlistBTagSFShapeDownEta0p7)
    
# use weight in single and double object analyzer modules for central values
modulelist   = process.analyzers_().keys()
genModules1  = process.kinFitGen.moduleNames()
genModules2  = process.kinFitGenPhaseSpace.moduleNames()
PUModuleList = set(['PUControlDistributionsDefault','PUControlDistributionsBeforeBtagging','PUControlDistributionsAfterBtagging'])

# a) Default configuration for event weights for all modules (gen + reco)
#    No re-weighting at all is to be applied
if(runningOnData=="MC" and not PUreweigthing and not effSFReweigthing and not BtagReweigthing):
    for module in modulelist:
        if module not in PUModuleList:
            if(MCweighting):
                getattr(process,module).weight=cms.InputTag("eventWeightMC")
            else:
                getattr(process,module).weight=cms.InputTag("")
# b) Re-weighting (PU and/or effSF) for reco-modules
if(runningOnData=="MC" and (PUreweigthing or effSFReweigthing)):
    if(PUreweigthing):
        print " All reco-modules will use PU event weights"
    if(effSFReweigthing):
        print " All reco-modules will use eff SF event weights"
    for module in modulelist:
        # exclude PU control modules and gen-modules
        if(module not in PUModuleList and module not in genModules1 and module not in genModules2):
            # 'eventWeightNoBtagSFWeight' treats the correct combination of 'PUreweigthing' and 'effSFReweigthing'
            getattr(process,module).weight=cms.InputTag("eventWeightNoBtagSFWeight")
        
# c) Re-weighting (b-tag SF) only for reco-modules after btagging
if(runningOnData=="MC" and BtagReweigthing):
    btagModules1 = process.monitorKinematicsAfterBtagging.moduleNames()
    btagModules2 = process.kinFit.moduleNames()
    btagModules3 = process.monitorElectronKinematicsAfterBtagging.moduleNames()
    print
    print " The following reco-modules will use additionally the b-tag event weights:"
    print btagModules1
    print btagModules2
    print btagModules3
    # 'eventWeightFinal' combines the correct combination of 'PUreweigthing' and 'effSFReweigthing' with the b-tag SF
    for module1 in btagModules1:
        getattr(process,module1).weight=cms.InputTag("eventWeightFinal")
    for module2 in btagModules2:
        if(not module2=='makeGenLevelBJets'):
            getattr(process,module2).weight=cms.InputTag("eventWeightFinal")
    for module3 in btagModules3:
        getattr(process,module3).weight=cms.InputTag("eventWeightFinal")
    print
        
# d) Re-weighting of gen-modules - PU reweighting only
if(runningOnData=="MC" and PUreweigthing):
    print
    print " The following gen-modules will only use the PU reweighting:"
    if(MCweighting):
        print "(and MC weight!)"
    #print genModules1
    #print genModules2
    for module1 in genModules1:
        if(not module1=='makeGenLevelBJets'):
            print module1
            if(MCweighting):
                getattr(process,module1).weight=cms.InputTag("eventWeightMCandPUweight")
            else:            
                getattr(process,module1).weight=PUweight
    for module2 in genModules2:
        if(not module2=='makeGenLevelBJets'):
            print module2
            if(MCweighting):
                getattr(process,module2).weight=cms.InputTag("eventWeightMCandPUweight")
            else:            
                getattr(process,module2).weight=PUweight
    genModules3 = process.kinFitGenPhaseSpaceHad.moduleNames()
    #print genModules3
    for module3 in genModules3:
        if(not module3=='makeGenLevelBJets'):
            print module3
            if(MCweighting):
                getattr(process,module3).weight=cms.InputTag("eventWeightMCandPUweight")
            else:            
                getattr(process,module3).weight=PUweight
    genModules4 = process.hadLvObjectMonitoring.moduleNames()
    #print genModules4
    for module4 in genModules4:
        if(not module4=='makeGenLevelBJets'):
            print module4
            if(MCweighting):
                getattr(process,module4).weight=cms.InputTag("eventWeightMCandPUweight")
            else:            
                getattr(process,module4).weight=PUweight
    if(additionalEventWeights and eventFilter=='signal only'):
        print "those gen modules are also cloned in order to also use NoPU, PUup and PUdown event weights "

        process.GENFULLanalyzers     =cms.Sequence(process.dummy)
        process.GENpartonPSanalyzers =cms.Sequence(process.dummy)
        process.GENhadronPSanalyzers =cms.Sequence(process.dummy)
        process.GENFULLbanalyzers    =cms.Sequence(process.dummy)
        process.GENpartonPSbanalyzers=cms.Sequence(process.dummy)
        process.GENhadronPSbanalyzers=cms.Sequence(process.dummy)
        genSystExt=["NoPUWeight", "PUup", "PUdown"]
        for sys in genSystExt:
            # get correct weight
            if(MCweighting):
                weightTagName=cms.InputTag("eventWeightMC")
            else:
                weightTagName=cms.InputTag("")
            if(not sys.find("PUup") == -1):
                if(MCweighting):
                    weightTagName=cms.InputTag("eventWeightMCandPUweightUp")
                else:
                    weightTagName=PUweightUp
            elif(not sys.find("PUdown") == -1):
                if(MCweighting):
                    weightTagName=cms.InputTag("eventWeightMCandPUweightDown")
                else:
                    weightTagName=PUweightDown
            # create plots for full PS
            setattr(process,"analyzeTopPartonLevelKinematics"+sys, process.analyzeTopPartonLevelKinematics.clone(weight=weightTagName))
            getattr(process,"analyzeTopPartonLevelKinematics"+sys).analyze.useTree = False
            #setattr(process,"analyzeTopPartonLevelKinematicsBjets"+sys, process.analyzeTopPartonLevelKinematicsBjets.clone(weight=weightTagName))
            #getattr(process,"analyzeTopPartonLevelKinematicsBjets"+sys).useTree = False
            #setattr(process,"analyzeTopPartonLevelKinematicsLepton"+sys, process.analyzeTopPartonLevelKinematicsLepton.clone(weight=weightTagName))
            #getattr(process,"analyzeTopPartonLevelKinematicsLepton"+sys).useTree = False
            #setattr(process,"compositedPartonGen"+sys, process.compositedPartonGen.clone(weight=weightTagName))
            #getattr(process,"compositedPartonGen"+sys).useTree = False 
            
            # create plots for parton level PS
            setattr(process,"analyzeTopPartonLevelKinematicsPhaseSpace"+sys, process.analyzeTopPartonLevelKinematicsPhaseSpace.clone(weight=weightTagName))
            getattr(process,"analyzeTopPartonLevelKinematicsPhaseSpace"+sys).analyze.useTree = False
            #setattr(process,"analyzeTopPartonLevelKinematicsBjetsPhaseSpace"+sys, process.analyzeTopPartonLevelKinematicsBjetsPhaseSpace.clone(weight=weightTagName))
            #getattr(process,"analyzeTopPartonLevelKinematicsBjetsPhaseSpace"+sys).useTree = False            
            #setattr(process,"analyzeTopPartonLevelKinematicsLeptonPhaseSpace"+sys, process.analyzeTopPartonLevelKinematicsLeptonPhaseSpace.clone(weight=weightTagName))
            #getattr(process,"analyzeTopPartonLevelKinematicsLeptonPhaseSpace"+sys).useTree = False
            #setattr(process,"compositedPartonGenPhaseSpace"+sys, process.compositedPartonGenPhaseSpace.clone(weight=weightTagName))
            #getattr(process,"compositedPartonGenPhaseSpace"+sys).useTree = False 

            # create plots for hadron level PS
            #setattr(process,"analyzeTopHadronLevelKinematicsPhaseSpace"+sys, process.analyzeTopHadronLevelKinematicsPhaseSpace.clone(weight=weightTagName))
            #getattr(process,"analyzeTopHadronLevelKinematicsPhaseSpace"+sys).analyze.useTree = False
            setattr(process,"analyzeTopHadronLevelKinematicsBjetsPhaseSpace"+sys, process.analyzeTopHadronLevelKinematicsBjetsPhaseSpace.clone(weight=weightTagName))
            getattr(process,"analyzeTopHadronLevelKinematicsBjetsPhaseSpace"+sys).useTree = False           
            setattr(process,"analyzeTopHadronLevelKinematicsLeptonPhaseSpace"+sys, process.analyzeTopHadronLevelKinematicsLeptonPhaseSpace.clone(weight=weightTagName))
            getattr(process,"analyzeTopHadronLevelKinematicsLeptonPhaseSpace"+sys).useTree = False
            setattr(process,"compositedHadronGenPhaseSpace"+sys, process.compositedHadronGenPhaseSpace.clone(weight=weightTagName))
            #getattr(process,"compositedHadronGenPhaseSpace"+sys).useTree = False
            
            # analyzer to be added to sequence
            process.GENFULLanalyzers    *=getattr(process,"analyzeTopPartonLevelKinematics"+sys)
            process.GENpartonPSanalyzers*=getattr(process,"analyzeTopPartonLevelKinematicsPhaseSpace"+sys)
            #process.GENhadronPSanalyzers*=getattr(process,"analyzeTopHadronLevelKinematicsPhaseSpace"+sys)
            #process.GENFULLbanalyzers    *=getattr(process,"analyzeTopPartonLevelKinematicsBjets"+sys)
            #process.GENpartonPSbanalyzers*=getattr(process,"analyzeTopPartonLevelKinematicsBjetsPhaseSpace"+sys)
            process.GENhadronPSbanalyzers*=getattr(process,"analyzeTopHadronLevelKinematicsBjetsPhaseSpace"+sys)
            #process.GENFULLbanalyzers    *=getattr(process,"analyzeTopPartonLevelKinematicsLepton"+sys)
            #process.GENpartonPSbanalyzers*=getattr(process,"analyzeTopPartonLevelKinematicsLeptonPhaseSpace"+sys)
            process.GENhadronPSbanalyzers*=getattr(process,"analyzeTopHadronLevelKinematicsLeptonPhaseSpace"+sys)
            #process.GENFULLbanalyzers    *=getattr(process,"compositedPartonGen"+sys)
            #process.GENpartonPSbanalyzers*=getattr(process,"compositedPartonGenPhaseSpace"+sys)
            process.GENhadronPSbanalyzers*=getattr(process,"compositedHadronGenPhaseSpace"+sys)
            

        process.GENFULLanalyzers.remove(process.dummy)
        process.GENpartonPSanalyzers.remove(process.dummy)
        process.GENhadronPSanalyzers.remove(process.dummy)
        process.GENFULLbanalyzers.remove(process.dummy)
        process.GENpartonPSbanalyzers.remove(process.dummy)
        process.GENhadronPSbanalyzers.remove(process.dummy)
        
        # add them to sequence
        process.kinFitGen              *= (process.GENFULLanalyzers *
                                           process.GENFULLbanalyzers)
        
        process.kinFitGenPhaseSpace    *= (process.GENpartonPSanalyzers *
                                           process.GENpartonPSbanalyzers)
        
        process.kinFitGenPhaseSpaceHad *= (process.GENhadronPSanalyzers *
                                           process.GENhadronPSbanalyzers)
        
elif(runningOnData=="MC" and not PUreweigthing):
    for module1 in genModules1:
        if(not module1=='makeGenLevelBJets'):
            if(MCweighting):
                getattr(process,module1).weight=cms.InputTag("eventWeightMC")
            else:
                getattr(process,module1).weight=cms.InputTag("")
            
    for module2 in genModules2:
        if(not module2=='makeGenLevelBJets'):
            if(MCweighting):
                getattr(process,module2).weight=cms.InputTag("eventWeightMC")
            else:
                getattr(process,module2).weight=cms.InputTag("")
##     for module3 in genModules3:
##        if(not module3=='makeGenLevelBJets'):
##            if(MCweighting):
##                getattr(process,module3).weight=cms.InputTag("eventWeightMC")
##            else:
##                getattr(process,module3).weight=cms.InputTag("")
##     for module4 in genModules4:
##        if(not module4=='makeGenLevelBJets'):
##            if(MCweighting):
##                getattr(process,module4).weight=cms.InputTag("eventWeightMC")
##            else:
##                getattr(process,module4).weight=cms.InputTag("")

# e) PU Modules - special configuration because of different possible PU weights that are handled in the modules themselves
if(runningOnData=="MC"):
    if(effSFReweigthing and decayChannel=="muon"):
        if(MCweighting):
            weightlistMCandSFMuon=cms.VInputTag()
            weightlistMCandSFMuon.append("eventWeightMC")
            weightlistMCandSFMuon.append("effSFMuonEventWeight")
            process.eventWeightMCandSFMuon=process.eventWeightMultiplier.clone(eventWeightTags = weightlistMCandSFMuon)
            process.PUControlDistributionsDefault.DefEventWeight        = cms.InputTag("eventWeightMCandSFMuon")
            process.PUControlDistributionsBeforeBtagging.DefEventWeight = cms.InputTag("eventWeightMCandSFMuon")
        else:
            process.PUControlDistributionsDefault.DefEventWeight        = cms.InputTag("effSFMuonEventWeight")
            process.PUControlDistributionsBeforeBtagging.DefEventWeight = cms.InputTag("effSFMuonEventWeight")
    elif(effSFReweigthing and decayChannel=="electron"):
        if(MCweighting):
            weightlistMCandSFElectron=cms.VInputTag()
            weightlistMCandSFElectron.append("eventWeightMC")
            weightlistMCandSFElectron.append("effSFElectronEventWeight")
            process.eventWeightMCandSFElectron=process.eventWeightMultiplier.clone(eventWeightTags = weightlistMCandSFElectron)
            process.PUControlDistributionsDefault.DefEventWeight        = cms.InputTag("eventWeightMCandSFElectron")
            process.PUControlDistributionsBeforeBtagging.DefEventWeight = cms.InputTag("eventWeightMCandSFElectron")
        else:
            process.PUControlDistributionsDefault.DefEventWeight        = cms.InputTag("effSFElectronEventWeight")
            process.PUControlDistributionsBeforeBtagging.DefEventWeight = cms.InputTag("effSFElectronEventWeight")
    else:
        if(MCweighting):
            process.PUControlDistributionsDefault.DefEventWeight        = cms.InputTag("eventWeightMC")
            process.PUControlDistributionsBeforeBtagging.DefEventWeight = cms.InputTag("eventWeightMC")
        else:
            process.PUControlDistributionsDefault.DefEventWeight        = cms.InputTag("")
            process.PUControlDistributionsBeforeBtagging.DefEventWeight = cms.InputTag("")
    # 'eventWeightNoPUWeight' combines the correct combination of 'effSFReweigthing' and 'BtagReweigthing'
    process.PUControlDistributionsAfterBtagging.DefEventWeight  = cms.InputTag("eventWeightNoPUWeight")    
    ## copies of TopRecoKinematicsKinFit analyzers with varied weights for monitoring and systematic unc.
if(runningOnData=="MC" and applyKinFit==True and additionalEventWeights):

    ## list of systematic variations via event weights
    process.weights=cms.Sequence(process.dummy)
    process.myAnalyzers=cms.Sequence(process.dummy)
    process.myBanalyzers=cms.Sequence(process.dummy)
    process.myLepanalyzers=cms.Sequence(process.dummy)
    process.myMixAnalyzers=cms.Sequence(process.dummy)
    systExt=["NoWeight", "OnlyPUWeight", "NoBtagSFWeight",
             "PUup", "PUdown",
             "FlatEffSF", "EffSFNormUpStat", "EffSFNormDownStat", "EffSFShapeUpEta", "EffSFShapeDownEta",
             "EffSFShapeUpPt", "EffSFShapeDownPt", "EffSFShapeUpPt40", "EffSFShapeDownPt40",
             "EffSFNormUpSys", "EffSFNormDownSys",
             "BtagSFup", "BtagSFdown", "MisTagSFup", "MisTagSFdown", "BTagSFShapeUpPt65", "BTagSFShapeDownPt65", "BTagSFShapeUpEta1p2", "BTagSFShapeDownEta1p2",
             "BTagSFFullShapeUpPt65", "BTagSFFullShapeDownPt65", "BTagSFFullShapeUpEta0p7", "BTagSFFullShapeDownEta0p7",
             "BTagSFShapeUpPt100", "BTagSFShapeDownPt100", "BTagSFShapeUpEta0p7", "BTagSFShapeDownEta0p7"]
    # loop them
    for sys in systExt:
        # take care of correct name and non existing weights
        weightTagName=cms.InputTag("eventWeight"+sys)
        if(not sys.find("OnlyPUWeight") == -1):
            if(MCweighting):
                weightTagName=cms.InputTag("eventWeightMCandPUweight")
            else:
                weightTagName=PUweight
        noweight=False
        if(not PUreweigthing and not sys.find("PU") == -1):
            noweight=True
        if(not PUreweigthing and not effSFReweigthing and not sys.find("NoBtagSFWeight") == -1):
            noweight=True
        if(not sys.find("NoWeight") == -1):
            noweight=True
        if(noweight):
            if(MCweighting):
                weightTagName=cms.InputTag("eventWeightMC")
            else:
                weightTagName=cms.InputTag("")
        # create plots for standard analyzer
        setattr(process,"analyzeTopRecoKinematicsKinFit"+sys, process.analyzeTopRecoKinematicsKinFit.clone(weight=weightTagName))
        getattr(process,"analyzeTopRecoKinematicsKinFit"+sys).analyze.useTree = False
        # create plots for top/antitop analyzer
        setattr(process,"analyzeTopRecoKinematicsKinFitTopAntitop"+sys, process.analyzeTopRecoKinematicsKinFitTopAntitop.clone(weight=weightTagName))
        getattr(process,"analyzeTopRecoKinematicsKinFitTopAntitop"+sys).analyze.useTree = False
        # create plots for lepton/bjets analyzer
        setattr(process,"analyzeTopRecoKinematicsBjets"+sys, process.analyzeTopRecoKinematicsBjets.clone(weight=weightTagName, useTree = False))
        setattr(process,"analyzeTopRecoKinematicsLepton"+sys, process.analyzeTopRecoKinematicsLepton.clone(weight=weightTagName, useTree = False))
        # create plots for mixed objects analyzer  
        setattr(process,"compositedKinematicsKinFit"+sys, process.compositedKinematicsKinFit.clone(weight=weightTagName))        
        # weights to be added to sequence
        if(sys.find("NoWeight") == -1 and sys.find("OnlyPUWeight") == -1 and sys.find("NoBtagSFWeight") == -1):
            process.weights*=getattr(process,"eventWeight"+sys)
        # analyzer to be added to sequence
        process.myAnalyzers*=getattr(process,"analyzeTopRecoKinematicsKinFit"+sys)
        process.myMixAnalyzers*=getattr(process,"compositedKinematicsKinFit"+sys)
        #analyzer*=getattr(process,"analyzeTopRecoKinematicsKinFitTopAntitop"+sys)
        if(sys.find("NoWeight") == -1 and sys.find("OnlyPUWeight") == -1 and sys.find("NoBtagSFWeight") == -1):
            process.myBanalyzers*=getattr(process,"analyzeTopRecoKinematicsBjets"+sys)
            process.myLepanalyzers*=getattr(process,"analyzeTopRecoKinematicsLepton"+sys)
        
    process.weights   .remove(process.dummy)
    process.myAnalyzers .remove(process.dummy)
    process.myBanalyzers.remove(process.dummy)
    process.myLepanalyzers.remove(process.dummy)
    process.myMixAnalyzers.remove(process.dummy)
        
    ## add to Sequence
    # a) event weights 
    process.kinFit.replace(process.analyzeTopRecoKinematicsKinFit, 
                           process.analyzeTopRecoKinematicsKinFit*                            
                           process.bTagSFEventWeightBTagSFUp    *
                           process.bTagSFEventWeightBTagSFDown  *
                           process.bTagSFEventWeightMisTagSFUp  *
                           process.bTagSFEventWeightMisTagSFDown*
                           process.bTagSFEventWeightBTagSFShapeUpPt65    *
                           process.bTagSFEventWeightBTagSFShapeDownPt65  *
                           process.bTagSFEventWeightBTagSFShapeUpEta1p2    *
                           process.bTagSFEventWeightBTagSFShapeDownEta1p2  *
                           process.bTagSFEventWeightBTagSFFullShapeUpPt65    *
                           process.bTagSFEventWeightBTagSFFullShapeDownPt65  *
                           process.bTagSFEventWeightBTagSFFullShapeUpEta0p7    *
                           process.bTagSFEventWeightBTagSFFullShapeDownEta0p7  *
                           process.bTagSFEventWeightBTagSFShapeUpPt100    *
                           process.bTagSFEventWeightBTagSFShapeDownPt100  *
                           process.bTagSFEventWeightBTagSFShapeUpEta0p7    *
                           process.bTagSFEventWeightBTagSFShapeDownEta0p7
                           )
    if(decayChannel=="muon"):
	process.kinFit.replace(process.analyzeTopRecoKinematicsKinFit, 
                               process.analyzeTopRecoKinematicsKinFit*
			       process.effSFMuonEventWeightPUup*
			       process.effSFMuonEventWeightPUdown*
                               process.effSFMuonEventWeightFlatEffSF*           
                               process.effSFMuonEventWeightEffSFNormUpStat*      
                               process.effSFMuonEventWeightEffSFNormDownStat*    
                               process.effSFMuonEventWeightEffSFShapeUpEta*  
                               process.effSFMuonEventWeightEffSFShapeDownEta*
                               process.effSFMuonEventWeightEffSFShapeUpPt*  
                               process.effSFMuonEventWeightEffSFShapeDownPt*
			       process.effSFMuonEventWeightEffSFShapeUpPt40*  
                               process.effSFMuonEventWeightEffSFShapeDownPt40*
                               process.effSFMuonEventWeightEffSFNormUpSys    *
                               process.effSFMuonEventWeightEffSFNormDownSys
                               )
    elif(decayChannel=="electron"):
        process.kinFit.replace(process.analyzeTopRecoKinematicsKinFit, 
                               process.analyzeTopRecoKinematicsKinFit*
                               process.effSFElectronEventWeightPUup*      
                               process.effSFElectronEventWeightPUdown*
			       process.effSFElectronEventWeightFlatEffSF*           
                               process.effSFElectronEventWeightEffSFNormUpStat*      
                               process.effSFElectronEventWeightEffSFNormDownStat*    
                               process.effSFElectronEventWeightEffSFShapeUpEta*  
                               process.effSFElectronEventWeightEffSFShapeDownEta*
                               process.effSFElectronEventWeightEffSFShapeUpPt*  
                               process.effSFElectronEventWeightEffSFShapeDownPt*
			       process.effSFElectronEventWeightEffSFShapeUpPt40*  
                               process.effSFElectronEventWeightEffSFShapeDownPt40*
                               process.effSFElectronEventWeightEffSFNormUpSys    *
                               process.effSFElectronEventWeightEffSFNormDownSys 
                               )
    process.kinFit.replace(process.bTagSFEventWeightBTagSFShapeDownEta0p7,
                           process.bTagSFEventWeightBTagSFShapeDownEta0p7      *
                           ## additional weights
                           process.weights *
                           ## additional std top analyzers
                           process.myAnalyzers                         
                           )
    process.kinFit.replace(process.analyzeTopRecoKinematicsBjets,
                           process.analyzeTopRecoKinematicsBjets *
                           ## additional b-hadron jet analyzers
                           process.myBanalyzers                     
                           )
    process.kinFit.replace(process.analyzeTopRecoKinematicsLepton,
                           process.analyzeTopRecoKinematicsLepton *
                           ## additional b-hadron jet analyzers
                           process.myLepanalyzers                     
                           )

    process.kinFit.replace(process.compositedKinematicsKinFit,
                           process.compositedKinematicsKinFit *
                           ## additional mixed objects analyzers
                           process.myMixAnalyzers                     
                           )
    
## SSV event weights
if(runningOnData=="MC" and BtagReweigthing):
    SSVModules = process.SSVMonitoring.moduleNames()
    print
    print " The following reco-modules will use the SSV b-tag event weights:"
    print SSVModules
    # 'eventWeightFinal' combines the correct combination of 'PUreweigthing' and 'effSFReweigthing' with the SSV b-tag SF
    for module in SSVModules:
        getattr(process,module).weight=cms.InputTag("eventWeightFinalSSV")
        
## test isolation
process.newvertexSelectedMuons=process.vertexSelectedMuons.clone(src="noCutPatMuons")
process.newtrackMuons=process.trackMuons.clone(src="newvertexSelectedMuons")
#process.testIsoMuons=process.goldenMuons.clone(muons = "newtrackMuons")
process.testIsoMuons=process.newtrackMuons.clone()

process.testIsoMuonSelection= process.muonSelection.clone (src = 'testIsoMuons', minNumber = 1, maxNumber = 99999999)
process.testIsoMuonQuality  = process.tightMuonQualityTagged.clone(src = 'testIsoMuons')

process.newvertexSelectedElectrons=process.vertexSelectedElectrons.clone(src="noCutPatElectrons")
process.testIsoElectrons=process.tightElectronsEJ.clone(
    src = "newvertexSelectedElectrons",
    cut = 'et > 30. &abs(eta) <  2.1  &( abs(superCluster.eta) < 1.4442   |  abs(superCluster.eta) > 1.5660 ) &abs(dB)  <  0.02 &test_bit( electronID("eidHyperTight1MC"), 0 )'
    )

process.testIsoElectronSelection= process.convElecRejection.clone (src = 'testIsoElectrons', minNumber = 1, maxNumber = 99999999)
process.testIsoElectronQuality  = process.tightElectronQualityTagged.clone(src = 'testIsoElectrons')


## test event particle listing
## using ParticleListDrawer from PhysicsTools/HepMCCandAlgos
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint  = cms.untracked.int32(10)
) 


## ---
##    run the final sequences
## ---
## standard sequence for cross section analyis and detailed cut monitoring
process.p1 = cms.Path(## gen event selection (decay channel) and the trigger selection (hltFilter)
                      process.filterSequence                        *
                      ## PV event selection
                      #process.PVSelection                           *
                      ## introduce some collections
                      process.semiLeptonicSelection                 *
                      process.isolatedGenLeptons                    *
                      process.semiLeptGenCollections                *
                      process.cleanedGenJetCollection               *
                      process.makeGenLevelBJets                     *
                      process.noOverlapGenJetCollection             *
                      process.altermakeGenLevelBJets                *
                      process.bjetGenJets                           *
                      process.visibleHadronPS                       *
                      process.additionalGenJet                      *
                      ## create PU event weights
                      process.makeEventWeightsPU                    *
                      ## create shape distortion event weights
                      process.eventWeightDileptonModelVariation     *
                      process.eventWeightPUDistort                  *
                      process.eventWeightPUupDistort                *
                      process.eventWeightPUdownDistort              *
		      ## create effSF eventWeight
		      process.effSFMuonEventWeight                  *
		      ## multiply event weights
		      process.eventWeightNoBtagSFWeight             *
                      ## muon selection
                      process.muonCuts                              *
                      ## veto on additional leptons
                      process.secondMuonVeto                        *
                      process.electronVeto                          *
                      ## jet selection and monitoring
                      process.jetSelection                          *
                      ## monitoring before b-tagging
                      process.monitorKinematicsBeforeBtagging       *
                      process.PUControlDistributionsBeforeBtagging  *                      
                      ## b-tagging
                      process.btagSelection1                        *
                      process.btagSelection                         *
                      ## create PU event weights
                      process.bTagSFEventWeight                     *
                      ## create combined weights
                      process.eventWeightNoPUWeight                 *
                      ## create combined weight
                      process.eventWeightFinal                      *
                      ## monitor after b-tagging
                      process.monitorKinematicsAfterBtagging        *
                      process.PUControlDistributionsAfterBtagging   *
                      ## apply kinematic fit
                      process.kinFit                                *
                      ## monitor after kinematic fit
                      process.monitorJetKinematicsAfterKinFit
                      )

if(applyKinFit==False or eventFilter!="signal only"):
    # rm gen collection- only needed when processing signal for gen-reco correlations
    process.p1.remove(process.makeGenLevelBJets            )
    process.p1.remove(process.noOverlapGenJetCollection    )
    process.p1.remove(process.isolatedGenLeptons           )
    process.p1.remove(process.semiLeptGenCollections       )
    process.p1.remove(process.bjetGenJets                  )
    process.p1.remove(process.visibleHadronPS              )
    process.p1.remove(process.cleanedGenJetCollection      )
    process.p1.remove(process.altermakeGenLevelBJets       )
    process.p1.remove(process.additionalGenJet             )
    process.p1.remove(process.dummy)

if(runningOnData=="data"):
    process.p1.remove(process.isolatedGenMuons)
    process.p1.remove(process.semiLeptGenCollections)

## analysis with SSV btagger
process.p2 = cms.Path(process.dummy)
## if(not cutflowSynch):
##     process.p2 = cms.Path(## gen event selection (decay channel) and the trigger selection (hltFilter)
##                           process.filterSequence                        *
##                           ## PV event selection
##                           #process.PVSelection                           *
##                           ## introduce some collections
##                           process.semiLeptonicSelection                 *
##                           ## create PU event weights
##                           process.makeEventWeightsPU                    *
##                           ## create shape distortion event weights
##                           process.eventWeightDileptonModelVariation     *
##                           process.eventWeightPUDistort                  *
##                           process.eventWeightPUupDistort                *
##                           process.eventWeightPUdownDistort              *
##                           ## create effSF eventWeight
##                           process.effSFMuonEventWeight                  *
##                           ## multiply event weights
##                           process.eventWeightNoBtagSFWeight             *
##                           ## monitoring of PU reweighting and PV
##                           process.PUControlDistributionsDefault         *
##                           ## std sel wo b-tagging
##                           process.muonCuts                              *
##                           process.secondMuonVeto                        *
##                           process.electronVeto                          *
##                           process.leadingJetSelectionNjets4             *
##                           ## SSV btag-selection
##                           process.btagSelectionSSV                      *
##                           ## create PU event weights
##                           process.bTagSFEventWeightSSV                  *
##                           ## create combined weight
##                           process.eventWeightFinalSSV                   *
##                           ## monitoring after SSV
##                           process.SSVMonitoring                         
##                           )
## else:
##     process.p2 = cms.Path(process.dummy)

## std analysis with generator objects as input for efficiency determination
## no phase space cuts and phase space cuts for hadron level
if(runningOnData=="MC" and not cutflowSynch):
    print "running on Monte Carlo, gen-plots produced"
    process.p3 = cms.Path(## gen event selection: semileptonic (lepton &
                          ## tau->lepton if tau==True), lepton=muon/electron
                          process.genFilterSequence                     *
                          process.printTree                             *
                          ## introduce some collections
                          process.isolatedGenLeptons                    *
                          process.semiLeptGenCollections                *
                          process.cleanedGenJetCollection               *
                          process.makeGenLevelBJets                     *
                          process.noOverlapGenJetCollection             *
                          process.altermakeGenLevelBJets                *
                          process.bjetGenJets                           *
                          process.visibleHadronPS                       *
                          process.additionalGenJet                      *  
                          ## create PU event weights
                          process.makeEventWeightsPU                    *
                          ## create shape distortion event weights
                          process.eventWeightDileptonModelVariation     *
                          process.eventWeightPUDistort                  *
                          process.eventWeightPUupDistort                *
                          process.eventWeightPUdownDistort              *
                          ## investigate top reconstruction full PS
                          process.kinFitGen                             *
                          ## monitoring of hadron level kinematics
                          process.hadLvObjectMonitoring                 *
                          ## hadron level muon selection
                          process.genMuonSelection                      *
                          ## hadron level gen jet selection
                          process.genJetCuts                            *
                          ## investigate top reconstruction hadron level PS
                          process.kinFitGenPhaseSpaceHad                *
                          ## parton level phase space cuts on the basis of genTtbarEvent
                          process.filterLeptonPhaseSpace                *
                          process.filterGenPhaseSpace   
                          )
    ## delete gen filter
    if(removeGenTtbar==True):    
        process.p3.remove(process.genFilterSequence)        
        process.p3.remove(process.filterGenPhaseSpace)
        process.p3.remove(process.filterLeptonPhaseSpace)
        process.p3.remove(process.genMuonSelection)
        process.p3.remove(process.genJetCuts)
        process.p3.remove(process.noOverlapGenJetCollection)
        process.p3.remove(process.bjetGenJets)
        process.p3.remove(process.visibleHadronPS)
        process.p3.remove(process.cleanedGenJetCollection)
        process.p3.remove(process.altermakeGenLevelBJets)
        process.p3.remove(process.makeGenLevelBJets)
        process.p3.remove(process.additionalGenJet)
    if(eventFilter=='background only'):
        process.p3.remove(process.filterGenPhaseSpace)
        process.p3.remove(process.filterLeptonPhaseSpace)
        process.p3.remove(process.genMuonSelection)
        process.p3.remove(process.genJetCuts)
        process.p3.remove(process.noOverlapGenJetCollection)
        process.p3.remove(process.bjetGenJets)
        process.p3.remove(process.visibleHadronPS)
        process.p3.remove(process.cleanedGenJetCollection)
        process.p3.remove(process.altermakeGenLevelBJets)
        process.p3.remove(process.makeGenLevelBJets)
        process.p3.remove(process.additionalGenJet)
        process.p3.remove(process.dummy)
    if(not printGenParticleList):
        process.p3.remove(process.printTree)
    ## delete dummy sequence
    if(applyKinFit==False or eventFilter!="signal only"):
        process.p3.remove(process.dummy)
else:
    process.p3 = cms.Path(process.dummy)
    
## std analysis with generator objects as input for efficiency determination
## phase space cuts for parton level
if(runningOnData=="MC" and not cutflowSynch):
    process.s4 = cms.Sequence(## introduce some collections
                              process.isolatedGenLeptons                    *
                              process.semiLeptGenCollections                *
                              process.noOverlapGenJetCollection             *
                              process.altermakeGenLevelBJets                *
                              process.bjetGenJets                           *
                              process.visibleHadronPS                       *
                              process.additionalGenJet                      *
                              ## create PU event weights
                              process.makeEventWeightsPU                    *
                              ## create shape distortion event weights
                              process.eventWeightDileptonModelVariation     *
                              process.eventWeightPUDistort                  *
                              process.eventWeightPUupDistort                *
                              process.eventWeightPUdownDistort              *
			      ## new phase space cuts on the basis of genTtbarEvent
                              process.filterLeptonPhaseSpace                *  
			      process.filterGenPhaseSpace                   *
                              ## investigate top reconstruction parton level PS
                              process.kinFitGenPhaseSpace                   *
                              ## hadron level muon selection
                              process.genMuonSelection                      *
                              ## hadron level gen jet selection
                              process.genJetCuts                            
                              )
    process.p4 = cms.Path(## gen event selection: semileptonic (muon & tau->lepton)
                          ## tau->Mu if eventFilter=='background only' and
                          ## process.ttSemiLeptonicFilter.invert = True
                          process.genFilterSequence                      *
                          ## sequence with gen selection and histograms
                          process.s4
                          )			   
    ## delete gen filter
    if(removeGenTtbar==True):
        process.p4.remove(process.genFilterSequence)
        process.p4.remove(process.filterLeptonPhaseSpace)
	process.p4.remove(process.filterGenPhaseSpace)
        process.p4.remove(process.genMuonSelection)
        process.p4.remove(process.genJetCuts)
        process.p4.remove(process.noOverlapGenJetCollection)
        process.p4.remove(process.altermakeGenLevelBJets)
        process.p4.remove(process.bjetGenJets     )
        process.p4.remove(process.visibleHadronPS)
        process.p4.remove(process.additionalGenJet)
    if(eventFilter=='background only'):
        process.p4.remove(process.filterGenPhaseSpace)
        process.p4.remove(process.filterLeptonPhaseSpace)
        process.p4.remove(process.genMuonSelection)
        process.p4.remove(process.genJetCuts)
        process.p4.remove(process.noOverlapGenJetCollection)
        process.p4.remove(process.altermakeGenLevelBJets)
        process.p4.remove(process.bjetGenJets     )
        process.p4.remove(process.visibleHadronPS)
        process.p4.remove(process.additionalGenJet)
    ## delete dummy sequence
    if(applyKinFit==False or eventFilter!="signal only"):
        process.p4.remove(process.dummy)
elif(runningOnData=="data" or cutflowSynch):
    process.p4 = cms.Path(process.dummy)
    print "running on data, no gen-plots"
else:
    print "choose runningOnData= data or MC, creating no gen-plots"
    
#process.p5 = cms.Path(## gen event selection (decay channel) and the trigger selection (hltFilter)
#		  process.filterSequence                        *
#		  ## PV event selection
#		  #process.PVSelection                           *
#		  ## introduce some collections
#		  process.semiLeptonicSelection                 *
#		  ## create PU event weights
#		  process.makeEventWeightsPU                    *
#		  ## create effSF eventWeight
#		  process.effSFMuonEventWeight                  *
#		  ## multiply event weights
#		  process.eventWeightNoBtagSFWeight             *
#		  ## jet selection and monitoring
#		  process.leadingJetSelectionNjets1             *
#		  process.leadingJetSelectionNjets2             *
#		  process.leadingJetSelectionNjets3             *
#		  process.leadingJetSelectionNjets4             *
#		  ## b-tagging
#		  process.btagSelection                         *
#		  ## mod. muon selection (>0 mu with all but isolation)
#		  process.newvertexSelectedMuons                *
#		  process.newtrackMuons                         *
#		  process.testIsoMuons                          *
#		  process.testIsoMuonSelection                  *
#		  ## create PU event weights
#		  process.bTagSFEventWeight                     *
#		  ## create combined weight
#		  process.eventWeightFinal                      *
#		  process.testIsoMuonQuality
#		  )

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag

## switch to PF objects
if(jetType=="particleFlow"):
    pathlist = [process.p1, process.p2, process.p3, process.p4]#, *process.p5]
    for path in pathlist:  
        massSearchReplaceAnyInputTag(path, 'tightLeadingJets', 'tightLeadingPFJets')
        massSearchReplaceAnyInputTag(path, 'tightBottomJets' , 'tightBottomPFJets' )
        massSearchReplaceAnyInputTag(path, 'goodJets'        , 'goodJetsPF30'      )
        massSearchReplaceAnyInputTag(path, 'centralJets'     , 'centralJetsPF30'   )
        massSearchReplaceAnyInputTag(path, 'reliableJets'    , 'reliableJetsPF30'  )
        massSearchReplaceAnyInputTag(path, 'noEtaJets'       , 'noEtaJetsPF30'     )
        massSearchReplaceAnyInputTag(path, 'noPtJets'        , 'noPtJetsPF'        )
        massSearchReplaceAnyInputTag(path, 'patMETs'         , 'patMETsPF'         )
        path.remove(process.centralJets)
        path.remove(process.reliableJets)
        path.remove(process.goodJets)
        path.remove(process.tightLeadingJets)
        path.remove(process.tightBottomJets)

## switch to from muon to electron collections
if(decayChannel=="electron"):
    # adpat trigger
    if(options.mctag=="Summer11"):
      process.hltFilter.HLTPaths=["HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v*"]
      process.dummy.HLTPaths=["HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v*"]
    elif(options.mctag=="Fall11"):
      process.hltFilter.HLTPaths=["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v*"]
      process.dummy.HLTPaths=["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v*"]
    elif((options.mctag=="Summer12") or ("data" in options.sample)):
      process.hltFilter.HLTPaths=["HLT_Ele27_WP80_v*"]
      process.dummy.HLTPaths=["HLT_Ele27_WP80_v*"]
    # adapt gen filter
    process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.electron = True
    process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon= False
    
    ## lepton-jet veto - REMOVED as Top Projections take care of it
    #from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
    #from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import *
    #process.noOverlapJetsPFelec = cleanPatJets.clone(
        #src = cms.InputTag("selectedPatJetsAK5PF"),
        #preselection = cms.string(''),
        #checkOverlaps = cms.PSet(
          #electrons = cms.PSet(
            #src       = cms.InputTag("tightElectronsEJ"),
            #algorithm = cms.string("byDeltaR"),
            #preselection        = cms.string(''),
            #deltaR              = cms.double(0.3),
            #checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
            #pairCut             = cms.string(""),
            #requireNoOverlaps   = cms.bool(True), # overlaps don't cause the jet to be discared
            #)
          #),
        #finalCut = cms.string(''),
        #)
    #process.goodJetsPF20.src  ='noOverlapJetsPFelec'
    #process.centralJetsPF.src ='noOverlapJetsPFelec'
    #process.reliableJetsPF.src='noOverlapJetsPFelec'
    #process.noEtaJetsPF.src   ='noOverlapJetsPFelec'
    #process.noPtJetsPF.src    ='noOverlapJetsPFelec'
    #process.noConstJetsPF.src ='noOverlapJetsPFelec'
    #process.noCEFJetsPF.src   ='noOverlapJetsPFelec'
    #process.noNHFJetsPF.src   ='noOverlapJetsPFelec'
    #process.noNEFJetsPF.src   ='noOverlapJetsPFelec'
    #process.noCHFJetsPF.src   ='noOverlapJetsPFelec'
    #process.noNCHJetsPF.src   ='noOverlapJetsPFelec'
    #process.noKinJetsPF.src   ='noOverlapJetsPFelec' 
    
    # gen selection
    process.p3.replace(process.genMuonSelection, process.genElectronSelection)
    process.p4.replace(process.genMuonSelection, process.genElectronSelection)
    process.cleanedGenJetCollection.checkOverlaps.muons.src=cms.InputTag("selectedGenElectronCollection")
    pathlist = [process.p1, process.p2, process.p3, process.p4]#, process.p5]
    for path in pathlist:
        # replace jet lepton veto
        #path.replace(process.noOverlapJetsPF, process.noOverlapJetsPFelec)
        # replace muon selection
        path.remove(process.muonCuts)
        path.remove(process.secondMuonVeto)
        path.replace( process.electronVeto, process.electronSelection)
	## replace effSF
        path.replace(process.effSFMuonEventWeight, process.effSFElectronEventWeight)
        ## replace gen object kinematics
        if(runningOnData=="MC"):
            process.genAllElectronKinematics.weight=cms.InputTag("")
        # remove muon monitoring
        path.remove(process.tightMuontightJetsKinematics)
        path.remove(process.tightMuonKinematics)
        path.remove(process.tightMuonQuality)
        path.remove(process.tightMuontightJetsKinematicsTagged)
        path.remove(process.tightMuonKinematicsTagged         )
        path.remove(process.tightMuonQualityTagged            )
        path.remove(process.tightMuontightJetsKinematicsSSV)
        # add electron monitoring
        path.replace(process.tightLead_0_JetKinematics      , process.tightElectronKinematics * process.tightElectronQuality * process.tightLead_0_JetKinematics)
        path.replace(process.tightLead_0_JetKinematicsTagged, process.tightElectronKinematicsTagged * process.tightElectronQualityTagged  * process.tightLead_0_JetKinematicsTagged)
        path.replace(process.tightMuonKinematicsSSV, process.tightElectronKinematicsSSV)
        path.replace(process.tightMuonQualitySSV   , process.tightElectronQualitySSV   )
        path.replace(process.newvertexSelectedMuons , process.newvertexSelectedElectrons)
        path.remove(process.newtrackMuons)
        path.replace(process.testIsoMuons, process.testIsoElectrons)
        path.replace(process.testIsoMuonSelection, process.testIsoElectronSelection)
        path.replace(process.testIsoMuonQuality, process.testIsoElectronQuality)
        # replace muon by electron in (remaining) kinfit analyzers
        massSearchReplaceAnyInputTag(path, 'tightMuons', 'goodElectronsEJ')
        # take care of replacements you do NOT want to do!
        process.compositedKinematics.MuonSrc='tightMuons'
        process.compositedKinematicsKinFit.MuonSrc='tightMuons'
        process.compositedKinematicsProbSel.MuonSrc='tightMuons'
        #if(eventFilter=='signal only'):
        process.compositedPartonGen.MuonSrc='tightMuons'
        process.compositedPartonGenPhaseSpace.MuonSrc='tightMuons'
        process.compositedHadronGenPhaseSpace.MuonSrc='tightMuons'

allpaths  = process.paths_().keys()

# switch to PF2PAT
if(pfToPAT):
    from TopAnalysis.TopUtils.usePatTupleWithParticleFlow_cff import prependPF2PATSequence
    if(not cutflowSynch):
        recoPaths=['p1','p2']#,'p5']
    else:
        recoPaths=['p1']
    # define general options
    PFoptions = {
        'runOnMC': True,
        'runOnAOD': True,
        'switchOffEmbedding': False,
        'addResolutions': True,
        'resolutionsVersion': 'fall11',
        'runOnOLDcfg': True,
        'cutsMuon': 'pt > 10. & abs(eta) < 2.5',
        'cutsElec': 'et > 20. & abs(eta) < 2.5',
        'cutsJets': 'pt > 10 & abs(eta) < 5.0', 
        'electronIDs': ['CiC','classical','MVA'],
        'pfIsoConeMuon': 0.4,
        'pfIsoConeElec': 0.3,
        'pfIsoValMuon': 0.2,
        'pfIsoValElec': 0.15,
        'doDeltaBetaCorrMuon' : True,
        'doDeltaBetaCorrElec' : True,
        'skipIfNoPFMuon': False,
        'skipIfNoPFElec': False,
        'addNoCutPFMuon': False,
        'addNoCutPFElec': False,
        'noMuonTopProjection': False,
        'noElecTopProjection': False,
        'analyzersBeforeMuonIso':cms.Sequence(),
        'analyzersBeforeElecIso':cms.Sequence(),
        'excludeElectronsFromWsFromGenJets': True,
        'METCorrectionLevel': 1
        }
    # adaptions when running on data
    if(runningOnData=="data"):
        PFoptions['runOnMC']=False
    if(decayChannel=="electron"):
    # skip events (and jet calculation) if no lepton is found
    # only done in data, as in MC you need the events for parton truth plots
        #PFoptions['skipIfNoPFElec']=True
    # collection without cuts is added
        PFoptions['addNoCutPFElec']=True
    # project no other leptons than the selected ones
        #PFoptions['noMuonTopProjection']=True
    elif(decayChannel=="muon"):
        #PFoptions['skipIfNoPFMuon']=True
        PFoptions['addNoCutPFMuon']=True
        #PFoptions['noElecTopProjection']=True
	## option to change PF2PAT settings for cutflow excercise:
	if(cutflowSynch):
	    PFoptions['pfIsoValMuon'] = 0.2
	    #PFoptions['cutsElec'] = 'pt>5 && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits<2'
    prependPF2PATSequence(process, recoPaths, PFoptions)

    # remove electron collections as long as id does not exist in the tuples
    for path in recoPaths:
        #getattr(process,path).remove( process.looseElectronsEJ )
        #getattr(process,path).remove( process.tightElectronsEJ )
        #getattr(process,path).remove( process.unconvTightElectronsEJ )
        #getattr(process,path).remove( process.goodElectronsEJ )
        # replace object consistently with names from PF2PAT
        #massSearchReplaceAnyInputTag(getattr(process,path), 'patMETsPF', 'patMETs')
        #massSearchReplaceAnyInputTag(getattr(process,path), 'selectedPatJetsAK5PF', 'selectedPatJets')        
        # run trigger at the beginning to save a lot of time
        getattr(process,path).insert(0,process.hltFilter)
        # investigate specific events
        #if(cutflowSynch):
        #    if(decayChannel=="electron"):
        #        getattr(process,path).insert(0,~process.multiEventFilterElec)
        #    elif(decayChannel=="muon"):
        #        getattr(process,path).insert(0,~process.multiEventFilterMuon)
## change decay subset to parton level (ME)
#process.decaySubset.fillMode = cms.string("kME")

## generator kinematics (check E-p conservation on gen level)
if (PythiaSample=="False" and runningOnData=="MC"):
    if (options.mctag=='Fall11' or options.mctag=='Summer11'):
        print
        print "Applying MC Kinematics Filter"
        for path in allpaths:
            getattr(process,path).insert(0,process.totalKinematicsFilterDefault)

# include MC event weights (needed for MCatNLO)
if(MCweighting):
    MCpathlist = [process.p1, process.p2, process.p3, process.p4]#, *process.p5]
    for path in MCpathlist:
        path.insert(0,process.eventWeightMC)
        path.replace(process.makeEventWeightsPU,
                     process.makeEventWeightsPU  * process.combineEventWeightsMCandPU)

## Output Module Configuration
if(writeOutput):
    from PhysicsTools.PatAlgos.patEventContent_cff import *
    process.out = cms.OutputModule("PoolOutputModule",
                                   fileName = cms.untracked.string('patTuple_selectedNjets4.root'),
                                   # save only events passing the full path
                                   SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p1') ),
                                   #save output (comment to keep everything...)
                                   outputCommands = cms.untracked.vstring('drop *') 
                                   )
    process.out.outputCommands += patEventContentNoCleaning
    process.out.outputCommands += patExtraAodEventContent
    process.outpath = cms.EndPath(process.out)


## possibly remove event reweighting

# define allpaths if not done yet
if(not pfToPAT):
    allpaths  = process.paths_().keys()
    
# Pile up
if(not PUreweigthing or runningOnData=="data"):
    for path in allpaths:
        getattr(process,path).remove( process.makeEventWeightsPU )
        if(MCweighting):
            getattr(process,path).remove( process.combineEventWeightsMCandPU )         
# Eff SF
if(not effSFReweigthing or runningOnData=="data"):
    for path in allpaths:
        if(decayChannel=="muon"):
            getattr(process,path).remove( process.effSFMuonEventWeight )
	elif(decayChannel=="electron"):
            getattr(process,path).remove( process.effSFElectronEventWeight )

# Btag scale factor
if(not BtagReweigthing or runningOnData=="data"):
    for path in allpaths:
        getattr(process,path).remove( process.bTagSFEventWeight )

# combined scale factor
if(runningOnData=="data" or (not PUreweigthing and not effSFReweigthing) ):
    for path in allpaths:
        getattr(process,path).remove( process.eventWeightNoBtagSFWeight )

# combined scale factor   
if(runningOnData=="data" or (not BtagReweigthing and not effSFReweigthing) ):
    for path in allpaths:
        getattr(process,path).remove( process.eventWeightNoPUWeight )

# combined scale factor
if(runningOnData=="data" or (not PUreweigthing and not BtagReweigthing and not effSFReweigthing) ):
    for path in allpaths:
        getattr(process,path).remove( process.eventWeightFinal )

# shape distortion scale factor
if(sysDistort==''):
    for path in allpaths:
        getattr(process,path).remove( process.eventWeightDileptonModelVariation )
        getattr(process,path).remove( process.eventWeightPUDistort              )
        getattr(process,path).remove( process.eventWeightPUupDistort            )
        getattr(process,path).remove( process.eventWeightPUdownDistort          )

# remove reco paths
if(genFull):
    process.p1=cms.Path(process.dummy)
    process.p2=cms.Path(process.dummy)
        
## use ecal driven momentum for electron collection
process.load("TopAnalysis.TopUtils.JetEnergyScale_cfi")
from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
for path in pathlist:
    massSearchReplaceAnyInputTag(path, 'selectedPatElectrons', 'ecalMomentum:selectedPatElectrons')
process.ecalMomentum=process.scaledJetEnergy.clone()
process.ecalMomentum.scaleFactor         = cms.double(1.0)
process.ecalMomentum.scaleFactorB        = cms.double(1.0)
process.ecalMomentum.scaleType           = cms.string("abs")
process.ecalMomentum.resolutionFactors   = cms.vdouble(1.0)
process.ecalMomentum.resolutionEtaRanges = cms.vdouble(0, -1)
process.ecalMomentum.inputElectrons      = cms.InputTag("selectedPatElectrons")
process.ecalMomentum.inputJets           = cms.InputTag("selectedPatJetsAK5PF")
for path in pathlist:
    path.replace(process.pf2pat,
                 process.pf2pat * process.ecalMomentum)
if(runningOnData=="data"):
    process.ecalMomentum.inputJets = cms.InputTag("selectedPatJets")

# reduce module content if reduced==True
if(reduced):
    for path in allpaths:
        #if(hasattr(process,                      'eventWeightNoPUWeight')):  # still needed for PU control distribution analyzer!!!
        #    getattr(process,path).remove( process.eventWeightNoPUWeight    )
        if(hasattr(process,                      'eventWeightFlatEffSF')):
            getattr(process,path).remove( process.eventWeightFlatEffSF     )
        if(hasattr(process,                      'analyzeTopRecoKinematicsKinFitNoWeight')):
            getattr(process,path).remove( process.analyzeTopRecoKinematicsKinFitNoWeight      )
        if(hasattr(process,                      'analyzeTopRecoKinematicsKinFitOnlyPUWeight')):
            getattr(process,path).remove( process.analyzeTopRecoKinematicsKinFitOnlyPUWeight  )
        if(hasattr(process,                      'analyzeTopRecoKinematicsKinFitNoBtagSFWeight')):
            getattr(process,path).remove( process.analyzeTopRecoKinematicsKinFitNoBtagSFWeight)
        if(hasattr(process,                      'analyzeTopRecoKinematicsKinFitFlatEffSF')):
            getattr(process,path).remove( process.analyzeTopRecoKinematicsKinFitFlatEffSF  )
        if(hasattr(process,                      'analyzeTopRecoKinematicsBjetsFlatEffSF')):
            getattr(process,path).remove( process.analyzeTopRecoKinematicsBjetsFlatEffSF   )
        if(hasattr(process,                      'analyzeTopRecoKinematicsLeptonFlatEffSF')):
            getattr(process,path).remove( process.analyzeTopRecoKinematicsLeptonFlatEffSF  )
        if(hasattr(process,                      'analyzeTopPartonLevelKinematicsNoPUWeight')):
            getattr(process,path).remove( process.analyzeTopPartonLevelKinematicsNoPUWeight      )
        if(hasattr(process,                      'analyzeTopPartonLevelKinematicsBjetsNoPUWeight')):
            getattr(process,path).remove( process.analyzeTopPartonLevelKinematicsBjetsNoPUWeight )
        if(hasattr(process,                      'analyzeTopPartonLevelKinematicsLeptonNoPUWeight')):
            getattr(process,path).remove( process.analyzeTopPartonLevelKinematicsLeptonNoPUWeight)
        if(hasattr(process,                      'analyzeTopHadronLevelKinematicsPhaseSpaceNoPUWeight')):
            getattr(process,path).remove( process.analyzeTopHadronLevelKinematicsPhaseSpaceNoPUWeight      )
        if(hasattr(process,                      'analyzeTopHadronLevelKinematicsBjetsPhaseSpaceNoPUWeight')):
            getattr(process,path).remove( process.analyzeTopHadronLevelKinematicsBjetsPhaseSpaceNoPUWeight )
        if(hasattr(process,                      'analyzeTopHadronLevelKinematicsLeptonPhaseSpaceNoPUWeight')):
            getattr(process,path).remove( process.analyzeTopHadronLevelKinematicsLeptonPhaseSpaceNoPUWeight)
        if(hasattr(process,                      'analyzeTopPartonLevelKinematicsBjetsPhaseSpaceNoPUWeight')):
            getattr(process,path).remove( process.analyzeTopPartonLevelKinematicsBjetsPhaseSpaceNoPUWeight )
        if(hasattr(process,                      'analyzeTopPartonLevelKinematicsLeptonPhaseSpaceNoPUWeight')):
            getattr(process,path).remove( process.analyzeTopPartonLevelKinematicsLeptonPhaseSpaceNoPUWeight)
        if(hasattr(process,                      'tightMuonKinematicsNjets1')):
            getattr(process,path).remove( process.tightMuonKinematicsNjets1)
        if(hasattr(process,                      'tightMuonQualityNjets1')):
            getattr(process,path).remove( process.tightMuonQualityNjets1   )
        if(hasattr(process,                      'tightMuonKinematicsNjets2')):
            getattr(process,path).remove( process.tightMuonKinematicsNjets2)
        if(hasattr(process,                      'tightMuonQualityNjets2')):
            getattr(process,path).remove( process.tightMuonQualityNjets2   )
        if(hasattr(process,                      'tightMuonKinematicsNjets3')):
            getattr(process,path).remove( process.tightMuonKinematicsNjets3)
        if(hasattr(process,                      'tightMuonQualityNjets3')):
            getattr(process,path).remove( process.tightMuonQualityNjets3   )
        if(hasattr(process,                      'tightJetKinematicsNjets1')):
            getattr(process,path).remove( process.tightJetKinematicsNjets1)
        if(hasattr(process,                      'tightJetQualityNjets1')):
            getattr(process,path).remove( process.tightJetQualityNjets1   )
        if(hasattr(process,                      'tightJetKinematicsNjets2')):
            getattr(process,path).remove( process.tightJetKinematicsNjets2)
        if(hasattr(process,                      'tightJetQualityNjets2')):
            getattr(process,path).remove( process.tightJetQualityNjets2   )
        if(hasattr(process,                      'tightJetKinematicsNjets3')):
            getattr(process,path).remove( process.tightJetKinematicsNjets3)
        if(hasattr(process,                      'tightJetQualityNjets3')):
            getattr(process,path).remove( process.tightJetQualityNjets3   )
        if(hasattr(process,                      'tightElectronKinematicsNjets1')):
            getattr(process,path).remove( process.tightElectronKinematicsNjets1)
        if(hasattr(process,                      'tightElectronQualityNjets1')):
            getattr(process,path).remove( process.tightElectronQualityNjets1   )
        if(hasattr(process,                      'tightElectronKinematicsNjets2')):
            getattr(process,path).remove( process.tightElectronKinematicsNjets2)
        if(hasattr(process,                      'tightElectronQualityNjets2')):
            getattr(process,path).remove( process.tightElectronQualityNjets2   )
        if(hasattr(process,                      'tightElectronKinematicsNjets3')):
            getattr(process,path).remove( process.tightElectronKinematicsNjets3)
        if(hasattr(process,                      'tightElectronQualityNjets3')):
            getattr(process,path).remove( process.tightElectronQualityNjets3   )
            
if(runningOnData=="MC"):
    ## change input collections to JER-shifted collections (for reco paths only)
    #from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    pathlist = [process.p1, process.p2]#, process.p5]
    for path in pathlist:
        if(jetType=="particleFlow"):
            massSearchReplaceAnyInputTag(path, 'selectedPatJetsAK5PF', 'scaledJetEnergy:selectedPatJetsAK5PF')
            massSearchReplaceAnyInputTag(path, 'patMETsPF'           , 'scaledJetEnergy:patMETsPF'           )
        elif(jetType=="Calo"):
            massSearchReplaceAnyInputTag(path, 'selectedPatJets', 'scaledJetEnergy:selectedPatJets') 
            massSearchReplaceAnyInputTag(path, 'patMETs'        , 'scaledJetEnergy:patMETs'        )
        if(pfToPAT==True):
            massSearchReplaceAnyInputTag(path, 'selectedPatJets', 'scaledJetEnergy:selectedPatJets')
            massSearchReplaceAnyInputTag(path, 'patMETs'        , 'scaledJetEnergy:patMETs'        )
        else:
            print "unknown jetType"
    # keep uncorrected jets to avoid out of range error in this module
    process.ecalMomentum.inputJets = "selectedPatJets"
    # eta-dependent smearing of the jet energy
    process.scaledJetEnergy.resolutionFactors   = cms.vdouble( 1.052 , 1.057 , 1.096 , 1.134 , 1.288 )
    process.scaledJetEnergy.resolutionEtaRanges = cms.vdouble(0.0,0.5,0.5,1.1,1.1,1.7,1.7,2.3,2.3,-1.)
    if(options.JEtag=="JERUp"):
        process.scaledJetEnergy.resolutionFactors   = cms.vdouble( 1.115 , 1.114 , 1.161 , 1.228 , 1.488 )
    elif(options.JEtag=="JERDown"):
        process.scaledJetEnergy.resolutionFactors   = cms.vdouble( 0.990 , 1.001 , 1.032 , 1.042 , 1.089 )

    # up and down scaling of the jet energy
    if(jetType=="particleFlow"):
        process.scaledJetEnergy.inputJets = "selectedPatJetsAK5PF"
        process.scaledJetEnergy.inputMETs = "patMETsPF"
        process.scaledJetEnergy.payload   = "AK5PF"
    elif(jetType=="Calo"):
        process.scaledJetEnergy.inputJets = "selectedPatJets"
        process.scaledJetEnergy.inputMETs = "patMETs"
        process.scaledJetEnergy.payload   = "AK5Calo"
    else:
        print "unknown jetType"
    if(pfToPAT==True):
        process.scaledJetEnergy.inputJets = "selectedPatJets"
        process.scaledJetEnergy.inputMETs = "patMETs"

    process.scaledJetEnergy.JECUncSrcFile = "TopAnalysis/TopUtils/data/Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt"
    if(options.JEtag=="JESUp"):
        process.scaledJetEnergy.scaleType   = "jes:up"
    elif(options.JEtag=="JESDown"):
        process.scaledJetEnergy.scaleType   = "jes:down"

    ## include module to create JES-shifted collection
    for path in pathlist:
        path.replace(process.selectedPatJets,
                     process.selectedPatJets * process.scaledJetEnergy)

    if(applyKinFit==True):
        # use status 3 particles (!)
        process.decaySubset.fillMode = cms.string("kME")


## adjust lepton pT cut
import re
ptval="33" # adjust cut value here! [GeV]
relevantLeptonCollections = [process.tightElectronsEJ, process.goldenMuons, process.selectedGenMuonCollection, process.selectedGenElectronCollection]
exp = re.compile('(?:t\s?>\s?30)') 
for lep in relevantLeptonCollections:
    if(exp.search(lep.cut.pythonValue())!=None):
        tmpExp=exp.sub("t > "+ptval, str(lep.cut.pythonValue()))
        lep.cut=tmpExp
        if(tmpExp.find("\'")>-1):
            lep.cut=tmpExp.strip("'")

