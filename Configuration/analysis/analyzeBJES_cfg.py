import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os
import sys
options = VarParsing.VarParsing ('standard')

options.register('mcversion', 'Summer12', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "MC campaign or data")
options.register('mcWeight', 1.0 , VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.float, "MC sample event weight")
options.register('lepton', 'muon', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Lepton+jets channel")
options.register('metcl', 1, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "MET correction level")

options.register('scaleType', 'abs', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "JES scale type")
options.register('jessource', '', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Uncertainty source for JES variation for source:up/down")
options.register('lJesFactor', 1.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "JES")
options.register('bJesFactor', 1.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "bJES")
options.register('resolution', 'nominal', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "JER")
options.register('mesShift', 0.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "Muon energy scale shift in sigma")
options.register('eesShift', 0.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "Electron energy scale shift in sigma")
options.register('uncFactor', 1.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "Unclustered energy factor")

options.register('csvm', 0.679, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "CSVM working point")
options.register('nbjets', 2, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "Minimum number of bjets")

options.register('brCorrection', True, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "Do BR correction (MadGraph)")

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

print options

if (options.mcversion == "data"): data = True
else:                             data = False

## top mass measurement
process = cms.Process("bJES")

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.categories.append('TtSemiLeptonicEvent')
process.MessageLogger.cerr.TtSemiLeptonicEvent = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)

## define input
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles,
                              dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
)

if (options.mcversion == "Summer12GEN"):
  readFiles.extend( [ 'file:TOP-Summer12-00190_step1.root' ] )
elif (options.mcversion == "Sherpa12"):
  readFiles.extend( [
         '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_0_fastreco_FASTSIM_HLT_PU.root',
         '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_2_fastreco_FASTSIM_HLT_PU.root',
         '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_4_fastreco_FASTSIM_HLT_PU.root',
         '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_6_fastreco_FASTSIM_HLT_PU.root',
         '/store/user/mseidel/TT_cluster_8TeV-sherpa/job_8_fastreco_FASTSIM_HLT_PU.root',
         '/store/user/mseidel/TT_lund_8TeV-sherpa/job_1_fastreco_FASTSIM_HLT_PU.root',
         '/store/user/mseidel/TT_lund_8TeV-sherpa/job_3_fastreco_FASTSIM_HLT_PU.root',
         '/store/user/mseidel/TT_lund_8TeV-sherpa/job_5_fastreco_FASTSIM_HLT_PU.root',
         '/store/user/mseidel/TT_lund_8TeV-sherpa/job_7_fastreco_FASTSIM_HLT_PU.root',
         '/store/user/mseidel/TT_lund_8TeV-sherpa/job_9_fastreco_FASTSIM_HLT_PU.root',
         '/store/user/mseidel/TT_lund_8TeV-sherpa/job_13_fastreco_FASTSIM_HLT_PU.root',
         
  ] )
elif (options.mcversion == "Summer12"):
  readFiles.extend( [
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/F4BE4558-BE42-E311-8359-7845C4F9321B.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/74171050-BF42-E311-9AE0-7845C4FC3A1C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/0E010BFF-1D43-E311-8CBE-848F69FD287A.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/7609E0BD-1E43-E311-8A6C-00A0D1EE9274.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/40632113-CE42-E311-9580-7845C4FC375E.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/1E4FA40B-1F43-E311-A539-00266CF25878.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/8EC57142-EE42-E311-B789-00A0D1EEE660.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/DA4BA806-EF42-E311-A0B2-848F69FD288F.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/AC384B68-EF42-E311-AC76-008CFA00879C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/CA199456-F042-E311-BF97-008CFA000744.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B888280F-F142-E311-866E-7845C4FC3983.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/F41C8366-F142-E311-BAF4-00266CF25B24.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/E4432BB2-0943-E311-9B90-008CFA001EE4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/BED49457-1743-E311-861D-00266CFAE4A0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/20FF45D4-2A43-E311-B50F-848F69FD2988.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/60D5B384-BA42-E311-8FB6-00266CF9BCAC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D872F65D-C042-E311-97A2-00A0D1EE8B08.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/86BE67E9-0043-E311-969F-00A0D1EEE3B0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/42A618CB-0143-E311-9E61-00266CFAE468.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/BA5B6F83-0143-E311-A9E0-00A0D1EE2A50.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/4447F70B-0343-E311-B2F1-008CFA001B18.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/8AD92272-0243-E311-81B4-848F69FD28C2.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/660B1F76-0D43-E311-B069-008CFA001334.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/5E18D20E-0D43-E311-A9FF-00266CF25B24.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/7689EB6B-0D43-E311-B78E-00266CF9B254.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/121B1FB5-0D43-E311-811E-7845C4FC39C5.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/341AD85C-0E43-E311-9CE7-008CFA001334.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/E84D8AD1-0E43-E311-93E8-00A0D1EE95CC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/98D59F3E-1843-E311-84BB-00266CF9B04C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/5A750928-1843-E311-902D-848F69FD47A5.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/0A13566A-1843-E311-9DA6-00266CFAE464.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/AC9B6868-1943-E311-8AA9-848F69FD4C9A.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/58B74B1A-2643-E311-A357-00266CFAE788.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/50B86656-2643-E311-8803-00266CFAE318.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/7C53F186-2643-E311-90EB-00266CFAE6A4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/9CB437B0-2643-E311-9473-00266CF25C20.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/20A61701-3543-E311-BE09-008CFA002028.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/402D3753-3543-E311-9DBA-00266CF9B1C4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/161CB867-3543-E311-BD2D-00A0D1EE95FC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A2D8DABF-3543-E311-9F5A-00A0D1EE2A50.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/4ECFFCE1-3543-E311-A69D-7845C4FC3A67.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/186034D2-3543-E311-98BC-00A0D1EE8C64.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/CA90D5DA-3543-E311-916D-00266CFADE34.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/825400FD-3543-E311-87E5-00266CFAED7C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/06AA1CF3-3543-E311-BF8F-180373FF93CE.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/CC3DDCF7-3543-E311-90C1-848F69FD5027.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/0AE26C0C-3643-E311-B258-00266CF25D18.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/42218B5D-3643-E311-90CE-848F69FD28C8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/68117468-3643-E311-8B3A-00266CF2506C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/226FAB55-3643-E311-A294-00266CF25DC4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/00F72AF2-0443-E311-A2D7-7845C4FC362F.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/702A8328-0F43-E311-B543-00266CF25218.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/CA8CA892-0E43-E311-AD50-180373FF8446.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/26F373A4-0F43-E311-9796-848F69FD2A30.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/92D3C7EF-0F43-E311-A4AD-00A0D1EE95CC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/BE76F30D-1043-E311-A9E0-008CFA001334.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/5058EDD8-1843-E311-8E6D-00266CF27130.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/E630D75D-3643-E311-AA1B-848F69FD294F.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D2A96353-3643-E311-9038-848F69FD2A27.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/FC1D717A-3643-E311-93F3-180373FF8CD4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/5E6C2034-2743-E311-951C-00266CFAE074.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/0CD7B7A1-2743-E311-95D6-848F69FD2931.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/467F656B-2843-E311-8BD1-00266CF9B1C4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/BA0A3F98-2743-E311-95C7-848F69FD2928.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/5C4FCF01-1A43-E311-AA71-848F69FD47A5.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/56F20362-1A43-E311-9F38-848F69FD2817.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/584D7F45-4043-E311-B9B8-00A0D1EE95AC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/8605F860-BC42-E311-BD64-7845C4FC3983.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D49FFA1E-4043-E311-A004-008CFA0014EC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/FCDED6AD-4043-E311-94E2-00266CFAED7C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/E0B90FB1-4043-E311-8804-7845C4FC3641.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/704BC1B0-4043-E311-8C45-7845C4FC39FB.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/3836DFB0-B642-E311-A049-00266CFAE740.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D44B55C9-1B43-E311-A405-00266CF252D4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/7EBB0FB2-4043-E311-8E49-848F69FD24D2.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/68FF02CD-4043-E311-90FD-7845C4FC39CB.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/FA52EE92-4043-E311-B5E7-00A0D1EE89E0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/96F4FFBB-4043-E311-9641-00266CF9ADA0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/08D5724B-4143-E311-AA52-00266CF9B878.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/724CA04B-4143-E311-8DBD-00266CFAEBA0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A6D43027-4143-E311-A006-848F69FD0884.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/3EB74436-F942-E311-A04C-7845C4FC346A.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D4DB47F2-0643-E311-AE15-008CFA00206C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/56B25CA2-0743-E311-8DFA-848F69FD2955.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/621E7B24-DB42-E311-A015-00266CF25DC4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/724C1B0F-F942-E311-A238-00266CFAEBF8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/6418D6DB-F942-E311-9D2E-7845C4F92EC7.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/32177EE5-F842-E311-B156-7845C4FC3A61.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/E8E9A8AA-0743-E311-9FDB-008CFA001FDC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/E87E57F3-0743-E311-8E65-008CFA001EE4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/02B81357-0843-E311-8F5E-001D09FDD7D4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/1021E6DA-1C43-E311-B12A-848F69FD29D3.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/9ADF421C-1D43-E311-8F3D-00266CF9B1C4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B2299732-1E43-E311-9F6A-001D09FDD7D4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/4A5B13F2-1D43-E311-9C11-7845C4FC397A.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/64E6122D-CB42-E311-9182-008CFA00018C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B0ADB255-CC42-E311-9B17-848F69FD46BB.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/90EB7662-D342-E311-8CC4-00266CF25DC4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/3421BEC2-D442-E311-9890-008CFA00210C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/68A440E2-1143-E311-96DB-008CFA001334.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/885EC857-1243-E311-8E14-00A0D1EE8AF0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/30E51962-0043-E311-8011-00266CFAE468.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/303276EA-0043-E311-8871-00A0D1EE29B8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/AADBFD63-1343-E311-AFAE-7845C4F92ECD.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/48C8DFAB-1343-E311-A5B6-00266CF9B5D0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/E2BE6B6C-1343-E311-AF36-008CFA000744.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/4416C4FF-1343-E311-90F7-00A0D1EE2A50.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/66B55BA2-0943-E311-B40E-008CFA001FDC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/DA818EEC-E242-E311-9A51-00266CF9B404.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/26F928CF-E342-E311-A5C4-848F69FD28E3.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A0EB0B02-FD42-E311-BF7B-7845C4FC3A61.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A82CC59F-FD42-E311-921C-008CFA008D4C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/5C57E4ED-FC42-E311-B6E7-00266CF9BEF8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/74ACCC43-0B43-E311-9593-848F69FD44B1.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A8A38A5F-0B43-E311-AE09-00266CF25B24.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/94BDBEA5-2B43-E311-BC72-00A0D1EEAA00.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/2A546D56-2C43-E311-88C8-00266CFAE6E0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/6E36C253-2C43-E311-B36B-008CFA008D0C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A07D8EB8-2C43-E311-826D-00266CF9B1C4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D8429D36-0C43-E311-83A8-00A0D1EE2990.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A249E2F3-0B43-E311-A51A-008CFA002E80.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A6247E28-0C43-E311-9D2F-7845C4FC39C5.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/46144784-2F43-E311-888D-848F69FD4586.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B2FF9029-2F43-E311-9FD8-848F69FD4568.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/5491A7A0-3143-E311-8203-00266CF9B034.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/1E99689A-3143-E311-9ACF-848F69FD4FC1.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/420A3AC7-1B43-E311-9EA8-7845C4F9162D.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/7ADB5BE6-1C43-E311-83F8-7845C4F91495.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/8894544F-1C43-E311-B7D3-008CFA000744.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/025537DD-1B43-E311-BAA5-008CFA0016A4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/5EA0A5CC-1B43-E311-B4A4-848F69FD4C94.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B0AF5F55-4A43-E311-B779-00266CF9AFF0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/F04CE795-4A43-E311-939F-7845C4FC3B6F.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/565FD670-4943-E311-A9F1-848F69FD2853.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/FA2B5634-4A43-E311-B43C-7845C4F932D8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/AAC88B93-4A43-E311-A78B-00266CF279F8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/8A5117FD-C442-E311-AA92-848F69FD29CA.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/5C40A70A-C642-E311-94E0-7845C4FC3998.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/ECA2FEA8-1443-E311-BF23-00A0D1EE8E64.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/7A0789AE-1443-E311-91B2-00A0D1EE95AC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/AE840A4C-1C43-E311-A1E5-7845C4FC3653.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/84765FDF-1B43-E311-B1C5-7845C4FC3B00.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/38EFC3DC-4643-E311-8909-008CFA008F50.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/AC2BD9DA-4643-E311-9302-008CFA00038C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/BAE5789A-4743-E311-B4C1-008CFA002E80.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/84B704AF-4743-E311-B3F7-00A0D1EC3950.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/BADBEF5F-FC42-E311-96EB-7845C4FC39C5.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D082F39D-FD42-E311-B539-848F69FD2952.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B83AAA0F-FE42-E311-924E-00A0D1EE8D7C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B62D5D5A-0543-E311-8DD7-00266CF25C88.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/FC832A50-0543-E311-9B21-008CFA00206C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/06AA3A27-2A43-E311-8510-008CFA007F18.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/5C11C746-2A43-E311-A7DC-00A0D1EE9238.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/08C72E24-2E43-E311-AFC2-848F69FD29B2.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/2C785A3C-2E43-E311-A2CB-00A0D1EEE3B0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/7A897E54-2E43-E311-8399-848F69FD28E3.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/E8B7E652-1243-E311-8ED6-848F69FD2A30.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/9A07057E-1243-E311-989A-00266CF279F8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D6DE9A7E-1343-E311-BAB5-7845C4FC3AE5.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D4FF582E-1343-E311-963E-00A0D1EE8E64.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/1CBDEA44-F742-E311-9C9D-008CFA007F18.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/2E6304CA-F542-E311-A2AE-00266CFAEBF8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/62BCEBB2-F742-E311-AEE4-00266CFACC38.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/2E93E772-1643-E311-BEC0-00266CFAE464.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D437D2CE-1643-E311-A69A-00A0D1EE8E78.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/2E357871-FE42-E311-8E8A-00266CFAE818.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/4CEB73F2-FD42-E311-95B1-7845C4FC39C8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/2C1D297B-FE42-E311-B848-00A0D1EE8F40.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/DE9DA89B-C742-E311-B657-7845C4F91633.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/EA3D4895-2243-E311-9CDB-00266CFAE6A4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/2A140EF9-2343-E311-80F9-7845C4FC3C11.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/40885F1D-2443-E311-93A1-7845C4FC3B00.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/6EC3DB19-2443-E311-995A-7845C4FC3A7F.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/C6B1416D-CF42-E311-B8ED-7845C4FC39FB.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/EAAEAB12-D142-E311-A154-008CFA001334.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/CED67F13-D242-E311-9D8F-00A0D1EE8E14.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D45614C4-0143-E311-A0B0-00266CF9AED8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/9E471201-2543-E311-8B66-00266CFAE6A4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/4A560E05-2443-E311-9063-848F69FD2823.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/3E578818-2443-E311-BFFE-00266CF9AFF0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/DC724A03-0443-E311-8050-848F69FD2892.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/00A0ED45-0443-E311-9EC4-848F69FD2952.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/3E2F034F-0443-E311-9C35-180373FFCE1E.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/F6188549-0443-E311-8256-848F69FD2955.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/DC08EE58-F242-E311-A4B4-00266CF20468.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/BC706152-F242-E311-8913-00A0D1EEAA00.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/6C38A248-F342-E311-946B-00266CF9B884.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/120BD795-F342-E311-A36B-00266CF20468.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/00BB7AF4-F342-E311-A83A-00266CF9B868.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/BA285EDE-F442-E311-BADC-008CFA0025E8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/70ABA580-F542-E311-BE41-00266CFAC810.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/DE4A2A75-F642-E311-81AD-008CFA0025E8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/0C6BE286-1A43-E311-AF92-00266CF27130.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/EA5BAB43-1B43-E311-8236-848F69FD4FC1.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/90927D3F-1B43-E311-A4C4-00266CF9B1C4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/227E7D8A-2F43-E311-9EEC-848F69FD28E3.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/16AA081A-3B43-E311-B8D2-00266CFADEC0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/26A7A02E-3B43-E311-BCEF-00266CFADFCC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/0A429F74-3B43-E311-9E22-008CFA001FDC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/72CE987C-3B43-E311-9BEB-848F69FD3FAA.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/545DBF7B-3B43-E311-8F37-848F69FD294C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/4A939678-3B43-E311-AF5A-7845C4FC36AD.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/E8BF8591-3B43-E311-9AE5-848F69FD2904.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/CE055AA5-3B43-E311-94A7-008CFA001444.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B04F85BC-3B43-E311-A160-00266CF27130.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/6E291BC1-3B43-E311-83A4-848F69FD29D3.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/CEB430DE-3B43-E311-B943-00266CF258D8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B47C39B9-3B43-E311-87E5-7845C4F91450.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/12062649-3C43-E311-BFFC-00266CF9BDFC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/767DE478-3B43-E311-885B-7845C4FC373A.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A08F3B77-3B43-E311-8E5A-848F69FD4EDD.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/4EB96E76-3B43-E311-8499-848F69FD283E.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/5E49D1B8-4B43-E311-B429-7845C4FC3B6F.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/ACC8FCC3-4C43-E311-8DAC-00266CF9B274.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/365C9EF5-F742-E311-9459-7845C4F92EC7.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/3EFD1967-F842-E311-9E4D-00A0D1EE95AC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/7E6AFED2-4043-E311-9790-00266CF9ABA8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D85AD4C9-4243-E311-957D-7845C4FC39FB.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/3AF835CD-4243-E311-9FA4-7845C4FC3641.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/4A7276F5-FF42-E311-BC8A-00A0D1EEE3B0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/F2DC6869-FF42-E311-B49D-7845C4FC3893.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/CCD4397B-FF42-E311-A949-00A0D1EE8D7C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/F0631108-0643-E311-A361-848F69FD2955.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/26370BEC-0543-E311-A17D-00A0D1EE9238.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/86E34DCB-DC42-E311-862F-00266CF25DC4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/2C0DB27C-DE42-E311-9E02-00A0D1EE8F34.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/225D3EC7-DF42-E311-B9C2-00A0D1EE2990.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/CC48CD7E-E042-E311-BF01-7845C4FC3614.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/C41D2360-E142-E311-B9D3-848F69FD2931.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/AA8FFA70-2843-E311-945F-00266CFAEA68.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/86D7936D-2843-E311-B98B-00266CF9B86C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/58555D3C-FF42-E311-93EB-00A0D1EE29B8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/7A712D58-2043-E311-BA29-7845C4F91621.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/3876AC87-FA42-E311-A656-180373FF93CE.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/CA19A4C3-FA42-E311-A854-00A0D1EE8E78.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/883FFA06-FB42-E311-A2E8-7845C4FC346A.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/EA4E6A64-2043-E311-B30A-008CFA002ED8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/700257AC-2043-E311-8A0F-00266CF25878.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D835FBBC-2043-E311-8B98-008CFA002FA4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/94C31502-0943-E311-89A5-7845C4F92F87.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/E4BC625B-0A43-E311-8AF8-00A0D1EE2990.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/9434D4F3-0A43-E311-8A02-7845C4FC3B48.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/721341BF-2443-E311-8E68-7845C4FC39AD.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/882C252E-2543-E311-B0B9-00266CF9B630.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/26A6B017-2543-E311-8D96-7845C4FC3C02.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/927B47C7-2443-E311-AD71-00266CFAE074.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B42E6839-2643-E311-B296-00266CF9B1C4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/FCA06A18-4443-E311-8807-848F69FD2817.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/AE2D3717-4443-E311-8645-848F69FD2970.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/9E5F231D-4443-E311-BFC5-00266CFADE34.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/16BC6B41-4443-E311-BE78-00266CFAE268.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B65B5842-4443-E311-B750-848F69FD28AA.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/627A54C6-4443-E311-AC8C-00266CF9ABA8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/4676EFFF-4443-E311-A490-00266CF9B414.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/92DBB54B-4543-E311-AE1B-00266CF271E8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/3240FBF1-4543-E311-A7DD-848F69FD45A4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/28AADC8E-2143-E311-B474-00A0D1EEF204.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D02F9BA1-2243-E311-A013-00A0D1EEE68C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/BC736E3F-2243-E311-B843-00266CF9AFF0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/9084E150-2243-E311-AF63-7845C4F91621.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/BA436C5F-2243-E311-8F59-008CFA002FA4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A436C82A-2D43-E311-9560-7845C4FC3A19.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/FA68052E-2D43-E311-9B51-7845C4FC3A7F.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/1284D07A-2D43-E311-93DD-848F69FD2823.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/68EE9993-C242-E311-BD3C-848F69FD29CA.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/5ACC0DC2-C342-E311-A303-848F69FD29CA.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/F289A0A9-D542-E311-9E6B-008CFA000744.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/26011160-D642-E311-9E99-848F69FD290A.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/2475F963-D742-E311-AFF3-848F69FD28AD.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/C4DDD7A8-D842-E311-BF4A-180373FF8DE0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/68434E4C-D842-E311-9117-7845C4FC365C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/F27A80A7-D942-E311-880F-00266CFAE318.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/F4348C74-E542-E311-A395-00266CF9B318.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/4C363610-E742-E311-9F09-001D09FDD7D4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/7E665327-E742-E311-8AB6-00266CF9B318.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D095C469-E942-E311-8467-00266CF253C4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/00C5D534-EA42-E311-8B42-008CFA001EE4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/9C53F5B1-EB42-E311-81D9-848F69FD4DCC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/4C7D9BF4-EC42-E311-B033-7845C4F932D8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/BC491589-ED42-E311-AB93-00266CFAEBA0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/447A5670-3943-E311-A45B-00266CFAE890.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B45F027E-3943-E311-B3A2-00266CFADEC0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B438098A-3943-E311-82CB-00266CFADFCC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/FE85CE8C-3943-E311-8F68-008CFA001FDC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/AAC716DE-3943-E311-9389-008CFA001444.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/105744FC-3943-E311-B8A4-7845C4F91450.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/78C45D17-3A43-E311-8797-00266CF27130.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B27A2066-2943-E311-B4FE-7845C4FC3998.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/BC014B73-2943-E311-B288-00266CF9B1C4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/50936944-2843-E311-ABEB-00266CF9ACB0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A6A9B20B-2943-E311-B938-848F69FD29D0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/E85E59AE-2943-E311-8A56-848F69FD2988.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/6448F4AD-4743-E311-A923-00266CF9BCC4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B6D16D3F-4843-E311-B95E-7845C4FC37AF.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/C6B6732A-4943-E311-8100-008CFA002E80.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/74F484A6-4843-E311-BD1D-848F69FD46C1.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/C0806F29-4943-E311-B8A0-008CFA000898.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A216FE27-4943-E311-8980-008CFA0027B4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/22BAF73C-2143-E311-975D-00266CFAE788.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/5292CB2C-3C43-E311-8EE0-848F69FD29DF.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B8395ECD-3C43-E311-B9E9-00266CFADFCC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D257EE32-3C43-E311-80CA-7845C4FC359F.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/C0AF68B8-3C43-E311-9505-00266CFADEC0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/9E654818-FB42-E311-91D6-00266CFAEBF8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/4220EC88-FB42-E311-98BC-7845C4FC3A61.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/24BCAE06-FC42-E311-B0E1-7845C4FC370A.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/F08C8E75-1043-E311-875E-00266CF25218.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/C8878466-1043-E311-998A-180373FF8446.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/E6844052-1143-E311-9D9B-848F69FD2A30.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/641D50D4-2C43-E311-AFD6-00266CFAEBF8.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/46AF6934-B342-E311-BC91-001D09FDD831.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A6797ADB-1143-E311-A337-848F69FD45A4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/FA3784AB-1143-E311-8BC5-00A0D1EE95CC.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A0A197FE-FF42-E311-816D-7845C4F9329C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/9E3AFB90-0E43-E311-97D1-00266CF97FF4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B64BA480-1543-E311-976B-848F69FD4ED4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B0D70486-1543-E311-8A82-00266CF9B874.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/5413D4D6-1543-E311-97D9-00266CFAE464.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/62634F60-0343-E311-8B6F-00266CFAE838.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/80FF35C4-0343-E311-8C83-7845C4FC3CA1.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/DE98DD11-0343-E311-995C-00A0D1EE2A50.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/16B752F7-0643-E311-B995-00266CF25C88.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/48C4F074-2E43-E311-8C10-848F69FD2955.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/A8C8A475-2E43-E311-B23F-008CFA008C94.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/92254507-2E43-E311-B380-848F69FD47A5.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/92308093-9C43-E311-8C4A-848F69FD2910.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/C415800F-2343-E311-BD00-00A0D1EE9274.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/9421022F-C942-E311-A637-848F69FD28C2.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B633E345-CA42-E311-9089-7845C4FC370D.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/C004E1C2-1E43-E311-BBD2-848F69FD28FE.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/548019A7-1F43-E311-8CC7-001D09FDD7D4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/C2F808BD-1F43-E311-83BB-00A0D1EEF204.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/8E9105C1-9B43-E311-8B39-848F69FD28FE.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/AC405E17-9D43-E311-89A7-7845C4FC3C65.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/F473C043-9D43-E311-9826-00A0D1EE8A0C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/0E2B44C4-6843-E311-BB9C-00266CF9C210.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/70B4E205-7543-E311-BA85-001D09FDD7DA.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/7A428241-7543-E311-83C5-008CFA0104E0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/0AC4C04F-A843-E311-97B3-00266CF275E0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/803B7F5D-A843-E311-828E-00A0D1EE8A10.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/D8D93470-A843-E311-B003-848F69FD29BB.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/B6FB9573-9E43-E311-B3C6-00A0D1EE8A0C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/7E252DCB-A943-E311-8523-00266CFAE69C.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/90412937-AA43-E311-81D5-00A0D1EE8A10.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/8E5BB61F-C743-E311-BD2B-7845C4FC3653.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/64678910-C743-E311-8D99-7845C4FC3668.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/2A147F57-C743-E311-9A4E-008CFA0008C4.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/E23CC40E-C743-E311-9CC0-00266CF9ACB0.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/768B7264-C743-E311-8559-001D09FDD80D.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/C05A7268-7843-E311-BECF-00266CF25590.root',
          '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/20000/F0A31DB7-7843-E311-B65A-7845C4FC3C98.root',
  ] )
elif (data and options.lepton == "muon"):
  readFiles.extend( [
         '/store/data/Run2012A/SingleMu/AOD/22Jan2013-v1/30000/FEDCB8E2-5270-E211-8FD6-00266CFFBC38.root'
  ] )
elif (data and options.lepton == "electron"):
  readFiles.extend( [
         '/store/data/Run2012A/SingleElectron/AOD/22Jan2013-v1/30002/72352080-9372-E211-B431-00266CFFA1AC.root'
  ] )
  process.source.skipEvents=cms.untracked.uint32(1800)
else:
  readFiles.extend( [
          '/store/user/mgosseli/mc/TT_TuneZ2_7TeV_madgraph_FASTSIM_172_5GeV_matchingdown_v1/TT_TuneZ2_7TeV_madgraph_FASTSIM_172_5GeV_matchingdown_v1_139.root',
  ] )

secFiles.extend( [
               ] )

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

## configure geometry & conditions
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if os.getenv('CMSSW_VERSION').startswith('CMSSW_4_1_'):
  process.GlobalTag.globaltag = cms.string('START41_V0::All')
elif os.getenv('CMSSW_VERSION').startswith('CMSSW_4_1_'):
  process.GlobalTag.globaltag = cms.string('START42_V17::All')
elif os.getenv('CMSSW_VERSION').startswith('CMSSW_5_3_'):
  process.GlobalTag.globaltag = cms.string('START53_V23::All')
if data:
  process.GlobalTag.globaltag = cms.string('FT_53_V21_AN5::All')

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")

if not data:
    ## MC weights
    process.load("TopAnalysis.TopUtils.EventWeightMC_cfi")
    
    ## B-JES
    process.load("TopAnalysis.TopUtils.EventWeightBJES_cfi")
    


# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeBJES.root')
)

# register TreeRegistryService
process.load("TopMass.TopEventTree.TreeRegistryService_cfi")
process.TreeRegistryService.treeName  = "eventTree"
process.TreeRegistryService.treeTitle = ""

#process.MessageLogger = cms.Service("MessageLogger")
process.content = cms.EDAnalyzer("EventContentAnalyzer")

## end path   
process.path = cms.Path(
                        process.genJetParticles *
                        process.ak5GenJets *
                        process.EventWeightBJES
                        )


