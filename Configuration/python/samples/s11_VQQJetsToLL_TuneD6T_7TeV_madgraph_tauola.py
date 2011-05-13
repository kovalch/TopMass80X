import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_9_1_DJs.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_8_1_Nf4.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_7_1_TnA.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_75_1_sad.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_74_1_jQM.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_73_1_rxR.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_72_1_Nqg.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_71_1_DUP.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_70_1_d4J.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_6_1_0Z5.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_69_1_RmS.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_68_1_HZp.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_67_1_exv.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_66_1_0iV.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_65_1_WVj.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_64_1_BIr.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_63_1_dDU.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_62_1_Lt9.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_61_1_dll.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_60_1_e0B.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_5_1_qGS.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_59_1_6wZ.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_58_1_Fci.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_57_1_COr.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_56_1_4Eu.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_55_1_pCq.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_54_1_5ZS.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_53_1_B4l.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_52_1_LQ6.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_51_1_H4d.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_50_1_27o.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_4_1_Bk6.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_49_1_wd6.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_48_1_CdA.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_47_1_ONu.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_46_1_bJz.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_45_1_992.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_44_1_cnf.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_43_1_PKj.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_42_1_6r5.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_41_1_j3s.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_40_1_jVb.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_3_1_bpi.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_39_1_1kF.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_38_1_owM.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_37_1_3qP.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_36_1_OAO.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_35_1_yG5.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_34_1_q3z.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_33_1_iC6.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_32_1_5x7.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_31_1_g7M.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_30_1_MqE.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_2_1_9rz.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_29_1_R8T.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_28_1_b8B.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_27_1_ZFu.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_26_1_Cxa.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_25_1_Kvd.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_24_1_2LY.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_23_1_EHo.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_22_1_tod.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_21_1_5De.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_20_1_0Ss.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_1_1_usn.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_19_1_DRz.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_18_1_X2T.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_17_1_Uws.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_16_1_0oE.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_15_1_C8U.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_14_1_BwT.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_13_1_5Sx.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_12_1_PuQ.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_11_1_WUR.root',
       '/store/user/wbehrenh/VQQJetsToLL_TuneD6T_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_10_1_97r.root' ] );


secFiles.extend( [
               ] )

