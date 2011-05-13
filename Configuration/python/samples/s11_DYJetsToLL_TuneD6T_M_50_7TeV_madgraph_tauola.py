import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_9_1_XKc.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_99_1_JED.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_98_1_C45.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_97_1_4Hr.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_96_1_L6F.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_95_1_l9n.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_94_1_fJy.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_93_1_cKn.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_92_1_DeU.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_91_1_YN4.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_90_1_zBN.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_8_1_gGm.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_89_1_1Bq.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_88_1_Mkt.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_87_1_dfJ.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_86_1_XKO.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_85_1_MDX.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_84_1_oYv.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_83_1_Mvr.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_82_1_Nch.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_81_1_Qxc.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_80_1_6jc.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_7_1_7dP.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_79_1_U0w.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_78_1_voj.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_77_1_0Jn.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_76_1_A2n.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_75_1_807.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_74_1_ZiJ.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_73_1_gjB.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_72_1_baA.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_71_1_Cik.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_70_1_Szd.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_6_1_70u.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_69_1_81F.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_68_1_37q.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_67_1_nPU.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_66_1_Sc7.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_65_1_2uU.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_64_1_2AV.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_63_1_eV5.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_62_1_SNc.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_61_1_Cs1.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_60_1_BRF.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_5_1_Pcr.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_59_1_qmZ.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_58_1_NkK.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_57_1_Hh4.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_56_1_4d6.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_55_1_OLn.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_54_1_D5i.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_53_1_024.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_52_1_ysW.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_51_1_8Av.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_50_1_80N.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_4_1_wW8.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_49_1_v3H.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_48_1_wfq.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_47_1_yAj.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_46_1_dsQ.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_45_1_lRk.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_44_1_q8N.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_43_1_X9s.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_42_1_vTL.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_41_1_kCO.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_40_1_0T4.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_3_1_MDr.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_39_1_5jI.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_38_1_Kdq.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_37_1_m4H.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_36_1_p3R.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_35_1_fSf.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_34_1_QdW.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_33_1_u7s.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_32_1_NoB.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_31_1_qyy.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_30_1_L2l.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_2_1_bSf.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_29_1_Vwc.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_28_1_ali.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_27_1_rw0.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_26_1_1jM.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_25_1_cs8.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_256_1_jEd.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_255_1_1ql.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_254_1_Q2f.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_253_1_sYM.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_252_1_gjF.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_251_1_FvL.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_250_1_1eo.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_24_1_oIa.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_249_1_Yk3.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_248_1_UAC.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_247_1_CKW.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_246_1_97G.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_245_1_1X5.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_244_1_tzP.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_243_1_AYS.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_242_1_Umu.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_241_1_iAx.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_240_1_jhq.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_23_1_aYh.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_239_1_DMD.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_238_1_qZe.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_237_1_bUR.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_236_1_rDE.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_235_1_3YC.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_234_1_2ry.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_233_1_IVu.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_232_1_sqP.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_231_1_5lL.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_230_1_4sA.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_22_1_e2B.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_229_1_SIa.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_228_1_L5g.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_227_1_QLz.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_226_1_s17.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_225_1_Ky9.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_224_1_PQZ.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_223_1_XaO.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_222_1_Frf.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_221_1_ktO.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_220_1_2Rp.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_21_1_Ah7.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_219_1_3Il.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_218_1_ICw.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_217_1_mGM.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_216_1_VtS.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_215_1_wp5.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_214_1_WIV.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_213_1_cY2.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_212_1_E7H.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_211_1_jdt.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_210_1_XN2.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_20_1_3zb.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_209_1_c6P.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_208_1_IJv.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_207_1_Ofv.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_206_1_XnV.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_205_1_TcM.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_204_1_tfb.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_203_1_43f.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_202_1_vOJ.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_201_1_9mG.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_200_1_ggs.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_1_1_ntG.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_19_1_DmY.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_199_1_3Ia.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_198_1_ijQ.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_197_1_h14.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_196_1_zKc.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_195_1_UvZ.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_194_1_o2O.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_193_1_52t.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_192_1_Dgd.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_191_1_ydg.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_190_1_7zZ.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_18_1_Ptb.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_189_1_o31.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_188_1_0Rd.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_187_1_ppw.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_186_1_13x.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_185_1_VRW.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_184_1_sGp.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_183_1_Ewl.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_182_1_HGr.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_181_1_HKG.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_180_1_Udi.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_17_1_skX.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_179_1_4SY.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_178_1_0Px.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_177_1_NUs.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_176_1_5Vg.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_175_1_Y6T.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_174_1_jut.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_173_1_jrZ.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_172_1_DOr.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_171_1_MwS.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_170_1_auj.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_16_1_pmC.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_169_1_Opc.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_168_1_AXP.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_167_1_U31.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_166_1_xiy.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_165_1_dlf.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_164_1_ORv.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_163_1_rCF.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_162_1_FTE.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_161_1_GMB.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_160_1_IDT.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_15_1_eNA.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_159_1_RUH.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_158_1_1lu.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_157_1_3cg.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_156_1_Nzg.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_155_1_dxe.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_154_1_IQ5.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_153_1_Rf7.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_152_1_CPB.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_151_1_nhw.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_150_1_xRr.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_14_1_T4Y.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_149_1_DlW.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_148_1_jbc.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_147_1_ViQ.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_146_1_IiM.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_145_1_DaR.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_144_1_qA3.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_143_1_ZH4.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_142_1_SWs.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_141_1_NRS.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_140_1_Ml6.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_13_1_R7h.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_139_1_YIq.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_138_1_VDZ.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_137_1_D7n.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_136_1_Stw.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_135_1_ycx.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_134_1_C0h.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_133_1_uxK.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_132_1_UZJ.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_131_1_Xmi.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_130_1_7zR.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_12_1_6R3.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_129_1_Fku.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_128_1_EJB.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_127_1_Ohv.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_126_1_JzH.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_125_1_aFm.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_124_1_diT.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_123_1_6gH.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_122_1_DKu.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_121_1_Xku.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_120_1_aFN.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_11_1_xlV.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_119_1_C7R.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_118_1_r7n.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_117_1_HNQ.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_116_1_3Ht.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_115_1_22M.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_114_1_qpL.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_113_1_TUD.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_112_1_1Bp.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_111_1_xiF.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_110_1_4tA.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_10_1_fhH.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_109_1_Cgz.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_108_1_jaj.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_107_1_rs9.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_106_1_ktc.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_105_1_uct.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_104_1_aWg.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_103_1_AsF.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_102_1_d0G.root',
       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_101_1_JlS.root' ] );
readFiles.extend( [

       '/store/user/wbehrenh/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola/Spring11-PAT/6e6559812e09b52af172f27db20ae337/mc2pat_100_1_F0e.root' ] );


secFiles.extend( [
               ] )

