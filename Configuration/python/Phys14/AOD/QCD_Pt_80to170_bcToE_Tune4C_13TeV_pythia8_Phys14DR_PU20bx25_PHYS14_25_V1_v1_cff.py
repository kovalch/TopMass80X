import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00B25D46-816F-E411-BF35-0025905B85E8.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00BB4A60-DB6E-E411-90FF-002618943904.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00F801BF-5970-E411-9C3A-002590593878.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0201236B-E26E-E411-B30E-002618943900.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/02356941-E76E-E411-A736-0025905A612C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/049A5D3E-CB6E-E411-A15F-00261894397F.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/04C6A341-DD6E-E411-A75F-0025905A6118.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/06834BB7-C96E-E411-AE96-00248C0BE018.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0AB82CF3-816F-E411-9D10-003048FFD75C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0C11B73C-E06E-E411-8086-003048FFD7C2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0CA68853-CB6F-E411-8E13-003048FFCC2C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0E887E06-EA6E-E411-86C5-003048FFD740.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/10D83B90-8C6F-E411-86AE-0025905A6126.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/12980AAD-E46E-E411-83F6-0025905B8562.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/12D5C27D-736F-E411-BE80-002618943982.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/144174DE-7A6F-E411-AFF3-003048FFD75C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/16342435-7A6F-E411-99E2-0025905A607E.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/18066A9B-7C6F-E411-8770-002618FDA204.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/18CB4961-946F-E411-B43D-0025905A48E4.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/18EE07DC-CF6E-E411-89C3-0025905B860C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1A640B45-D26E-E411-8DC0-002618943882.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1C8B1DA2-DF6E-E411-878E-0025905938B4.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1CC0E61D-DE6E-E411-B4D0-003048FFCC1E.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1CF603F5-C26E-E411-8D90-0025905B8576.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/202F5EBF-DE6E-E411-96EE-0025905A60F8.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/22F67361-DB6E-E411-A7CC-0025905B858A.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/24036B4B-8B6F-E411-A711-0025905A48F0.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/24A9D7FF-E26E-E411-A7C2-0026189437F8.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2622E2B0-E36E-E411-AEDC-0025905A6092.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/28242FA6-836F-E411-BC0F-0025905A60A8.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2C4DECCF-6F6F-E411-A98C-0025905A6076.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2E0103FF-876F-E411-9CE8-0025905964A6.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/306E620B-C26E-E411-AA0C-0025905B8562.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3289077D-886F-E411-9B60-0025905A60F4.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/34CC7895-DD6E-E411-8DE5-0025905A6056.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3A894478-D06E-E411-9EA5-002618943800.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3C1FCFF2-BA6E-E411-971F-0026189438A5.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3C31E4AC-E36E-E411-A232-002618943838.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3C59D46C-E86E-E411-B9FF-0025905B85F6.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3E240A71-CD6E-E411-BD40-002618FDA248.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3E2E94C8-6D6F-E411-93BF-0025905A60B4.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/405D23E6-716F-E411-BA34-0025905822B6.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/425A2D5B-D16E-E411-A5A7-0025905A6122.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/440F7884-8C6F-E411-A5A4-0026189438A2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/441386CD-6D6F-E411-B1B1-0025905A60B8.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4469EC52-FA6E-E411-AF0E-0025905B8596.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/44CB7BB9-DE6E-E411-AAA2-0025905A60B6.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/46228B7B-866F-E411-B3A3-0026189437EC.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/468F1ADA-806F-E411-9A09-0025905A4964.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/46A5E107-BD6E-E411-82AB-003048FFD76E.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4A2AC777-926F-E411-ACC7-00261894393C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4C33BC89-856F-E411-AACC-003048FFD7BE.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4E22372D-8A6F-E411-8CD6-002618943957.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5047EACE-E16E-E411-B6F1-0026189438A0.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/50A51739-E76E-E411-A8C3-00261894397B.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/50B92A85-826F-E411-A02F-0025905B855C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5200ABBF-6F6F-E411-909C-0025905B857C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/523A2E46-D26E-E411-8310-002618943857.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/52B5397C-926F-E411-8B73-0025905938B4.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5A22FE3B-716F-E411-AC53-003048FFD75C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5C200D07-E36E-E411-A306-0025905A48EC.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5C681775-776F-E411-823B-0025905A6082.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5EA847FF-CB6E-E411-889B-0025905A60E4.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5EF45D9F-CE6E-E411-803F-0026189438F2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/60AD6203-9B6F-E411-9CF9-00248C55CC40.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/60E061DF-7D6F-E411-9343-0025905B860C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/60F17BB1-E56E-E411-8BD4-0025905A60D6.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/62BF2727-DB6E-E411-9DA8-0025905A605E.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/62F5BD3D-7A6F-E411-B313-0025905A48BA.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/661C32AC-E56E-E411-88B2-003048FFCB74.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/661D73E8-616F-E411-9064-0025905B85E8.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/661DFFD3-B86E-E411-AA63-0026189438B5.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6669C4D6-CC6E-E411-B1BD-00259059642E.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/66709C37-D86E-E411-92AC-0025905938A8.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/66DAB6AB-E46E-E411-98EA-0025905B85A2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/68109DA3-CE6E-E411-9C42-0025905AA9F0.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6876E9D9-E06E-E411-9B42-003048FFCB9E.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/68E53130-736F-E411-8D89-0025905A60A6.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6A62D180-E16E-E411-A190-0026189438D8.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6C9E9AC7-D96E-E411-9599-0025905A48D6.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6E74E2C0-746F-E411-9683-0025905A48D8.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/70528D5D-D56E-E411-8E90-0025905A60BC.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/70E2F579-976F-E411-9FE7-002618943972.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/70E74FF1-CB6E-E411-841E-002618FDA248.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/724C029C-896F-E411-9336-0025905B8606.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/768C4839-E06E-E411-BAD5-0025905964BA.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7CF8FB72-5D6F-E411-AA83-00261894390E.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/80932172-6C6F-E411-8B0D-0025905A6122.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/80A53F83-856F-E411-AA38-002618943832.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/80F99FAC-E46E-E411-A0AF-0025905A6138.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/822CE4DD-846F-E411-BA42-00261894393E.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/826B8FD3-E16E-E411-854C-0025905B85D0.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/82834668-E26E-E411-B109-002618943898.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/843D8CD4-D36E-E411-A8D9-0026189438A2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/863270ED-9D6F-E411-99E0-002618943874.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/88A3CCAD-906F-E411-BACC-0025905A48F2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/88D6DD27-876F-E411-9B74-0025905A6082.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/88EE3725-C86E-E411-A95A-003048FFCBFC.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8AAB9068-D16E-E411-B7DC-002618943800.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8ACCB9B6-8D6F-E411-9F6D-0026189438E7.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8C4B26B0-C56E-E411-9C10-0026189437FE.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8ED8C97D-976F-E411-8F8D-0025905A609E.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/90BDC791-9E6F-E411-9692-002618943949.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/925DFD60-C66E-E411-B090-0025905B8562.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9414948C-EB6E-E411-B878-0025905A6082.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/94588F88-8C6F-E411-96C0-0026189438E4.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/96912DBE-5970-E411-867A-0025905A48B2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/96A0A2EC-7F6F-E411-90B2-0025905A60B2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9AC7F538-CF6E-E411-9378-00261894394F.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9CB83883-766F-E411-9710-002618FDA248.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9E4AA995-DD6E-E411-BD88-0025905938D4.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9E5A173E-CB6E-E411-8BA3-0025905B860C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9E671392-DC6E-E411-8664-0025905A60B2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9E9E33A3-CE6E-E411-92CD-0025905B857C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A2078E4D-756F-E411-ABB4-003048FFD7C2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A27A0343-CF6E-E411-B619-0025905A48D6.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A45B27C8-D96E-E411-8BFB-0025905A48D0.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A83D0603-896F-E411-82D7-0025905A6084.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AC7DEA05-906F-E411-BF4D-0025905A48BA.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/ACE7314A-8B6F-E411-AF1F-0025905A6090.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B64EF7BA-D06E-E411-8E14-0025905A60C6.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B6B5CAE3-946F-E411-8BFD-002618943937.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B6D77BA1-7C6F-E411-91B9-003048FFCC2C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B833DCB8-8D6F-E411-844D-002618943981.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B861313C-916F-E411-8F72-00248C55CC97.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BA04355F-D56E-E411-B6C2-0025905A60DA.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BE919C6E-CD6E-E411-BC68-0025905A60B4.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BEC535E6-BF6E-E411-B48B-003048FFCBFC.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BEDB68D4-D36E-E411-AFB8-0025905A607E.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C2B1EF93-E96F-E411-97C8-0025905B8598.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C2B81D01-E46E-E411-9B4B-0025905A60A0.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C47C125C-786F-E411-AEBF-0025905A60D2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C4D13E22-C86E-E411-B828-002618FDA248.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C4FEBEB2-8E6F-E411-911F-0025905AA9CC.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C671D2ED-756F-E411-9906-0025905B8582.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C698ADF4-CC6E-E411-916E-0025905A60DA.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C6AD8CF5-E66E-E411-A442-0025905A60DE.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C6BA3968-786F-E411-B1C3-0025905A606A.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C85F9152-E56E-E411-8BE1-002618943953.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CA7DF93B-916F-E411-B88A-0025905B85E8.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CCF025C5-8E6F-E411-A8CA-0025905B8572.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CE0CF803-7E6F-E411-8DC0-0025905A6066.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CE4F2AF2-CA6E-E411-AEBC-0026189438AD.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CE4F3E0D-786F-E411-8941-0025905B8596.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D0819B0E-906F-E411-A012-003048FF9AA6.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D208D934-8A6F-E411-B8F7-0025905B860E.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D24319A7-836F-E411-99F9-0025905A60AA.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D2EB5E6E-E86E-E411-B533-0025905A606A.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D48E0FCF-E16E-E411-9226-002618FDA204.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D49D6640-DD6E-E411-8551-002618B27F8A.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D4ACB625-7F6F-E411-BC49-0025905A60FE.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D852DA4E-886F-E411-9178-0026189437F5.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D879EA7D-866F-E411-ADB5-003048FFCC2C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DAAA7F18-746F-E411-BFBD-0026189438E7.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DAED01CA-D06E-E411-9190-0025905B8598.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DC0BD609-886F-E411-B2DB-0025905B85AA.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DCEBCA03-9B6F-E411-A8E3-003048FFD736.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E05B10B7-746F-E411-A50A-0025905B8606.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E0E12635-D86E-E411-85D0-0025905B85EE.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E2D91A19-836F-E411-A463-0025905A60B8.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E45D2A9B-896F-E411-BC92-0025905A611C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E81582DF-8B6F-E411-8966-0025905A60F8.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E830B9B9-C96E-E411-B107-002618943857.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E873ED7C-D06E-E411-9D12-0025905A6064.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E88A633A-716F-E411-A945-002618FDA28E.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/EEDAC9B5-936F-E411-9610-00261894393F.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F014241C-DE6E-E411-A468-0025905B860C.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F021A75D-DF6E-E411-9547-0025905A605E.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F0D758A0-DF6E-E411-AEE3-0025905A60D2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F247F1DA-E06E-E411-9C89-0025905A60BE.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FA425D2E-736F-E411-9F66-0025905938D4.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FA899F55-FA6E-E411-B2C5-0025905B85B2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FE9DCFDD-CF6E-E411-B825-0025905A48F2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/085EDE9B-886F-E411-84D7-0025905A60DA.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/0A28D0E5-946F-E411-A47B-0025905938A8.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/30687DD4-8D6F-E411-B9E6-0025905A60B0.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/48864F4D-806F-E411-9FA3-0025905A48B2.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/64CF82B4-766F-E411-AE83-0025905A608A.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/72F476EA-786F-E411-BAB6-0025905A6132.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/8A6DA12B-8B6F-E411-AF3F-0025905A60D0.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/A6228612-DF6F-E411-94CA-0025905A60D0.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/ACA99A11-DF6F-E411-AC9F-002618943935.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/B40FFE2B-716F-E411-8F13-0025905A608A.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/EAD4A470-8A6F-E411-B8B8-0025905A609A.root',
       '/store/mc/Phys14DR/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/FA9CE298-886F-E411-BB02-0025905B8572.root' ] );


secFiles.extend( [
               ] )

