import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/FE81BD5F-30EF-E111-8060-1CC1DE046FB0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/FE41F845-68EF-E111-808B-00266CFFC7E0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/FCCC8B0C-03EF-E111-AA99-AC162DACC3F0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/F8E6BF18-11EF-E111-8296-AC162DA8C2B0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/F82B15C7-1DEF-E111-A229-D8D3855B79E8.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/F694C03B-46EF-E111-BE9F-00266CFFCD00.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/F0EDFE5F-27EF-E111-8A69-00266CFFB7D0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/F0343F10-EBEE-E111-AC11-00266CFFB7D0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/EE6E628F-6AEF-E111-AF14-00266CFFCC54.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/EC44CF4C-10F1-E111-85DA-00237DA1AC24.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/EC116F60-79EE-E111-81B5-00266CFFBCB0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/EA060AE7-11F1-E111-ADA4-00266CFFCD14.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/E8026A57-6FF1-E111-9222-0017A4770C1C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/E67928F1-01F1-E111-8726-0017A477102C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/E6727F3B-F6EE-E111-BD51-1CC1DE047F98.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/E43C6307-27F1-E111-B138-0017A477082C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/E21738CF-44EF-E111-AE8B-00266CFFCD00.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/DED52B7E-24EF-E111-8D77-1CC1DE050110.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/DEA68D49-D6EE-E111-A67C-00266CFFBF50.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/DCEB764D-F6EE-E111-9DEB-1CC1DE1CEAEE.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/DC87BB4F-EBEE-E111-AC91-0025B3E0216C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/DC246062-F6EE-E111-A400-00237DA1DDE4.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/DAE435E6-1FEF-E111-B00B-00266CFFBF68.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/DA52C33C-5DF1-E111-A894-AC162DACB208.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/DA113228-4EF1-E111-B665-0025B3E0228C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/D882CEA3-08F1-E111-B10D-1CC1DE0437C8.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/D4B03C19-28F1-E111-BC97-1CC1DE051118.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/D21EB4C0-C7EE-E111-AE38-00266CFF0ACC.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/D02DCB8E-F8EE-E111-B7C9-0017A4770010.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/CE5F1346-32EF-E111-A2AE-00237DA12CBE.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/CCA0412B-77EE-E111-9776-1CC1DE046FB0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/C84D836C-FAEE-E111-966D-0017A4770828.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/C837738B-09F1-E111-8B9B-1CC1DE051118.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/C4545B28-EBEE-E111-AF79-00266CFF0234.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/C2D21BF3-2EEF-E111-818B-00266CFFCC50.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/BC7F089D-3CEF-E111-ADA6-0017A4771020.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/B8F946C1-20F1-E111-A729-00266CFFBF38.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/B83F3945-14F1-E111-BC7B-00266CFFCD6C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/B6D66D36-42EF-E111-A874-0017A4770838.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/B4A06285-09EF-E111-A1BB-0017A4770438.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/B2C7985B-0FEF-E111-BD6A-00266CFFBF64.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/B0446854-1FF1-E111-A2F3-00237DA13C2E.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/ACA3FB79-F9EE-E111-BACB-0017A4771020.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/AAFA6C61-FBEE-E111-A47F-AC162DAC3428.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/AAF2DE65-46EF-E111-87E1-00266CFFC848.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/AA4B51C1-B6EE-E111-9CEB-00266CFFCAC0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/A4C397B6-12EF-E111-8918-1CC1DE1CDCAE.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/A4296676-03F1-E111-8DCE-1CC1DE1D1F80.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/A29805A1-B0EE-E111-A0EE-00266CFFC13C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/A26A4CDA-17EF-E111-AACF-0017A4770034.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/A0AD8B8A-76F1-E111-B0A4-00237DA1AC2A.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/9A66A132-35EF-E111-8319-1CC1DE046F00.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/98D6DC18-F6EE-E111-AD83-00266CFFC9EC.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/98AE18C7-22F1-E111-B3CD-0017A4770024.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/92D3525C-31EF-E111-B51B-00266CFFCD6C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/92B44793-36F1-E111-8200-00266CFFCA1C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/90A6EC18-30F1-E111-8506-1CC1DE046F18.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/9095BE3D-04F1-E111-9408-00237DA1B34E.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/8ED8B748-75EF-E111-985E-00266CFFBE14.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/8EA78267-D1EE-E111-B39E-00266CFFC198.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/8CABC4F0-D6EE-E111-AC25-00266CFFBC3C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/8AB4413E-F6EE-E111-92FB-AC162DACE1B8.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/86D9E024-9BEE-E111-8FAC-00266CFFC43C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/862E982A-79EF-E111-B445-0017A477041C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/84E25702-DBEE-E111-A890-0017A477043C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/82C7D7FF-06EF-E111-ACB2-1CC1DE0530F8.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/80DCC6FF-18EF-E111-B87C-00266CFFBDB4.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/80A8EA1A-6DEF-E111-8A31-AC162DACC3F0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/80674BBE-2CEF-E111-A9C4-1CC1DE040FA0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/801D3610-0DEF-E111-9E73-1CC1DE041FD8.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/7E3802DD-ABEE-E111-8A7C-1CC1DE04FF50.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/7C72321B-2DEF-E111-92DA-00237DA0F456.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/7ADF7B7F-0DF1-E111-B6F7-78E7D1E4B6E8.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/7A559345-BCEE-E111-A5D5-AC162DA87230.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/7603C836-1CF1-E111-870D-0017A4770034.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/72C278FB-FCEE-E111-9FDD-1CC1DE0500F0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/722B0067-F6EE-E111-A5A4-00266CFF0ACC.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/70BA9FE1-60F1-E111-9A48-00266CFFC948.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/70727D9B-A5EE-E111-9F8E-AC162DACB208.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/7069F722-43EF-E111-BD82-1CC1DE1CED22.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/6C8E36D2-D4EE-E111-AE2C-00266CFFCB80.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/68C7D20E-81EE-E111-947E-00266CFFBF34.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/6860ED86-02F1-E111-A0F3-1CC1DE1D1F80.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/68226CAA-0CF1-E111-8EA5-0017A4770410.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/649380B2-05F1-E111-AC60-00266CFFBDE8.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/644CA1D2-06F1-E111-AC32-00266CFFB7D0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/62C1DC2E-42F1-E111-8FB3-0017A4770C18.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/62BDF18A-0EEF-E111-B424-0025B3E020D0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/620141D8-32EF-E111-A357-AC162DAB0B08.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/60DA8914-4AEF-E111-B62D-0017A4771028.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/60555211-55EF-E111-95B7-D485645C8BC8.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/5EC5B2C7-2FEF-E111-81F3-AC162DACC3F8.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/5E916C9F-5EEF-E111-A943-1CC1DE1D03EA.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/5CF789D2-25F1-E111-9B2C-AC162DA87230.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/5C6B0CBF-2BEF-E111-BD50-1CC1DE1D03DE.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/5AC5A2BC-65F1-E111-9DE2-0025B3E020D0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/58072B56-0EF1-E111-9075-0017A4770410.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/54E1EF2A-4AF0-E111-B118-1CC1DE052048.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/54C7851D-41EF-E111-9548-0017A4770838.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/549D883F-F6EE-E111-8D81-AC162DA87230.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/542BEB8C-4FEE-E111-8797-1CC1DE04DF20.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/50F985CA-01EF-E111-925C-1CC1DE1CED22.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/4CD37D10-32F1-E111-BEBD-78E7D1E4B874.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/4A213E63-3EEF-E111-93E3-00266CFFBC38.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/4A12D97F-FFEE-E111-BB54-00266CFF0034.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/4A0D5FD4-4CEF-E111-BFDA-0017A4770C28.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/46AE4517-1BF1-E111-A359-00266CFF0608.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/466E1216-74EE-E111-9697-00266CFFCCBC.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/446B22C2-0AF1-E111-A084-1CC1DE0503C0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/44111D9A-16EF-E111-97EF-1CC1DE1CDF2A.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/40F29C2D-FEEE-E111-A205-0025B3E0216C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/4051841B-15EF-E111-A37D-00266CFFCD50.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/4050C77C-A9EE-E111-B61F-AC162DABAF78.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/3E8C523D-F6EE-E111-A249-1CC1DE052030.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/3E5D7CC4-D8EE-E111-ADCC-AC162DA87230.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/3C991743-2AEF-E111-9FB2-00266CFF0608.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/3A178D03-00EF-E111-86B5-00266CFF0034.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/3A0415D7-7BEE-E111-9A75-00266CFFBF84.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/38AAC437-3DF1-E111-B682-0017A477102C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/368E2884-34EF-E111-8491-1CC1DE046F00.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/363994AF-F7EE-E111-B5A9-1CC1DE1CDF30.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/36214013-D8EE-E111-BD1B-0025B3E020D0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/3448B925-36EF-E111-87EF-00266CFFCC50.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/326CDDA8-E0EE-E111-AC53-00237DA1A548.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/2EBA0579-0CEF-E111-B425-00266CFEFDE0.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/2CD2538F-49F1-E111-868B-0017A4770C0C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/284262A6-71F1-E111-83F3-0017A4770C1C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/26E32228-5BEF-E111-9C2C-00266CFFBDAC.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/22984E34-96EE-E111-91BD-00266CFFBC60.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/2069460C-1DF1-E111-9B84-00266CFFBF38.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/203C8D20-F6EE-E111-B5FE-AC162DACC3E8.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/1C4B39FF-37EF-E111-855D-00266CFEFCE8.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/1C0E5A49-F6EE-E111-AB75-00266CFFBF90.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/1676DD35-8EEE-E111-9820-00266CFFCD60.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/14882B81-00F1-E111-931B-1CC1DE05D2F8.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/0CDE3810-3BEF-E111-859C-AC162DACB208.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/0897BF94-19F1-E111-B115-00266CFFCAC8.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/083870AF-77EF-E111-9F5E-0025B3E02292.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/0827BF65-1AEF-E111-83C8-1CC1DE0530F8.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/081F9B91-2EEF-E111-B580-00266CFFCB14.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/066AED18-CBEE-E111-A42E-AC162DACB208.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/04AE6C47-08EF-E111-B25D-0025B3E0216C.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/049BE89B-04EF-E111-8E77-00266CFFC9C4.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/02C89B02-F7EE-E111-801D-1CC1DE1CEAEE.root',
       '/store/data/Run2012C/MuEG/AOD/24Aug2012-v1/00000/00C6DE04-D9EE-E111-80B5-0017A4770010.root' ] );


secFiles.extend( [
               ] )


