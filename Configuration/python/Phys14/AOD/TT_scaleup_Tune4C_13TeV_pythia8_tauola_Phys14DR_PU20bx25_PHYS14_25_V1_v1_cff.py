import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/04B85240-2E6F-E411-B7CA-0025904B1026.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/04CFAAE0-466F-E411-9ACA-002590DB9152.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/04F00AF0-8D6F-E411-B7B9-0025901D4C18.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0CA5E64C-7A6F-E411-BBFA-002590AC4CCC.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0CDF014B-276F-E411-A965-002481E0DC7C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0CDF3185-716F-E411-81A0-002590AC4CC8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0E297521-506F-E411-BD1A-0025901D4936.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0E47C34F-886F-E411-BF56-002590AC505C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0EF532F1-416F-E411-8F80-0025907FD3CC.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/106DF197-716F-E411-8618-002590AC4C56.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/121CFEC7-606F-E411-AB75-0025901D493C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/127000E3-4E6F-E411-922B-00266CF2ABA8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/14B4E9F7-466F-E411-8E1C-002590A2CCE6.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/16AD3144-5D6F-E411-AE38-002590DB052A.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1866C87F-576F-E411-BB4F-00266CF33118.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/187A4FCD-466F-E411-93F3-00266CF32CD0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1A17099A-716F-E411-AD9F-00266CF33340.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1A32A475-0770-E411-B461-00266CFFA5C4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/202E6879-256F-E411-A642-0025904B0F96.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/242A8F90-286F-E411-BB9F-0025904E32F2.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/24A83AAD-326F-E411-B460-00266CF32BC4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/26C5DED6-5E6F-E411-AFC7-002590494C94.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2A982648-7E6F-E411-BBCE-002590DB9232.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2CD82042-E36F-E411-ABAF-0025901D481E.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2E2D106A-346F-E411-A63F-002590AC4B54.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2ECDB8B2-7A6F-E411-891E-003048D4DFBA.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3032C1C6-4C6F-E411-9224-002481E0DC82.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/342EA035-596F-E411-B58E-00259074B2BE.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3812CAC8-4A6F-E411-880E-002481E0D66C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3AE2C85E-7A6F-E411-AE0B-00266CFFA7A8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3ECFBB9C-716F-E411-90A1-002590AC4CCC.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4018C8AA-F170-E411-92C1-0025907DC9D6.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/40D2D3A5-2D6F-E411-BA03-0025907FD3CC.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/40D39401-426F-E411-ABF8-002590A52B4A.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4A012675-8B70-E411-9309-00266CF2E2C8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4AFC3559-A26F-E411-81E5-003048D4DFA4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4CBE820D-426F-E411-ACAB-00237DDC5E96.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/50551A88-716F-E411-A1FE-002590DB92C4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/50B9B5E0-876F-E411-AC67-0025907DC9D6.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5420B5D5-466F-E411-9C27-002590DB923E.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/54BD629C-2B6F-E411-94F0-002590AC4CCC.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/56B53C0B-576F-E411-902E-002590AC4CAA.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5A93CBEF-8D6F-E411-B8C4-0025901D47AA.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5AE7B553-376F-E411-A7DA-003048F0E82C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/62F80C51-516F-E411-AFDC-002590AB3A70.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/64FA6408-576F-E411-AC56-0025904B130A.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6C74660A-726F-E411-8ECE-00266CFFA754.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6CF8865C-FE6E-E411-A404-0025901D42C0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6EDC67B3-716F-E411-A123-003048F0E838.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/707DCB81-746F-E411-B68A-003048F0E780.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/72C2FA04-4E6F-E411-BBAB-002590AC4C7C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/76867C7C-CD6F-E411-9224-002590DBDFE0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/785BD7B5-716F-E411-9834-002590A36084.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/78DA6087-6A6F-E411-B6DB-002590AC4C48.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7A9E2C71-306F-E411-886A-00266CFFA120.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7C1EFB57-EE6F-E411-B548-003048D439C0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7C742F33-256F-E411-A958-002590AC4B5C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/820A19F0-416F-E411-AD49-00266CFFA344.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/846E164B-296F-E411-B438-002481E94B26.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/849B9B05-A86F-E411-B095-0025904B5FB8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/88AF22A2-356F-E411-AD6F-002481E7636A.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/942D5E0B-936F-E411-A051-002590AC4B5C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/96A8F9C6-5A6F-E411-8389-0025904B12FC.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9C106047-296F-E411-8029-002590DB91C6.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A64D7FED-D970-E411-8664-003048CF6334.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A84E45E9-F36E-E411-9DDC-0025904B578E.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A87F46B7-9D6F-E411-9BBA-003048F0E82A.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A8E95E73-CD6F-E411-8600-003048F0E5A4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AA0DC8E3-466F-E411-B302-002590AC4C7C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AABC8011-0670-E411-A867-00266CFFA184.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/ACB91B86-986F-E411-9F12-0025901D4C46.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B089D304-EE6E-E411-A005-002590DB91F0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B229A6BF-4A6F-E411-B205-0025907FD424.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B2E83B4B-296F-E411-BE67-002590DB9294.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B4369A1D-2B6F-E411-A0C0-0025904B0F96.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B67BDFD3-746F-E411-9EE4-002590DB924E.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B8E8D14A-7E6F-E411-9BDC-0025904B12A4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BEFED0AB-F170-E411-8F60-002481E0D66C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C02D2FA2-816F-E411-AD53-002590AC4C66.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C05528EC-986F-E411-8A8F-008CFA104E64.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C0D86106-776F-E411-918C-003048D4DEAC.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C0F9C3D7-5E6F-E411-93B6-002590AC4C56.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C47D9CE8-336F-E411-805C-0025901D4938.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C6BEEF9D-AB6F-E411-8880-002590494CBC.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C6D15EA8-836F-E411-87F8-00266CF32AF4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C80AB58B-376F-E411-B10D-002590DBDFE0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C84EAE0B-756F-E411-B42D-002590AC4E28.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C857A393-BC6F-E411-B875-0025907DC9C4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C89A2E2E-816F-E411-8767-0025907DBA06.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CA32BDD1-466F-E411-B802-002590AC4C80.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CECBB7FE-4A6F-E411-BA3C-0025907FD2B6.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D448C84F-936F-E411-B3EC-002590AC4C56.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DEBA1617-576F-E411-9998-003048D437E4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DED75DCF-B26F-E411-93C3-003048D437C4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E01B1CC4-5A6F-E411-9A90-00266CFFA7C0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E0246290-266F-E411-994E-002590A2CDA8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E272F936-436F-E411-83CD-003048D439A8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E4F2E2BB-9D6F-E411-9D29-003048D4DEA6.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E6972BA6-326F-E411-878D-002590AC4C52.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E820876D-4D6F-E411-A58D-003048D3CD6C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E824FFAB-6A6F-E411-A276-003048D462C8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E8D9480C-576F-E411-B97D-003048F0E836.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/EC3B6E21-A26F-E411-B1C8-002590AC4C7C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/EC5D98A8-386F-E411-8778-002590DB918C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F0E83063-426F-E411-94FB-003048F0E3BE.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F0EB2597-526F-E411-9B59-00266CF2718C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F2D1903F-E36F-E411-B13D-0025907FD2BA.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F43DD5D5-4A6F-E411-A218-003048F0E83C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F4646A3D-596F-E411-BDBD-00266CF327E0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F4987275-716F-E411-BB06-0025904B0FE4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F4AC4135-426F-E411-9D4A-00266CFFA750.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F6CC04C0-4A6F-E411-BEAB-003048F02CBA.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F8E4607D-716F-E411-BB20-002590DB9258.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FA06FD02-776F-E411-9AC9-002590AC5074.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FC60230D-B36F-E411-8A72-002590AC5074.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/00336D7D-B46F-E411-BF9C-00266CFFA184.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/022A8464-9F6F-E411-A7A1-0025901D4B06.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/0274F5D5-9A6F-E411-8D09-0025901D4932.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/04FD813D-7B6F-E411-98A8-003048F0E188.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/0624732D-766F-E411-839A-002590494C94.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/067EB46D-FE6F-E411-B320-00266CF327C4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/0895B3D7-316F-E411-A0C5-002481E0D500.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/0C682895-DC6F-E411-909F-00266CFFA754.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/0E3F8847-476F-E411-8DAA-002481E0D6EE.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/10033A15-3E6F-E411-88B2-00266CF2AE10.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/101A24FF-6B6F-E411-85C7-002590A7A2E0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/1278FC79-A46F-E411-82D1-00266CF33100.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/142C1DC2-216F-E411-9F70-002590AC4C08.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/14B70CD5-256F-E411-B611-0025904B1026.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/1668C25B-776F-E411-846F-0025904B141E.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/18869EA1-ED70-E411-9AEE-003048F0E526.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/1A5893DF-246F-E411-84D2-002590AC5082.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/1AA4DBF3-356F-E411-9521-0025904B12E0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/1AE3B3DB-316F-E411-AAD9-002590DB9214.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/1AECB7AE-336F-E411-87C7-002481E14F2A.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/1C0AADAF-A46F-E411-8B1C-002590AC4BF0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/1E99B452-A06F-E411-9C31-002590DB917C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/24E8126E-A96F-E411-9C7D-00266CFFA7BC.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/26E9507C-B46F-E411-8FE2-0025901D4D64.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/286EEC54-646F-E411-BE9E-003048D4DEAE.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/287612E5-F66E-E411-A0E1-003048D4DFA6.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/2A0BC184-576F-E411-80E8-00266CFFA7A8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/38840BA2-226F-E411-A89B-0025904B1026.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/3CCC1EDC-566F-E411-B8AF-002590A2CD44.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/3E7407DF-756F-E411-8726-003048F0E194.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/3EAF6DB2-426F-E411-A745-0025907DC9B4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/40DB6B39-356F-E411-BF4A-00266CF32930.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/40DF8BA3-CE6F-E411-91C1-002590A2CDA8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/40F8CE3F-956F-E411-9677-002590AC505C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/421DEE7E-246F-E411-A91B-0025907DCA64.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/4626B80F-7B6F-E411-BD2D-0025904B141E.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/4A48D859-726F-E411-93DA-0025901D4124.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/4CF88A36-816F-E411-84DD-002590DB05DA.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/52684D9C-426F-E411-A464-00266CF327C4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/52A72EE9-4D6F-E411-A48C-002590DB91C4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/561ADDBC-4F6F-E411-AC42-00266CFFB390.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/5673AC6F-5B6F-E411-8EB3-002590AC4CAA.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/56BA6FBA-296F-E411-BDA4-002590AC5082.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/5A8B284D-6D6F-E411-ACAA-002481E0D500.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/5E375C22-286F-E411-B384-00266CF326A8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/5E61279D-8F6F-E411-8973-00266CF32BC4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/60390276-576F-E411-9E7F-002590AC4BF8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/665665F9-6B6F-E411-A156-0025904B0FB6.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/68F123B7-4D6F-E411-B9B9-002590DBDFE0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/6AF0CC33-816F-E411-97FB-002590DB91C4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/6C23FCF1-8B6F-E411-988E-00266CFFA5E0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/6CE2CB5F-216F-E411-A3AB-00266CF32EAC.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/70E27047-726F-E411-B938-00266CFFA678.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/7AD444A0-8A6F-E411-86E8-002481E0D974.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/7C5AEBD6-3E6F-E411-89E5-0025907FD442.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/80167160-BA6F-E411-A9F5-0025901D4124.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/805E4858-946F-E411-B666-00266CF32CD0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/8291237F-216F-E411-AA9D-002590AC4CEC.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/84594137-7C6F-E411-8EBA-002590DB05F4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/8A704F47-DC6F-E411-9D34-0025907DBA06.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/8E21D227-C96F-E411-ABAB-003048F0E5B4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/8E273CCA-4D6F-E411-BF68-00266CF32F00.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/8E4BF6ED-646F-E411-82D6-003048D43996.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/90511C49-726F-E411-AAA2-00266CF32F14.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/90DF6540-766F-E411-9766-00266CFFA5C4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/922741E0-566F-E411-B899-0025907FD40C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/92900669-FE6F-E411-B4D1-00266CF32F14.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/96E13BA0-426F-E411-87A2-00266CF3338C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/9AEC5374-3D6F-E411-8A88-0025904B130A.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/9C29D9C6-8A6F-E411-AB3A-002590DB91CE.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/9E705705-7B6F-E411-BAC9-0025901D4932.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/9EE7725B-476F-E411-9967-0025907FD424.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/A226595B-2D6F-E411-A0E7-002590AC52CA.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/A4C57AFF-6B6F-E411-8849-002590494E36.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/A4E9BF4C-476F-E411-AC92-002590AC4B50.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/A61D4579-826F-E411-8557-003048D43960.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/AA7013A6-ED70-E411-A109-00266CF33100.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/AC6444FA-506F-E411-B44F-003048D439C0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/AEBBBCC3-426F-E411-B6B5-00266CFFA6A8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/B0DCC109-436F-E411-852B-002590A2CCE6.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/B0E4AF8E-9F6F-E411-A5B7-002590AC4E28.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/B2815D37-646F-E411-A2C2-003048D3CB3C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/B8B12208-BB6F-E411-8C43-002590494CB2.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/C2109555-ED6E-E411-A205-0025901D4D54.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/C2F495FB-A170-E411-9278-002481E0DBE0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/C42EBA61-646F-E411-A50F-003048D47976.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/C49903DB-9A6F-E411-9A22-0025901D4932.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/C6C1CD0C-576F-E411-8DAE-002590DB9188.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/CA2D0ACC-276F-E411-B354-003048D4DFB8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/CA416E98-426F-E411-AF7C-002590DB923E.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/CC9D0001-6C6F-E411-B0E8-0025901D4B08.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/CEA7E575-5B6F-E411-891F-0025904B12B2.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/D44FD3C2-A96F-E411-8886-003048F0E522.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/D4806B9D-386F-E411-B7DF-0025907DCA9C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/D4E01FEC-4D6F-E411-A4A1-002590DB05F4.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/DCED2908-226F-E411-8DE7-0025907DCA0C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/DE2AFA5C-726F-E411-9E26-00266CFFB7B8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/E09E3832-2F6F-E411-B244-002481E76372.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/E2C75110-7B6F-E411-BAB4-0025901D4C46.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/E411BEC1-4D6F-E411-B029-002481E0D6EE.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/E4EB5A37-2F6F-E411-978D-002590DB91C8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/E6B7D062-336F-E411-A4B1-002590DB91C8.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/EC87D7AB-8F6F-E411-BC09-002590AC4CEC.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/EE78F09A-8A6F-E411-9B5C-003048D4DF08.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/EEF4ADCC-426F-E411-89C2-002590A7A2E0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/F25B5547-476F-E411-8FF7-002590AC4BF0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/F2E950E8-2B6F-E411-9269-003048D47746.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/F41B1C11-816F-E411-BCFC-002590AC504A.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/F443F35B-C36F-E411-98E9-002481E76372.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/F6BD9272-F66E-E411-ABF2-003048D436C6.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/F88B12D9-AC70-E411-8BEB-002590AC5062.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/F890FBC6-726F-E411-9C95-002590DB91E0.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/FC1432AC-426F-E411-A143-00266CF2718C.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/FCBDD065-646F-E411-BF59-003048D4399E.root',
       
'/store/mc/Phys14DR/TT_scaleup_Tune4C_13TeV-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/FE2B97DD-CE6F-E411-8E43-002481E101DC.root' 
] );


secFiles.extend( [
               ] )

