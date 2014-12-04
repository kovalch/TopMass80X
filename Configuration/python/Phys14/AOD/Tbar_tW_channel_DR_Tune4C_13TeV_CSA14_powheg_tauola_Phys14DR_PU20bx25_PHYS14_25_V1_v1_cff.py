import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/000179EE-B771-E411-97E4-00259074AE28.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/004EC757-1272-E411-BC24-002590747D94.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/04D8B2B2-A771-E411-BB91-E0CB4E5536BB.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/04E8B657-0772-E411-AA61-E0CB4E29C4C1.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0892D3E6-BC71-E411-95A4-002590D0B054.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/08BB6F64-BB71-E411-8657-00259073E4DA.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0A8B6EA4-AA71-E411-B37F-002590D0B06E.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0AAD543C-F171-E411-A32F-E0CB4E29C4E6.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0E8EFBC2-2E72-E411-8A78-002590D0B0C4.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/14988092-3372-E411-B0A6-00259073E3AC.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/16CEAA79-2B72-E411-B63D-002590D0B054.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/18A463AA-F671-E411-9800-00259073E370.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/18D3C748-5B72-E411-B223-002590D0B096.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1E700417-3072-E411-B9F2-20CF3027A5C9.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1E7AFB90-A671-E411-92A9-20CF305B04D1.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/247BDEC3-3872-E411-99A2-E0CB4E4408F9.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/28CBEF45-DE71-E411-B2C7-00259074AE38.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2A8E2426-3672-E411-8BC9-00259073E4E6.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2AFB2D7B-5672-E411-AA22-002590747D90.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2EBCC076-C771-E411-99A3-00259073E536.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/300D8DE4-A871-E411-86D4-00259073E506.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/32D6E9AB-BE71-E411-95A6-20CF3027A5A6.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/34CF84C9-9771-E411-9B14-0025907277BE.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/369682F5-4E72-E411-B01C-E0CB4EA0A917.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3A1902ED-D371-E411-A035-00259073E3A8.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3AFCFFBC-A571-E411-B2B0-485B3989720C.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/40154F7D-9071-E411-A96B-00259074AE3C.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/444AA4C2-9971-E411-8B7A-0025907B4E6C.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/44595495-E971-E411-B3D6-002590D0AF76.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/46991524-E471-E411-9466-002590D0B04E.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/48B5AA7E-B871-E411-B48B-E0CB4E29C511.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4AD1C751-B971-E411-A167-20CF305B0512.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4C9A29D5-3F72-E411-8E6B-002590D0B0C4.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4EC5C814-AD71-E411-A459-002590D0B0BA.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/506BCAD2-3A72-E411-9252-E0CB4E19F9B4.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/523D5CCF-E271-E411-9431-E0CB4E19F961.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/54ED692F-9671-E411-B008-00259073E474.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/580AFD08-EA71-E411-B25C-485B39800C08.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/58310C7E-F371-E411-A880-E0CB4E5536BE.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5A452192-EB71-E411-B89A-002590D0B000.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5C968794-9E71-E411-93D8-E0CB4E55364C.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5CBAD5D2-A271-E411-AFB4-E0CB4E4408EE.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6061FB69-E071-E411-93D8-00259074AE38.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6282FD90-2A72-E411-BF1E-20CF3027A5C9.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/64779468-A771-E411-87DF-20CF305616EC.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/64C79D87-ED71-E411-907C-002590D0B03A.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6631F5D7-5B72-E411-8468-00259074AE3C.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/663F84ED-A471-E411-8BA2-E0CB4EA0A8E0.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/66C11F1A-D371-E411-9C18-00259073E466.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/68343433-F971-E411-B307-002590D0AFC8.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/689AF987-FB71-E411-A891-485B39800C3B.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/68DA179A-E071-E411-B18E-E0CB4E19F976.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7029CB37-B671-E411-A7A2-0025907B4F16.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7271086D-D171-E411-A161-00259073E388.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/74BB7BBA-E771-E411-8585-0025907B4F1C.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7C0AB575-AB71-E411-A1EB-002590D0B0BA.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/802BDCFC-8C71-E411-BFD5-00259074AE3C.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8264EC9F-B771-E411-86D7-0025907B4F44.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/82F7C4BE-B471-E411-9B3B-485B3989720C.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/869F0BC3-EF71-E411-BCA1-002590D0B050.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8A39EAB9-B871-E411-AC63-E0CB4EA0A936.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8C8930ED-0072-E411-AAF6-002590D0B02A.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8E332B19-9371-E411-83AE-E0CB4E0ED8C0.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8EFEFC13-2D72-E411-AF8C-00259073E4F6.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/925F58D6-5E72-E411-BCB7-0025907277A0.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/96827614-AF71-E411-8C08-E0CB4E19F987.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A6282BA1-A871-E411-8286-002590D0B0B6.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A679546C-B371-E411-B2C2-002590D6009C.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A6928882-F571-E411-9706-E0CB4E29C4F6.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AE1BB86E-DF71-E411-BC1A-E0CB4EA0A8FA.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AEBE8D1B-BA71-E411-9FDB-00259074AE28.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B00E1C1C-B671-E411-A220-002590D0B09C.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B0581C82-6972-E411-968C-485B39897219.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B2D0A47C-2172-E411-9B4E-00259073E53C.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B409377D-5872-E411-AB53-002590747D90.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B8230AC0-B171-E411-9413-00261834B586.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BC50219F-DF71-E411-A495-E0CB4E19F961.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BCE801D2-C071-E411-841B-00259073E35A.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C024BD7E-9B71-E411-8D39-002590D0B09E.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C03CB71E-FD71-E411-9E98-485B39800C3B.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C2F43D55-E571-E411-B1C9-00259073E35A.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C4E1DC8D-8871-E411-B842-00259074AE3C.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C83C0D30-CA71-E411-89B3-00259073E3FA.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CA6FCE7A-2972-E411-9FB0-485B39800BBF.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CAF35970-C571-E411-8982-0025907B4FBC.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CE6443D7-3172-E411-8338-90E6BA19A214.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D80A9CAA-9C71-E411-895E-00259074AE38.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E2B2D75F-C271-E411-B699-E0CB4E55366A.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E689BA75-DF71-E411-A644-00259022277E.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/EEEB98F4-A971-E411-8EB0-002590D0B07A.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F29E0243-D171-E411-AF6B-002590747D94.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F6134D44-B271-E411-966B-002590D0B0B4.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FC00EB9C-E171-E411-BC79-E0CB4EA0A8FA.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/08DFB949-9A71-E411-96CD-20CF3027A598.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/1E4F54E1-9671-E411-A513-E0CB4E5536F2.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/5AE30939-9B71-E411-B5B8-002590D0AFAC.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/64641E15-9D71-E411-B176-002590D0B056.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/6A0C6609-A971-E411-96A8-0025907B4FD6.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/6A78545C-AB71-E411-8BF6-0025907B4F72.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/6C9A410E-9871-E411-9DD0-E0CB4E553665.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/6C9F870B-A071-E411-B4FB-E0CB4E1A114B.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/8EADB99F-B371-E411-9BDD-00259073E398.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/ACB35EB0-A771-E411-82F6-00259073E51A.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/C04459AC-9D71-E411-A2D4-002590D0B052.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/CEF6CF41-9A71-E411-9FCE-002590D0B014.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/D4BBDF44-9971-E411-9BBA-00259074AE3C.root',
       
'/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/F6CC7F71-8571-E411-B924-00259073E41E.root' 
] );


secFiles.extend( [
               ] )

