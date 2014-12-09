import FWCore.ParameterSet.Config as cms

mumuTriggers = cms.vstring(
                            'HLT_Mu17_Mu8_v*',
                            'HLT_Mu17_TkMu8_v*',
                            'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*',
                            'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*',
                            'HLT_DoubleMu33NoFiltersNoVtx_v*',
                            )

emuTriggers  = cms.vstring(
                            'HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v*',
                            'HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v*',
                            )

eeTriggers   = cms.vstring(
                            'HLT_Ele23_Ele12_CaloId_TrackId_Iso_v*',
                            )




#mumuTriggers = cms.vstring(
                            #'HLT_Mu17_Mu8_v*',
                            #'HLT_Mu17_TkMu8_v*',
                            #)

#emuTriggers  = cms.vstring( 
                            #'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
                            #'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
                            #)

#eeTriggers   = cms.vstring( 
                            #'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
                            #)



#mumuTriggers = cms.vstring( 'HLT_DoubleMu6_v*',
                            #'HLT_DoubleMu7_v*',
                            #'HLT_Mu13_Mu8_v*',
                            #'HLT_Mu17_Mu8_v*',
                            #'HLT_DoubleMu45_v*'
			    #'HLT_Mu17_TkMu8_v*',	# data
                            #'HLT_Mu17_Mu8_v16',		# only for MC
			    #'HLT_Mu17_TkMu8_v9',	# only for MC
                            #)
#emuTriggers  = cms.vstring( 'HLT_Mu10_Ele10_CaloIdL_*',
                            #'HLT_Mu8_Ele17_CaloIdL_*',
                            #'HLT_Mu17_Ele8_CaloIdL_*',
                            #'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*',
                            #'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*',
                            #'HLT_Mu10_Ele10_CaloIdL_v*',
                            #'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',	# data
			    #'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',	# data
			    #'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6',	# MC with GT=START52_V9::All  #To be updated for CMSSW53X 
			    #'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6',	# MC with GT=START52_V9::All  #To be updated for CMSSW53X 
                            #)
#eeTriggers   = cms.vstring( 'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*',
                            #'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*',
                            #'HLT_DoubleEle45_CaloIdL_v*',
			    #'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',	# data
			    #'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17',	# MC with GT=START52_V9::All  #To be updated for CMSSW53X 
                            #)

#METTriggers  = cms.vstring( 'HLT_CentralJet80_MET100_v*',
#                            'HLT_CentralJet80_MET110_v*',
#                            'HLT_CentralJet80_MET160_v*',
#                            'HLT_CentralJet80_MET65_v*',
#                            'HLT_CentralJet80_MET80_v*',
#                            'HLT_CentralJet80_MET95_v*',
#                            'HLT_DiCentralJet20_MET100_HBHENoiseFiltered_v*',
#                            'HLT_DiCentralJet20_MET80_v*',
#                            'HLT_DiJet60_MET45_v*',
#                            'HLT_MET100_v*',
#                            'HLT_MET120_v*',
#                            'HLT_MET200_v*',
#                            'HLT_MET400_v*',
#                            'HLT_MET65_v*',
#                            )


