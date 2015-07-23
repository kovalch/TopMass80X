import FWCore.ParameterSet.Config as cms

## created GenJet collection with selected jets

from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu

genParticlesForJets = cms.EDFilter("CandViewSelector",
                       src = cms.InputTag("genParticles"),
                       cut = cms.string("((abs(pdgId) = 11 | abs(pdgId) = 13) & pt > 33 & abs(eta) < 2.5)")
                       )

newAk5GenJets = ak5GenJets.clone(src = "genParticlesForJetsNoNu")

selectedAk5GenJets = cms.EDFilter("GenJetSelector",
                                src = cms.InputTag("newAk5GenJets"),
                                cut = cms.string("pt > 20 & abs(eta) < 2.4"),
                            )

myAk5GenJets = cms.EDProducer("MyGenJetsProducer",
    genJets = cms.InputTag("selectedAk5GenJets"),
    genParticles = cms.InputTag("genParticlesForJets")
  )

inputGenJetCollection = "myAk5GenJets"


## convert GenJets to PatJets
from TopAnalysis.TopUtils.MyBTagProducer_cfi import createMyBTags
from TopAnalysis.TopUtils.MyPatJetCorrFactorsProducer_cfi import createMyPatJetCorrFactors
from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import patJets
from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import patJets
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets as mySelectedPatJets

myPatJetCorrFactors = createMyPatJetCorrFactors.clone(jetSrc = inputGenJetCollection)

combinedSecondaryVertexBJetTags      = createMyBTags.clone(jetSrc = inputGenJetCollection, bDiscValue = 1.0)
combinedSecondaryVertexMVABJetTags   = createMyBTags.clone(jetSrc = inputGenJetCollection, bDiscValue = 1.0)
jetProbabilityBJetTags               = createMyBTags.clone(jetSrc = inputGenJetCollection, bDiscValue = 1.0)
jetBProbabilityBJetTags              = createMyBTags.clone(jetSrc = inputGenJetCollection, bDiscValue = 5.0)
simpleSecondaryVertexHighEffBJetTags = createMyBTags.clone(jetSrc = inputGenJetCollection, bDiscValue = 5.0)
simpleSecondaryVertexHighPurBJetTags = createMyBTags.clone(jetSrc = inputGenJetCollection, bDiscValue = 5.0)
#softElectronByPtBJetTags             = createMyBTags.clone(jetSrc = inputGenJetCollection, bDiscValue = 1.0)
#softElectronByIP3dBJetTags           = createMyBTags.clone(jetSrc = inputGenJetCollection, bDiscValue = 1.0)
#softMuonBJetTags                     = createMyBTags.clone(jetSrc = inputGenJetCollection, bDiscValue = 1.0)
#softMuonByPtBJetTags                 = createMyBTags.clone(jetSrc = inputGenJetCollection, bDiscValue = 1.0)
#softMuonByIP3dBJetTags               = createMyBTags.clone(jetSrc = inputGenJetCollection, bDiscValue = 1.0)
trackCountingHighEffBJetTags         = createMyBTags.clone(jetSrc = inputGenJetCollection, bDiscValue = 20.0)
trackCountingHighPurBJetTags         = createMyBTags.clone(jetSrc = inputGenJetCollection, bDiscValue = 20.0)

myBTaggingSequence = cms.Sequence( combinedSecondaryVertexBJetTags
                                   * combinedSecondaryVertexMVABJetTags  
                                   * jetBProbabilityBJetTags             
                                   * jetProbabilityBJetTags              
                                   * simpleSecondaryVertexHighEffBJetTags
                                   * simpleSecondaryVertexHighPurBJetTags
                                   #* softElectronByPtBJetTags            
                                   #* softElectronByIP3dBJetTags          
                                   #* softMuonBJetTags                    
                                   #* softMuonByPtBJetTags                
                                   #* softMuonByIP3dBJetTags              
                                   * trackCountingHighEffBJetTags        
                                   * trackCountingHighPurBJetTags        
                                   )

myPatJets = patJets.clone(jetSource = inputGenJetCollection,
                                #embedCaloTowers = False,
                                embedPFCandidates = False, 
                                addJetCorrFactors    = True,
                                jetCorrFactorsSource = cms.VInputTag(cms.InputTag("myPatJetCorrFactors") ),
                                addBTagInfo          = True,
                                addDiscriminators    = True,
                                discriminatorSources = cms.VInputTag(cms.InputTag("combinedSecondaryVertexBJetTags"),
                                                                     cms.InputTag("combinedSecondaryVertexMVABJetTags"),     
                                                                     cms.InputTag("jetBProbabilityBJetTags"),                
                                                                     cms.InputTag("jetProbabilityBJetTags"),                 
                                                                     cms.InputTag("simpleSecondaryVertexHighEffBJetTags"),   
                                                                     cms.InputTag("simpleSecondaryVertexHighPurBJetTags"),   
                                                                     #cms.InputTag("softElectronByPtBJetTags"),               
                                                                     #cms.InputTag("softElectronByIP3dBJetTags"),             
                                                                     #cms.InputTag("softMuonBJetTags"),                       
                                                                     #cms.InputTag("softMuonByPtBJetTags"),                   
                                                                     #cms.InputTag("softMuonByIP3dBJetTags"),                 
                                                                     cms.InputTag("trackCountingHighEffBJetTags"),           
                                                                     cms.InputTag("trackCountingHighPurBJetTags"),           
                                                                     ),                                          
                                addTagInfos         = False,
                                addAssociatedTracks = False,
                                addJetCharge        = False,
                                addJetID            = False,
                                addGenPartonMatch   = False,
                                embedGenPartonMatch = False,
                                addGenJetMatch      = False,
                                embedGenJetMatch    = False,
                                addPartonJetMatch   = False,
                                getJetMCFlavour     = False,
                                addEfficiencies     = False,
                                addResolutions      = False
                                )

selectedPatJets = mySelectedPatJets.clone(src = 'myPatJets',
                                        cut = ''
                                        )
                                     
## process everything
convertGenToPatJets = cms.Sequence(genParticlesForJets*
                                   genParticlesForJetsNoNu*
                                   newAk5GenJets*
                                   selectedAk5GenJets*
                                   myAk5GenJets *
                                   #mygenMuons *
                                   myPatJetCorrFactors *
                                   myBTaggingSequence *
                                   selectedPatJets
                                   )

gen4JetSelector = cms.EDFilter("CandViewSelector",
                               src = cms.InputTag("ak5GenJets"),
                               cut = cms.string("pt > 60 & abs(eta) < 2.4")
                               )

gen4JetFilter   = cms.EDFilter("CandViewCountFilter",
                               src = cms.InputTag("gen4JetSelector"),
                               minNumber = cms.uint32(4)
                               )

gen5JetSelector = cms.EDFilter("CandViewSelector",
                               src = cms.InputTag("ak5GenJets"),
                               cut = cms.string("pt > 50 & abs(eta) < 2.4")
                               )

gen5JetFilter   = cms.EDFilter("CandViewCountFilter",
                               src = cms.InputTag("gen5JetSelector"),
                               minNumber = cms.uint32(5)
                               )

gen6JetSelector = cms.EDFilter("CandViewSelector",
                               src = cms.InputTag("ak5GenJets"),
                               cut = cms.string("pt > 40 & abs(eta) < 2.4")
                               )

gen6JetFilter   = cms.EDFilter("CandViewCountFilter",
                               src = cms.InputTag("gen6JetSelector"),
                               minNumber = cms.uint32(6)
                               )

fullHadGenJetSelection =  cms.Sequence(gen4JetSelector*gen4JetFilter*
                                       gen5JetSelector*gen5JetFilter*
                                       gen6JetSelector*gen6JetFilter
                                       )
