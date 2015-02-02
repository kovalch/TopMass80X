--- src/HiggsAnalysis.cc.txt	2015-02-02 12:21:10.000000000 +0100
+++ src/HiggsAnalysis.cc	2015-02-02 12:26:28.000000000 +0100
@@ -242,6 +242,217 @@
                                                weightGenerator, weightTopPt,
                                                trueLevelWeightNoPileup, trueLevelWeight);
     
+    // Reading necessary branches for generator level signal histograms (before the fill statement of the step0b)
+    if(this->isTopSignal() && additionalBjetMode_/100 > 0) {
+        // === FULL RECO OBJECT SELECTION === (can thus be used at each selection step)
+    
+        // Access reco objects, and common generator objects
+        const RecoObjects& recoObjects = this->getRecoObjects(entry);
+        const CommonGenObjects& commonGenObjects = this->getCommonGenObjects(entry);
+        
+        // Get allLepton indices, apply selection cuts and order them by pt (beginning with the highest value)
+        const VLV& allLeptons = *recoObjects.allLeptons_;
+        const std::vector<int>& lepPdgId = *recoObjects.lepPdgId_;
+        std::vector<int> allLeptonIndices = initialiseIndices(allLeptons);
+        selectIndices(allLeptonIndices, allLeptons, LVeta, LeptonEtaCUT, false);
+        selectIndices(allLeptonIndices, allLeptons, LVeta, -LeptonEtaCUT);
+        selectIndices(allLeptonIndices, allLeptons, LVpt, LeptonPtCut);
+        orderIndices(allLeptonIndices, allLeptons, LVpt);
+        //const int numberOfAllLeptons = allLeptonIndices.size();
+        
+        // Get indices of leptons and antiLeptons separated by charge, and get the leading ones if they exist
+        std::vector<int> leptonIndices = allLeptonIndices;
+        std::vector<int> antiLeptonIndices = allLeptonIndices;
+        selectIndices(leptonIndices, lepPdgId, 0);
+        selectIndices(antiLeptonIndices, lepPdgId, 0, false);
+        const int numberOfLeptons = leptonIndices.size();
+        const int numberOfAntiLeptons = antiLeptonIndices.size();
+        const int leptonIndex = numberOfLeptons>0 ? leptonIndices.at(0) : -1;
+        const int antiLeptonIndex = numberOfAntiLeptons>0 ? antiLeptonIndices.at(0) : -1;
+        
+        // In case of an existing opposite-charge dilepton system,
+        // get their indices for leading and next-to-leading lepton
+        int leadingLeptonIndex(-1);
+        int nLeadingLeptonIndex(-1);
+        if(numberOfLeptons>0 && numberOfAntiLeptons>0){
+            leadingLeptonIndex = leptonIndex;
+            nLeadingLeptonIndex = antiLeptonIndex;
+            orderIndices(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, LVpt);
+        }
+        const bool hasLeptonPair = this->hasLeptonPair(leadingLeptonIndex, nLeadingLeptonIndex, lepPdgId);
+        
+        // Get two indices of the two leptons in the right order for trigger scale factor, if existing
+        int leptonXIndex(leadingLeptonIndex);
+        int leptonYIndex(nLeadingLeptonIndex);
+        if(hasLeptonPair){
+            //in ee and mumu channel leptonX must be the highest pt lepton, i.e. this is already correct
+            // in emu channel leptonX must be electron
+            if(std::abs(lepPdgId.at(leptonXIndex)) != std::abs(lepPdgId.at(leptonYIndex))){
+                orderIndices(leptonYIndex, leptonXIndex, lepPdgId, true);
+            }
+        }
+        
+        // Get dilepton system, if existing
+        const LV dummyLV(0.,0.,0.,0.);
+        const LV dilepton(hasLeptonPair ? allLeptons.at(leadingLeptonIndex)+allLeptons.at(nLeadingLeptonIndex) : dummyLV);
+        
+        // Get jet indices, apply selection cuts and order them by pt (beginning with the highest value)
+        const VLV& jets = *recoObjects.jets_;
+        std::vector<int> jetIndices = initialiseIndices(jets);
+        selectIndices(jetIndices, jets, LVeta, JetEtaCUT, false);
+        selectIndices(jetIndices, jets, LVeta, -JetEtaCUT);
+        selectIndices(jetIndices, jets, LVpt, JetPtCUT);
+        orderIndices(jetIndices, jets, LVpt);
+        const int numberOfJets = jetIndices.size();
+        
+        // Fill a vector with all jet pair indices, while sorting each pair by the jet charge:
+        // first entry is antiBIndex i.e. with higher jet charge, second entry is bIndex
+        //const std::vector<double>& jetChargeGlobalPtWeighted = *recoObjects.jetChargeGlobalPtWeighted_;
+        const std::vector<double>& jetChargeRelativePtWeighted = *recoObjects.jetChargeRelativePtWeighted_;
+        const tth::IndexPairs& jetIndexPairs = this->chargeOrderedJetPairIndices(jetIndices, jetChargeRelativePtWeighted);
+        
+        // Get b-jet indices, apply selection cuts
+        // and apply b-tag efficiency MC correction using random number based tag flipping (if requested correction mode is applied)
+        // and order b-jets by btag discriminator (beginning with the highest value)
+        const std::vector<double>& jetBTagCSV = *recoObjects.jetBTagCSV_;
+        const std::vector<int>& jetPartonFlavour = *commonGenObjects.jetPartonFlavour_;
+        std::vector<int> bjetIndices = jetIndices;
+        selectIndices(bjetIndices, jetBTagCSV, this->btagCutValue());
+        this->retagJets(bjetIndices, jetIndices, jets, jetPartonFlavour, jetBTagCSV);
+        orderIndices(bjetIndices, jetBTagCSV);
+        
+        // Get MET, and in case of MVA MET apply recoil correction for Drell-Yan sample
+        this->correctMvaMet(dilepton, numberOfJets, entry);
+        
+        const tth::RecoObjectIndices recoObjectIndices(allLeptonIndices,
+                                                       leptonIndices, antiLeptonIndices,
+                                                       leptonIndex, antiLeptonIndex,
+                                                       leadingLeptonIndex, nLeadingLeptonIndex,
+                                                       leptonXIndex, leptonYIndex,
+                                                       jetIndices, jetIndexPairs,
+                                                       bjetIndices);
+
+        
+        
+        // === FULL GEN OBJECT SELECTION ===
+    
+        // Access top generator object struct, and higgs generator object struct
+        const TopGenObjects& topGenObjects = this->getTopGenObjects(entry);
+        const HiggsGenObjects& higgsGenObjects = this->getHiggsGenObjects(entry);
+        
+        // Generated jets
+        const VLV& allGenJets =  topGenObjects.valuesSet_ ? *topGenObjects.allGenJets_ : VLV();
+        std::vector<int> allGenJetIndices = initialiseIndices(allGenJets);
+        std::vector<int> genJetIndices = allGenJetIndices;
+        selectIndices(genJetIndices, allGenJets, LVeta, GenJetEtaCUT, false);
+        selectIndices(genJetIndices, allGenJets, LVeta, -GenJetEtaCUT);
+        selectIndices(genJetIndices, allGenJets, LVpt, GenJetPtCUT);
+        if(GenDeltaRLeptonJetCUT > 0.){
+            // Vector of genLeptons from which genJets need to be separated in deltaR
+            VLV allGenLeptons;
+            if(topGenObjects.valuesSet_){
+                if(topGenObjects.GenLepton_) allGenLeptons.push_back(*topGenObjects.GenLepton_);
+                if(topGenObjects.GenAntiLepton_) allGenLeptons.push_back(*topGenObjects.GenAntiLepton_);
+            }
+            this->leptonCleanedJetIndices(genJetIndices, allGenJets, allGenLeptons, GenDeltaRLeptonJetCUT);
+        }
+        
+        // Match for all genJets all B hadrons
+        std::vector<std::vector<int> > allGenJetBhadronIndices;
+        std::vector<std::vector<int> > genJetBhadronIndices;
+        std::vector<int> allGenBjetIndices;
+        std::vector<int> genBjetIndices;
+        std::vector<int> genJetMatchedRecoBjetIndices;
+        if(topGenObjects.valuesSet_){
+            allGenJetBhadronIndices = this->matchHadronsToGenJets(allGenJetIndices, allGenJets, *topGenObjects.genBHadJetIndex_);
+            genJetBhadronIndices = this->matchHadronsToGenJets(genJetIndices, allGenJets, *topGenObjects.genBHadJetIndex_);
+            allGenBjetIndices = this->genBjetIndices(allGenJetBhadronIndices);
+            genBjetIndices = this->genBjetIndices(genJetBhadronIndices);
+            genJetMatchedRecoBjetIndices = this->matchRecoToGenJets(jetIndices, jets, allGenBjetIndices, allGenJets);
+        }
+        
+        // Match for all genJets all C hadrons
+        std::vector<std::vector<int> > allGenJetChadronIndices;
+        std::vector<std::vector<int> > genJetChadronIndices;
+        std::vector<int> allGenCjetIndices;
+        std::vector<int> genCjetIndices;
+        std::vector<int> genJetMatchedRecoCjetIndices;
+        if(topGenObjects.valuesSet_){
+            allGenJetChadronIndices = this->matchHadronsToGenJets(allGenJetIndices, allGenJets, *topGenObjects.genCHadJetIndex_);
+            genJetChadronIndices = this->matchHadronsToGenJets(genJetIndices, allGenJets, *topGenObjects.genCHadJetIndex_);
+            allGenCjetIndices = this->genCjetIndices(allGenJetBhadronIndices, allGenJetChadronIndices);
+            genCjetIndices = this->genCjetIndices(genJetBhadronIndices, genJetChadronIndices);
+            genJetMatchedRecoCjetIndices = this->matchRecoToGenJets(jetIndices, jets, allGenCjetIndices, allGenJets);
+        }
+        
+        // Jet matchings for ttbar system
+        int genBjetFromTopIndex(-1);
+        int genAntiBjetFromTopIndex(-1);
+        int matchedBjetFromTopIndex(-1);
+        int matchedAntiBjetFromTopIndex(-1);
+        if(topGenObjects.valuesSet_){
+            genBjetFromTopIndex = this->genBjetIndex(topGenObjects, 6);
+            genAntiBjetFromTopIndex = this->genBjetIndex(topGenObjects, -6);
+//             matchedBjetFromTopIndex = genBjetFromTopIndex>=0 ? genJetMatchedRecoBjetIndices.at(genBjetFromTopIndex) : -1;
+//             matchedAntiBjetFromTopIndex = genAntiBjetFromTopIndex>=0 ? genJetMatchedRecoBjetIndices.at(genAntiBjetFromTopIndex) : -1;
+        }
+        
+        // Jet matchings for Higgs system
+        int genBjetFromHiggsIndex(-1);
+        int genAntiBjetFromHiggsIndex(-1);
+        int matchedBjetFromHiggsIndex(-1);
+        int matchedAntiBjetFromHiggsIndex(-1);
+        if(topGenObjects.valuesSet_ && higgsDecayMode == 5){
+            genBjetFromHiggsIndex = this->genBjetIndex(topGenObjects, 25);
+            genAntiBjetFromHiggsIndex = this->genBjetIndex(topGenObjects, -25);
+            matchedBjetFromHiggsIndex = genBjetFromHiggsIndex>=0 ? genJetMatchedRecoBjetIndices.at(genBjetFromHiggsIndex) : -1;
+            matchedAntiBjetFromHiggsIndex = genAntiBjetFromHiggsIndex>=0 ? genJetMatchedRecoBjetIndices.at(genAntiBjetFromHiggsIndex) : -1;
+        }
+        
+        
+        const tth::GenObjectIndices genObjectIndices(genJetIndices,
+                                                     allGenBjetIndices,
+                                                     genBjetIndices,
+                                                     genJetBhadronIndices,
+                                                     genJetMatchedRecoBjetIndices,
+                                                     allGenCjetIndices,
+                                                     genCjetIndices,
+                                                     genJetChadronIndices,
+                                                     genJetMatchedRecoCjetIndices,
+                                                     genBjetFromTopIndex, genAntiBjetFromTopIndex,
+                                                     matchedBjetFromTopIndex, matchedAntiBjetFromTopIndex,
+                                                     genBjetFromHiggsIndex, genAntiBjetFromHiggsIndex,
+                                                     matchedBjetFromHiggsIndex, matchedAntiBjetFromHiggsIndex);
+        
+        // Create dummies for objects, non-dummies are created only as soon as needed
+        const EventMetadata eventMetadataDummy;
+        const RecoObjects recoObjectsDummy;
+        const CommonGenObjects commonGenObjectsDummy;
+        const TopGenObjects topGenObjectsDummy;
+        const HiggsGenObjects higgsGenObjectsDummy;
+        const KinematicReconstructionSolutions kinematicReconstructionSolutionsDummy;
+
+        // Set up dummies for weights and indices, as needed for generic functions
+        const tth::GenObjectIndices genObjectIndicesDummy({}, {}, {}, {}, {}, {}, {}, {}, {}, -1, -1, -1, -1, -1, -1, -1, -1);
+        const tth::RecoObjectIndices recoObjectIndicesDummy({}, {}, {}, -1, -1, -1, -1, -1, -1, {}, {}, {});
+        const tth::GenLevelWeights genLevelWeightsDummy(0., 0., 0., 0., 0., 0.);
+        const tth::RecoLevelWeights recoLevelWeightsDummy(0., 0., 0., 0., 0., 0.);
+        
+        // Applying the reweighting
+        const double trueLevelWeight_reweighted = trueLevelWeight*reweightingWeight(topGenObjects, genObjectIndices);
+        
+        // ++++ Control Plots ++++
+        
+        this->fillAll(selectionStep,
+                      eventMetadataDummy,
+                      recoObjectsDummy, commonGenObjects,
+                      topGenObjects, higgsGenObjects,
+                      kinematicReconstructionSolutionsDummy,
+                      genObjectIndices, recoObjectIndicesDummy,
+                      genLevelWeights, recoLevelWeightsDummy,
+                      trueLevelWeight_reweighted);
+    } else
+    
     // ++++ Control Plots ++++
     
     this->fillAll(selectionStep,
