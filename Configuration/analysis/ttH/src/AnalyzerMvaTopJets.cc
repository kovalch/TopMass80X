#include <iostream>

#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TString.h>
#include <TFile.h>
#include <TList.h>
#include <TKey.h>
#include <TObjArray.h>
#include <TObjString.h>

#include "AnalyzerMvaTopJets.h"
#include "MvaReaderBase.h"
#include "MvaReaderTopJets.h"
#include "MvaVariablesTopJets.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
#include "../../common/include/RootFileReader.h"
#include "../../common/include/KinematicReconstructionSolution.h"





AnalyzerMvaTopJets::AnalyzerMvaTopJets(const char* mva2dWeightsFile,
                                       const std::vector<TString>& selectionStepsNoCategories,
                                       const std::vector<TString>& stepsForCategories,
                                       const JetCategories* jetCategories):
AnalyzerBase("mvaA_", selectionStepsNoCategories, stepsForCategories, jetCategories)
{
    std::cout<<"--- Beginning setting up MVA validation\n";
    
    // Open the root file
    TFile* weightsFile = new TFile(mva2dWeightsFile, "READ");
    // Generate the list of all objects stored in the root file, to find steps and training names
    TList* list = weightsFile->GetListOfKeys();
    for(int keyNum = 0; keyNum < list->GetSize(); ++keyNum){
        TKey* key = (TKey*)list->At(keyNum);
        
        TString objType(key->GetClassName());
        if(!objType.EqualTo("TObjArray")) continue;
        
        TString objName(key->GetName());
        if(!objName.BeginsWith("correct_")) continue;
        //if(!objName.Contains("cate0")) continue;
        
        // Extract name of the step (substring from 8-th character to the end)
        // Synchronized with the length of "correct_"
        std::string stepName(objName.Data(), 8, objName.Length());
        
        TObjArray* trainingsList;
        // Extract names of all correct trainings for the current step
        std::vector<std::string> v_correct;
        trainingsList = (TObjArray*)(weightsFile->Get(TString("correct_").Append(stepName))->Clone());
        for(int trainId = 0; trainId < trainingsList->GetEntries(); ++trainId){
            v_correct.push_back(std::string(((TObjString*)trainingsList->At(trainId))->String()));
        }
        
        // Extract names of all swapped trainings for the current step
        std::vector<std::string> v_swapped;
        trainingsList = (TObjArray*)(weightsFile->Get(TString("swapped_").Append(stepName))->Clone());
        for(int trainId = 0; trainId < trainingsList->GetEntries(); ++trainId){
            v_swapped.push_back(std::string(((TObjString*)trainingsList->At(trainId))->String()));
        }
        
        // Create struct of correct and swapped weights for the current step
        MvaWeightsStruct mvaWeightsStruct = MvaWeightsStruct(stepName, v_correct, v_swapped, mva2dWeightsFile);
        v_mvaWeightsStruct_.push_back(mvaWeightsStruct);
    }

    weightsFile->Close();

    std::cout<<"=== Finishing setting up MVA validation\n\n";
}



void AnalyzerMvaTopJets::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    for(const MvaWeightsStruct& mvaWeightsStruct : v_mvaWeightsStruct_){
        // FIXME: treatment of different MVA sets neglected for now, could be treated via folders in root file?
        //const std::string& mvaSetStep = mvaWeightsStruct.stepName();
        
        this->bookHistosPerSet(step, m_histogram, mvaWeightsStruct);
    }
}



void AnalyzerMvaTopJets::bookHistosPerSet(const TString& step, std::map<TString, TH1*>& m_histogram, const MvaWeightsStruct& mvaWeightsStruct)
{
    TString name;
    
    
    name = "dijet_mass";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "dijetMass;m_{jj} [GeV];# jet pairs", 100, 0, 500));
    
    // Book histograms for MVA weights for correct trainings
    for(const auto& weights : mvaWeightsStruct.correctWeights()){
        const std::string& name = weights.first;
        this->bookMvaSpecificHistos(step, m_histogram, "_" + name, "correct");
        
        // Book histograms for MVA weights for combined trainings
        if(mvaWeightsStruct.combinedWeights().find(name) != mvaWeightsStruct.combinedWeights().end()){
            const auto& weights1 = mvaWeightsStruct.combinedWeights().at(name);
            for(const auto weights2 : weights1){
                const std::string& name2 = weights2.first;
                const std::string& nameCombined = name + "_" + name2;
                this->bookMvaSpecificHistos(step, m_histogram, "_" + nameCombined, "combined");
            }
        }
    }
    
    // Book histograms for MVA weights for swapped trainings
    for(const auto& weights : mvaWeightsStruct.swappedWeights()){
        const std::string& name = weights.first;
        this->bookMvaSpecificHistos(step, m_histogram, "_" + name, "swapped");
    }
}



void AnalyzerMvaTopJets::bookMvaSpecificHistos(const TString& step, std::map<TString, TH1*>& m_histogram,
                                          const TString& mvaConfigName, const TString& mvaType)
{
    TString mvaTypeC;
    if(mvaType == "correct") mvaTypeC = "Correct";
    else if(mvaType == "swapped") mvaTypeC = "Swapped";
    else if(mvaType == "combined") mvaTypeC = "Combined";
    else{
        std::cerr<<"Error in bookMvaSpecificHistos()! Undefined mvaType: "<<mvaType<<"\n...break\n"<<std::endl;
        exit(439);
    }
    
    TString name;
    
    
    name = "dijet_mass_best"+mvaTypeC+mvaConfigName;
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "dijetMass_{"+mvaType+"}^{best};m_{jj (best)}^{"+mvaType+"} [GeV];# events", 100, 0, 500));
    
    
    name = "jetsFromTop_best"+mvaTypeC+mvaConfigName;
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "jetsFromTop_{"+mvaType+"}^{best};;# events",4,0,4));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "matched Top jets");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "pair from Top");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "correct pair");
    m_histogram[name]->GetXaxis()->SetBinLabel(4, "swapped pair");
    
    
    name = "jetsFromHiggs_best"+mvaTypeC+mvaConfigName;
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "jetsFromHiggs_{"+mvaType+"}^{best};;# events",4,0,4));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "matched Higgs jets");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "pair from Higgs");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "overlap with Top");
    m_histogram[name]->GetXaxis()->SetBinLabel(4, "fully solved");
    
    
    name = "jetsFromBoth_best"+mvaTypeC+mvaConfigName;
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "jetsFromBoth_{"+mvaType+"}^{best};;# events",4,0,4));
    m_histogram[name]->GetXaxis()->SetBinLabel(1, "matched jets");
    m_histogram[name]->GetXaxis()->SetBinLabel(2, "pair from top");
    m_histogram[name]->GetXaxis()->SetBinLabel(3, "pair from higgs");
    m_histogram[name]->GetXaxis()->SetBinLabel(4, "fully solved");
    
    
    const double mvaWeightLow = mvaType=="combined" ? 0. : -1.2;
    const double mvaWeightHigh = mvaType=="combined" ? 1. : 0.2;
    
    name = "mvaWeight"+mvaTypeC+mvaConfigName;
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, "mvaWeight_{"+mvaType+"};w_{MVA}^{"+mvaType+"};# jet pairs", 80, mvaWeightLow, mvaWeightHigh);
    
    name = "mvaWeight"+mvaTypeC+"_best"+mvaTypeC+mvaConfigName;
    this->bookHistosInclExcl(m_histogram, prefix_, step, name, "best mvaWeight_{"+mvaType+"};w_{MVA,1}^{"+mvaType+"};# events", 80, mvaWeightLow, mvaWeightHigh);
    
    name = "mvaWeight"+mvaTypeC+"_TopPairMinusHiggsPair"+mvaConfigName;
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "mvaWeight_{"+mvaType+"}^{top} - mvaWeight_{"+mvaType+"}^{higgs};w_{top}^{"+mvaType+"} - w_{higgs}^{"+mvaType+"};# events", 80, -1.2, 1.2));
    
    name = "mvaWeight"+mvaTypeC+"_TopPairVsHiggsPair"+mvaConfigName;
    m_histogram[name] = this->store(new TH2D(prefix_+name+step, "mvaWeightTopPairVsHiggsPair;w_{MVA}^{top};w_{MVA}^{higgs}", 40, mvaWeightLow, mvaWeightHigh, 40, mvaWeightLow, mvaWeightHigh));
    
    if(mvaType == "combined"){
        name = "mvaWeightCorrectVsSwapped"+mvaConfigName;
        this->bookHistosInclExcl2D(m_histogram, prefix_, step, name, "mvaWeightCorrectVsSwapped;w_{MVA}^{correct};w_{MVA}^{swapped}", 40, -1.2, 0.2, 40, -1.2, 0.2);
        
        name = "mvaWeightCorrectVsSwapped_bestCorrectBestSwapped"+mvaConfigName;
        this->bookHistosInclExcl2D(m_histogram, prefix_, step, name, "mvaWeightCorrectVsSwapped for best Correct pair;w_{MVA,1}^{correct};w_{MVA,1}^{swapped}", 40, -1.2, 0.2, 40, -1.2, 0.2);
    }
    
    
//     name = "dijet_massVsMvaWeightHigh";
//     m_histogram[name] = this->store(new TProfile(prefix_+name+step, "dijetMassVsMvaWeightHigh;w_{MVA,1};m_{jj} [GeV]",20, -1.2, 0.2, "s"));
//     name = "dijet_massVsMvaWeightDiff";
//     m_histogram[name] = this->store(new TProfile(prefix_+name+step, "dijetMassVsMvaWeightDiff;w_{MVA,1} - w_{MVA,2};m_{jj} [GeV]",20, 0, 3, "s"));
//     
//     name = "dijet_massVsMvaWeightHigh_bestCorrect";
//     m_histogram[name] = this->store(new TProfile(prefix_+name+step, "dijetMassVsMvaWeightHigh_{correct}^{best};w_{MVA,1};m_{jj (best)}^{correct} [GeV]",20, -1.2, 0.2, "s"));
//     name = "dijet_massVsMvaWeightHigh_bestSwapped";
//     m_histogram[name] = this->store(new TProfile(prefix_+name+step, "dijetMassVsMvaWeightHigh_{swapped}^{best};w_{MVA,1};m_{jj (best)}^{swapped} [GeV]",20, -1.2, 0.2, "s"));
//     name = "dijet_massVsMvaWeightHigh_bestCombined";
//     m_histogram[name] = this->store(new TProfile(prefix_+name+step, "dijetMassVsMvaWeightHigh_{combined}^{best};w_{MVA,1};m_{jj (best)}^{combined} [GeV]",20, -1.2, 0.2, "s"));
//     
//     name = "dijet_massVsMvaWeightDiff_bestCorrect";
//     m_histogram[name] = this->store(new TProfile(prefix_+name+step, "dijetMassVsMvaWeightDiff_{correct}^{best};w_{MVA,1} - w_{MVA,2};m_{jj (best)}^{correct} [GeV]",20, 0, 3, "s"));
//     name = "dijet_massVsMvaWeightDiff_bestSwapped";
//     m_histogram[name] = this->store(new TProfile(prefix_+name+step, "dijetMassVsMvaWeightDiff_{swapped}^{best};w_{MVA,1} - w_{MVA,2};m_{jj (best)}^{swapped} [GeV]",20, 0, 3, "s"));
//     name = "dijet_massVsMvaWeightDiff_bestCombined";
//     m_histogram[name] = this->store(new TProfile(prefix_+name+step, "dijetMassVsMvaWeightDiff_{combined}^{best};w_{MVA,1} - w_{MVA,2};m_{jj (best)}^{combined} [GeV]",20, 0, 3, "s"));
}



void AnalyzerMvaTopJets::fillHistos(const EventMetadata&,
                                    const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                    const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                                    const KinematicReconstructionSolutions& kinRecoObjects,
                                    const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                                    const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                                    const double& weight, const TString& step,
                                    std::map<TString, TH1*>& m_histogram)
{
    for(const MvaWeightsStruct& mvaWeightsStruct : v_mvaWeightsStruct_){
        // FIXME: treatment of different MVA sets neglected for now, could be treated via folders in root file?
        //const std::string& mvaSetStep = mvaWeightsStruct.stepName();
        
        this->fillHistosPerSet(recoObjects, commonGenObjects, topGenObjects, higgsGenObjects, kinRecoObjects,
                               recoObjectIndices, genObjectIndices, genLevelWeights, recoLevelWeights, weight, step, m_histogram,
                               mvaWeightsStruct);
    }
}



void AnalyzerMvaTopJets::fillHistosPerSet(const RecoObjects& recoObjects, const CommonGenObjects&,
                                          const TopGenObjects&, const HiggsGenObjects&,
                                          const KinematicReconstructionSolutions&,
                                          const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                                          const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                          const double& weight, const TString&,
                                          std::map<TString, TH1*>& m_histogram,
                                          const AnalyzerMvaTopJets::MvaWeightsStruct& mvaWeightsStruct)
{
    TString name;
    
    const tth::IndexPairs& jetIndexPairs = recoObjectIndices.jetIndexPairs_;
    for(const auto& jetIndexPair : jetIndexPairs){
        const LV dijet = recoObjects.jets_->at(jetIndexPair.first) + recoObjects.jets_->at(jetIndexPair.second);
        name = "dijet_mass";
        m_histogram.at(name)->Fill(dijet.M(), weight);
    }
    
    // Loop over all jet combinations and get MVA input variables
    std::vector<MvaVariablesBase*> v_mvaVariables =
            MvaVariablesTopJets::fillVariables(recoObjectIndices, genObjectIndices, recoObjects, weight);
    
    // Get the MVA weights from file, one entry per jet pair
    std::map<std::string, std::vector<float> > m_weightsCorrect;
    for(const auto& weights : mvaWeightsStruct.correctWeights()){
        m_weightsCorrect[weights.first] = weights.second->mvaWeights(v_mvaVariables);
    }
    std::map<std::string, std::vector<float> > m_weightsSwapped;
    for(const auto& weights : mvaWeightsStruct.swappedWeights()){
        m_weightsSwapped[weights.first] = weights.second->mvaWeights(v_mvaVariables);
    }
    std::map<std::string, std::map<std::string, std::vector<float> > > m_weightsCombined;
    for(const auto& weights1 : mvaWeightsStruct.combinedWeights()){
        for(const auto& weights2 : weights1.second){
            m_weightsCombined[weights1.first][weights2.first] = MvaReaderBase::combinedWeights(weights2.second,
                                                                                               m_weightsCorrect.at(weights1.first),
                                                                                               m_weightsSwapped.at(weights2.first));
        }
    }
    
    // Fill all variables and associated MVA weights per event
    const MvaVariablesTopJetsPerEvent mvaVariablesTopJetsPerEvent(v_mvaVariables, m_weightsCorrect, m_weightsSwapped, m_weightsCombined);
    
    // Fill correct trainings
    for(const auto& mvaWeights : mvaVariablesTopJetsPerEvent.mvaWeightsCorrectMap()){
        const std::string& mvaConfigName = mvaWeights.first;
        const std::vector<float>& v_mvaWeight = mvaWeights.second;
        const size_t maxWeightIndex = mvaVariablesTopJetsPerEvent.maxWeightCorrectIndex(mvaConfigName);
        
        this->fillWeightHistos(mvaVariablesTopJetsPerEvent, v_mvaWeight, maxWeightIndex, weight, recoObjectIndices, genObjectIndices, m_histogram,
                               "Correct", mvaConfigName);
        this->fillBestWeightHistos(v_mvaWeight, recoObjects, recoObjectIndices, genObjectIndices, weight, m_histogram,
                                  "Correct", "_"+mvaConfigName);
    }
    
    // Fill swapped trainings
    for(const auto& mvaWeights : mvaVariablesTopJetsPerEvent.mvaWeightsSwappedMap()){
        const std::string& mvaConfigName = mvaWeights.first;
        const std::vector<float>& v_mvaWeight = mvaWeights.second;
        const size_t maxWeightIndex = mvaVariablesTopJetsPerEvent.maxWeightSwappedIndex(mvaConfigName);
        
        this->fillWeightHistos(mvaVariablesTopJetsPerEvent, v_mvaWeight, maxWeightIndex, weight, recoObjectIndices, genObjectIndices, m_histogram,
                               "Swapped", mvaConfigName);
        this->fillBestWeightHistos(v_mvaWeight, recoObjects, recoObjectIndices, genObjectIndices, weight, m_histogram,
                                   "Swapped", "_"+mvaConfigName);
    }
    
    // Fill combined trainings
    for(const auto& mvaWeights1 : mvaVariablesTopJetsPerEvent.mvaWeightsCombinedMap()){
        for(const auto& mvaWeights2 : mvaWeights1.second){
            const std::string& mvaConfigName1 = mvaWeights1.first;
            const std::string& mvaConfigName2 = mvaWeights2.first;
            const std::string mvaConfigName = mvaConfigName1 + "_" + mvaConfigName2;
            const std::vector<float>& v_mvaWeight = mvaWeights2.second;
            const size_t maxWeightIndex = mvaVariablesTopJetsPerEvent.maxWeightCombinedIndex(mvaConfigName1, mvaConfigName2);
            
            this->fillWeightHistos(mvaVariablesTopJetsPerEvent, v_mvaWeight, maxWeightIndex, weight, recoObjectIndices, genObjectIndices, m_histogram,
                                   "Combined", mvaConfigName1, mvaConfigName2);
            this->fillBestWeightHistos(v_mvaWeight, recoObjects, recoObjectIndices, genObjectIndices, weight, m_histogram,
                                       "Combined", "_"+mvaConfigName1+"_"+mvaConfigName2);
        }
    }
    
    // Cleanup
    MvaVariablesTopJets::clearVariables(v_mvaVariables);
}



void AnalyzerMvaTopJets::fillWeightHistos(const MvaVariablesTopJetsPerEvent& mvaVariablesTopJetsPerEvent,
                                          const std::vector<float>& v_mvaWeight, const size_t maxWeightIndex,
                                          const double& weight, const tth::RecoObjectIndices& recoObjectIndices, 
                                          const tth::GenObjectIndices& genObjectIndices, std::map<TString, TH1*>& m_histogram,
                                          const TString& mvaType, const std::string& mvaConfigName1, const std::string& mvaConfigName2)
{
    TString mvaConfigName = "_" + mvaConfigName1;
    if(mvaConfigName2 != "") mvaConfigName.Append("_").Append(mvaConfigName2);
    
    TString name;
    
    // Get the indices of the jet pairs and order them by MVA weights, biggest value first
    const tth::IndexPairs& jetIndexPairs = recoObjectIndices.jetIndexPairs_;
    
    double weightTop = -100.;
    double weightHiggs = -100.;
    for(size_t index = 0; index<mvaVariablesTopJetsPerEvent.variables().size(); ++index){
        const MvaVariablesTopJets* mvaVariablesTopJets = dynamic_cast<const MvaVariablesTopJets*>(mvaVariablesTopJetsPerEvent.variables().at(index));
        if(!mvaVariablesTopJets){
            std::cerr<<"ERROR in AnalyzerMvaTopJets::fillWeightHistos()! MvaVariables are of wrong type, cannot typecast\n...break\n"<<std::endl;
            exit(395);
        }
        
        const bool isMaxWeight = index == maxWeightIndex;
        
        
        name = "mvaWeight"+mvaType+mvaConfigName;
        this->fillHistosInclExcl(m_histogram, name, v_mvaWeight.at(index), mvaVariablesTopJets, weight);
        
        name = "mvaWeight"+mvaType+"_best"+mvaType+mvaConfigName;
        if(isMaxWeight) this->fillHistosInclExcl(m_histogram, name, v_mvaWeight.at(index), mvaVariablesTopJets, weight);
        
        if(mvaType == "Combined"){
            const float weightCorrect = mvaVariablesTopJetsPerEvent.mvaWeightsCorrect(mvaConfigName1).at(index);
            const float weightSwapped = mvaVariablesTopJetsPerEvent.mvaWeightsSwapped(mvaConfigName2).at(index);
            name = "mvaWeightCorrectVsSwapped"+mvaConfigName;
            this->fillHistosInclExcl2D(m_histogram, name, weightCorrect, weightSwapped, mvaVariablesTopJets, weight);
            
            if(index==0){
                const float maxWeightCorrect = mvaVariablesTopJetsPerEvent.maxWeightCorrect(mvaConfigName1);
                const float maxWeightSwapped = mvaVariablesTopJetsPerEvent.maxWeightSwapped(mvaConfigName2);
                name = "mvaWeightCorrectVsSwapped_bestCorrectBestSwapped"+mvaConfigName;
                this->fillHistosInclExcl2D(m_histogram, name, maxWeightCorrect, maxWeightSwapped, mvaVariablesTopJets, weight);
            }
        }
        if(genObjectIndices.isPairFromTop(jetIndexPairs.at(index).first, jetIndexPairs.at(index).second)) weightTop=v_mvaWeight.at(index);
        if(genObjectIndices.isPairFromHiggs(jetIndexPairs.at(index).first, jetIndexPairs.at(index).second)) weightHiggs=v_mvaWeight.at(index);
    }
    
    if(weightTop>-99. && weightHiggs>-99.){
        name = "mvaWeight"+mvaType+"_TopPairVsHiggsPair"+mvaConfigName;
        ((TH2*)m_histogram[name])->Fill(weightTop, weightHiggs, weight);
        
        name = "mvaWeight"+mvaType+"_TopPairMinusHiggsPair"+mvaConfigName;
        m_histogram.at(name)->Fill(weightTop - weightHiggs, weight);
    }
}



void AnalyzerMvaTopJets::fillBestWeightHistos(const std::vector<float>& v_mvaWeights,
                                              const RecoObjects& recoObjects,
                                              const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                                              const double& weight,
                                              std::map<TString, TH1*>& m_histogram,
                                              const std::string& mvaType, const std::string& mvaConfigName)
{
    TString name;
    
    // Get the indices of the jet pairs and order them by MVA weights, biggest value first
    const tth::IndexPairs& jetIndexPairs = recoObjectIndices.jetIndexPairs_;
    std::vector<int> jetIndexPairsIndices = common::initialiseIndices(jetIndexPairs);
    common::orderIndices(jetIndexPairsIndices, v_mvaWeights);
    
    // Get jet pair leading in MVA weight, and extract bIndex and antiBIndex
    const std::pair<int, int>& leadingJetIndexPair = jetIndexPairs.at(jetIndexPairsIndices.at(0));
    const int antiBFromTopIndex = leadingJetIndexPair.first;
    const int bFromTopIndex = leadingJetIndexPair.second;
    
    // Check whether the two jets correspond to the b's from tops, and if the two are correct or swapped
    const bool isSwappedPairFromTop(genObjectIndices.isSwappedPairFromTop(bFromTopIndex, antiBFromTopIndex));
    const bool isCorrectPairFromTop(genObjectIndices.isCorrectPairFromTop(bFromTopIndex, antiBFromTopIndex));
    const bool isPairFromTop(genObjectIndices.isPairFromTop(bFromTopIndex, antiBFromTopIndex));
    
    
    // Histograms for the top system
    name = "jetsFromTop_best"+mvaType+mvaConfigName;
    if(genObjectIndices.uniqueRecoTopMatching()){
        m_histogram.at(name)->Fill(0., weight);
        if(isCorrectPairFromTop){
            m_histogram.at(name)->Fill(1., weight);
            m_histogram.at(name)->Fill(2., weight);
        }
        else if(isSwappedPairFromTop){
            m_histogram.at(name)->Fill(1., weight);
            m_histogram.at(name)->Fill(3., weight);
        }
    }
    else{
        m_histogram.at(name)->Fill(-999., weight);
    }
    
    
    // Get all jets except the ones assigned to ttbar system, and order them by b-tag discriminator value
    std::vector<int> remainingJetIndices;
    for(const int index : recoObjectIndices.jetIndices_){
        if(index==antiBFromTopIndex || index==bFromTopIndex) continue;
        remainingJetIndices.push_back(index);
    }
    common::orderIndices(remainingJetIndices, *recoObjects.jetBTagCSV_);
    
    
    // Calculations requiring the presence of at least 4 jets
    if(recoObjectIndices.jetIndices_.size()>=4){
        // Get the two jets assigned to Higgs
        const int jet1FromHiggsIndex = remainingJetIndices.at(0);
        const int jet2FromHiggsIndex = remainingJetIndices.at(1);
        
        // Check whether the two jets correspond to the b's from Higgs
        const bool isPairFromHiggs(genObjectIndices.isPairFromHiggs(jet1FromHiggsIndex, jet2FromHiggsIndex));
        
        const LV dijet = recoObjects.jets_->at(jet1FromHiggsIndex) + recoObjects.jets_->at(jet2FromHiggsIndex);
        name = "dijet_mass_best"+mvaType+mvaConfigName;
        m_histogram.at(name)->Fill(dijet.M(), weight);
        
        name = "jetsFromHiggs_best"+mvaType+mvaConfigName;
        if(genObjectIndices.uniqueRecoHiggsMatching()){
            m_histogram.at(name)->Fill(0., weight);
            if(isPairFromHiggs){
                m_histogram.at(name)->Fill(1., weight);
                if(!genObjectIndices.uniqueRecoMatching()) m_histogram.at(name)->Fill(2., weight);
                else if(isPairFromTop) m_histogram.at(name)->Fill(3., weight);
            }
        }
        else{
            m_histogram.at(name)->Fill(-999., weight);
        }
        
        name = "jetsFromBoth_best"+mvaType+mvaConfigName;
        if(genObjectIndices.uniqueRecoMatching()){
            m_histogram.at(name)->Fill(0., weight);
            if(isPairFromTop) m_histogram.at(name)->Fill(1., weight);
            if(isPairFromHiggs){
                m_histogram.at(name)->Fill(2., weight);
                if(isPairFromTop) m_histogram.at(name)->Fill(3., weight);
            }
        }
    }
}



void AnalyzerMvaTopJets::bookHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& prefix, const TString& step,
                                            const TString& name, const TString& title,
                                            const int& nBinX, const double& xMin, const double& xMax)
{
    const TString correct("_correct");
    const TString swapped("_swapped");
    const TString wrong("_wrong");
    
    m_histogram[name] = this->store(new TH1D(prefix+name+step, title, nBinX, xMin, xMax));
    
    m_histogram[name+correct] = this->store(new TH1D(prefix+name+correct+step, title, nBinX, xMin, xMax));
    m_histogram[name+swapped] = this->store(new TH1D(prefix+name+swapped+step, title, nBinX, xMin, xMax));
    m_histogram[name+wrong] = this->store(new TH1D(prefix+name+wrong+step, title, nBinX, xMin, xMax));
}



void AnalyzerMvaTopJets::bookHistosInclExcl2D(std::map<TString, TH1*>& m_histogram, const TString& prefix, const TString& step,
                                              const TString& name, const TString& title,
                                              const int& nBinX, const double& xMin, const double& xMax,
                                              const int& nBinY, const double& yMin, const double& yMax)
{
    const TString correct("_correct");
    const TString swapped("_swapped");
    const TString wrong("_wrong");
    
    m_histogram[name] = this->store(new TH2D(prefix+name+step, title, nBinX, xMin, xMax, nBinY, yMin, yMax));
    
    m_histogram[name+correct] = this->store(new TH2D(prefix+name+correct+step, title, nBinX, xMin, xMax, nBinY, yMin, yMax));
    m_histogram[name+swapped] = this->store(new TH2D(prefix+name+swapped+step, title, nBinX, xMin, xMax, nBinY, yMin, yMax));
    m_histogram[name+wrong] = this->store(new TH2D(prefix+name+wrong+step, title, nBinX, xMin, xMax, nBinY, yMin, yMax));
}



void AnalyzerMvaTopJets::fillHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& name,
                                            const float& variable,
                                            const MvaVariablesTopJets* mvaTopJetsVariables, const double& weight)
{
    const TString correct("_correct");
    const TString swapped("_swapped");
    const TString wrong("_wrong");
    
    m_histogram.at(name)->Fill(variable, weight);
    
    if(mvaTopJetsVariables->correctCombination_.value_) m_histogram.at(name+correct)->Fill(variable, weight);
    else if(mvaTopJetsVariables->swappedCombination_.value_) m_histogram.at(name+swapped)->Fill(variable, weight);
    else m_histogram.at(name+wrong)->Fill(variable, weight);
}



void AnalyzerMvaTopJets::fillHistosInclExcl2D(std::map<TString, TH1*>& m_histogram, const TString& name,
                                              const float& variable1, const float& variable2,
                                              const MvaVariablesTopJets* mvaTopJetsVariables, const double& weight)
{
    const TString correct("_correct");
    const TString swapped("_swapped");
    const TString wrong("_wrong");
    
    ((TH2*)m_histogram.at(name))->Fill(variable1, variable2, weight);
    
    if(mvaTopJetsVariables->correctCombination_.value_) ((TH2*)m_histogram.at(name+correct))->Fill(variable1, variable2, weight);
    else if(mvaTopJetsVariables->swappedCombination_.value_) ((TH2*)m_histogram.at(name+swapped))->Fill(variable1, variable2, weight);
    else ((TH2*)m_histogram.at(name+wrong))->Fill(variable1, variable2, weight);
}










AnalyzerMvaTopJets::MvaWeightsStruct::MvaWeightsStruct(const std::string& stepName,
                                                       const std::vector<std::string>& v_nameCorrect,
                                                       const std::vector<std::string>& v_nameSwapped,
                                                       const char* mva2dWeightsFile):
stepName_(stepName)
{
    // Extract the path to the folder containing the root file and xml files with weights
    TString mva2dWeightsFolder(mva2dWeightsFile);
    mva2dWeightsFolder.Remove(mva2dWeightsFolder.Last('/')+1);
    
    // Access correct weights
    for(const auto& nameCorrect : v_nameCorrect){
        TString weightsCorrectFilename(mva2dWeightsFolder);
        weightsCorrectFilename.Append("correct_").Append(stepName).Append("_").Append(nameCorrect).Append(".weights.xml");
        
        m_correct_[nameCorrect] = new MvaReaderTopJets("BDT method");
        m_correct_.at(nameCorrect)->book(weightsCorrectFilename);
        
        // Access combined weights
        for(const auto& nameSwapped : v_nameSwapped){
            TString weights2dHistoname("combined_");
            weights2dHistoname.Append(stepName).Append("_").Append(nameCorrect).Append("_").Append(nameSwapped);
            
            RootFileReader* fileReader(RootFileReader::getInstance());
            // FIXME: do normalisation already during storage
            TH2D* weightsCombined = fileReader->GetClone<TH2D>(mva2dWeightsFile, weights2dHistoname);
            const double integral = weightsCombined->Integral(0, weightsCombined->GetNbinsX()+1, 0, weightsCombined->GetNbinsY()+1);
            weightsCombined->Scale(1./integral);
            m_combined_[nameCorrect][nameSwapped] = weightsCombined;
        }
    }
    
    // Access swapped weights
    for(const auto& nameSwapped : v_nameSwapped){
        TString weightsSwappedFilename(mva2dWeightsFolder);
        weightsSwappedFilename.Append("swapped_").Append(stepName).Append("_").Append(nameSwapped).Append(".weights.xml");
        
        m_swapped_[nameSwapped] = new MvaReaderTopJets("BDT method");
        m_swapped_.at(nameSwapped)->book(weightsSwappedFilename);
    }
}








