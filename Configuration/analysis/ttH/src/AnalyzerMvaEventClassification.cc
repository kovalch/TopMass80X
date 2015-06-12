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

#include "AnalyzerMvaEventClassification.h"
#include "MvaReaderBase.h"
#include "MvaReaderEventClassification.h"
#include "MvaVariablesEventClassification.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "higgsUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
#include "../../common/include/RootFileReader.h"
#include "../../common/include/KinematicReconstructionSolution.h"
#include "../../common/include/sampleHelpers.h"





AnalyzerMvaEventClassification::AnalyzerMvaEventClassification(const std::vector<TString>& selectionStepsNoCategories,
                                                               const std::vector<TString>& stepsForCategories,
                                                               const JetCategories* jetCategories):
AnalyzerBase("mvaEventA_", selectionStepsNoCategories, stepsForCategories, jetCategories)
{
    std::cout<<"--- Beginning setting up MVA validation\n";
    
    // Sanity check
    if(!stepsForCategories.size() || !jetCategories){
        std::cerr<<"ERROR in Constructor of AnalyzerMvaEventClassification! No steps for categories defined, but required\n...break\n"<<std::endl;
        exit(912);
    }
    
    // Access training names from text file
    TString filename = common::accessFolder("mvaOutput_mvaEvent", Channel::combined, Systematic::nominalSystematic());
    filename.Append("trainingNames.txt");
    const std::vector<TString> v_trainingName = common::readFile(filename);
    
    // Loop over steps, and if training is valid for this step, set up reader
    for(const auto& stepShort : stepsForCategories){
        for(int category = 0; category<jetCategories->numberOfCategories(); ++category){
            const TString step = tth::stepName(stepShort, category);
            for(const TString& training : v_trainingName){
                if(!training.Contains(step)) continue;
                std::cout<<"\n\ntraining: "<<step<<" , "<<training<<"\n\n";
                m_m_mvaWeight_[step][training] = new MvaReaderEventClassification("BDT method");
                m_m_mvaWeight_.at(step).at(training)->book(training);
            }
        }
    }

    std::cout<<"=== Finishing setting up MVA validation\n\n";
}



void AnalyzerMvaEventClassification::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    // Exclude steps where no valid training exists
    if(m_m_mvaWeight_.find(step) == m_m_mvaWeight_.end()) return;
    
    TString name;
    
    // Loop over trainings for given step and book histograms
    const auto& m_mvaWeight = m_m_mvaWeight_.at(step);
    for(const auto& mvaWeight : m_mvaWeight){
        const TString& trainingName = mvaWeight.first;
        
        name = "mvaWeight_"+trainingName;
        m_histogram[name] = this->store(new TH1D(prefix_+name+step, "mvaWeight;w_{MVA};# events", 80, -1.2, 0.2));
    }
}



void AnalyzerMvaEventClassification::fillHistos(const EventMetadata&,
                                                const RecoObjects& recoObjects, const CommonGenObjects&,
                                                const TopGenObjects&, const HiggsGenObjects&,
                                                const KinematicReconstructionSolutions&,
                                                const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices&,
                                                const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                                const double& weight, const TString& step,
                                                std::map<TString, TH1*>& m_histogram)
{
    // Exclude steps where no valid training exists
    if(m_m_mvaWeight_.find(step) == m_m_mvaWeight_.end()) return;
    
    TString name;
    
    // Get MVA input variables
    MvaVariablesBase* mvaVariables = MvaVariablesEventClassification::fillVariables(recoObjectIndices, recoObjects, weight);
    
    // Loop over trainings for given step and fill histograms
    const auto& m_mvaWeight = m_m_mvaWeight_.at(step);
    for(const auto& mvaWeight : m_mvaWeight){
        const TString& trainingName = mvaWeight.first;
        const float trainingWeight = mvaWeight.second->mvaWeight(mvaVariables);
        
        name = "mvaWeight_"+trainingName;
        m_histogram.at(name)->Fill(trainingWeight, weight);
    }
}



/*
AnalyzerMvaEventClassification::MvaWeightsStruct::MvaWeightsStruct(const std::string& stepName,
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
        
    }
}
*/







