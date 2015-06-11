#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCut.h>
#include <TMVA/Tools.h>
#include <TMVA/Config.h>
#include <TMVA/Factory.h>
#include <TMVA/MethodBase.h>

#include "MvaFactoryBase.h"
#include "MvaVariablesBase.h"
#include "mvaSetup.h"
#include "higgsUtils.h"
#include "Samples.h"
#include "../../common/include/RootFileReader.h"





MvaFactoryBase::MvaFactoryBase(const TString& mvaOutputDir, const TString& weightFileDir,
                               const Samples& samples,
                               const bool inOneFactory):
mvaOutputDir_(mvaOutputDir),
weightFileDir_(weightFileDir),
samples_(samples),
inOneFactory_(inOneFactory),
fileReader_(RootFileReader::getInstance())
{
    
}



MvaFactoryBase::MvaFactoryBase(const TString& mvaOutputDir, const TString& weightFileDir,
                               const TString& treeFileName,
                               const bool inOneFactory):
mvaOutputDir_(mvaOutputDir),
weightFileDir_(weightFileDir),
samples_(Samples()),
treeFileName_(treeFileName),
inOneFactory_(inOneFactory),
fileReader_(0)
{
    // Check if input file exists
    std::ifstream inputFileStream(treeFileName);
    if(!inputFileStream.is_open()){
        std::cerr<<"ERROR in constructor of MvaFactoryBase! File containing the MVA input trees not found: "<<treeFileName
                 <<"\n...break\n"<<std::endl;
        exit(1);
    }
    inputFileStream.close();
    
    if(inOneFactory_) std::cout<<"Running all training configurations in one single factory\n";
    else std::cout<<"Running each training configuration in own factory\n";
}



void MvaFactoryBase::train(const std::vector<mvaSetup::MvaSet>& v_mvaSet)const
{
    std::cout<<"--- Beginning MVA training\n";
    
    // Open input file
    TFile* treeFile = TFile::Open(treeFileName_);
    
    // Loop over the steps/categories and train MVAs
    for(const auto& mvaSet : v_mvaSet){
        const TString mergedStepName = tth::stepName(mvaSet.step_, mvaSet.v_category_);
        this->trainMva(treeFile, mvaSet, mergedStepName);
    }
    
    // Cleanup
    treeFile->Close();
    
    std::cout<<"=== Finishing MVA training\n\n";
}



void MvaFactoryBase::train2(const std::vector<mvaSetup::MvaSet>& v_mvaSet)const
{
    std::cout<<"--- Beginning MVA training\n";
    
    // Fill map containing scale factors for all steps
    std::map<TString, SystematicChannelFactors> m_stepFactors;
    for(const auto& mvaSet : v_mvaSet){
        const TString step = tth::stepName(mvaSet.step_, mvaSet.v_category_);
        m_stepFactors[step] = samples_.globalWeights(step).first;
    }
    
    // Loop over all channels and systematics
    const SystematicChannelSamples& m_systematicChannelSample(samples_.getSystematicChannelSamples());
    for(const auto& systematicChannelSamples : m_systematicChannelSample){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        for(const auto& channelSample : systematicChannelSamples.second){
            const Channel::Channel& channel(channelSample.first);
            const std::vector<Sample>& v_sample(channelSample.second);
            const TString outputFolder = common::assignFolder(mvaOutputDir_, channel, systematic);
            
            // Create text file for writing training names
            TString trainingNameFilename(outputFolder);
            trainingNameFilename.Append("trainingNames.txt");
            std::ofstream trainingNameFile;
            trainingNameFile.open(trainingNameFilename);
            
            // Loop over the steps/categories and train MVAs
            for(const auto& mvaSet : v_mvaSet){
                // Check if training should be applied for channel
                if(std::find(mvaSet.v_channel_.begin(), mvaSet.v_channel_.end(), channel) == mvaSet.v_channel_.end()) continue;
                
                const TString step = tth::stepName(mvaSet.step_, mvaSet.v_category_);
                const SystematicChannelFactors globalWeights = m_stepFactors.at(step);
                const std::vector<double>& v_weight(globalWeights.at(systematic).at(channel));
                this->trainMva2(mvaSet, outputFolder, v_sample, v_weight, trainingNameFile, step);
            }
            
            // Cleanup
            trainingNameFile.close();
        }
    }
    
    std::cout<<"=== Finishing MVA training\n\n";
}



void MvaFactoryBase::trainMva(TFile* const, const mvaSetup::MvaSet&, const TString&)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method trainMva() in MvaFactoryBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void MvaFactoryBase::trainMva2(const mvaSetup::MvaSet&, const TString&,
                               const std::vector<Sample>&, const std::vector<double>&,
                               std::ofstream&, const TString&)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method trainMva() in MvaFactoryBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void MvaFactoryBase::runMva(const char* const methodPrefix, const TCut& cutSignal, const TCut& cutBackground,
                            TTree* const treeTraining, TTree* const treeTesting,
                            const std::vector<mvaSetup::MvaConfig>& v_mvaConfig,
                            const TString& stepName)const
{
    // Output directory for the weights
    TString mvaOutputWeightsDirectory(mvaOutputDir_);
    mvaOutputWeightsDirectory.Append(weightFileDir_);
    
    // Output filename bases for root files and weight files
    const TString filenameBase = ((TString)methodPrefix).Append(stepName);
    const TString mvaOutputFilenameBase = ((TString)mvaOutputDir_).Append(filenameBase);
    
    // Train all configs in one factory, or each in individual factory
    if(inOneFactory_){
        // Get a TMVA instance
        TMVA::Tools::Instance();
        
        // Set the output directory for the weights (if not specified, default is "weights")
        (TMVA::gConfig().GetIONames()).fWeightFileDir = mvaOutputWeightsDirectory;
        
        // Create ROOT output file for TMVA
        TString mvaOutputFilename(mvaOutputFilenameBase);
        mvaOutputFilename.Append(".root");
        TFile* outputFile = TFile::Open(mvaOutputFilename, "RECREATE");
        
        // Create and configure factory
        TMVA::Factory* factory = new TMVA::Factory(filenameBase, outputFile, "!V:!Silent");
        this->configureFactory(factory, cutSignal, cutBackground, treeTraining, treeTesting);
        
        // Book the MVA methods (e.g. boosted decision tree with specific setup)
        for(const mvaSetup::MvaConfig& mvaConfig : v_mvaConfig)
            factory->BookMethod(mvaConfig.mvaType_, mvaConfig.methodAppendix_, mvaConfig.options_);
        
        // Run factory
        factory->TrainAllMethods();
        factory->TestAllMethods();
        factory->EvaluateAllMethods();
        
        // Cleanup
        outputFile->Close();
        delete factory;
    }
    else{
        for(const mvaSetup::MvaConfig& mvaConfig : v_mvaConfig){
            TMVA::Tools::Instance();
            
            // Set the output directory for the weights (if not specified, default is "weights")
            (TMVA::gConfig().GetIONames()).fWeightFileDir = mvaOutputWeightsDirectory;
            
            // Create ROOT output file for TMVA
            TDirectory* oldDir = gDirectory;
            TString mvaOutputFilename(mvaOutputFilenameBase);
            mvaOutputFilename.Append("_").Append(mvaConfig.methodAppendix_).Append(".root");
            TFile* outputFile = TFile::Open(mvaOutputFilename, "RECREATE");
            oldDir->cd();
            
            // Create and configure factory
            TMVA::Factory* factory = new TMVA::Factory(filenameBase, outputFile, "!V:!Silent");
            this->configureFactory(factory, cutSignal, cutBackground, treeTraining, treeTesting);
            
            // Book the MVA method (e.g. boosted decision tree with specific setup)
            //TMVA::MethodBase* method = factory->BookMethod(mvaConfig.mvaType_, mvaConfig.methodAppendix_, mvaConfig.options_);
            factory->BookMethod(mvaConfig.mvaType_, mvaConfig.methodAppendix_, mvaConfig.options_);
            
            // Run factory
            factory->TrainAllMethods();
            factory->TestAllMethods();
            factory->EvaluateAllMethods();
            
            // Trials of cleanup to avoid memory leak in TMVA (not successful)
            //for(int i = 0; i < 1000000; ++i) if(!(i%5)) std::cout<<"Something to show for testing: "<<mvaConfig.methodAppendix_<<"\n";
            //method->Clear();
            //method->Delete();
            //factory->DeleteAllMethods();
            //factory->Clear();
            //factory->Delete();
            //TMVA::Tools::DestroyInstance();
            
            // Cleanup
            outputFile->Close();
            delete factory;
        }
    }
}



void MvaFactoryBase::runMva2(const std::vector<mvaSetup::MvaConfig>& v_mvaConfig, const TString& outputFolder,
                            const std::vector<Sample>& v_sample, const std::vector<double>& v_weight,
                            std::ofstream& trainingNameFile, const TString& stepName)const
{
    const TCut cutSignal = "";
    const TCut cutBackground = "";
    
    // Output directory for the weights
    TString mvaOutputWeightsDirectory(outputFolder);
    mvaOutputWeightsDirectory.Append(weightFileDir_);
    
    // Output filename bases for root files and weight files
    const TString filenameBase = TString("classification").Append(stepName);
    const TString mvaOutputFilenameBase = TString(outputFolder).Append(filenameBase);
    
    // Train all configs in one factory, or each in individual factory
    if(inOneFactory_){
        // Get a TMVA instance
        TMVA::Tools::Instance();
        
        // Set the output directory for the weights (if not specified, default is "weights")
        (TMVA::gConfig().GetIONames()).fWeightFileDir = mvaOutputWeightsDirectory;
        
        // Create ROOT output file for TMVA
        TString mvaOutputFilename(mvaOutputFilenameBase);
        mvaOutputFilename.Append(".root");
        TFile* outputFile = TFile::Open(mvaOutputFilename, "RECREATE");
        
        // Create and configure factory
        TMVA::Factory* factory = new TMVA::Factory(filenameBase, outputFile, "!V:!Silent");
        this->configureFactory2(factory, cutSignal, cutBackground, v_sample, v_weight, stepName);
        
        // Book the MVA methods (e.g. boosted decision tree with specific setup)
        for(const mvaSetup::MvaConfig& mvaConfig : v_mvaConfig){
            factory->BookMethod(mvaConfig.mvaType_, mvaConfig.methodAppendix_, mvaConfig.options_);
            
            // Write training names to text file
            TString trainingName(mvaOutputWeightsDirectory);
            trainingName.Append(filenameBase).Append("_").Append(mvaConfig.methodAppendix_).Append(".weights.xml");
            trainingNameFile<<trainingName<<std::endl;
        }
        
        // Run factory
        factory->TrainAllMethods();
        factory->TestAllMethods();
        factory->EvaluateAllMethods();
        
        // Cleanup
        outputFile->Close();
        delete factory;
    }
    else{
        for(const mvaSetup::MvaConfig& mvaConfig : v_mvaConfig){
            TMVA::Tools::Instance();
            
            // Set the output directory for the weights (if not specified, default is "weights")
            (TMVA::gConfig().GetIONames()).fWeightFileDir = mvaOutputWeightsDirectory;
            
            // Create ROOT output file for TMVA
            TDirectory* oldDir = gDirectory;
            TString mvaOutputFilename(mvaOutputFilenameBase);
            mvaOutputFilename.Append("_").Append(mvaConfig.methodAppendix_).Append(".root");
            TFile* outputFile = TFile::Open(mvaOutputFilename, "RECREATE");
            oldDir->cd();
            
            // Create and configure factory
            TMVA::Factory* factory = new TMVA::Factory(filenameBase, outputFile, "!V:!Silent");
            this->configureFactory2(factory, cutSignal, cutBackground, v_sample, v_weight, stepName);
            
            // Book the MVA method (e.g. boosted decision tree with specific setup)
            //TMVA::MethodBase* method = factory->BookMethod(mvaConfig.mvaType_, mvaConfig.methodAppendix_, mvaConfig.options_);
            factory->BookMethod(mvaConfig.mvaType_, mvaConfig.methodAppendix_, mvaConfig.options_);
            
            // Write training names to text file
            TString trainingName(mvaOutputWeightsDirectory);
            trainingName.Append(filenameBase).Append("_").Append(mvaConfig.methodAppendix_).Append(".weights.xml");
            trainingNameFile<<trainingName<<std::endl;
            
            // Run factory
            factory->TrainAllMethods();
            factory->TestAllMethods();
            factory->EvaluateAllMethods();
            
            // Trials of cleanup to avoid memory leak in TMVA (not successful)
            //for(int i = 0; i < 1000000; ++i) if(!(i%5)) std::cout<<"Something to show for testing: "<<mvaConfig.methodAppendix_<<"\n";
            //method->Clear();
            //method->Delete();
            //factory->DeleteAllMethods();
            //factory->Clear();
            //factory->Delete();
            //TMVA::Tools::DestroyInstance();
            
            // Cleanup
            outputFile->Close();
            delete factory;
        }
    }
}



void MvaFactoryBase::configureFactory(TMVA::Factory* const,
                                      const TCut&, const TCut&,
                                      TTree* const, TTree* const)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method configureFactory() in MvaFactoryBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void MvaFactoryBase::configureFactory2(TMVA::Factory* const,
                                      const TCut&, const TCut&,
                                      const std::vector<Sample>&, const std::vector<double>&,
                                      const TString&)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method configureFactory() in MvaFactoryBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



TTree* MvaFactoryBase::readTree(const TString&, const TString&)
{
    //TTree* const tree = fileReader_->Get<TTree*>(filename, treename);
    return 0;
}



void MvaFactoryBase::addVariable(TMVA::Factory* const factory, const MvaVariableInt& variable)const
{
    factory->AddVariable(variable.name().data(), *variable.type());
}



void MvaFactoryBase::addVariable(TMVA::Factory* const factory, const MvaVariableFloat& variable)const
{
    factory->AddVariable(variable.name().data(), *variable.type());
}



void MvaFactoryBase::addSpectator(TMVA::Factory* const factory, const MvaVariableInt& variable)const
{
    factory->AddSpectator(variable.name().data(), *variable.type());
}



void MvaFactoryBase::addSpectator(TMVA::Factory* const factory, const MvaVariableFloat& variable)const
{
    factory->AddSpectator(variable.name().data(), *variable.type());
}










