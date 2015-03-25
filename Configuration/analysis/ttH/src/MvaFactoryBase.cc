#include <vector>
#include <iostream>
#include <fstream>

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





MvaFactoryBase::MvaFactoryBase(const TString& mvaOutputDir, const TString& weightFileDir,
                               const TString& treeFileName,
                               const bool inOneFactory):
mvaOutputDir_(mvaOutputDir),
weightFileDir_(weightFileDir),
treeFileName_(treeFileName),
inOneFactory_(inOneFactory)
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



void MvaFactoryBase::trainMva(TFile* const, const mvaSetup::MvaSet&, const TString&)const
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



void MvaFactoryBase::configureFactory(TMVA::Factory* const,
                                      const TCut&, const TCut&,
                                      TTree* const, TTree* const)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method configureFactory() in MvaFactoryBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
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










