#include <iostream>
#include <cstdlib>
#include <fstream>
#include <utility>

#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TSelectorList.h>
#include <TIterator.h>
#include <TObject.h>
#include <TList.h>
#include <Rtypes.h>

#include "MvaTreeHandlerBase.h"
#include "MvaVariablesBase.h"
#include "JetCategories.h"
#include "analysisStructs.h"
#include "higgsUtils.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/analysisObjectStructs.h"





MvaTreeHandlerBase::MvaTreeHandlerBase(const TString& prefix,
                                       const char* mvaInputDir,
                                       const std::vector<TString>& selectionStepsNoCategories,
                                       const std::vector<TString>& stepsForCategories,
                                       const JetCategories* jetCategories):
prefix_(prefix),
selectorList_(0),
selectionSteps_(selectionStepsNoCategories),
stepsForCategories_(stepsForCategories),
jetCategories_(jetCategories),
mvaInputDir_(mvaInputDir)
{
    if(!jetCategories_ && stepsForCategories_.size()>0){
        std::cerr<<"ERROR in constructor for MvaTreeHandlerBase! "
                 <<"No jet categories passed, but request for category-wise selection steps\n...break\n"<<std::endl;
        exit(236);
    }
    if(jetCategories_ && stepsForCategories_.size()==0){
        std::cerr<<"ERROR in constructor for MvaTreeHandlerBase! "
                 <<"Jet categories passed, but no category-wise selection steps defined\n...break\n"<<std::endl;
        exit(237);
    }
}



void MvaTreeHandlerBase::book()
{
    for(const auto& stepShort : selectionSteps_){
        const TString step = tth::stepName(stepShort);
        this->addStep(step);
    }
    
    for(const auto& stepShort : stepsForCategories_){
        for(int category = 0; category<jetCategories_->numberOfCategories(); ++category){
            const TString step = tth::stepName(stepShort, category);
            this->addStep(step);
        }
    }
}



void MvaTreeHandlerBase::addStep(const TString& step)
{
    // Check whether step already exists
    if(this->checkExistence(step)){
        std::cout<<"Warning in addStep()! Selection step already contained: "<<step
                 <<"\n...skip this one\n";
        return;
    }

    // Book the step
    m_stepMvaVariables_[step] = std::vector<MvaVariablesBase*>();
}



bool MvaTreeHandlerBase::checkExistence(const TString& step)const
{
    return m_stepMvaVariables_.find(step) != m_stepMvaVariables_.end();
}



void MvaTreeHandlerBase::fill(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                              const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                              const KinRecoObjects& kinRecoObjects,
                              const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                              const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                              const double& weight, const TString& stepShort)
{
    // Number of selected jets and bjets
    const int numberOfJets = recoObjectIndices.jetIndices_.size();
    const int numberOfBjets = recoObjectIndices.bjetIndices_.size();
    
    // Set up step name and check if step exists
    const bool stepInCategory = stepShort.Contains("_cate");
    const TString step = stepInCategory ? stepShort : tth::stepName(stepShort);
    const bool stepExists(this->checkExistence(step));
    if(!stepInCategory && jetCategories_){
        const int categoryId = jetCategories_->categoryId(numberOfJets, numberOfBjets);
        
        // Here check the individual jet categories
        const TString fullStepName = tth::stepName(stepShort, categoryId);
        this->fill(recoObjects, commonGenObjects,
                   topGenObjects, higgsGenObjects,
                   kinRecoObjects,
                   recoObjectIndices, genObjectIndices,
                   genLevelWeights, recoLevelWeights,
                   weight, fullStepName);
    }
    if(!stepExists) return;
    
    // Fill the MVA variables of the specific treeHandler
    std::vector<MvaVariablesBase*>& mvaVariables = m_stepMvaVariables_.at(step);
    this->fillVariables(recoObjects, commonGenObjects,
                        topGenObjects, higgsGenObjects,
                        kinRecoObjects,
                        recoObjectIndices, genObjectIndices,
                        genLevelWeights, recoLevelWeights,
                        weight, step,
                        mvaVariables);
}



void MvaTreeHandlerBase::writeTrees(const std::string& outputFilename,
                                    const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    // Create output file for MVA tree
    std::string f_savename = static_cast<std::string>(common::assignFolder(mvaInputDir_, channel, systematic));
    f_savename.append(outputFilename);
    TFile outputFile(f_savename.c_str(),"RECREATE");
    std::cout<<"\nOutput file for MVA input trees: "<<f_savename<<"\n";
    
    // Produce MVA input TTree and store it in output
    TSelectorList* output = new TSelectorList();
    this->writeTrees(output);
    
    // Write file and cleanup
    TIterator* it = output->MakeIterator();
    while(TObject* obj = it->Next()){
        obj->Write();
    }
    outputFile.Close();
    output->SetOwner();
    output->Clear();
}



void MvaTreeHandlerBase::writeTrees(TSelectorList* output)
{
    std::cout<<"--- Beginning production of MVA input trees\n";
    
    // Set pointer to output, so that TTrees are owned by it
    selectorList_ = output;
    
    std::map<TString, TTree*> m_stepTree;
    for(const auto& stepMvaVariables : m_stepMvaVariables_){
        const TString& step = stepMvaVariables.first;
        const std::vector<MvaVariablesBase*>& v_mvaVariables = stepMvaVariables.second;
        TTree* tree = m_stepTree[step];
        tree = this->store(new TTree(prefix_+"mvaVariables"+step, prefix_+"mvaVariables"));
        this->createAndFillBranches(tree, v_mvaVariables);
    }
    
    std::cout<<"MVA variables multiplicity per step (step, no. of entries):\n";
    for(auto vars : m_stepMvaVariables_){
        std::cout<<"\t"<<vars.first<<" , "<<vars.second.size()<<"\n";
    }
    
    std::cout<<"=== Finishing production of MVA input trees\n\n";
}



void MvaTreeHandlerBase::createAndFillBranches(TTree*, const std::vector<MvaVariablesBase*>&)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method createAndFillBranches() in MvaTreeHandlerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void MvaTreeHandlerBase::createBranch(TTree* tree, const MvaVariableInt& variable)const
{
    std::string name(variable.name());
    std::string nameType(name);
    nameType.append("/").append(variable.type());
    tree->Branch(name.data(), (Long_t)&variable.value_, nameType.data());
}



void MvaTreeHandlerBase::createBranch(TTree* tree, const MvaVariableFloat& variable)const
{
    std::string name(variable.name());
    std::string nameType(name);
    nameType.append("/").append(variable.type());
    tree->Branch(name.c_str(), (Long_t)&variable.value_, nameType.data());
}



void MvaTreeHandlerBase::importTrees(const std::string& f_savename, const std::string& prefix)
{
    std::cout<<"--- Beginning import of TTrees with MVA variables\n";
    
    // Open input file
    TFile* inputFile = TFile::Open(f_savename.c_str());
    if(inputFile->IsZombie()){
        std::cerr<<"ERROR in importTrees()! Cannot open input file to import TTrees, filename is: "
                 <<f_savename<<"\n...break\n"<<std::endl;
        exit(77);
    }
    
    // Find all trees of all steps/categories containing MVA input variables
    const TString treeName = prefix;
    const std::vector<std::pair<TString, TString> > v_nameStepPair =
        tth::nameStepPairs(f_savename.c_str(), treeName);
    
    // Loop over steps and import trees
    this->clear();
    for(const auto& nameStepPair : v_nameStepPair){
        TTree* tree(0);
        tree = dynamic_cast<TTree*>(inputFile->Get(nameStepPair.first));
        if(!tree){
            std::cerr<<"ERROR in importTrees()! TTree not found in file, tree name is: "
                     <<treeName<<"\n...break\n"<<std::endl;
            exit(78);
        }
        this->importBranches(tree, m_stepMvaVariables_[nameStepPair.second]);
        tree = 0;
        
        std::cout<<"Step, Number of entries: "<<nameStepPair.second<<" , "<<m_stepMvaVariables_[nameStepPair.second].size()<<"\n";
    }
    inputFile->Close();
    
    std::cout<<"=== Finishing import of TTrees with MVA variables\n\n";
}



void MvaTreeHandlerBase::importBranches(TTree*, std::vector<MvaVariablesBase*>&)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method importBranches() in MvaTreeHandlerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void MvaTreeHandlerBase::importBranch(TTree* tree, MvaVariableInt& variable)const
{
    tree->SetBranchAddress(variable.name().data(), &variable.value_, &variable.branch_);
}



void MvaTreeHandlerBase::importBranch(TTree* tree, MvaVariableFloat& variable)const
{
    tree->SetBranchAddress(variable.name().data(), &variable.value_, &variable.branch_);
}



const std::map<TString, std::vector<MvaVariablesBase*> >& MvaTreeHandlerBase::stepMvaVariablesMap()const
{
    return m_stepMvaVariables_;
}



void MvaTreeHandlerBase::clear()
{
    selectorList_ = 0;
    
    for(auto& stepMvaVariables : m_stepMvaVariables_){
        MvaVariablesBase::clearVariables(stepMvaVariables.second);
    }
    m_stepMvaVariables_.clear();
}










tth::mvaHelpers::SystematicChannelFileNames tth::mvaHelpers::systematicChannelFileNames(const char* fileListBase,
                                                                                        const std::vector<Channel::Channel>& v_channel,
                                                                                        const std::vector< Systematic::Systematic >& v_systematic,
                                                                                        const bool forTraining)
{
    SystematicChannelFileNames m_systematicChannelFileNames;
    
    for(const auto& systematic : v_systematic){
        std::map<Channel::Channel, std::vector<TString> >& m_channelFileNames = m_systematicChannelFileNames[systematic];
        for(const auto& channel : v_channel){
            std::vector<TString>& v_inputFileName = m_channelFileNames[channel];
            // FIXME: for now systematic is not used to study systematic variations which modify Higgs samples,
            // FIXME: but running on all Higgs masses
            for(const auto& systematicMass : Systematic::allowedSystematicsHiggsPlotting){
                
                // Access FileList containing list of input root files
                // FIXME: almost same functionality as in Samples.cc, unify after MVA training is established
                const TString histoListName(fileListBase + Systematic::convertSystematic(systematicMass) + "_" + Channel::convertChannel(channel) + ".txt");
                //std::cout << "Reading file: " << histoListName << std::endl;
                ifstream fileList(histoListName);
                if(fileList.fail()){
                    std::cerr<<"Error reading file: "<<histoListName<<std::endl;
                    exit(1);
                }
                while(!fileList.eof()){
                    TString filename;
                    fileList>>filename;
                    if(filename==""){continue;} // Skip empty lines
                    if(filename.BeginsWith("#")){continue;} // Comment lines in FileList with '#'
                    
                    if(forTraining && filename.Contains("ttbarH") && filename.Contains("inclusiveBbbar"))
                        v_inputFileName.push_back(filename);
                    else if(!forTraining && filename.Contains("ttbarH") && filename.Contains("tobbbar"))
                        v_inputFileName.push_back(filename);
                    else continue;
                }
            }
        }
    }
    
    return m_systematicChannelFileNames;
}



void MvaTreeHandlerBase::fillVariables(const RecoObjects&, const CommonGenObjects&,
                                       const TopGenObjects&, const HiggsGenObjects&,
                                       const KinRecoObjects&,
                                       const tth::RecoObjectIndices&, const tth::GenObjectIndices&,
                                       const tth::GenLevelWeights&, const tth::RecoLevelWeights&,
                                       const double&, const TString&,
                                       std::vector<MvaVariablesBase*>&)const
{
    // WARNING: this is empty template method, overwrite for inherited class
    
    std::cerr<<"ERROR! Dummy method fillVariables() in MvaTreeHandlerBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



















tth::mvaHelpers::SystematicChannelFileNames tth::mvaHelpers::mergeTrees(
                    const char* mvaInputDir,
                    const tth::mvaHelpers::SystematicChannelFileNames& m_systematicChannelFileNamesTraining,
                    const tth::mvaHelpers::SystematicChannelFileNames& m_systematicChannelFileNamesTesting,
                    const std::vector<std::pair<TString, TString> >& v_nameStepPair)
{
    SystematicChannelFileNames result;
    
    // Loop over all channels and systematics
    for(const auto& systematicChannelFileNamesTraining : m_systematicChannelFileNamesTraining){
        const Systematic::Systematic& systematic = systematicChannelFileNamesTraining.first;
        for(const auto& channelFileNamesTraining : systematicChannelFileNamesTraining.second){
            const Channel::Channel& channel = channelFileNamesTraining.first;
            const std::vector<TString>& v_fileNameTraining = channelFileNamesTraining.second;
            const std::vector<TString>& v_fileNameTesting = m_systematicChannelFileNamesTesting.at(systematic).at(channel);
            std::cout<<"\nProcessing (Channel, Systematic): "<<Channel::convertChannel(channel)<<" , "<<Systematic::convertSystematic(systematic)<<"\n\n";
            
            // Open the input files and access the MVA input training trees
            std::map<TString, TList*> m_stepListTraining;
            for(const auto& fileName : v_fileNameTraining){
                std::cout<<"File for training: "<<fileName<<std::endl;
                // FIXME: need to check whether input file really exists
                TFile* inputFile(0);
                inputFile = TFile::Open(fileName);
                for(const auto& nameStepPair : v_nameStepPair){
                    //std::cout<<"Tree and step: "<<nameStepPair.first<<" , "<<nameStepPair.second<<"\n\n";
                    // FIXME: need to check whether input tree really exists
                    TTree* inputTree = (TTree*)inputFile->Get(nameStepPair.first);
                    if(m_stepListTraining.find(nameStepPair.second) == m_stepListTraining.end()){
                        m_stepListTraining[nameStepPair.second] = new TList;
                    }
                    m_stepListTraining.at(nameStepPair.second)->Add(inputTree);
                }
            }
            std::cout<<std::endl;
            
            // Open the input files and access the MVA input testing trees
            std::map<TString, TList*> m_stepListTesting;
            for(const auto& fileName : v_fileNameTesting){
                std::cout<<"File for testing: "<<fileName<<std::endl;
                // FIXME: need to check whether input file really exists
                TFile* inputFile(0);
                inputFile = TFile::Open(fileName);
                for(const auto& nameStepPair : v_nameStepPair){
                    //std::cout<<"Tree and step: "<<nameStepPair.first<<" , "<<nameStepPair.second<<"\n\n";
                    // FIXME: need to check whether input tree really exists
                    TTree* inputTree = (TTree*)inputFile->Get(nameStepPair.first);
                    if(m_stepListTesting.find(nameStepPair.second) == m_stepListTesting.end()){
                        m_stepListTesting[nameStepPair.second] = new TList;
                    }
                    m_stepListTesting.at(nameStepPair.second)->Add(inputTree);
                }
            }
            std::cout<<std::endl;
            
            // Unfortunately this output file is needed to prevent from strange ROOT message
            TString mergedTreesFileName = common::assignFolder(mvaInputDir, channel, systematic);
            mergedTreesFileName.Append("/");
            mergedTreesFileName.Append("mergedTrees.root");
            TFile* mergedTrees = new TFile(mergedTreesFileName, "RECREATE");
            for(const auto& nameStepPair : v_nameStepPair){
                TTree* treeTraining = TTree::MergeTrees(m_stepListTraining.at(nameStepPair.second));
                treeTraining->SetName("training"+nameStepPair.first);
                TTree* treeTesting = TTree::MergeTrees(m_stepListTesting.at(nameStepPair.second));
                treeTesting->SetName("testing"+nameStepPair.first);
                treeTraining->Write();
                treeTesting->Write();
            }
            mergedTrees->Close();
            result[systematic][channel].push_back(mergedTreesFileName);
        }
    }
    
    return result;
}














