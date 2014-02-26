#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <fstream>

#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TSelectorList.h>
#include <TIterator.h>
#include <TObject.h>
#include <TList.h>
#include <Rtypes.h>

#include "MvaTreeHandler.h"
#include "MvaVariablesTopJets.h"
#include "analysisStructs.h"
#include "JetCategories.h"
#include "higgsUtils.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/analysisObjectStructs.h"





MvaTreeHandler::MvaTreeHandler(const char* mvaInputDir,
                               const std::vector<TString>& selectionStepsNoCategories,
                               const std::vector<TString>& stepsForCategories,
                               const JetCategories* jetCategories):
selectorList_(0),
selectionSteps_(selectionStepsNoCategories),
stepsForCategories_(stepsForCategories),
jetCategories_(jetCategories),
mvaInputDir_(mvaInputDir)
{
    std::cout<<"--- Beginning setting up MVA input tree handler\n";
    
    if(!jetCategories_ && stepsForCategories_.size()>0){
        std::cerr<<"ERROR in constructor for MvaTreeHandler! "
                 <<"No jet categories passed, but request for category-wise selection steps\n...break\n"<<std::endl;
        exit(236);
    }
    if(jetCategories_ && stepsForCategories_.size()==0){
        std::cerr<<"ERROR in constructor for MvaTreeHandler! "
                 <<"Jet categories passed, but no category-wise selection steps defined\n...break\n"<<std::endl;
        exit(237);
    }
    
    std::cout<<"=== Finishing setting up MVA input tree handler\n\n";
}



void MvaTreeHandler::book()
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



void MvaTreeHandler::addStep(const TString& step)
{
    // Check whether step already exists
    if(this->checkExistence(step)){
        std::cout<<"Warning in addStep()! Selection step already contained: "<<step
                 <<"\n...skip this one\n";
        return;
    }

    // Book the step
    m_stepMvaVariables_[step] = std::vector<MvaVariablesTopJets*>();
}



bool MvaTreeHandler::checkExistence(const TString& step)const
{
    return m_stepMvaVariables_.find(step) != m_stepMvaVariables_.end();
}



void MvaTreeHandler::fill(const RecoObjects& recoObjects,
                          const tth::GenObjectIndices& genObjectIndices, const tth::RecoObjectIndices& recoObjectIndices,
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
        // Here check the individual jet categories
        const int categoryId = jetCategories_->categoryId(numberOfJets, numberOfBjets);
        const TString fullStepName = tth::stepName(stepShort, categoryId);
        this->fill(recoObjects, genObjectIndices, recoObjectIndices, weight, fullStepName);
    }
    if(!stepExists) return;
    
    // Loop over all jet combinations and get MVA input variables
    const std::vector<MvaVariablesTopJets*> v_mvaTopJetsVariables = 
            MvaVariablesTopJets::fillVariables(recoObjectIndices, genObjectIndices, recoObjects, weight);
    
    // Fill the MVA variables
    m_stepMvaVariables_.at(step).insert(m_stepMvaVariables_.at(step).end(), v_mvaTopJetsVariables.begin(), v_mvaTopJetsVariables.end());
}



void MvaTreeHandler::writeTrees(const std::string& outputFilename,
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



void MvaTreeHandler::writeTrees(TSelectorList* output)
{
    std::cout<<"--- Beginning production of MVA input trees\n";
    
    // Set pointer to output, so that TTrees are owned by it
    selectorList_ = output;
    
    std::map<TString, TTree*> m_stepTree;
    for(const auto& stepMvaVariables : m_stepMvaVariables_){
        const TString& step = stepMvaVariables.first;
        const std::vector<MvaVariablesTopJets*>& v_mvaVariables = stepMvaVariables.second;
        TTree* tree = m_stepTree[step];
        tree = this->store(new TTree("mvaInputTopJets"+step, "mvaInputTopJets"));
        
        // Set up struct for branch handling
        MvaVariablesTopJets* const mvaVariablesTopJets = new MvaVariablesTopJets();
        this->createBranches(tree, mvaVariablesTopJets);
        this->fillBranches(tree, mvaVariablesTopJets, v_mvaVariables);
        delete mvaVariablesTopJets;
    }
    
    std::cout<<"Dijet combinations per step (step, no. of combinations):\n";
    for(auto vars : m_stepMvaVariables_){
        std::cout<<"\t"<<vars.first<<" , "<<vars.second.size()<<"\n";
    }
    
    std::cout<<"=== Finishing production of MVA input trees\n\n";
}



void MvaTreeHandler::createBranches(TTree* tree, const MvaVariablesTopJets* mvaTopJetsVariables)const
{
    this->createBranch(tree, mvaTopJetsVariables->lastInEvent_);
    this->createBranch(tree, mvaTopJetsVariables->eventWeight_);
    this->createBranch(tree, mvaTopJetsVariables->bQuarkRecoJetMatched_);
    this->createBranch(tree, mvaTopJetsVariables->correctCombination_);
    this->createBranch(tree, mvaTopJetsVariables->swappedCombination_);
    this->createBranch(tree, mvaTopJetsVariables->jetChargeDiff_);
    this->createBranch(tree, mvaTopJetsVariables->meanDeltaPhi_b_met_);
    this->createBranch(tree, mvaTopJetsVariables->massDiff_recoil_bbbar_);
    this->createBranch(tree, mvaTopJetsVariables->pt_b_antiLepton_);
    this->createBranch(tree, mvaTopJetsVariables->pt_antiB_lepton_);
    this->createBranch(tree, mvaTopJetsVariables->deltaR_b_antiLepton_);
    this->createBranch(tree, mvaTopJetsVariables->deltaR_antiB_lepton_);
    this->createBranch(tree, mvaTopJetsVariables->btagDiscriminatorSum_);
    this->createBranch(tree, mvaTopJetsVariables->deltaPhi_antiBLepton_bAntiLepton_);
    this->createBranch(tree, mvaTopJetsVariables->massDiff_fullBLepton_bbbar_);
    this->createBranch(tree, mvaTopJetsVariables->meanMt_b_met_);
    this->createBranch(tree, mvaTopJetsVariables->massSum_antiBLepton_bAntiLepton_);
    this->createBranch(tree, mvaTopJetsVariables->massDiff_antiBLepton_bAntiLepton_);
}



void MvaTreeHandler::createBranch(TTree* tree, const MvaVariableInt& variable)const
{
    std::string name(variable.name());
    std::string nameType(name);
    nameType.append("/").append(variable.type());
    tree->Branch(name.data(), (Long_t)&variable.value_, nameType.data());
}



void MvaTreeHandler::createBranch(TTree* tree, const MvaVariableFloat& variable)const
{
    std::string name(variable.name());
    std::string nameType(name);
    nameType.append("/").append(variable.type());
    tree->Branch(name.c_str(), (Long_t)&variable.value_, nameType.data());
}



void MvaTreeHandler::fillBranches(TTree* tree, MvaVariablesTopJets* const mvaTopJetsVariables,
                                  const std::vector<MvaVariablesTopJets*>& v_mvaTopJetsVariables)const
{
    for(const MvaVariablesTopJets* mvaVariables : v_mvaTopJetsVariables){
        *mvaTopJetsVariables = *mvaVariables;
        tree->Fill();
    }
}



void MvaTreeHandler::clear()
{
    selectorList_ = 0;
    
    for(auto& stepMvaVariables : m_stepMvaVariables_){
        MvaVariablesTopJets::clearVariables(stepMvaVariables.second);
    }
    m_stepMvaVariables_.clear();
}



void MvaTreeHandler::importTrees(const std::string& f_savename, const std::string& prefix)
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
    m_stepMvaVariables_.clear();
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



void MvaTreeHandler::importBranches(TTree* tree, std::vector<MvaVariablesTopJets*>& v_mvaTopJetsVariables)
{
    // Set up variables struct
    MvaVariablesTopJets mvaTopJetsVariables;
    
    // Set branch addresses
    this->importBranch(tree, mvaTopJetsVariables.lastInEvent_);
    this->importBranch(tree, mvaTopJetsVariables.eventWeight_);
    this->importBranch(tree, mvaTopJetsVariables.bQuarkRecoJetMatched_);
    this->importBranch(tree, mvaTopJetsVariables.correctCombination_);
    this->importBranch(tree, mvaTopJetsVariables.swappedCombination_);
    this->importBranch(tree, mvaTopJetsVariables.jetChargeDiff_);
    this->importBranch(tree, mvaTopJetsVariables.meanDeltaPhi_b_met_);
    this->importBranch(tree, mvaTopJetsVariables.massDiff_recoil_bbbar_);
    this->importBranch(tree, mvaTopJetsVariables.pt_b_antiLepton_);
    this->importBranch(tree, mvaTopJetsVariables.pt_antiB_lepton_);
    this->importBranch(tree, mvaTopJetsVariables.deltaR_b_antiLepton_);
    this->importBranch(tree, mvaTopJetsVariables.deltaR_antiB_lepton_);
    this->importBranch(tree, mvaTopJetsVariables.btagDiscriminatorSum_);
    this->importBranch(tree, mvaTopJetsVariables.deltaPhi_antiBLepton_bAntiLepton_);
    this->importBranch(tree, mvaTopJetsVariables.massDiff_fullBLepton_bbbar_);
    this->importBranch(tree, mvaTopJetsVariables.meanMt_b_met_);
    this->importBranch(tree, mvaTopJetsVariables.massSum_antiBLepton_bAntiLepton_);
    this->importBranch(tree, mvaTopJetsVariables.massDiff_antiBLepton_bAntiLepton_);
    
    // Loop over all tree entries and fill vector of structs
    for(Long64_t iEntry = 0; iEntry < tree->GetEntries(); ++iEntry){
        tree->GetEntry(iEntry);
        MvaVariablesTopJets* const mvaVariablesTopJetsPtr = new MvaVariablesTopJets();
        *mvaVariablesTopJetsPtr = mvaTopJetsVariables;
        v_mvaTopJetsVariables.push_back(mvaVariablesTopJetsPtr);
    }
}



void MvaTreeHandler::importBranch(TTree* tree, MvaVariableInt& variable)
{
    tree->SetBranchAddress(variable.name().data(), &variable.value_, &variable.branch_);
}



void MvaTreeHandler::importBranch(TTree* tree, MvaVariableFloat& variable)
{
    tree->SetBranchAddress(variable.name().data(), &variable.value_, &variable.branch_);
}



const std::map<TString, std::vector<MvaVariablesTopJets*> >& MvaTreeHandler::stepMvaVariablesMap()const
{
    return m_stepMvaVariables_;
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














