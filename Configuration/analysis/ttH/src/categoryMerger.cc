#include <string>
#include <vector>
#include <map>
#include <fstream>

#include <TString.h>
#include <TFile.h>
#include <TObject.h>
#include <TH1.h>
#include <TTree.h>
#include <TList.h>

#include "higgsUtils.h"
#include "../../common/include/utils.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/CommandLineParameters.h"




void categoryMerger(const std::vector<Channel::Channel>& v_channel,
                    const std::vector<Systematic::Systematic>& v_systematic,
                    const TString& filelistDir,
                    const std::vector<int>& v_jetCategory,
                    std::vector<TString>& v_pattern)
{
    for(const auto& systematic : v_systematic){
        for(const auto& channel : v_channel){
    
            // Read the full input filenames from the filelist
            const std::vector<TString> v_filename = common::readFilelist(filelistDir, channel, systematic, v_pattern);
            
            // Loop over all selected files and merge categories
            for(const TString& filename : v_filename){
                
                // Check if file exists
                std::ifstream filestream;
                filestream.open(filename);
                if(filestream.is_open()) filestream.close();
                else{
                    std::cerr<<"Error in categoryMerger()! Root file not found: "<<filename<<"\n...break\n"<<std::endl;
                    exit(28);
                }
                
                // Create step name with merged categories, and check if merged categories already exist
                const TString mergedCategoriesName = tth::categoryName(v_jetCategory);
                const std::vector<std::pair<TString, TString> > v_nameStepPairMerged = tth::nameStepPairs(filename, mergedCategoriesName);
                for(const auto& nameStepPair : v_nameStepPairMerged){
                    if(static_cast<int>(v_jetCategory.size()) == tth::numberOfCategories(nameStepPair.second)){
                        std::cout<<"Merged categories already exist in file: "<<filename
                                 <<"\n...not merging again, but stopping\n"<<std::endl;
                        exit(0);
                    }
                }
                
                // Check if all requested categories exist
                std::map<int, std::vector<std::pair<TString, TString> > > m_categoryNameStepPair;
                for(const int category : v_jetCategory){
                    m_categoryNameStepPair[category] = std::vector<std::pair<TString, TString> >();
                    std::vector<std::pair<TString, TString> > result;
                    const TString categoryName = tth::categoryName(category);
                    const std::vector<std::pair<TString, TString> > v_nameStepPair = tth::nameStepPairs(filename, categoryName);
                    for(const auto& nameStepPair : v_nameStepPair){
                        if(tth::numberOfCategories(nameStepPair.second) != 1) continue;
                        m_categoryNameStepPair.at(category).push_back(nameStepPair);
                    }
                    if(!m_categoryNameStepPair.at(category).size()){
                        std::cerr<<"Error in categoryMerger()! Requested category not found: "
                                 <<category<<"\n...break\n"<<std::endl;
                        exit(913);
                    }
                }
                
                // Check if all categories have same amount of objects in root file
                size_t numberOfObjects = m_categoryNameStepPair.begin()->second.size();
                for(auto i_pair = ++m_categoryNameStepPair.begin(); i_pair != m_categoryNameStepPair.end(); ++i_pair){
                    if(i_pair->second.size() != numberOfObjects){
                        std::cerr<<"Error in categoryMerger()! Different numbers of objects for different categories found"
                                 <<"\n...break\n"<<std::endl;
                        exit(914);
                    }
                }
                
                // Open file with write permissions, create merged histos and trees, and write to file
                TFile file(filename, "UPDATE");
                std::vector<TObject*> v_object;
                std::vector<std::pair<TString, TString> > v_nameStepPair = m_categoryNameStepPair.at(v_jetCategory.at(0));
                for(size_t iObject = 0; iObject < v_nameStepPair.size(); ++iObject){
                    TObject* object = file.Get(v_nameStepPair.at(iObject).first);
                    const TString categoryName = tth::extractJetCategory(v_nameStepPair.at(iObject).second);
                    TString objectName = object->GetName();
                    objectName.ReplaceAll(categoryName, mergedCategoriesName);
                    TH1* hist = dynamic_cast<TH1*>(object);
                    if(hist){
                        for(size_t iCategory = 1; iCategory < v_jetCategory.size(); ++iCategory){
                            TH1* tmpHist = (TH1*)file.Get(m_categoryNameStepPair.at(v_jetCategory.at(iCategory)).at(iObject).first);
                            hist->Add(tmpHist);
                        }
                        hist->SetName(objectName);
                        hist->Write();
                        continue;
                    }
                    TTree* tree = dynamic_cast<TTree*>(object);
                    if(tree){
                        TList* list = new TList;
                        list->Add(tree);
                        for(size_t iCategory = 1; iCategory < v_jetCategory.size(); ++iCategory){
                            TTree* tmpTree = (TTree*)file.Get(m_categoryNameStepPair.at(v_jetCategory.at(iCategory)).at(iObject).first);
                            list->Add(tmpTree);
                        }
                        TTree* outputTree = TTree::MergeTrees(list);
                        outputTree->SetName(objectName);
                        outputTree->Write();
                        list->Delete();
                        continue;
                    }
                }
                
                // Close file
                file.Close();
            }
        }
    }
}



/// All systematics allowed for merging, i.e. all which come out of load_Analysis
namespace Systematic{
    const std::vector<Type> allowedSystematics = {
        nominal, all, allAvailable,
        pu, lept, trig,
        jer, jes,
        btag, 
        btagPt, btagEta,
        btagLjet, 
        btagLjetPt, btagLjetEta,
        btagDiscrBstat1, btagDiscrBstat2,
        btagDiscrLstat1, btagDiscrLstat2,
        btagDiscrPurity,
        kin,
        topPt,
        mass, match, scale,
        powheg, powhegHerwig, mcatnlo, perugia11, perugia11NoCR,
        pdf
    };
}



int main(int argc, char** argv){
    
    // Get and check configuration parameters
    CLParameter<std::string> opt_fileDir("d", "Directory where filelist to be processed is contained", true, 1, 1);
    CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu. Default: all", false, 1, 3,
        common::makeStringCheck(Channel::convert(Channel::allowedChannelsAnalysis)));
    CLParameter<std::string> opt_systematic("s", "Systematic variation. Default: Nominal", false, 1, 100,
        common::makeStringCheckBegin(Systematic::convertType(Systematic::allowedSystematics)));
    CLParameter<int> opt_jetCategory("j", "Specify categories to be joined (at least two)", true, 2, 100,
        [](int id){return id >= 0;});
    CLParameter<std::string> opt_filelistPattern("f", "Name (pattern) of samples in filelist to be processed", false, 1, 100);
    CLAnalyser::interpretGlobal(argc, argv);
    
    // Set up filelist to read
    const TString filelistDir = opt_fileDir[0];
    std::cout<<"Processing file list in directory: "<<filelistDir<<"\n\n";
    
    // Set up patterns of root files for samples which should be processed
    std::vector<TString> v_pattern;
    if(opt_filelistPattern.isSet()){
        for(auto argument : opt_filelistPattern.getArguments()) v_pattern.push_back(argument);
        std::cout<<"Restrict to samples containing in name: ";
        for(auto pattern : v_pattern) std::cout<<pattern<<" ";
        std::cout<<"\n\n";
    }
    
    // Set up channels
    std::vector<Channel::Channel> v_channel(Channel::allowedChannelsPlotting);
    if(opt_channel.isSet()) v_channel = Channel::convert(opt_channel.getArguments());
    std::cout<<"Processing channels: ";
    for(auto channel : v_channel) std::cout<<Channel::convert(channel)<<" ";
    std::cout<<"\n\n";
    
    // Set up systematics
    std::vector<Systematic::Systematic> v_systematic;
    if(opt_systematic.isSet()){
        if(opt_systematic[0] == Systematic::convertType(Systematic::all))
            v_systematic = Systematic::allowedSystematicsAnalysis(Systematic::allowedSystematics);
        else if(opt_systematic[0] == Systematic::convertType(Systematic::allAvailable))
            v_systematic = common::findSystematicsFromFilelists(filelistDir, v_channel, v_systematic);
        else
            v_systematic = Systematic::setSystematics(opt_systematic.getArguments());
    }
    else{
        v_systematic.clear();
        v_systematic.push_back(Systematic::nominalSystematic());
    }
    std::cout <<"Processing systematics: "; 
    for(auto systematic : v_systematic) std::cout<<systematic.name()<<" ";
    std::cout<<"\n\n";
    
    const std::vector<int> v_jetCategory = opt_jetCategory.getArguments();
    std::cout<<"Categories to merge: ";
    for(auto jetCategory : v_jetCategory) std::cout<<jetCategory<<" ";
    std::cout<<"\n\n";
    
    categoryMerger(v_channel, v_systematic, filelistDir, v_jetCategory, v_pattern);
}