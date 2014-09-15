#define JetCategogies_cxx

#include <sstream>
#include <string>
#include <iostream>
#include <cstdlib>
#include <algorithm>

#include <TString.h>

#include "JetCategories.h"







// --------------------------------- Methods of struct JetCategory ------------------------------------------


JetCategories::JetCategory::JetCategory():
numberOfJets_(0),
numberOfBjets_(0),
greaterEqualJetNumber_(0),
greaterEqualBjetNumber_(0)
{}



JetCategories::JetCategory::JetCategory(const int numberOfJets, const int numberOfBjets,
                                        const bool greaterEqualJetNumber, const bool greaterEqualBjetNumber):
numberOfJets_(numberOfJets),
numberOfBjets_(numberOfBjets),
greaterEqualJetNumber_(greaterEqualJetNumber),
greaterEqualBjetNumber_(greaterEqualBjetNumber)
{}



TString JetCategories::JetCategory::binLabel(const bool rootLabel)const
{
    const std::string greaterEqual(rootLabel ? "#geq" : ">=");
    
    std::stringstream ss_label;
    if(greaterEqualJetNumber_)ss_label<<greaterEqual;
    ss_label<<numberOfJets_<<",";
    if(greaterEqualBjetNumber_)ss_label<<greaterEqual;
    ss_label<<numberOfBjets_;
    
    return ss_label.str().c_str();
}







// --------------------------------- Methods of class JetCategories ------------------------------------------



JetCategories::JetCategories()
{
    this->clear();
}



JetCategories::JetCategories(const int minNumberOfJets, const int maxNumberOfJets,
                             const int minNumberOfBjets, const int maxNumberOfBjets,
                             const bool greaterEqualJetNumber, const bool greaterEqualBjetNumber)
{
    std::cout<<"--- Beginning setting up jet categories\n";
    this->clear();
    
    // Check sanity of input
    if(minNumberOfJets > maxNumberOfJets){
        std::cerr<<"ERROR in JetCategories! Minimum number of jets bigger than maximum number (min, max): "
                 <<minNumberOfJets<<" , "<<maxNumberOfJets<<"...break\n"<<std::endl;
        exit(307);
    }
    if(minNumberOfBjets > maxNumberOfBjets){
        std::cerr<<"ERROR in JetCategories! Minimum number of b-jets bigger than maximum number (min, max): "
                 <<minNumberOfBjets<<" , "<<maxNumberOfBjets<<"...break\n"<<std::endl;
        exit(308);
    }
    if(minNumberOfJets < minNumberOfBjets){
        std::cerr<<"ERROR in JetCategories! Minimum number of jets smaller than minimum number of b-jets (jets, b-jets): "
                 <<minNumberOfJets<<" , "<<minNumberOfBjets<<"...break\n"<<std::endl;
        exit(309);
    }
    if(maxNumberOfJets < maxNumberOfBjets && !greaterEqualJetNumber){
        std::cerr<<"ERROR in JetCategories! Maximum number of jets smaller than maximum number of b-jets (jets, b-jets): "
                 <<maxNumberOfJets<<" , "<<maxNumberOfBjets<<"...break\n"<<std::endl;
        exit(310);
    }
    
    // Add all categories by looping over requested number of jets and number of b-jets
    for(int nBjets = minNumberOfBjets; nBjets <= maxNumberOfBjets; ++nBjets){
        const int maxNumberOfJetsForNumberOfBjets = greaterEqualJetNumber ? std::max(maxNumberOfJets, nBjets) : maxNumberOfJets;
        const int minNumberOfJetsForNumberOfBjets = greaterEqualJetNumber ? std::max(minNumberOfJets, nBjets) : minNumberOfJets;
        for(int nJets = minNumberOfJetsForNumberOfBjets; nJets <= maxNumberOfJetsForNumberOfBjets; ++nJets){
            if(nJets<maxNumberOfJets && nBjets<maxNumberOfBjets){
                this->addCategory(nJets, nBjets, false, false);
            }
            else if(nJets<maxNumberOfJets && nBjets==maxNumberOfBjets){
                const bool greaterEqual = nJets>maxNumberOfBjets ? greaterEqualBjetNumber : false;
                this->addCategory(nJets, nBjets, false, greaterEqual);
            }
            else if(nJets>=maxNumberOfJets && nBjets<maxNumberOfBjets){
                this->addCategory(nJets, nBjets, greaterEqualJetNumber, false);
            }
            else if(nJets>=maxNumberOfJets && nBjets>=maxNumberOfBjets){
                this->addCategory(nJets, nBjets, greaterEqualJetNumber, greaterEqualBjetNumber);
            }
        }
    }
    std::cout<<"Number of categories: "<<this->numberOfCategories()<<"\n";
    std::cout<<"Labels:\t";
    for(const auto label : this->binLabels(false)) std::cout<<label<<"  ;  ";
    std::cout<<std::endl;
        
    std::cout<<"=== Finishing setting up jet categories\n\n";
}



JetCategories::~JetCategories(){}



void JetCategories::addCategory(const int numberOfJets, const int numberOfBjets,
                                const bool greaterEqualJetNumber, const bool greaterEqualBjetNumber)
{
    const JetCategory jetCategory(numberOfJets, numberOfBjets, greaterEqualJetNumber, greaterEqualBjetNumber);
    if(ambiguousCategory(jetCategory)){
        std::cout<<"... Category is not added\n\n";
        return;
    }
    for(std::vector<JetCategory>::iterator i_jetCategory = v_jetCategory_.begin(); i_jetCategory != v_jetCategory_.end(); ++i_jetCategory){
        if(numberOfJets<i_jetCategory->numberOfJets_){
            v_jetCategory_.insert(i_jetCategory, jetCategory);
            return;
        }
        else if(numberOfJets==i_jetCategory->numberOfJets_ && numberOfBjets<i_jetCategory->numberOfBjets_){
            v_jetCategory_.insert(i_jetCategory, jetCategory);
            return;
        }
    }
    v_jetCategory_.push_back(jetCategory);
}



bool JetCategories::ambiguousCategory(const JetCategories::JetCategory& jetCategory)const
{
    if(jetCategory.numberOfJets_<0 || jetCategory.numberOfBjets_<0){
        std::cerr<<"\nBad Jet Category, negative number of (b-)jets: "<<jetCategory.binLabel()
                 <<"\n...break\n\n";
        exit(311);
    }
    if(jetCategory.numberOfJets_<jetCategory.numberOfBjets_){
        std::cerr<<"\nBad Jet Category, number of jets smaller than number of b-jets: "<<jetCategory.binLabel()
                 <<"\n...break\n\n";
        exit(312);
    }
    if(jetCategory.numberOfJets_==jetCategory.numberOfBjets_ && !jetCategory.greaterEqualJetNumber_ && jetCategory.greaterEqualBjetNumber_){
        std::cerr<<"\nBad Jet Category, setting >= for number of b-jets makes no sense: "<<jetCategory.binLabel()
                 <<"\n...break\n\n";
        exit(313);
    }
    
    for(std::vector<JetCategory>::const_iterator i_jetCategory = v_jetCategory_.begin(); i_jetCategory != v_jetCategory_.end(); ++i_jetCategory){
        if(hasOverlap(jetCategory.numberOfJets_, jetCategory.greaterEqualJetNumber_, 
                      i_jetCategory->numberOfJets_, i_jetCategory->greaterEqualJetNumber_) && 
           hasOverlap(jetCategory.numberOfBjets_, jetCategory.greaterEqualBjetNumber_, 
                      i_jetCategory->numberOfBjets_, i_jetCategory->greaterEqualBjetNumber_)){
            std::cout<<"Ambiguity while adding Jet Category: "<<jetCategory.binLabel(false)
                     <<"\nOverlap with existing category: "<<i_jetCategory->binLabel(false)
                     <<"\n";
            return true;
        }
    }
    return false;
}



bool JetCategories::hasOverlap(const int numberOfJets1, const bool greaterEqualJetNumber1,
                               const int numberOfJets2, const bool greaterEqualJetNumber2)const
{
    if(numberOfJets1 == numberOfJets2)return true;
    if(greaterEqualJetNumber1 && numberOfJets1<numberOfJets2)return true;
    if(greaterEqualJetNumber2 && numberOfJets1>numberOfJets2)return true;
    return false;
}



void JetCategories::clear()
{
    v_jetCategory_.clear();
}



int JetCategories::numberOfCategories()const
{
    return v_jetCategory_.size();
}



int JetCategories::categoryId(const int numberOfJets, const int numberOfBjets)const
{
    // Events belonging to specific category should get value >=0
    // If the category definition is bad, it should end up in the underflow (return value <0, depending on error)
    // If the category definition is fine, but the event does not belong to any category,
    // it should end up in the overflow (return value +999)
    
    if(v_jetCategory_.size()==0)return -1;
    bool alreadyInCategory(false);
    int position(-99);
    std::vector<JetCategory>::const_iterator i_jetCategory;
    for(i_jetCategory = v_jetCategory_.begin(); i_jetCategory != v_jetCategory_.end(); ++i_jetCategory){
        if(numberOfJets<i_jetCategory->numberOfJets_)continue;
        if(!i_jetCategory->greaterEqualJetNumber_ && numberOfJets>i_jetCategory->numberOfJets_)continue;
        if(numberOfBjets<i_jetCategory->numberOfBjets_)continue;
        if(!i_jetCategory->greaterEqualBjetNumber_ && numberOfBjets>i_jetCategory->numberOfBjets_)continue;
        if(alreadyInCategory){
            std::cerr<<"Bad jet category definition, event with (# jets = "<<numberOfJets
                     <<", # bjets = "<<numberOfBjets<<") fits in more than one category!!!\n\n";
            exit(314);
        }
        position = std::distance(v_jetCategory_.begin(), i_jetCategory);
        alreadyInCategory = true;
    }
    // is event contained in any category?
    if(!alreadyInCategory)return +999;
    return position;
}



std::vector<TString> JetCategories::binLabels(const bool rootLabel)const
{
    std::vector<TString> v_label;
    for(std::vector<JetCategory>::const_iterator i_jetCategory = v_jetCategory_.begin(); i_jetCategory != v_jetCategory_.end(); ++i_jetCategory){
        v_label.push_back(i_jetCategory->binLabel(rootLabel));
    }
    return v_label;
}














