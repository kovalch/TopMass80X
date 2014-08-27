#include <iostream>
#include <iomanip>

#include "KinematicReconstructionSolution.h"







// -------------------------------------- Methods for KinematicReconstructionSolution --------------------------------------


KinematicReconstructionSolution::KinematicReconstructionSolution():
leptonIndex_(-1),
antiLeptonIndex_(-1),
bjetIndex_(-1),
antiBjetIndex_(-1),
wPlus_(LV(0.,0.,0.,0.)),
wMinus_(LV(0.,0.,0.,0.)),
top_(LV(0.,0.,0.,0.)),
antiTop_(LV(0.,0.,0.,0.)),
neutrino_(LV(0.,0.,0.,0.)),
antiNeutrino_(LV(0.,0.,0.,0.)),
reconstructedTopMass_(-999.),
numberOfBtags_(-1),
weight_averaged_sumSmearings_mlbBased_(-999.)
{}



KinematicReconstructionSolution::KinematicReconstructionSolution(const KinematicReconstructionSolution& kinematicReconstructionSolution):
leptonIndex_(kinematicReconstructionSolution.leptonIndex_),
antiLeptonIndex_(kinematicReconstructionSolution.antiLeptonIndex_),
bjetIndex_(kinematicReconstructionSolution.bjetIndex_),
antiBjetIndex_(kinematicReconstructionSolution.antiBjetIndex_),
wPlus_(kinematicReconstructionSolution.wPlus_),
wMinus_(kinematicReconstructionSolution.wMinus_),
top_(kinematicReconstructionSolution.top_),
antiTop_(kinematicReconstructionSolution.antiTop_),
neutrino_(kinematicReconstructionSolution.neutrino_),
antiNeutrino_(kinematicReconstructionSolution.antiNeutrino_),
reconstructedTopMass_(kinematicReconstructionSolution.reconstructedTopMass_),
numberOfBtags_(kinematicReconstructionSolution.numberOfBtags_),
weight_averaged_sumSmearings_mlbBased_(kinematicReconstructionSolution.weight_averaged_sumSmearings_mlbBased_)
{}



KinematicReconstructionSolution::KinematicReconstructionSolution(const int leptonIndex, const int antiLeptonIndex,
                                                                 const int bjetIndex, const int antiBjetIndex,
                                                                 const LV& wPlus, const LV& wMinus,
                                                                 const LV& top, const LV& antiTop,
                                                                 const LV& neutrino, const LV& antiNeutrino,
                                                                 const double& reconstructedTopMass,
                                                                 const int numberOfBtags,
                                                                 const double& weight_averaged_sumSmearings_mlbBased):
leptonIndex_(leptonIndex),
antiLeptonIndex_(antiLeptonIndex),
bjetIndex_(bjetIndex),
antiBjetIndex_(antiBjetIndex),
wPlus_(wPlus),
wMinus_(wMinus),
top_(top),
antiTop_(antiTop),
neutrino_(neutrino),
antiNeutrino_(antiNeutrino),
reconstructedTopMass_(reconstructedTopMass),
numberOfBtags_(numberOfBtags),
weight_averaged_sumSmearings_mlbBased_(weight_averaged_sumSmearings_mlbBased)
{}



void KinematicReconstructionSolution::print()const
{
    const LV ttbar = this->ttbar();
    std::cout<<"Index (lepton, antilepton):         "<<leptonIndex_<<" , "<<antiLeptonIndex_<<"\n"
             <<"Index (b jet, anti-b jet):          "<<bjetIndex_<<" , "<<antiBjetIndex_<<"\n"
             <<"W+ (pt, eta, phi, mass):            "<<std::fixed<<std::setprecision(3)<<wPlus_.pt()<<" , "<<wPlus_.eta()<<" , "<<wPlus_.phi()<<" , "<<wPlus_.mass()<<"\n"
             <<"W- (pt, eta, phi, mass):            "<<wMinus_.pt()<<" , "<<wMinus_.eta()<<" , "<<wMinus_.phi()<<" , "<<wMinus_.mass()<<"\n"
             <<"Top (pt, eta, phi, mass):           "<<top_.pt()<<" , "<<top_.eta()<<" , "<<top_.phi()<<" , "<<top_.mass()<<"\n"
             <<"Anti-top (pt, eta, phi, mass):      "<<antiTop_.pt()<<" , "<<antiTop_.eta()<<" , "<<antiTop_.phi()<<" , "<<antiTop_.mass()<<"\n"
             <<"tt system(pt, eta, phi, mass):      "<<ttbar.pt()<<" , "<<ttbar.eta()<<" , "<<ttbar.phi()<<" , "<<ttbar.mass()<<"\n"
             <<"Neutrino (pt, eta, phi, mass):      "<<neutrino_.pt()<<" , "<<neutrino_.eta()<<" , "<<neutrino_.phi()<<" , "<<neutrino_.mass()<<"\n"
             <<"Anti-neutrino (pt, eta, phi, mass): "<<antiNeutrino_.pt()<<" , "<<antiNeutrino_.eta()<<" , "<<antiNeutrino_.phi()<<" , "<<antiNeutrino_.mass()<<"\n"
             <<"Reconstructed top mass:             "<<reconstructedTopMass_<<"\n"
             <<"Number of b-tags:                   "<<numberOfBtags_<<"\n"
             <<"Weight:                             "<<weight_averaged_sumSmearings_mlbBased_<<"\n";
}










// -------------------------------------- Methods for KinematicReconstructionSolutions --------------------------------------


KinematicReconstructionSolutions::KinematicReconstructionSolutions()
{}




void KinematicReconstructionSolutions::addSolutionTwoBtags(const KinematicReconstructionSolution& solution)
{
    // Add solution
    v_solution_.push_back(solution);
    const std::vector<KinematicReconstructionSolution>::const_iterator i_solution = v_solution_.end() - 1;
    v_solutionTwoBtags_.push_back(&*i_solution);
    
    // Get positions in vectors
    const int solutionIndex = std::distance(v_solution_.begin(), v_solution_.end()) - 1;
    const int solutionCateogryIndex = std::distance(v_solutionTwoBtags_.begin(), v_solutionTwoBtags_.end()) - 1;
    
    // Add index of solution, ordered by weight
    this->insertIndex(solutionIndex, v_solution_, solution.weight_averaged_sumSmearings_mlbBased_,
                      v_index_weight_averaged_sumSmearings_mlbBased_);
    this->insertIndex(solutionCateogryIndex, v_solutionTwoBtags_, solution.weight_averaged_sumSmearings_mlbBased_,
                      v_indexTwoBtags_weight_averaged_sumSmearings_mlbBased_);
}



void KinematicReconstructionSolutions::addSolutionOneBtag(const KinematicReconstructionSolution& solution)
{
    v_solution_.push_back(solution);
    const std::vector<KinematicReconstructionSolution>::const_iterator i_solution = v_solution_.end() - 1;
    v_solutionOneBtag_.push_back(&*i_solution);
    
    // Get positions in vectors
    const int solutionIndex = std::distance(v_solution_.begin(), v_solution_.end()) - 1;
    const int solutionCateogryIndex = std::distance(v_solutionOneBtag_.begin(), v_solutionOneBtag_.end()) - 1;
    
    // Add index of solution, ordered by weight
    this->insertIndex(solutionIndex, v_solution_, solution.weight_averaged_sumSmearings_mlbBased_,
                      v_index_weight_averaged_sumSmearings_mlbBased_);
    this->insertIndex(solutionCateogryIndex, v_solutionOneBtag_, solution.weight_averaged_sumSmearings_mlbBased_,
                      v_indexOneBtag_weight_averaged_sumSmearings_mlbBased_);
}



void KinematicReconstructionSolutions::addSolutionNoBtags(const KinematicReconstructionSolution& solution)
{
    v_solution_.push_back(solution);
    const std::vector<KinematicReconstructionSolution>::const_iterator i_solution = v_solution_.end() - 1;
    v_solutionNoBtags_.push_back(&*i_solution);
    
    // Get positions in vectors
    const int solutionIndex = std::distance(v_solution_.begin(), v_solution_.end()) - 1;
    const int solutionCateogryIndex = std::distance(v_solutionNoBtags_.begin(), v_solutionNoBtags_.end()) - 1;
    
    // Add index of solution, ordered by weight
    this->insertIndex(solutionIndex, v_solution_, solution.weight_averaged_sumSmearings_mlbBased_,
                      v_index_weight_averaged_sumSmearings_mlbBased_);
    this->insertIndex(solutionCateogryIndex, v_solutionNoBtags_, solution.weight_averaged_sumSmearings_mlbBased_,
                      v_indexNoBtags_weight_averaged_sumSmearings_mlbBased_);
}



void KinematicReconstructionSolutions::insertIndex(const size_t index, const std::vector<KinematicReconstructionSolution>& v_solution,
                                                   const double weight, std::vector<size_t>& v_index)
{
    if(!v_index.size()){
        v_index.push_back(index);
        return;
    }
    
    bool isInserted(false);
    for(std::vector<size_t>::iterator i_index = v_index.begin(); i_index != v_index.end(); ++i_index){
        const double& iWeight = v_solution.at(*i_index).weight_averaged_sumSmearings_mlbBased_;
        if(iWeight < weight){
            v_index.insert(i_index, index);
            isInserted = true;
            break;
        }
    }
    if(!isInserted) v_index.push_back(index);
}



void KinematicReconstructionSolutions::insertIndex(const size_t index, const std::vector<const KinematicReconstructionSolution*>& v_solution,
                                                   const double weight, std::vector<size_t>& v_index)
{
    if(!v_index.size()){
        v_index.push_back(index);
        return;
    }
    
    bool isInserted(false);
    for(std::vector<size_t>::iterator i_index = v_index.begin(); i_index != v_index.end(); ++i_index){
        const double& iWeight = v_solution.at(*i_index)->weight_averaged_sumSmearings_mlbBased_;
        if(iWeight < weight){
            v_index.insert(i_index, index);
            isInserted = true;
            break;
        }
    }
    if(!isInserted) v_index.push_back(index);
}



const KinematicReconstructionSolution& KinematicReconstructionSolutions::solution_averaged_sumSmearings_mlbBased(const size_t solutionNumber)const
{
    const size_t index = v_index_weight_averaged_sumSmearings_mlbBased_.at(solutionNumber);
    return v_solution_.at(index);
}



const KinematicReconstructionSolution& KinematicReconstructionSolutions::solutionTwoBtags_averaged_sumSmearings_mlbBased(const size_t solutionNumber)const
{
    const size_t index = v_indexTwoBtags_weight_averaged_sumSmearings_mlbBased_.at(solutionNumber);
    return *(v_solutionTwoBtags_.at(index));
}



const KinematicReconstructionSolution& KinematicReconstructionSolutions::solutionOneBtag_averaged_sumSmearings_mlbBased(const size_t solutionNumber)const
{
    const size_t index = v_indexOneBtag_weight_averaged_sumSmearings_mlbBased_.at(solutionNumber);
    return *(v_solutionOneBtag_.at(index));
}


const KinematicReconstructionSolution& KinematicReconstructionSolutions::solutionNoBtags_averaged_sumSmearings_mlbBased(const size_t solutionNumber)const
{
    const size_t index = v_indexNoBtags_weight_averaged_sumSmearings_mlbBased_.at(solutionNumber);
    return *(v_solutionNoBtags_.at(index));
}





