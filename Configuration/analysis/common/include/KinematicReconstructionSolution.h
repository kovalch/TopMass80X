#ifndef KinematicReconstructionSolution_h
#define KinematicReconstructionSolution_h

#include "classes.h"









struct KinematicReconstructionSolution{
    
    /// Default constructor
    KinematicReconstructionSolution();
    
    /// Constructor for setting values of the solution
    KinematicReconstructionSolution(const int leptonIndex, const int antiLeptonIndex,
                                    const int bjetIndex, const int antiBjetIndex,
                                    const LV& wPlus, const LV& wMinus,
                                    const LV& top, const LV& antiTop,
                                    const LV& neutrino, const LV& antiNeutrino,
                                    const double& reconstructedTopMass,
                                    const int numberOfBtags,
                                    const double& weight_averaged_sumSmearings_mlbBased);
    
    /// Copy constructor
    KinematicReconstructionSolution(const KinematicReconstructionSolution& kinematicReconstructionSolution);
    
    /// Destructor
    ~KinematicReconstructionSolution(){}
    
    /// Assignment operator (needed to allow std::vector<KinematicReconstructionSolution> with const data members)
    const KinematicReconstructionSolution& operator=(const KinematicReconstructionSolution& rhs){return rhs;}
    
    /// Print all stored quantities
    void print()const;
    
    LV ttbar()const{return top_ + antiTop_;}
    
    // Indices of input quantities
    const int leptonIndex_;
    const int antiLeptonIndex_;
    const int bjetIndex_;
    const int antiBjetIndex_;
    
    // Reconstructed quantities
    const LV wPlus_;
    const LV wMinus_;
    const LV top_;
    const LV antiTop_;
    const LV neutrino_;
    const LV antiNeutrino_;
    
    const double reconstructedTopMass_;
    
    // Number of b-tagged jets for solution
    const int numberOfBtags_;
    
    // FIXME: need to define enums for specific weights and use them consistently throughout the code
    // FIXME: give useful name to weight, and description
    // Weights
    const double weight_averaged_sumSmearings_mlbBased_;
};





class KinematicReconstructionSolutions{
    
public:
    
    /// Empty constructor
    KinematicReconstructionSolutions();
    
    /// Destructor
    ~KinematicReconstructionSolutions(){}
    
    
    
    // FIXME: could be made one single function, since numberOfBtags is known for solution
    void addSolutionTwoBtags(const KinematicReconstructionSolution& solution);
    void addSolutionOneBtag(const KinematicReconstructionSolution& solution);
    void addSolutionNoBtags(const KinematicReconstructionSolution& solution);
    
    
    
    /// Number of all solutions
    size_t numberOfSolutions()const{return v_solution_.size();}
    
    /// Number of 2 b-tag solutions
    size_t numberOfSolutionsTwoBtags()const{return v_solutionTwoBtags_.size();}
    
    /// Number of 1 b-tag solutions
    size_t numberOfSolutionsOneBtag()const{return v_solutionOneBtag_.size();}
    
    /// Number of 0 b-tag solutions
    size_t numberOfSolutionsNoBtags()const{return v_solutionNoBtags_.size();}
    
    
    
    /// Access from all solutions the one selected with solutionNumber, ranked by decreasing specific weight
    const KinematicReconstructionSolution& solution_averaged_sumSmearings_mlbBased(const size_t solutionNumber =0)const;
    
    /// Access from solutions with 2 b-tags the one selected with solutionNumber, ranked by decreasing specific weight
    const KinematicReconstructionSolution& solutionTwoBtags_averaged_sumSmearings_mlbBased(const size_t solutionNumber =0)const;
    
    /// Access from solutions with 1 b-tag the one selected with solutionNumber, ranked by decreasing specific weight
    const KinematicReconstructionSolution& solutionOneBtag_averaged_sumSmearings_mlbBased(const size_t solutionNumber =0)const;
    
    /// Access from solutions with 0 b-tags the one selected with solutionNumber, ranked by decreasing specific weight
    const KinematicReconstructionSolution& solutionNoBtags_averaged_sumSmearings_mlbBased(const size_t solutionNumber =0)const;
    
    
    
private:
    
    /// Insert solution index for specific weight in vector for all solutions, ranked by weight
    void insertIndex(const size_t index, const std::vector<KinematicReconstructionSolution>& v_solution,
                     const double weight, std::vector<size_t>& v_index);
    
    /// Insert solution index for specific weight in vector for b-tag categorised solutions, ranked by weight
    void insertIndex(const size_t index, const std::vector<const KinematicReconstructionSolution*>& v_solution,
                     const double weight, std::vector<size_t>& v_index);
    
    
    
    /// Vector containing all solutions
    std::vector<KinematicReconstructionSolution> v_solution_;
    
    /// Vector containing pointers to the solutions with 2 b-tags stored in v_solution_
    std::vector<const KinematicReconstructionSolution*> v_solutionTwoBtags_;
    
    /// Vector containing pointers to the solutions with 1 b-tag stored in v_solution_
    std::vector<const KinematicReconstructionSolution*> v_solutionOneBtag_;
    
    /// Vector containing pointers to the solutions with 0 b-tags stored in v_solution_
    std::vector<const KinematicReconstructionSolution*> v_solutionNoBtags_;
    
    
    
    /// Vector containing indices of all solutions, ordered for the specific weight
    std::vector<size_t> v_index_weight_averaged_sumSmearings_mlbBased_;
    
    /// Vector containing indices of 2 b-tag solutions, ordered for the specific weight
    std::vector<size_t> v_indexTwoBtags_weight_averaged_sumSmearings_mlbBased_;
    
    /// Vector containing indices of 1 b-tag solutions, ordered for the specific weight
    std::vector<size_t> v_indexOneBtag_weight_averaged_sumSmearings_mlbBased_;
    
    /// Vector containing indices of 0 b-tag solutions, ordered for the specific weight
    std::vector<size_t> v_indexNoBtags_weight_averaged_sumSmearings_mlbBased_;
};





#endif


