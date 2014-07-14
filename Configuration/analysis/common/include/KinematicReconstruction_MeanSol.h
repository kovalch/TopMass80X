#ifndef KinematicReconstruction_MeanSol_h
#define KinematicReconstruction_MeanSol_h

#include <vector>

#include <TLorentzVector.h>


using namespace std;

class KinematicReconstruction_MeanSol {

public:
    
  KinematicReconstruction_MeanSol(const double& topm);
 ~KinematicReconstruction_MeanSol();  
  void add(const TLorentzVector& top, const TLorentzVector& topbar, const TLorentzVector& n, const TLorentzVector& nbar,const double& weight,const double& mbl_weight);
  void add(const TLorentzVector& top, const TLorentzVector& topbar, const TLorentzVector& n, const TLorentzVector& nbar,const double& weight);
  
  void getMeanVect(TLorentzVector& lv, const vector<TLorentzVector>& vlv,double mass);
  void getMeanSol(TLorentzVector& top,TLorentzVector& topbar,TLorentzVector& n,TLorentzVector& nbar);
  double getSumWeight();
  int getNsol();
  void clear();
   
  private:
    vector<TLorentzVector> v_top_;
    vector<TLorentzVector> v_topbar_;
    vector<TLorentzVector> v_n_;
    vector<TLorentzVector> v_nbar_;
    
    vector<double> v_weight_;
    double sum_weight_;  
    double max_sum_weight_;

    double mass_top_;
    
};

#endif