#ifndef DoubleObject_h
#define DoubleObject_h

#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

/**
   \class   DoubleObject DoubleObject.h "TopAnalysis/TopAnalyzer/interface/DoubleObject.h"

   \brief   Interface class to analyze the relation of two object types within full framework of fwlite

   This class is an interface to a common analyer or a corresponding fwlite mimic
   in TopAnalysis/TopAnalzyer/bin. It takes a two object classes as argument. An 
   interface to a common book, fill & process function is expected by the common 
   analyzer (mimic) class. Each derived class has to provide these functions. The 
   functions, which need to be implemented include afill function for the first 
   template argument class only which needs to be present as the parameter for a
   second incout collection _srcB_ is optional. If this is not applicable to your 
   module the implementation of this function should throw an exception to prevent 
   missconfiguration. 

   The module provides containers for one and 2-dimensional histogram management 
   and corresponding fill functions. The fill functions are of type:

   fill("name", valueA, weight)         (to fill 1-dimensional histograms)
   fill("name", valueA, valueB, weight) (to fill 2-dimensional histograms)

   where _name_ is expectec to be to be the name of the histogram under which it
   has been booked for module internal histogram management (i.e. the std::string
   in the map container), _valueA_ and _valueB_ are the values of the quantities 
   to be filled (e.g. muon->pt()) and weight is the weight with which the histo-
   grams is to be filled.
*/

template <typename CollectionA, typename CollectionB> 
class DoubleObject{

 protected:
  /// histogram container
  std::map<std::string, TH1*>  hists_;
  std::map<std::string, TH2*>  hists2D_;
  std::map<std::string, float> treeVars_;
  TTree * tree;
  /// weights
  std::vector<edm::InputTag> wgts_;
  std::vector<double> weights;
    
 public:
  /// default constructor for fw lite
  explicit DoubleObject(){};
  /// default constructor for full fw
  explicit DoubleObject(const edm::ParameterSet& configFile){};
  /// default destructor
  virtual ~DoubleObject(){};
  /// write histograms to file for fwlite
  void write(TFile& file, const char* directory);
  /// book histograms or tree variables
  void bookVariable(edm::Service<TFileService>& fs, const char * variable,
		    unsigned int binsX, float lowX, float upX, unsigned int binsY, float lowY, float upY,
		    bool useTree=false);
  void bookVariable(edm::Service<TFileService>& fs, const char * variable, unsigned int binsX, float lowX, float upX, bool useTree=false);
  void bookVariable(edm::Service<TFileService>& fs, const char * variable);
  /// fill values into map for histograms or tree
  void fillValue(std::string variable, float value1, float value2, const double& weight=1.);
  void fillValue(std::string variable, float value1, const double& weight=1.);
  void initializeTrees(float value, const double& weight=1.);

  /**
     The following functions have to be implemented for any class
     derived from DoubleObject<CollectionA, CollectionB>
  **/
  /// histogramm booking for fwlite 
  virtual void book() = 0;
  /// histogramm booking for fwfull
  virtual void book(edm::Service<TFileService>& fileService) = 0;
  /// histogram filling for a single object collection. This 
  /// function is needed as the input of the second collection 
  /// can be optional. If this doea not apply to your module 
  /// throw an exception here.
  void fill2(const CollectionA& inputCollectionA, const double& runNumber, const double& luminosityBlockNumber, const double& eventNumber, const double& weight=1., std::vector<double> weights=0);
  virtual void fill(const CollectionA& inputCollectionA, const double& weight=1.) = 0;
  /// histogram filling for fwlite and for fwfull
  void fill2(const CollectionA& inputCollectionA, const CollectionB& inputCollectionB, const double& runNumber, const double& luminosityBlockNumber, const double& eventNumber, const double& weight=1., std::vector<double> weights=0);
  /// books the weight tags
  void book2(std::vector<edm::InputTag> wgts2_) {wgts_ = wgts2_;};
  virtual void fill(const CollectionA& inputCollectionA, const CollectionB& inputCollectionB, const double& weight=1.) = 0;
  /// everything which needs to be done after the event loop
  virtual void process() = 0;
};

template <typename CollectionA, typename CollectionB>
  void DoubleObject<CollectionA, CollectionB>::fill2(const CollectionA& inputCollectionA, const double& runNumber, const double& luminosityBlockNumber, const double& eventNumber, const double& weight, std::vector<double> weights)
{
  fillValue("runNumber", runNumber, weight);
  fillValue("luminosityBlockNumber", luminosityBlockNumber, weight);
  fillValue("eventNumber", eventNumber, weight);
  // if more than one weight is supposed to be stored
  for(unsigned int iWeight=0; iWeight < wgts_.size(); iWeight++){
    std::string weightName = wgts_[iWeight].label()+wgts_[iWeight].instance();
    if(wgts_.size() == weights.size()) {
      fillValue(weightName, weights[iWeight], weight);
    }
    else std::cout<< "ERROR!!! Size of weights ("<<weights.size()<<") != size of weight tags ("<<wgts_.size()<<")!!! No weights are filled in tree!" <<std::endl;
  }
  fill(inputCollectionA, weight);
}

template <typename CollectionA, typename CollectionB>
    void DoubleObject<CollectionA, CollectionB>::fill2(const CollectionA& inputCollectionA, const CollectionB& inputCollectionB, const double& runNumber, const double& luminosityBlockNumber, const double& eventNumber, const double& weight, std::vector<double> weights)
{
  fillValue("runNumber", runNumber, weight);
  fillValue("luminosityBlockNumber", luminosityBlockNumber, weight);
  fillValue("eventNumber", eventNumber, weight);
  // if more than one weight is supposed to be stored
  for(unsigned int iWeight=0; iWeight < wgts_.size(); iWeight++){
    std::string weightName = wgts_[iWeight].label()+wgts_[iWeight].instance();
    if(wgts_.size() == weights.size()) {
      fillValue(weightName, weights[iWeight], weight);
    }
    else std::cout<< "ERROR!!! Size of weights ("<<weights.size()<<") != size of weight tags ("<<wgts_.size()<<")!!! No weights are filled in tree!" <<std::endl;
  }
  fill(inputCollectionA, inputCollectionB, weight);
}

/// book histograms or tree variables
template <typename CollectionA, typename CollectionB>
  void DoubleObject<CollectionA, CollectionB>::bookVariable(edm::Service<TFileService>& fs, const char * variable,
							    unsigned int binsX, float lowX, float upX, unsigned int binsY, float lowY, float upY,
							    bool useTree)
{
  if(useTree && !binsY){
    //std::cout << "Adding *" << variable << "* to TTree" << std::endl;
    if(!tree){
      tree = fs->make<TTree>("tree","tree",0);
      treeVars_["weight"];
      tree->Branch("weight", &treeVars_["weight"], (std::string("weight") + "/F").c_str());
      // if more than one weight is supposed to be stored in tree
      for(unsigned iWeight=0; iWeight < wgts_.size(); iWeight++){
	if(iWeight==0) std::cout << wgts_.size() <<" additional event weights are stored: " <<std::endl;
	std::string weightName = wgts_[iWeight].label()+wgts_[iWeight].instance();
	std::cout << "weightName = " <<weightName <<std::endl;
	treeVars_[weightName];
	tree->Branch(weightName.c_str(), &treeVars_[weightName], (std::string(weightName) + "/F").c_str());
      }
      treeVars_["runNumber"];
      tree->Branch("runNumber", &treeVars_["runNumber"], (std::string("runNumber") + "/F").c_str());
      treeVars_["luminosityBlockNumber"];
      tree->Branch("luminosityBlockNumber", &treeVars_["luminosityBlockNumber"], (std::string("luminosityBlockNumber") + "/F").c_str());
      treeVars_["eventNumber"];
      tree->Branch("eventNumber", &treeVars_["eventNumber"], (std::string("eventNumber") + "/F").c_str());
    }
    treeVars_[variable];
    tree->Branch(variable, &treeVars_[variable], (std::string(variable) + "/F").c_str());
  }
  else{
    //std::cout << "Adding *" << variable << "* to Histograms" << std::endl;
    if     (!binsY &&  !lowY && !upY )  hists_  [variable] = fs->make<TH1F>( variable, variable, binsX, lowX, upX );
    else if( binsY && ( lowY ||  upY )) hists2D_[variable] = fs->make<TH2F>( variable, variable, binsX, lowX, upX, binsY, lowY, upY );
  }
}

template <typename CollectionA, typename CollectionB>
  void DoubleObject<CollectionA, CollectionB>::bookVariable(edm::Service<TFileService>& fs, const char * variable,
							    unsigned int binsX, float lowX, float upX, bool useTree)
{
  bookVariable( fs, variable, binsX, lowX, upX, 0, 0, 0, useTree);
}

template <typename CollectionA, typename CollectionB>
  void DoubleObject<CollectionA, CollectionB>::bookVariable(edm::Service<TFileService>& fs, const char * variable)
{
  bookVariable( fs, variable, 0, 0, 0, 0, 0, 0, true);
}

/// fill values into map for histograms or tree
template <typename CollectionA, typename CollectionB>
  void DoubleObject<CollectionA, CollectionB>::fillValue(std::string variable, float value1, float value2, const double& weight)
{
  if(hists2D_.find(variable) != hists2D_.end()) hists2D_.find(variable)->second->Fill(value1, value2, weight);
}

template <typename CollectionA, typename CollectionB>
  void DoubleObject<CollectionA, CollectionB>::fillValue(std::string variable, float value, const double& weight)
{
  if(treeVars_.find(variable) != treeVars_.end()){
    treeVars_.find(variable)->second = value;
    treeVars_.find("weight")->second = weight;
  }
  if(hists_.find(variable) != hists_.end()){
    hists_.find(variable)->second->Fill(value, weight);
  }
}

/// writing histograms to file in fwlite
template <typename CollectionA, typename CollectionB> 
void DoubleObject<CollectionA, CollectionB>::write(TFile& file, const char* directory)
{
  file.cd( directory );
  for(std::map<std::string, TH1*>::const_iterator hist = hists_.begin(); hist !=hists_.end(); ++hist){
    hist->second->Write( );
  }
  for(std::map<std::string, TH2*>::const_iterator hist = hists2D_.begin(); hist !=hists2D_.end(); ++hist){
    hist->second->Write( );
  }
}

/// initialize all branches with default value (can be called in all events)
template <typename CollectionA, typename CollectionB> 
void DoubleObject<CollectionA, CollectionB>::initializeTrees(float value, const double& weight)
{
  // loop all branches in the tree
  for(std::map<std::string, float>::iterator treeEntry=treeVars_.begin(); treeEntry!=treeVars_.end(); ++treeEntry){
    // skip initialising if treeEntry is one of the weights or run numbers etc.
    for(unsigned iWeight=0; iWeight < wgts_.size(); iWeight++) {
      if(treeEntry->first==(wgts_[iWeight].label()+wgts_[iWeight].instance())) return;
    }
    if(treeEntry->first!="weight"&&treeEntry->first!="runNumber"&&treeEntry->first!="eventNumber"&&treeEntry->first!="luminosityBlockNumber") treeEntry->second = value;
    else if (treeEntry->first=="weight") treeEntry->second = weight;
  }
}

#endif
