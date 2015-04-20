#ifndef SingleObject_h
#define SingleObject_h

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
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

/**
   \class   SingleObject SingleObject.h "TopAnalysis/TopAnalyzer/interface/SingleObject.h"

   \brief   Interface class to analyze single objects within full framework of fwlite

   Ths class is an interface to a common analyer or a corresponding fwlite mimic
   in TopAnalysis/TopAnalzyer/bin. It takes a single object class as argument. 
   An interface to a common book, fill & process function is expected by the common 
   analyzer (mimic) class. Each derived class has to provide these functions. It 
   provides containers for one and 2-dimensional histogram management and corres-
   ponding fill functions. The fill functions are of type:

   fill("name", valueA, weight)         (to fill 1-dimensional histograms)
   fill("name", valueA, valueB, weight) (to fill 2-dimensional histograms)

   where _name_ is expectec to be to be the name of the histogram under which it
   has been booked for module internal histogram management (i.e. the std::string
   in the map container), _valueA_ and _valueB_ are the values of the quantities 
   to be filled (e.g. muon->pt()) and weight is the weight with which the histo-
   grams is to be filled.
*/

template <typename Collection> 
class SingleObject{

 protected:
  /// histogram container
  std::map<std::string, TH1*> hists_;
  std::map<std::string, TH2*> hists2D_;
  std::map<std::string, float> treeVars_;
  std::map<std::string, unsigned int> treeVarsUI_;
  std::map<std::string, int> treeVarsI_;
  std::map<std::string, double> treeVarsD_;
  TTree * tree;
  /// weights
  std::vector<edm::InputTag> wgts_;
  std::vector<double> weights;

 public:
  /// default constructor for fwlite
  explicit SingleObject(){};
  /// default constructor for fwfull
  explicit SingleObject(const edm::ParameterSet& configFile){};
  /// default destructor
  virtual ~SingleObject(){};
  /// write histograms to file for fwlite
  void write(TFile& file, const char* directory);
  /// book histograms or tree variables
  void bookVariable(edm::Service<TFileService>& fs, const char * variable,
		    unsigned int binsX, float lowX, float upX, unsigned int binsY, float lowY, float upY,
		    bool useTree=false);
  void bookVariable(edm::Service<TFileService>& fs, const char * variable, unsigned int binsX, float lowX, float upX, bool useTree=false);
  void bookVariable(edm::Service<TFileService>& fs, const char * variable);
  /// fill values into map for histograms or tree
  void fillValue(std::string variable, float  value1, float value2, const double& weight=1.);
  void fillValue(std::string variable, float  value1, const double& weight=1.);
  void fillValue(std::string variable, int    value1, const double& weight=1.);
  void fillValue(std::string variable, unsigned int    value1, const double& weight=1.);
  void fillValue(std::string variable, double value1, const double& weight=1.);
  void initializeTrees(float value, const double& weight=1.);

  /**
     The following functions have to be implemented for any class
     derived from SingleObject<Collection>
  **/
  /// histogramm booking for fwlite 
  virtual void book() = 0;
  /// histogramm booking for fwfull
  virtual void book(edm::Service<TFileService>& fileService) = 0;
  /// histogram filling for fwlite and for fwfull
  virtual void fill(const Collection& inputCollection, const double& weight=1.) = 0;
  void fill2(const Collection& inputCollection, const double& runNumber, const double& luminosityBlockNumber, const double& eventNumber, const double& x1=0, const double& x2=0, const float& Q=0, const int& id1=-42, const int& id2=-42, const double& weight=1., std::vector<double> weights=0) ;
  /// books the weight tags
  void book2(std::vector<edm::InputTag> wgts2_) {wgts_ = wgts2_; };
  /// everything which needs to be done after the event loop
  virtual void process() = 0;

 protected:
  /// check if histogram was booked in the corresponding map
  bool booked(std::map<std::string, TH1*> map, const std::string histName) const { return map.find(histName.c_str())!=map.end(); };
  /// fill 1-dimensional histogram if it had been booked before
  void fill(const std::string histName, double valueA, double weight) const { if(booked(hists_, histName)) hists_.find(histName)->second->Fill(valueA, weight); };
  /// fill 2-dimensional histogram if it had been booked before
  void fill(const std::string histName, double valueA, double valueB, double weight) const { if(booked(hists2D_, histName)) hists2D_.find(histName)->second->Fill(valueA, valueB, weight); };

};

template <typename Collection>
    void SingleObject<Collection>::fill2(const Collection& inputCollection, const double& runNumber, const double& luminosityBlockNumber, const double& eventNumber, const double& x1, const double& x2, const float& Q, const int& id1, const int& id2, const double& weight, std::vector<double> weights)
{
  fillValue("runNumber", runNumber, weight);
  fillValue("luminosityBlockNumber", luminosityBlockNumber, weight);
  fillValue("eventNumber", eventNumber, weight);
  fillValue("id1", id1, weight);
  fillValue("id2", id2, weight);
  fillValue("x1" , x1 , weight);
  fillValue("x2" , x2 , weight);
  fillValue("Q"  , Q  , weight);
  // if more than one weight is supposed to be stored
  for(unsigned int iWeight=0; iWeight < wgts_.size(); iWeight++){
    std::string weightName = wgts_[iWeight].label()+wgts_[iWeight].instance();
    if(wgts_.size() == weights.size()) {
      fillValue(weightName, weights[iWeight], weight);
    }
    else std::cout<< "ERROR!!! Size of weights ("<<weights.size()<<") != size of weight tags ("<<wgts_.size()<<")!!! No weights are filled in tree!" <<std::endl;
  }
  fill(inputCollection, weight);
}

/// book histograms or tree variables
template <typename Collection>
  void SingleObject<Collection>::bookVariable(edm::Service<TFileService>& fs, const char * variable,
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
	tree->Branch(weightName.c_str(), &treeVars_[weightName], (weightName + "/F").c_str());
      }
      treeVars_["runNumber"];
      tree->Branch("runNumber", &treeVars_["runNumber"], (std::string("runNumber") + "/F").c_str());
      treeVars_["luminosityBlockNumber"];
      tree->Branch("luminosityBlockNumber", &treeVars_["luminosityBlockNumber"], (std::string("luminosityBlockNumber") + "/F").c_str());
      treeVars_["eventNumber"];
      tree->Branch("eventNumber", &treeVars_["eventNumber"], (std::string("eventNumber") + "/F").c_str());
      treeVars_  ["Q"   ];
      tree->Branch("Q"  , &treeVars_ ["Q"    ], (std::string("Q"  ) + "/F").c_str());
      treeVarsI_ ["id1" ];
      tree->Branch("id1", &treeVarsI_["id1"  ], (std::string("id1") + "/I").c_str());
      treeVarsI_ ["id2" ];
      tree->Branch("id2", &treeVarsI_["id2"  ], (std::string("id2") + "/I").c_str());
      treeVarsD_ ["x1"  ];
      tree->Branch("x1" , &treeVarsD_["x1"   ], (std::string("x1" ) + "/D").c_str());
      treeVarsD_ ["x2"  ];
      tree->Branch("x2" , &treeVarsD_["x2"   ], (std::string("x2" ) + "/D").c_str());
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

template <typename Collection>
  void SingleObject<Collection>::bookVariable(edm::Service<TFileService>& fs, const char * variable,
					      unsigned int binsX, float lowX, float upX, bool useTree)
{
  bookVariable( fs, variable, binsX, lowX, upX, 0, 0, 0, useTree);
}

template <typename Collection>
  void SingleObject<Collection>::bookVariable(edm::Service<TFileService>& fs, const char * variable)
{
  bookVariable( fs, variable, 0, 0, 0, 0, 0, 0, true);
}

/// fill values into map for histograms or tree
template <typename Collection>
  void SingleObject<Collection>::fillValue(std::string variable, float value1, float value2, const double& weight)
{
  if(hists2D_.find(variable) != hists2D_.end()) hists2D_.find(variable)->second->Fill(value1, value2, weight);
}

template <typename Collection>
  void SingleObject<Collection>::fillValue(std::string variable, int value, const double& weight)
{
  if(treeVarsI_.find(variable) != treeVarsI_.end()){
    treeVarsI_.find(variable)->second = value;
    treeVars_.find("weight")->second = weight;
  }
  else if(treeVars_.find(variable) != treeVars_.end()){
    treeVars_.find(variable)->second = value;
    treeVars_.find("weight")->second = weight;
  }
  if(hists_.find(variable) != hists_.end()){
    hists_.find(variable)->second->Fill(value, weight);
  }
}

template <typename Collection>
  void SingleObject<Collection>::fillValue(std::string variable, unsigned int value, const double& weight)
{
  if(treeVarsUI_.find(variable) != treeVarsUI_.end()){
    treeVarsUI_.find(variable)->second = value;
    treeVars_.find("weight")->second = weight;
  }
  else if(treeVars_.find(variable) != treeVars_.end()){
    treeVars_.find(variable)->second = value;
    treeVars_.find("weight")->second = weight;
  }
  if(hists_.find(variable) != hists_.end()){
    hists_.find(variable)->second->Fill(value, weight);
  }
}

template <typename Collection>
  void SingleObject<Collection>::fillValue(std::string variable, float value, const double& weight)
{
  if(treeVars_.find(variable) != treeVars_.end()){
    treeVars_.find(variable)->second = value;
    treeVars_.find("weight")->second = weight;
  }
  else if(treeVars_.find(variable) != treeVars_.end()){
    treeVars_.find(variable)->second = value;
    treeVars_.find("weight")->second = weight;
  }
  if(hists_.find(variable) != hists_.end()){
    hists_.find(variable)->second->Fill(value, weight);
  }
}

template <typename Collection>
  void SingleObject<Collection>::fillValue(std::string variable, double value, const double& weight)
{
  if(treeVarsD_.find(variable) != treeVarsD_.end()){
    treeVarsD_.find(variable)->second = value;
    treeVars_.find("weight")->second = weight;
  }
  else if(treeVars_.find(variable) != treeVars_.end()){
    treeVars_.find(variable)->second = value;
    treeVars_.find("weight")->second = weight;
  }
  if(hists_.find(variable) != hists_.end()){
    hists_.find(variable)->second->Fill(value, weight);
  }
}

/// writing histograms to file in fwlite
template <typename Collection> 
void SingleObject<Collection>::write(TFile& file, const char* directory)
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
template <typename Collection>
void SingleObject<Collection>::initializeTrees(float value, const double& weight)
{
  // loop all branches in the tree
  for(std::map<std::string, float>::iterator treeEntry=treeVars_.begin(); treeEntry!=treeVars_.end(); ++treeEntry){
    // skip initialising if treeEntry is one of the weights or run numbers etc.
    for(unsigned iWeight=0; iWeight < wgts_.size(); iWeight++) {
      if(treeEntry->first==(wgts_[iWeight].label()+wgts_[iWeight].instance())) return;
    }
    if(treeEntry->first!="weight"&&treeEntry->first!="runNumber"&&treeEntry->first!="eventNumber"&&treeEntry->first!="luminosityBlockNumber"&&treeEntry->first!="x1"&&treeEntry->first!="x2"&&treeEntry->first!="Q"&&treeEntry->first!="id1"&&treeEntry->first!="id2") treeEntry->second = value;
    else if (treeEntry->first=="weight") treeEntry->second = weight;
  }
  //for(std::map<std::string, int>::iterator treeEntry=treeVarsI_.begin(); treeEntry!=treeVarsI_.end(); ++treeEntry){
  //  if(treeEntry->first!="weight"&&treeEntry->first!="runNumber"&&treeEntry->first!="eventNumber"&&treeEntry->first!="luminosityBlockNumber"&&treeEntry->first!="x1"&&treeEntry->first!="x2"&&treeEntry->first!="Q"&&treeEntry->first!="id1"&&treeEntry->first!="id2") treeEntry->second = -42;
  //}
  //for(std::map<std::string, unsigned int>::iterator treeEntry=treeVarsUI_.begin(); treeEntry!=treeVarsUI_.end(); ++treeEntry){
  //  if(treeEntry->first!="weight"&&treeEntry->first!="runNumber"&&treeEntry->first!="eventNumber"&&treeEntry->first!="luminosityBlockNumber"&&treeEntry->first!="x1"&&treeEntry->first!="x2"&&treeEntry->first!="Q"&&treeEntry->first!="id1"&&treeEntry->first!="id2") treeEntry->second =  42;
  //}
  //for(std::map<std::string, double>::iterator treeEntry=treeVarsD_.begin(); treeEntry!=treeVarsD_.end(); ++treeEntry){
  //  if(treeEntry->first!="weight"&&treeEntry->first!="runNumber"&&treeEntry->first!="eventNumber"&&treeEntry->first!="luminosityBlockNumber"&&treeEntry->first!="x1"&&treeEntry->first!="x2"&&treeEntry->first!="Q"&&treeEntry->first!="id1"&&treeEntry->first!="id2") treeEntry->second = -42.0;
  //}
}

#endif
