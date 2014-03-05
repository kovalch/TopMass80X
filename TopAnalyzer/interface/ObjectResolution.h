#ifndef ObjectResolution_h
#define ObjectResolution_h

#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

/**
   \class   ObjectResolution ObjectResolution.h "TopAnalysis/TopAnalyzer/interface/ObjectResolution.h"

   \brief   Interface class to determinethe resolution of a single object within full framework of fwlite

   The class is an interface to a common analyer or a corresponding fwlite mimic
   in TopAnalysis/TopAnalzyer/bin in analogy to the SingleObject<Collection> class 
   to determine object resolutions based on generator level information. It takes 
   a single object class as argument. An interface to a common book, fill & process 
   function is expected by the common analyzer (mimic) class. Each derived class 
   has to provide these. It needs a dedicated contructor for fwlite to obtain the 
   information of histogram binnings in pt, eta and phi of the top quark.
*/

template <typename Collection> 
class ObjectResolution{

 public:
  /// default constructor for fw lite
  explicit ObjectResolution(double matchDR, std::vector<double> binsPt, std::vector<double> binsEta, std::vector<double> binsPhi);
  /// default constructor for fwfull
  explicit ObjectResolution(const edm::ParameterSet& configFile);
  /// default destructor
  virtual ~ObjectResolution(){};
  /// histogramm booking for fwlite 
  virtual void book();
  /// histogramm booking for fwfull
  virtual void book(edm::Service<TFileService>& fileService);
  /// dummy function
  void book2(std::vector<edm::InputTag> wgts2_){};
  /// everything which needs to be done after the event loop
  virtual void process(){};
  /// write histograms to file for fwlite
  void write(TFile& file, const char* directory);

  /**
     The following functions have to be implemented for any class
     derived from ObjectResolution<Collection>
  **/
  /// histogram filling for fwlite and for fwfull
  virtual void fill(const Collection& inputCollection, const double& weight=1.) = 0;
  void fill2(const Collection& inputCollection, const double& a, const double& b, const double& c, const double& x1, const double& x2, const float& Q, const int& id1, const int& id2, const double& weight, std::vector<double> weights){
    //fillValue("id1", id1, weight);
    //fillValue("id2", id2, weight);
    //fillValue("x1" , x1 , weight);
    //fillValue("x2" , x2 , weight);
    //fillValue("Q"  , Q  , weight);
    return fill(inputCollection, weight);
};

 protected:
  /**
     Helper function for special histogram management
  **/
  /// book a vector of histograms according to a given binning
  void book(std::vector<TH1F*>& hists, const char* kin, std::vector<double> binning, unsigned int nBins, double boundary);
  /// book a vector of histograms according to a given binning
  void book(edm::Service<TFileService>& fs, std::vector<TH1F*>& hists, const char* kin, std::vector<double> binning, unsigned int nBins, double boundary);

 protected:
  /// maximally acceptable distance in deltaR (for matching)
  double matchDR_;
  /// binning of resolution histograms
  std::vector<double> binsPt_;
  std::vector<double> binsEta_;
  std::vector<double> binsPhi_;

  /// resolution histograms
  std::vector<TH1F*> pt_;
  std::vector<TH1F*> eta_;
  std::vector<TH1F*> phi_;
};

/// default constructor for fw lite
template <typename Collection> 
ObjectResolution<Collection>::ObjectResolution(double matchDR, std::vector<double> binsPt, std::vector<double> binsEta, std::vector<double> binsPhi):
matchDR_( matchDR ), binsPt_ ( binsPt  ), binsEta_( binsEta ), binsPhi_( binsPhi )
{
}

/// default constructor for full fw
template <typename Collection> 
ObjectResolution<Collection>::ObjectResolution(const edm::ParameterSet& cfg):
  matchDR_( cfg.getParameter<double>( "matchDR" ) ),
  binsPt_ ( cfg.getParameter<std::vector<double> >( "binsPt" ) ),
  binsEta_( cfg.getParameter<std::vector<double> >( "binsEta") ),
  binsPhi_( cfg.getParameter<std::vector<double> >( "binsPhi") )
{
}

/// book a vector of histograms according to a given binning
template <typename Collection> 
void ObjectResolution<Collection>::book(edm::Service<TFileService>& fs, std::vector<TH1F*>& hists, const char* kin, std::vector<double> binning, unsigned int nBins, double boundary)
{
  NameScheme res("res");
  for(unsigned int idx=0; idx<binning.size()-1; ++idx){
    hists.push_back( fs->make<TH1F>(res(kin, binning[idx], binning[idx+1]), res.name(kin,binning[idx], binning[idx+1]), nBins, -boundary, boundary) );
  }  
}

/// book a vector of histograms according to a given binning
template <typename Collection> 
void ObjectResolution<Collection>::book(std::vector<TH1F*>& hists, const char* kin, std::vector<double> binning, unsigned int nBins, double boundary)
{
  NameScheme res("res");
  for(unsigned int idx=0; idx<binning.size()-1; ++idx){
    hists.push_back( new TH1F(res(kin, binning[idx], binning[idx]+1), res.name(kin,binning[idx], binning[idx+1]), nBins, -boundary, boundary) );
  }  
}

/// histogramm booking for fwlite 
template <typename Collection> 
void ObjectResolution<Collection>::book()
{
  /** 
      Resolution Histograms (with a default binning)
  **/
  book( pt_  ,  "pt" ,  binsPt_  ,  100 ,    0.5 );
  book( eta_ , "eta" ,  binsEta_ ,  100 ,  0.005 );
  book( phi_ , "phi" ,  binsPhi_ ,  100 ,  0.005 );
}

/// histogramm booking for fwlite 
template <typename Collection> 
void ObjectResolution<Collection>::book(edm::Service<TFileService>& fs)
{
  /** 
      Resolution Histograms (with a default binning)
  **/
  book( fs , pt_  ,  "pt"  , binsPt_  ,  100 ,    0.5 );
  book( fs , eta_ ,  "eta" , binsEta_ ,  100 ,  0.005 );
  book( fs , phi_ ,  "phi" , binsPhi_ ,  100 ,  0.005 );
}

/// writing histograms to file in fwlite
template <typename Collection> 
void ObjectResolution<Collection>::write(TFile& file, const char* directory)
{
  file.cd( directory );
  // write resolution histograms binned in pt to the file
  for(std::vector<TH1F*>::const_iterator hist=pt_ .begin(); hist!=pt_ .end(); ++hist){
    const TH1F* out = *hist;
    out->Write();
  }
  // write resolution histograms binned in eta to the file
  for(std::vector<TH1F*>::const_iterator hist=eta_.begin(); hist!=eta_.end(); ++hist){
     const TH1F* out = *hist;
    out->Write();
  }
  // write resolution histograms binned in phi to the file
  for(std::vector<TH1F*>::const_iterator hist=phi_.begin(); hist!=phi_.end(); ++hist){
     const TH1F* out = *hist;
    out->Write();
  }
}

#endif
