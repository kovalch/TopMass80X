#include <iostream>
#include <sstream>
#include "boost/program_options.hpp"
#include "boost/variant.hpp"

#include "TCanvas.h"
#include "TChain.h"
#include "TH2F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TFile.h"

#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "Helper.h"
#include "GenMatchAnalyzer.h"
#include "MVAAnalyzer.h"
#include "IdeogramAnalyzer.h"
#include "RooFitTemplateAnalyzer.h"

#include "tinyxml2.h"

namespace po = boost::program_options;
namespace xml = tinyxml2;

class Analysis {
 private:
  TString _samplePath;
  TString _fIdentifier;
  TString _fFile;
  TString _fMethod;

  TString _selection;

  int _fBins;
  double _fLumi;
    
  TChain* _fChain;
  TTree* _fTree;
  TTree* _tTree;
  TTree* _tTreeBkg;
  TH2F* _hEntries;

  TString _tempFilePath;
  TFile* _tempFile;

  TH2F* _hMass;
  TH2F* _hMassError;
  TH2F* _hJES;
  TH2F* _hJESError;
    
  TH2F* _hMassConstJES;
  TH2F* _hMassConstJESError;

  TH2F* _hFSig;
  TH2F* _hFSigError;
  TH2F* _hMassfSig;
  TH2F* _hMassfSigError;
  TH2F* _hJESfSig;
  TH2F* _hJESfSigError;

  TH2F* _bTagEff;
  TH2F* _cTagEff;
  TH2F* _lTagEff;

  /*
  std::map<TString, short>        _shortVariables;
  std::map<TString, int>          _intVariables;
  std::map<TString, float>        _floatVariables;
  std::map<TString, double>       _doubleVariables;
  std::map<TString, unsigned int> _uintVariables;

  std::map<TString, unsigned short*> _ushortArrayVariables;
  std::map<TString, short*>          _shortArrayVariables;
  std::map<TString, double*>         _doubleArrayVariables;

  std::map<TString, TClonesArray*> _TClonesArrayVariables;
  */

  typedef boost::variant< short, int, float, double, unsigned int, unsigned short*, short*, double*, TClonesArray*> variableTypes;

  //std::map<TString, variableTypes> _variables;
  std::map<TString, std::map<TString, variableTypes> > _variables;

  static xml::XMLDocument* _config;

  void CreateHistos();
  void CreateRandomSubset();
  void ReadConfigFromXMLFile();

  void SetBranchStatuses(TTree* tree);
  void SetBranchAddresses(TTree* tree);

  // return the PU weights for the different samples
  enum enumForPUWeights {kSummer11, kSummer11Plus05, kSummer11Minus05, kFall11, kFall11Plus05, kFall11Minus05, kFall10};
  double calcPUWeight_(enum enumForPUWeights sample, short nPU);
  double calcPDFWeight_(int whichPDFUncertainty, bool upVariation, double x1, int id1, float Q, double x2, int id2);
  double eventBTagProbability_(std::vector<double> &oneMinusBEffies, std::vector<double> &oneMinusBMistags, bool verbose = false);
  double calcBTagWeight_(int Njet=0, short* pdgId=0, TClonesArray* jets=0);

 public:
  Analysis(po::variables_map vm);
  Analysis(TString identifier, TString file, TString method, int bins, double lumi);
  ~Analysis();
    
  void Analyze(po::variables_map vm);

  void AddWeights(TTree* tempTree, bool isData=false);
    
  TH2F* GetH2Mass();
  TH2F* GetH2MassError();
  TH2F* GetH2JES();
  TH2F* GetH2JESError();
    
  TH2F* GetH2MassConstJES();
  TH2F* GetH2MassConstJESError();
    
  TH2F* GetH2FSig();
  TH2F* GetH2FSigError();
  TH2F* GetH2MassfSig();
  TH2F* GetH2MassfSigError();
  TH2F* GetH2JESfSig();
  TH2F* GetH2JESfSigError();
    
  TString GetIdentifier();

  /*
  template <class objectType>
  class variable_visitor : public boost::static_visitor<objectType>
  {
    public:
      objectType operator()(objectType& operand) const
      {
        return operand;
      }
      objectType* operator()(objectType* operand) const
      {
        return operand;
      }
  };
  */
};

