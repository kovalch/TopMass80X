#ifndef homelessFunctions_h
#define homelessFunctions_h

#include <vector>
#include <set>

class TH1;
class TH2;
class TString;
class RootFileReader;


class homelessFunctions{

public:
    
    /// Constructor for producing Drell-Yan scale factors
    homelessFunctions(RootFileReader*  rootFileReader,bool doClosureTest,bool doDYScale);
    
    /// Default destructor
    ~homelessFunctions(){}
    
    void setLumi(double newLumi, double xSec);
    
    //IVAN's Scaling Code
    double SampleXSection(const TString& filename);
    double CalcLumiWeight(const TString& WhichSample);
    std::vector<TString> InputFileList(TString mode, TString Systematic);
    void ApplyFlatWeights(TH1* varhists,   const double weight);
    void ApplyFlatWeights(TH2* varhists,   const double weight);
    
    void DYScaleFactor(TString SpetialComment,std::vector<double>& DYScale,TString name);
    
    
   static void fillSetListOfSystematics(std::set<TString>& SetOfValidSystematics);
   static void fillVectorOfValidSystematics(std::vector<const char*>& vect);
    
    std::set<TString> ListOfSyst;
    double lumi, topxsec;
    
private:
    
    /// File reader for accessing specific histogram from given file
    RootFileReader* fileReader;
    
    bool doClosureTest ;
    bool doDYScale;
    TString name;
    
};







#endif





