#ifndef UsefulTools_h
#define UsefulTools_h

#include <vector>
#include <set>

class TH1;
class TH2;
class TString;
class RootFileReader;


class UsefulTools{

public:
    
    static constexpr double lumi_ = 19712.; // data luminosity in pb-1
    
    /// Constructor for producing Drell-Yan scale factors
    UsefulTools(RootFileReader*  rootFileReader,bool doClosureTest,bool doDYScale);
    
    /// Default destructor
    ~UsefulTools(){} 
    
    //IVAN's Scaling Code
    double SampleXSection(const TString& filename);
    double CalcLumiWeight(const TString& WhichSample);
    std::vector<TString> InputFileList(TString mode, TString Systematic);
    void ApplyFlatWeights(TH1* varhists,   const double weight);
    void ApplyFlatWeights(TH2* varhists,   const double weight);
    
    void DYScaleFactor(TString SpetialComment,std::vector<double>& DYScale,TString name);
    
    // Draw official labels (CMS Preliminary, luminosity and CM energy) above plot
    static void DrawDecayChLabel(TString decaychannel="", double textSize=0.04);
    static void DrawCMSLabels(double lumi,int cmsprelim, double energy =8 , double textSize = 0.04);
    void setStyle(TH1 *hist, TString Axis);
    
   static void fillSetListOfSystematics(std::set<TString>& SetOfValidSystematics);
   static void fillVectorOfValidSystematics(std::vector<const char*>& vect);
   static void fillLegendColorDataset(const TString& fileListName, std::vector<TString>& legends, std::vector<int>& colors, std::vector<TString>& dataset);
    
    std::set<TString> ListOfSyst;
    double topxsec;
    
private:
    
    /// File reader for accessing specific histogram from given file
    RootFileReader* fileReader;
    
    bool doClosureTest ;
    bool doDYScale;
    TString name;
    
};







#endif





