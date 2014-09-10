#ifndef HistoListReader_h
#define HistoListReader_h

#include <vector>
#include <map>

#include <TString.h>

class TH1;



struct PlotProperties{
    
public:
    
    PlotProperties();
    ~PlotProperties();
    //return a histogram with the binning and labels as defined in the properties
    TH1* getHistogram();
    //like getHistogram, but returns a clone, i.e. the caller must delete the histogram
    TH1* getClonedHistogram();
    
    TString name;
    TString specialComment;
    TString ytitle;
    TString xtitle;
    int rebin;
    bool do_dyscale;
    bool logX;
    bool logY;
    double ymin;
    double ymax;
    double xmin;
    double xmax;
    int bins;
    std::vector<double> xbinbounds;
    std::vector<double> bincenters;
    
private:
    
    void MakeHisto();
    
    TH1* histo_;
};



class HistoListReader{
    
public:
    
    HistoListReader(const char* filename);
    bool isZombie()const;
    const PlotProperties& getPlotProperties(const TString& name)const;
    //rootcint does not like decltype! :-(
    //auto begin() -> decltype(plots_.begin()) { return plots_.begin(); }
    //auto end() -> decltype(plots_.end()) { return plots_.end(); }
    std::map<TString, PlotProperties>::const_iterator begin()const{return plots_.begin();}
    std::map<TString, PlotProperties>::const_iterator end()const{return plots_.end();}
    
private:
    
    /// Read string and keep blank spaces if string starts with quotes  "... ..."
    void readString(std::stringstream &input, TString &output)const;
    
    const char* filename_;
    bool isZombie_;
    std::map<TString, PlotProperties> plots_;
};



#endif


