#ifndef Output_h
#define Output_h

#include <vector>
#include <string>
#include <map>

#include <TString.h>

typedef std::vector<std::pair<TString, std::vector<double> > > mapDouble;
// typedef std::map<TString, std::vector<double> >  mapDouble;
typedef  std::map<TString, std::vector<TString> >  mapString;

class Output{
    
public:
    
    /// Constructor
    Output(const char* );
    
    /// Destructor
    //~Output(){};
    
    void addLine(TString);
    void add(const TString&, TString);
    void add(const TString&, std::vector<double>);
    void save(const TString&);
   
    
private:
    
    std::vector<TString> headerLine_;
    mapDouble mainBodyD_;
    mapString mainBodyS_;
    int nLines;
    const char* type;
    
};





#endif







