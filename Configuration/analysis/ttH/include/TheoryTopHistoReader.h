#ifndef TheoryTopHistoReader_h
#define TheoryTopHistoReader_h

#include <vector>
#include <map>

#include <TString.h>

class TH1;



struct BinLine{
    
public:
    
    BinLine(std::string string);
    
    bool hasBinInfo;
    float x;
    float val;
    float err;
};



class TheoryTopHistoReader{
    
public:    
    TH1* getHisto1D(const char* fileName, const char* histoName);
    
private:
    
    std::vector<BinLine> cleanedBinList(const std::vector<BinLine>& binLines)const;
    
    const char* fileName_;
    const char* histoName_;
    bool isZombie_;
};



#endif


