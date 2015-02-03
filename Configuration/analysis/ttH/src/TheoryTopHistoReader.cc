#include "TheoryTopHistoReader.h"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include <TH1.h>
#include <TString.h>





TH1* TheoryTopHistoReader::getHisto1D(const char* fileName, const char* histoName)
{
    TH1* histo = 0;
    
    std::ifstream controlHistStream(fileName, std::ifstream::in);
    if (!controlHistStream.good()) {
        std::cerr << "TheoryTopHistoReader: cannot read " << fileName << std::endl;
        return histo;
    }
    
    bool histoName_found = false;
    std::vector<BinLine> binLines;
    
    // Reading each line and storing each bin line to the list
    while(controlHistStream.good()){
        // Read HistoList-File
        std::string line;
        getline(controlHistStream, line);
        //remove leading whitespace
        line.erase(0, line.find_first_not_of(" \t"));    
        
        //skip empty lines and/or comments
        if (line.size() == 0 || line[0] == '#') continue;
        
        
        // Checking if the line contains name of the histogram
        if(line.find(histoName) != std::string::npos) {
            // Ensuring that no other part of the name follows after at least 2 whitespaces
            std::string lineEnd(line);
            lineEnd.replace(0,line.find(histoName)+std::string(histoName).length(), "");
            if(lineEnd.length() < 1 || lineEnd.find_first_not_of(" \t") - lineEnd.length() > 2) histoName_found = true;
        }
        
        // Doing nothing unless the block with the proper histogram is found
        if(!histoName_found) continue;
        
        
        BinLine binLine(line);
        // Skipping the lines that do not contain bin information
        if(binLines.size() < 1 && !binLine.hasBinInfo) continue;
        // Stopping if the last line with bin information has been analysed
        if(binLines.size() > 2 && !binLine.hasBinInfo) break;
        
        binLines.push_back(binLine);
    }
    
    // Cleaning up the list of bin lines: in case bin centers and bin boundaries are printed separately
    std::vector<BinLine> binLines_clean = cleanedBinList(binLines);
    // Extracting the list of bins for the histogram
    std::vector<double> xBins;
    for(BinLine binLine : binLines_clean) {
        xBins.push_back(binLine.x);
    }
    histo = new TH1D(TString(histoName).ReplaceAll(" ","_"), histoName, xBins.size()-1,&xBins[0]);
    
    // Setting content and error for each bin
    for(size_t binId = 0; binId < xBins.size()-1; ++binId) {
        histo->SetBinContent(binId+1, binLines_clean.at(binId).val);
        histo->SetBinError(binId+1, binLines_clean.at(binId).err);
    }
    
    return histo;
}


std::vector<BinLine> TheoryTopHistoReader::cleanedBinList(const std::vector<BinLine>& binLines)const
{
    std::vector<BinLine> binLines_clean;
    
    for(size_t lineId = 0; lineId < binLines.size(); ++lineId) {
        bool isLastLine = lineId == binLines.size()-1;
        BinLine binLine = binLines.at(lineId);
        // Filling only lines that have non-zero error
        if(binLines.at(lineId).err == 0. && !isLastLine) continue;
        // Updating the left bin boundary if the current one points to the bin center
        if(lineId > 0 && binLines.at(lineId-1).err ==0. && !isLastLine) binLine.x = binLines.at(lineId-1).x;
        
        binLines_clean.push_back(binLine);
    }
    
    return binLines_clean;
}



BinLine::BinLine(std::string string):
hasBinInfo(false),
x(-1.),
val(-1.),
err(-1.)
{
    const int nMatched = sscanf(string.c_str(), "%f %E %E", &x, &val, &err);
    hasBinInfo = nMatched == 3;
}

