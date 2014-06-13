#include "HistoListReader.h"

#include <TH1.h>
#include <TString.h>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>



HistoListReader::HistoListReader(const char* filename) : 
    filename_ ( filename ),
    isZombie_ ( false )
{
    std::ifstream controlHistStream(filename, std::ifstream::in);
    if (!controlHistStream.good()) {
        isZombie_ = true;
        std::cerr << "HistoListReader: cannot read " << filename << std::endl;
        return;
    }
    plots_.clear();
    
    while(controlHistStream.good()){
        // Read HistoList-File
        std::string line;
        getline(controlHistStream, line);
        //remove leading whitespace
        line.erase(0, line.find_first_not_of(" \t"));    
        
        //skip empty lines and/or comments
        if (line.size() == 0 || line[0] == '#') continue;
        
        PlotProperties m;
        std::string dummy;
        std::ostringstream combined;
        std::stringstream linestream(line);
//        # Name, Extra, axis labels (y,x), rebin, do_dyscale, logx, logy, ymin, ymax, xmin, xmax, nbins, xbins, bcs
        linestream 
            >> m.name
            >> m.specialComment;
            readString(linestream, m.ytitle);
            readString(linestream, m.xtitle);
            linestream >> m.rebin
            >> m.do_dyscale
            >> m.logX
            >> m.logY
            >> m.ymin
            >> m.ymax
            >> m.xmin
            >> m.xmax
            >> m.bins;

        m.xbinbounds.clear();
        m.bincenters.clear();

        for(int i = 0; i <= m.bins; ++i){
            double temp;
            linestream>>temp;
            m.xbinbounds.push_back(temp);
        }
        for(int i = 0; i < m.bins; i++){//only until bincenter code is finalized
            double temp;
            linestream>>temp;
            m.bincenters.push_back(temp);
        }
        plots_[m.name] = m;
        
        if (linestream.fail()) {
            std::cerr 
            << "********************************************************************\n"
            << "Error reading file (too few entries?)\n" 
            << "File: '" << filename << "'\n"
            << "Line: '" << line << "'\n"
            << "********************************************************************\n"
            << std::endl; 
            exit(981);
        }
        
        if (!linestream.eof()) {
            std::string rest;
            getline(linestream, rest);
            std::cerr 
            << "********************************************************************\n"
            << "Too many entries!\n" 
            << "File: '" << filename << "'\n"
            << "Line: '" << line << "'\n"
            << "Not used: '" << rest << "'\n"
            << "********************************************************************\n"
            << std::endl; 
            exit(981);
        }
    }
}

void HistoListReader::readString(std::stringstream &input, TString &output)const
{
    std::string dummy;
    std::stringstream combined;

    input >> dummy;

    if(dummy.size() > 0 && dummy[0] == '"'){
        dummy.erase(0,1);
        if(dummy.size() ==0) input >> dummy;
        int nTokens =0;
        while(dummy[dummy.size()-1] != '"' && !input.eof())
        {
            if (nTokens++ > 0) combined << " ";
            combined << dummy;
            input >> dummy;
        }
        if(dummy[dummy.size()-1] == '"')
        {
            dummy.erase(dummy.size()-1);
        }
        if (nTokens++ > 0) combined << " ";
        combined << dummy;
        output = combined.str();
    } else {
        output = dummy;
    }
}

bool HistoListReader::IsZombie() const
{
    return isZombie_;
}

/////////////////////////////////////////

PlotProperties& HistoListReader::getPlotProperties(TString name)
{
    std::map <TString, PlotProperties >::iterator it = plots_.find(name);
    if (it == plots_.end()) { std::cerr << "no such plot in " << filename_ << ": ``" << name << "''" <<std::endl; exit(671); }
    return it->second;
}

PlotProperties::PlotProperties() :
    histo_ (0)
{
}

TH1* PlotProperties::getHistogram()
{
    if (!histo_) MakeHisto();
    return histo_;
}

TH1* PlotProperties::getClonedHistogram()
{
    bool old = TH1::AddDirectoryStatus();
    TH1::AddDirectory(false);
    TH1* clone = static_cast<TH1*>(getHistogram()->Clone());
    TH1::AddDirectory(old);
    return clone;
}


void PlotProperties::MakeHisto()
{
    bool old = TH1::AddDirectoryStatus();
    TH1::AddDirectory(false);
    histo_ = new TH1D(name, name, bins, &xbinbounds[0]);
    TH1::AddDirectory(old);
    histo_->GetXaxis()->SetTitle(xtitle);
    histo_->GetYaxis()->SetTitle(ytitle);
}

PlotProperties::~PlotProperties()
{
    delete histo_;
}



////                                  !!!  2d - part of HistoListReader  !!!                                      ////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

HistoListReader_2d::HistoListReader_2d(const char* filename) : 
    filename_ ( filename ),
    isZombie_ ( false )
{
    std::ifstream controlHistStream(filename, std::ifstream::in);
    if (!controlHistStream.good()) {
        isZombie_ = true;
        std::cerr << "HistoListReader: cannot read " << filename << std::endl;
        return;
    }
    plots_.clear();
    
    while(controlHistStream.good()){
        // Read HistoList-File
        std::string line;
        getline(controlHistStream, line);
        //remove leading whitespace
        line.erase(0, line.find_first_not_of(" \t"));
        
        //skip empty lines and/or comments
        if (line.size() == 0 || line[0] == '#') continue;
        
        PlotProperties_2d m;
        //std::string dummy;
        //std::ostringstream combined;
        std::stringstream linestream(line);
        std::string trash;
//        # Name, yAxisName_, nbinY, y_binsArr_, # xAxisName_, nbinX, x_binsArr_, # fileNameY, binWidthY, ymin, ymax # fileNameX, binWidthX, xmin, xmax # Rmin , Rmax
        
        m.y_binsArr.clear();
        m.x_binsArr.clear();
        
        linestream 
            >> m.name;
            
            readString(linestream, m.ytitle);
            linestream >> m.nbinY;
            for(int i = 0; i <= m.nbinY; ++i){
            Double_t temp;
            linestream>>temp;
            m.y_binsArr.push_back(temp);
            }
            linestream >> trash;
            
            readString(linestream, m.xtitle);
            linestream >> m.nbinX;
            for(int i = 0; i <= m.nbinX; ++i){
            Double_t temp;
            linestream>>temp;
            m.x_binsArr.push_back(temp);
            }
            linestream >> trash;
            
            linestream >> m.fileNameY
            >> m.binWidthY
            >> m.ymin
            >> m.ymax
            >> trash
            >> m.fileNameX
            >> m.binWidthX
            >> m.xmin
            >> m.xmax
            >> trash
            >> m.Rmin
            >> m.Rmax;
            
        plots_[m.name] = m;
        
        if (linestream.fail()) {
            std::cerr 
            << "********************************************************************\n"
            << "Error reading file (too few entries?)\n" 
            << "File: '" << filename << "'\n"
            << "Line: '" << line << "'\n"
            << "********************************************************************\n"
            << std::endl; 
            exit(981);
        }
        
        if (!linestream.eof()) {
            std::string rest;
            getline(linestream, rest);
            std::cerr 
            << "********************************************************************\n"
            << "Too many entries!\n" 
            << "File: '" << filename << "'\n"
            << "Line: '" << line << "'\n"
            << "Not used: '" << rest << "'\n"
            << "********************************************************************\n"
            << std::endl; 
            exit(981);
        }
    }
}

void HistoListReader_2d::readString(std::stringstream &input, TString &output)const
{
    std::string dummy;
    std::stringstream combined;

    input >> dummy;

    if(dummy.size() > 0 && dummy[0] == '"'){
        dummy.erase(0,1);
        if(dummy.size() ==0) input >> dummy;
        int nTokens =0;
        while(dummy[dummy.size()-1] != '"' && !input.eof())
        {
            if (nTokens++ > 0) combined << " ";
            combined << dummy;
            input >> dummy;
        }
        if(dummy[dummy.size()-1] == '"')
        {
            dummy.erase(dummy.size()-1);
        }
        if (nTokens++ > 0) combined << " ";
        combined << dummy;
        output = combined.str();
    } else {
        output = dummy;
    }
}

bool HistoListReader_2d::IsZombie() const
{
    return isZombie_;
}

/////////////////////////////////////////

PlotProperties_2d& HistoListReader_2d::getPlotProperties(TString name)
{
    std::map <TString, PlotProperties_2d >::iterator it = plots_.find(name);
    if (it == plots_.end()) { std::cerr << "no such plot in " << filename_ << ": ``" << name << "''" <<std::endl; exit(671); }
    return it->second;
}


PlotProperties_2d::PlotProperties_2d(){}

PlotProperties_2d::~PlotProperties_2d(){}








