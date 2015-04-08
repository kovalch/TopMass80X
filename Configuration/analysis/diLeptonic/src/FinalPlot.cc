
#include <fstream>
#include <iostream>
// #include <cstdio>
#include <sstream>
#include <cmath>
// #include <iomanip>
#include <stdlib.h>

#include <TCanvas.h>
//#include <TLegend.h>
//#include <TExec.h>
//#include <TStyle.h>
#include <TSystem.h>
//#include <TMath.h>
//#include <TROOT.h>
//#include <THStack.h>
//#include <TDirectory.h>
//#include <TFile.h>
// #include <TString.h>
//#include <TTree.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
//#include <TH2.h>
//#include <TGraph.h>

#include "utils.h"
#include "Output.h"
#include "FinalPlot.h"
#include "Plotter.h"
//#include "ttbarUtils.h"
#include "Samples.h"
//#include "Sample.h"
#include "../../common/include/utils.h"
//#include "../../common/include/plotterUtils.h"
#include "UsefulTools.h"


//#include "TUnfold.h"
#include <TLine.h>
//#include <TMath.h>
//#include "TUnfoldBinning.h"
//#include "TUnfoldDensity.h"





FinalPlot::FinalPlot(const std::vector<Channel::Channel>& v_channel,
                     const std::vector<Systematic::Systematic>& v_systematic,
                     const double lumi, const double topxsec):

v_channel_(v_channel),
v_systematic_(v_systematic),
lumi_(lumi),
topxsec_(topxsec),
nD_(-999),
v_plotName_(std::vector<TString>(0)),
v_plotTitle_(std::vector<TString>(0)),
v_plotUnits_(std::vector<TString>(0))
{
    
}



void FinalPlot::setOptions(const std::vector<TString> v_plotName)
{
    
    //std::cout<<"[FinalPlot]:--- Beginning of plot options settings\n\n";
    
    //Clearing previous declarations  
    this->clearMemory();
    
    v_plotName_ = v_plotName;
    nD_ = (int)v_plotName_.size();
    TString nDflag;
    if(nD_==1) nDflag="1d";
    else if(nD_==2) nDflag="2d";
    else{
         std::cout<<"ERROR in FinalPlot! nD_ is not equal to 1 or 2.\n...break\n"<<std::endl;
         exit(237);
     }
    
        // reading Control Cards
    for(auto plotName : v_plotName_){
        const std::string ccFileName(common::CMSSW_BASE() + "/src/TopAnalysis/Configuration/analysis/diLeptonic/ControlCards/" + plotName + ".txt");
        std::ifstream ccFileStream(ccFileName.data(), std::ifstream::in);
        if (!ccFileStream.good()) {
          std::cerr<<"Error in FinalPlot! Cannot find file with name: "<< ccFileName <<"\n...break\n"<<std::endl;
          exit(12);
        }

       // Loop over all lines in ccFile
       bool isInBlock = false;
       while(ccFileStream.good()){
        // Read cc-File
            std::string line;
            getline(ccFileStream, line);
            line.erase(0, line.find_first_not_of(" \t"));
            if (line.size() == 0 || line[0] == '#') continue;

            // Loop over words in cc-File line and fill vWord
             std::vector<TString> vWord;
             std::string word;
             for (std::stringstream ss(line); ss >> word; ){
                vWord.push_back(word);
             }
             
             if(vWord.at(0) == "title")v_plotTitle_.push_back(vWord.at(1));
             
             if(vWord.at(0) == "units" ){
                if(vWord.size()>1) v_plotUnits_.push_back(vWord.at(1)); 
                if(vWord.size()<2) v_plotUnits_.push_back("");
             }
             
             //Pass to dimensional block
             if( (int)vWord.size()==1 && vWord.at(0) == nDflag ){
                 if(isInBlock) isInBlock = false;
                 if(!isInBlock) isInBlock = true;
            }
             
             // Reading only from dimensional block we are interesting
             if(isInBlock){
                 
             }//read from dimensional block
             else continue;
             
       }//stream line
       
    }//dimension name
    
    
    
    
    //std::cout<<"[FinalPlot]:--- Finishing of plot options settings\n\n";
}


void FinalPlot::clearMemory()
{
    v_plotName_.clear();
    v_plotTitle_.clear();
    v_plotUnits_.clear();
}


void FinalPlot::producePlots(const TString& fileType)
{
    
    bool allSyst = ((int)v_systematic_.size() > 1);
        //std::cout<<"[FinalPlot]:--- Beginning of plot production\n\n";
    
    // Loop over all channels , bins and systematics
    for(const auto& channel : v_channel_){//channel loop
        TString outPlotsFolder = "./Plots/FinalPlot/"+Channel::convert(channel)+"/"+v_plotName_.at(0)+"-"+v_plotName_.at(1)+"/";
            gSystem->mkdir(outPlotsFolder, true);
        TString inFolder = "./Plots/Nominal/"+Channel::convert(channel)+"/"+v_plotName_.at(0)+"-"+v_plotName_.at(1)+"/";
        
        std::vector<std::vector<double> > vvBinBinSystStat;
        std::vector<std::vector<TString> > vvLinksXsecCp;
        
        for(int nD=0;nD<nD_;nD++){//dimension loop
            TString firstPartName=v_plotName_.at(nD);
            TString secondPartName=v_plotName_.at(nD_-1-nD);
            TString nominalFile = inFolder + "xSec_"+firstPartName+"_IN_"+secondPartName+"_0.txt";
            std::vector<double> vBins;
                utils::readLineToVector(nominalFile,"bins" ,vBins);
            TString nbinsFile = inFolder + "xSec_"+secondPartName+"_IN_"+firstPartName+"_0.txt";
            std::vector<double> v_nBins;
                utils::readLineToVector(nbinsFile,"bins" ,v_nBins);
            int nBins = (int)(v_nBins.size())-1;
            if(nD==0)nBinsVar2_=nBins;
            
            for(int ibin=0;ibin<nBins;ibin++){//bin loop
                TString file = inFolder + "xSec_"+firstPartName+"_IN_"+secondPartName+"_" + utils::numToString(ibin) + ".txt";
                std::vector<double> vBin;
                    utils::readLineToVector(file,"bin" ,vBin);
                std::vector<double> vXsec;
                std::vector<double> vMC;
                std::vector<double> vMCpowheg;
                std::vector<double> vMCpowhegherwig;
                std::vector<double> vMCmcatnlo;
                std::vector<double> vStatErr;
                    utils::readLineToVector(file,"xsec",vXsec);
                    utils::readLineToVector(file,"stat" ,vStatErr);
                    utils::readLineToVector(file,"mc",vMC);
                std::vector<double> vSystErrUp((int)vXsec.size(),0);
                std::vector<double> vSystErrDown((int)vXsec.size(),0);
                std::vector<double> vTotalErrUp;
                std::vector<double> vTotalErrDown;
                Output systInfo("syst");
                    systInfo.add(secondPartName+"-bin",vBin);
                    systInfo.add(firstPartName+"-bins",vBins);
                    systInfo.add("Nominal-xsec",vXsec);
                    systInfo.add("mc-xsec",vMC);
                    systInfo.add("stat-err [%]",utils::divideVect(vStatErr,vXsec,100,1));
                    
                //fill bc info
                std::vector<double> vBcXsec;
                TString bcFile = file;
                    bcFile.ReplaceAll("Plots","Plots-bc");
                    utils::readLineToVector(bcFile,"xsec",vBcXsec);
                    systInfo.add("bc [%]",utils::divideVect(utils::diffVect(vBcXsec,vXsec),vXsec,100,2));
                // ...
                
                for(const auto& systematic : v_systematic_){//systematics loop
                            std::cout << "channel: " << Channel::convert(channel) << " " << systematic.name() << std::endl;
                        TString tempFile = file;
                            tempFile.ReplaceAll("Nominal",systematic.name());
                            
                        if(systematic.name() == "MCATNLO"){ utils::readLineToVector(tempFile,"mc",vMCmcatnlo); continue;}
                        if(systematic.name() == "POWHEG") {utils::readLineToVector(tempFile,"mc",vMCpowheg); continue;}
                        if(systematic.name() == "POWHEGHERWIG") {utils::readLineToVector(tempFile,"mc",vMCpowhegherwig); continue;}
                        
                        std::vector<double> vSyst;
                        std::vector<double> vDiff;
                            utils::readLineToVector(tempFile,"xsec",vSyst);
                            vDiff = utils::diffVect(vSyst,vXsec);
                        for(size_t i=0;i<vXsec.size();i++){
                            vSystErrUp.at(i) = sqrt(pow(vSystErrUp.at(i),2) + (vDiff.at(i)>=0)*pow(vDiff.at(i),2));
                            vSystErrDown.at(i) = sqrt(pow(vSystErrDown.at(i),2) + (vDiff.at(i)<0)*pow(vDiff.at(i),2));
                        }
                        systInfo.add(systematic.name()+" [%]",utils::divideVect(vDiff,vXsec,100,1));
                }//systematic loop      
                
                if(nD==1){
                    for(size_t i=0;i<vXsec.size();i++)
                    {
                        std::vector<double> vBinBinSystStat;
                        vBinBinSystStat.push_back(vBin.at(0));
                        vBinBinSystStat.push_back(vBin.at(1));
                        vBinBinSystStat.push_back(vBins.at(i));
                        vBinBinSystStat.push_back(vBins.at(i+1));
                        vBinBinSystStat.push_back(vXsec.at(i)*1000000);
                        vBinBinSystStat.push_back(100*vSystErrDown.at(i)/vXsec.at(i));
                        vBinBinSystStat.push_back(100*vSystErrUp.at(i)/vXsec.at(i));
                        vBinBinSystStat.push_back(100*vStatErr.at(i)/vXsec.at(i));
                        utils::setPrecision( vBinBinSystStat,1);
                        vvBinBinSystStat.push_back(vBinBinSystStat);
                    }
                        std::vector<TString> vLinksXsecCp;
                           TString fileLink=file;
                               fileLink.ReplaceAll("Nominal","FinalPlot");
                               fileLink.ReplaceAll(".txt",".pdf");
                          vLinksXsecCp.push_back(fileLink);
                               fileLink.ReplaceAll("FinalPlot","Nominal");
                               fileLink.ReplaceAll("xSec_","CP_");
                           vLinksXsecCp.push_back(fileLink);
                               fileLink.ReplaceAll("Nominal","FinalPlot");
                               fileLink.ReplaceAll(".pdf",".txt");
                               fileLink.ReplaceAll("CP_","systErr_");
                           vLinksXsecCp.push_back(fileLink);
                           vvLinksXsecCp.push_back(vLinksXsecCp);
                }
                
                    systInfo.save(outPlotsFolder+"systErr_"+firstPartName+"_IN_"+secondPartName+"_" + utils::numToString(ibin) + ".txt");
                    
                std::vector<double> vBinCentr;
                std::vector<double> xErr;
                TH1D histStat("histStat","Data",(int)vXsec.size(),vBins.data());
                TH1D histMC("histMC","MadGraph + Pythia",(int)vXsec.size(),vBins.data());
                TH1D histMCmcatnlo("histMCmcatnlo","MC@NLO + Herwig",(int)vXsec.size(),vBins.data());
                TH1D histMCpowheg("histMCpowheg","Powheg + Pythia",(int)vXsec.size(),vBins.data());
                TH1D histMCpowhegherwig("histMCpowhegherwig","Powheg + Herwig",(int)vXsec.size(),vBins.data());
                    
                    
                for(size_t i=0;i<vBins.size()-1;i++){
                    std::cout << vStatErr.at(i) << " " ;
                    vBinCentr.push_back((vBins.at(i+1)+vBins.at(i))*0.5);
                    xErr.push_back(0);
                    vTotalErrUp.push_back( sqrt(pow(vSystErrUp.at(i),2) + pow(vStatErr.at(i),2)) );
                    vTotalErrDown.push_back( sqrt(pow(vSystErrDown.at(i),2) + pow(vStatErr.at(i),2)) );
                    histStat.SetBinContent(i+1,vXsec.at(i));
                    histStat.SetBinError(i+1,vStatErr.at(i));
                    histMC.SetBinContent(i+1,vMC.at(i));
                    if (allSyst)
                    {
                        histMCmcatnlo.SetBinContent(i+1,vMCmcatnlo.at(i));
                        histMCpowheg.SetBinContent(i+1,vMCpowheg.at(i));
                        histMCpowhegherwig.SetBinContent(i+1,vMCpowhegherwig.at(i));
                    }
                    
                }
                std::cout << std::endl;
                
                
                TString plotTitle = utils::makeBinTitle(v_plotTitle_.at(nD_-1-nD),vBin.at(0),vBin.at(1));
                        plotTitle = plotTitle + ";" + v_plotTitle_.at(nD) + ", " + v_plotUnits_.at(nD);
                        //plotTitle = "";
                
                TGraphAsymmErrors graphTotal((int)vXsec.size(), vBinCentr.data(), vXsec.data(), xErr.data(), xErr.data(),vTotalErrDown.data(), vTotalErrUp.data());
                    
                
                    
                
                TCanvas* canvas = utils::setCanvas();
                double legX1=0.1,legY1=0.1,legDX=0.35,legDY=0.27;
                if(firstPartName == "top_pt"){
                    legX1=0.55;
                    legY1=0.55;
                }
                TLegend* legend = utils::setLegend(legX1,legY1,legX1+legDX,legY1+legDY);
                
                    histMC.SetLineColor(2);
                    histMC.SetTitle(plotTitle);
                    histMC.SetAxisRange(0,( histStat.GetMaximum() > histMC.GetMaximum() ? histStat.GetMaximum()*1.15 : histMC.GetMaximum()*1.15 ),"Y");
                    
                    histMCmcatnlo.SetLineColor(4);
                    histMCmcatnlo.SetLineStyle(4);
                    histMCmcatnlo.SetLineWidth(2);
                    histMCpowheg.SetLineColor(3);
                    histMCpowheg.SetLineStyle(2);
                    histMCpowheg.SetLineWidth(2);
                    histMCpowhegherwig.SetLineColor(1);
                    histMCpowhegherwig.SetLineStyle(4);
                    histMCpowhegherwig.SetLineWidth(2);
                    
                    histMC.SetStats(0);
                    histMC.Draw("HIST");
                    
                    if(allSyst){
                        histMCmcatnlo.DrawClone("HIST,SAME");
                        histMCpowheg.DrawClone("HIST,SAME");
                        histMCpowhegherwig.DrawClone("HIST,SAME");
                    }
                    
                    histStat.SetMarkerSize(histStat.GetMarkerSize()*0.5);
                    histStat.SetTitle(plotTitle);
                    histStat.Draw("e1, same");
                    
                    double endErrorSize =  gStyle->GetEndErrorSize();
                    gStyle->SetEndErrorSize(0);
                    graphTotal.SetMarkerSize(0);
                    graphTotal.Draw("P,SAME");
                    gStyle->SetEndErrorSize(endErrorSize);
                    
                    
                    UsefulTools::DrawCMSLabels(lumi_,2);
                    
                    legend->AddEntry(&histStat,"Data","p");
                    legend->AddEntry(&histMC,"MadGraph + Pythia","l");
                    if(allSyst)
                    {
                        legend->AddEntry(&histMCmcatnlo,histMCmcatnlo.GetTitle(),"l");
                        legend->AddEntry(&histMCpowheg,histMCpowheg.GetTitle(),"l");
                        legend->AddEntry(&histMCpowhegherwig,histMCpowhegherwig.GetTitle(),"l");
                    }
                    legend->Draw("same");
                    
                    utils::drawRatio(&histStat, &histMC, NULL,0.5,1.6,NULL,&graphTotal);
                    
                    if(allSyst){
                        histMCmcatnlo.Divide(&histMC);
                        histMCpowheg.Divide(&histMC);
                        histMCpowhegherwig.Divide(&histMC);
                        histMCmcatnlo.DrawClone("HIST,SAME");
                        histMCpowheg.DrawClone("HIST,SAME");
                        histMCpowhegherwig.DrawClone("HIST,SAME");
                    }
                    TLine line(histMC.GetXaxis()->GetXmin(),1,histMC.GetXaxis()->GetXmax(),1);
                    line.SetLineColor(kRed);
                    line.Draw("SAME");
                    
                    canvas->Print(outPlotsFolder+"xSec_"+firstPartName+"_IN_"+secondPartName+"_" + utils::numToString(ibin) + fileType);
                    canvas->Delete();
                    
            }//bin loop
        }//dimension loop
        
        printSystTableCore(v_plotTitle_.at(0)+","+v_plotUnits_.at(0), v_plotTitle_.at(1)+","+v_plotUnits_.at(1), nBinsVar2_,
                            vvBinBinSystStat,vvLinksXsecCp,"to_delete.tex");
        system("pdflatex to_delete.tex");
        system("cp to_delete.pdf "+outPlotsFolder+"/"+v_plotName_.at(0)+"-"+v_plotName_.at(1)+".pdf");
        system("rm to_delete.*");
        
    }//channel loop
    
    //std::cout<<"\n[FinalPlot]:=== Finishing of plot production\n\n";
    
}



void FinalPlot::printSystTableCore(const TString& name1, const TString& name2,const int& nBins, 
                                   const std::vector<std::vector<double> >& vv, const std::vector<std::vector<TString> >& vvLinks,
                                   const TString& saveFile)
{
    std::ofstream file_out(saveFile);
    file_out << "\\documentclass{article}\n";
    file_out << "\\usepackage{hyperref}\n";
    file_out << "\\usepackage{multirow}\n";
    file_out << "\\usepackage{rotating}\n";
    file_out << "\\usepackage{booktabs}\n";
    file_out << "\\usepackage{verbatim}\n";
    //\usepackage{graphicx}
    file_out << "\\begin{document}\n";
    file_out << "\\begin{table}[h]\n";
    file_out << "\\begin{tabular}{lllllll}\n";
    file_out << "\\hline\n";
    file_out << "$" << name1 << "$ & $" << name2 << "$ & " << "$xSec, *10^{-6}$" << " & " << "sys Down, [\\%]" << " & " << "sys Up, [\\%]" << "&" << "stat, [\\%]" << "&links" << " \\\\ \\hline\n";
    
    int pageCounter=2;
    double val1 = -999;
    for(auto v : vv)
    {
        if(val1 != v.at(0))
        {
            val1 = v.at(0);
            file_out << "\\hline\n";
            file_out << "\\multirow{" << nBins << "}{*}{" << v.at(0) << " " << v.at(1) << "}&"  << v.at(2) << " " << v.at(3)
                     << "&" << v.at(4) << "& -" << v.at(5) << "&" << v.at(6) << "&" << v.at(7) 
                     << "&" << "\\multirow{" << nBins << "}{*}{" <<"\\hyperlink{page." << pageCounter << "}{xsec}, \\hyperlink{page." << pageCounter+1 << "}{cp}, \\hyperlink{page." << pageCounter+2 << "}{err} " 
                     << "}\\\\\n";
            pageCounter+=3;
        }
        else file_out  << "&"  << v.at(2) << " " << v.at(3)
                       << "&" << v.at(4) << "& -" << v.at(5) << "&" << v.at(6) << "&" << v.at(7) 
                       << "\\\\\n";
    }
    
    file_out << "\\end{tabular}\n";
    file_out << "\\end{table}\n";
    
    for(int i=0;i<(int) vvLinks.size();i++){
    file_out << "\\newpage\n";
    file_out << "\\includegraphics[width=1.0\\textwidth]{" << vvLinks.at(i).at(0) << "}\\\\\n";
    file_out << "\\hyperlink{page.1}{go back}\n";
    
    file_out << "\\newpage\n";
    file_out << "\\includegraphics[width=1.0\\textwidth]{" << vvLinks.at(i).at(1) << "}\\\\\n";
    file_out << "\\hyperlink{page.1}{go back}\n";
    
    file_out << "\\newpage\n";
    file_out << "\\hyperlink{page.1}{go back}\n";
    file_out << "\\verbatiminput{" << vvLinks.at(i).at(2) << "}\n";
    file_out << "\\hyperlink{page.1}{go back}\n";
    
    }
    
    file_out << "\\end{document}\n";
    file_out.close();
}


