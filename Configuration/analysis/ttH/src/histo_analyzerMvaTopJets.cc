#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <iomanip>

#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TPaveStats.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TError.h"

#include "MvaFactory.h"
#include "MvaTreeHandler.h"
#include "MvaTreeAnalyzer.h"
#include "MvaWeights2d.h"
#include "mvaSetup.h"
#include "higgsUtils.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/CommandLineParameters.h"
#include "../../common/include/utils.h"
#include "../../common/include/RootFileReader.h"



/// The input base folder
constexpr const char* InputBaseDIR = "selectionRoot";

/// The output base folder
constexpr const char* OutputBaseDIR = "Plots_mvaValidation";






// FIXME: shift everything to own file
void plotWeight2dHistos(const std::vector<std::pair<TString, TString> >& v_nameStepPair,
                        const std::vector<std::pair<TString, std::vector<TString> > >& v_input,
                        const TString& outputFolder)
{
    gROOT->SetStyle("Plain");
    gROOT->ForceStyle();
    gStyle->SetOptStat(111110);
    gErrorIgnoreLevel = 1001;
    
    gStyle->SetPalette(1);      //Spektralpalette, Default: 0 resp. 50
    //gStyle->SetNumberContours(20);  // Default: 20

    double width = 600.;

    gStyle->SetCanvasDefW(width);
    gStyle->SetCanvasDefH(width);

    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.10);

    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.15);

    gStyle->SetTitleOffset(1.4,"Y");
    gStyle->SetTitleOffset(1.2,"X");


    TGaxis::SetMaxDigits(3);


    //TH1::StatOverflows(kTRUE);// compute mean etc including overflow!
    //gStyle->SetHistMinimumZero(kTRUE); // no zero-suppression on y-axis
    //gStyle->SetOptFit(222);         // 1: Fit-Parameter werden angezeigt
    //gStyle->SetCanvasDefX(400);     // canvas (default) upper edge in X
    //gStyle->SetCanvasDefY(200);     // canvas (default) upper edge in Y

    //gStyle->SetHistLineWidth(2);

    gStyle->SetTitleX(0.2);         // move upper left corner of title box to specified value
    //gStyle->SetTitleY(0.99);        // move upper left corner of title box to specified value

    gStyle->SetTitleXSize(0.05);
    gStyle->SetTitleYSize(0.05);
    gStyle->SetTitleSize(0.05,"XY");
    gStyle->SetLabelSize(0.05,"XY");
    
    
    RootFileReader* rootFileReader(RootFileReader::getInstance());
    
    for(const auto& input : v_input){
        const TString& fileShort = input.first;
        const std::vector<TString>& v_inputFilename = input.second;
        
        
        TCanvas* canvas(0);
        
        for(const auto& nameStepPair : v_nameStepPair){
            const TString& name = nameStepPair.first;
            
            TString histName = name;
            histName.ReplaceAll("mvaA_", "");
            
            // Set up canvas
            canvas = new TCanvas("canvas1");
            canvas->cd();
            
            TH2* hist(0);
            for(const auto& filename : v_inputFilename){
                if(!hist) hist = rootFileReader->GetClone<TH2D>(filename, name);
                else hist->Add(rootFileReader->GetClone<TH2D>(filename, name));
            }
            
            // Draw the histogram
            hist->SetTitle("");
            hist->Draw("box");
            canvas->Modified();
            canvas->Update();
            
            const TString plotName(outputFolder + fileShort + "/" + histName);
            
            canvas->Print(plotName + ".eps");
            canvas->Print(plotName + ".png");
            
            hist->Delete();
            canvas->Close();
        }
    }
}



void plotWeightHistos(const std::vector<std::pair<TString, TString> >& v_nameStepPairCorrect,
                      const std::vector<std::pair<TString, TString> >& v_nameStepPairSwapped,
                      const std::vector<std::pair<TString, TString> >& v_nameStepPairWrong,
                      const std::vector<std::pair<TString, std::vector<TString> > >& v_input,
                      const TString& outputFolder)
{
    gROOT->SetStyle("Plain");
    gROOT->ForceStyle();
    gStyle->SetOptStat(111110);
    gErrorIgnoreLevel = 1001;
    
    gStyle->SetPalette(1);      //Spektralpalette, Default: 0 resp. 50
    //gStyle->SetNumberContours(20);  // Default: 20

    const double width = 600.;

    gStyle->SetCanvasDefW(width);
    gStyle->SetCanvasDefH(width);

    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.10);

    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.15);

    gStyle->SetTitleOffset(1.4,"Y");
    gStyle->SetTitleOffset(1.2,"X");

    TGaxis::SetMaxDigits(3);
    
    
    RootFileReader* rootFileReader(RootFileReader::getInstance());
    
    for(const auto& input : v_input){
        const TString& fileShort = input.first;
        const std::vector<TString>& v_inputFilename = input.second;
        
        TCanvas* canvas(0);
        TLegend* legend(0);
        
        for(size_t iHist = 0; iHist < v_nameStepPairCorrect.size(); ++iHist){
            const TString& nameCorrect = v_nameStepPairCorrect.at(iHist).first;
            const TString& nameSwapped = v_nameStepPairSwapped.at(iHist).first;
            const TString& nameWrong = v_nameStepPairWrong.at(iHist).first;
            
            // 2D-histos need to be excluded here
            if(nameCorrect.Contains("_mvaWeightCorrectVsSwapped_")) continue;
            
            TString histName = nameCorrect;
            histName.ReplaceAll("_correct_", "_");
            histName.ReplaceAll("mvaA_", "");
            
            // Set up canvas
            canvas = new TCanvas("canvas1");
            canvas->cd();
            
            // Set up legend
            legend = new TLegend(0.41,0.84,0.76,0.99);
            legend->SetFillColor(0);
            legend->SetFillStyle(0);
            legend->SetTextSize(0.04);
            legend->SetMargin(0.12);
            legend->SetFillStyle(1001);
            
            // Set up histograms
            TH1* histCorrect(0);
            TH1* histSwapped(0);
            TH1* histWrong(0);
            for(const auto& filename : v_inputFilename){
                if(!histCorrect){
                    histCorrect = rootFileReader->GetClone<TH1D>(filename, nameCorrect);
                    histCorrect->SetLineColor(4);
                    legend->AddEntry(histCorrect, "correct", "l");
                }
                else histCorrect->Add(rootFileReader->GetClone<TH1D>(filename, nameCorrect));
                
                if(!histSwapped){
                    histSwapped = rootFileReader->GetClone<TH1D>(filename, nameSwapped);
                    histSwapped->SetLineColor(3);
                    legend->AddEntry(histSwapped, "swapped", "l");
                }
                else histSwapped->Add(rootFileReader->GetClone<TH1D>(filename, nameSwapped));
                
                if(!histWrong){
                    histWrong = rootFileReader->GetClone<TH1D>(filename, nameWrong);
                    histWrong->SetLineColor(2);
                    legend->AddEntry(histWrong, "wrong", "l");
                }
                else histWrong->Add(rootFileReader->GetClone<TH1D>(filename, nameWrong));
            }
            std::vector<TH1*> v_hist = {histCorrect, histSwapped, histWrong};
            
            // Scale histograms and get minimum and maximum value
            Double_t yMin(99999.);
            Double_t yMax(-99999.);
            for(std::vector<TH1*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
                TH1* hist = *i_hist;
                //hist->Scale(1./hist->Integral(0, hist->GetNbinsX()+1));
                
                const Double_t yMinHist = hist->GetBinContent(hist->GetMinimumBin());
                const Double_t yMaxHist = hist->GetBinContent(hist->GetMaximumBin());
                if(yMin > yMinHist) yMin = yMinHist;
                if(yMax < yMaxHist) yMax = yMaxHist;
            }
            
            // Draw the histograms
            for(std::vector<TH1*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
                TH1* hist = *i_hist;
                const bool firstHist = i_hist==v_hist.begin();
                hist->SetLineWidth(2);
                
                if(firstHist){
                    hist->SetTitle("");
                    hist->GetYaxis()->SetRangeUser(0.,1.4*yMax);
                    hist->Draw();
                }
                else{
                    hist->Draw("sameS");
                }
                canvas->Modified();
                canvas->Update();
            }
            
            
            // Adjust stats boxes, has to be done only after drawing
            std::vector<TPaveStats*> v_stats;
            for(size_t iHist = 0; iHist < v_hist.size(); ++iHist){
                TH1* hist = v_hist.at(iHist);
                
                TPaveStats* stats = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
                const double y2ndc = 0.99 - static_cast<double>(iHist)*0.21;
                stats->SetY2NDC(y2ndc);
                stats->SetY1NDC(y2ndc-0.2);
                stats->SetLineColor(hist->GetLineColor());
                
                canvas->Modified();
                canvas->Update();
            }
            
            legend->Draw("same");
            canvas->Modified();
            canvas->Update();
            
            const TString plotName(outputFolder + fileShort + "/" + histName);
            
            canvas->Print(plotName + ".eps");
            canvas->Print(plotName + ".png");
            
            legend->Delete();
            for(size_t iHist = 0; iHist < v_hist.size(); ++iHist) v_hist.at(iHist)->Delete();
            v_hist.clear();
            for(size_t iStats = 0; iStats < v_stats.size(); ++iStats) v_stats.at(iStats)->Delete();
            v_stats.clear();
            canvas->Close();
        }
    }
}



void plotBestWeightHistos(std::map<TString, std::vector<TString> >& m_mvaConfigCorrect,
                          std::map<TString, std::vector<TString> >& m_mvaConfigSwapped,
                          const std::vector<std::pair<TString, std::vector<TString> > >& v_input,
                          const TString& outputFolder,
                          const TString& histNameBase)
{
    gROOT->SetStyle("Plain");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = 1001;
    
    gStyle->SetPalette(1);      //Spektralpalette, Default: 0 resp. 50
    //gStyle->SetNumberContours(20);  // Default: 20
    
    const double width = 600.;

    gStyle->SetCanvasDefW(width);
    gStyle->SetCanvasDefH(width);

    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.10);

    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.15);

    gStyle->SetTitleOffset(1.4,"Y");
    gStyle->SetTitleOffset(1.2,"X");

    TGaxis::SetMaxDigits(3);
    
    
    
    RootFileReader* rootFileReader(RootFileReader::getInstance());
    
    for(const auto& input : v_input){
        const TString& fileShort = input.first;
        const std::vector<TString>& v_inputFilename = input.second;
        
        TCanvas* canvas(0);
        TLegend* legend(0);
        
        std::map<TString, std::map<TString, TH1*> > m_correct;
        for(const auto& mvaConfig : m_mvaConfigCorrect){
            const TString& step = mvaConfig.first;
            const std::vector<TString>& v_configName = mvaConfig.second;
            
            std::map<TString, TH1*> m_hist;
            for(const auto& configName : v_configName){
                TString histName = histNameBase;
                histName.Append("_best");
                histName.Append("Correct");
                histName.Append(configName);
                histName.Append(step);
                
                TH1* hist(0);
                for(const auto& filename : v_inputFilename){
                    if(!hist){
                        hist = rootFileReader->GetClone<TH1D>(filename, histName);
                    }
                    else hist->Add(rootFileReader->GetClone<TH1D>(filename, histName));
                }
                m_hist[configName] = hist;
            }
            m_correct[step] = m_hist;
        }
        
        std::map<TString, std::map<TString, TH1*> > m_swapped;
        for(const auto& mvaConfig : m_mvaConfigSwapped){
            const TString& step = mvaConfig.first;
            const std::vector<TString>& v_configName = mvaConfig.second;
            
            std::map<TString, TH1*> m_hist;
            for(const auto& configName : v_configName){
                TString histName = histNameBase;
                histName.Append("_best");
                histName.Append("Swapped");
                histName.Append(configName);
                histName.Append(step);
                
                TH1* hist(0);
                for(const auto& filename : v_inputFilename){
                    if(!hist){
                        hist = rootFileReader->GetClone<TH1D>(filename, histName);
                    }
                    else hist->Add(rootFileReader->GetClone<TH1D>(filename, histName));
                }
                m_hist[configName] = hist;
            }
            m_swapped[step] = m_hist;
        }
        
        std::map<TString, std::map<TString, std::map<TString, TH1*> > > m_combined;
        for(const auto& mvaConfig1 : m_mvaConfigCorrect){
            const TString& step = mvaConfig1.first;
            const std::vector<TString>& v_configName1 = mvaConfig1.second;
            if(m_swapped.find(step) != m_swapped.end()){
                const std::vector<TString>& v_configName2 = m_mvaConfigSwapped.at(step);
                
                for(const auto& configName1 : v_configName1){
                    for(const auto& configName2 : v_configName2){
                        TString histName = histNameBase;
                        histName.Append("_best");
                        histName.Append("Combined");
                        histName.Append(configName1);
                        histName.Append(configName2);
                        histName.Append(step);
                        
                        TH1* hist(0);
                        for(const auto& filename : v_inputFilename){
                            if(!hist){
                                hist = rootFileReader->GetClone<TH1D>(filename, histName);
                            }
                            else hist->Add(rootFileReader->GetClone<TH1D>(filename, histName));
                        }
                        m_combined[step][configName1][configName2] = hist;
                    }
                }
            }
        }
        
        // Plot all correct trainings overlaid
        for(const auto& mvaConfig : m_mvaConfigCorrect){
            const TString& step = mvaConfig.first;
            const std::vector<TString>& v_configName = mvaConfig.second;
            
            // Access all extrema in order to print them
            double minFraction1(999.);
            double maxFraction1(-999.);
            double minFraction2(999.);
            double maxFraction2(-999.);
            double minFraction3(999.);
            double maxFraction3(-999.);
            
            // Set up canvas
            canvas = new TCanvas("canvas1");
            canvas->cd();
            
            // Set up legend
            legend = new TLegend(0.84,0.30,0.99,0.99);
            legend->SetFillColor(0);
            legend->SetFillStyle(0);
            legend->SetTextSize(0.04);
            legend->SetMargin(0.12);
            legend->SetFillStyle(1001);
            
            TString histName = histNameBase;
            histName.ReplaceAll("mvaA_", "");
            histName.Append("_correct");
            histName.Append(step);
            
            int counter(1);
            std::vector<TH1*> v_hist;
            for(const auto& configName : v_configName){
                TH1* hist = m_correct.at(step).at(configName);
                hist->SetLineColor(counter);
                legend->AddEntry(hist, configName, "l");
                v_hist.push_back(hist);
                ++counter;
            }
            
            if(histNameBase == "mvaA_jetsFromTop"){
                std::ofstream outfile;
                TString fileString = outputFolder + fileShort + "/result_jetsFromTop";
                fileString.Append("_correct");
                fileString.Append(step).Append(".txt");
                outfile.open(fileString.Data());
                for(const auto iter : m_correct.at(step)){
                    const TString& configName = iter.first;
                    TH1* hist = iter.second;
                    const double fraction1 = hist->GetBinContent(2)/hist->GetBinContent(1)*100.;
                    const double fraction2 = hist->GetBinContent(3)/hist->GetBinContent(1)*100.;
                    const double fraction3 = hist->GetBinContent(4)/hist->GetBinContent(1)*100.;
                    if(fraction1 < minFraction1) minFraction1 = fraction1;
                    if(fraction1 > maxFraction1) maxFraction1 = fraction1;
                    if(fraction2 < minFraction2) minFraction2 = fraction2;
                    if(fraction2 > maxFraction2) maxFraction2 = fraction2;
                    if(fraction3 < minFraction3) minFraction3 = fraction3;
                    if(fraction3 > maxFraction3) maxFraction3 = fraction3;
                    outfile<<std::setw(10)<<configName<<std::setw(10)<<fraction1
                           <<std::setw(10)<<fraction2<<std::setw(10)<<fraction3<<std::endl;
                }
                outfile.close();
            }
            else if(histNameBase == "mvaA_jetsFromHiggs"){
                std::ofstream outfile;
                TString fileString = outputFolder + fileShort + "/result_jetsFromHiggs";
                fileString.Append("_correct");
                fileString.Append(step).Append(".txt");
                outfile.open(fileString.Data());
                for(const auto iter : m_correct.at(step)){
                    const TString& configName = iter.first;
                    TH1* hist = iter.second;
                    const double fraction1 = hist->GetBinContent(2)/hist->GetBinContent(1)*100.;
                    const double fraction2 = hist->GetBinContent(3)/hist->GetBinContent(1)*100.;
                    const double fraction3 = hist->GetBinContent(4)/hist->GetBinContent(1)*100.;
                    if(fraction1 < minFraction1) minFraction1 = fraction1;
                    if(fraction1 > maxFraction1) maxFraction1 = fraction1;
                    if(fraction2 < minFraction2) minFraction2 = fraction2;
                    if(fraction2 > maxFraction2) maxFraction2 = fraction2;
                    if(fraction3 < minFraction3) minFraction3 = fraction3;
                    if(fraction3 > maxFraction3) maxFraction3 = fraction3;
                    outfile<<std::setw(10)<<configName<<std::setw(10)<<fraction1
                           <<std::setw(10)<<fraction2<<std::setw(10)<<fraction3<<std::endl;
                }
                outfile.close();
            }
            else if(histNameBase == "mvaA_jetsFromBoth"){
                std::ofstream outfile;
                TString fileString = outputFolder + fileShort + "/result_jetsFromBoth";
                fileString.Append("_correct");
                fileString.Append(step).Append(".txt");
                outfile.open(fileString.Data());
                for(const auto iter : m_correct.at(step)){
                    const TString& configName = iter.first;
                    TH1* hist = iter.second;
                    const double fraction1 = hist->GetBinContent(2)/hist->GetBinContent(1)*100.;
                    const double fraction2 = hist->GetBinContent(3)/hist->GetBinContent(1)*100.;
                    const double fraction3 = hist->GetBinContent(4)/hist->GetBinContent(1)*100.;
                    if(fraction1 < minFraction1) minFraction1 = fraction1;
                    if(fraction1 > maxFraction1) maxFraction1 = fraction1;
                    if(fraction2 < minFraction2) minFraction2 = fraction2;
                    if(fraction2 > maxFraction2) maxFraction2 = fraction2;
                    if(fraction3 < minFraction3) minFraction3 = fraction3;
                    if(fraction3 > maxFraction3) maxFraction3 = fraction3;
                    outfile<<std::setw(10)<<configName<<std::setw(10)<<fraction1
                           <<std::setw(10)<<fraction2<<std::setw(10)<<fraction3<<std::endl;
                }
                outfile.close();
            }
            
            std::cout<<"\nExtreme results for correct trainings:\n"<<fileShort<<"\n"<<outputFolder<<" , "<<step<<" , "<<histNameBase<<"\n";
            std::cout<<"Fraction 1 (min, max): "<<minFraction1<<" , "<<maxFraction1<<"\n";
            std::cout<<"Fraction 2 (min, max): "<<minFraction2<<" , "<<maxFraction2<<"\n";
            std::cout<<"Fraction 3 (min, max): "<<minFraction3<<" , "<<maxFraction3<<"\n\n";
            
            // Scale histograms and get minimum and maximum value
            Double_t yMin(99999.);
            Double_t yMax(-99999.);
            for(std::vector<TH1*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
                TH1* hist = *i_hist;
                //hist->Scale(1./hist->Integral(0, hist->GetNbinsX()+1));
                
                const Double_t yMinHist = hist->GetBinContent(hist->GetMinimumBin());
                const Double_t yMaxHist = hist->GetBinContent(hist->GetMaximumBin());
                if(yMin > yMinHist) yMin = yMinHist;
                if(yMax < yMaxHist) yMax = yMaxHist;
            }
            
            // Draw the histograms
            for(std::vector<TH1*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
                TH1* hist = *i_hist;
                const bool firstHist = i_hist==v_hist.begin();
                hist->SetLineWidth(2);
                
                if(firstHist){
                    hist->SetTitle("");
                    hist->GetYaxis()->SetRangeUser(0.,1.4*yMax);
                    hist->Draw();
                }
                else{
                    hist->Draw("sameS");
                }
                canvas->Modified();
                canvas->Update();
            }
            
            
            legend->Draw("same");
            canvas->Modified();
            canvas->Update();
            
            const TString plotName(outputFolder + fileShort + "/" + histName);
            
            canvas->Print(plotName + ".eps");
            canvas->Print(plotName + ".png");
            
            legend->Delete();
            canvas->Close();
        }
        
        // Plot all swapped trainings overlaid
        for(const auto& mvaConfig : m_mvaConfigSwapped){
            const TString& step = mvaConfig.first;
            const std::vector<TString>& v_configName = mvaConfig.second;
            
            // Access all extrema in order to print them
            double minFraction1(999.);
            double maxFraction1(-999.);
            double minFraction2(999.);
            double maxFraction2(-999.);
            double minFraction3(999.);
            double maxFraction3(-999.);
            
            // Set up canvas
            canvas = new TCanvas("canvas1");
            canvas->cd();
            
            // Set up legend
            legend = new TLegend(0.84,0.30,0.99,0.99);
            legend->SetFillColor(0);
            legend->SetFillStyle(0);
            legend->SetTextSize(0.04);
            legend->SetMargin(0.12);
            legend->SetFillStyle(1001);
            
            TString histName = histNameBase;
            histName.ReplaceAll("mvaA_", "");
            histName.Append("_swapped");
            histName.Append(step);
            
            int counter(1);
            std::vector<TH1*> v_hist;
            for(const auto& configName : v_configName){
                TH1* hist = m_swapped.at(step).at(configName);
                hist->SetLineColor(counter);
                legend->AddEntry(hist, configName, "l");
                v_hist.push_back(hist);
                ++counter;
            }
            
            if(histNameBase == "mvaA_jetsFromTop"){
                std::ofstream outfile;
                TString fileString = outputFolder + fileShort + "/result_jetsFromTop";
                fileString.Append("_swapped");
                fileString.Append(step).Append(".txt");
                outfile.open(fileString.Data());
                for(const auto iter : m_swapped.at(step)){
                    const TString& configName = iter.first;
                    TH1* hist = iter.second;
                    const double fraction1 = hist->GetBinContent(2)/hist->GetBinContent(1)*100.;
                    const double fraction2 = hist->GetBinContent(3)/hist->GetBinContent(1)*100.;
                    const double fraction3 = hist->GetBinContent(4)/hist->GetBinContent(1)*100.;
                    if(fraction1 < minFraction1) minFraction1 = fraction1;
                    if(fraction1 > maxFraction1) maxFraction1 = fraction1;
                    if(fraction2 < minFraction2) minFraction2 = fraction2;
                    if(fraction2 > maxFraction2) maxFraction2 = fraction2;
                    if(fraction3 < minFraction3) minFraction3 = fraction3;
                    if(fraction3 > maxFraction3) maxFraction3 = fraction3;
                    outfile<<std::setw(10)<<configName<<std::setw(10)<<fraction1
                           <<std::setw(10)<<fraction2<<std::setw(10)<<fraction3<<std::endl;
                }
                outfile.close();
            }
            else if(histNameBase == "mvaA_jetsFromHiggs"){
                std::ofstream outfile;
                TString fileString = outputFolder + fileShort + "/result_jetsFromHiggs";
                fileString.Append("_swapped");
                fileString.Append(step).Append(".txt");
                outfile.open(fileString.Data());
                for(const auto iter : m_swapped.at(step)){
                    const TString& configName = iter.first;
                    TH1* hist = iter.second;
                    const double fraction1 = hist->GetBinContent(2)/hist->GetBinContent(1)*100.;
                    const double fraction2 = hist->GetBinContent(3)/hist->GetBinContent(1)*100.;
                    const double fraction3 = hist->GetBinContent(4)/hist->GetBinContent(1)*100.;
                    if(fraction1 < minFraction1) minFraction1 = fraction1;
                    if(fraction1 > maxFraction1) maxFraction1 = fraction1;
                    if(fraction2 < minFraction2) minFraction2 = fraction2;
                    if(fraction2 > maxFraction2) maxFraction2 = fraction2;
                    if(fraction3 < minFraction3) minFraction3 = fraction3;
                    if(fraction3 > maxFraction3) maxFraction3 = fraction3;
                    outfile<<std::setw(10)<<configName<<std::setw(10)<<fraction1
                           <<std::setw(10)<<fraction2<<std::setw(10)<<fraction3<<std::endl;
                }
                outfile.close();
            }
            else if(histNameBase == "mvaA_jetsFromBoth"){
                std::ofstream outfile;
                TString fileString = outputFolder + fileShort + "/result_jetsFromBoth";
                fileString.Append("_swapped");
                fileString.Append(step).Append(".txt");
                outfile.open(fileString.Data());
                for(const auto iter : m_swapped.at(step)){
                    const TString& configName = iter.first;
                    TH1* hist = iter.second;
                    const double fraction1 = hist->GetBinContent(2)/hist->GetBinContent(1)*100.;
                    const double fraction2 = hist->GetBinContent(3)/hist->GetBinContent(1)*100.;
                    const double fraction3 = hist->GetBinContent(4)/hist->GetBinContent(1)*100.;
                    if(fraction1 < minFraction1) minFraction1 = fraction1;
                    if(fraction1 > maxFraction1) maxFraction1 = fraction1;
                    if(fraction2 < minFraction2) minFraction2 = fraction2;
                    if(fraction2 > maxFraction2) maxFraction2 = fraction2;
                    if(fraction3 < minFraction3) minFraction3 = fraction3;
                    if(fraction3 > maxFraction3) maxFraction3 = fraction3;
                    outfile<<std::setw(10)<<configName<<std::setw(10)<<fraction1
                           <<std::setw(10)<<fraction2<<std::setw(10)<<fraction3<<std::endl;
                }
                outfile.close();
            }
            
            std::cout<<"\nExtreme results for swapped trainings:\n"<<fileShort<<"\n"<<outputFolder<<" , "<<step<<" , "<<histNameBase<<"\n";
            std::cout<<"Fraction 1 (min, max): "<<minFraction1<<" , "<<maxFraction1<<"\n";
            std::cout<<"Fraction 2 (min, max): "<<minFraction2<<" , "<<maxFraction2<<"\n";
            std::cout<<"Fraction 3 (min, max): "<<minFraction3<<" , "<<maxFraction3<<"\n\n";
            
            // Scale histograms and get minimum and maximum value
            Double_t yMin(99999.);
            Double_t yMax(-99999.);
            for(std::vector<TH1*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
                TH1* hist = *i_hist;
                //hist->Scale(1./hist->Integral(0, hist->GetNbinsX()+1));
                
                const Double_t yMinHist = hist->GetBinContent(hist->GetMinimumBin());
                const Double_t yMaxHist = hist->GetBinContent(hist->GetMaximumBin());
                if(yMin > yMinHist) yMin = yMinHist;
                if(yMax < yMaxHist) yMax = yMaxHist;
            }
            
            // Draw the histograms
            for(std::vector<TH1*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
                TH1* hist = *i_hist;
                const bool firstHist = i_hist==v_hist.begin();
                hist->SetLineWidth(2);
                
                if(firstHist){
                    hist->SetTitle("");
                    hist->GetYaxis()->SetRangeUser(0.,1.4*yMax);
                    hist->Draw();
                }
                else{
                    hist->Draw("sameS");
                }
                canvas->Modified();
                canvas->Update();
            }
            
            
            legend->Draw("same");
            canvas->Modified();
            canvas->Update();
            
            const TString plotName(outputFolder + fileShort + "/" + histName);
            
            canvas->Print(plotName + ".eps");
            canvas->Print(plotName + ".png");
            
            legend->Delete();
            canvas->Close();
        }
        
        // Plot all combined trainings overlaid
        for(const auto& mvaConfig1 : m_mvaConfigCorrect){
            const TString& step = mvaConfig1.first;
            const std::vector<TString>& v_configName1 = mvaConfig1.second;
            
            if(m_combined.find(step) == m_combined.end()) continue;
            
            // Access all extrema in order to print them
            double minFraction1(999.);
            double maxFraction1(-999.);
            double minFraction2(999.);
            double maxFraction2(-999.);
            double minFraction3(999.);
            double maxFraction3(-999.);
            
            const std::vector<TString>& v_configName2 = m_mvaConfigSwapped.at(step);
            
            // Set up canvas
            canvas = new TCanvas("canvas1");
            canvas->cd();
            
            // Set up legend
            legend = new TLegend(0.84,0.30,0.99,0.99);
            legend->SetFillColor(0);
            legend->SetFillStyle(0);
            legend->SetTextSize(0.04);
            legend->SetMargin(0.12);
            legend->SetFillStyle(1001);
            
            TString histName = histNameBase;
            histName.ReplaceAll("mvaA_", "");
            histName.Append("_combined");
            histName.Append(step);
            
            int counter(1);
            std::vector<TH1*> v_hist;
            for(const auto& configName1 : v_configName1){
                for(const auto&  configName2 : v_configName2){
                    const TString configName = configName1 + configName2;
                    TH1* hist = m_combined.at(step).at(configName1).at(configName2);
                    hist->SetLineColor(counter);
                    legend->AddEntry(hist, configName, "l");
                    v_hist.push_back(hist);
                    ++counter;
                }
            }
            
            if(histNameBase == "mvaA_jetsFromTop"){
                std::ofstream outfile;
                TString fileString = outputFolder + fileShort + "/result_jetsFromTop";
                fileString.Append("_combined");
                fileString.Append(step).Append(".txt");
                outfile.open(fileString.Data());
                for(const auto iter1 : m_combined.at(step)){
                    for(const auto iter2 : iter1.second){
                        const TString configName = iter1.first + iter2.first;
                        TH1* hist = iter2.second;
                        const double fraction1 = hist->GetBinContent(2)/hist->GetBinContent(1)*100.;
                        const double fraction2 = hist->GetBinContent(3)/hist->GetBinContent(1)*100.;
                        const double fraction3 = hist->GetBinContent(4)/hist->GetBinContent(1)*100.;
                        if(fraction1 < minFraction1) minFraction1 = fraction1;
                        if(fraction1 > maxFraction1) maxFraction1 = fraction1;
                        if(fraction2 < minFraction2) minFraction2 = fraction2;
                        if(fraction2 > maxFraction2) maxFraction2 = fraction2;
                        if(fraction3 < minFraction3) minFraction3 = fraction3;
                        if(fraction3 > maxFraction3) maxFraction3 = fraction3;
                        outfile<<std::setw(10)<<configName<<std::setw(10)<<fraction1
                               <<std::setw(10)<<fraction2<<std::setw(10)<<fraction3<<std::endl;
                    }
                }
                outfile.close();
            }
            else if(histNameBase == "mvaA_jetsFromHiggs"){
                std::ofstream outfile;
                TString fileString = outputFolder + fileShort + "/result_jetsFromHiggs";
                fileString.Append("_combined");
                fileString.Append(step).Append(".txt");
                outfile.open(fileString.Data());
                for(const auto iter1 : m_combined.at(step)){
                    for(const auto iter2 : iter1.second){
                        const TString configName = iter1.first + iter2.first;
                        TH1* hist = iter2.second;
                        const double fraction1 = hist->GetBinContent(2)/hist->GetBinContent(1)*100.;
                        const double fraction2 = hist->GetBinContent(3)/hist->GetBinContent(1)*100.;
                        const double fraction3 = hist->GetBinContent(4)/hist->GetBinContent(1)*100.;
                        if(fraction1 < minFraction1) minFraction1 = fraction1;
                        if(fraction1 > maxFraction1) maxFraction1 = fraction1;
                        if(fraction2 < minFraction2) minFraction2 = fraction2;
                        if(fraction2 > maxFraction2) maxFraction2 = fraction2;
                        if(fraction3 < minFraction3) minFraction3 = fraction3;
                        if(fraction3 > maxFraction3) maxFraction3 = fraction3;
                        outfile<<std::setw(10)<<configName<<std::setw(10)<<fraction1
                               <<std::setw(10)<<fraction2<<std::setw(10)<<fraction3<<std::endl;
                    }
                }
                outfile.close();
            }
            else if(histNameBase == "mvaA_jetsFromBoth"){
                std::ofstream outfile;
                TString fileString = outputFolder + fileShort + "/result_jetsFromBoth";
                fileString.Append("_combined");
                fileString.Append(step).Append(".txt");
                outfile.open(fileString.Data());
                for(const auto iter1 : m_combined.at(step)){
                    for(const auto iter2 : iter1.second){
                        const TString configName = iter1.first + iter2.first;
                        TH1* hist = iter2.second;
                        const double fraction1 = hist->GetBinContent(2)/hist->GetBinContent(1)*100.;
                        const double fraction2 = hist->GetBinContent(3)/hist->GetBinContent(1)*100.;
                        const double fraction3 = hist->GetBinContent(4)/hist->GetBinContent(1)*100.;
                        if(fraction1 < minFraction1) minFraction1 = fraction1;
                        if(fraction1 > maxFraction1) maxFraction1 = fraction1;
                        if(fraction2 < minFraction2) minFraction2 = fraction2;
                        if(fraction2 > maxFraction2) maxFraction2 = fraction2;
                        if(fraction3 < minFraction3) minFraction3 = fraction3;
                        if(fraction3 > maxFraction3) maxFraction3 = fraction3;
                        outfile<<std::setw(10)<<configName<<std::setw(10)<<fraction1
                               <<std::setw(10)<<fraction2<<std::setw(10)<<fraction3<<std::endl;
                    }
                }
                outfile.close();
            }
            
            std::cout<<"\nExtreme results for combined trainings:\n"<<fileShort<<"\n"<<outputFolder<<" , "<<step<<" , "<<histNameBase<<"\n";
            std::cout<<"Fraction 1 (min, max): "<<minFraction1<<" , "<<maxFraction1<<"\n";
            std::cout<<"Fraction 2 (min, max): "<<minFraction2<<" , "<<maxFraction2<<"\n";
            std::cout<<"Fraction 3 (min, max): "<<minFraction3<<" , "<<maxFraction3<<"\n\n";
            
            // Scale histograms and get minimum and maximum value
            Double_t yMin(99999.);
            Double_t yMax(-99999.);
            for(std::vector<TH1*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
                TH1* hist = *i_hist;
                //hist->Scale(1./hist->Integral(0, hist->GetNbinsX()+1));
                
                const Double_t yMinHist = hist->GetBinContent(hist->GetMinimumBin());
                const Double_t yMaxHist = hist->GetBinContent(hist->GetMaximumBin());
                if(yMin > yMinHist) yMin = yMinHist;
                if(yMax < yMaxHist) yMax = yMaxHist;
            }
            
            // Draw the histograms
            for(std::vector<TH1*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
                TH1* hist = *i_hist;
                const bool firstHist = i_hist==v_hist.begin();
                hist->SetLineWidth(2);
                
                if(firstHist){
                    hist->SetTitle("");
                    hist->GetYaxis()->SetRangeUser(0.,1.4*yMax);
                    hist->Draw();
                }
                else{
                    hist->Draw("sameS");
                }
                canvas->Modified();
                canvas->Update();
            }
            
            
            legend->Draw("same");
            canvas->Modified();
            canvas->Update();
            
            const TString plotName(outputFolder + fileShort + "/" + histName);
            
            canvas->Print(plotName + ".eps");
            canvas->Print(plotName + ".png");
            
            legend->Delete();
            canvas->Close();
        }
        
        // Plot set of combined and corresponding correct and swapped training
        for(const auto& stepConfigs : m_combined){
            const TString& step = stepConfigs.first;
            
            for(const auto& config1 : stepConfigs.second){
                const TString& configName1 = config1.first;
                
                for(const auto& config2 : config1.second){
                    const TString& configName2 = config2.first;
                    
                    const TString configName = configName1 + configName2;
                    
                    // Set up canvas
                    canvas = new TCanvas("canvas1");
                    canvas->cd();
                    
                    // Set up legend
                    legend = new TLegend(0.41,0.84,0.76,0.99);
                    legend->SetFillColor(0);
                    legend->SetFillStyle(0);
                    legend->SetTextSize(0.04);
                    legend->SetMargin(0.12);
                    legend->SetFillStyle(1001);
                    
                    TString histName = histNameBase;
                    histName.ReplaceAll("mvaA_", "");
                    histName.Append(configName);
                    histName.Append(step);
                    
                    std::vector<TH1*> v_hist;
                    {
                        TH1* histCorrect = m_correct.at(step).at(configName1);
                        histCorrect->SetLineColor(4);
                        legend->AddEntry(histCorrect, "Correct", "l");
                        v_hist.push_back(histCorrect);
                        TH1* histSwapped = m_swapped.at(step).at(configName2);
                        histSwapped->SetLineColor(3);
                        legend->AddEntry(histSwapped, "Swapped", "l");
                        v_hist.push_back(histSwapped);
                        TH1* histCombined = config2.second;
                        histCombined->SetLineColor(1);
                        legend->AddEntry(histCombined, "Combined", "l");
                        v_hist.push_back(histCombined);
                    }
                    
                    // Scale histograms and get minimum and maximum value
                    Double_t yMin(99999.);
                    Double_t yMax(-99999.);
                    for(std::vector<TH1*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
                        TH1* hist = *i_hist;
                        //hist->Scale(1./hist->Integral(0, hist->GetNbinsX()+1));
                        
                        const Double_t yMinHist = hist->GetBinContent(hist->GetMinimumBin());
                        const Double_t yMaxHist = hist->GetBinContent(hist->GetMaximumBin());
                        if(yMin > yMinHist) yMin = yMinHist;
                        if(yMax < yMaxHist) yMax = yMaxHist;
                    }
                    
                    // Draw the histograms
                    for(std::vector<TH1*>::iterator i_hist = v_hist.begin(); i_hist != v_hist.end(); ++i_hist){
                        TH1* hist = *i_hist;
                        const bool firstHist = i_hist==v_hist.begin();
                        hist->SetLineWidth(2);
                        
                        if(firstHist){
                            hist->SetTitle("");
                            hist->GetYaxis()->SetRangeUser(0.,1.4*yMax);
                            hist->Draw();
                        }
                        else{
                            hist->Draw("sameS");
                        }
                        canvas->Modified();
                        canvas->Update();
                    }
                    
                    
                    legend->Draw("same");
                    canvas->Modified();
                    canvas->Update();
                    
                    const TString plotName(outputFolder + fileShort + "/" + histName);
                    
                    canvas->Print(plotName + ".eps");
                    canvas->Print(plotName + ".png");
                    
                    legend->Delete();
                    canvas->Close();
                }
            }
        }
    }
}


void histoAnalyzerMvaTopJets(const std::vector<Channel::Channel>& v_channel,
                             const std::vector<Systematic::Systematic>& v_systematic,
                             const std::vector<std::string>& v_mode)
{
    // Hardcoded input files, since detailed validation makes sense only on samples containing ttbar system, and some only on real ttH events
    const std::vector<TString> v_inputFileTtbar =
        {
//            "ttbarW",
//            "ttbarZ",
//            "ttbarH125inclusiveBbbar",
//            "ttbarH125inclusiveOther",
            "ttbarH125tobbbar",
            "ttbarsignalPlusBbbar",
            "ttbarsignalPlusOther"
        };
    
    // Loop over all systematics and channels
    for(const auto& systematic : v_systematic){
        for(const auto& channel : v_channel){
            
            // Access input folders
            std::vector<TString> v_inputFolder;
            if(channel != Channel::combined){
                TString inputFolder = common::accessFolder(InputBaseDIR, channel, systematic);
                inputFolder.Append(Channel::convertChannel(channel)).Append("_");
                v_inputFolder.push_back(inputFolder);
            }
            else{
                for(const auto& realChannel : {Channel::ee, Channel::emu, Channel::mumu}){
                    TString inputFolder = common::accessFolder(InputBaseDIR, realChannel, systematic);
                    inputFolder.Append(Channel::convertChannel(realChannel)).Append("_");
                    v_inputFolder.push_back(inputFolder);
                }
            }
            
            // Access all filenames, and create output directories
            const TString outputFolder = common::assignFolder(OutputBaseDIR, channel, systematic);
            std::vector<std::pair<TString, std::vector<TString> > > v_input;
            for(const auto& inputFileTtbar : v_inputFileTtbar){
                TString path = outputFolder;
                path.Append(inputFileTtbar).Append("/");
                gSystem->MakeDirectory(path.Data());
                
                std::vector<TString> v_filename;
                for(const auto& inputFolder : v_inputFolder){
                    TString filename = inputFolder;
                    filename.Append(inputFileTtbar).Append(".root");
                    v_filename.push_back(filename);
                    std::cout<<"Filename: "<<filename<<"\n";
                }
                v_input.push_back(std::make_pair(inputFileTtbar, v_filename));
                //std::cout<<"Processing sample: "<<inputFileTtbar<<"\n";
            }
            
            
            if(std::find(v_mode.begin(), v_mode.end(), "weight") != v_mode.end()){
                // Find all histograms containing specific MVA weight from first file
                const std::vector<std::pair<TString, TString> > v_nameStepPair = 
                    tth::nameStepPairs(v_input.at(0).second.at(0), "mvaA_mvaWeight");
                
                // Access all histograms for all types of combinations
                std::vector<std::pair<TString, TString> > v_nameStepPairCorrect;
                std::vector<std::pair<TString, TString> > v_nameStepPairSwapped;
                std::vector<std::pair<TString, TString> > v_nameStepPairWrong;
                for(const auto& nameStepPair : v_nameStepPair){
                    const TString& name = nameStepPair.first;
                    if(name.Contains("_correct_")) v_nameStepPairCorrect.push_back(nameStepPair);
                    else if(name.Contains("_swapped_")) v_nameStepPairSwapped.push_back(nameStepPair);
                    else if(name.Contains("_wrong_")) v_nameStepPairWrong.push_back(nameStepPair);
                    //std::cout<<"Name, step: "<<nameStepPair.first<<" , "<<nameStepPair.second<<"\n";
                }
                
                // Plot the histograms
                plotWeightHistos(v_nameStepPairCorrect, v_nameStepPairSwapped, v_nameStepPairWrong, v_input, outputFolder);
            }
            
            
            if(std::find(v_mode.begin(), v_mode.end(), "weight2d") != v_mode.end()){
                // Find all histograms containing specific MVA weight from first file
                const std::vector<std::pair<TString, TString> > v_nameStepPair = 
                    tth::nameStepPairs(v_input.at(0).second.at(0), "mvaA_mvaWeightCorrectVsSwapped");
                
                // Plot the histograms
                plotWeight2dHistos(v_nameStepPair, v_input, outputFolder);
            }
            
            
            if(std::find(v_mode.begin(), v_mode.end(), "best") != v_mode.end()){
                
                // Find all correct and swapped trainings from histograms containing dijet mass with best MVA weight from first file
                const std::vector<std::pair<TString, TString> > v_nameStepPair1 =
                    tth::nameStepPairs(v_input.at(0).second.at(0), "mvaA_dijet_mass_best");
                std::map<TString, std::vector<TString> > m_mvaConfigCorrect;
                std::map<TString, std::vector<TString> > m_mvaConfigSwapped;
                for(const auto& nameStepPair : v_nameStepPair1){
                    const TString& name = nameStepPair.first;
                    const TString& step = nameStepPair.second;
                    //std::cout<<"Name, step: "<<nameStepPair.first<<" , "<<nameStepPair.second<<"\n";
                    if(name.Contains("_bestCorrect_")){
                        TString mvaConfigName = name;
                        mvaConfigName.ReplaceAll("mvaA_dijet_mass_bestCorrect", "");
                        mvaConfigName.ReplaceAll(nameStepPair.second, "");
                        if(m_mvaConfigCorrect.find(step) != m_mvaConfigCorrect.end()){
                            m_mvaConfigCorrect.at(step).push_back(mvaConfigName);
                        }
                        else{
                            m_mvaConfigCorrect[step].push_back(mvaConfigName);
                        }
                        //std::cout<<"Correct training: "<<step<<" , "<<mvaConfigName<<"\n";
                    }
                    else if(name.Contains("_bestSwapped_")){
                        TString mvaConfigName = name;
                        mvaConfigName.ReplaceAll("mvaA_dijet_mass_bestSwapped", "");
                        mvaConfigName.ReplaceAll(nameStepPair.second, "");
                        if(m_mvaConfigSwapped.find(step) != m_mvaConfigSwapped.end()){
                            m_mvaConfigSwapped.at(step).push_back(mvaConfigName);
                        }
                        else{
                            m_mvaConfigSwapped[step].push_back(mvaConfigName);
                        }
                        //std::cout<<"Swapped training: "<<step<<" , "<<mvaConfigName<<"\n";
                    }
                }
                
                // Plot the histograms
                plotBestWeightHistos(m_mvaConfigCorrect, m_mvaConfigSwapped, v_input, outputFolder, "mvaA_dijet_mass");
                plotBestWeightHistos(m_mvaConfigCorrect, m_mvaConfigSwapped, v_input, outputFolder, "mvaA_jetsFromTop");
                plotBestWeightHistos(m_mvaConfigCorrect, m_mvaConfigSwapped, v_input, outputFolder, "mvaA_jetsFromHiggs");
                plotBestWeightHistos(m_mvaConfigCorrect, m_mvaConfigSwapped, v_input, outputFolder, "mvaA_jetsFromBoth");
            }
        }
    }
    
    std::cout<<"MVA validation plotting successfully finished";
}



int main(int argc, char** argv)
{
    CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu, combined. Default: combined", false, 1, 4,
        common::makeStringCheck(Channel::convertChannels(Channel::allowedChannelsPlotting)));
    CLParameter<std::string> opt_systematic("s", "Systematic variation - default is Nominal", false, 1, 100,
        common::makeStringCheck({"Nominal", ""}));
    CLParameter<std::string> opt_mode("m", "Which modes to run. 1D weight histos (weight), 2D weight histos(weight2d), best weight validation [default] (best)", false, 1, 3,
        common::makeStringCheck({"weight", "weight2d", "best", ""}));
    CLAnalyser::interpretGlobal(argc, argv);
    
    // Set up channels
    std::vector<Channel::Channel> v_channel({Channel::combined});
    if(opt_channel.isSet()) v_channel = Channel::convertChannels(opt_channel.getArguments());
    std::cout << "Processing channels: ";
    for (auto channel : v_channel)std::cout << Channel::convertChannel(channel) << " ";
    std::cout << "\n\n";
    
    // Set up systematics
    std::vector<Systematic::Systematic> v_systematic({Systematic::nominal});
    if(opt_systematic.isSet()) v_systematic = Systematic::convertSystematics(opt_systematic.getArguments());
    std::cout << "Processing systematics: ";
    for (auto systematic : v_systematic) std::cout << Systematic::convertSystematic(systematic) << " ";
    std::cout << "\n\n";
    
    // Set up modes for drawing
    std::vector<std::string> v_mode({"best"});
    if(opt_mode.isSet()) v_mode = opt_mode.getArguments();
    std::cout << "Processing modes: ";
    for (auto mode : v_mode) std::cout << mode << " ";
    std::cout << "\n\n";
    
    histoAnalyzerMvaTopJets(v_channel, v_systematic, v_mode);
}





