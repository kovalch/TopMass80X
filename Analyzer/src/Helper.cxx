#include "Helper.h"

#include "ProgramOptionsReader.h"

//#include "TH1F.h"
#include "TStyle.h"
#include "TPaveLabel.h"

#include <iostream>

typedef ProgramOptionsReader po;

void Helper::init()
{
  channelID_ = channelIDFromString(po::GetOption<std::string>("channel"));
  energy_ = po::GetOption<int>("cmsenergy");
}

int Helper::channelIDFromString(const std::string& channel)
{
  if      (!strncmp(channel.c_str(), "alljets" , 7)) return kAllJets;
  else if (!strncmp(channel.c_str(), "muon"    , 4)) return kMuonJets;
  else if (!strncmp(channel.c_str(), "electron", 8)) return kElectronJets;
  else if (!strncmp(channel.c_str(), "lepton"  , 6)) return kLeptonJets;
  else {
    std::cerr << "Channel name *" << channel << "* not know! Aborting program execution!" << std::endl;
    exit(1);
    return kMaxChannels;
  }
}

int Helper::methodIDFromString(const std::string& method)
{
  if      (!strcmp(method.c_str(), "GenMatch"   )) return kGenMatch;
  else if (!strcmp(method.c_str(), "MVA"        )) return kMVA;
  else if (!strcmp(method.c_str(), "Ideogram"   )) return kIdeogram;
  else if (!strcmp(method.c_str(), "IdeogramNew")) return kIdeogramNew;
  else if (!strcmp(method.c_str(), "IdeogramMin")) return kIdeogramMin;
  else if (!strcmp(method.c_str(), "RooFit"     )) return kRooFit;
  else {
    std::cerr << "Stopping analysis! Specified analysis method *" << method << "* not known!" << std::endl;
    exit(1);
    return kMaxMethods;
  }
}

int Helper::getCMSEnergy()
{
  return po::GetOption<int>("cmsenergy");
}

TH1F* Helper::GetH1(const std::string& title) {
  const float* array = &vBinning[0];

  TH1F* hHelper = new TH1F(title.c_str(), title.c_str(), vBinning.size()-1, array);
  hHelper->SetStats(false);
  hHelper->SetXTitle(fBinning.c_str());
  hHelper->SetYTitle(title.c_str());

  return hHelper;
}

//// tdrGrid: Turns the grid lines on (true) or off (false)
//
//void tdrGrid(bool gridOn) {
//  tdrStyle->SetPadGridX(gridOn);
//  tdrStyle->SetPadGridY(gridOn);
//}
//
//// fixOverlay: Redraws the axis
//
//void fixOverlay() {
//  gPad->RedrawAxis();
//}

void Helper::SetTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadRightMargin(0.25);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(true);
  tdrStyle->SetPadGridY(true);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
//  tdrStyle->SetErrorMarker(20);
  //tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.16);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.3);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(1005, "X");
  tdrStyle->SetNdivisions(510,  "YZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);
  
  tdrStyle->SetPaintTextFormat(".2f");
  tdrStyle->SetPalette(1);

  tdrStyle->cd();

}

void Helper::DrawLabel(const std::string& text, const double x1, const double y1, const double x2, Color_t color)
{
  // function to directly draw a label into the active canvas
  double y2 = y1 + 0.05;
  double yOffset = 0.02;
  TPaveLabel *label = new TPaveLabel(x1, y1+yOffset, x2, y2+yOffset, text.c_str(), "br NDC");
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  label->SetTextSize(0.75);
  label->SetTextAlign(12);
  label->SetTextColor(color);
  label->Draw("same");
}

void Helper::DrawCMS(int channelID, int energy) {
  if(channelID < 0) channelID = channelID_;
  if(energy    < 0) energy    = energy_;
  if(energy == 7){
    if(po::GetOption<int>("preliminary")){
      if     (channelID == kAllJets     ) DrawLabel("CMS Preliminary, 3.54 fb^{-1},  #sqrt{s} = 7 TeV"         , 0.2, 0.93, 0.9);
      else if(channelID == kLeptonJets  ) DrawLabel("CMS Preliminary, 5.0 fb^{-1},  #sqrt{s} = 7 TeV, l+jets"  , 0.2, 0.93, 0.9);
      else if(channelID == kElectronJets) DrawLabel("CMS Preliminary, 5.0 fb^{-1},  #sqrt{s} = 7 TeV, e+jets"  , 0.2, 0.93, 0.9);
      else if(channelID == kMuonJets    ) DrawLabel("CMS Preliminary, 5.0 fb^{-1},  #sqrt{s} = 7 TeV, #mu+jets", 0.2, 0.93, 0.9);
    }
    else{
      if     (channelID == kAllJets     ) DrawLabel("CMS, 3.54 fb^{-1},  #sqrt{s} = 7 TeV"         , 0.2, 0.93, 0.9);
      else if(channelID == kLeptonJets  ) DrawLabel("CMS, 5.0 fb^{-1},  #sqrt{s} = 7 TeV, l+jets"  , 0.2, 0.93, 0.9);
      else if(channelID == kElectronJets) DrawLabel("CMS, 5.0 fb^{-1},  #sqrt{s} = 7 TeV, e+jets"  , 0.2, 0.93, 0.9);
      else if(channelID == kMuonJets    ) DrawLabel("CMS, 5.0 fb^{-1},  #sqrt{s} = 7 TeV, #mu+jets", 0.2, 0.93, 0.9);
    }
  }
  else if(energy == 8){
    if(po::GetOption<int>("preliminary") == 3){
      if     (channelID == kAllJets     ) DrawLabel("PRIVATE WORK,  #sqrt{s} = 8 TeV"         , 0.2, 0.93, 0.9);
      else if(channelID == kLeptonJets  ) DrawLabel("PRIVATE WORK,  #sqrt{s} = 8 TeV, l+jets"  , 0.2, 0.93, 0.9);
      else if(channelID == kElectronJets) DrawLabel("PRIVATE WORK,  #sqrt{s} = 8 TeV, e+jets"  , 0.2, 0.93, 0.9);
      else if(channelID == kMuonJets    ) DrawLabel("PRIVATE WORK,  #sqrt{s} = 8 TeV, #mu+jets", 0.2, 0.93, 0.9);
    }
    else if(po::GetOption<int>("preliminary") == 2){
      if     (channelID == kAllJets     ) DrawLabel("PRIVATE WORK, 18.2 fb^{-1},  #sqrt{s} = 8 TeV"         , 0.2, 0.93, 0.9);
      else if(channelID == kLeptonJets  ) DrawLabel("PRIVATE WORK, 19.7 fb^{-1},  #sqrt{s} = 8 TeV, l+jets"  , 0.2, 0.93, 0.9);
      else if(channelID == kElectronJets) DrawLabel("PRIVATE WORK, 19.7 fb^{-1},  #sqrt{s} = 8 TeV, e+jets"  , 0.2, 0.93, 0.9);
      else if(channelID == kMuonJets    ) DrawLabel("PRIVATE WORK, 19.7 fb^{-1},  #sqrt{s} = 8 TeV, #mu+jets", 0.2, 0.93, 0.9);
    }
    else if(po::GetOption<int>("preliminary") == 1){
      if     (channelID == kAllJets     ) DrawLabel("CMS Preliminary, 18.2 fb^{-1},  #sqrt{s} = 8 TeV"         , 0.2, 0.93, 0.9);
      else if(channelID == kLeptonJets  ) DrawLabel("CMS Preliminary, 19.7 fb^{-1},  #sqrt{s} = 8 TeV, l+jets"  , 0.2, 0.93, 0.9);
      else if(channelID == kElectronJets) DrawLabel("CMS Preliminary, 19.7 fb^{-1},  #sqrt{s} = 8 TeV, e+jets"  , 0.2, 0.93, 0.9);
      else if(channelID == kMuonJets    ) DrawLabel("CMS Preliminary, 19.7 fb^{-1},  #sqrt{s} = 8 TeV, #mu+jets", 0.2, 0.93, 0.9);
    }
    else{
      if     (channelID == kAllJets     ) DrawLabel("CMS, 18.2 fb^{-1},  #sqrt{s} = 8 TeV"         , 0.2, 0.93, 0.9);
      else if(channelID == kLeptonJets  ) DrawLabel("CMS, 19.7 fb^{-1},  #sqrt{s} = 8 TeV, l+jets"  , 0.2, 0.93, 0.9);
      else if(channelID == kElectronJets) DrawLabel("CMS, 19.7 fb^{-1},  #sqrt{s} = 8 TeV, e+jets"  , 0.2, 0.93, 0.9);
      else if(channelID == kMuonJets    ) DrawLabel("CMS, 19.7 fb^{-1},  #sqrt{s} = 8 TeV, #mu+jets", 0.2, 0.93, 0.9);
    }
  }
}

void Helper::DrawCMSSim(int energy) {
  if(energy < 0) energy = energy_;
  if(energy == 7){
    if(po::GetOption<int>("preliminary")){
      DrawLabel("CMS Simulation Preliminary,  #sqrt{s} = 7 TeV", 0.2, 0.93, 0.9);
    }
    else{
      DrawLabel("CMS Simulation,  #sqrt{s} = 7 TeV", 0.2, 0.93, 0.9);
    }
  }
  else if(energy == 8){
    if(po::GetOption<int>("preliminary") == 2){
      DrawLabel("PRIVATE WORK,  #sqrt{s} = 8 TeV", 0.2, 0.93, 0.9);
    }
    else if(po::GetOption<int>("preliminary") == 1){
      DrawLabel("CMS Simulation Preliminary,  #sqrt{s} = 8 TeV", 0.2, 0.93, 0.9);
    }
    else{
      DrawLabel("CMS Simulation,  #sqrt{s} = 8 TeV", 0.2, 0.93, 0.9);
    }
  }
}

int Helper::channelID()
{
  return channelIDFromString(po::GetOption<std::string>("channel"));
}

int Helper::methodID()
{
  return methodIDFromString(po::GetOption<std::string>("method"));
}

std::vector<double> Helper::readParameters(const char *whichParameter){
  std::vector<double> pars;
  std::vector<std::string> vsPars = readParametersString(whichParameter);
  for(auto var : vsPars)
    pars.push_back(std::atof(var.c_str()));
  return pars;
}

std::vector<std::string> Helper::readParametersString(const char *whichParameter){
  std::vector<std::string> vsPars;
  std::string sPars = po::GetOption<std::string>(whichParameter);
  boost::split(vsPars, sPars, boost::is_any_of("|"));
  return vsPars;
}




TH1* HelperFunctions::createRatioPlot(const TH1 *h1, const TH1 *h2, const std::string &yTitle){
  assert( h1->GetNbinsX() == h2->GetNbinsX() );
  std::string name = (std::string("Ratio_") + h1->GetName())+h2->GetName();
  TH1 *hRatio = static_cast<TH1*>(h1->Clone(name.c_str()));
  TH1 *h2Temp = static_cast<TH1*>(h2->Clone("h2Temp"));
  if(h1==h2){
	  std::cout << "both histograms are the same; setting demoninator errors to zero" << std::endl;
	  for(int i=1; i<h2Temp->GetNbinsX();i++){
		  h2Temp->SetBinError(i,0);
	  }
  }
  hRatio->SetMarkerStyle(h1->GetMarkerStyle());
  hRatio->SetMarkerColor(h1->GetLineColor());
  hRatio->SetYTitle(yTitle.c_str());
//  hRatio->Divide(h2);
  hRatio->Divide(h2Temp);
  hRatio->GetYaxis()->SetRangeUser(0,-1);
  delete h2Temp;
  return hRatio;
}


std::string HelperFunctions::cleanedName(std::string toBeCleaned){
    boost::replace_all(toBeCleaned,"(","_"      );
    boost::replace_all(toBeCleaned,")","_"      );
    boost::replace_all(toBeCleaned,"/","_"      );
    boost::replace_all(toBeCleaned,"#","_"      );
    boost::replace_all(toBeCleaned," ","_"      );
    boost::replace_all(toBeCleaned,"{","_"      );
    boost::replace_all(toBeCleaned,"}","_"      );
    boost::replace_all(toBeCleaned,"^","_"      );
    boost::replace_all(toBeCleaned,".","_"      );
    boost::replace_all(toBeCleaned,"*","_"      );
    boost::replace_all(toBeCleaned,",",""       );
    boost::replace_all(toBeCleaned,";",""       );
    boost::replace_all(toBeCleaned,":",""       );
    boost::replace_all(toBeCleaned,"|",""       );
    boost::replace_all(toBeCleaned,"+","_"      );
    boost::replace_all(toBeCleaned,"<","_st_"   );
    boost::replace_all(toBeCleaned,">","_gt_"   );
    boost::replace_all(toBeCleaned,"[","_"      );
    boost::replace_all(toBeCleaned,"]","_"      );
    boost::replace_all(toBeCleaned,"@","_"      );
    boost::replace_all(toBeCleaned,"$","_"      );
    return toBeCleaned;
}

// -------------------------------------------------------------------------------------
void HelperFunctions::findYRange(const TH1 *h, double& min, double& max) {
  min = 1E10;
  max = 0.;
  for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
    double val = h->GetBinContent(bin);
    if( val < min && val!=0) min = val;
    if( val > max ) max = val;
  }
  if( min > max ) {
    min = 1E-3;
    max = 1;
  }
}


// -------------------------------------------------------------------------------------
void HelperFunctions::setCommonYRange(std::vector <TH1 *> histos, double RelTopOffset, double logMin) {
  if(histos.size()>0){
    double min = 0.;
    double max = 0.;
    findYRange(histos.at(0),min,max);
    for(unsigned int i = 1; i < histos.size(); i++) {
      double minTmp = 0;
      double maxTmp = 0;
      findYRange(histos.at(i),minTmp,maxTmp);
      if( minTmp < min ) min = minTmp;
      if( maxTmp > max ) max = maxTmp;
    }
    for(unsigned int i = 0; i < histos.size(); i++) {
      if(logMin>0) histos.at(i)->GetYaxis()->SetRangeUser(logMin, logMin* pow(max/min, 1./ (1-RelTopOffset)));
      else histos.at(i)->GetYaxis()->SetRangeUser(min-((max-min)*0.1),(max-min)/(1-RelTopOffset)+min);
    }
  }
}

// --------------------------------------------------
//    static bool fitCoreWidth(const TH1* hist, double nSig, TF1* &gauss, double &width, double &widthErr, double &rms, double &rmsErr, double &mean) {
bool HelperFunctions::fitCoreWidth(const TH1* hist, double nSig, TF1* &gauss, double &width, double &widthErr, double &rms, double &rmsErr) {
  bool result  = false;

  TString name = hist->GetName();
  name += "_GaussFit";
  gauss = new TF1(name,"gaus",hist->GetXaxis()->GetBinLowEdge(1),hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()));
  gauss->SetLineWidth(1);
  gauss->SetLineColor(kRed);

  TH1* h = static_cast<TH1*>(hist->Clone("util::HistOps::fitCoreWidth::h"));
  double mean = h->GetMean();
  rms = h->GetRMS();
  rmsErr = h->GetRMSError();
  double sig = 2.5*rms;
  if( h->Fit(gauss,"0QIR","",mean-sig,mean+sig) == 0 ) {
mean = gauss->GetParameter(1);
sig = nSig*gauss->GetParameter(2);
if( h->Fit(gauss,"0QIR","",mean-sig,mean+sig) == 0 ) {
  result = true;
  width = gauss->GetParameter(2);
  widthErr = gauss->GetParError(2);

  mean = gauss->GetParameter(1);
  sig = nSig*width;
  gauss->SetRange(mean-sig,mean+sig);
} else {
  std::cerr << "WARNING in util::HistOps::fitCoreWidth: No convergence when fitting width of '" << h->GetName() << "'\n";
  width = 0.;
  widthErr = 10000.;
}
  }
  delete h;

  return result;
}

// -------------------------------------------------------------------------------------
bool HelperFunctions::equidistLogBins(std::vector<double>& binEdges, double min, double max, bool logarithm) {
  if( binEdges.size() < 2 ) return false;
  if( (min <= 0. && logarithm) || (max <= 0. && logarithm) || min >= max ) return false;

  binEdges.front() = min;
  binEdges.back() = max;
  const double minLog = log10(binEdges.front());
  const double maxLog = log10(binEdges.back());
  for(unsigned int i = 1; i < binEdges.size()-1; ++i) {
	  if(logarithm)binEdges.at(i) = pow(10., minLog + i*(maxLog-minLog)/(binEdges.size()-1));
	  else binEdges.at(i) = min + i*(max-min)/(binEdges.size()-1);
  }
  return true;
}


// -------------------------------------------------------------------------------------
std::string HelperFunctions::addProperArrayIndex(std::string inputexpression, std::string arrayIndex) {
  boost::replace_all(inputexpression,"/",arrayIndex + "/");
  inputexpression+=arrayIndex;
  return inputexpression;
}


