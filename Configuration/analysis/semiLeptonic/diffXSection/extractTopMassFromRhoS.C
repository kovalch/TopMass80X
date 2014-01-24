#include "basicFunctions.h"

std::map <TString, TH1F*> extraction(int verbose=0, double luminosity=19712., bool save=true, double minx=160., double maxx=185., double nbinsx=250, bool dataMassDependence=false, int binOfInterest=1, TString outputFolder="./", TString outputFile="diffXSecTopSemiLepHadronPhaseSpace.root");
std::map<TString, double> GetChi2Info(TH1F* chi2, bool draw=false, int verbose=0, double xmin=-1., double xmax=-1.);

void extractTopMassFromRhoS(int verbose=0, double luminosity=19712., bool save=true){

  // define some parameters
  // range of plot
  double minx2=130.;
  double maxx2=250.;
  double minx=160.;
  double maxx=190.;

  // fine bins
  double nbinsx=250.;
  // consider mass dependence of measurement
  bool dataMassDependence=true;  
  // where to store the output
  TString outputFolder="./diffXSecFromSignal/plots/combined/2012/rhos/";
  TString outputFile="diffXSecTopSemiLepHadronPhaseSpace.root";
  // comparison mass value (Tevatron)
  double mass   =173.18;
  double masserr=0.94;

  // ============================
  //  Set Root Style
  // ============================
  TStyle myStyle("HHStyle2","HHStyle2");
  setHHStyle(myStyle);
  TGaxis::SetMaxDigits(2);
  myStyle.cd();
  gROOT->SetStyle("HHStyle2");
  gStyle->SetEndErrorSize(8);

  // create empty mass from bin plot
  std::map< TString, std::vector<double> > binning_ = makeVariableBinning();
  int Nbins = std::abs(binning_["rhos"][binning_["rhos"].size()-1]-binning_["rhos"][0])*100;
  // get histogram in correct binning
  TH1F* extractedMass= new TH1F("extractedMass", "extractedMass", Nbins, binning_["rhos"][0], binning_["rhos"][binning_["rhos"].size()-1]);
  reBinTH1F(*extractedMass, binning_["rhos"], 0);
  TGraphAsymmErrors* extractedChi2Mass= new TGraphAsymmErrors(nbinsx);
  TH1F* TevatronMass= new TH1F("TevatronMass", "TevatronMass", 1, binning_["rhos"][0], binning_["rhos"][binning_["rhos"].size()-1]);
  histogramStyle(*TevatronMass, kSig, false, 0.1);
  TevatronMass->SetBinContent(1,mass);
  TGraphAsymmErrors* TevatronMassErr= new TGraphAsymmErrors(nbinsx);  
  for(int bin=1; bin<=extractedMass->GetNbinsX(); ++bin){
    double binWidth             = extractedMass->GetBinWidth(bin);
    double binCenter            = extractedMass->GetBinCenter(bin);
    TevatronMassErr->SetPoint     (bin-1, binCenter, mass);
    TevatronMassErr->SetPointError(bin-1, 0.5*binWidth, 0.5*binWidth, masserr, masserr);
  }
  TevatronMassErr->SetLineColor(color_[kSig]);
  TevatronMassErr->SetMarkerColor(color_[kSig]);
  TevatronMassErr->SetFillColor(color_[kBkg]);
  TevatronMassErr->SetFillStyle(1001);
  TevatronMassErr->SetLineWidth(2);
  TevatronMassErr->SetLineStyle(1);
  TevatronMassErr->SetMarkerSize(0.1);

  //loop all bins
  TH1F* globalchi2=0;
  std::map <TString, TH1F*> resultbin;
  for(int bin=1; bin<=4; ++bin){
    // !!! define bin of interest !!! 
    int binOfInterest=bin;
    // perform top mass extraction
    resultbin=extraction(verbose, luminosity, save, (binOfInterest==2||binOfInterest==4) ? minx : minx2, (binOfInterest==2||binOfInterest==4) ? maxx : maxx2, nbinsx, dataMassDependence, binOfInterest, outputFolder, outputFile);
    // fill plot for intersection method
    extractedMass->SetBinContent(bin, resultbin["massval"]->GetBinContent(1));
    extractedMass->SetBinError  (bin, std::abs(resultbin["massunc"]->GetBinContent(1)));
    if(resultbin["chi2"]){
      // collect all chi2 in global chi2
      if(globalchi2) globalchi2->Add((TH1F*)resultbin["chi2"]->Clone("globalChi2Bin"+getTStringFromInt(bin)));
      else globalchi2=(TH1F*)resultbin["chi2"]->Clone("globalChi2");
      // fill plot for chi2 method
      double binWidth             = extractedMass->GetBinWidth(bin);
      double binCenter            = extractedMass->GetBinCenter(bin);
      std::map<TString, double> chi2mass=GetChi2Info(resultbin["chi2"], false, verbose, (binOfInterest==2||binOfInterest==4) ? minx : minx2, (binOfInterest==2||binOfInterest==4) ? maxx : maxx2);
      extractedChi2Mass->SetPoint     (bin-1, binCenter, chi2mass["min"]);
      extractedChi2Mass->SetPointError(bin-1, 0.5*binWidth, 0.5*binWidth, chi2mass["dn"],  chi2mass["up"]);
    }
  }

  std::cout << "combined result:" << std::endl;
  // draw global chi2
  std::vector<TCanvas*> plotCanvas_;
  addCanvas(plotCanvas_);
  plotCanvas_[plotCanvas_.size()-1]->cd(0);
  TString titlechi2="globalChi2RhosAllBins";
  plotCanvas_[plotCanvas_.size()-1]->SetName (titlechi2);
  plotCanvas_[plotCanvas_.size()-1]->SetTitle(titlechi2);
  globalchi2->GetYaxis()->SetTitle("global  #chi^{2} = #Sigma_{#splitline{ all}{bins}} #chi^{2}(bin)");
  globalchi2->GetXaxis()->SetRangeUser(minx, maxx);
  globalchi2->Draw("p");
  std::map<TString, double> chi2Result=GetChi2Info(globalchi2, true, verbose, minx, maxx);
  DrawDecayChLabel("e/#mu + Jets Combined");
  DrawCMSLabels(prelim, 0.5*(constLumiMuon+constLumiElec), 0.04, false, false, false);

  std::cout << "mtop(global chi2)= " << chi2Result["min"] << "+" << chi2Result["up"] << "-" << chi2Result["dn"] << std::endl;
  
  // linear fit to extracted mass values
  addCanvas(plotCanvas_);
  plotCanvas_[plotCanvas_.size()-1]->cd(0);
  TF1* combfit=new TF1("combfitdata", "[0]", extractedMass->GetBinLowEdge(1), extractedMass->GetBinLowEdge(extractedMass->GetNbinsX()+1));
  combfit->SetLineColor(kBlue);
  combfit->SetLineWidth(3);
  TString option="Q";
  if(verbose>0)option="";
  if(verbose>1)option="V";
  option+="0";
  extractedMass->Fit(combfit, option,"",extractedMass->GetBinLowEdge(1), extractedMass->GetBinLowEdge(extractedMass->GetNbinsX()+1));
  double mcomb=combfit->GetParameter(0);
  double merr =combfit->GetParError(0);
  // histogram with global chi2 result
  TH1F* CombMass= new TH1F("CombMass", "CombMass", 1, binning_["rhos"][0], binning_["rhos"][binning_["rhos"].size()-1]);
  histogramStyle(*CombMass, kSig, false, 0.1, kBlue);
  CombMass->SetBinContent(1, chi2Result["min"]);
  TGraphAsymmErrors* massfinal= new TGraphAsymmErrors(nbinsx);  
  for(int bin=1; bin<=extractedMass->GetNbinsX(); ++bin){
    double binWidth             = extractedMass->GetBinWidth(bin);
    double binCenter            = extractedMass->GetBinCenter(bin);
    massfinal->SetPoint     (bin-1, binCenter, chi2Result["min"]);
    massfinal->SetPointError(bin-1, 0.5*binWidth, 0.5*binWidth, chi2Result["dn"], chi2Result["up"]);
  }
  massfinal->SetLineColor(kBlue);
  massfinal->SetMarkerColor(kBlue);
  massfinal->SetFillColor(kBlue-9);
  massfinal->SetFillStyle(1001);
  massfinal->SetLineWidth(2);
  massfinal->SetLineStyle(1);
  massfinal->SetMarkerSize(0.1);
  std::cout << "mtop(intersection)=" << mcomb << "+/-" << merr << "GeV (=" << 100*merr/mcomb << "%)" << std::endl;
  // create legend
  TLegend *leg = new TLegend(0.67, 0.55, 0.92, 0.88 );
  legendStyle(*leg,"");
  TH1F* legtev=(TH1F*)TevatronMass->Clone("legtevatron");
  legtev->SetFillColor  (TevatronMassErr->GetFillColor());
  legtev->SetFillStyle  (TevatronMassErr->GetFillStyle());
  TH1F* legcomb=(TH1F*)TevatronMass->Clone("legcomb");
  legcomb->SetLineColor(massfinal->GetLineColor());
  legcomb->SetLineWidth(massfinal->GetLineWidth());
  legcomb->SetFillColor(massfinal->GetFillColor());
  legcomb->SetFillStyle(massfinal->GetFillStyle());
  leg->AddEntry(legtev       , "Tevatron"               , "FLP");
  TH1F* legEmpty=(TH1F*)TevatronMass->Clone("legEmpty");
  legEmpty->SetMarkerColor(kWhite);
  legEmpty->SetMarkerSize(0.000001);
  leg->AddEntry(legEmpty, "#splitline{extracted from}{#rho_{S} in CMS data:}", "P"  );
  leg->AddEntry(extractedMass, "#splitline{#chi^{2} result of}{separate bins}", "LP" );
  //leg->AddEntry(legEmpty, " ", "P" );
  leg->AddEntry(legcomb      , "#splitline{global #chi^{2} result:}{"+getTStringFromDouble(chi2Result["min"], 1)+"^{+"+getTStringFromDouble(chi2Result["up"], 1)+"}_{- "+getTStringFromDouble(chi2Result["dn"], 1)+"} GeV}", "FLP");

  // draw extracted mass for all bins
  TString title="mtopExtractionFromRhosAllBins";
  if(dataMassDependence) title+="DataMassDependence";
  plotCanvas_[plotCanvas_.size()-1]->SetName (title);
  plotCanvas_[plotCanvas_.size()-1]->SetTitle(title);
  if(!dataMassDependence) plotCanvas_[plotCanvas_.size()-1]->SetLogy(1);
  axesStyle(*extractedMass, xSecLabelName("rhos"), "m^{top} [GeV]");
  histogramStyle(*extractedMass, kData);
  extractedMass->GetYaxis()->SetNoExponent(true);
  extractedMass->GetXaxis()->SetTitleOffset(1.1*extractedMass->GetXaxis()->GetTitleOffset());
  extractedMass->Draw("AXIS");
  massfinal->Draw("e2 same");
  CombMass->Draw("hist same");
  TevatronMassErr->Draw("e2 same");
  TevatronMass->Draw("hist same");
  //extractedMass->Draw("p e same");
  extractedChi2Mass->Draw("p e same");  
  leg->Draw("same");
  extractedMass->Draw("AXIS same");
  DrawCMSLabels(true, luminosity, 0.04, false, false, false);

  // saving 
  if(save){
    if(verbose>0) std::cout << "saving" << std::endl;
    int initWarningLV=gErrorIgnoreLevel;
    if(verbose<1) gErrorIgnoreLevel=kWarning;
    for(unsigned int plot=0; plot<=plotCanvas_.size()-1; ++plot){
      //if(verbose>1) std::cout << " - for " << outputFolder+plotCanvas_[plot]->GetName() << std::endl;
      plotCanvas_[plot]->Print(outputFolder+plotCanvas_[plot]->GetName()+".eps");
      plotCanvas_[plot]->Print(outputFolder+plotCanvas_[plot]->GetName()+".png");
      if(plot<plotCanvas_.size()-1) plotCanvas_[plot]->Print(outputFolder+"massFromRhoS.pdf");
      else plotCanvas_[plot]->Print(outputFolder+"massFromRhoS.pdf)");
      saveToRootFile(outputFile, plotCanvas_[plot], true, verbose-1, "massFromRhos");
    }
    gErrorIgnoreLevel=initWarningLV;
 } 

}

std::map <TString, TH1F*> extraction(int verbose, double luminosity, bool save, double minx, double maxx, double nbinsx, bool dataMassDependence, int binOfInterest, TString outputFolder, TString outputFile){

  // returned result
  std::map<TString, TH1F*> out_;
  // quadratic fit instead of linear for last bin
  bool quad = binOfInterest==4 ? true : false;
  // Draw lines for all intersection points
  bool DrawIntersection=false;
  
  // ============================  
  //  Set ROOT Style
  // ============================
  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  //TGaxis::SetMaxDigits(2);
  //myStyle.cd();
  //gROOT->SetStyle("HHStyle");
  //gROOT->ForceStyle();
  // A) get measurement from data
  // open data result file
  TFile* datafile = TFile::Open(outputFile, "READ");
  // dont associate new objects with file to be able to close it in the end
  gROOT->cd();
  // get canvas for norm rhos
  TCanvas* canvasNom = (TCanvas*)(datafile->Get("finalXSec/rhosNorm")->Clone());
  // get data result plot with final errors from canvas
  TGraphAsymmErrors* dataRaw  = (TGraphAsymmErrors*)((canvasNom->GetPrimitive("dataTotalError"))->Clone("dataTot" ));
  TGraphAsymmErrors* dataStat = (TGraphAsymmErrors*)((canvasNom->GetPrimitive("dataStatError" ))->Clone("dataStat"));
  // binning
  std::map< TString, std::vector<double> > binning_ = makeVariableBinning();
  int Nbins = std::abs(binning_["rhos"][binning_["rhos"].size()-1]-binning_["rhos"][0])*100;
  // get histogram in correct binning
  TH1F* data= new TH1F("dataRhoSNorm"+getTStringFromInt(binOfInterest), "dataRhoSNorm"+getTStringFromInt(binOfInterest), Nbins, binning_["rhos"][0], binning_["rhos"][binning_["rhos"].size()-1]);
  reBinTH1F(*data, binning_["rhos"], 0);
  TH1F* dataStati=(TH1F*)data->Clone("dataStati");
  //TH1F* dataStat=(TH1F*)(data->Clone("dataStat"));
  // refill TGraphAsymmErrors to rebinned histo
  for(int bin=1; bin<=data->GetNbinsX(); ++bin){
    //std::cout << dataRaw ->GetY()[bin] << std::endl;
    //std::cout << dataStat->GetY()[bin] << std::endl;
    data     ->SetBinContent(bin, dataRaw ->GetY()[bin]);
    dataStati->SetBinContent(bin, dataStat->GetY()[bin]);
    double err    =dataRaw ->GetErrorYhigh(bin);
    double errStat=dataStat->GetErrorYhigh(bin);
    if(err<dataRaw->GetErrorYlow(bin)) err=dataRaw->GetErrorYlow(bin);
    data     ->SetBinError(bin, err    );
    dataStati->SetBinError(bin, errStat);
  }
  axesStyle(*data, xSecLabelName("rhos"), "#frac{1}{#sigma} #frac{d#sigma}{d#rho_{S}}" , 0., 4.);
  data->SetTitleOffset(0.8*data->GetTitleOffset("Y"), "Y");
  std::cout << "bin: [" << data->GetBinLowEdge(binOfInterest) << ".." << data->GetBinLowEdge(binOfInterest+1) << "]" <<  std::endl; 
  // get values for different top masses
  std::vector<double > MassPointsData_;
  std::vector<TH1F*  > data_;
  if(dataMassDependence){
    MassPointsData_.push_back(161.5);
    MassPointsData_.push_back(163.5);
    MassPointsData_.push_back(166.5);
    MassPointsData_.push_back(169.5);
  }
  MassPointsData_.push_back(172.5);
  if(dataMassDependence){
    MassPointsData_.push_back(175.5);
    MassPointsData_.push_back(178.5);
    MassPointsData_.push_back(181.5);
    MassPointsData_.push_back(184.5);
  }
  // get result canvas for measurement from all masspoints
  std::vector<TCanvas* > canvasData_;
  canvasData_.push_back( (TCanvas*)(datafile->Get("xSec/sysTopMassDown4/rhosNorm")->Clone("dataCanvDn4")) );
  canvasData_.push_back( (TCanvas*)(datafile->Get("xSec/sysTopMassDown3/rhosNorm")->Clone("dataCanvDn3")) );
  canvasData_.push_back( (TCanvas*)(datafile->Get("xSec/sysTopMassDown2/rhosNorm")->Clone("dataCanvDn2")) );
  canvasData_.push_back( (TCanvas*)(datafile->Get("xSec/sysTopMassDown/rhosNorm" )->Clone("dataCanvDn" )) );
  canvasData_.push_back( (TCanvas*)(                                     canvasNom->Clone("dataCanvNom")) );
  canvasData_.push_back( (TCanvas*)(datafile->Get("xSec/sysTopMassUp/rhosNorm"   )->Clone("dataCanvUp" )) );
  canvasData_.push_back( (TCanvas*)(datafile->Get("xSec/sysTopMassUp2/rhosNorm"  )->Clone("dataCanvUp2")) );
  canvasData_.push_back( (TCanvas*)(datafile->Get("xSec/sysTopMassUp3/rhosNorm"  )->Clone("dataCanvUp3")) );
  canvasData_.push_back( (TCanvas*)(datafile->Get("xSec/sysTopMassUp4/rhosNorm"  )->Clone("dataCanvUp4")) );
  // get rhos cross section measurement curve from all masspoints
  TH1F* dataDn4= (TH1F*)(canvasData_[0]->GetPrimitive("rhoskData")->Clone("dataDn4"));
  TH1F* dataDn3= (TH1F*)(canvasData_[1]->GetPrimitive("rhoskData")->Clone("dataDn3"));
  TH1F* dataDn2= (TH1F*)(canvasData_[2]->GetPrimitive("rhoskData")->Clone("dataDn2"));
  TH1F* dataDn = (TH1F*)(canvasData_[3]->GetPrimitive("rhoskData")->Clone("dataDn" ));
  TH1F* dataUp = (TH1F*)(canvasData_[5]->GetPrimitive("rhoskData")->Clone("dataUp" ));
  TH1F* dataUp2= (TH1F*)(canvasData_[6]->GetPrimitive("rhoskData")->Clone("dataUp2"));
  TH1F* dataUp3= (TH1F*)(canvasData_[7]->GetPrimitive("rhoskData")->Clone("dataUp3"));
  TH1F* dataUp4= (TH1F*)(canvasData_[8]->GetPrimitive("rhoskData")->Clone("dataUp4"));
  // collect them in a vector
  if(dataMassDependence){
    data_.push_back((TH1F*)dataDn4->Clone("measurementDn4"    ));
    data_.push_back((TH1F*)dataDn3->Clone("measurementDn3"    ));
    data_.push_back((TH1F*)dataDn2->Clone("measurementDn2"    ));
    data_.push_back((TH1F*)dataDn ->Clone("measurementDn"     ));
  }
  int centralMP=data_.size();
  data_.push_back((TH1F*)data  ->Clone("measurementCentral"));
  if(dataMassDependence){
    data_.push_back((TH1F*)dataUp ->Clone("measurementUp"     ));
    data_.push_back((TH1F*)dataUp2->Clone("measurementUp2"    ));
    data_.push_back((TH1F*)dataUp3->Clone("measurementUp3"    ));
    data_.push_back((TH1F*)dataUp4->Clone("measurementUp4"    ));
  }
  // extract total uncertainty for bin of interest from central mass point 
  double central =data->GetBinContent(binOfInterest);
  double errorAbs=data->GetBinError  (binOfInterest);
  double errRel=std::abs(1.-(central+errorAbs)/central);
  if(verbose>0) std::cout << "central data result: " << central << "+/-" << errorAbs << " (=" << errRel*100 << "% unc.)" << std::endl;
  
  // extract statistical uncertainty for bin of interest from all mass points 
  std::map <double, double > statError_; // (masspoint, stat Err)
  for( unsigned int masspoint=0; masspoint<data_.size(); ++masspoint ){
    // get plot with stat. errors:
    TH1F* tempStat= (masspoint==centralMP ? dataStati : data_[masspoint]);    
    // save value in vector as pair (masspoint, stat Err)
    statError_[MassPointsData_[masspoint]]=tempStat->GetBinError(binOfInterest);
    if(verbose>1){
      std::cout << "bin: " << binOfInterest << std::endl;
      std::cout << MassPointsData_[masspoint] << ": " << tempStat->GetBinContent(binOfInterest) << "+/-" << tempStat->GetBinError(binOfInterest) << std::endl;
    }
  }

  // collect rhos data results from bin of interest for different masses in one separate plot
  TH1F* measurement= new TH1F("CMSdata"+getTStringFromInt(binOfInterest), "CMSdata"+getTStringFromInt(binOfInterest), nbinsx, minx, maxx);
  histogramStyle(*measurement, kData, false, 1.5, color_[0]);
  measurement->SetMarkerColor(color_[0]);
  measurement->SetMarkerStyle(33);
  TH1F* tempMeasStat=(TH1F*)measurement->Clone("tempMeasStat");
  // loop samples
  for( unsigned int masspoint=0; masspoint<data_.size(); ++masspoint ){
    double value=data_[masspoint]->GetBinContent(binOfInterest);
    if(verbose>0) std::cout << "data(" << MassPointsData_[masspoint] << "GeV): " << value << std::endl;
    measurement ->SetBinContent(measurement->FindBin(MassPointsData_[masspoint]), value   );
    measurement ->SetBinError  (measurement->FindBin(MassPointsData_[masspoint]), errorAbs);
    tempMeasStat->SetBinContent(measurement->FindBin(MassPointsData_[masspoint]), value   );
    tempMeasStat->SetBinError  (measurement->FindBin(MassPointsData_[masspoint]), statError_[MassPointsData_[masspoint]]);
  }

  // keep modified versions for drawing
  TH1F* measNomOnly   =((TH1F*)(measurement->Clone("measurementNominalOnly")));
  TH1F* measNonNomOnly=((TH1F*)(measurement->Clone("measurementNonNominal" )));
  TGraphAsymmErrors* measStatUncOnly = new TGraphAsymmErrors((TH1F*)tempMeasStat->Clone("temp"));
  measStatUncOnly->SetName ("measurementStatAll");
  measStatUncOnly->SetTitle(measStatUncOnly->GetName());
  measStatUncOnly->SetMarkerSize(0.01);
  
  for( unsigned int masspoint=0; masspoint<data_.size(); ++masspoint ){
    double value=data_[masspoint]->GetBinContent(binOfInterest);
    // measNomOnly: rm all but cwentral point
    if(value!=central) measNomOnly->SetBinContent(measNomOnly->FindBin(MassPointsData_[masspoint]), 0.);
    // measNonNomOnly: rm all errors
    measNonNomOnly->SetBinError(measNonNomOnly->FindBin(MassPointsData_[masspoint]), 0.); // FIXME
  }
  // do a linear or quadratic fit of the mass dependence
  TF1* lindata= quad ? new TF1("quadFitdata", "[2]*x*x+[0]*x+[1]", minx, maxx) : new TF1("linFitdata", "[0]*x+[1]", minx, maxx);
  
  lindata->SetLineColor(measurement->GetLineColor());
  lindata->SetLineWidth(2);
  TString option="Q";
  if(verbose>0)option="";
  if(verbose>1)option="V";
  option+="R";
  measStatUncOnly->Fit(lindata, option,"",MassPointsData_[0]-0.5, MassPointsData_[MassPointsData_.size()-1]+0.5);
  if(!dataMassDependence){
    lindata->SetParameter(0, 0.     );
    lindata->SetParameter(1, central);
    if(quad) lindata->SetParameter(2, 0.);
  }
  double adata=lindata->GetParameter(0);
  double bdata=lindata->GetParameter(1);
  double cdata= quad ? lindata->GetParameter(2) : 0.;
//   TF1* lindataUp=new TF1("linFitdataUp", "[0]*x+[1]", minx, maxx);
//   TF1* lindataDn=new TF1("linFitdataDn", "[0]*x+[1]", minx, maxx);
//   lindataUp->SetParameter(0, adata);
//   lindataDn->SetParameter(0, adata);
//   lindataUp->SetParameter(1, bdata+errorAbs);
//   lindataDn->SetParameter(1, bdata-errorAbs);

  // create data band plot
  TH1F* measurementCentral= new TH1F  ("CMSdataNom"+getTStringFromInt(binOfInterest), "CMSdataNom"+getTStringFromInt(binOfInterest)    , nbinsx, minx, maxx);
  measurementCentral->Add(lindata);
  histogramStyle(*measurementCentral, kSig, false, 0.1);
  TH1F* measurementUp= new TH1F("CMSdataUp"+getTStringFromInt(binOfInterest), "CMSdataUp"+getTStringFromInt(binOfInterest), nbinsx, minx, maxx);
  //measurement->Add(lindataUp);
  histogramStyle(*measurementUp, kBkg, true, 0.1); 
  measurementUp->SetFillColor(color_[1]);
  measurementUp->SetFillStyle(1001);
  measurementUp->SetLineColor(10);
  TH1F* measurementDn= new TH1F("CMSdataDn"+getTStringFromInt(binOfInterest), "CMSdataDn"+getTStringFromInt(binOfInterest), nbinsx, minx, maxx);
  //measurementDn->Add(lindataDn);
  histogramStyle(*measurementDn, kBkg, true, 0.1, 10); 
  for(int bin=1; bin<=measurementCentral->GetNbinsX(); ++bin){
    measurementUp->SetBinContent(bin, measurementCentral->GetBinContent(bin)+errorAbs);
    measurementDn->SetBinContent(bin, measurementCentral->GetBinContent(bin)-errorAbs);
  }

  // B) get prediction with top mass dependence
  //  load rootfiles
  std::vector<TFile* > file_;
  file_.push_back(TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/TopMassDown/combinedDiffXSecSigTopMassDown4Summer12PF.root", "Open"));
  file_.push_back(TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/TopMassDown/combinedDiffXSecSigTopMassDown3Summer12PF.root", "Open"));
  file_.push_back(TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/TopMassDown/combinedDiffXSecSigTopMassDown2Summer12PF.root", "Open"));
  file_.push_back(TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/TopMassDown/combinedDiffXSecSigTopMassDownSummer12PF.root" , "Open"));
  file_.push_back(TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/combinedDiffXSecSigSummer12PF.root"                        , "Open"));
  file_.push_back(TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/TopMassUp/combinedDiffXSecSigTopMassUpSummer12PF.root"     , "Open"));
  file_.push_back(TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/TopMassUp/combinedDiffXSecSigTopMassUp2Summer12PF.root"    , "Open"));
  file_.push_back(TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/TopMassUp/combinedDiffXSecSigTopMassUp3Summer12PF.root"    , "Open"));
  file_.push_back(TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/TopMassUp/combinedDiffXSecSigTopMassUp4Summer12PF.root"    , "Open"));
  std::vector<double > mtop_;
  // fill mass indication
  mtop_.push_back(161.5);
  mtop_.push_back(163.5);
  mtop_.push_back(166.5);
  mtop_.push_back(169.5);
  mtop_.push_back(172.5);
  mtop_.push_back(175.5);
  mtop_.push_back(178.5);
  mtop_.push_back(181.5);
  mtop_.push_back(184.5);

  // C) get prediction with theory parameter dependence
  std::vector<TFile* > fileSys_;
  //  load rootfiles
  fileSys_.push_back(TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/ScaleUp/combinedDiffXSecSigScaleUpSummer12PF.root"     , "Open"));
  fileSys_.push_back(TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/ScaleDown/combinedDiffXSecSigScaleDownSummer12PF.root" , "Open"));
  fileSys_.push_back(TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/MatchUp/combinedDiffXSecSigMatchUpSummer12PF.root"     , "Open"));
  fileSys_.push_back(TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/MatchDown/combinedDiffXSecSigMatchDownSummer12PF.root" , "Open"));

  // get rhos(mtop) plots
  std::vector<TH1F* > rhosMC_;
  // loop samples
  for(int sample=0; sample<(int)file_.size(); ++sample ){
    // get plot
    TH1F* temp=(TH1F*)((file_[sample]->Get("compositedHadronGenPhaseSpace/rhosGen"))->Clone("temp"));
    // rebin it
    reBinTH1F(*temp, binning_["rhos"], 0);
    // divide by integral and binwidth
    double N=temp->Integral(1, temp->GetNbinsX());
    for(int i=1; i<= temp->GetNbinsX(); i++){
      temp->SetBinContent(i, ( (double)(temp->GetBinContent(i)) / (N*temp->GetBinWidth(i)) ) );
    }
    // keep it in vector
    rhosMC_.push_back((TH1F*)temp->Clone(TString("rhosMC")+getTStringFromInt(sample)));  
  }

  // get rhos(parameter) plots
  std::vector<TH1F* > rhosMCsys_;
  // loop samples
  for(int sample=0; sample<(int)fileSys_.size(); ++sample ){
    // get plot
    TH1F* temp=(TH1F*)((fileSys_[sample]->Get("compositedHadronGenPhaseSpace/rhosGen"))->Clone("tempSys"));
    // rebin it
    reBinTH1F(*temp, binning_["rhos"], 0);
    // divide by integral and binwidth
    double N=temp->Integral(1, temp->GetNbinsX());
    for(int i=1; i<= temp->GetNbinsX(); i++){
      temp->SetBinContent(i, ( (double)(temp->GetBinContent(i)) / (N*temp->GetBinWidth(i)) ) );
    }
    // keep it in vector
    rhosMCsys_.push_back((TH1F*)temp->Clone(TString("rhosMCsys")+getTStringFromInt(sample)));  
  }

  // collect rhos results from bin of interest for different masses is separate plot
  double MCnom=0.;
  TH1F* MC= new TH1F("CMSMC", "CMSMC", nbinsx, minx, maxx);
  //reBinTH1F(*data, mtop_,0);
  // loop samples
  for( unsigned int sample=0; sample<file_.size(); ++sample ){
    double value=rhosMC_[sample]->GetBinContent(binOfInterest);
    if(verbose>0) std::cout << "MC(" << mtop_[sample] << "GeV): " << value << std::endl;
    MC->SetBinContent(MC->FindBin(mtop_[sample]), value);
    if(mtop_[sample]==172.5) MCnom=value;
  }
  histogramStyle(*MC, kData, false, 1.2, kGreen-2);
  

  // extract values for parameter dependence
  double MCnomUnc=0.;
  for( unsigned int sample=0; sample<fileSys_.size(); sample=sample+2){
    double valueUp=rhosMCsys_[sample  ]->GetBinContent(binOfInterest);
    double valueDn=rhosMCsys_[sample+1]->GetBinContent(binOfInterest);
    double unc=0.5*(std::abs(valueUp-MCnom)+std::abs(valueDn-MCnom));
    MCnomUnc+=unc;
    //std::cout << "c : " << MCnom << std::endl;
    //std::cout << "up=" << fileSys_[sample  ]->GetName() << ": " << valueUp << std::endl;
    //std::cout << "dn=" << fileSys_[sample+1]->GetName() << ": " << valueDn << std::endl;
  }
  
  // do a linear fit for the mass dependence of the simulation 
  TF1* lin=new TF1("linFit", "[0]*x+[1]", minx, maxx);
  lin->SetLineColor(MC->GetLineColor());
  lin->SetLineWidth(2);
  MC->Fit(lin,option,"",mtop_[0]-0.5, mtop_[mtop_.size()-1]+0.5);
  double a=lin->GetParameter(0);
  double b=lin->GetParameter(1);
  // C) calculate intersection
  // c1) linear
  double mresultfit=(b-bdata)/(adata-a);
  double centralfit=adata*mresultfit+bdata;
  double mup=(b-MCnomUnc-(bdata+errorAbs))/(adata-a);
  double mdn=(b+MCnomUnc-(bdata-errorAbs))/(adata-a);
  double up=adata*mup+(bdata+errorAbs);
  double dn=adata*mdn+(bdata-errorAbs);
  // c2 quadratic
  if(quad){
    if(verbose>1){
      std::cout << "quadratic intersection" << std::endl;
      std::cout << "a=" << a << std::endl;
      std::cout << "b=" << b << std::endl;
      std::cout << "adata=" << adata << std::endl;
      std::cout << "bdata=" << bdata << std::endl;
      std::cout << "cdata=" << cdata << std::endl;
    }
    // c2.1 central value
    double p=(a-adata)/(2*cdata); // -p/2 from p-q
    double q=(bdata-b)/cdata;     // q from p-q
    if(verbose>1){
      std::cout << "-p/2=" << p << std::endl;
      std::cout << "q="    << q << std::endl;
    }
    // solutions from p-q
    double sol1=p+sqrt(p*p-q);
    double sol2=p-sqrt(p*p-q);
    if(verbose>1){
      std::cout << "sol1=" << sol1 << std::endl;
      std::cout << "sol2=" << sol2 << std::endl;
    }
    // take smaller solution for central value
    mresultfit=sol1;
    if(sol1<0.||(sol2>=0.&&sol2<sol1)) mresultfit=sol2;
    if(verbose>1) std::cout << "mresultfit=" << mresultfit << std::endl;
    centralfit=adata*mresultfit+bdata+cdata*mresultfit*mresultfit;
    if(verbose>1) std::cout << "centralfit=" << centralfit << std::endl;
    // c2.2 up (b->b-MCnommUnc, bdata->bdata+errorAbs)
    double qup=(bdata+errorAbs-(b-MCnomUnc))/cdata; // q from p-q
    double sol1up=p+sqrt(p*p-qup);
    double sol2up=p-sqrt(p*p-qup);
    mup=sol1up;
    if(sol1up<0.||(sol2up>=0.&&sol2up<sol1up)) mup=sol2up;
    up=adata*mup+(bdata+errorAbs)+cdata*mup*mup;
    // c2.3 dn (b->b+MCnommUnc, bdata->bdata-errorAbs)
    double qdn=(bdata-errorAbs-(b+MCnomUnc))/cdata; // q from p-q
    double sol1dn=p+sqrt(p*p-qdn);
    double sol2dn=p-sqrt(p*p-qdn);
    mdn=sol1dn;
    if(sol1dn<0.||(sol2dn>=0.&&sol2dn<sol1dn)) mdn=sol2dn;
    dn=adata*mdn+(bdata-errorAbs)+cdata*mdn*mdn;
  }
  if(verbose>0){
    std::cout << "mtop= (" << mdn << ".." << mresultfit << ".." << mup << ")" << std::endl; 
    if(verbose>0) std::cout << "mtop= " << mresultfit << " +" << mresultfit-mup << " -" << mdn-mresultfit << std::endl; 
  }
  TGraphAsymmErrors* MCunc= new TGraphAsymmErrors(nbinsx);  
  for(int bin=1; bin<=measurementCentral->GetNbinsX(); ++bin){
    //double binEdgeDown          = measurementCentral->GetBinLowEdge(bin);
    //double binEdgeUp            = measurementCentral->GetBinLowEdge(bin+1);
    double binWidth             = measurementCentral->GetBinWidth(bin);
    double binCenter            = measurementCentral->GetBinCenter(bin);
    double value=a*binCenter+b;
    MCunc->SetPoint     (bin-1, binCenter, value);
    MCunc->SetPointError(bin-1, 0.5*binWidth, 0.5*binWidth, MCnomUnc, MCnomUnc);
  }
  MCunc->SetLineColor(kGreen);
  MCunc->SetMarkerColor(kGreen);
  MCunc->SetFillColor(kGreen);
  MCunc->SetFillStyle(3344);
  MCunc->SetLineWidth(1);
  MCunc->SetLineStyle(1);
  MCunc->SetMarkerSize(0.1);
  
  // chi2 approach for binOfInterest
  double minxchi2=120;
  double maxxchi2=250;
  TH1F* chi2CMS= new TH1F  ("chi2Bin"+getTStringFromInt(binOfInterest), "chi2Bin"+getTStringFromInt(binOfInterest), 1300*(maxxchi2-minxchi2), minxchi2, maxxchi2);
  histogramStyle(*chi2CMS, kData, false, 1.3);
  axesStyle(*chi2CMS, "m^{top}", "#chi^{2} ("+getTStringFromDouble(data->GetBinLowEdge(binOfInterest))+" #leq#rho_{S}#leq "+getTStringFromDouble(data->GetBinLowEdge(binOfInterest+1))+")");
  chi2CMS->GetXaxis()->SetNoExponent(true);
  // quadratic combination of uncertainty of MC pred.+ data measurement from central masspoint
  double totUnc=sqrt(MCnomUnc*MCnomUnc+errorAbs*errorAbs);
  // loop all chi2 bins
  for(int masspoint=1; masspoint<chi2CMS->GetNbinsX(); ++masspoint ){
    // mass value
    double m=chi2CMS->GetBinCenter(masspoint);
    // data and prediction from 
    double valuedata=adata*m+bdata+m*m*cdata;
    double valueMC  =a*m+b;
    // chi2
    double chi2=(std::abs(valuedata-valueMC)*std::abs(valuedata-valueMC))/(totUnc*totUnc);
    chi2CMS->SetBinContent(masspoint, chi2);
  }
  out_["chi2"   ]=(TH1F*)chi2CMS->Clone(TString("chi2Bin")+getTStringFromInt(binOfInterest));
  chi2CMS->GetXaxis()->SetRangeUser(minx, maxx);
  // quadratic fit
  //TF1* quadchi= new TF1("quadchi", "[2]*x*x+[0]*x+[1]", minx, maxx);
  //quadchi->SetLineColor(kBlack);
  //quadchi->SetLineWidth(2);
  //TString optionchi="Q";
  //if(verbose>0)optionchi="";
  //if(verbose>1)optionchi="V";
  //optionchi+="R";
  //chi2CMS->Fit(quadchi, optionchi,"",MassPointsData_[0]-0.5, MassPointsData_[MassPointsData_.size()-1]+0.5);

  // ---
  //    legend
  // ---
  TLegend *leg = new TLegend(0.61, (binOfInterest!=3) ? 0.7 : 0.25, 0.96, (binOfInterest!=3) ? 0.86 : 0.41);
  legendStyle(*leg,"");
  TH1F* legmeasurement=(TH1F*)measurementCentral->Clone("legdata");
  legmeasurement->SetFillColor  (measurementUp->GetFillColor());
  legmeasurement->SetFillStyle  (measurementUp->GetFillStyle());
  legmeasurement->SetMarkerStyle(measurement->GetMarkerStyle());
  legmeasurement->SetMarkerColor(measurement->GetMarkerColor());
  legmeasurement->SetMarkerSize (measurement->GetMarkerSize() );
  TH1F* legMC=(TH1F*)MC->Clone("legMC");
  legMC->SetFillColor(MCunc->GetFillColor());
  legMC->SetFillStyle(MCunc->GetFillStyle());
  leg->AddEntry(legmeasurement, "Data"             , "PL" );
  leg->AddEntry(legmeasurement, "Fit to Data"      , "FL" );
  leg->AddEntry(legMC         , "Fit to Simulation", "PLF");
  //leg->AddEntry(lin           , "linear fit"    , "L" );
  
  // D) create final plots
  std::vector<TCanvas*> plotCanvas_;
  // intersection method
  addCanvas(plotCanvas_);
  plotCanvas_[plotCanvas_.size()-1]->cd(0);
  TString title="mtopExtractionFromRhosBin"+getTStringFromInt(binOfInterest);
  if(dataMassDependence) title+="DataMassDependence";
  plotCanvas_[plotCanvas_.size()-1]->SetName (title);
  plotCanvas_[plotCanvas_.size()-1]->SetTitle(title);
  double min=0.2;
  double max=0.8;
  if(binOfInterest==1){
    min=0.0005;
    max=0.06;
  }
  else if(binOfInterest==2){
    min=1.0;
    max=2.2;
  }
  else if(binOfInterest==3){
    min=1.8;
    max=3.4;
  }
  else if(binOfInterest==4){
    min=0.2;
    max=0.8;
  }
  MC->SetMaximum(max);
  MC->SetMinimum(min);
  MC->GetXaxis()->SetNoExponent(true);
  axesStyle(*MC, "m^{top}", "#frac{1}{#sigma} #frac{d#sigma}{d#rho_{S}} ("+getTStringFromDouble(data->GetBinLowEdge(binOfInterest))+" #leq#rho_{S}#leq "+getTStringFromDouble(data->GetBinLowEdge(binOfInterest+1))+")");
  MC->Draw("AXIS");
  measurementUp     ->Draw("hist same");
  measurementDn     ->Draw("hist same");
  measurementCentral->Draw("hist same");
  MCunc->Draw("p e2 same");
  MC   ->Draw("p same");
  lin->Draw("same");
  measNonNomOnly->Draw("p same");
  measNomOnly->Draw("p same");
  measStatUncOnly->Draw("e same");
  MC->Draw("AXIS same");
  if(DrawIntersection){
    if(mresultfit<=maxx&&mresultfit>=minx) drawLine(mresultfit, min, mresultfit, centralfit, kBlack, 2, 2);
    if(mup<=maxx       &&mup>=minx       ) drawLine(mup       , min, mup       , up        , kBlack, 1, 2);
    if(mdn<=maxx       &&mdn>=minx       ) drawLine(mdn       , min, mdn       , dn        , kBlack, 1, 2);
  }
  leg->Draw("same");
  DrawCMSLabels(true, luminosity, 0.04, false, false, false);
  // chi2 method
  addCanvas(plotCanvas_);
  plotCanvas_[plotCanvas_.size()-1]->cd(0);
  TString title2="mtopChi2ScanRhosBin"+getTStringFromInt(binOfInterest);
  if(dataMassDependence) title2+="DataMassDependence";
  plotCanvas_[plotCanvas_.size()-1]->SetName (title2);
  plotCanvas_[plotCanvas_.size()-1]->SetTitle(title2);
  chi2CMS->Draw("p");
  std::map<TString, double> chi2Result=GetChi2Info(chi2CMS, true, verbose, minx, maxx);  
  //if(binOfInterest<=2) quadchi->Draw("hist same");
  DrawDecayChLabel("e/#mu + Jets Combined");
  DrawCMSLabels(prelim, 0.5*(constLumiMuon+constLumiElec), 0.04, false, false, false);

  // E) Canvas with different mass predictions
  //    duplicate canvas from nominal data measurement
  if(binOfInterest==1){
    gROOT->cd();
    std::cout << "create cross section for different masses plot" << std::endl;
    //plotCanvas_.push_back((TCanvas*)(canvasNom->Clone()));
    addCanvas(plotCanvas_);  
    plotCanvas_[plotCanvas_.size()-1]->cd(0);
    plotCanvas_[plotCanvas_.size()-1]->SetName ("finalXSecrhosNormCombinedHadronLvPSMassCurves");
    plotCanvas_[plotCanvas_.size()-1]->SetTitle("finalXSecrhosNormCombinedHadronLvPSMassCurves");
    // remove other MC curves
    data->Draw("AXIS");
    // style for mass samples
    std::vector <int > colorMC_;
    colorMC_.push_back(6         );
    colorMC_.push_back(kBlue     );
    colorMC_.push_back(kGreen+1  );
    colorMC_.push_back(kGreen+3  );
    colorMC_.push_back(kRed      );
    colorMC_.push_back(kAzure+8  );
    colorMC_.push_back(kOrange+3 ); 
    colorMC_.push_back(kMagenta+2); 
    colorMC_.push_back(kOrange+8 );
    std::vector <int > style_;
    style_.push_back(2);
    style_.push_back(10);
    style_.push_back(4);
    style_.push_back(6);
    style_.push_back(1);
    style_.push_back(7);
    style_.push_back(8);
    style_.push_back(9);
    style_.push_back(3);
    // create legend
    TLegend *legMasses = new TLegend(0.7, 0.4 , 1.2, 0.8 );
    TLegend *legDat    = new TLegend(0.6, 0.75, 1.2, 0.95);
    legDat->AddEntry(data, "#splitline{Data#scale[0.8]{, unfolded}}{#scale[0.8]{with m_{top}=172.5GeV}}", "LP");
    legendStyle(*legMasses,"#scale[0.9]{#splitline{MadGraph+Pythia}{with m_{top}[GeV]=}}");
    legendStyle(*legDat   , "");
    // loop mass points
    for(int sample=0; sample<(int)file_.size(); ++sample ){
      // get and draw curves
      rhosMC_[sample]->SetLineColor(  colorMC_[sample]);
      rhosMC_[sample]->SetMarkerColor(colorMC_[sample]);
      rhosMC_[sample]->SetLineWidth(3);
      rhosMC_[sample]->SetLineStyle(style_[sample]);
      rhosMC_[sample]->Draw("hist same");
      // add legend entry
      double val=mtop_[sample];
      legMasses->AddEntry(rhosMC_[sample], Form("%2.1f", val), "L" );
    }
    // redraw data
    dataStat->Draw("e p same");
    data    ->Draw("e p same");
    // draw legend
    legMasses->Draw("same");
    legDat   ->Draw("same");
    // draw labels
    DrawDecayChLabel("e/#mu + Jets Combined");
    DrawCMSLabels(prelim, 0.5*(constLumiMuon+constLumiElec), 0.04, false, false, false);
    // saving
    if(save){
      TCanvas* tempCanv=(TCanvas*)(plotCanvas_[plotCanvas_.size()-1]->Clone());
      int initWarningLV=gErrorIgnoreLevel;
      if(verbose<1) gErrorIgnoreLevel=kWarning;   
      tempCanv->Print(outputFolder+tempCanv->GetName()+".eps");
      tempCanv->Print(outputFolder+tempCanv->GetName()+".png");
      tempCanv->Print(outputFolder+"massFromRhoS.pdf(");
      saveToRootFile(outputFile, tempCanv, true, verbose-1, "massFromRhos");
      gErrorIgnoreLevel=initWarningLV;
    }
    // canvas with ratio
    gROOT->cd();
    bool drawRat =true;
    bool drawRat2=false;
    std::cout << "create ratios for cross section for different masses" << std::endl;
    if(drawRat){
      // adapt some stylings
      legMasses->SetY1(1.12*legMasses->GetY1());
      data->SetLabelSize(0.);
      data->SetTitleSize(0.);
      data->SetMaximum(5.);
      rhosMC_[0]->SetLabelSize(0.85*rhosMC_[0]->GetLabelSize());
      rhosMC_[0]->SetTitleSize(1.3*rhosMC_[0]->GetTitleSize());
      rhosMC_[0]->SetLabelSize(1.3*rhosMC_[0]->GetLabelSize("y"), "Y");
      rhosMC_[0]->SetTitleSize(1.4*rhosMC_[0]->GetTitleSize("y"), "Y");

      plotCanvas_[plotCanvas_.size()-1]->Update();
      // ratio function
      plotCanvas_.push_back(drawFinalResultRatio((TH1F*)data->Clone("dataTotRat"), 0.35, 1.75, myStyle, verbose, rhosMC_, (TCanvas*)(plotCanvas_[plotCanvas_.size()-1]->Clone()), -1., -1., dataStat, true));
      plotCanvas_[plotCanvas_.size()-1]->cd();
      plotCanvas_[plotCanvas_.size()-1]->SetName ("ratioFinalXSecrhosNormCombinedHadronPhaseSpaceMassCurves");
      plotCanvas_[plotCanvas_.size()-1]->SetTitle("ratioFinalXSecrhosNormCombinedHadronPhaseSpaceMassCurves");
      plotCanvas_[plotCanvas_.size()-1]->Draw();
      plotCanvas_[plotCanvas_.size()-1]->Update();
      if(save){    
	TCanvas* tempCanv=(TCanvas*)(plotCanvas_[plotCanvas_.size()-1]->Clone());
	int initWarningLV=gErrorIgnoreLevel;
	if(verbose<1) gErrorIgnoreLevel=kWarning;   
	tempCanv->Print(outputFolder+tempCanv->GetName()+".eps");
	tempCanv->Print(outputFolder+tempCanv->GetName()+".png");
	tempCanv->Print(outputFolder+"massFromRhoS.pdf");
	saveToRootFile(outputFile, tempCanv, true, verbose-1, "massFromRhos");
	gErrorIgnoreLevel=initWarningLV;
      }
    }
    if(drawRat2){
     gPad->cd();
      //plotCanvas_[plotCanvas_.size()-1]->Draw();
      //plotCanvas_[plotCanvas_.size()-1]->Update();
      //data->Draw("hist same");
      for(int bin=1; bin<=data->GetNbinsX(); ++bin) std::cout << data->GetBinContent(bin) << std::endl;
      // error vector dummy
      std::vector<double> zeroerr_;
      for(int bin=1; bin<=rhosMC_[0]->GetNbinsX(); ++bin) zeroerr_.push_back(0.);
      std::cout << "- 1" << std::endl;
      // axis
      drawRatio(data, data, 0.6, 1.6, myStyle, 0, zeroerr_, "prediction", "data", "axis", kBlack, false);
      std::cout << "- 2" << std::endl;
      // MCs
      for(int sample=0; sample<(int)file_.size(); ++sample ){
	drawRatio(rhosMC_[sample], data, 0.6, 1.6, myStyle, 0, zeroerr_, "prediction", "data", "hist same", rhosMC_[sample]->GetLineColor(), false);
      }
      std::cout << "- 3" << std::endl;
      // data error vector
      std::vector<double> err_;
      for(int bin=1; bin<=data->GetNbinsX(); ++bin) err_.push_back(data->GetBinError(bin)/data->GetBinContent(bin));
      // axis
      drawRatio(data, data, 0.6, 1.6, myStyle, 0, err_, "prediction", "data", "e p same", kBlack, false);
    }
  }

  // F) do saving 
  if(save){
    gROOT->cd();
    int initWarningLV=gErrorIgnoreLevel;
    if(verbose<1) gErrorIgnoreLevel=kWarning;
    // loop canvases
    for(unsigned int plot=0; plot<plotCanvas_.size(); ++plot){
      // skip last two canvas for the first bin as different mass predictions have been already saved
      if(binOfInterest!=1||(binOfInterest==1&&plot<plotCanvas_.size()-2)){
	if(verbose>1) std::cout << plot+1 << "/" << plotCanvas_.size()-1 << std::endl;
	if(!plotCanvas_[plot]) std::cout << "ERROR: canvas #" << plot << " not found" << std::endl; 
	else if(verbose>1) std::cout << plotCanvas_[plot]->GetName() << std::endl;
	if(!((TString(plotCanvas_[plot]->GetName())).Contains("final"))){
	  plotCanvas_[plot]->Print(outputFolder+plotCanvas_[plot]->GetName()+".eps");
	  plotCanvas_[plot]->Print(outputFolder+plotCanvas_[plot]->GetName()+".png");
	  plotCanvas_[plot]->Print(outputFolder+"massFromRhoS.pdf");
	  saveToRootFile(outputFile, plotCanvas_[plot], true, verbose-1, "massFromRhos");
	  gROOT->cd();
	}
      }
    }
    gErrorIgnoreLevel=initWarningLV;
  }
  double errIup=std::abs(mresultfit-mup);
  double errIdn=std::abs(mresultfit-mdn);
  double largerUnc= errIup>errIdn ? errIup : errIdn;
  std::cout << "intersection method: " << mresultfit << "+" << errIup<< "-" << errIdn << std::endl;
  std::cout << "chi2 method:         " << chi2Result["min"] << "+" << chi2Result["up"] << "-" << chi2Result["dn"] << std::endl;
  
  // collect output
  TString binext=TString("Bin")+getTStringFromInt(binOfInterest);
  TH1F* massvaltemp=new TH1F("massvaltemp"+binext, "massvaltemp"+binext, 1, 0., 1.);
  massvaltemp->SetBinContent(1, mresultfit);
  TH1F* massunctemp=new TH1F("massunctemp"+binext, "massunctemp"+binext, 1, 0., 1.);
  massunctemp->SetBinContent(1, largerUnc);
  out_["massval"]=massvaltemp;
  out_["massunc"]=massunctemp;
  return out_;
}

std::map<TString, double> GetChi2Info(TH1F* chi2, bool draw, int verbose, double xmin, double xmax){
  // this function scans the provided "chi2", searching for the minimum and also for the uncertainties, coming from Dchi2=1
  // output: the returned map contains 3 entries: "min", "up", "dn" for the minumum, and the deltaChi2+/-1 values
  // draw: do you want to draw lines, indicating the minimum and the values for dChi2=1?
  // verbose: detail of output, 0=no,1=info,2+=degbugging in certain steps
  // xmin/xmax: plotted range of the chi2 function, only needed for draw==true, -1 means first to last bin
 
  // loop bins to find the minimum
  double minval =-1.;
  double minchi2=10000.;
  for(int bin=1; bin<chi2->GetNbinsX(); ++bin){
    double binval =chi2->GetBinCenter(bin);
    double chi2val=chi2->GetBinContent(bin);
    if(verbose>2) std::cout << "chi2(" << binval << ")=" << chi2val << std::endl;
    // check also that minimum is invisible range to avoid second minima somewhere
    if(chi2val<minchi2&&(xmax==-1.||binval<=xmax)&&(xmin==-1||binval>=xmin)){
      if(verbose>2) std::cout << "-> new Minimum!" << std::endl;
      minchi2=chi2val;
      minval=binval;
    }
    // stop at the first minimum
    //else break;
  }
  // check if valid point is found
  if(minval==-1.){
    std::cout << "ERROR: no minimum can be found in GetChi2Info for plot " << chi2->GetTitle() << std::endl;
    exit(0);
  }
  else if(verbose>0) std::cout << "minimum in chi2 found for plot " << chi2->GetTitle() << ": chi2(" << minval << ")=" <<   minchi2 << std::endl;
  // search for deltaChi2=1 bins
  double delChval=minchi2+1.0;
  double upval   =-1.;
  double downval =-1.;
  for(int bin=1; bin<chi2->GetNbinsX()-1; ++bin){
    double chi2val  =chi2->GetBinContent(bin  );
    double chi2valup=chi2->GetBinContent(bin+1);
    // skip bin if the last one was the one for DelChi2=1
    if(downval==-1||(upval==-1.&&bin-1!=chi2->FindBin(downval))||(upval!=-1.&&bin-1!=chi2->FindBin(upval))){
      double tempval=-1.;
      // check if this bin has exactly the searched chi2 value
      if(delChval==chi2val) tempval=chi2->GetBinCenter(bin);
      // check if DelChi2=1 is between this and the next bin
      else if((chi2val<chi2valup&&delChval>chi2val&&delChval<chi2valup)||(chi2valup<chi2val&&delChval<chi2val&&delChval>chi2valup)){
	// take average of those bins, weighted with the difference between their chi2 and the one searched
	tempval= (chi2->GetBinCenter(bin)*std::abs(chi2val-delChval)+chi2->GetBinCenter(bin+1)*std::abs(chi2valup-delChval))/(std::abs(chi2val-delChval)+std::abs(chi2valup-delChval));
      }
      // assign value as down or up
      if(tempval!=-1.){
	if(verbose>1) std::cout << "delchi2=1 for chi2(" << chi2->GetBinCenter(bin) << ")=" << chi2val << ", mtop=" << tempval << std::endl;
      
	if(     downval ==-1.) downval=tempval;
	else if(upval   ==-1.) upval  =tempval;
	else{
	  if(verbose>0) std::cout << "WARNING: more than 2 values are found for deltaChi2=1 in GetChi2Info for plot " << chi2->GetTitle() << ": " << downval << ", "  << upval << ", " << tempval << " (will keep the first two)" << std::endl;
	}
      }
    }
  }
  // check if two valid points have been found
  if(upval==-1.){
    std::cout << "ERROR: no values for deltaChi2=1 wrt. the minimum can be found in GetChi2Info for plot " << chi2->GetTitle() << std::endl;
    exit(0);
  }
  else if(verbose>0) std::cout << "delchi2=1 found for " << downval << " & " << upval << std::endl;
  // do the drawing if wanted
  // expecting a open canvas, where chi2 is drawn already
  if(draw){
    // line at minimum 
    double binF=xmin==-1 ? chi2->GetBinContent(1) : chi2->GetBinContent(chi2->FindBin(xmin));
    double binL=xmax==-1 ? chi2->GetBinContent(chi2->GetNbinsX()-1) : chi2->GetBinContent(chi2->FindBin(xmax));
    double height=binF>binL ? binF : binL;
    drawLine(minval , 0., minval ,  0.25*height                               , kBlack, 2, 2);
    // line at deltaChi2=-1
    drawLine(downval, 0., downval, chi2->GetBinContent(chi2->FindBin(downval)), kBlack, 2, 2);
    // line at deltaChi2=+1
    drawLine(upval  , 0., upval  , chi2->GetBinContent(chi2->FindBin(upval  )), kBlack, 2, 2);
    // line at 1 from left  to deltaChi2=-1  
    if(xmin==-1.) xmin=chi2->GetBinLowEdge(1);
    drawLine(xmin , delChval, downval, delChval, kBlack, 2, 2);
    // line at 1 from right to deltaChi2=+1
    if(xmax==-1.) xmax=chi2->GetBinLowEdge(chi2->GetNbinsX()+1);
    drawLine(upval, delChval, xmax   , delChval, kBlack, 2, 2);
    // draw values as label
    DrawLabel(TString("#splitline{minimum: }{")+getTStringFromDouble(minval)+"^{ + "+getTStringFromDouble(upval -minval)+"}_{  - "+getTStringFromDouble(minval-downval)+"} GeV}", 0.5, 0.4, 0.75, 0.6, 12, 0.03);
  }
  // prepare output 
  std::map <TString, double> out_;
  out_["min"]=minval;
  out_["up" ]=(upval -minval );
  out_["dn" ]=(minval-downval);
  return out_;
}
