/**
*
*   Exectue this macro from the 
*       $CMSSW_BASE/src/TopAnalysis/Configuration/analysis/semiLeptonic/diffXSection
*   via:
*       root -l -b -q sevenToEightTeVComparison.C++g
*
*   Needed input:
*       Dilepton DiffXS 8TeV results in TGraphAssymErrors in a root file:
*           TString line= readLineFromFile(p+4, "/nfs/dust/cms/user/asincruz/CMSSW5314p1_Development/CMSSW_5_3_14_patch1/src/TopAnalysis/Configuration/analysis/diLeptonic/analysis_postCWR/Plots/FinalResults/combined/Hyp"+DilepFileName(name)+"LaTeX.txt");
*       l+jets 8TeV DiffXS results stored on the NAF:
*           TString line= readLineFromFile(p+1, groupSpace+"CommonFiles/topPtInputForReweighting/diffXSecTopSemiLepParton"+name+".txt");
*
*       The 7TeV l+jets results are hardcoded in the function 'seven'
*       The dilepton 7TeV results, theory and BCC are also hardcoded in the function 'dilepton7TeVCurves'
*       The dilepton 8TeV Madgraph+Pythia values are hardcoded in the function 'DilepMCMadGraph'
*       The BCCs for 8TeV for dilepton and l+jets are hardcoded in the function 'BCCvalues'
*
*   If wanted ratio: theory/MadGraph+Pythia.
*       Uncomment the L.535.
*   the function 'drawTheoryLines' will draw the theory/MG+PY ratio curves in all 8TeV ratio plots.
*
*   Needed input: the theory curves as TH1D in single ROOT (per variable)
*       TString filename =  "/nfs/dust/cms/user/asincruz/CMSSW5314p1_Development/CMSSW_5_3_14_patch1/src/TopAnalysis/Configuration/analysis/diLeptonic/theoryHistograms/combined/VisGen"+variable+".root";
*   the theory curves should be stored in the ROOT file using the following naming convention
*       Madgrpah+Pythia:    Madgraph
*       Powheg+Pythia:      Powheg
*       Powheg+Herwig:      PowhegHerwig
*       MC@NLO+Herwig:      Mcatnlo
*   the ROOT file can be obtained by running
*       cd $CMSSW_BASE/src/TopAnalysis/Configuration/analysis/diLeptonic
*       root -l -b -q macros/obtainTheoryCurves.C++g
*/





#include "basicFunctions.h"
#include "TSystem.h"

void processQuantity(TString name, TString func, bool channel);
void addTAE(TGraphAsymmErrors* in, TGraphAsymmErrors* out);
TGraphAsymmErrors* seven(TString quantity, bool BCC);
TString DilepFileName(TString name);
int NbinsLL(TString name);
std::vector<double> BCCvalues(TString channel, TString name);
std::vector<double> DilepMCMadGraph(TString name);

void sevenToEightTeVComparison(){

    std::cout<<"\n\n\n\033[1;31m******************************************************************\n"
             <<"******************************************************************\033[1;m\n\n"
             <<"\033[1;34m Macro to obtain the dilepton vs l+jets differential cross section result comparison"
             <<"  if wanted 7TeV vs 8TeV results are also obtained\n"
             <<"      processQuantity(xSecVariables_[i], 'exp', false);\n\n"
             <<"  If wanted also the theory-to-Magraphd+Pythia ratio can be obtained\033[1;m\n\n\n" 
             <<" Exectue this macro from\n" 
             <<"     $CMSSW_BASE/src/TopAnalysis/Configuration/analysis/semiLeptonic/diffXSection\n"
             <<" via:\n"
             <<"     root -l -b -q sevenToEightTeVComparison.C++g\n\n"
             <<"\033[1;34m Output stored in\n"
             <<"      $CMSSW_BASE/src/TopAnalysis/Configuration/analysis/semiLeptonic/diffXSection/diffXSecFromSignal/plots/combined/2012/xSec/\n"
             <<"  which MUST exist a priori\033[1;m\n\n"
             <<"\033[1;34m Needed input:\033[1;m\n"
             <<"     Dilepton DiffXS 8TeV results in TGraphAssymErrors in a root file:\n"
             <<"         TString line= readLineFromFile(p+4,'/nfs/dust/cms/user/asincruz/CMSSW5314p1_Development/CMSSW_5_3_14_patch1/src/TopAnalysis/Configuration/analysis/diLeptonic/analysis_postCWR/Plots/FinalResults/combined/Hyp'+DilepFileName(name)+'LaTeX.txt)';\n"
             <<"     l+jets 8TeV DiffXS results stored on the NAF:\n"
             <<"         TString line= readLineFromFile(p+1,groupSpace+'CommonFiles/topPtInputForReweighting/diffXSecTopSemiLepParton'+name+'.txt');\n"
             <<"     The 7TeV l+jets results are hardcoded in the function 'seven'\n"
             <<"     The dilepton 7TeV results, theory and BCC are also hardcoded in the function 'dilepton7TeVCurves'\n"
             <<"     The dilepton 8TeV Madgraph+Pythia values are hardcoded in the function 'DilepMCMadGraph'\n"
             <<"     The BCCs for 8TeV for dilepton and l+jets are hardcoded in the function 'BCCvalues'\n"
             <<"\n If wanted ratio: theory/MadGraph+Pythia.\n"
             <<"     Uncomment the L.535.\n"
             <<" the function 'drawTheoryLines' will draw the theory/MG+PY ratio curves in all 8TeV ratio plots.\n"
             <<"   Needed input: the theory curves as TH1D in single ROOT (per variable)\n"
             <<"       TString filename = '/nfs/dust/cms/user/asincruz/CMSSW5314p1_Development/CMSSW_5_3_14_patch1/src/TopAnalysis/Configuration/analysis/diLeptonic/theoryHistograms/combined/VisGen'+variable+'.root';\n"
             <<"   the theory curves should be stored in the ROOT file using the following naming convention\n"
             <<"       Madgrpah+Pythia:    Madgraph\n"
             <<"       Powheg+Pythia:      Powheg\n"
             <<"       Powheg+Herwig:      PowhegHerwig\n"
             <<"       MC@NLO+Herwig:      Mcatnlo\n"
             <<"   the ROOT file can be obtained by running\n"
             <<"       cd $CMSSW_BASE/src/TopAnalysis/Configuration/analysis/diLeptonic\n"
             <<"       root -l -b -q macros/obtainTheoryCurves.C++g\n"
             <<"\033[1;31m******************************************************************\n"
             <<"******************************************************************\033[1;m\n\n\n\n\n"<<std::endl;




  // test quatity "" means none
  //TString TESTME="topPtTtbarSys";
  TString TESTME="";
  // collect full PS quantities
  std::vector<TString> xSecVariables_;
  xSecVariables_.insert(xSecVariables_.end(), xSecVariablesKinFit, xSecVariablesKinFit + sizeof(xSecVariablesKinFit)/sizeof(TString));
  // process all variations
  for(unsigned i=0; i<xSecVariables_.size(); ++i){
    if((TESTME==""||xSecVariables_[i]==TESTME)&&xSecVariables_[i]!="ttbarPhiStar"){
      if(xSecVariables_[i]!="ttbarPhiStar"){
	// decay channel
	processQuantity(xSecVariables_[i], "exp", true);   // dilepton vs ljets only at 8TeV
    processQuantity(xSecVariables_[i], "exp", false);  // dilepton vs ljets and 7TeV vs 8TeV 
	// sqrt(s)
	//if(xSecVariables_[i]=="topPt"||xSecVariables_[i]=="topY"||xSecVariables_[i]=="ttbarPt"||xSecVariables_[i]=="ttbarY"||xSecVariables_[i]=="ttbarMass") processQuantity(xSecVariables_[i], "exp", false);
      }
    }
  }
}


TGraphAsymmErrors* dilepton7TeVCurves(TString name)
{
    if (name == ""){
        std::cout <<" seventToEightTeVComparison::dilepton7TeVCurves: you didn't provide a variable name.\nCannot continue" << std::endl;
        exit(11);
    }
    std::cout << __LINE__<< " Variable name = " << name << std::endl;
        
    int nbins = 5;
    if(      name=="ttbarMass") nbins = 7;
    else if (name=="ttbarPt"  ) nbins = 4;
    else if (name=="ttbarY"   ) nbins = 6;
    else if (name=="topPt"    ) nbins = 5;
    else if (name=="topY"     ) nbins = 8;

    // create asymmerrors
    TGraphAsymmErrors* out = new TGraphAsymmErrors(nbins);
    
    std::vector<double> v_data, v_binCenter, v_theory, v_error;
    
    if (name == "topPt")
    {
        double data[] = {0.00509572, 0.00626002, 0.00296467, .000701592, 0.00012036};
        double error[] = {0.0601381, 0.0469906, 0.0555114, 0.071274, 0.0924826};
        double theory[] = {0.00453114, 0.00600115, 0.00321705, 0.000931674, 0.000191065};
        double binCenter[] = {34.7, 107., 162., 242., 343.};
        v_data.assign(data, data+nbins);
        v_error.assign(error, error+nbins);
        v_theory.assign(theory, theory+nbins);
        v_binCenter.assign(binCenter, binCenter+nbins);
    }
    if (name == "topY")
    {
        double data[] = {0.091, 0.255, 0.302, 0.351, 0.371, 0.306, 0.241, 0.090};
        double error[] = {0.080, 0.058, 0.052, 0.050, 0.049, 0.053, 0.059, 0.079};
        double theory[] = {0.08688615, 0.2418012, 0.3199502, 0.3587744, 0.3588348, 0.3205259, 0.2424955, 0.08722365};
        double binCenter[] = {-1.85, -1.05, -0.61, -0.23, 0.23, 0.61, 1.05, 1.85};
        v_data.assign(data, data+nbins);
        v_error.assign(error, error+nbins);
        v_theory.assign(theory, theory+nbins);
        v_binCenter.assign(binCenter, binCenter+nbins);
    }
    if (name == "ttbarPt")
    {
        double data[] = {0.0160, 0.0097, 0.0032, 0.0005};
        double error[] = {0.250, 0.109, 0.135, 0.079};
        double theory[] = {0.01516833, 0.009621917, 0.003224087, 0.0005992394};
        double binCenter[] = {5., 39., 86., 231.};
        v_data.assign(data, data+nbins);
        v_error.assign(error, error+nbins);
        v_theory.assign(theory, theory+nbins);
        v_binCenter.assign(binCenter, binCenter+nbins);
    }
    if (name == "ttbarY")
    {
        double data[] = {0.030, 0.219, 0.418, 0.393, 0.218, 0.040};
        double error[] = {0.201, 0.054, 0.043, 0.044, 0.055, 0.180};
        double theory[] = {0.04102893, 0.2246687, 0.3976315, 0.3982489, 0.2258003, 0.04113086};
        double binCenter[] = {-1.95, -1.1, -0.4, 0.4, 1.1, 1.95};
        v_data.assign(data, data+nbins);
        v_error.assign(error, error+nbins);
        v_theory.assign(theory, theory+nbins);
        v_binCenter.assign(binCenter, binCenter+nbins);
    }
    if (name == "ttbarMass")
    {
        double data[] = {0.00526, 0.00458, 0.00246, 0.00107, 0.00039, 0.00008, 0.00001};
        double error[] = {0.117, 0.056, 0.090, 0.072, 0.134, 0.301, 0.490};
        double theory[] = {0.004905427, 0.004427047, 0.002457022, 0.0012, 0.00042, 0.000091, 0.0000088};
        double binCenter[] = {364., 439., 515., 616., 660., 857., 1275.0};
        v_data.assign(data, data+nbins);
        v_error.assign(error, error+nbins);
        v_theory.assign(theory, theory+nbins);
        v_binCenter.assign(binCenter, binCenter+nbins);
    }
    
    for (int iter = 0; iter < nbins; iter++)
    {
        double ratio = v_data.at(iter) / v_theory.at(iter);
        double rel_error  = 1. * v_error.at(iter) * ratio;
        out->SetPoint( iter, v_binCenter.at(iter),  ratio );
        out->SetPointError( iter, 0., 0., rel_error, rel_error);
    }
    return out;

}

/// Check existence file. Return 1 is exists, 0 if not. 
bool checkFileExistence(TString file)
{
    std::ifstream f(file);
    if (f.good())
    {
        f.close();
        return true;
    } else {
        std::cout<<"File '"<<file<<"'\nDoes not exist.\nEXIT!"<<std::endl;
        f.close();
        return false;
    }
}



void setStyle(TH1* histo, TString theoryName, TLegend *leg = 0)
{
    if(!histo) return;
    histo->SetLineWidth(3);
    histo->Scale(1./histo->Integral("width"));
    if(theoryName == "madgraph"){
        histo->SetLineColor(kRed+1);
        histo->SetLineStyle(1);
//         if(leg) leg->AddEntry(histo, "MadGraph+Pythia",  "l");
    }
    if(theoryName == "powheg"){
        histo->SetLineColor(kGreen+1);
        histo->SetLineStyle(7);
        if(leg) leg->AddEntry(histo, "Powheg+Pythia",  "l");
    }
    if(theoryName == "powhegherwig"){
        histo->SetLineColor(kGreen+3);
        histo->SetLineStyle(9);
        if(leg) leg->AddEntry(histo, "Powheg+Herwig",  "l");
    }
    if(theoryName == "mcatnlo"){
        histo->SetLineColor(kBlue);
        histo->SetLineStyle(5);
        if(leg) leg->AddEntry(histo, "MC@NLO+Herwig",  "l");
    }
}


void drawTheoryLines(TString name, float xmin=-1000, float xmax = -1000, int nrebin = 1)
{
    TString variable = "";

    if (name.Contains("topPtSubLead")) variable = TString("VisGenToppTNLead");
    if (name.Contains("topPtLead"))    variable = TString("VisGenToppTLead");
    if (name.Contains("topPtTtbarSys"))variable = TString("VisGenToppTTTRestFrame");
    if (name.Contains("topPt"))        variable = TString("VisGenToppT");
    if (name.Contains("topY"))         variable = TString("VisGenTopRapidity");
    if (name.Contains("ttbarDelPhi"))  variable = TString("VisGenTTBarDeltaPhi");
    if (name.Contains("ttbarMass"))    variable = TString("VisGenTTBarMass");
    if (name.Contains("ttbarPt"))      variable = TString("VisGenTTBarpT");
    if (name.Contains("ttbarY"))       variable = TString("VisGenTTBarRapidity");


    TString filename =  "/nfs/dust/cms/user/asincruz/CMSSW5314p1_Development/CMSSW_5_3_14_patch1/src/TopAnalysis/Configuration/analysis/diLeptonic/theoryHistograms/combined/"+variable+".root";
    if(!checkFileExistence(filename))
    {
        std::cout<<"*-*-*-*-*-*-File '"<<filename<<"'\nDoes not exist.\nEXIT"<<std::endl;
        exit(23);
    }
    TFile *curves_file = new TFile(filename);
    
    TH1D *madgraph     = (TH1D *)curves_file->Get("Madgraph")->Clone("madgraph");
    TH1D *powheg       = (TH1D *)curves_file->Get("Powheg")->Clone("powheg");
    TH1D *powhegherwig = (TH1D *)curves_file->Get("PowhegHerwig")->Clone("powhegherwig");
    TH1D *mcatnlo      = (TH1D *)curves_file->Get("Mcatnlo")->Clone("mcatnlo");

    if(!madgraph)
    {
        std::cout<<"Madrgraph theory curve not available\nEXIT"<<std::endl;
        exit(22); 
    }

    madgraph->Rebin(nrebin);
    powheg->Rebin(nrebin);
    powhegherwig->Rebin(nrebin);
    mcatnlo->Rebin(nrebin);

    if(xmin != -1000 && xmax != -1000)
    {
        madgraph->GetXaxis()->SetRangeUser(xmin, xmax);
        powheg->GetXaxis()->SetRangeUser(xmin, xmax);
        powhegherwig->GetXaxis()->SetRangeUser(xmin, xmax);
        mcatnlo->GetXaxis()->SetRangeUser(xmin, xmax);
    }
    
    TLegend *leg = new TLegend();
    setStyle(madgraph, "madgraph", leg);
    setStyle(powheg, "powheg", leg);
    setStyle(powhegherwig, "powhegherwig", leg);
    setStyle(mcatnlo, "mcatnlo", leg);

    leg->SetX1NDC(0.225);    leg->SetX2NDC(0.45);
    leg->SetY1NDC(0.17);    leg->SetY2NDC(0.30);
    leg->SetTextFont(42);   leg->SetTextAlign(12);
    leg->SetTextSize(0.04);
    leg->SetFillStyle(0);   leg->SetBorderSize(0);


    powheg->Divide(madgraph);
    powhegherwig->Divide(madgraph);
    mcatnlo->Divide(madgraph);
    
    powheg->Draw("same");
    powhegherwig->Draw("same");
    mcatnlo->Draw("same");
    leg->Draw("same");
}

void processQuantity(TString name, TString func, bool channel){
  // channel == true  -> 8TeV ljets vs dilepton plots
  // channel == false -> ljets 7TeV vs 8TeV
  // some parameters
  bool BCC=true;//name=="topPt"||!channel ? true : false;
  bool grey=false;
  bool plotsepfit   =true;
  bool plotfiterrors=false;
  TString optD= plotsepfit ? "" : "0";
  double ratmax=channel ? 1.8 : (name=="ttbarMass"  ? 2.75 : 1.5);
  double ratmin=channel ? 0.5 : (name.Contains("Y") ? 0.7 : (name == "ttbarMass" ? 0.375 : 0.6));

  bool isAvailableAt7TeV = 0;
  if(!channel && (name == "topPt" || name == "topY" || name == "ttbarPt" || name == "ttbarMass" || name == "ttbarY")) isAvailableAt7TeV = 1;


  // colors
  int color7     = kRed+1; //kRed-4
  int color8     = kGreen; //kBlue+2 
  int ljets7color= color7;
  int dilep7color= color7;
  int ljets8color= channel ? kBlue+2 : color8;
  int dilep8color= channel ? kRed-4 : color8;
  int fit7color  = color7;
  int fit8color  = color8;
  int colorband  = kCyan-7;
  int nrebinTheoryCurves = 1;
  // ---
  //    canvas style 
  // ---
  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  myStyle.SetErrorX(0.5);
  myStyle.cd();
  gROOT->SetStyle("HHStyle");
  gStyle->SetEndErrorSize(10);
  gStyle->SetOptFit(0);
  
  // ---
  //    collect all curves
  // ---
  // 7 TeV
  //int NbinsLjets7=7;
  //TGraphAsymmErrors* SFljets = new TGraphAsymmErrors(NbinsLjets7);
  TGraphAsymmErrors* SFljets =seven(name,BCC);
  int NbinsDilep7=5;
  TGraphAsymmErrors* SFdilep = isAvailableAt7TeV ? dilepton7TeVCurves(name) : new TGraphAsymmErrors(NbinsDilep7);
  //int Nbins7=NbinsLjets7+NbinsDilep7;
  TGraphAsymmErrors* SF7 = new TGraphAsymmErrors(0);
  // 8 TeV
  std::vector<double> ljetsBinning_=makeVariableBinning(false)[name];
  int NbinsLjets8=ljetsBinning_.size()-1;
  double xmin=ljetsBinning_.at(0);
  double xmax=ljetsBinning_.at(ljetsBinning_.size()-1);
  TString Txmax=getTStringFromDouble(xmax);
  TGraphAsymmErrors* SFljets8 = new TGraphAsymmErrors(NbinsLjets8);
  int NbinsDilep8=NbinsLL(name);
  TGraphAsymmErrors* SFdilep8 = new TGraphAsymmErrors(NbinsDilep8);
  //int Nbins8=NbinsLjets8+NbinsDilep8;
  TGraphAsymmErrors* SF8 = new TGraphAsymmErrors(0);
  // combined
  //int Nbins=NbinsLjets7+NbinsDilep7+NbinsLjets8+NbinsDilep8;
  TGraphAsymmErrors* SF = new TGraphAsymmErrors(0);

  // extension for different modes
  TString nameext= channel ? "DecayChannel" : "Energy";
  std::cout << name << "(" << nameext << ")" << std::endl;
  // style of 7TeV ratio
  SFljets->SetLineWidth(3.);
  SFljets->SetMarkerSize(1.5);
  SFljets->SetMarkerStyle(22);
  SFljets->SetMarkerColor(ljets7color);
  SFljets->SetLineColor(ljets7color);

/*  // b) dilepton 7TeV data points
  //           bin x(BCC)    data  / Madgraph               // BCCNNLO // BCC MG
  SFdilep->SetPoint( 0, 33.7,  (0.00509572 / 0.00453114 )  );// 33.7    // 34 
  SFdilep->SetPoint( 1, 107 ,  (0.00626002 / 0.00600115 )  );// 106     // 107
  SFdilep->SetPoint( 2, 162 ,  (0.00296467 / 0.00321705 )  );// 162     // 163
  SFdilep->SetPoint( 3, 242 ,  (0.000701592/ 0.000931674)  );// 242     // 247
  SFdilep->SetPoint( 4, 343 ,  (0.00012036 / 0.000191065)  );// 343     // 350
  //                   x errors   rel.err(data) *( data  / Madgraph)
  SFdilep->SetPointError( 0, 0., 0., 0.0601381*(0.00509572 / 0.00453114 ), 0.0601381*(0.00509572 / 0.00453114 ) );
  SFdilep->SetPointError( 1, 0., 0., 0.0469906*(0.00626002 / 0.00600115 ), 0.0469906*(0.00626002 / 0.00600115 ) );
  SFdilep->SetPointError( 2, 0., 0., 0.0555114*(0.00296467 / 0.00321705 ), 0.0555114*(0.00296467 / 0.00321705 ) );
  SFdilep->SetPointError( 3, 0., 0., 0.071274* (0.000701592/ 0.000931674), 0.071274* (0.000701592/ 0.000931674) );
  SFdilep->SetPointError( 4, 0., 0., 0.0924826*(0.00012036 / 0.000191065), 0.0924826*(0.00012036 / 0.000191065) );*/
  //style of ratio
  SFdilep->SetLineWidth(3.);
  SFdilep->SetMarkerSize(1.5);
//   SFdilep->SetMarkerStyle(22);
  SFdilep->SetMarkerStyle(26);
  SFdilep->SetMarkerColor(dilep7color);
  SFdilep->SetLineColor(dilep7color);

  // collect 8 TeV BCC x values for analysis binning 
  std::vector<double> xBCCljets_;
  std::vector<double> xBCCdilep_;
  if(BCC){
    xBCCljets_=BCCvalues("ljets"   , name);
    xBCCdilep_=BCCvalues("dilepton", name);
  }
  if(!BCC||xBCCljets_.size()==0){
    for(unsigned int i=0; i<ljetsBinning_.size()-1; ++i){
      double value=ljetsBinning_[i]+0.5*(ljetsBinning_[i+1]-ljetsBinning_[i]);
      //std::cout << value << std::endl;
      xBCCljets_.push_back(value);
    }
  }
  // c) l+jets 8TeV data points
  for(int p=0; p<NbinsLjets8; ++p){
    // get line with all informations
    TString line= readLineFromFile(p+1, groupSpace+"CommonFiles/topPtInputForReweighting/diffXSecTopSemiLepParton"+name+".txt");
    //std::cout << line << std::endl;
    // data value
    TString temp = getStringEntry(line, 3 , "&");
    temp.ReplaceAll(" ","");
    double data=atof(temp.Data());
    temp = getStringEntry(line, 2 , "&");
    temp.ReplaceAll(" ","");
    double MC  =atof(temp.Data());
    SFljets8->SetPoint( p,  xBCCljets_.at(p) , data/MC ); 
    temp = getStringEntry(line, 6 , "&");
    temp.ReplaceAll(" ","");
    temp.ReplaceAll("\\","");
    double unc=atof(temp.Data());
    SFljets8->SetPointError( p, 0., 0., (unc/100.)*(data/MC), (unc /100.)*(data/MC) );
  }
  whipEmptyBinsAway(SFljets8, 0);
  //style of ratio
  SFljets8->SetLineWidth(3.);
  SFljets8->SetMarkerSize(1.5);
//   SFljets8->SetMarkerStyle(24);
  SFljets8->SetMarkerStyle(channel ? 24 : 20);
  //SFljets8->SetLineStyle(2);
  SFljets8->SetMarkerColor(ljets8color);
  SFljets8->SetLineColor(ljets8color);

  // d) dilepton 8TeV data points
  // MC prediction point (as not in provided table)
  std::vector<double> MCdilep_=DilepMCMadGraph(name);
  if(channel || 1){
    for(int p=0; p<NbinsDilep8; ++p){
      // get line with all informations
//       if(p==0) std::cout << groupSpace+"CommonFiles/topPtInputForReweighting/Hyp"+DilepFileName(name)+"LaTeX.txt" << std::endl;
//       TString line= readLineFromFile(p+4, groupSpace+"CommonFiles/topPtInputForReweighting/Hyp"+DilepFileName(name)+"LaTeX.txt");
      if(p==0) std::cout << "/nfs/dust/cms/user/asincruz/CMSSW5314p1_Development/CMSSW_5_3_14_patch1/src/TopAnalysis/Configuration/analysis/diLeptonic/analysis_postCWR/Plots/FinalResults/combined/Hyp"+DilepFileName(name)+"LaTeX.txt" << std::endl;
      TString line= readLineFromFile(p+4, "/nfs/dust/cms/user/asincruz/CMSSW5314p1_Development/CMSSW_5_3_14_patch1/src/TopAnalysis/Configuration/analysis/diLeptonic/analysis_postCWR/Plots/FinalResults/combined/Hyp"+DilepFileName(name)+"LaTeX.txt");
      //std::cout << line << std::endl;
      // data value
      TString datatemp= getStringEntry(line, 3, "&");
      TString xctemp  = getStringEntry(line, 1, "$");
      TString x1temp  = getStringEntry(getStringEntry(line, 2, "&"), 2, "$");
      TString x2temp  = getStringEntry(getStringEntry(line, 2, "&"), 4, "$");
      double data   =atof(datatemp.Data());
      double xcenter=atof(xctemp.Data());
      double x1D    =atof(x1temp .Data());
      double x2D    =atof(x2temp .Data());
      double xBCC= ((int)xBCCdilep_.size()==NbinsDilep8 ? xBCCdilep_.at(p) : x1D+0.5*(x2D-x1D));
      //std::cout << "->x1D :" << x1D << std::endl;
      //std::cout << "->x2D :" << x2D << std::endl;
      //std::cout << "->xBCC :" << BCCtemp << "->" << xBCC << std::endl;
      //datatemp.ReplaceAll(" ","");
      //temp = getStringEntry(line, 2 , "&");
      //temp.ReplaceAll(" ","");
      double MC  = MCdilep_[p];
      SFdilep8->SetPoint( p,  BCC ? xBCC : xcenter, data/MC ); 
      TString unctemp = getStringEntry(line, 6 , "&");
      double unc=atof(unctemp.Data());
      SFdilep8->SetPointError( p, BCC ? 0. : 0.5*(x2D-x1D), BCC ? 0. : 0.5*(x2D-x1D), (unc/100.)*(data/MC), (unc /100.)*(data/MC) );
    }
  }
  //style of ratio
  SFdilep8->SetLineWidth(3.);
  SFdilep8->SetMarkerSize(1.5);
//   SFdilep8->SetMarkerStyle(22);
  SFdilep8->SetMarkerStyle(channel ? 22 : 24);
  SFdilep8->SetMarkerColor(dilep8color);
  SFdilep8->SetLineColor(dilep8color);

  // e) combined 7 TeV data points
  addTAE(SFdilep , SF7);
  addTAE(SFljets , SF7);
  //style of ratio
  SF7->SetLineWidth(3.);
  SF7->SetMarkerSize(0.1);
  SF7->SetMarkerStyle(20);
  SF7->SetMarkerColor(kWhite);
  SF7->SetLineColor(kWhite);

  // f) combined 8 TeV data points
  if(channel){
    addTAE(SFdilep8, SF8);
    addTAE(SFljets8, SF8);
  }
  //style of ratio
  SF8->SetLineWidth(3.);
  SF8->SetMarkerSize(0.1);
  SF8->SetMarkerStyle(20);
  SF8->SetMarkerColor(kWhite);
  SF8->SetLineColor(kWhite);

  // g) combined 7+8TeV data points
  //addTAE(SF7, SF);
  //addTAE(SF8, SF);
  //style of ratio
  SF->SetLineWidth(3.);
  SF->SetMarkerSize(0.1);
  SF->SetMarkerStyle(20);
  SF->SetMarkerColor(kWhite);
  SF->SetLineColor(kWhite);

  // ---
  //    dummy plots for axis
  // ---
  TH1F* dummy= new TH1F("","",1,xmin,xmax);
  histogramStyle(*dummy, kSig);
  TString xlabel="";
  for(unsigned int i=0; i<sizeof(xSecVariablesKinFit)/sizeof(TString);++i){
    if(xSecVariablesKinFit[i]==name) xlabel=xSecLabelKinFit[i];    
  }
  //std::cout << "TEST: " << xlabel << std::endl;
  dummy->GetXaxis()->SetTitle(getStringEntry(xlabel, 1, "/")+" "+getStringEntry(xlabel, 2, "/"));
  dummy->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+getStringEntry(xlabel, 1, "/")+"} Ratio: (Data / Simulation)");
  dummy->GetYaxis()->SetTitleOffset(0.9*dummy->GetYaxis()->GetTitleOffset());
  dummy->SetMaximum(ratmax);
  dummy->SetMinimum(ratmin);

  // ---
  //    legends
  // ---
  double x1=0.25;
  double x2=0.88;
    
  TLegend *leg0 = new TLegend(x1, 0.67, x2, 0.87);
  leg0->SetFillStyle(0);
  leg0->SetTextSize(0.035);
  leg0->SetBorderSize(0);
//   leg0->SetHeader("#font[22]{Data / MadGraph+Pythia(CTEQ6L1)}");
  leg0->SetHeader("#font[22]{Data / MadGraph+Pythia}");

  TLegend *leg1 = new TLegend(x1, 0.57, x2, 0.69);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.035);
  leg1->SetBorderSize(0);
  leg1->SetHeader("#font[22]{Fit: exp(a+b#upointx)}");

  if(plotsepfit) leg1->SetY1(leg1->GetY1()-0.2);

  // canvas
  std::vector<TCanvas*> plotCanvas_;
  addCanvas(plotCanvas_);
  plotCanvas_[plotCanvas_.size()-1]->cd(0);
  plotCanvas_[plotCanvas_.size()-1]->SetTitle("data/MC "+name+"("+nameext+") ratio");
  // drawing
  dummy->Draw("axis");
  drawLine(xmin, 1.0, xmax, 1.0, kBlack, 4, 2); // line at one
  //SF->Draw("p e1 same");
  //SF7->Draw("p e1 same");
  //SF8->Draw("p e1 same");
//   if(channel) drawTheoryLines(name, xmin, xmax, nrebinTheoryCurves);
  if(isAvailableAt7TeV) SFljets->Draw("p e1 same");
  if(isAvailableAt7TeV) SFdilep->Draw("p e1 same");
//   if( channel) SFdilep8->Draw("p e1 same");
  SFdilep8->Draw("p e1 same");
  SFljets8->Draw("p e1 same");
  // fit polynomial or exponential function
  TString def = "";
  if(func=="pol2")def="[0]*x*x+[1]*x+[2]";
  if(func=="exp" )def="exp([0]+[1]*x)";
  double fitLowEdge=0.;
  double fitHighEdge=xmax;
  // a) to all 8 and 7 TeV points
//   TF1* function=new TF1("function",def,fitLowEdge, fitHighEdge);
//   function->SetLineColor(kMagenta+2);
//   SF->Fit(function,"R","same",fitLowEdge, fitHighEdge);
//   for(int i=0; i<function->GetNumberFreeParameters(); i++){
//     function->SetParameter(i,round(function->GetParameter(i),3));
//   }
//   TString fitEntry="#splitline{}{#splitline{}{#splitline{}{#splitline{combined fit: ";
//   fitEntry+=function->GetExpFormula("p")+",}{                          #chi^{2}/ndof=";
//   fitEntry+=getTStringFromDouble(function->GetChisquare())+"/"+getTStringFromInt(function->GetNDF())+"}}}}";
//   fitEntry.ReplaceAll("+(","");
//   fitEntry.ReplaceAll("))",")");
//   leg0->AddEntry( function, fitEntry, "L");
  // b) to all 7 TeV points
  TF1* function7=new TF1("function7",def,fitLowEdge, fitHighEdge);
  function7->SetLineColor(fit7color);
  function7->SetLineWidth(6);
  function7->SetLineStyle(2);
  //SF7->Fit(function7,"R","same",fitLowEdge, fitHighEdge);
  //for(int i=0; i<function7->GetNumberFreeParameters(); i++){
  //  function7->SetParameter(i,round(function7->GetParameter(i),3));
  //}
  //TString fitEntry7="fit 7 TeV: ";
  //fitEntry7+=function7->GetExpFormula("p");
  //fitEntry7+=",  #chi^{2}/ndof=";
  //fitEntry7+=getTStringFromDouble(function7->GetChisquare())+"/"+getTStringFromInt(function7->GetNDF());
  //fitEntry7.ReplaceAll("+(","");
  //fitEntry7.ReplaceAll("))",")");
  TString fitEntry7="7 TeV: ";
  if(plotfiterrors) fitEntry7+="              ";
  fitEntry7+="a=";
  fitEntry7+=getTStringFromDouble(function7->GetParameter(0), 3);
  if(plotfiterrors){
    fitEntry7+="#pm";
    fitEntry7+=getTStringFromDouble(function7->GetParError(0) , 3);
  }
  fitEntry7+=", b=";
  fitEntry7+=getTStringFromDouble(function7->GetParameter(1), 5);
  if(plotfiterrors){
    fitEntry7+="#pm";
    fitEntry7+=getTStringFromDouble(function7->GetParError(1) , 5);
  }

  // b1) to l+jets 7 TeV points
  TF1* functionljets7=new TF1("functionljets7",def,fitLowEdge, fitHighEdge);
  functionljets7->SetLineColor(kRed+1);
  functionljets7->SetLineWidth(2);
  //SFljets->Fit(functionljets7,"R"+optD,"same",fitLowEdge, fitHighEdge);
  //for(int i=0; i<functionljets7->GetNumberFreeParameters(); i++){
  //  functionljets7->SetParameter(i,round(functionljets7->GetParameter(i),3));
  //}
  //TString fitEntryljets7="fit 7 TeV l+jets: ";
  //fitEntryljets7+=functionljets7->GetExpFormula("p");
  //fitEntryljets7+=",  #chi^{2}/ndof=";
  //fitEntryljets7+=getTStringFromDouble(functionljets7->GetChisquare())+"/"+getTStringFromInt(functionljets7->GetNDF());
  //fitEntryljets7.ReplaceAll("+(","");
  //fitEntryljets7.ReplaceAll("))",")");
  TString fitEntryljets7="7 TeV l+jets:     ";
  fitEntryljets7+="a=";
  fitEntryljets7+=getTStringFromDouble(functionljets7->GetParameter(0), 3);
  if(plotfiterrors){
    fitEntryljets7+="#pm";
    fitEntryljets7+=getTStringFromDouble(functionljets7->GetParError(0) , 3);
  }
  fitEntryljets7+=", b=";
  fitEntryljets7+=getTStringFromDouble(functionljets7->GetParameter(1), 5);
  if(plotfiterrors){
    fitEntryljets7+="#pm";
    fitEntryljets7+=getTStringFromDouble(functionljets7->GetParError(1) , 5);
  }

  // b2) to dilepton 7 TeV points
  TF1* functiondilep7=new TF1("functiondilep7",def,fitLowEdge, fitHighEdge);
  functiondilep7->SetLineColor(kOrange+7);
  functiondilep7->SetLineWidth(2);
  //SFdilep->Fit(functiondilep7,"R"+optD,"same",fitLowEdge, fitHighEdge);
  //for(int i=0; i<functiondilep7->GetNumberFreeParameters(); i++){
  //  functiondilep7->SetParameter(i,round(functiondilep7->GetParameter(i),3));
  //}
  //TString fitEntrydilep7="fit 7 TeV dilepton: ";
  //fitEntrydilep7+=functiondilep7->GetExpFormula("p");
  //fitEntrydilep7+=",  #chi^{2}/ndof=";
  //fitEntrydilep7+=getTStringFromDouble(functiondilep7->GetChisquare())+"/"+getTStringFromInt(functiondilep7->GetNDF());
  //fitEntrydilep7.ReplaceAll("+(","");
  //fitEntrydilep7.ReplaceAll("))",")");
  TString fitEntrydilep7="7 TeV dilepton: ";
  fitEntrydilep7+="a=";
  fitEntrydilep7+=getTStringFromDouble(functiondilep7->GetParameter(0), 3);
  if(plotfiterrors){
    fitEntrydilep7+="#pm";
    fitEntrydilep7+=getTStringFromDouble(functiondilep7->GetParError(0) , 3);
  }
  fitEntrydilep7+=", b=";
  if(plotfiterrors){
    fitEntrydilep7+=getTStringFromDouble(functiondilep7->GetParameter(1), 5);
    fitEntrydilep7+="#pm";
  }
  fitEntrydilep7+=getTStringFromDouble(functiondilep7->GetParError(1) , 5);

  // c1) to l+jets 8 TeV points
  TF1* functionljets8=new TF1("functionljets8",def,fitLowEdge, fitHighEdge);
  functionljets8->SetLineColor(kBlue);
  functionljets8->SetLineWidth(2);
  //SFljets8->Fit(functionljets8,"R"+optD,"same",fitLowEdge, fitHighEdge);
  //for(int i=0; i<functionljets8->GetNumberFreeParameters(); i++){
  //  functionljets8->SetParameter(i,round(functionljets8->GetParameter(i),3));
  //}
  //TString fitEntryljets8="fit 8 TeV l+jets: ";
  //fitEntryljets8+=functionljets8->GetExpFormula("p");
  //fitEntryljets8+=",  #chi^{2}/ndof=";
  //fitEntryljets8+=getTStringFromDouble(functionljets8->GetChisquare())+"/"+getTStringFromInt(functionljets8->GetNDF());
  //fitEntryljets8.ReplaceAll("+(","");
  //fitEntryljets8.ReplaceAll("))",")");
  TString fitEntryljets8="8 TeV l+jets:     ";
  fitEntryljets8+="a=";
  fitEntryljets8+=getTStringFromDouble(functionljets8->GetParameter(0), 3);
  if(plotfiterrors){
    fitEntryljets8+="#pm";
    fitEntryljets8+=getTStringFromDouble(functionljets8->GetParError(0) , 3);
  }
  fitEntryljets8+=", b=";
  fitEntryljets8+=getTStringFromDouble(functionljets8->GetParameter(1), 5);
  if(plotfiterrors){
    fitEntryljets8+="#pm";
    fitEntryljets8+=getTStringFromDouble(functionljets8->GetParError(1) , 5);
  }

  // c2) to dilepton 8 TeV points
  TF1* functiondilep8=new TF1("functiondilep8",def,fitLowEdge, fitHighEdge);
  functiondilep8->SetLineColor(kAzure+6);
  functiondilep8->SetLineWidth(2);
  
  //SFdilep8->Fit(functiondilep8,"R"+optD,"same",fitLowEdge, fitHighEdge);
  //for(int i=0; i<functiondilep8->GetNumberFreeParameters(); i++){
  //  functiondilep8->SetParameter(i,round(functiondilep8->GetParameter(i),3));
  //}
  //TString fitEntrydilep8="fit 8 TeV dilepton: ";
  //fitEntrydilep8+=functiondilep8->GetExpFormula("p");
  //fitEntrydilep8+=",  #chi^{2}/ndof=";
  //fitEntrydilep8+=getTStringFromDouble(functiondilep8->GetChisquare())+"/"+getTStringFromInt(functiondilep8->GetNDF());
  //fitEntrydilep8.ReplaceAll("+(","");
  //fitEntrydilep8.ReplaceAll("))",")");
  TString fitEntrydilep8="8 TeV dilepton: ";
  fitEntrydilep8+="a=";
  fitEntrydilep8+=getTStringFromDouble(functiondilep8->GetParameter(0), 3);
  if(plotfiterrors){
    fitEntrydilep8+="#pm";
    fitEntrydilep8+=getTStringFromDouble(functiondilep8->GetParError(0) , 3);
  }
  fitEntrydilep8+=", b=";
  fitEntrydilep8+=getTStringFromDouble(functiondilep8->GetParameter(1), 5);
  if(plotfiterrors){
    fitEntrydilep8+="#pm";
    fitEntrydilep8+=getTStringFromDouble(functiondilep8->GetParError(1) , 5);
  }

  // Draw legend
  //leg0->AddEntry(SFljets8,"8 TeV, e/#mu+jets "+TString(PHD ? "(this thesis)" : "(CMS-PAPER-TOP-12-028)"), "P");
  leg0->AddEntry(SFljets8,"e/#mu + Jets Combined"+TString(PHD ? "(this thesis)" : "") + TString(!channel ? " (8 TeV)" : ""), "P");
//   if(!channel) leg0->AddEntry(SFljets, "#splitline{7 TeV, e/#mu+jets}{(Eur. Phys. J. C73 (2013) 2339)}" , "P");
  if(!channel && isAvailableAt7TeV) leg0->AddEntry(SFljets, "e/#mu + Jets Combined (7 TeV)" , "P");
//   //leg0->AddEntry(SFdilep, "7 TeV ee/e#mu/#mu#mu (TOP-11-013)", "P");
  leg0->AddEntry(SFdilep8,"Dilepton Combined" + TString(!channel ? " (8 TeV)" : ""), "P");
  if(!channel && isAvailableAt7TeV) leg0->AddEntry(SFdilep, "Dilepton Combined (7 TeV)", "P");
  //if(channel) leg0->AddEntry(SFdilep8,"#splitline{8 TeV, ee/e#mu/#mu#mu}{(CMS-PAPER-TOP-12-028)}", "P");
    leg0->Draw("same");
  //leg1->AddEntry( function7, fitEntry7, "L");
  //if(plotsepfit) leg1->AddEntry( functiondilep7, fitEntrydilep7, "L");
  //if(plotsepfit) leg1->AddEntry( functionljets7, fitEntryljets7, "L");
  //leg1->AddEntry( function8, fitEntry8, "L");
  //if(plotsepfit) leg1->AddEntry( functiondilep8, fitEntrydilep8, "L");
  //if(plotsepfit) leg1->AddEntry( functionljets8, fitEntryljets8, "L");
  //leg1->Draw("same");
  // Draw cms label
  TPaveText *label = new TPaveText();
  label -> SetX1NDC(gStyle->GetPadLeftMargin());
  label -> SetY1NDC(1.0-gStyle->GetPadTopMargin());
  label -> SetX2NDC(1.0-gStyle->GetPadRightMargin());
  label -> SetY2NDC(1.0);
  label -> SetTextFont(42);
  TString CMSlab="";
  //if(!PHD) CMSlab+="CMS Preliminary, ";  
  if(!PHD) CMSlab+="CMS, ";  
  if(channel || !isAvailableAt7TeV) CMSlab+="19.7 fb^{-1} at #sqrt{s} = 8 TeV";
  else if (isAvailableAt7TeV) CMSlab+="5.0/19.7 fb^{-1} at #sqrt{s} = 7/8 TeV";
  label -> AddText(CMSlab);
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  label->SetTextSize(0.04);
  label->SetTextAlign(32);
  label-> Draw("same");
  // BCC label
  double positionX=xmax+0.045*(xmax-xmin)*(gStyle->GetCanvasDefW()/600.);
  double positionY=ratmin;
  TLatex *bcclabel = new TLatex(positionX,positionY, " (horizontal BCC wrt. MadGraph+Pythia)");
  bcclabel->SetTextAlign(11);
  bcclabel->SetTextAngle(90);
  bcclabel->SetTextSize(0.035);
  bcclabel->Draw("same");
  if(grey) plotCanvas_[0]->SetGrayscale();
  //saving
  plotCanvas_[0]->Print("diffXSecFromSignal/plots/combined/2012/xSec/dataVsMadgraph"+name+nameext+"7and8TeV.eps");
  plotCanvas_[0]->Print("diffXSecFromSignal/plots/combined/2012/xSec/dataVsMadgraph"+name+nameext+"7and8TeV.png");
  // ---
  // ERROR band plot
  // ---
  if(channel&&name=="topPtTtbarSys"){
    addCanvas(plotCanvas_);
    plotCanvas_[plotCanvas_.size()-1]->cd(0);
    plotCanvas_[plotCanvas_.size()-1]->SetTitle("data/MC top Pt errorband");
    // drawing
    // - axis
    dummy->GetYaxis()->SetTitle("(Data / Simulation) SF (#sqrt{s}=8 TeV)");
    dummy->GetYaxis()->SetRangeUser(0.35, 1.65);
    dummy->SetFillColor(10);
    dummy->SetLineColor(10);
    dummy->Draw("axis");
    // - 8 TeV data points
    SF8->Draw("p same");
    // perform fit to all 8 TeV points
    TF1* function8=new TF1("function8",def,fitLowEdge, fitHighEdge);
    function8->SetLineWidth(6);
    function8->SetLineColor(fit8color);
    function8->SetLineStyle(2);
    SF8->Fit(function8,"R","same",fitLowEdge, fitHighEdge);
    for(int i=0; i<function8->GetNumberFreeParameters(); i++){
      function8->SetParameter(i,round(function8->GetParameter(i),3));
    }
    TString fitEntry8="8 TeV: ";
    if(plotfiterrors) fitEntry8+="              ";
    fitEntry8+="a=";
    fitEntry8+=getTStringFromDouble(function8->GetParameter(0), 3);
    if(plotfiterrors){
      fitEntry8+="#pm";
      fitEntry8+=getTStringFromDouble(function8->GetParError(0) , 3);
    }
    fitEntry8+=", b=";
    fitEntry8+=getTStringFromDouble(function8->GetParameter(1), 5);
    if(plotfiterrors){
      fitEntry8+="#pm";
      fitEntry8+=getTStringFromDouble(function8->GetParError(1) , 5);
    }
    // extract parameters
    double a=function8->GetParameter(0);
    double b=function8->GetParameter(1);
    // turning point
    double min=0;
    double max=500;
    double TP=-a/b;
    // get functions for high, low and central
    TF1* centralErr=new TF1("centralErr",def,min, max);
    centralErr->SetParameter(0, a);
    centralErr->SetParameter(1, b);
    TF1* upErr=new TF1("upErr",def,min, max);
    upErr->SetParameter(0, 2*a);
    upErr->SetParameter(1, 2*b);
    TF1* dnErr=new TF1("upErr",def,min, max);
    dnErr->SetParameter(0, 0.);
    dnErr->SetParameter(1, 0.);
    // - draw errorbands
    upErr->SetFillStyle(1001);
    dnErr->SetFillStyle(1001);
    upErr->SetLineColor(10);
    dnErr->SetLineColor(10);
    upErr->SetFillColor(colorband);
    upErr->SetRange(min,TP);
    upErr->DrawClone("hist same");
    dnErr->SetFillColor(10);
    dnErr->SetLineColor(10);
    dnErr->SetRange(min,TP);
    dnErr->DrawClone("hist same");
    dnErr->SetFillColor(colorband);
    dnErr->SetLineColor(colorband);
    dnErr->SetRange(TP, max);
    dnErr->DrawClone("hist same");
    upErr->SetFillColor(10);
    upErr->SetLineColor(10);
    upErr->SetRange(TP, max);
    upErr->DrawClone("hist same");
    drawLine(TP, 0.35, TP, 1.05, 10, 2, 1);
    // - draw central fit
    centralErr->SetFillStyle(0);
    centralErr->SetFillColor(0);
    centralErr->SetLineColor(kBlue);
    centralErr->SetLineWidth(6);
    centralErr->SetLineColor(fit8color);
    centralErr->SetLineStyle(2);
    centralErr->Draw("hist same");
    // - re-draw 8 TeV data points
    SFdilep8->Draw("p same");
    SFljets8->Draw("p same");
    // legend and labels
    dummy->Draw("axis same");
    TPaveText *label2 = new TPaveText();
    label2 -> SetX1NDC(gStyle->GetPadLeftMargin());
    label2 -> SetY1NDC(1.0-gStyle->GetPadTopMargin());
    label2 -> SetX2NDC(1.0-gStyle->GetPadRightMargin());
    label2 -> SetY2NDC(1.0);
    label2 -> SetTextFont(42);
    TString CMSlab2="";
    if(!PHD) CMSlab2+="CMS Preliminary, ";  
    CMSlab2+="19.7 fb^{-1} at #sqrt{s} = 8 TeV";
    label2->AddText(CMSlab2);
    label2->SetFillStyle(0);
    label2->SetBorderSize(0);
    label2->SetTextSize(0.04);
    label2->SetTextAlign(32);
    label2->Draw("same");
    TLegend *leg3 = new TLegend(0.22, 0.17, 0.85, 0.33);
    leg3->SetFillStyle(0);
    leg3->SetTextSize(0.035);
    leg3->SetBorderSize(0);
    leg3->SetHeader("#font[22]{Parametrisation: exp(a+b#upointx)}");
    TString entryErr=fitEntry8;
    entryErr.ReplaceAll("8 TeV: ", "");
    leg3->AddEntry(centralErr, entryErr , "L");
    leg3->AddEntry(dnErr     , "a,b #pm 100%", "F");
    leg0->Draw("same");
    leg3->Draw("same");
    bcclabel->Draw("same");
    //saving
    if(grey) plotCanvas_[1]->SetGrayscale();
    plotCanvas_[1]->Print("diffXSecFromSignal/plots/combined/2012/xSec/comp"+name+nameext+"7vs8TeVunc.eps");
    plotCanvas_[1]->Print("diffXSecFromSignal/plots/combined/2012/xSec/comp"+name+nameext+"7vs8TeVunc.png");
  }
}


void addTAE(TGraphAsymmErrors* in, TGraphAsymmErrors* out){
  // this functions adds the points of "in" to out

  // calculate new number of points
  int Nin =in ->GetMaxSize();
  int Nout=out->GetMaxSize();
  int Nsum=Nin+Nout;
  out->Expand(Nsum);

  // add points
  for(int p=0; p<Nin; ++p){
    out->SetPoint     ( Nout+p, in->GetX()[p], in->GetY()[p] );
    out->SetPointError( Nout+p, in->GetEXlow()[p], in->GetEXhigh()[p], in->GetEYlow()[p], in->GetEYhigh()[p] );
  }
}

TGraphAsymmErrors* seven(TString quantity, bool BCC){
  // number of bins
  int Nbins7=7;
  if(      quantity=="ttbarMass") Nbins7=7;
  else if (quantity=="ttbarPt"  ) Nbins7=6;
  else if (quantity=="ttbarY"   ) Nbins7=10;
  else if (quantity=="topPt"    ) Nbins7=7;
  else if (quantity=="topY"     ) Nbins7=10;

  // create asymmerrors
  TGraphAsymmErrors* out = new TGraphAsymmErrors(Nbins7);
  
  if(quantity=="ttbarMass"){
    // BCC BIN MC data stat syst tot
    // 345 to   400 & 0.004899  & 0.004812 &   5.2 &   9.7 &  11.1
    // 400 to   470 & 0.004477  & 0.004603 &   5.0 &   8.4 &   9.8
    // 470 to   550 & 0.002516  & 0.002462 &   5.2 &  10.2 &  11.4
    // 550 to   650 & 0.001175  & 0.001144 &   5.6 &  10.6 &  12.0
    // 650 to   800 & 0.000435  & 0.000432 &   6.2 &   8.3 &  10.3
    // 800 to  1100 & 0.000095  & 0.000099 &   7.1 &  20.0 &  21.2
    //1100 to  1600 & 0.000009  & 0.000014 &  13.5 &  19.4 &  23.7
    out->SetPoint( 0, BCC ?  361.5 :  345+0.5*( 400- 345), 0.004812/0.004899 );
    out->SetPoint( 1, BCC ?  435.5 :  400+0.5*( 470- 400), 0.004603/0.004477 );
    out->SetPoint( 2, BCC ?  508.5 :  470+0.5*( 550- 470), 0.002462/0.002516 );
    out->SetPoint( 3, BCC ?  596.5 :  550+0.5*( 650- 550), 0.001144/0.001175 );
    out->SetPoint( 4, BCC ?  715.5 :  650+0.5*( 800- 650), 0.000432/0.000435 );
    out->SetPoint( 5, BCC ?  928.5 :  800+0.5*(1100- 800), 0.000099/0.000095 );
    out->SetPoint( 6, BCC ? 1290.0 : 1100+0.5*(1600-1100), 0.000014/0.000009 );
    out->SetPointError( 0, BCC ? 0. : 0.5*( 400- 345), BCC ? 0. : 0.5*( 400- 345), (0.004812/0.004899)*(11.1/100), (0.004812/0.004899)*(11.1/100) );
    out->SetPointError( 1, BCC ? 0. : 0.5*( 470- 400), BCC ? 0. : 0.5*( 470- 400), (0.004603/0.004477)*( 9.8/100), (0.004603/0.004477)*( 9.8/100) );
    out->SetPointError( 2, BCC ? 0. : 0.5*( 550- 470), BCC ? 0. : 0.5*( 550- 470), (0.002462/0.002516)*(11.4/100), (0.002462/0.002516)*(11.4/100) );
    out->SetPointError( 3, BCC ? 0. : 0.5*( 650- 550), BCC ? 0. : 0.5*( 650- 550), (0.001144/0.001175)*(12.0/100), (0.001144/0.001175)*(12.0/100) );
    out->SetPointError( 4, BCC ? 0. : 0.5*( 800- 650), BCC ? 0. : 0.5*( 800- 650), (0.000432/0.000435)*(10.3/100), (0.000432/0.000435)*(10.3/100) );
    out->SetPointError( 5, BCC ? 0. : 0.5*(1100- 800), BCC ? 0. : 0.5*(1100- 800), (0.000099/0.000095)*(21.2/100), (0.000099/0.000095)*(21.2/100) );
    out->SetPointError( 6, BCC ? 0. : 0.5*(1600-1100), BCC ? 0. : 0.5*(1600-1100), (0.000014/0.000009)*(23.7/100), (0.000014/0.000009)*(23.7/100) );
  }
  else if (quantity=="ttbarPt"  ){
    // BIN MC data stat syst tot
    //   0 to   20 & 0.015187  & 0.015011 &  4.1 &  11.8 &  12.5
    //  20 to   45 & 0.011405  & 0.012083 &  3.5 &   7.0 &   7.8
    //  45 to   75 & 0.005723  & 0.005808 &  3.8 &   9.2 &  10.0
    //  75 to  120 & 0.002706  & 0.002578 &  4.3 &  14.0 &  14.6
    // 120 to  190 & 0.001077  & 0.001032 &  4.5 &   7.8 &   8.9
    // 190 to  300 & 0.000291  & 0.000247 &  6.3 &  18.0 &  19.1
    out->SetPoint( 0, BCC ?   4 :   0+0.5*(20 -  0), 0.015011/0.015187 );
    out->SetPoint( 1, BCC ?  32 :  20+0.5*(45 - 20), 0.012083/0.011405 );
    out->SetPoint( 2, BCC ?  58 :  45+0.5*(75 - 45), 0.005808/0.005723 );
    out->SetPoint( 3, BCC ?  96 :  75+0.5*(120- 75), 0.002578/0.002706 );
    out->SetPoint( 4, BCC ? 152 : 120+0.5*(190-120), 0.001032/0.001077 );
    out->SetPoint( 5, BCC ? 238 : 190+0.5*(300-190), 0.000247/0.000291 );

    out->SetPointError( 0, BCC ? 0. : 0.5*(20 -  0), BCC ? 0. : 0.5*(20 -  0), (0.015011/0.015187)*(12.5/100), (0.015011/0.015187)*(12.5/100) );
    out->SetPointError( 1, BCC ? 0. : 0.5*(45 - 20), BCC ? 0. : 0.5*(45 - 20), (0.012083/0.011405)*( 7.8/100), (0.012083/0.011405)*( 7.8/100) );
    out->SetPointError( 2, BCC ? 0. : 0.5*(75 - 45), BCC ? 0. : 0.5*(75 - 45), (0.005808/0.005723)*(10.0/100), (0.005808/0.005723)*(10.0/100) );
    out->SetPointError( 3, BCC ? 0. : 0.5*(120- 75), BCC ? 0. : 0.5*(120- 75), (0.002578/0.002706)*(14.6/100), (0.002578/0.002706)*(14.6/100) );
    out->SetPointError( 4, BCC ? 0. : 0.5*(190-120), BCC ? 0. : 0.5*(190-120), (0.001032/0.001077)*( 8.9/100), (0.001032/0.001077)*( 8.9/100) );
    out->SetPointError( 5, BCC ? 0. : 0.5*(300-190), BCC ? 0. : 0.5*(300-190), (0.000247/0.000291)*(19.1/100), (0.000247/0.000291)*(19.1/100) );
  }
  else if (quantity=="ttbarY"   ){
    // BIN MC data stat syst tot
    // -2.5 to  -1.3 & 0.058060  & 0.055296 &  6.4 &  10.8 &  12.5  
    // -1.3 to  -0.9 & 0.224558  & 0.217093 &  3.4 &   5.8 &   6.7  
    // -0.9 to  -0.6 & 0.321032  & 0.312239 &  3.6 &   4.4 &   5.7  
    // -0.6 to  -0.3 & 0.388325  & 0.400025 &  3.1 &   3.3 &   4.5  
    // -0.3 to   0.0 & 0.424172  & 0.434507 &  3.1 &   4.1 &   5.1  
    //  0.0 to   0.3 & 0.425370  & 0.468782 &  2.8 &   3.8 &   4.8  
    //  0.3 to   0.6 & 0.388678  & 0.393596 &  3.1 &   5.9 &   6.7  
    //  0.6 to   0.9 & 0.321439  & 0.316880 &  3.4 &   4.7 &   5.8  
    //  0.9 to   1.3 & 0.225316  & 0.222247 &  3.3 &   5.8 &   6.6  
    //  1.3 to   2.5 & 0.057744  & 0.049806 &  6.8 &   9.7 &  11.9  
    out->SetPoint( 0, BCC ? -1.806 : -2.5+0.5*(2.5-1.3), 0.055296/0.058060  );
    out->SetPoint( 1, BCC ? -1.094 : -1.3+0.5*(1.3-0.9), 0.217093/0.224558  );
    out->SetPoint( 2, BCC ? -0.756 : -0.9+0.5*(0.9-0.6), 0.312239/0.321032  );
    out->SetPoint( 3, BCC ? -0.456 : -0.6+0.5*(0.6-0.3), 0.400025/0.388325  );
    out->SetPoint( 4, BCC ? -0.181 : -0.3+0.5*(0.3-0.0), 0.434507/0.424172  );
    out->SetPoint( 5, BCC ?  0.194 :  0.0+0.5*(0.3-0.0), 0.468782/0.425370  );
    out->SetPoint( 6, BCC ?  0.444 :  0.3+0.5*(0.6-0.3), 0.393596/0.388678  );
    out->SetPoint( 7, BCC ?  0.756 :  0.6+0.5*(0.9-0.6), 0.316880/0.321439  );
    out->SetPoint( 8, BCC ?  1.094 :  0.9+0.5*(1.3-0.9), 0.222247/0.225316  );
    out->SetPoint( 9, BCC ?  1.806 :  1.3+0.5*(2.5-1.3), 0.049806/0.057744  );
    out->SetPointError( 0, BCC ? 0. : 0.5*(2.5-1.3), BCC ? 0. : 0.5*(2.5-1.3), (0.055296/0.058060)*(12.5/100), (0.055296/0.058060)*(12.5/100) );
    out->SetPointError( 1, BCC ? 0. : 0.5*(1.3-0.9), BCC ? 0. : 0.5*(1.3-0.9), (0.217093/0.224558)*( 6.7/100), (0.217093/0.224558)*( 6.7/100) );
    out->SetPointError( 2, BCC ? 0. : 0.5*(0.9-0.6), BCC ? 0. : 0.5*(0.9-0.6), (0.312239/0.321032)*( 5.7/100), (0.312239/0.321032)*( 5.7/100) );
    out->SetPointError( 3, BCC ? 0. : 0.5*(0.6-0.3), BCC ? 0. : 0.5*(0.6-0.3), (0.400025/0.388325)*( 4.5/100), (0.400025/0.388325)*( 4.5/100) );
    out->SetPointError( 4, BCC ? 0. : 0.5*(0.3-0.0), BCC ? 0. : 0.5*(0.3-0.0), (0.434507/0.424172)*( 5.1/100), (0.434507/0.424172)*( 5.1/100) );
    out->SetPointError( 5, BCC ? 0. : 0.5*(0.3-0.0), BCC ? 0. : 0.5*(0.3-0.0), (0.468782/0.425370)*( 4.8/100), (0.468782/0.425370)*( 4.8/100) );
    out->SetPointError( 6, BCC ? 0. : 0.5*(0.6-0.3), BCC ? 0. : 0.5*(0.6-0.3), (0.393596/0.388678)*( 6.7/100), (0.393596/0.388678)*( 6.7/100) );
    out->SetPointError( 7, BCC ? 0. : 0.5*(0.9-0.6), BCC ? 0. : 0.5*(0.9-0.6), (0.316880/0.321439)*( 5.8/100), (0.316880/0.321439)*( 5.8/100) );
    out->SetPointError( 8, BCC ? 0. : 0.5*(1.3-0.9), BCC ? 0. : 0.5*(1.3-0.9), (0.222247/0.225316)*( 6.6/100), (0.222247/0.225316)*( 6.6/100) );
    out->SetPointError( 9, BCC ? 0. : 0.5*(2.5-1.3), BCC ? 0. : 0.5*(2.5-1.3), (0.049806/0.057744)*(11.9/100), (0.049806/0.057744)*(11.9/100) );
  }

  else if (quantity=="topPt"    ){
    // BIN MC data stat syst tot
    // 0  to   60 & 0.003806  & 0.004536 &  2.5 &  3.6 &   4.4  
    // 60 to  100 & 0.006574  & 0.006658 &  2.4 &  4.9 &   5.5  
    //100 to  150 & 0.005077  & 0.004740 &  2.4 &  3.2 &   4.0  
    //150 to  200 & 0.002748  & 0.002501 &  2.6 &  5.1 &   5.8  
    //200 to  260 & 0.001195  & 0.001042 &  2.9 &  5.5 &   6.2  
    //260 to  320 & 0.000454  & 0.000378 &  3.7 &  8.2 &   9.0  
    //320 to  400 & 0.000154  & 0.000120 &  5.8 &  9.5 &  11.1  
    out->SetPoint( 0, BCC ?  26.2 :   0+0.5*( 60-  0), 0.004536 / 0.003806 );
    out->SetPoint( 1, BCC ?  88.8 :  60+0.5*(100- 60), 0.006658 / 0.006574 ); 
    out->SetPoint( 2, BCC ? 126.2 : 100+0.5*(150-100), 0.004740 / 0.004740 ); 
    out->SetPoint( 3, BCC ? 173.8 : 150+0.5*(200-150), 0.002501 / 0.002748 ); 
    out->SetPoint( 4, BCC ? 228.8 : 200+0.5*(260-200), 0.001042 / 0.001195 ); 
    out->SetPoint( 5, BCC ? 288.8 : 260+0.5*(320-260), 0.000378 / 0.000454 ); 
    out->SetPoint( 6, BCC ? 356.2 : 320+0.5*(400-320), 0.000120 / 0.000154 ); 
    out->SetPointError( 0, BCC ? 0. : 0.5*( 60-  0), BCC ? 0. : 0.5*( 60-  0), (4.4 /100.)*(0.004536 / 0.003806), (4.4 /100.)*(0.004536 / 0.003806) );
    out->SetPointError( 1, BCC ? 0. : 0.5*(100- 60), BCC ? 0. : 0.5*(100- 60), (5.5 /100.)*(0.006658 / 0.006574), (5.5 /100.)*(0.006658 / 0.006574) );
    out->SetPointError( 2, BCC ? 0. : 0.5*(150-100), BCC ? 0. : 0.5*(150-100), (4.0 /100.)*(0.004740 / 0.004740), (4.0 /100.)*(0.004740 / 0.004740) );
    out->SetPointError( 3, BCC ? 0. : 0.5*(200-150), BCC ? 0. : 0.5*(200-150), (5.8 /100.)*(0.002501 / 0.002748), (5.8 /100.)*(0.002501 / 0.002748) );
    out->SetPointError( 4, BCC ? 0. : 0.5*(260-200), BCC ? 0. : 0.5*(260-200), (6.2 /100.)*(0.001042 / 0.001195), (6.2 /100.)*(0.001042 / 0.001195) );
    out->SetPointError( 5, BCC ? 0. : 0.5*(320-260), BCC ? 0. : 0.5*(320-260), (9.0 /100.)*(0.000378 / 0.000454), (9.0 /100.)*(0.000378 / 0.000454) );
    out->SetPointError( 6, BCC ? 0. : 0.5*(400-320), BCC ? 0. : 0.5*(400-320), (11.1/100.)*(0.000120 / 0.000154), (11.1/100.)*(0.000120 / 0.000154) );
  }
  else if (quantity=="topY"     ){
    // BIN MC data stat syst tot
    // -2.5 to  -1.6 & 0.062394  & 0.065109 &  5.1 &  10.3 &  11.5 
    // -1.6 to  -1.2 & 0.171139  & 0.172594 &  2.9 &   5.9 &   6.6 
    // -1.2 to  -0.8 & 0.252392  & 0.262370 &  2.8 &   4.1 &   5.0 
    // -0.8 to  -0.4 & 0.320539  & 0.316197 &  2.6 &   3.8 &   4.6 
    // -0.4 to   0.0 & 0.358094  & 0.333846 &  2.7 &   4.8 &   5.5 
    //  0.0 to   0.4 & 0.359286  & 0.357691 &  2.5 &   2.6 &   3.6 
    //  0.4 to   0.8 & 0.320713  & 0.327341 &  2.5 &   5.2 &   5.8 
    //  0.8 to   1.2 & 0.252567  & 0.255667 &  2.7 &   5.0 &   5.7 
    //  1.2 to   1.6 & 0.170898  & 0.168034 &  3.0 &   5.7 &   6.4 
    //  1.6 to   2.5 & 0.062206  & 0.064393 &  5.0 &   7.1 &   8.7 
    out->SetPoint( 0, BCC ? -2.012 : -2.5+0.5*(2.5-1.6), 0.065109/0.062394 );
    out->SetPoint( 1, BCC ? -1.387 : -1.6+0.5*(1.6-1.2), 0.172594/0.171139 );
    out->SetPoint( 2, BCC ? -1.012 : -1.2+0.5*(1.2-0.8), 0.262370/0.252392 );
    out->SetPoint( 3, BCC ? -0.613 : -0.8+0.5*(0.8-0.4), 0.316197/0.320539 );
    out->SetPoint( 4, BCC ? -0.237 : -0.4+0.5*(0.4-0.0), 0.333846/0.358094 );
    out->SetPoint( 5, BCC ?  0.213 :  0.0+0.5*(0.4-0.0), 0.357691/0.359286 );
    out->SetPoint( 6, BCC ?  0.613 :  0.4+0.5*(0.8-0.4), 0.327341/0.320713 );
    out->SetPoint( 7, BCC ?  1.012 :  0.8+0.5*(1.2-0.8), 0.255667/0.252567 );
    out->SetPoint( 8, BCC ?  1.387 :  1.2+0.5*(1.6-1.2), 0.168034/0.170898 );
    out->SetPoint( 9, BCC ?  2.013 :  1.6+0.5*(2.5-1.6), 0.064393/0.062206 );
    out->SetPointError( 0, BCC ? 0. : 0.5*(2.5-1.6), BCC ? 0. : 0.5*(2.5-1.6), (0.065109/0.062394)*(11.5/100), (0.065109/0.062394)*(11.5/100) );
    out->SetPointError( 1, BCC ? 0. : 0.5*(1.6-1.2), BCC ? 0. : 0.5*(1.6-1.2), (0.172594/0.171139)*( 6.6/100), (0.172594/0.171139)*( 6.6/100) );
    out->SetPointError( 2, BCC ? 0. : 0.5*(1.2-0.8), BCC ? 0. : 0.5*(1.2-0.8), (0.262370/0.252392)*( 5.0/100), (0.262370/0.252392)*( 5.0/100) );
    out->SetPointError( 3, BCC ? 0. : 0.5*(0.8-0.4), BCC ? 0. : 0.5*(0.8-0.4), (0.316197/0.320539)*( 4.6/100), (0.316197/0.320539)*( 4.6/100) );
    out->SetPointError( 4, BCC ? 0. : 0.5*(0.4-0.0), BCC ? 0. : 0.5*(0.4-0.0), (0.333846/0.358094)*( 5.5/100), (0.333846/0.358094)*( 5.5/100) );
    out->SetPointError( 5, BCC ? 0. : 0.5*(0.4-0.0), BCC ? 0. : 0.5*(0.4-0.0), (0.357691/0.359286)*( 3.6/100), (0.357691/0.359286)*( 3.6/100) );
    out->SetPointError( 6, BCC ? 0. : 0.5*(0.8-0.4), BCC ? 0. : 0.5*(0.8-0.4), (0.327341/0.320713)*( 5.8/100), (0.327341/0.320713)*( 5.8/100) );
    out->SetPointError( 7, BCC ? 0. : 0.5*(1.2-0.8), BCC ? 0. : 0.5*(1.2-0.8), (0.255667/0.252567)*( 5.7/100), (0.255667/0.252567)*( 5.7/100) );
    out->SetPointError( 8, BCC ? 0. : 0.5*(1.6-1.2), BCC ? 0. : 0.5*(1.6-1.2), (0.168034/0.170898)*( 6.4/100), (0.168034/0.170898)*( 6.4/100) );
    out->SetPointError( 9, BCC ? 0. : 0.5*(2.5-1.6), BCC ? 0. : 0.5*(2.5-1.6), (0.064393/0.062206)*( 8.7/100), (0.064393/0.062206)*( 8.7/100) );
  }
  return out;
}

TString DilepFileName(TString name){
  
  if(name=="topPt") return "ToppT";
  else if(name=="ttbarDelPhi"  ) return "TTBarDeltaPhi";
  else if(name=="ttbarMass"    ) return "TTBarMass";
  else if(name=="ttbarY"       ) return "TTBarRapidity";
  else if(name=="ttbarPt"      ) return "TTBarpT";
  else if(name=="topY"         ) return "TopRapidity";
  else if(name=="topPtLead"    ) return "ToppTLead";
  else if(name=="topPtSubLead" ) return "ToppTNLead";
  else if(name=="topPtTtbarSys") return "ToppTTTRestFrame";
  return "UNKNOWN";
}

int NbinsLL(TString name){
  // get number of bins from number of lines in dilepton .tex files
//   TString temp= getStringEntry(gSystem->GetFromPipe("wc -l "+groupSpace+"CommonFiles/topPtInputForReweighting/Hyp"+DilepFileName(name)+"LaTeX.txt"), 1, " ");
  TString temp= getStringEntry(gSystem->GetFromPipe("wc -l /nfs/dust/cms/user/asincruz/CMSSW5314p1_Development/CMSSW_5_3_14_patch1/src/TopAnalysis/Configuration/analysis/diLeptonic/analysis_postCWR/Plots/FinalResults/combined/Hyp"+DilepFileName(name)+"LaTeX.txt"), 1, " ");
  int out=atof(temp.Data())-3;
  return out;
}


std::vector<double> BCCvalues(TString channel, TString name){
  // hardcoded BCC values for the 8TeV results
  std::vector<double> x_;
  if(channel=="ljets"){
    if(name=="topPt"){
      x_.push_back( 26.25 ); //Bin   0.00 ....  60.00 
      x_.push_back( 88.75 ); //Bin  60.00 .... 100.00 
      x_.push_back(126.25 ); //Bin 100.00 .... 150.00 
      x_.push_back(173.75 ); //Bin 150.00 .... 200.00 
      x_.push_back(228.75 ); //Bin 200.00 .... 260.00 
      x_.push_back(286.25 ); //Bin 260.00 .... 320.00 
      x_.push_back(356.25 ); //Bin 320.00 .... 400.00 
      x_.push_back(446.25 ); //Bin 400.00 .... 500.00 
    }
    else if(name=="ttbarDelPhi"  ){
      x_.push_back(1.26);//Bin 0.00 .... 2.00 
      x_.push_back(2.44);//Bin 2.00 .... 2.75 
      x_.push_back(2.89);//Bin 2.75 .... 3.00 
      x_.push_back(3.14);//Bin 3.00 .... 3.15 
    }
    else if(name=="ttbarMass"    ){
      x_.push_back( 362.50);//Bin  345.00 ....  400.00 
      x_.push_back( 435.50);//Bin  400.00 ....  470.00 
      x_.push_back( 508.50);//Bin  470.00 ....  550.00 
      x_.push_back( 595.50);//Bin  550.00 ....  650.00 
      x_.push_back( 717.50);//Bin  650.00 ....  800.00 
      x_.push_back( 927.50);//Bin  800.00 .... 1100.00 
      x_.push_back(1328.50);//Bin 1100.00 .... 1600.00 
    }
    else if(name=="ttbarY"       ){
      x_.push_back(-1.82);// Bin -2.50 .... -1.30 
      x_.push_back(-1.09);// Bin -1.30 .... -0.90 
      x_.push_back(-0.76);// Bin -0.90 .... -0.60 
      x_.push_back(-0.46);// Bin -0.60 .... -0.30 
      x_.push_back(-0.14);// Bin -0.30 ....  0.00 
      x_.push_back( 0.14);// Bin  0.00 ....  0.30 
      x_.push_back( 0.47);// Bin  0.30 ....  0.60 
      x_.push_back( 0.76);// Bin  0.60 ....  0.90 
      x_.push_back( 1.09);// Bin  0.90 ....  1.30 
      x_.push_back( 1.82);// Bin  1.30 ....  2.50 
    }
    else if(name=="ttbarPt"      ){
      x_.push_back(  4.50);//Bin   0.00 ....  20.00
      x_.push_back( 32.50);//Bin  20.00 ....  45.00
      x_.push_back( 58.50);//Bin  45.00 ....  75.00
      x_.push_back( 95.50);//Bin  75.00 .... 120.00
      x_.push_back(152.50);//Bin 120.00 .... 190.00
      x_.push_back(237.50);//Bin 190.00 .... 300.00
    }
    else if(name=="topY"         ){
      x_.push_back(-2.01);//Bin -2.50 .... -1.60 
      x_.push_back(-1.39);//Bin -1.60 .... -1.20 
      x_.push_back(-1.01);//Bin -1.20 .... -0.80 
      x_.push_back(-0.61);//Bin -0.80 .... -0.40 
      x_.push_back(-0.24);//Bin -0.40 ....  0.00 
      x_.push_back( 0.24);//Bin  0.00 ....  0.40 
      x_.push_back( 0.61);//Bin  0.40 ....  0.80 
      x_.push_back( 1.01);//Bin  0.80 ....  1.20 
      x_.push_back( 1.41);//Bin  1.20 ....  1.60 
      x_.push_back( 2.01);//Bin  1.60 ....  2.50 
    }
    else if(name=="topPtLead"    ){
      x_.push_back( 31.25 ); // Bin   0.00 ....  60.00
      x_.push_back( 76.25 ); // Bin  60.00 .... 100.00
      x_.push_back(126.25 ); // Bin 100.00 .... 150.00
      x_.push_back(173.75 ); // Bin 150.00 .... 200.00
      x_.push_back(228.75 ); // Bin 200.00 .... 260.00
      x_.push_back(286.25 ); // Bin 260.00 .... 320.00
      x_.push_back(356.25 ); // Bin 320.00 .... 400.00
      x_.push_back(446.25 ); // Bin 400.00 .... 500.00
    }
    else if(name=="topPtSubLead" ){
      x_.push_back( 23.75 );//Bin   0.00 ....  60.00 
      x_.push_back( 83.75 );//Bin  60.00 .... 100.00 
      x_.push_back(123.75 );//Bin 100.00 .... 150.00 
      x_.push_back(173.75 );//Bin 150.00 .... 200.00 
      x_.push_back(228.75 );//Bin 200.00 .... 260.00 
      x_.push_back(286.25 );//Bin 260.00 .... 320.00 
      x_.push_back(356.25 );//Bin 320.00 .... 400.00 
      x_.push_back(443.75 );//Bin 400.00 .... 500.00 
    }
    else if(name=="topPtTtbarSys"){ 
      x_.push_back( 26.25);//Bin   0.00 ....  60.00 
      x_.push_back( 63.75);//Bin  60.00 .... 100.00 
      x_.push_back(126.25);//Bin 100.00 .... 150.00 
      x_.push_back(173.75);//Bin 150.00 .... 200.00 
      x_.push_back(226.25);//Bin 200.00 .... 260.00 
      x_.push_back(286.25);//Bin 260.00 .... 320.00 
      x_.push_back(356.25);//Bin 320.00 .... 400.00 
      x_.push_back(443.75);//Bin 400.00 .... 500.00 
    }
  }
  else if(channel=="dilepton"){
    if(name=="topPt"){
      x_.push_back( 28.75 ); //Bin   0.00 ....  65.00 
      x_.push_back(101.25 ); //Bin  65.00 .... 125.00 
      x_.push_back(161.25 ); //Bin 125.00 .... 200.00 
      x_.push_back(238.75 ); //Bin 200.00 .... 290.00 
      x_.push_back(336.25 ); //Bin 290.00 .... 400.00 
    }
    else if(name=="ttbarDelPhi"  ){
      x_.push_back(1.19);// Bin 0.00 .... 1.89 
      x_.push_back(2.44);// Bin 1.89 .... 2.77 
      x_.push_back(2.94);// Bin 2.77 .... 3.04 
      x_.push_back(3.09);// Bin 3.04 .... 3.15 
    }
    else if(name=="ttbarMass"    ){
      x_.push_back( 354.50);//Bin  340.00 ....  380.00 
      x_.push_back( 428.50);//Bin  380.00 ....  470.00 
      x_.push_back( 537.50);//Bin  470.00 ....  620.00 
      x_.push_back( 705.50);//Bin  620.00 ....  820.00 
      x_.push_back( 940.50);//Bin  820.00 .... 1100.00 
      x_.push_back(1328.50);//Bin 1100.00 .... 1600.00 
    }
    else if(name=="ttbarY"       ){
      x_.push_back(-1.93);//Bin -2.50 .... -1.50 
      x_.push_back(-1.24);//Bin -1.50 .... -1.00 
      x_.push_back(-0.76);//Bin -1.00 .... -0.50 
      x_.push_back(-0.31);//Bin -0.50 ....  0.00 
      x_.push_back( 0.29);//Bin  0.00 ....  0.50 
      x_.push_back( 0.76);//Bin  0.50 ....  1.00 
      x_.push_back( 1.24);//Bin  1.00 ....  1.50 
      x_.push_back( 1.93);//Bin  1.50 ....  2.50 
    }
    else if(name=="ttbarPt"      ){
      x_.push_back(  4.50);//Bin   0.00 ....  30.00 
      x_.push_back( 51.50);//Bin  30.00 ....  80.00 
      x_.push_back(118.50);//Bin  80.00 .... 170.00 
      x_.push_back(223.50);//Bin 170.00 .... 300.00 
    }
    else if(name=="topY"         ){
      x_.push_back(-2.01);//Bin -2.50 .... -1.60 
      x_.push_back(-1.31);//Bin -1.60 .... -1.00 
      x_.push_back(-0.76);//Bin -1.00 .... -0.50 
      x_.push_back(-0.29);//Bin -0.50 ....  0.00 
      x_.push_back( 0.29);//Bin  0.00 ....  0.50 
      x_.push_back( 0.76);//Bin  0.50 ....  1.00 
      x_.push_back( 1.31);//Bin  1.00 ....  1.60 
      x_.push_back( 2.01);//Bin  1.60 ....  2.50 
    }
    else if(name=="topPtLead"    ){
      x_.push_back( 36.25);//Bin   0.00 ....  75.00
      x_.push_back(111.25);//Bin  75.00 .... 130.00
      x_.push_back(163.75);//Bin 130.00 .... 200.00
      x_.push_back(241.25);//Bin 200.00 .... 290.00
      x_.push_back(336.25);//Bin 290.00 .... 400.00
    }
    else if(name=="topPtSubLead" ){
      x_.push_back( 23.75);// Bin   0.00 ....  55.00 
      x_.push_back( 91.25);// Bin  55.00 .... 120.00 
      x_.push_back(156.25);// Bin 120.00 .... 200.00 
      x_.push_back(238.75);// Bin 200.00 .... 290.00 
      x_.push_back(338.75);// Bin 290.00 .... 400.00 
    }
    else if(name=="topPtTtbarSys"){
      x_.push_back( 26.25);// Bin   0.00 ....  60.00 
      x_.push_back( 93.75);// Bin  60.00 .... 115.00 
      x_.push_back(151.25);// Bin 115.00 .... 190.00 
      x_.push_back(226.25);// Bin 190.00 .... 275.00 
      x_.push_back(318.75);// Bin 275.00 .... 380.00 
      x_.push_back(428.75);// Bin 380.00 .... 500.00 
    }
  }
  return x_;
}
std::vector<double> DilepMCMadGraph(TString name){
  std::vector<double> x_;
  if(name=="topPt"){
    x_.push_back(0.00396076);  //Bin   0.00 ....  65.00 
    x_.push_back(0.00620269);  //Bin  65.00 .... 125.00 
    x_.push_back(0.00336987);  //Bin 125.00 .... 200.00 
    x_.push_back(0.00102834);  //Bin 200.00 .... 290.00 
    x_.push_back(0.000228163); //Bin 290.00 .... 400.00 
  }
  else if(name=="ttbarDelPhi"  ){
    x_.push_back(0.0607185);// Bin 0.00 .... 1.89 
    x_.push_back(0.283538 );// Bin 1.89 .... 2.77 
    x_.push_back(1.24641  );// Bin 2.77 .... 3.04 
    x_.push_back(2.74429  );// Bin 3.04 .... 3.15 
  }
  else if(name=="ttbarMass"    ){
    x_.push_back(0.00375959 );//Bin  340.00 ....  380.00 
    x_.push_back(0.00464199 );//Bin  380.00 ....  470.00 
    x_.push_back(0.0019876  );//Bin  470.00 ....  620.00 
    x_.push_back(0.000505343);//Bin  620.00 ....  820.00 
    x_.push_back(9.67355e-05);//Bin  820.00 .... 1100.00 
    x_.push_back(1.10849e-05);//Bin 1100.00 .... 1600.00 
  }
  else if(name=="ttbarY"       ){
    x_.push_back(0.048721 );//Bin -2.50 .... -1.50 
    x_.push_back(0.189206 );//Bin -1.50 .... -1.00 
    x_.push_back(0.313871 );//Bin -1.00 .... -0.50 
    x_.push_back(0.398119 );//Bin -0.50 ....  0.00 
    x_.push_back(0.398312 );//Bin  0.00 ....  0.50 
    x_.push_back(0.314823 );//Bin  0.50 ....  1.00 
    x_.push_back(0.190316 );//Bin  1.00 ....  1.50 
    x_.push_back(0.0489551);//Bin  1.50 ....  2.50 
  }
  else if(name=="ttbarPt"      ){
    x_.push_back(0.0142172  );//Bin   0.00 ....  30.00 
    x_.push_back(0.00702319 );//Bin  30.00 ....  80.00 
    x_.push_back(0.00192234 );//Bin  80.00 .... 170.00 
    x_.push_back(0.000379327);//Bin 170.00 .... 300.00 
  }
  else if(name=="topY"         ){
    x_.push_back(0.0698754);//Bin -2.50 .... -1.60 
    x_.push_back(0.195331 );//Bin -1.60 .... -1.00 
    x_.push_back(0.292762 );//Bin -1.00 .... -0.50 
    x_.push_back(0.345366 );//Bin -0.50 ....  0.00 
    x_.push_back(0.346687 );//Bin  0.00 ....  0.50 
    x_.push_back(0.293418 );//Bin  0.50 ....  1.00 
    x_.push_back(0.196023 );//Bin  1.00 ....  1.60 
    x_.push_back(0.0702028);//Bin  1.60 ....  2.50 
  }
  else if(name=="topPtLead"    ){
    x_.push_back(0.0030416  );//Bin   0.00 ....  75.00
    x_.push_back(0.00618785 );//Bin  75.00 .... 130.00
    x_.push_back(0.00391269 );//Bin 130.00 .... 200.00
    x_.push_back(0.00136176 );//Bin 200.00 .... 290.00
    x_.push_back(0.000319097);//Bin 290.00 .... 400.00
  }
  else if(name=="topPtSubLead" ){
    x_.push_back(0.00492779 );// Bin   0.00 ....  55.00 
    x_.push_back(0.00652115 );// Bin  55.00 .... 120.00 
    x_.push_back(0.00284016 );// Bin 120.00 .... 200.00 
    x_.push_back(0.000696973);// Bin 200.00 .... 290.00 
    x_.push_back(0.000137789);// Bin 290.00 .... 400.00 
  }
  else if(name=="topPtTtbarSys"){
    x_.push_back(0.00403399 );// Bin   0.00 ....  60.00 
    x_.push_back(0.00662169 );// Bin  60.00 .... 115.00 
    x_.push_back(0.00365223 );// Bin 115.00 .... 190.00 
    x_.push_back(0.00106479 );// Bin 190.00 .... 275.00 
    x_.push_back(0.000229744);// Bin 275.00 .... 380.00 
    x_.push_back(4.35003e-05);// Bin 380.00 .... 500.00 
  }
  return x_;
}
