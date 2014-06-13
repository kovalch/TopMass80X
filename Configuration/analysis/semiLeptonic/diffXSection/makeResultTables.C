#include "basicFunctions.h"
#include <stdlib.h>

void makeResultTables(std::string decayChannel = "combined", bool extrapolate=false, bool hadron=true, bool addCrossCheckVariables=false, bool useBCC=false, int verbose=1){
  
  // ============================
  //  Set Root Style
  // ============================

  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  TGaxis::SetMaxDigits(2);
  myStyle.cd();
  gROOT->SetStyle("HHStyle");

  // ============================
  //  Open file
  // ============================
  if(extrapolate==true) hadron=false; 
  TString filename="diffXSecTopSemi";
  if(decayChannel=="combined") filename+="Lep";
  else if(decayChannel=="electron") filename+="Elec";
  else if(decayChannel=="muon"    ) filename+="Mu";
  if(extrapolate) filename+="Parton";
  else{
    if(hadron) filename+="Hadron";
    else filename+="Parton";
    filename+="PhaseSpace";
  }
  filename+=".root";
  if(verbose>0) std::cout << "opening file " << filename << std::endl;
  TFile* file = TFile::Open(filename, "READ");
  if(!file){
    std::cout << "ERROR: can not open file " << filename << std::endl;
    exit(0);
  }
  // dont associate new objects with file to be able to close it in the end
  gROOT->cd();
  // ============================
  //  Get plots
  // ============================
  // variables to be processed
  std::vector<TString> xSecVariables_;
  // a) top and ttbar quantities
  if(!hadron){
    xSecVariables_.insert(xSecVariables_.end(), xSecVariablesKinFit, xSecVariablesKinFit + sizeof(xSecVariablesKinFit)/sizeof(TString));
  }
  // b) lepton and b-jet quantities
  if(hadron||!extrapolate){
    xSecVariables_.insert(xSecVariables_.end(), xSecVariablesFinalState    , xSecVariablesFinalState     + sizeof(xSecVariablesFinalState    )/sizeof(TString));
  }
  // c) cross check variables presently only available for parton level cross-sections
  if (addCrossCheckVariables && !hadron){
    xSecVariables_.insert( xSecVariables_.end(),   xSecVariablesCCVar,     xSecVariablesCCVar     + sizeof(xSecVariablesCCVar    )/sizeof(TString)    );
    xSecVariables_.insert( xSecVariables_.end(),   xSecVariablesCCVarNorm, xSecVariablesCCVarNorm + sizeof(xSecVariablesCCVarNorm)/sizeof(TString));
  }

  for(unsigned int i=0; i<xSecVariables_.size(); ++i){
    TString plotName=xSecVariables_[i];
    if(verbose>0) std::cout << std::endl << "variable: " << plotName << std::endl;
    // get canvas for chosen cross section
    TCanvas* canvas = (TCanvas*)(file->Get("finalXSec/"+plotName+"Norm")->Clone());
    if(!canvas){
      std::cout << "ERROR: can not load canvas finalXSec/"+plotName+"Norm" << std::endl;
      exit(0);
    }
    // GET DATA: with final errors from canvas
    TGraphAsymmErrors* dataTot  = (TGraphAsymmErrors*)canvas->GetPrimitive("dataTotalError");
    TGraphAsymmErrors* dataStat = (TGraphAsymmErrors*)canvas->GetPrimitive("dataStatError" );
    TH1F* binned = (TH1F*)canvas->GetPrimitive(plotName);
    TH1F* binnedMCatNLO = (TH1F*)canvas->GetPrimitive(plotName+"MC@NLO");
    TH1F* binnedPowheg  = (TH1F*)canvas->GetPrimitive(plotName+"POWHEG");
    TH1F* binnedPowhegHerwig  = (TH1F*)canvas->GetPrimitive(plotName+"POWHEGHERWIG");
    TH1F* binnedNNLO    = (TH1F*)canvas->GetPrimitive(plotName+"nnlo");
    if(!dataTot){
      std::cout << "ERROR: can not load TGraphAsymmErrors dataTotalError in canvas finalXSec/"+plotName+"Norm" << std::endl;
      exit(0);
    }
    if(!dataStat){
      std::cout << "ERROR: can not load TGraphAsymmErrors dataStatError in canvas finalXSec/"+plotName+"Norm" << std::endl;
      exit(0);
    }
    if(!binned){
      std::cout << "ERROR: can not load TH1F topPt in canvas finalXSec/"+plotName << std::endl;
      exit(0);
    }
    // define range of relevant plots
    // INFO: keep this consistent with the range as defined in setXAxisRange
    //       and makevariableBinning in basicFunctions.h
    double xMin=-999999999;
    double xMax= 999999999;
    if(plotName.Contains     ("topPt"    )){ xMin=0.   ; xMax=401. ;} 
    else if(plotName.Contains("topY"     )){ xMin=-2.51; xMax=2.51 ;}
    else if(plotName.Contains("ttbarY"   )){ xMin=-2.51; xMax=2.51 ;}
    else if(plotName.Contains("ttbarMass")){ xMin=344. ; xMax=1601.;}
    else if(plotName.Contains("ttbarPt"  )){ xMin=0.   ; xMax=301. ;}
    else if(plotName.Contains("lepPt"    )){ xMin=29   ; xMax=201. ;}
    else if(plotName.Contains("lepEta"   )){ xMin=-2.11; xMax=2.11 ;}
    else if(plotName.Contains("bqPt"     )){ xMin=29.  ; xMax=401. ;}
    else if(plotName.Contains("bqEta"    )){ xMin=-2.41; xMax=2.41 ;}
    else if(plotName.Contains("bbbarPt"  )){ xMin=0.   ; xMax=800. ;}
    else if(plotName.Contains("bbbarMass")){ xMin=0.   ; xMax=1200.;}
    else if(plotName.Contains("ttbarDelPhi" )){ xMin=0. ; xMax=3.16;}
    else if(plotName.Contains("ttbarPhiStar")){ xMin=0. ; xMax=2.01;}
    else if(plotName.Contains("lbMass"      )){ xMin=0. ; xMax=501.;}
    else if(plotName.Contains("Njets"       )){ xMin=3. ; xMax=10. ;}
    else if(plotName.Contains("rhos"        )){ xMin=0. ; xMax=1.1 ;}

    // initialize ndof counter
    int ndof=0;    
    // initialize global chi2
    double chi2=0;
    double chi2Mc=0;
    double chi2Po=0;
    double chi2PoHer=0;
    double chi2NN=0;
    //  loop all bins
    for(int bin=1; bin<=binned->GetNbinsX(); ++bin){
      if(verbose>1) std::cout << "bin #" << bin;
      // collect information
      double MCxSec  =binned       ->GetBinContent(bin);
      double MCxSecMc=binnedMCatNLO ? binnedMCatNLO->GetBinContent(bin) : 0;
      double MCxSecPo=binnedPowheg  ? binnedPowheg ->GetBinContent(bin) : 0;
      double MCxSecPoHer=binnedPowhegHerwig ? binnedPowhegHerwig ->GetBinContent(bin) : 0;
      double MCxSecNN=binnedNNLO    ? binnedNNLO   ->GetBinContent(bin) : 0;
      // FIXME: current topY NNLO prediction is shifted by one! make sure this is still the case if you update the new prediction
      if(plotName.Contains("topY")) MCxSecNN=binnedNNLO ? binnedNNLO->GetBinContent(bin+1) : 0;
      double xSec=dataTot->GetY()[bin];
      double totError=dataTot->GetErrorYhigh(bin);
      double statError=dataStat->GetErrorYhigh(bin);
      double sysError=sqrt(totError*totError-statError*statError);
      double BCCxValue=dataTot->GetX()[bin];
      double xValueUp=binned->GetBinLowEdge(bin+1);
      double xValueDn=binned->GetBinLowEdge(bin);      
      if(verbose>1) std::cout << std::setprecision(2) << std::fixed << ", xvalue: " << BCCxValue << " (" << xValueDn << ".." << xValueUp << ")" << std::endl;
      // combine information in Latex line style in one TString
      int precXSec=6;
      int precErr=1;
      int precXBCC=2;
      int precX=1;
      if(plotName.Contains("Pt")){
	precXBCC=1;
	if(plotName.Contains("Lep"  )) precXBCC=2;
	if(plotName.Contains("ttbar")) precXBCC=0;
	precX=0;	
      }
      if(plotName.Contains("Eta")||plotName.Contains("Y")){
	precXBCC=3;
	precX=1;
      }
      if(plotName.Contains("Mass")){
	precXBCC=1;
	precX=0;
      }
      if(plotName.Contains("ttbarDelPhi")||plotName.Contains("rhos")) precX=2;


//       if(plotName.Contains("ttbarY")){
// 	std::cout << std::endl << "xValueDn=" << xValueDn << std::endl;
// 	std::cout << "precX=" << precX << std::endl;
// 	std::cout << "xValueDn +5./(pow(10,precX+1))=" << xValueDn +5./(pow(10,precX   +1)) << std::endl;
// 	TString help=getTStringFromDouble(xValueDn +5./(pow(10,precX+1)), precX, true);
// 	std::cout << "rounded number=" << help << std::endl;
//       }
      TString out= "";
      if(useBCC){
	out+=getTStringFromDouble(BCCxValue, precXBCC);
	out+=" & ";
      }
      out+=fillspace(xValueDn, getBigitFromDouble(binned->GetBinLowEdge(binned->GetNbinsX()+1)));	  
      out+=getTStringFromDouble(xValueDn, precX);
      out+=" to ";	
      out+=fillspace(xValueUp, getBigitFromDouble(binned->GetBinLowEdge(binned->GetNbinsX()+1)));
      out+=getTStringFromDouble(xValueUp, precX);
      out+=" & ";
      out+=getTStringFromDouble(MCxSec, precXSec);
      out+=" & ";
      out+=getTStringFromDouble(xSec  , precXSec);
      out+=" & ";
      out+=fillspace(100*(statError/xSec), 2);
      out+=getTStringFromDouble(100*(statError/xSec),  precErr);
      out+=" & ";
      out+=fillspace(100*(sysError/xSec), 2);
      out+=getTStringFromDouble(100*(sysError/xSec ),  precErr);
      out+=" & ";	
      out+=fillspace(100*(totError/xSec), 2);		    
      out+=getTStringFromDouble(100*(totError/xSec ),  precErr);
      out+=" \\\\ ";
      bool append= (bin==1 ? false : true);
      TString txtfile="./diffXSecFromSignal/plots/"+TString(decayChannel)+"/2012/"+filename;
      txtfile.ReplaceAll(".root", plotName+".txt");
      TString txtfile2=txtfile;
      txtfile2.ReplaceAll(plotName, plotName+"chi2");
      writeToFile(out, txtfile, append);
      // chi2 for this distribution
      if(xValueDn>=xMin&&xValueUp<=xMax){
	++ndof;	
	chi2+=                 ((std::abs(MCxSec  -xSec)/totError)*(std::abs(MCxSec  -xSec)/totError));
	if(MCxSecMc!=0)chi2Mc+=((std::abs(MCxSecMc-xSec)/totError)*(std::abs(MCxSecMc-xSec)/totError));
	if(MCxSecPo!=0)chi2Po+=((std::abs(MCxSecPo-xSec)/totError)*(std::abs(MCxSecPo-xSec)/totError));
	if(MCxSecPoHer!=0)chi2PoHer+=((std::abs(MCxSecPoHer-xSec)/totError)*(std::abs(MCxSecPoHer-xSec)/totError));
	if(MCxSecNN!=0)chi2NN+=((std::abs(MCxSecNN-xSec)/totError)*(std::abs(MCxSecNN-xSec)/totError));
	if(verbose>1) std::cout << "-> considered for chi2" << std::endl;
      }
      //std::cout << out << std::endl;
      //std::cout << BCCxValue << " &  " << xValueDn << " to  " << xValueUp << " & " << MCxSec << "  & " << xSec << " &  " << statError/xSec << " &  " << sysError/xSec << " &  " << totError/xSec << " \\\\ " << std::endl;
      if(verbose>1){
	std::cout << std::setprecision(7) << std::fixed << "data:     " << xSec << "+/-" << statError << "+/-" << sysError << std::endl;
	std::cout << std::setprecision(7) << std::fixed << "MadGraph+Pythia: " << MCxSec; 
	std::cout << std::setprecision(2) << std::fixed << " (" << std::abs(MCxSec  -xSec)/totError << " std variations)" << std::endl;
	if(MCxSecMc!=0){
	std::cout << std::setprecision(7) << std::fixed << "MC@NLO+Herwig:   " << MCxSecMc;
	std::cout << std::setprecision(2) << std::fixed << " (" << std::abs(MCxSecMc-xSec)/totError << " std variations)" << std::endl;
	}
	if(MCxSecPo!=0){
	  std::cout << std::setprecision(7) << std::fixed << "Powheg+Pythia:   " << MCxSecPo;
	  std::cout << std::setprecision(2) << std::fixed << " (" << std::abs(MCxSecPo-xSec)/totError << " std variations)" << std::endl;
	}
	if(MCxSecPoHer!=0){
	  std::cout << std::setprecision(7) << std::fixed << "Powheg+Herwig:   " << MCxSecPoHer;
	  std::cout << std::setprecision(2) << std::fixed << " (" << std::abs(MCxSecPoHer-xSec)/totError << " std variations)" << std::endl;
	}
	if(MCxSecNN!=0){
	   std::cout << std::setprecision(7) << std::fixed << "NNLO:     " << MCxSecNN;
	   std::cout << std::setprecision(2) << std::fixed << " (" << std::abs(MCxSecNN-xSec)/totError << " std variations)" << std::endl;
	}
      }
      int precchi=2;
      if(bin==binned->GetNbinsX()&&ndof!=0){
	TString out2=" ";//"$"+xSecLabelName(plotName)+"$";
	//out2.ReplaceAll("GeV", "\\GeV");
	//out2.ReplaceAll("p_{T}^{", "{\\rm p}_{\\rm T}^{\\rm ");
	//out2.ReplaceAll("y^{", "{\\rm y}^{\\rm ");
	//out2.ReplaceAll("m^{", "m^{\\rm ");
	//out2.ReplaceAll("#eta^{", "#eta^{\\rm ");
	//out2.ReplaceAll("N_{jets} (p_{T}>30 GeV, |#eta|<2.4)", "{\\rm N}_{\\rm jets}");
	//out2.ReplaceAll("#rho_{S}= #frac{#scale[0.55]{2#upoint170 GeV}}{#scale[0.55]{m(t#bar{t}+1jet)}}", "#rho_{\\rm S}");
	//out2.ReplaceAll("#", "\\");

	ndof--;// normalisation eats up one dof
	chi2  /=ndof;
	chi2Mc/=ndof;
	chi2Po/=ndof;
	chi2PoHer/=ndof;
	chi2NN/=ndof;
	if(verbose>1) std::cout << std::endl;
	if(chi2  !=0   ){
	  writeToFile("%chi2(MadGraph+Pythia): "+getTStringFromDouble(chi2), txtfile, true); 
	  if(verbose>0){std::cout << "chi2(MadGraph+Pythia): " << chi2      << std::endl;}
	  out2+=" & "+getTStringFromDouble(chi2, precchi);
	}
	else out2+=" & - ";
	if(chi2Mc!=0   ){
	  writeToFile("%chi2(MC@NLO+Herwig  ): "+getTStringFromDouble(chi2Mc   ), txtfile, true); 
	  if(verbose>0){std::cout << "chi2(MC@NLO+Herwig  ): " << chi2Mc    << std::endl;}
	   out2+=" & "+getTStringFromDouble(chi2Mc, precchi);
	}
	else out2+=" & - ";
	if(chi2Po!=0   ){
	  writeToFile("%chi2(Powheg+Pythia  ): "+getTStringFromDouble(chi2Po   ), txtfile, true); 
	  if(verbose>0){std::cout << "chi2(Powheg+Pythia  ): " << chi2Po    << std::endl;}
	   out2+=" & "+getTStringFromDouble(chi2Po, precchi);
	}
	else out2+=" & - ";
	if(chi2PoHer!=0){
	  writeToFile("%chi2(Powheg+Herwig  ): "+getTStringFromDouble(chi2PoHer), txtfile, true); 
	  if(verbose>0){std::cout << "chi2(Powheg+Herwig  ): " << chi2PoHer << std::endl;}
	  out2+=" & "+getTStringFromDouble(chi2PoHer, precchi);
	}
	//else out2+=" & - ";
	if(chi2NN!=0   ){
	  writeToFile("%chi2(NNLO/NLO+NNLL  ): "+getTStringFromDouble(chi2NN   ), txtfile, true); 
	  if(verbose>0){std::cout << "chi2(NNLO/NLO+NNLL  ): " << chi2NN    << std::endl;}
	   out2+=" & "+getTStringFromDouble(chi2NN, precchi);
	}
	else{
	  if(extrapolate) out2+=" & - ";
	}
	out2+=" \\\\ ";	
	writeToFile(out2, txtfile2, false); 
      }
    }
  }
}
