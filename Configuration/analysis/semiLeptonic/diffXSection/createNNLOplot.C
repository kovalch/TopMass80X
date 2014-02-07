#include "basicFunctions.h"

void createNNLOplot(TString theory="ahrens")
{
  // theory="ahrens", "kidonakis"
  TString outputRootFile="test.root";
  // NB: add new datset name here
  if(theory.Contains("ahrens")   ){
    outputRootFile="AhrensNLONNLL";
    //if(theory.Contains("mtt")    ) outputRootFile+="mttbar" ;
    //else if(theory.Contains("pt")) outputRootFile+="pTttbar";
    outputRootFile+=".root";
  }
  else if(theory=="kidonakis") outputRootFile="KidonakisAproxNNLO.root";
  // general options
  gStyle->SetOptStat(0);
  bool usequad=true;
  bool divideByBinwidth=true;
  // list of variables
  std::vector<TString> xSecVariables_, xSecLabel_;
  // NB: add variables for new datset name here
  TString xSecVariablesUsedAhrens[]    ={"ttbarMass", "ttbarPt"};
  TString xSecVariablesUsedKidonakis[] ={"topPt"    , "topY"   };
  if(     theory.Contains("ahrens")   ) xSecVariables_ .insert( xSecVariables_.begin(), xSecVariablesUsedAhrens   , xSecVariablesUsedAhrens    + sizeof(xSecVariablesUsedAhrens   )/sizeof(TString) );
  else if(theory=="kidonakis") xSecVariables_ .insert( xSecVariables_.begin(), xSecVariablesUsedKidonakis, xSecVariablesUsedKidonakis + sizeof(xSecVariablesUsedKidonakis)/sizeof(TString) );
  // get variable binning used for final cross sections 
  std::map<TString, std::vector<double> > bins_=makeVariableBinning();
  //std::vector<double> tempbins_;
  //double ttbarMassBins[]={345,445,545,745, 1045};
  //tempbins_.insert( tempbins_.begin(), ttbarMassBins, ttbarMassBins + sizeof(ttbarMassBins)/sizeof(double) );
  //bins_["ttbarMass"]=tempbins_;

  // loop all variables
  for(unsigned var=0; var<xSecVariables_.size(); ++var){
    TString variable=xSecVariables_[var];
    std::cout << "----------" << std::endl;
    std::cout << theory << ": " << variable << std::endl;
    // get bin boundaries
    double low =bins_[variable][0];
    double high=bins_[variable][bins_[variable].size()-1];
    // --------------------------
    // get raw nnlo theory points
    // --------------------------
    // NB: add new datset file names here
    TGraph * rawHist;
    TString predictionfolder="/afs/naf.desy.de/group/cms/scratch/tophh/CommonFiles/";
    if(theory=="ahrens"   ){
      if(variable.Contains("ttbarMass")) rawHist= new TGraph(predictionfolder+"AhrensTheory_Mttbar_8000_172.5_NLONNLL_norm.dat");//AhrensTheory_Mtt_7000_172.5_Mtt_fin.dat
      if(variable.Contains("ttbarPt"  )) rawHist= new TGraph(predictionfolder+"AhrensTheory_pTttbar_8000_172.5_NLONNLL_abs.dat");//AhrensTheory_pTttbar_7000_172.5_NLONNLL_abs.dat
    }
    else if(theory=="kidonakis"){
      if(variable.Contains("topPt")) rawHist= new TGraph(predictionfolder+"pttopnnloapprox8lhc173m.dat");//"ptnormalNNLO7lhc173m.dat" //"pttopnnloapprox8lhc173m.dat"
      if(variable.Contains("topY" )) rawHist= new TGraph(predictionfolder+"ytopnnloapprox8lhc173m.dat" );//"ynormalNNLO7lhc173m.dat"  //"ytopnnloapprox8lhc173m.dat"
    }
    // NB: say if points should be interpreted as single points with 
    //     nothing in between or as integrated value over the range
    bool points=true;
    if(theory.Contains("ahrens")) points=false;
    else if(theory=="kidonakis")  points=true;
    std::cout << "input: " << rawHist->GetTitle() << std::endl;
    std::cout << "interprete values as points?: " << points << std::endl;

    // --------------------
    // convertion to TH1F*
    // --------------------
    double *Xvalues = rawHist->GetX();
    int Nbins=rawHist->GetN();
    //double *Yvalues = rawHist->GetY(); // not working
    double xmax=-1;
    double xmin=-1;
    double binwidth=-1;
    //
    TH1F* hist;
    // NB: add additional binning for new theory here
    // NB: if loaded data should be interpreted as points with 
    //     nothing in between (like kidonakis), make suree you 
    //     choose a binning that is fine enough for the 
    //     data points loaded
    if(theory=="ahrens"){
      if(variable.Contains("ttbarMass")){
	xmin= 345.;
	xmax=2720.;
	binwidth=25.;// 5 for 8TeV, 25 for 7TeV
	if(TString(rawHist->GetTitle()).Contains("8000")) binwidth=5.;
      }
      else if(variable.Contains("ttbarPt")){
	xmin= 0.;
	xmax=1300.;
	binwidth=5.;
      }
    }
    else if(theory=="kidonakis"){ 
      if(variable.Contains("topPt")){
	xmin=   0.;
	xmax=1500.;
	binwidth=1.;
      }
      else if(variable.Contains("topY")){
	xmin=-3.8;
	xmax= 3.8;
	binwidth=0.01;
      }
    }
    // fill data in binned TH1F
    hist= new TH1F ( variable, variable, (int)((xmax-xmin)/binwidth), xmin, xmax);
    TH1F* ori=(TH1F*)hist->Clone("original points");
    std::cout << "fine binned theory prediction has " << hist->GetNbinsX() << " bins" << std::endl;
    std::cout << "loaded values from .dat file: " << std::endl;
    // list all data values loaded and the corresponding bin
    for(int bin=1; bin<=Nbins; ++bin){
      double x=-999;
      double y=-999;
      // NB: choose if loaded data is interpreted as points with nothing 
      //     between (like kidonakis) or as integrated over the bin range (like ahrens)
      if(points){
	// check if you are still inside the array
	//std::cout << "data point " << bin-1 << "/" << sizeof(Xvalues)/sizeof(double) << std::endl;
	if(rawHist->GetPoint(bin-1, x, y)!=-1){
	  //x=Xvalues[bin-1]; // get value from data points
	  x+=0.5*binwidth;  // add half the binwidth to get the center of the bin
	}
      }
      else{
	x=hist->GetBinCenter(bin);  // get bin center
	y=rawHist->GetY()[bin-1];
      }
      if(x!=-999){
	std::cout << "data point: " << bin;
	std::cout << "(<x>=" << x << ")-> bin";
	// get bin in target (fine binned) plot
	int bin2=bin;
	if(points) bin2=hist->FindBin(x);
	//double y=Yvalues[bin2-1];
	std::cout << bin2 << " (";
	std::cout << hist->GetBinLowEdge(bin2) << ".." << hist->GetBinLowEdge(bin2+1);
	std::cout << "): " << y << std::endl;
	// fill target (fine binned) plot
	if(!points) hist->SetBinContent(bin2, y);
	// mark bins coming from the original prediction
	ori->SetBinContent(bin2, 1.);
	// -------------------------------
	// fit range without data entries
	// -------------------------------
	// NB: needed if loaded data is interpreted as points with nothing 
	//     between (like kidonakis)
	if(points){
	  // perform a linear fit wrt previous point
	  // get the two points (this bin and the previous one)
	  int binPrev= (bin==1&&variable=="topPt") ? 0 : hist->FindBin(Xvalues[bin-2]+0.5*binwidth);
	  double x2=-1;
	  double x1= 0;
	  double y2=-1;
	  double y1= 0.;
	  rawHist->GetPoint(bin-1, x2, y2);
	  x2+=0.5*binwidth;
	  if(bin==1&&variable=="topPt"){
	    y1=0;
	    x1=0;
	  }
	  else{
	    rawHist->GetPoint(bin-2, x1, y1);
	    x1+=0.5*binwidth;
	  }
	  // calculate linear funtion
	  double a=(y2-y1)/(x2-x1);
	  double b=y1-a*x1;
	  TF1* linInterpol=new TF1("linInterpol"+getTStringFromInt(bin), "[0]*x+[1]", x1, x2);
	  linInterpol->SetParameter(0,a);
	  linInterpol->SetParameter(1,b);
	  double xlow  = (bin==1&&variable=="topPt") ? 0. : hist->GetBinLowEdge(binPrev+1);
	  double xhigh = hist->GetBinLowEdge(bin2+1);
	  std::cout << " lin. interpolation [ (" << x1 << "," << y1 << ") .. (" << x2 << "," << y2 << ") ]: " << a << "*x+" << b << std::endl;
	  linInterpol->SetRange(xlow, xhigh);
	  hist->Add(linInterpol);
	}
      }
    }
    // check theory curve
    std::cout << std::endl << "analyze theory curve:" << std::endl;
    double integralRawTheory=hist->Integral(0.,hist->GetNbinsX()+1);
    std::cout << "Integral: " << integralRawTheory << std::endl;
    if(integralRawTheory==0.){
      std::cout << "ERROR: Integral can not be 0!" << std::endl;
      exit(0);
    }
    std::cout << "binwidth: " << binwidth << std::endl;
    std::cout << "Integral*binwidth: " << integralRawTheory*binwidth << std::endl;
    std::cout << "binRange: " << xmin << ".." << xmax << std::endl;
    std::cout << "   -> reco range: " << low << ".." << high << std::endl;
    std::vector<double> recoBins_=bins_[variable];
    for(unsigned int i=0; i<recoBins_.size(); ++i){
      i==0 ? std::cout << "(" : std::cout << ", ";
      std::cout << recoBins_[i];
      if(i==recoBins_.size()-1) std::cout << ")" << std::endl;
    }
    // check if you need to divide by binwidth
    if(std::abs(1.-integralRawTheory)<0.03){
      std::cout << "Integral is approx. 1 -> need to divide by binwidth!" << std::endl;
      hist->Scale(1./binwidth);
    }
    else if(std::abs(1-integralRawTheory*binwidth)<0.03){
      std::cout << "Integral*binwidth is approx. 1 -> no scaling needed!" << std::endl;
    }
    else{
      std::cout << "need to normalize and divide by binwidth";
      hist->Scale(1./(binwidth*integralRawTheory));
    }
    // styling
    hist->GetXaxis()->SetRangeUser(low,high);
    hist->SetMarkerColor(kMagenta);
    hist->SetLineColor(  kMagenta);
    hist->SetMarkerStyle(24);
    // create canvas
    std::cout << std::endl << "create canvas " << std::endl;
    TCanvas *canv = new TCanvas(variable,variable,800,600);
    canv->cd();
    std::cout << "draw original theory curve " << std::endl;
    //temp->Draw("axis");
    hist->Draw("p");
    //hist->Draw("hist same");
    // draw original points
    if(points){
      for(int bin=1; bin<ori->GetNbinsX(); ++bin){
	// create copy of original data points with the normalized values
	if(ori->GetBinContent(bin)!=0) ori->SetBinContent(bin, hist->GetBinContent(bin));
	ori->SetMarkerColor(kBlack);
	ori->SetMarkerStyle(29);
	ori->SetMarkerSize(1);
	ori->Draw("p same");
      }
    }
    // --------------------
    // create rebinned plot
    // --------------------
    std::cout << std::endl << "create rebinned curve:" << std::endl;
    TString name="";
    // NB: add name for dataset here
    name=variable;
    if(    theory=="kidonakis") name+="approxnnlo";
    else if(theory=="ahrens"  ) name+="nlonnll"   ;
    TH1F* binnedPlot=new TH1F(name, name, bins_[variable].size()-1, &bins_[variable][0]);
    for(int bin=1; bin<=hist->GetNbinsX(); ++bin){
      double y=hist->GetBinContent(bin)*hist->GetBinWidth(bin);
      double xlow =hist->GetBinLowEdge(bin);
      double xhigh=hist->GetBinLowEdge(bin+1);
      // search corresponding bin in rebinned curve
      bool found=false;
      //std::cout << "xlow: " << xlow << ", xhigh: " << xhigh << std::endl; // FIXME
      for(int binRebinned=0; binRebinned<=binnedPlot->GetNbinsX()+1; ++binRebinned){
	if(binnedPlot->GetBinLowEdge(binRebinned)<=xlow&&binnedPlot->GetBinLowEdge(binRebinned+1)>=xhigh){
	  found=true;
	  binnedPlot->SetBinContent(binRebinned, binnedPlot->GetBinContent(binRebinned)+y);
	  break;
	}
	//else{
	//  std::cout << "not in: " << binnedPlot->GetBinLowEdge(binRebinned) << ".." << binnedPlot->GetBinLowEdge(binRebinned+1)<< std::endl;
	//}
      }
      // --------------------
      // use linear interpolation for edge bins
      // --------------------
      if(hist->GetBinCenter(bin)<high&&!found){
	std::cout << "need interpolation for bin" << bin << "(<x>=" << hist->GetBinCenter(bin) << ")!"<< std::endl;
	// search for the two bins involved
	double binLow=0;
	double binHigh=0;
	for(int binRebinned=1; binRebinned<=binnedPlot->GetNbinsX(); ++binRebinned){
	  // search for bin in binned histo where upper border of is close to lower border of unbinned histo
	  if(std::abs(binnedPlot->GetBinLowEdge(binRebinned+1)-xlow)<binwidth){
	    binLow =binRebinned; 
	    binHigh=binRebinned+1;
	    break;
	  }
	}
	std::cout << "theory bin " << xlow << ".." << xhigh << "-> reco bins ";
	std::cout << binnedPlot->GetBinLowEdge(binLow) << ".." << binnedPlot->GetBinLowEdge(binLow+1) << " & ";
	std::cout << binnedPlot->GetBinLowEdge(binHigh) << ".." << binnedPlot->GetBinLowEdge(binHigh+1) << std::endl;
	// get the two points (this bin and the previous one)
	double x2=hist->GetBinCenter (bin  );
	double x1=hist->GetBinCenter (bin-1);
	double x3=hist->GetBinCenter (bin+1);
	double y2=hist->GetBinContent(bin  );
	double y1=hist->GetBinContent(bin-1);
	double y3=hist->GetBinContent(bin+1);	
	// calculate linear funtion
	double a=(y2-y1)/(x2-x1);
	double b=y1-a*x1;
	TF1* linInterpol=new TF1("linInterpol"+getTStringFromInt(bin), "[0]*x+[1]", x1, x2);
	linInterpol->SetParameter(0,a);
	linInterpol->SetParameter(1,b);
	// calculate the corresponding area of linear function to binned curve
	double contributionLowerBin=linInterpol->Integral(xlow,binnedPlot->GetBinLowEdge(binHigh));
	double contributionUpperBin=linInterpol->Integral(binnedPlot->GetBinLowEdge(binHigh),xhigh);
	// draw interpolation function for checking
	linInterpol->SetRange(xlow,xhigh);
	linInterpol->SetLineColor(kMagenta);
	linInterpol->DrawClone("same");
	// eventually use quadratic interpolation for the first harens m(ttbar) bin
	if(theory=="ahrens"&&variable.Contains("ttbarMass")&&usequad&&x2<450){
	  // calculate quadratic funtion
	  double a2=(y2-((y3-y1))*(x2-x1)/(x3-x1))/(x2*x2-x1*x1-(x2-x1)*(x3*x3+x1*x1)/(x3-x1));
	  double b2=((y3-y1)-a2*(x3*x3-x1*x1))/(x3-x1);
	  double c2=y1-a2*x1*x1-b*x1;
	  TF1* quadInterpol=new TF1("quadInterpol"+getTStringFromInt(bin), "[0]*x*x+[1]*x+[2]", x1, x2);
	  quadInterpol->SetParameter(0,a2);
	  quadInterpol->SetParameter(1,b2);
	  quadInterpol->SetParameter(2,c2);
	  // draw quad interpolation function for checking
	  quadInterpol->SetRange(xlow,xhigh);
	  quadInterpol->SetLineColor(kGreen);
	  quadInterpol->SetLineStyle(2);
	  hist->Fit(quadInterpol, "", "same", x1, x3);
	  // calculate the corresponding area of linear function to binned curve
	  quadInterpol->DrawClone("same");
	  double areaLow =quadInterpol->Integral(xlow,binnedPlot->GetBinLowEdge(binHigh) );
	  double areaHigh=quadInterpol->Integral(binnedPlot->GetBinLowEdge(binHigh),xhigh);
	  std::cout << "ratio(high/low) linear/quadratic: " << contributionUpperBin/contributionLowerBin << "/"<< areaHigh/areaLow << std::endl;
	  contributionLowerBin=y2*binwidth*areaLow /(areaLow+areaHigh);
	  contributionUpperBin=y2*binwidth*areaHigh/(areaLow+areaHigh);
	}
	// add fitted result
	binnedPlot->SetBinContent(binLow , binnedPlot->GetBinContent(binLow )+contributionLowerBin);
	binnedPlot->SetBinContent(binHigh, binnedPlot->GetBinContent(binHigh)+contributionUpperBin);
      }
    }
    // ensure over/underflow is 0
    binnedPlot->SetBinContent(0, 0.);
    binnedPlot->SetBinContent(binnedPlot->GetNbinsX()+1, 0.);
    // ensure normalization
    binnedPlot->Scale(1./(binnedPlot->Integral(0.,binnedPlot->GetNbinsX()+1)));
    // divide by binwidth
    if(divideByBinwidth){
      for(int bin=1; bin<=binnedPlot->GetNbinsX(); ++bin){
	binnedPlot->SetBinContent(bin, binnedPlot->GetBinContent(bin)/binnedPlot->GetBinWidth(bin));
      }
    }
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "result: binned output histo" << std::endl;
    for(int bin=1; bin<=binnedPlot->GetNbinsX(); ++bin){
      std::cout << "content bin " << bin << " (" << binnedPlot->GetBinLowEdge(bin) << ".." << binnedPlot->GetBinLowEdge(bin+1) << ")= " << binnedPlot->GetBinContent(bin) << std::endl;
    }

    // styling
    binnedPlot->SetLineColor(kBlue);
    binnedPlot->SetLineWidth(2);
    std::cout << "draw rebinned theory curve " << std::endl;
    binnedPlot->Draw("hist same");
    // draw bin boundaries
    std::cout << "draw bin boundaries " << std::endl;
    int binColor=kRed;
    int binWidth=2;
    int binStyle=6;
    for(int bin=0; bin<(int)bins_.size(); ++bin){
      if(!variable.Contains("ttbarMass")||bins_[variable][bin]>=345.) drawLine(bins_[variable][bin], 0, bins_[variable][bin], hist->GetMaximum(), binColor, binWidth, binStyle);
    }
    TH1F* line=(TH1F*)hist->Clone("line");
    line->SetLineColor(binColor);
    line->SetLineWidth(binWidth);
    line->SetLineStyle(binStyle);
    // legend
    TLegend *leg = new TLegend(0.7, 0.6, 0.95, 0.9);
    legendStyle(*leg,theory);
    if(points) leg ->AddEntry(ori, "original data points","P");
    leg ->AddEntry(hist      , "theory prediction" ,"P");
    leg ->AddEntry(line      , "reco binning"      ,"L");
    leg ->AddEntry(hist      , "linear interpolation" ,"L");
    leg ->AddEntry(binnedPlot, "rebinned curve"    ,"L");
    leg->Draw("same");
    std::cout << "done" << std::endl;
    // save in png and rootfile
    std::cout << std::endl << "do saving..." << std::endl;
    canv->SaveAs(variable+"Norm_Theory.png");
    TH1F* out=(TH1F*)binnedPlot->Clone();
    out->SetTitle(variable);
    out->SetName (variable);
    out->GetXaxis()->SetTitle(xSecLabelName(variable));
    TString yTile="#frac{1}{#sigma} #frac{d#sigma}{d";
    if(variable=="ttbarMass") yTile+="m^{t#bar{t}}} [GeV^{-1}]";
    if(variable=="ttbarPt") yTile+="p_{T}^{t#bar{t}}} [GeV^{-1}]";
    if(variable=="topPt") yTile+="p_{T}^{t}} [GeV^{-1}]";
    if(variable=="topY" ) yTile+="y^{t}}";
    out->GetYaxis()->SetTitle(yTile);
    out->SetLineColor(kOrange+2);
    out->SetLineStyle(2);
    std::cout << std::endl << "draw final result " << std::endl;
    TCanvas *canv2 = new TCanvas(variable+"Rebinned",variable+"Rebinned",800,600);
    canv2->cd();
    out->Draw();
    saveToRootFile(outputRootFile, out    , true, 0,""      );
    saveToRootFile(outputRootFile, rawHist, true, 0,"graph" );
    saveToRootFile(outputRootFile, canv   , true, 0,"detail");
    std::cout << "done!" << std::endl;
  }
}
