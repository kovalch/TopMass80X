#ifndef BCC_CPP_
#define BCC_CPP_

// =============================================================================
//
//  Class for calculating Bin Center Corrections
// 
//   a) - Calculating correction in x as center of gravity for each bin
//      - Calculating statistical error in x
//   b) - Calculating correction in x as intersection between fine binned
//        histogram and flat bin, code copied and adapted from [2]
//        
//        References:
//          [1] G. Lafferty and T. Wyatt, 
// 		"Where to Stick Data Points: 
//               The Treatment of Measurements within Wide Bins"
//              Nucl. Instr. and Methods, A355 (1995) 541
//          [2] https://hypernews.cern.ch/HyperNews/CMS/get/ewk/610.html
//  
//      - error in x is symmetric and half a bin width
//   
//   c) - In case of errors default values are geometrical bin center
//   
//  Input: 
//   
//   - ROOT-file with data tree
//   - std::vector<TString> with data branches to access
//   - Information if top and antitop branches shall be merged
// 
//  Output (accessible by public get-functions):
//  
//   - Corrected bin centers:
//     std::map<TString, std::vector<double> >
//   - Error on corrected bin centers:
//     std::map<TString, std::vector<double> >
//   - error code (definiton below)
//
//  Error Codes:
// 
//    0: No error found
//   -1: Input file not found
//   -2: Error in configuration file: 
//       number of binning-vectors and variables not corresponding
//   -3: Error while data filling into histograms 
//       number of binning-vectors and histograms not corresponding
//
// =============================================================================

#include "BCC.h"
#include "basicFunctions.h"

// ==============
//  Constructors
// ==============

BCC::BCC(){}

BCC::BCC(TString filename, std::vector<TString> vecPathnames, std::vector<TString> vecBranches, bool mergeLepHad, bool addCrossCheckVariables)
{
    error_code = 0;
    
    // Please provide your own function 'makeVariableBinning' producing a std::map with your binning
    
    mapBinning  = makeVariableBinning(addCrossCheckVariables);
    
    // Set variable 'largeMGFile', open file and get tree

    largeMGFile = (filename.Contains("Large")) ? 1 : 0;

    datafile = new TFile(filename,"READ");
    
    // Local variables

    float tempData = 0;
    std::vector<double> vecBranchData;
    BCCMap mapTempData;

    if (datafile)
    { 
        // fill data from tree into temporary map
       
        for (std::vector<TString>::iterator iterBranch = vecBranches.begin(); iterBranch != vecBranches.end(); iterBranch++){
	  
	    datafile -> cd(vecPathnames[iterBranch-vecBranches.begin()]);
	    TTree *treeData = (TTree*)gDirectory->Get("tree");
	    vecBranchData.clear();

	    std::cout << " Processing Data for .... " << filename << ": ";
	    std::cout << vecPathnames.at(iterBranch-vecBranches.begin()) << "/" << *iterBranch << std::endl;
	    TString activeBranchName = "";
	    activeBranchName = (*iterBranch);	 
	    // for b-jet quantities: treat thm like had/lep altough their name does not contain these appendices
	    activeBranchName.ReplaceAll("GenLep","Gen");
	    activeBranchName.ReplaceAll("GenHad","Gen");

	    if( activeBranchName.Contains("lepEtaMinus") || activeBranchName.Contains("lepEtaPlus")) activeBranchName="lepEta";

	    std::cout << " --> " << activeBranchName << std::endl;

	    treeData->SetBranchStatus("*",0);
	    treeData->SetBranchStatus(activeBranchName,1);
	    treeData->SetBranchAddress(activeBranchName,(&tempData));
	    for (int j=0; j<treeData->GetEntries(); j++)
	    {	
		treeData->GetEntry(j);
		vecBranchData.push_back(tempData);		
	    }
	    TString tempmapname=*iterBranch;
	    tempmapname.ReplaceAll("bbarq","bq");
	    mapTempData[tempmapname] = vecBranchData; 	    
	}
	vecBranchData.clear();

	// Merge data for top- and antitop-quarks if no distinction between them is required

	if (mergeLepHad)
	{
	    TString strNew, strHad, strLep;
	    for (BCCMap::iterator iter = mapTempData.begin(); iter != mapTempData.end(); iter++)
	    {
		vecBranchData.clear();

		// Default values

		strHad = "DefaultHad";
		strLep = "DefaultLep";
		strNew = iter->first;
		vecBranchData = iter->second; 
		//std::cout << "variable: " << strNew << std::endl;
		// Check if data belong to top -> bqq'
		
		if (iter != mapTempData.end() && iter->first.EndsWith("Had"))
		{ 
		    strHad = (iter->first); 
		    strHad.ReplaceAll("Gen"  ,"");
		    //std::cout << "had identified: " << strHad << std::endl;
		    iter++;
		}
		
		// Check if immediately following entry in map belongs to top -> blnu 
		// This works because the map is sorted alphabetically -> must be reviewed and fixed if new a branch is added between topxxHad and topxxLep
		//std::cout << iter->first << std::endl;
		//std::cout << "if (iter != mapTempData.end())"     << (iter != mapTempData.end()) << std::endl;
		//std::cout << "if (iter->first.EndsWith(\"Lep\")"  << iter->first.EndsWith("Lep") << std::endl;
		if (iter != mapTempData.end() && iter->first.EndsWith("Lep")){
		  strLep = (iter->first); 
		  strLep.ReplaceAll("Gen"  ,"");
		  //std::cout << "lep identified: " << strLep << std::endl;

		}
		if (strHad != "DefaultHad" && strLep != "DefaultLep" && strHad(0,strHad.Length()-3) == strLep(0,strLep.Length()-3) )
		{
		    strNew = strHad(0,strHad.Length()-3);
		    vecBranchData.insert(vecBranchData.end(),iter->second.begin(),iter->second.end());	
	    	    //std::cout << "merging vectors " << std::endl;
		}

		mapData[strNew] = vecBranchData;
	    }
	}
	if (mapBinning.size() != mapData.size() )
	{
	    error_code = -2;
	    std::cout << " WARNING ---- Size of maps for binning and data do not correspond!"     << std::endl;
	    std::cout << " WARNING ---- Calculation of corrected bin centres will be suppressed." << std::endl;
	    std::cout << " WARNING ---- Switch to geometrical centres as default values."         << std::endl;
	    std::cout << "        data entries: "     << std::endl;
	    for (BCCMap::iterator iter = mapData.begin(); iter != mapData.end(); iter++){
	      std::cout << "        " << iter->first << std::endl;
	    }
	    std::cout << "        binning entries: "     << std::endl;
	    for (BCCMap::iterator iter = mapBinning.begin(); iter != mapBinning.end(); iter++){
	      std::cout << "        " << iter->first << std::endl;
	    }
	}
	else
	  {
	    for (BCCMap::iterator iter = mapData.begin(); iter != mapData.end(); iter++) 

	      std::cout << " Creating new branch .... " << iter->first << " with " << iter->second.size() << " entries." << std::endl;

	    MakeHistos();	    
	  }	
    }
    else 
    {
	std::cout << "Input Data File not found. " << std::endl;
	error_code = -1;
    }
}

// ============
//  Destructor
// ============

BCC::~BCC(){}

// ============
//  Error Code
// ============

int BCC::getErrorCode()
{
  return error_code;
}

// =================================================
//  Set Default Values to Maps (error handling)
// =================================================

void BCC::setDefValues()
{
    std::cout << "Accessing Default Function." << std::endl;
    
    std::vector<double> vecTempBCCs;
    std::vector<double> vecTempBCCErrors;

    for (BCCMap::iterator iter = mapBinning.begin(); iter !=  mapBinning.end(); iter++ )
    {  
	std::vector<double>& refVecTempBinning = mapBinning[iter->first];  // required for ROOT purposes

	for (unsigned int bin = 1; bin < refVecTempBinning.size(); bin++)
	{
	    vecTempBCCs.push_back(refVecTempBinning[bin-1]+0.5*(refVecTempBinning[bin]-refVecTempBinning[bin-1])); 
	    vecTempBCCErrors.push_back(0.5*(refVecTempBinning[bin]-refVecTempBinning[bin-1]));
	}	
	TString plotname=iter->first;
	plotname.ReplaceAll("Gen","");
	mapCorrBinCentersInX[plotname]      = vecTempBCCs;
	mapCorrBinCenterErrorsInX[plotname] = vecTempBCCErrors;
	
	vecTempBCCs.clear();
	vecTempBCCErrors.clear();       
    }
}
// =================================================
//  Set Theory Histos
// =================================================

void BCC::MakeHistos()
{
  TString hname = "";
  double binmin = 0.0;
  double binmax = 1.0;
	
  // Can be replaced by a non-static way of finding an optimal binning
	
  std::map<TString,double> BinWidthScaler;
  std::map<TString,int>    HistoBins; 
  std::map<TString,int>    HistoSmoothingValue;

  BinWidthScaler["lepEta"]      = 80.0;
  BinWidthScaler["lepEtaPlus"]  = 80.0;
  BinWidthScaler["lepEtaMinus"] = 80.0;
  BinWidthScaler["lepPt"]       =    2;
  BinWidthScaler["topPt"]       =  0.4;
  BinWidthScaler["topPtLead"]   =  0.4;
  BinWidthScaler["topPtSubLead"]=  0.4;
  BinWidthScaler["topPtTtbarSys"]= 0.4;
  BinWidthScaler["topPtPlus"]   =  0.4;
  BinWidthScaler["topPtminus"]  =  0.4;
  BinWidthScaler["topY"]        = 40.0;
  BinWidthScaler["ttbarDelPhi"] = 40.0;
  BinWidthScaler["topYPlus"]    = 40.0;
  BinWidthScaler["topYMinus"]   = 40.0;
  BinWidthScaler["ttbarMass"]   =  1.0;
  BinWidthScaler["ttbarPt"]     =    1;
  BinWidthScaler["ttbarY"]      = 80.0;
  BinWidthScaler["bqPt"]        =    2;
  BinWidthScaler["bqEta"]       = 40.0;

  HistoBins["lepEta"]      = 200;
  HistoBins["lepEtaPlus"]  = 200;
  HistoBins["lepEtaMinus"] = 200;
  HistoBins["lepPt"]       = 160;
  HistoBins["topPt"]       = 160;
  HistoBins["topPtLead"]   = 160;
  HistoBins["topPtSubLead"]= 160;
  HistoBins["topPtTtbarSys"]= 160;
  HistoBins["topPtPlus"]   = 160;
  HistoBins["topPtMinus"]  = 160;
  HistoBins["topY"]        = 200;
  HistoBins["topYPlus"]    = 200;
  HistoBins["ttbarDelPhi"] = 200;
  HistoBins["topYMinus"]   = 200;
  HistoBins["ttbarMass"]   = 240;
  HistoBins["ttbarPt"]     = 300;
  HistoBins["ttbarY"]      = 200;
  HistoBins["bqPt"]        = 160;
  HistoBins["bqEta"]       = 200;

  HistoSmoothingValue["lepEta"]      = (largeMGFile) ? 2 :  4;
  HistoSmoothingValue["lepEtaPlus"]  = (largeMGFile) ? 2 :  4;
  HistoSmoothingValue["lepEtaMinus"] = (largeMGFile) ? 2 :  4;
  HistoSmoothingValue["lepPt"]       = (largeMGFile) ? 0 :  0;
  HistoSmoothingValue["topPtPlus"]   = (largeMGFile) ? 1 : 10;
  HistoSmoothingValue["topPtMinus"]  = (largeMGFile) ? 1 : 10; 
  HistoSmoothingValue["topPt"]       = (largeMGFile) ? 1 : 10;
  HistoSmoothingValue["topY"]        = (largeMGFile) ? 1 :  5; 
  HistoSmoothingValue["topYPlus"]    = (largeMGFile) ? 1 :  5;
  HistoSmoothingValue["topYMinus"]   = (largeMGFile) ? 1 :  5;
  HistoSmoothingValue["ttbarMass"]   = (largeMGFile) ? 1 : 10;
  HistoSmoothingValue["ttbarPt"]     = (largeMGFile) ? 1 :  5;
  HistoSmoothingValue["ttbarY"]      = (largeMGFile) ? 2 :  2;
  HistoSmoothingValue["bqPt"]        = (largeMGFile) ? 0 :  0;
  HistoSmoothingValue["bqEta"]       = (largeMGFile) ? 2 :  2;

  for (BCCMap::iterator iter = mapData.begin(); iter != mapData.end(); iter++){

    //TString& refname= iter->first;
    TString plotname=iter->first;
    //TString plotname=(TString)(refname);
    plotname.ReplaceAll("Gen","");
    const std::vector<double>& refVecTempBinning = mapBinning[plotname];  // required for ROOT purposes
    const std::vector<double>& refVecTempData    = iter->second;	       // required for ROOT purposes
    
    std::vector<double>::const_iterator iterVecData;
    
    hname  = iter->first;
    hname.ReplaceAll("Gen","");
    plotname.ReplaceAll("Gen","");
    std::cout << "refVecTempBinning.size(): " << refVecTempBinning.size() << std::endl;
    std::cout << "mapData entry valid? " << iter->first << std::endl;
    binmin = refVecTempBinning[0];
    binmax = refVecTempBinning[refVecTempBinning.size()-1];

    TH1F temphisto(hname,hname, BinWidthScaler[hname]*(binmax-binmin),binmin,binmax);
    //TH1F temphisto(hname,hname, HistoBins[hname], binmin, binmax);

    
    for (iterVecData = refVecTempData.begin(); iterVecData != refVecTempData.end(); iterVecData++ )
      temphisto.Fill(*iterVecData);

    //temphisto.Smooth(HistoSmoothingValue[hname]);

    mapHistos[hname] = temphisto;
  }
  
  if (mapBinning.size() != mapHistos.size() ){

    error_code = -3;
    std::cout << " WARNING ---- Size of maps for binning and histograms do not correspond!" << std::endl;
    std::cout << " WARNING ---- Calculation of corrected bin centres will be suppressed."   << std::endl;
    std::cout << " WARNING ---- Switch to geometrical centres as default values."           << std::endl;
  }
  else{

    for (std::map<TString, TH1F>::iterator iterHistos = mapHistos.begin(); iterHistos != mapHistos.end(); iterHistos++)
          
      std::cout << " Creating new histo .... " << (iterHistos->first) << " with " << (iterHistos->second).GetEntries() << " entries." << std::endl;
  }	
}

// =======================================================
//  Set Beam Centre Corrections in x where coarsely 
//  binned histogram intersects smooth curve
// =======================================================

void BCC::setBCCinX_IntersectionInBin()
{
  for (std::map<TString, TH1F>::iterator iterHistos = mapHistos.begin(); iterHistos != mapHistos.end(); iterHistos++)
  {
    TH1F& refHisto                     = iterHistos->second;	         // required for ROOT purposes
    std::vector<double>& refVecBinning = mapBinning[iterHistos->first];  // required for ROOT purposes
    
    //Note that the bins of the theory histogram should be sub-bins of the corr. data histogram
    
    const int NBinsTheory = refHisto.GetNbinsX();
    const int NBinsData   = refVecBinning.size()-1;  // -1 because refVecBinning contains bin boundaries
    
    std::cout << " Processing BCC for " << iterHistos->first << " .... NBinsTheory: " << NBinsTheory << " and NBinsData: " << NBinsData << std::endl;
    
    double weight = 0.0;
    double sumWeights[NBinsData];
    double integral[NBinsData];
    double average[NBinsData];
    double avXvalues[NBinsData]; 
    
    std::vector<double> vecTempBCCs;
    std::vector<double> vecTempBCCErrors;
    
    int NTheoryBinsPerDataBin[NBinsData];
    
    // Reset variables:
    
    for(int jBinData=0; jBinData < NBinsData; jBinData++)
    {
      sumWeights[jBinData] = 0.0;
      integral[jBinData]   = 0.0;
      average[jBinData]    = 0.0;
      NTheoryBinsPerDataBin[jBinData] = 0;
      avXvalues[jBinData]  = 0.0;
    }
    
    // Loop over bins of theory histogram and check into which data this falls
    // Calculate expressions needed for the weighted average
    
    for (int jBinTheory = 1; jBinTheory <= NBinsTheory; jBinTheory++)
    {
      double xTheory = refHisto.GetBinCenter(jBinTheory);
      
      for(int jBinData = 0; jBinData < NBinsData; jBinData++)
      {
	if (xTheory > refVecBinning[jBinData] && xTheory < refVecBinning[jBinData+1])
	{
	  weight                = refHisto.GetBinContent(jBinTheory);  // y-value of theory bin
	  integral[jBinData]   += weight * xTheory;                    // integral over theory histogram 
	  sumWeights[jBinData] += weight;			       // will later be used for mean value over data bin	
	  NTheoryBinsPerDataBin[jBinData]++;
	}	
      }			
    }
    
    // Determine x-values at which mean y-entry in data bin intersects y-value of theory histogram
    
    double min  = 0.0;
    
    for (int jBinData = 0; jBinData < NBinsData; jBinData++)
    {
      if(sumWeights[jBinData] > 0.0)
      {
	min = (double)LONG_MAX;
	avXvalues[jBinData] =0.0;
	average[jBinData] = sumWeights[jBinData]/NTheoryBinsPerDataBin[jBinData];
	for(int jBinTheory = 1; jBinTheory <= NBinsTheory; jBinTheory++)
	{
	  double xTheory     = refHisto.GetBinCenter(jBinTheory);
	  double diff        = std::abs(average[jBinData] - refHisto.GetBinContent(jBinTheory));
	  
	  if ( (xTheory > refVecBinning[jBinData] && xTheory < refVecBinning[jBinData+1]) && diff <= min) 
	  {
	    // use first intersection point in ttbarPt bin1
	    if(!((iterHistos->first).Contains("ttbarPt")&&jBinData==0)||(refHisto.GetBinCenter(jBinTheory)<15)){
	      // use no intersection points near to bin borders
	      //if(std::abs((refHisto.GetBinCenter(jBinTheory)-refVecBinning[jBinData])/refVecBinning[jBinData])>0.05&&std::abs((refHisto.GetBinCenter(jBinTheory)-refVecBinning[jBinData+1])/refVecBinning[jBinData+1])<0.05){
	      //std::cout << iterHistos->first << ": theory bin " << jBinTheory << "y average: " << average[jBinData] << ", y-diff(bin) " << diff << "<" << min << ", xval "<< refHisto.GetBinCenter(jBinTheory) << std::endl;
	      min = diff;
	      avXvalues[jBinData] = refHisto.GetBinCenter(jBinTheory);	      
	      if(((iterHistos->first).Contains("ttbarMass"))&&avXvalues[jBinData]>1270&&avXvalues[jBinData]<1310) avXvalues[jBinData] =1290;
	    }
	  }
	}
      }
      
      std::cout << Form(" Bin %1.2f .... %1.2f : mean data = %1.4f and mean x = %1.2f",
      			  refVecBinning[jBinData], refVecBinning[jBinData+1], average[jBinData], avXvalues[jBinData]) << std::endl;			
      
      vecTempBCCs.push_back(avXvalues[jBinData]);
      vecTempBCCErrors.push_back(0.5*refHisto.GetXaxis()->GetBinWidth(1));
    }
    
    mapCorrBinCentersInX[iterHistos->first]      = vecTempBCCs;
    mapCorrBinCenterErrorsInX[iterHistos->first] = vecTempBCCErrors;
    
    vecTempBCCs.clear();
    vecTempBCCErrors.clear();
  }	
}

// =======================================================
//  Set Beam Centre Corrections in x as mean over bin
// =======================================================

void BCC::setBCCinX_MeanOverBin()
{		
  double sum, error;
  int NEVE;

  for (BCCMap::iterator iter = mapData.begin(); iter != mapData.end(); iter++)
  {	
    std::vector<double>& refVecTempBinning = mapBinning[iter->first];  // required for ROOT purposes
    std::vector<double>& refVecTempData    = iter->second;	           // required for ROOT purposes
    std::vector<double>  vecTempBCCs;
    std::vector<double>  vecTempBCCErrors;
    
    std::cout << iter->first                << ": ";
    std::cout << refVecTempBinning.size()-1 << " bins, ";
    std::cout << refVecTempData.size()      << " events." << std::endl;
    
    // Sort vector to ease comparision to bin boundaries (think of 'binning' the vector)
    
    std::sort(refVecTempData.begin(),refVecTempData.end());
    
    // ======================================================
    //  Get Corrected Bin Centre - center of gravity per bin
    // ======================================================
    
    // Get Mean Value
    
    unsigned int indexData = 0;
    
    for (indexData = 0; indexData < refVecTempData.size();)
    {
      for (unsigned int bin = 1; bin < refVecTempBinning.size(); bin++)
      {
	std::cout << "Upper Bin Boundary: " << refVecTempBinning[bin]; 
	
	sum  = 0;
	NEVE = 0;
	
	while (refVecTempData[indexData] < refVecTempBinning[bin] && indexData < refVecTempData.size() )
	{
	  sum += refVecTempData[indexData];
	  indexData++;
	  NEVE++;
	}
	
	// std::cout << " " <<(sum/(double)NEVE) << " " << NEVE << std::endl;
	
	(NEVE > 0) ? vecTempBCCs.push_back(sum/(double)NEVE) : vecTempBCCs.push_back(refVecTempBinning[bin-1]+0.5*(refVecTempBinning[bin]-refVecTempBinning[bin-1]));  // fill vector with corrected bin centers 
	
	//  compensation of endless loop for data which would be assigned to the overflow bin
	
	while (refVecTempData[indexData] > refVecTempBinning[refVecTempBinning.size()-1] && indexData < refVecTempData.size() ) indexData++; 
      }	    
    }
    
    // Get Statistical Error
    
    indexData = 0;
	
    for (unsigned int bin = 1; bin < refVecTempBinning.size(); bin++)
    {
      error = 0;
      NEVE  = 0;
      
      while (refVecTempData[indexData] < refVecTempBinning[bin] && indexData < refVecTempData.size() )
      {
	error += pow(refVecTempData[indexData]-vecTempBCCs[bin-1],2);
	indexData++;
	NEVE++;
      }
      
      // std::cout << " " << sqrt(error)/(double)NEVE << " " << NEVE << std::endl; 
      
      (NEVE != 0) ? vecTempBCCErrors.push_back(sqrt(error)/(double)NEVE) : vecTempBCCErrors.push_back(0);  // fill vector with statistical errors on corrected bin centers
    }
    
    // fill maps
    
    TString plotname=iter->first;
    plotname.ReplaceAll("Gen","");
    mapCorrBinCentersInX[plotname]      = vecTempBCCs;
    mapCorrBinCenterErrorsInX[plotname] = vecTempBCCErrors;
    
    vecTempBCCs.clear();
    vecTempBCCErrors.clear();       
  }
}

// ================================================
//  Functions to pass results to external calls
// ================================================

// Corrected bin centers

BCCMap BCC::getMapWithCorrectedCentersInX()
{
    return mapCorrBinCentersInX;  
}

// Error on corrected bin centers

BCCMap BCC::getMapWithCenterErrorsInX()
{
    return mapCorrBinCenterErrorsInX;  
}

// ================================================
//  Main test routine
// ================================================

void BCC::runBCCCalculation()
{
  if (error_code == 0)
  { 
    setBCCinX_IntersectionInBin();
    //setBCCinX_MeanOverBin(); 
  }
  else setDefValues();    
}

#endif
