#ifndef _NNLO_H_
#define _NNLO_H_

#include <iostream>
#include <vector>
#include <cmath>
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveText.h"

class NNLO {

	private: 
	
		TString parPlotName;
		unsigned int parVerbose;
	
		std::vector<double> vecBins;
		double* binArray;	
	
		TGraph *graphUnnormalizedNNLO, *graphUnnormalizedNNLO_ScaleDn, *graphUnnormalizedNNLO_ScaleUp;
		TGraph *graphNormalizedNNLO,   *graphNormalizedNNLO_ScaleDn,   *graphNormalizedNNLO_ScaleUp;

		TH1F *histoRebinnedNormalizedNNLO, *histoRebinnedNormalizedNNLO_ScaleDn, *histoRebinnedNormalizedNNLO_ScaleUp;
		TGraphAsymmErrors *graphRebinnedNormalizedNNLO_ErrorBand;
	
		TLegend *legUnnormalizedNNLO, *legNormalizedNNLO, *legNormalizedReBinnedNNLO;

		TCanvas *c1_graph, *c2_graph, *c3_histo;
	
		// Basic Functions
	
		void drawControlPlots();
		void printControlOutput();
		void setObjectStyles(TString xAxisTitle="", TString xAxisUnit="", 
							 TString yAxisTitle="", TString yAxisUnit="", TString yAxisTitleNorm="", TString yAxisUnitNorm="", 
							 bool logXScale=0, double ymax=1.0);
		void fillLegends();
	
		void writeOutputFile(TString file_open_option);
	
		// Auxiliary functions
	
		inline double getLocalIntegral(double xLeft, double xRight, double yLeft, double yRight);
		inline double getNormalization(const TGraph* graph);
		inline double getNormalization(const TH1* histo);
		inline void   rescaleGraph(TGraph* graph, const double scaleFactor);
		void   createRebinnedHisto(const TGraph* graph, TH1* histo, const std::vector<double>& vecBinning);
		void   setBinErrors(TGraphAsymmErrors* graph, const TH1* histoCentral, const TH1* histoDn, const TH1* histoUp);
	
		// Style functions
	
		void makeGraphStyle(TGraph* graph, int color, int markerStyle, TString name, TString title,TString xAxisTitle, TString xAxisUnit, TString yAxisTitle, TString yAxisUnit);
		void makeHistoStyle(TH1* histo, int color, int markerStyle, TString name, TString title,TString xAxisTitle, TString xAxisUnit, TString yAxisTitle, TString yAxisUnit);
		void makeLegendStyle(TLegend *leg);
		void makeCanvasStyle(TCanvas *canvas, bool logXScale=0);

	public:
	
		NNLO();
		NNLO(TString plotname, bool verbose=0);
		~NNLO();
		
		int runNNLO(TString file_open_option="");	
	
};

#endif
