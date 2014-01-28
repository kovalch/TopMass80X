/*
HOW TO:

mkdir -p ControlPlots
mkdir -p InputFiles
mkdir -p ROOTFiles
cp data/pttopnnloapprox8lhc173m.txt InputFiles/topPt.txt
cp data/ytopnnloapprox8lhc173m.txt InputFiles/topY.txt
cp data/Mtt_8000_172.5_Mtt_fin.dat InputFiles/ttbarM.txt
cp data/Ptt_8000_172.5_Ptt_fin.dat InputFiles/ttbarPt.txt
cd InputFiles
ln -s topPt.txt topPt_ScaleUp.txt
ln -s topPt.txt topPt_ScaleDown.txt
ln -s topY.txt topY_ScaleUp.txt
ln -s topY.txt topY_ScaleDown.txt
ln -s ttbarM.txt ttbarM_ScaleUp.txt
ln -s ttbarM.txt ttbarM_ScaleDown.txt
ln -s ttbarPt.txt ttbarPt_ScaleUp.txt
ln -s ttbarPt.txt ttbarPt_ScaleDown.txt
cd -

if you want to create kidonakis curves execute:

    root -l
    .L macros/NNLO.cpp++
    NNLO pt=NNLO("kidonakis", "topPt",2)
    pt.runNNLO("UPDATE")
    NNLO y=NNLO("kidonakis", "topY",2)
    y.runNNLO("UPDATE")

instead if you want to create Ahrens curves execute:

    root -l
    .L macros/NNLO.cpp++
    NNLO ptttbar=NNLO("ahrens", "ttbarPt",2)
    ptttbar.runNNLO("UPDATE")
    NNLO mttbar=NNLO("ahrens", "ttbarM",2)
    mttbar.runNNLO("UPDATE")



-> result is in:
        ROOTFiles/kidonakisNNLO_8TeV.root
        ROOTFiles/ahrensNNLL_8TeV.root


*/


#ifndef _NNLO_CPP_
#define _NNLO_CPP_

#include "NNLO.h"

// ==============
//  Constructor
// ==============

NNLO::NNLO(){}

NNLO::NNLO(TString theory, TString plotname, bool verbose, TString Energy){

    // General seetings

    gStyle -> SetOptStat(0);
    gStyle -> SetPadLeftMargin(0.15);

    Theory = theory;
    parPlotName = plotname;
    parVerbose  = verbose;

    if(Theory != "kidonakis" && Theory != "ahrens")
    {
        std::cout<<"Not valid theory prediction. I don't know to treat the binning"<<std::endl;
        exit(88);
    } else {
        std::cout<<"Making plots alla '"<<Theory<<"'\n"<<std::endl;
        outputFilename = Theory;
        if(Theory == "ahrens") outputFilename.Append("NNLL_");
        else if(Theory == "kidonakis") outputFilename.Append("NNLO_");
        outputFilename.Append(Energy+"TeV");
    }

    // l+jets
//     const double topPtBins[] = { 1.0,60.0,100.0,150.0,200.0,260.0,320.0,400.0,800.0,1500.};
//     const double topYBins[]  = {-5.0,-2.5,-1.6,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.6,2.5,5.0};
//    const double ttbarMBins[] = {345.,400.,470.,550.,650.,800.,1100.,1600.};
//     const double ttbarPtBins[] = {0.,20.,60.,120.,300.,1300.};

    // dilepton

    const double topPtBins[] = { 0.0,80.0,130.0,200.0,300.0,400.0,1500.};
    const double topYBins[]  = {-5.0,-2.5,-1.3,-0.8,-0.4,0.0,0.4,0.8,1.3,2.5,5.0};
    const double ttbarMBins[] = {340.,400.,470.,550.,650.,800.,1100.,1600.,2720.};
    const double ttbarPtBins[] = {0.,20.,60.,120.,300.,1300.};

    // Cannel specific variables

    TString xAxisTitle, yAxisTitle, yAxisTitleNorm;
    TString xAxisUnit,  yAxisUnit,  yAxisUnitNorm;

    bool logXScale = 0;
    double ymax=1.0;

    if (plotname=="topPt"){
        vecBins.insert(vecBins.begin(), topPtBins, topPtBins + sizeof(topPtBins)/sizeof(double));
        xAxisTitle="p_{T}";
        xAxisUnit="#left[GeV#right]";
        yAxisTitle="d#sigma/dp_{T}";
        yAxisUnit="#left[pb/GeV#right]";
        yAxisTitleNorm="1/#sigma d#sigma/dp_{T}";
        yAxisUnitNorm="#left[1/GeV#right]";
        logXScale=1;
        ymax=0.01;
    } else if (plotname=="topY"){
        vecBins.insert(vecBins.begin(), topYBins,  topYBins  + sizeof(topYBins)/sizeof(double));
        xAxisTitle="y";
        xAxisUnit="";
        yAxisTitle="d#sigma/dy";
        yAxisUnit="#left[pb#right]";
        yAxisTitleNorm="1/#sigma d#sigma/dy";
        yAxisUnitNorm="";
        logXScale=0;
        ymax=0.5;
    } else if (plotname=="ttbarM"){
        vecBins.insert(vecBins.begin(), ttbarMBins,  ttbarMBins  + sizeof(ttbarMBins)/sizeof(double));
        xAxisTitle="m ";
        xAxisUnit="#left[GeV#right]";
        yAxisTitle="d#sigma/dm";
        yAxisUnit="#left[pb#right]";
        yAxisTitleNorm="1/#sigma d#sigma/dm";
        yAxisUnitNorm="#left[1/GeV#right]";
        logXScale=0;
        ymax=0.006;
    } else if (plotname=="ttbarPt"){
        vecBins.insert(vecBins.begin(), ttbarPtBins,  ttbarPtBins  + sizeof(ttbarPtBins)/sizeof(double));
        xAxisTitle="p_{T} ";
        xAxisUnit="#left[GeV#right]";
        yAxisTitle="d#sigma/dp_{T}";
        yAxisUnit="#left[pb/GeV#right]";
        yAxisTitleNorm="1/#sigma d#sigma/dp_{T}";
        yAxisUnitNorm="#left[1/GeV#right]";
        logXScale=1;
        ymax=0.03;
    } else {
        std::cout<<"Variable '"<<plotname<<"' not valid"<<std::endl;
        std::cout<<"  Choose between topPt, topY, ttbarM, ttbarPt"<<std::endl;
        exit(89);
    }

    // Prepare array with new binning (for the sake of common naming - no essential difference to topPtBins[] or topYBins[])

    const int binArraySize = vecBins.size();
    binArray = new double[binArraySize];

    for(int j=0;j<binArraySize;j++){
        binArray[j]=vecBins[j];
        if (parVerbose>0) std::cout << "Binning array entry " << j << ": " << binArray[j] << std::endl;
    }

    // ============================
    //  Initialize all objects
    // ============================

    // Canvases

    c1_graph = new TCanvas("c1_graph","c1_graph",800,600);
    c2_graph = new TCanvas("c2_graph","c2_graph",800,600);
    c3_histo = new TCanvas("c3_histo","c3_histo",800,600);

    makeCanvasStyle(c1_graph,logXScale);
    makeCanvasStyle(c2_graph,logXScale);
    makeCanvasStyle(c3_histo,logXScale);

    // Legends 

    legUnnormalizedNNLO       = new TLegend();
    legNormalizedNNLO         = new TLegend();
    legNormalizedReBinnedNNLO = new TLegend();

    // unnormalized graphs

    graphUnnormalizedNNLO         = new TGraph("InputFiles/"+plotname+".txt");
    graphUnnormalizedNNLO_ScaleDn = new TGraph("InputFiles/"+plotname+"_ScaleDown.txt");
    graphUnnormalizedNNLO_ScaleUp = new TGraph("InputFiles/"+plotname+"_ScaleUp.txt");

    // normalized graphs

    graphNormalizedNNLO         = (TGraph*)graphUnnormalizedNNLO         -> Clone("graphNormalizedNNLO");
    graphNormalizedNNLO_ScaleDn = (TGraph*)graphUnnormalizedNNLO_ScaleDn -> Clone("graphNormalizedNNLO_ScaleDown");
    graphNormalizedNNLO_ScaleUp = (TGraph*)graphUnnormalizedNNLO_ScaleUp -> Clone("graphNormalizedNNLO_ScaleUp");

    // rebinned, normalized histograms

    histoRebinnedNormalizedNNLO         = new TH1F("histoRebinnedNormalizedNNLO","histoRebinnedNormalizedNNLO",binArraySize-1,binArray);
    histoRebinnedNormalizedNNLO_ScaleDn = new TH1F("histoRebinnedNormalizedNNLO_ScaleDn","histoRebinnedNormalizedNNLO_ScaleDn",binArraySize-1,binArray);
    histoRebinnedNormalizedNNLO_ScaleUp = new TH1F("histoRebinnedNormalizedNNLO_ScaleUp","histoRebinnedNormalizedNNLO_ScaleUp",binArraySize-1,binArray);

    // rebinned, normalized graph with asymmetric error

    //graphRebinnedNormalizedNNLO_ErrorBand = new TGraphAsymmErrors();

    // object styles

    setObjectStyles(xAxisTitle, xAxisUnit, yAxisTitle, yAxisUnit, yAxisTitleNorm, yAxisUnitNorm, logXScale,ymax);
}

// ==============
//  Destructor
// ==============

NNLO::~NNLO(){}

// ============================
//  Draw basic control plots
// ============================

void NNLO::drawControlPlots(){

    // Unnormalized graphs

    c1_graph->cd();
    graphUnnormalizedNNLO         -> Draw("AP");
    graphUnnormalizedNNLO_ScaleDn -> Draw("P");
    graphUnnormalizedNNLO_ScaleUp -> Draw("P");
    legUnnormalizedNNLO           -> Draw("SAME");
    c1_graph->Print("ControlPlots/"+parPlotName+"_UnnormalizedGraph.eps");
    c1_graph->Print("ControlPlots/"+parPlotName+"_UnnormalizedGraph.C");

    // Normalized graphs

    c2_graph->cd();
    graphNormalizedNNLO         -> Draw("AP");
    graphNormalizedNNLO_ScaleDn -> Draw("P");
    graphNormalizedNNLO_ScaleUp -> Draw("P");
    legNormalizedNNLO           -> Draw("SAME");
    c2_graph->Print("ControlPlots/"+parPlotName+"_NormalizedGraph.eps");
    c2_graph->Print("ControlPlots/"+parPlotName+"_NormalizedGraph.C");

    // Rebinned histograms

    c3_histo->cd();
    graphRebinnedNormalizedNNLO_ErrorBand->Draw("AF");
    //histoRebinnedNormalizedNNLO->Draw("HIST SAME");
    legNormalizedReBinnedNNLO->Draw("SAME");
    c3_histo->Print("ControlPlots/"+parPlotName+"_RebinnedFromGraph.eps");
    c3_histo->Print("ControlPlots/"+parPlotName+"_RebinnedFromGraph.C");
}

// ==========================================
//  Print control output (normalization)
// ==========================================

void NNLO::printControlOutput(){

    if(parVerbose == 0) return;
    std::cout << std::endl;
    std::cout << "UNNORMALIZED GRAPHS" << std::endl;

    std::cout << "Graph NNLO, unnormalized, no sys     - integral: " << getNormalization(graphUnnormalizedNNLO)         << std::endl;
    std::cout << "Graph NNLO, unnormalized, scale down - integral: " << getNormalization(graphUnnormalizedNNLO_ScaleDn) << std::endl;
    std::cout << "Graph NNLO, unnormalized, scale up   - integral: " << getNormalization(graphUnnormalizedNNLO_ScaleUp) << std::endl;

    // ============================

    std::cout << std::endl;
    std::cout << "NORMALIZED GRAPH" << std::endl;

    std::cout << "Graph NNLO, normalized, no sys     - integral: " << getNormalization(graphNormalizedNNLO)         << std::endl;
    std::cout << "Graph NNLO, normalized, scale down - integral: " << getNormalization(graphNormalizedNNLO_ScaleDn) << std::endl;
    std::cout << "Graph NNLO, normalized, scale up   - integral: " << getNormalization(graphNormalizedNNLO_ScaleUp) << std::endl;

    // ============================

    std::cout << std::endl;
    std::cout << "REBINNED NORMALIZED HISTOGRAM - from Graph" << std::endl;

    std::cout << "Histo NNLO, rebinned and normalized, no sys     - integral: " << getNormalization(histoRebinnedNormalizedNNLO)         << std::endl;
    std::cout << "Histo NNLO, rebinned and normalized, scale down - integral: " << getNormalization(histoRebinnedNormalizedNNLO_ScaleDn) << std::endl;
    std::cout << "Histo NNLO, rebinned and normalized, scale up   - integral: " << getNormalization(histoRebinnedNormalizedNNLO_ScaleUp) << std::endl;

}

// ============================
//  Write output to file
// ============================

void NNLO::fillLegends(){

    legUnnormalizedNNLO -> AddEntry(graphUnnormalizedNNLO_ScaleDn, "Scale Down",   "L");
    legUnnormalizedNNLO -> AddEntry(graphUnnormalizedNNLO,         "No Variation", "L");
    legUnnormalizedNNLO -> AddEntry(graphUnnormalizedNNLO_ScaleUp, "Scale Up",     "L");

    legNormalizedNNLO -> AddEntry(graphNormalizedNNLO_ScaleDn, "Scale Down",  "L");
    legNormalizedNNLO -> AddEntry(graphNormalizedNNLO,         "No Variation","L");
    legNormalizedNNLO -> AddEntry(graphNormalizedNNLO_ScaleUp, "Scale Up",    "L");

    legNormalizedReBinnedNNLO -> AddEntry(histoRebinnedNormalizedNNLO_ScaleDn,   "Scale Down"  , "L" );
    legNormalizedReBinnedNNLO -> AddEntry(histoRebinnedNormalizedNNLO,           "No Variation", "FL");
    legNormalizedReBinnedNNLO -> AddEntry(histoRebinnedNormalizedNNLO_ScaleUp,   "Scale Up"    , "L" );
    legNormalizedReBinnedNNLO -> AddEntry(graphRebinnedNormalizedNNLO_ErrorBand, "Error Band"  , "FL");
}

// ============================
//  Write output to file
// ============================

void NNLO::writeOutputFile(TString file_open_option){

    TFile *f = new TFile("ROOTFiles/"+outputFilename+".root",file_open_option);

    graphNormalizedNNLO->SetName(parPlotName+"_graph");
    histoRebinnedNormalizedNNLO->SetName(parPlotName);

    graphNormalizedNNLO->Write();
    histoRebinnedNormalizedNNLO->Write();

    f->Close();
}

// ============================
//  Format class objects
// ============================

void NNLO::setObjectStyles(TString xAxisTitle, TString xAxisUnit, 
                           TString yAxisTitle, TString yAxisUnit, TString yAxisTitleNorm, TString yAxisUnitNorm, 
                           bool logXScale, double ymax){

    makeCanvasStyle(c1_graph,logXScale);
    makeCanvasStyle(c2_graph,logXScale);
    makeCanvasStyle(c3_histo,logXScale);

    makeLegendStyle(legUnnormalizedNNLO);
    makeLegendStyle(legNormalizedNNLO);
    makeLegendStyle(legNormalizedReBinnedNNLO);

    makeGraphStyle(graphUnnormalizedNNLO,         38, 7, "Graph - Unnormalized, No Sys  ", "", xAxisTitle, xAxisUnit, yAxisTitle, yAxisUnit);
    makeGraphStyle(graphUnnormalizedNNLO_ScaleDn, 36, 7, "Graph - Unnormalized, Scale Dn", "", xAxisTitle, xAxisUnit, yAxisTitle, yAxisUnit);
    makeGraphStyle(graphUnnormalizedNNLO_ScaleUp, 40, 7, "Graph - Unnormalized, Scale Up", "", xAxisTitle, xAxisUnit, yAxisTitle, yAxisUnit);

    makeGraphStyle(graphNormalizedNNLO,         38, 7, "Graph - Normalized, No Sys  ", "", xAxisTitle, xAxisUnit, yAxisTitleNorm, yAxisUnitNorm);
    makeGraphStyle(graphNormalizedNNLO_ScaleDn, 36, 7, "Graph - Normalized, Scale Dn", "", xAxisTitle, xAxisUnit, yAxisTitleNorm, yAxisUnitNorm);
    makeGraphStyle(graphNormalizedNNLO_ScaleUp, 40, 7, "Graph - Normalized, Scale Up", "", xAxisTitle, xAxisUnit, yAxisTitleNorm, yAxisUnitNorm);

    makeHistoStyle(histoRebinnedNormalizedNNLO,         38, 5, "Histo - Rebinned from Graph, No Sys  ", "", xAxisTitle, xAxisUnit, yAxisTitleNorm, yAxisUnitNorm);
    makeHistoStyle(histoRebinnedNormalizedNNLO_ScaleDn, 36, 5, "Histo - Rebinned from Graph, Scale Dn", "", xAxisTitle, xAxisUnit, yAxisTitleNorm, yAxisUnitNorm);
    makeHistoStyle(histoRebinnedNormalizedNNLO_ScaleUp, 40, 5, "Histo - Rebinned from Graph, Scale Up", "", xAxisTitle, xAxisUnit, yAxisTitleNorm, yAxisUnitNorm);

    histoRebinnedNormalizedNNLO->SetMaximum(ymax);
    histoRebinnedNormalizedNNLO->SetMarkerSize(0);
    histoRebinnedNormalizedNNLO_ScaleDn->SetMaximum(ymax);
    histoRebinnedNormalizedNNLO_ScaleUp->SetMaximum(ymax);
}

// ===================================
//  Integral between two data points
// ===================================

inline double NNLO::getLocalIntegral(double xLeft, double xRight, double yLeft, double yRight){

    if(Theory == "kidonakis"){// Kidonakis give the value 'y' in a given point 'x', eg: f(x) = y
        return 0.5*(yLeft+yRight)*std::abs(xRight-xLeft);
    } else if(Theory == "ahrens"){// Ahrens, instead, gives the integrated value 'y' in a bin with bin center 'x'
        return yLeft*std::abs(xRight-xLeft);
    }
    return 0;
}

// ===============================================================
//  Normalization of TGraph-object (overloaded with TH1 version)
// ===============================================================

inline double NNLO::getNormalization(const TGraph* graph){

    double integral = 0.0;

    double xLeft, xRight;
    double yLeft, yRight;
    int j;

    int NPoints = 0;
    if(Theory == "kidonakis") NPoints = graph->GetN()-1;
    else if(Theory == "ahrens") NPoints = graph->GetN();

    for (j=0;j<NPoints;j++){
        graph->GetPoint(j,  xLeft, yLeft);
        graph->GetPoint(j+1,xRight,yRight);
        integral += getLocalIntegral(xLeft,xRight,yLeft,yRight);
    }
    return integral;
}

// ===============================================================
//  Normalization of TH1-object (overloaded with TGraph version)
// ===============================================================

inline double NNLO::getNormalization(const TH1* histo){

    double integral = 0.0;

    for(int j=0;j<1+histo->GetNbinsX();j++){
        integral += (histo->GetBinContent(j+1) * histo->GetBinWidth(j+1));
    }

    return integral;
}

// ========================================================
//  Rescaling of TGraph-object with constant scale factor
// ========================================================

inline void NNLO::rescaleGraph(TGraph* graph, const double scaleFactor){

    double x,y;
    int j;

    for (j=0;j<graph->GetN();j++){
        graph->GetPoint(j,x,y);
        graph->SetPoint(j,x,y/scaleFactor);
    }
}

// =================================================
//  Create rebinned histogram from normalized graph
// =================================================

void NNLO::createRebinnedHisto(const TGraph* graph, TH1* histo, const std::vector<double>& vecBinning){

    double xLeft, xRight;
    double yLeft, yRight;

    for (std::vector<double>::size_type vecIndex = 1;vecIndex != vecBinning.size();vecIndex++){

        double lowEdge  = vecBinning[vecIndex-1];
        double highEdge = vecBinning[vecIndex];
        double newBin   = 0.5*(highEdge+lowEdge);

        double binSum_graph1 = 0.0;
        double binSum_graph2 = 0.0;

        double area     = 0.0;
        double binratio = 0.0;

        // distribute bin integral from bin where new bin center falls between two data points 
        // according to ratio of left bin edge, bin center, right bin edge 
        int NPoints = 0;
        if(Theory == "kidonakis") NPoints = graph->GetN()-1;
        else if (Theory == "ahrens") NPoints = graph->GetN();
        for (int j=0;j<NPoints;j++){

            graph->GetPoint(j,  xLeft, yLeft);
            graph->GetPoint(j+1,xRight,yRight);
            area = getLocalIntegral(xLeft,xRight,yLeft,yRight);

            if(Theory == "ahrens"){
                // Modify the binWidth in graph to accomodate Ahrens description
                //  eg: 'x' = center of a bin
                //      'y' = integral over that given bin
                double binWidth = std::fabs(xRight - xLeft);
                xRight = xLeft + 0.5 * binWidth;
                xLeft -= 0.5* binWidth;
            }

            if (xLeft>=lowEdge && xRight<=highEdge){
                binSum_graph1 += area;
            }
            else if (xLeft<highEdge && xRight>highEdge){
                binratio = fabs((highEdge-xLeft)/(xRight-xLeft));// taking case into account when new bin center falls between two data points
                binSum_graph2 += binratio * area;
            }
            else if (xLeft<lowEdge && xRight>lowEdge){
                binratio = fabs((xRight-lowEdge)/(xRight-xLeft));// taking case into account when new bin center falls between two data points
                binSum_graph2 += binratio * area;
            }
        }

        if (parVerbose>0){
            std::cout << std::endl;
            std::cout << "Bin Content, Graph: " << newBin << " ---- " << binSum_graph1 << ", " << binSum_graph1/(highEdge-lowEdge) << std::endl;
            std::cout << "Bin Content, Graph: " << newBin << " ---- " << binSum_graph1+binSum_graph2 << ", " << binSum_graph1+binSum_graph2/(highEdge-lowEdge) << std::endl;
        }
        histo->SetBinContent(vecIndex,(binSum_graph1+binSum_graph2)/(highEdge-lowEdge));
    }
    histo->Scale(1./histo->Integral("width"));
}

// ========================================================
//  Calculate errors in y and assign them to TH1-object
// ========================================================

void NNLO::setBinErrors(TGraphAsymmErrors* graph, const TH1* histoCentral, const TH1* histoDn, const TH1* histoUp){

    double ey_dn = 0.0;
    double ey_up = 0.0;
    //double ey_total = 0.0;

    double error_min = 0.0;
    double error_max = 0.0;

    for (int j=0;j<histoCentral->GetNbinsX();j++){

        ey_dn = histoDn->GetBinContent(j+1) - histoCentral->GetBinContent(j+1);
        ey_up = histoUp->GetBinContent(j+1) - histoCentral->GetBinContent(j+1);

        if(parVerbose>0) std::cout << j+1 << ": " << histoCentral->GetBinCenter(j+1) << " .... " << ey_dn << " .... " << ey_up << std::endl;

        error_min=TMath::Min( TMath::Min(ey_dn,ey_up), 0.0);
        error_max=TMath::Max( TMath::Max(ey_dn,ey_up), 0.0);

        graph -> SetPointEYlow(j,error_min);
        graph -> SetPointEYhigh(j,error_max);

        //ey_total = sqrt(pow(ey_dn,2)+pow(ey_up,2));
        //histoCentral->SetBinError(j+1,ey_total);
    }
}

// =================================================
//  Generic style function for TGraph-objects
// =================================================

void NNLO::makeGraphStyle(TGraph* graph, int color, int markerStyle, TString name, TString title,
                          TString xAxisTitle, TString xAxisUnit, TString yAxisTitle, TString yAxisUnit)
{
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(markerStyle);
    graph->SetLineWidth(2);
    graph->SetName(name);
    graph->SetTitle(title);
    graph->GetXaxis()->SetTitle(xAxisTitle+" "+xAxisUnit);
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->SetTitle(yAxisTitle+" "+yAxisUnit);
    graph->GetYaxis()->SetTitleOffset(1.5);
}

// ==========================================
//  Generic style function for TH1-objects
// ==========================================

void NNLO::makeHistoStyle(TH1* histo, int color, int markerStyle, TString name, TString title, 
                          TString xAxisTitle, TString xAxisUnit, TString yAxisTitle, TString yAxisUnit)
{
    histo->SetLineColor(color);
    histo->SetMarkerColor(color);
    histo->SetMarkerStyle(markerStyle);
    histo->SetLineWidth(2);
    histo->SetName(name);
    histo->SetTitle(title);
    histo->GetXaxis()->SetTitle(xAxisTitle+" "+xAxisUnit);
    histo->GetXaxis()->CenterTitle();
    histo->GetYaxis()->SetTitle(yAxisTitle+" "+yAxisUnit);
    histo->GetYaxis()->SetTitleOffset(1.5);
}

// =================================================
//  Generic style function for TLegend-objects
// =================================================

void NNLO::makeLegendStyle(TLegend *leg){

    leg->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength());
    leg->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.1);
    leg->SetX2NDC(   gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.2);
    leg->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength());

    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextAlign(12);
    leg->SetTextSize(0.03);
}

// ===================================
//  generic style of TCanvas-objects
// ===================================

void NNLO::makeCanvasStyle(TCanvas *canvas, bool logXScale){

    canvas->cd();
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->SetLogx(logXScale);
}

// =====================
//  Main Function
// =====================

int NNLO::runNNLO(TString file_open_option){

    rescaleGraph(graphNormalizedNNLO,         getNormalization(graphUnnormalizedNNLO));
    rescaleGraph(graphNormalizedNNLO_ScaleDn, getNormalization(graphUnnormalizedNNLO_ScaleDn));
    rescaleGraph(graphNormalizedNNLO_ScaleUp, getNormalization(graphUnnormalizedNNLO_ScaleUp));

    createRebinnedHisto(graphNormalizedNNLO,         histoRebinnedNormalizedNNLO,         vecBins);
    createRebinnedHisto(graphNormalizedNNLO_ScaleDn, histoRebinnedNormalizedNNLO_ScaleDn, vecBins);
    createRebinnedHisto(graphNormalizedNNLO_ScaleUp, histoRebinnedNormalizedNNLO_ScaleUp, vecBins);

    graphRebinnedNormalizedNNLO_ErrorBand =new TGraphAsymmErrors((TH1F*)histoRebinnedNormalizedNNLO->Clone("histoRebinnedNormalizedNNLO_ErrorBand"));
    graphRebinnedNormalizedNNLO_ErrorBand->SetLineColor(2);
    graphRebinnedNormalizedNNLO_ErrorBand->SetFillStyle(3001);
    graphRebinnedNormalizedNNLO_ErrorBand->SetFillColor(kGray);

    setBinErrors(graphRebinnedNormalizedNNLO_ErrorBand,histoRebinnedNormalizedNNLO,histoRebinnedNormalizedNNLO_ScaleDn,histoRebinnedNormalizedNNLO_ScaleUp);

    fillLegends();

    printControlOutput();
    drawControlPlots();
    if (file_open_option!="") writeOutputFile(file_open_option);

    return 0;
}

#endif
