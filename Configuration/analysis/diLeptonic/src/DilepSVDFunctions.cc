# include "DilepSVDFunctions.h"
#include <TSystem.h>





// Namespaces
using namespace std;




////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////    Function Definitions  ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////


DilepSVDFunctions::DilepSVDFunctions()
{
    outpath = "";
}

    
// SVD Helper Function
// Return String with Steering
TString DilepSVDFunctions::SVD_GetSteering(TString channel, TString particle, TString quantity, TString special, TString syst)
{
    
    
    
    // Determine the Steering 
    int flag_regmode = 4;            // Always read in tau from file (4) 
    int flag_regpar = 2;             // Always use Tau (2)
    int flag_scan = 1;               // Turn scan on (2) or off (1) 
    int flag_ps = 4;                 // Write everything to PS (4)
    int flag_root = 1;               // Do not make a ROOT file (1)
    int flag_text = 2;               // Do write Text file with optimal tau (2)
    int flag_verbose = 2;            // Medium Verbosity (2)
    int flag_scanpoints = 3;         // Use 125 Scan Points (3)
    int flag_scanrange = 2;          // Use Factor 100 for Scan Range (2)
    int flag_lowersidebin = 3;       // Cut out side bin on rec, not on gen level (3)
    int flag_uppersidebin = 3;       // Cut out side bin on rec, not on gen level (3)
    int flag_matrixorientation = 1;  // Transpose matrix prior to unfolding (1)
    int flag_norm = 3;               // Extrinsic Normalization (1)
    int flag_closure = 1;            // Turn off closure test (1) 
    int flag_preweighting = 1;       // Reweighting prior to unfolding off (1)
    
    // Cut out the lower gen level side bin for the PT quantities
    TString concatenation = DilepSVDFunctions::SVD_CPQSS(channel, particle, quantity, special, ""); 
    if ( concatenation.EndsWith("LepPair_Pt") == true ) flag_lowersidebin = 5;
    if ( concatenation.EndsWith("TopQuarks_Pt") == true ) flag_lowersidebin = 5;
    if ( concatenation.EndsWith("TtBar_Pt") == true ) flag_lowersidebin = 5;
    
    
    // Make preweighting iterations
    //if ( concatenation.EndsWith("LepPair_Pt") == true ) flag_preweighting = 2;
    //if ( concatenation.EndsWith("TopQuarks_Pt") == true ) flag_preweighting = 2;
    //if ( concatenation.EndsWith("TtBar_Pt") == true ) flag_preweighting = 2;
    //if ( concatenation.EndsWith("BJets_Pt") == true ) flag_preweighting = 4;
    
    
    // Fix Side bins for Jet Quantities
    // if ( concatenation.EndsWith("BJets_Pt") == true ) flag_lowersidebin = 4;
    // if ( concatenation.EndsWith("BJets_Pt") == true ) flag_uppersidebin = 4;
    // if ( concatenation.EndsWith("BJets_Rapidity") == true ) flag_lowersidebin = 4;
    // if ( concatenation.EndsWith("BJets_Rapidity") == true ) flag_uppersidebin = 4;
    
    // Turn off scan if this is a systematic sample 
    if ( syst.CompareTo("") != 0 ) {
    	flag_scan = 1 ;
    } 
    
     
    // Calculate Steering Parameter
    TString steering = "";
    steering.Prepend(TString::Format("%i", flag_regmode));                     //1
    steering.Prepend(TString::Format("%i", flag_regpar));                      //2
    steering.Prepend(TString::Format("%i", flag_scan));                        //3
    steering.Prepend(TString::Format("%i", flag_ps));                          //4
    steering.Prepend(TString::Format("%i", flag_root));                        //5
    steering.Prepend(TString::Format("%i", flag_text));                        //6
    steering.Prepend(TString::Format("%i", flag_verbose));                     //7
    steering.Prepend(TString::Format("%i", flag_scanpoints));                  //8
    steering.Prepend(TString::Format("%i", flag_scanrange));                   //9
    steering.Prepend(TString::Format("%i", flag_lowersidebin));                //10
    steering.Prepend(TString::Format("%i", flag_uppersidebin));                //11
    steering.Prepend(TString::Format("%i", flag_matrixorientation));           //12 
    steering.Prepend(TString::Format("%i", flag_norm));                        //13 
    steering.Prepend(TString::Format("%i", flag_closure));                     //14 
    steering.Prepend(TString::Format("%i", flag_preweighting));                //15 
    
    // Return it
    return steering;
}

    

// DAVID
// SVD Unfolding Helper Function
// Create Tex Snippets for the histograms
void DilepSVDFunctions::SVD_Tex(TString channel, TString particle, TString quantity, TString special, TString& channelTex, TString& particleTex, TString& quantityTex, TString& specialTex)
{

    channelTex = "";
    if ( strcmp(channel, "mumu") == 0 ) {
        channelTex.Append("#mu^{+}#mu^{-}");
    } else if ( strcmp(channel, "emu") == 0 ) {
        channelTex.Append("e^{+}#mu^{-}");
    } else if ( strcmp(channel, "ee") == 0 ) {
        channelTex.Append("e^{+}e^{-}");
    } else if ( strcmp(channel, "combined") == 0 ) {
        channelTex.Append("cb.");
    } else {
        channelTex.Append(channel);
    }
    
    particleTex = "";
    if ( strcmp(particle, "Leptons") == 0 ) {
        particleTex.Append("l^{+} / l^{-}");
    } else if ( strcmp(particle, "LepPair") == 0 ) {
        particleTex.Append("(l^{+}l^{-})");
    } else if ( strcmp(particle, "Jets") == 0 ){
        particleTex.Append("Jets");
    } else if ( strcmp(particle, "BJets") == 0 ){
        particleTex.Append("BJets");
    } else if ( strcmp(particle, "BBbar") == 0 ){
        particleTex.Append("b#bar{b}");
    } else if ( strcmp(particle, "TopQuarks") == 0 ) {
        particleTex.Append("t / #bar{t}");
    } else if ( strcmp(particle, "TtBar") == 0 ) {
        particleTex.Append("(t#bar{t})");
    } else {
        particleTex.Append(particle);
    } 
    
    quantityTex = "";
    if ( strcmp(quantity, "Pt") == 0 ) {
        quantityTex.Append("p_{T}");
    } else if ( strcmp(quantity, "Eta") == 0 ) {
        quantityTex.Append("#eta");
    } else if ( strcmp(quantity, "Rapidity") == 0 ) {
        quantityTex.Append("y");
    } else if ( strcmp(quantity, "Mass") == 0 ) {
        quantityTex.Append("M_{t #bar{t}}");
    } else {
        quantityTex.Append(quantity);
    } 
    
    specialTex = "";
    if ( strcmp(special, "MoreBins") == 0 ) {
        specialTex.Append("(MB)"); 
    } else {
        specialTex.Append(particle);
    }

} 

// DAVID
// SVD Unfolding Helper Function
TString DilepSVDFunctions::SVD_GetOutputPath()
{

    TString outputfilename = "";
    outputfilename.Append(outpath);
    outputfilename.Append("SVD"); 

    return outputfilename;

}

// DAVID
// SVD Unfolding Helper Function
// Get complete Output File Name for Txt
TString DilepSVDFunctions::SVD_GetOutputFileNameTxt()
{

    TString outputfilename = SVD_GetOutputPath();
    outputfilename.Append("/optimalTau"); 
    outputfilename.Append(".txt"); 
    
    return outputfilename;
     
}

// DAVID
// SVD Unfolding Helper Function
// Get complete Output File Name for Ps
TString DilepSVDFunctions::SVD_GetOutputFileNameEps(TString channel, TString particle, TString quantity, TString syst, TString special)
{

    TString outputfilename = SVD_GetOutputPath();
    outputfilename.Append("/Unfolding_");
    outputfilename.Append(SVD_CPQSS(channel, particle, quantity, special, syst));  
    outputfilename.Append(".eps"); 
    return outputfilename;

}

// DAVID
// SVD Unfolding Helper Function
// Get complete Output File Name for Ps
TString DilepSVDFunctions::SVD_GetOutputFileNamePs(TString channel, TString particle, TString quantity, TString syst, TString special)
{

    TString outputfilename = SVD_GetOutputPath();
    outputfilename.Append("/Unfolding_");
    outputfilename.Append(SVD_CPQSS(channel, particle, quantity, special, syst));  
    outputfilename.Append(".ps"); 
    return outputfilename;

}


// DAVID
// SVD Unfolding Helper Function
// Get complete Output File Name for Root
TString DilepSVDFunctions::SVD_GetOutputFileNameRoot(TString channel, TString particle, TString quantity, TString syst, TString special)
{
    TString outputfilename = SVD_GetOutputPath();
    outputfilename.Append("/Unfolding_");
    outputfilename.Append(SVD_CPQSS(channel, particle, quantity, special, syst)); 
    outputfilename.Append(".root"); 
    return outputfilename;

}
 
 

 
// Unfolding with systematics.
// This takes "3 of each kind": Nominal, Up and Down
// The return value will contain absolute shifts !!!
// For more information see function SVD_DoUnfold() 
double DilepSVDFunctions::SVD_DoUnfoldSys(
        TH1D* dataInputHist,                              // DATA.
        TH1D* bgrNom, TH1D* bgrUp, TH1D* bgrDown,         // Background.
        TH1D* ttbgrNom, TH1D* ttbgrUp, TH1D* ttbgrDown,   // Tt-Background.
        TH1D* genNom, TH1D* genUp, TH1D* genDown,         // Gen. Note: if genUp==NULL then genUp=genNom, etc. 
        TH1D* recNom, TH1D* recUp, TH1D* recDown,         // Rec
        TH2D* respNom, TH2D* respUp, TH2D* respDown,      // Response Matrices
        double totalDataEventsNom,
        double totalBgrEventsNom,  double totalBgrEventsUp,  double totalBgrEventsDown, 
        double totalTtBgrEventsNom,  double totalTtBgrEventsUp,  double totalTtBgrEventsDown, 
        double totalRecEventsNom,  double totalRecEventsUp,  double totalRecEventsDown, 
        double totalGenEventsNom,  double totalGenEventsUp,  double totalGenEventsDown,     
        const double thebins[],                 
        const int numbins,                       
        TH1D*& shifts,                        
        TString channel,                       
        TString particle,                      
        TString quantity,                      
        TString special,    
        TString syst                           // Specify Name for the Down-Systematic
)
{           
    // Align Data
    TH1D* arrData = new TH1D[1];
    *(arrData+0) = *dataInputHist; 
     
    
    // Align Background
    TH1D* arrBgr = new TH1D[3];
    *(arrBgr+0) = *bgrNom;
    *(arrBgr+1) = *bgrUp;
    *(arrBgr+2) = *bgrDown;
     
    
    // Align Background
    TH1D* arrTtBgr = new TH1D[3];
    *(arrTtBgr+0) = *ttbgrNom;
    *(arrTtBgr+1) = *ttbgrUp;
    *(arrTtBgr+2) = *ttbgrDown;
    
    
    // Align Gen Dists
    TH1D* arrGen =  new TH1D[3];
    *(arrGen+0) = *genNom;
    *(arrGen+1) = *genUp;
    *(arrGen+2) = *genDown; 
         
    
    // Align Rec Dists
    TH1D* arrRec = new TH1D[3];
    *(arrRec+0) = *recNom;
    *(arrRec+1) = *recUp;
    *(arrRec+2) = *recDown;
    
    
    // Align Response Matrices
    TH2D* arrResp = new TH2D[3];
    *(arrResp+0) = *respNom;
    *(arrResp+1) = *respUp;
    *(arrResp+2) = *respDown;
    
    // Align Normalization Input
    double totalDataEvents[1]  = {totalDataEventsNom};
    double totalBgrEvents[3]   = {totalBgrEventsNom, totalBgrEventsUp, totalBgrEventsDown};
    double totalTtBgrEvents[3] = {totalTtBgrEventsNom, totalTtBgrEventsUp, totalTtBgrEventsDown};
    double totalRecEvents[3]   = {totalRecEventsNom, totalRecEventsUp, totalRecEventsDown};
    double totalGenEvents[3]   = {totalGenEventsNom, totalGenEventsUp, totalGenEventsDown};
         
    // Create Pointer for Output 
    TH1D* arrUnfolded = NULL;
    TH1D* arrUnfoldedNorm = NULL;
    
    
    // Call the unfolding function
    int numSys = 2;
    double bestTau = SVD_DoUnfold(
        arrData,               
        arrBgr,               
        arrTtBgr,                       
        arrGen,                      
        arrRec,                    
        arrResp, 
        totalDataEvents, 
        totalBgrEvents,  
        totalTtBgrEvents,  
        totalRecEvents, 
        totalGenEvents,                    
        thebins,                 
        numbins,                       
        arrUnfolded,      // Return value 
        arrUnfoldedNorm,  // Return value
        numSys,          
        channel,                     
        particle,                      
        quantity,                      
        special,
        syst);                    
     
     
    // Get the shifts histograms 
    // Get them from the NORMALIZED distributions!
    TH1D* unfNom = arrUnfoldedNorm;  
    TH1D* unfUp = arrUnfoldedNorm + 1; 
    TH1D* unfDown = arrUnfoldedNorm + 2; 
     
     
    TCanvas *c = new TCanvas("Ivan", "Ivan");
    TLegend * leg =  new TLegend(0.70,0.55,0.98,0.85);
    TH1 *Nominal= (TH1*)unfNom->Clone("Nominal");
    TH1 *NominalUp= (TH1*)unfUp->Clone("NominalUp");
    TH1 *NominalDown= (TH1*)unfDown->Clone("NominalDown");
    Nominal->SetLineColor(kBlack);
    Nominal->SetMarkerStyle(20);
    Nominal->SetMarkerColor(kBlack);
    NominalUp->SetLineColor(kRed);
    NominalUp->SetMarkerStyle(22);
    NominalUp->SetMarkerColor(kRed);
    NominalDown->SetLineColor(kBlue);
    NominalDown->SetLineStyle(2);
    NominalDown->SetMarkerColor(kBlue);
    NominalDown->SetMarkerStyle(23);
    leg->AddEntry(Nominal, "Nominal", "p");
    leg->AddEntry(NominalUp, "Up", "lp");
    leg->AddEntry(NominalDown, "Down", "lp");
    Nominal->DrawNormalized();
    NominalUp->DrawNormalized("h,same");
    NominalDown->DrawNormalized("h,same");
    leg->Draw("same");
    c->Print(TString("Plots/").Copy().Append(channel).Append("/SystComparison_").Append(channel).Append("_").Append(particle).Append("_").Append(quantity).Append("_").Append(syst).Append(".eps"));
    leg->Clear();
    Nominal->Clear(); delete Nominal;
    NominalUp->Clear(); delete NominalUp;
    NominalDown->Clear(); delete NominalDown;
    c->Clear();
    
    
    // Create the output histogram for this function
    if ( shifts == NULL ) {
        shifts =  (TH1D*) unfNom->Clone("shifts");
    }
     
    
    // Calculate the bin wise differences
    for ( int i = 1 ; i <= numbins+2 ; i++ ) {
        double valunfNom = unfNom->GetBinContent(i);
        double valunfUp = unfUp->GetBinContent(i);
        double valunfDown = unfDown->GetBinContent(i);  
        
        // Calculate average of up and down shift
        double Sys_Error_Up   = std::fabs(valunfUp   - valunfNom);
        double Sys_Error_Down = std::fabs(valunfDown - valunfNom);
        
        double Sys_Error = 0.5* (Sys_Error_Up + Sys_Error_Down) ;
       
        // Make a relative uncertainty from that
        double Sys_Error_Relative = SVD_Divide(Sys_Error, valunfNom);

        // Save it
        shifts->SetBinContent(i, Sys_Error_Relative);
        shifts->SetBinError(i, 0.);
    }
      
    delete[] arrData;
    delete[] arrBgr;
    delete[] arrGen;
    delete[] arrRec;
    delete[] arrResp;
    delete[] arrTtBgr;
    delete[] arrUnfolded;
    delete[] arrUnfoldedNorm;
    delete leg;
    delete c;
    
    // return best tau value
    return bestTau;
        
}    
 
 
// PERFORM UNFOLDING  
double DilepSVDFunctions::SVD_DoUnfold(
        TH1D* dataInputHist,                     // Data Input (RAW Data including the background)
        TH1D* bgrInputHist,                      // Background (will be substracted from data)
        TH1D* ttbgrInputHist,                    // Tt-Background 
        TH1D* genInputHist,                      // Generated MC
        TH1D* recInputHist,                      // Reconstructed MC
        TH2D* respInputHist,                     // Response Matrix 
        double* totalDataEvents,                 // For normalization. If zet to NULL, will be calculated from input distributions
        double* totalBgrEvents,                  // For normalization. If zet to NULL, will be calculated from input distributions 
        double* totalTtBgrEvents,                // For normalization. If zet to NULL, will be calculated from input distributions 
        double* totalRecEvents,                  // For normalization. If zet to NULL, will be calculated from input distributions
        double* totalGenEvents,                  // For normalization. If zet to NULL, will be calculated from input distributions
        const double thebins[],                  // Binning for the unfolding
        const int numbins,                       // Number of bins for unfolding (not counting OF bins !!!) 
        TH1D*& unfolded,                         // Returned: Unfolded Distribution.
        TH1D*& unfoldedNorm,                         // Returned: Unfolded Distribution.
        const int numSys,                        // Number of used systematics (standard = 0 )
        TString channel,                         // Specify Name for the Channel ("mumu", "emu", "ee" ...)
        TString particle,                        // Specify Name for the Physics Object ("Top", "BJets", "Leptons")
        TString quantity,                        // Specify Name for the Quantity ("Eta", "Pt", or "Mass");
        TString special,                         // Special Postfix to be used for names of files and plots
        TString syst                             // Specify Name for the Up-Systematic
)   
{ 

    
    
    
    // Determine the Steering 
    TString steering = SVD_GetSteering(channel, particle, quantity, special, syst);  
     
    
    // Filenames for output 
    TString rootFile = SVD_GetOutputFileNameRoot(channel, particle, quantity, syst, special);
    TString psFile = SVD_GetOutputFileNamePs(channel, particle, quantity, syst, special);
    TString epsFile = SVD_GetOutputFileNameEps(channel, particle, quantity, syst, special);
    TString txtFile =  SVD_GetOutputFileNameTxt(); 
    
    // Filenames for input
    TString regParFile(gSystem->pwd());
    regParFile.Append("/optimalTau.txt");

  
    // TexStrings for Formatting
    TString channelTex = "";
    TString particleTex = "";
    TString quantityTex = ""; 
    TString specialTex = "";
    SVD_Tex(channel, particle, quantity, special, channelTex, particleTex, quantityTex, specialTex);


    // Tex for Systematics
    TString systTex = syst; 
    
    
    // Dummy Tau Parameter
    // will be ignored, since we read regPar from file 
    double regPar = 1.; 


    // Do the Unfolding 
    double optimalTauX = TopSVDFunctions::SVD_Unfold(dataInputHist, bgrInputHist, ttbgrInputHist, genInputHist, recInputHist, respInputHist,
        unfolded, unfoldedNorm,
        totalDataEvents, totalBgrEvents, totalTtBgrEvents, totalRecEvents, totalGenEvents,
        thebins, numbins, regPar, steering, numSys,
        channel, particle, quantity, special, syst,
        channelTex, particleTex, quantityTex, specialTex, systTex,
        rootFile, psFile, epsFile, txtFile, regParFile);
  
    // Return
    return optimalTauX;
    
}


 
    
