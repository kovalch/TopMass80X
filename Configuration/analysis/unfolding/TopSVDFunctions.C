#include <cstdio>
#include <fstream>
#include <TROOT.h>
#include <TSystem.h>
#include "TopSVDFunctions.h"



#ifdef STANDALONE_COMPILE
#include "../common/include/plotterUtils.h"
#else
#include "TopSVDUnfold.C"
#include "../semiLeptonic/diffXSection/basicFunctions.h"
#include "../semiLeptonic/diffXSection/HHStyle.h"
#endif





// Namespaces
using namespace std;



////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////    Function Definitions  ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////



// SVD Unfolding Helper Function
TString TopSVDFunctions::SVD_CPQSS(TString channel, TString particle, TString quantity, TString special, TString syst)
{
    TString cpq = ""; 
    cpq.Append(channel);
    cpq.Append("_");
    cpq.Append(particle);
    cpq.Append("_");
    cpq.Append(quantity);
    if ( special.CompareTo("") != 0 ) {
        if ( special.BeginsWith("_") == true ) special.Remove(0,1); 
        cpq.Append("_");
        cpq.Append(special);
    }
    if ( syst.CompareTo("") != 0 ) {
        if ( syst.BeginsWith("_") == true ) syst.Remove(0,1); 
        cpq.Append("_");
        cpq.Append(syst);
    }
    return cpq;
}


// Create a name for a plot according to a certain scheme
// This name will be used for saving the histogram in a root file
TString TopSVDFunctions::SVD_PlotName(TString channel, TString particle, TString quantity, TString special, TString syst, TString plotName)
{
    TString rawStr("SVD_");
    rawStr.Append(SVD_CPQSS(channel, particle, quantity, special, syst));
    rawStr.Append("_");
    rawStr.Append(plotName); 

    return rawStr;
}

 
// Division of positive doubles
// If an argument is smaller or equal to zero,
// the result is zero.
double TopSVDFunctions::SVD_Divide(double numerator, double denominator)
{
    double result = 0.;
    if ( numerator > 0. && denominator > 0. ) {
        result = numerator / denominator;
    }
    return result;
}


// Calculates binomial errors                  
double TopSVDFunctions::SVD_BinomialError(double numerator, double denominator)
{
    if (numerator>denominator) {
        return 0.;
        cout << "WARNING in BinomialError: efficiency is larger than one!" << endl;
    }

    double  error;
    error=TMath::Sqrt((1-numerator/denominator)
                      *numerator/denominator/denominator);

    return error;
}



// Sqrt of positive doubles
// If an argument is smaller or equal to zero,
// the result is zero.
double TopSVDFunctions::SVD_Sqrt(double square)
{
    double result = 0.;
    if ( square > 0. ) {
        result = TMath::Sqrt(square);
    }

    return result;
}





// Delete 1D histograms safely
// ATTENTION: May only be called on arrays  
// (works of course also on 1-entry-arrays.)
void TopSVDFunctions::SVD_DeleteHists1D(TH1D*& hist, int numHist)
{
    // Existence of Objects
    if ( hist == NULL ) return;  
 
    // Full Array 
    for ( int h = 0 ; h < numHist ; h++ ) {
        (hist+h)->SetDirectory(NULL);
    }  
    delete[] hist; 
    hist = NULL;
}


// Delete 2D histograms safely
// ATTENTION: May only be called on arrays  
// (works of course also on 1-entry-arrays.)
void TopSVDFunctions::SVD_DeleteHists2D(TH2D*& hist, int numHist)
{ 
    // Existence of Objects
    if ( hist == NULL ) return;  
 
    // Full Array 
    for ( int h = 0 ; h < numHist ; h++ ) {
        (hist+h)->SetDirectory(NULL);
    }  
    delete[] hist; 
    hist = NULL;
}

 

// Helper Function to place legend
void TopSVDFunctions::SVD_FormatLegendStandard(TLegend* leg) {
 
   
    leg->SetX1NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength()-0.25);
    leg->SetY1NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength()-0.15);
    leg->SetX2NDC(1.0-gStyle->GetPadRightMargin()-gStyle->GetTickLength());
    leg->SetY2NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength());
    leg->SetTextFont(42);
    leg->SetTextSize(0.03); 
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextAlign(12);   
} 


// Helper Function to place legend 
void TopSVDFunctions::SVD_FormatLegendKScan(TLegend* leg) {
 
    leg->SetX1NDC(1.0-gStyle->GetPadRightMargin());
    leg->SetY1NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength()-0.5);
    leg->SetX2NDC(1.0);
    leg->SetY2NDC(1.0-gStyle->GetPadTopMargin()-gStyle->GetTickLength());
    leg->SetTextFont(42);
    leg->SetTextSize(0.03); // 0.04
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextAlign(12);
 

} 


// Helper Function
// This clears a given stack and a given legend, 
// and resets all the formating to the standard values.
// ATTENTION:
// (1) This deletes ALL histograms in the stack! 
//     This is OK, as long as you use the helper function
//     SVD_Array2Stack() to add things to your stack, because
//     this function makes copies first and adds them afterwards.
// (2) If NULL pointers are given for the stack or the leg,
//     new objects will be created 
void TopSVDFunctions::SVD_ClearStackLeg(THStack*& stack, TLegend*& leg, TString header1, TString header2, TString header3) {

    
    // Create new Objects if needed
    if ( leg   == NULL ) leg   = new TLegend();
    if ( stack == NULL ) stack = new THStack("stack","stack");  



    // Clear the legend
    leg->Clear();


    // Standard Format
    SVD_FormatLegendStandard(leg); 
    
    
    // Header, 1. line
    TString theHeader1 = header1;
    if ( theHeader1.CompareTo("") != 0 ) {
        leg->AddEntry((TObject*) NULL, theHeader1, "");
    } 
    
    
    // Header, 2. line
    TString theHeader2 = header2;
    if ( theHeader2.CompareTo("") != 0 ) {
        leg->AddEntry((TObject*) NULL, theHeader2, "");
    }
    
    
    // Header, 3. line
    TString theHeader3 = header3;
    if ( theHeader3.CompareTo("") != 0 ) {
        leg->AddEntry((TObject*) NULL, theHeader3, "");
    }
    
    
    // Clear the list of histograms in the stack
    // Set the "kIsOwner"-Flag to really remove all the histograms 
    TList* theHists = stack->GetHists();
    if ( theHists != NULL ) {
        theHists->SetOwner(true); 
        theHists->Delete();
    }  
    
    // Clear the list of SUM of histograms in the stack
    // Set the "kIsOwner"-Flag to really remove all the histograms
    TObjArray* theStack = stack->GetStack();
    if ( theStack != NULL ) {
        theStack->SetOwner(true);
        theStack->Delete(); 
    } 
    
     
    // Reset the histogram
    TH1* theHistogram = stack->GetHistogram();
    if ( theHistogram != NULL ) {
        theHistogram->Reset("M"); 
    }
    
    // Set Standard Values 
    stack->SetMaximum(-1111.);
    stack->SetMinimum(-1111.);
    stack->SetNameTitle("","");
     

} 

// Clone Array of 1D Hists
// Attention: New objects will be created on the heap.
// Do not fortget to delete them somewhere
TH1D* TopSVDFunctions::SVD_CloneHists1D(TH1D* hist, int numHist)
{
     
    // Existence of Objects
    if ( hist == NULL ) return NULL;  
 
    // Create new objects
    TH1D* outputhists = new TH1D[numHist];   
    for ( int h = 0 ; h < numHist ; h++ ) { 
        
        // Copy the histograms  
        *(outputhists+h) = *(hist+h);    
        
        // Make sure the new histogram is not on a TDirectory
        // You will most likely run into garbage collector problems,
        (outputhists+h)->SetDirectory(NULL);
    }  
    
    // return new array
    return outputhists;
}

// Clone Array of 2D Hists
// Attention: New objects will be created on the heap.
// Do not fortget to delete them somewhere
TH2D* TopSVDFunctions::SVD_CloneHists2D(TH2D* hist, int numHist)
{
     
    // Existence of Objects
    if ( hist == NULL ) return NULL;  
    
    // Create new objects
    TH2D* outputhists = new TH2D[numHist];
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Copy the histograms
        (hist+h)->Copy(*(outputhists+h)); 
        
        // Make sure the new histogram is not on a TDirectory
        // You will most likely run into garbage collector problems,
        (outputhists+h)->SetDirectory(NULL);
    }
    
    // return new array
    return outputhists;
}

// Move Events in OF/UF bins in dedicated side bins.
// Attention: Sidebins need to be created first
// Error recalculated in Poisson manner.
// Can handle arrays of histograms
void TopSVDFunctions::SVD_MoveOFBins1D(TH1D* hist, int numHist)
{
     
    // Existence of Objects
    if ( hist == NULL ) return;  
    
    int nBins = hist->GetNbinsX();
    
    // Loop over all histograms in the array
    for ( int h = 0 ; h < numHist ; h++ ) {

        // Count OF Events 
        double ufBin = (hist+h)->GetBinContent(0);
        ufBin += (hist+h)->GetBinContent(1);
        double ofBin = (hist+h)->GetBinContent(nBins+1);
        ofBin += (hist+h)->GetBinContent(nBins);
    
        // Set Side Bins
        (hist+h)->SetBinContent(0, 0.);
        (hist+h)->SetBinContent(1, ufBin);
        (hist+h)->SetBinContent(nBins, ofBin);
        (hist+h)->SetBinContent(nBins+1, 0.);
    
        // Set Side Bins / Error
        (hist+h)->SetBinError(0, 0.);
        (hist+h)->SetBinError(1, 0.);
        (hist+h)->SetBinError(nBins, 0.);
        (hist+h)->SetBinError(nBins+1, 0.);
        if ( ufBin > 0. ) {
            (hist+h)->SetBinError(1, TMath::Sqrt(ufBin));
        }
        if ( ofBin > 0. ) {
            (hist+h)->SetBinError(nBins, TMath::Sqrt(ofBin));
        }
    }
}




// Move Events in OF/UF bins in dedicated side bins
// Attention: Sidebins need to be created first
// Error recalculated in Possion manner
// Can handle arrays of histograms
void TopSVDFunctions::SVD_MoveOFBins2D(TH2D* hist, int numHist){
    
    // Existence of Objects
    if ( hist == NULL ) return;  
    
    int nBinsx = hist->GetNbinsX();
    int nBinsy = hist->GetNbinsY();
    
    // Loop over all histograms in the array 
    for ( int h = 0 ; h < numHist ; h++ ) {
    
    
        // Loop over bins in X direction
        for ( int i = 0 ; i < nBinsx+2 ; i++ ) {
    
             // Count OF Events
            double ufBin = (hist+h)->GetBinContent(i,0);
            ufBin += (hist+h)->GetBinContent(i,1);
            double ofBin = (hist+h)->GetBinContent(i,nBinsy+1);
            ofBin += (hist+h)->GetBinContent(i,nBinsy);
        
            // Set Side Bins
            (hist+h)->SetBinContent(i, 0, 0.);
            (hist+h)->SetBinContent(i, 1, ufBin);
            (hist+h)->SetBinContent(i, nBinsy, ofBin);
            (hist+h)->SetBinContent(i, nBinsy+1, 0.);
        
            // Set Side Bins / Error
            (hist+h)->SetBinError(i, 0, 0.);
            (hist+h)->SetBinError(i, 1, 0.);
            (hist+h)->SetBinError(i, nBinsy, 0.);
            (hist+h)->SetBinError(i, nBinsy+1, 0.);
            if ( ufBin > 0. ) {
                (hist+h)->SetBinError(i, 1, TMath::Sqrt(ufBin));
            }
            if ( ofBin > 0. ) {
                (hist+h)->SetBinError(i, nBinsy, TMath::Sqrt(ofBin));
            }
        }
    
        // Now, loop over bins in Y direction
        // The first and last row should be empty already
        // But this does not hurt
        for ( int i = 0 ; i < nBinsy+2 ; i++ ) {
    
             // Count OF Events
            double ufBin = (hist+h)->GetBinContent(0,i);
            ufBin += (hist+h)->GetBinContent(1,i);
            double ofBin = (hist+h)->GetBinContent(nBinsx+1,i);
            ofBin += (hist+h)->GetBinContent(nBinsx,i);
    
            // Set Side Bins
            (hist+h)->SetBinContent(0, i, 0.);
            (hist+h)->SetBinContent(1, i, ufBin);
            (hist+h)->SetBinContent(nBinsx,   i, ofBin);
            (hist+h)->SetBinContent(nBinsx+1, i, 0.);
    
            // Set Side Bins / Error
            (hist+h)->SetBinError(0, i, 0.);
            (hist+h)->SetBinError(1, i, 0.);
            (hist+h)->SetBinError(nBinsx, i, 0.);
            (hist+h)->SetBinError(nBinsx+1, i, 0.);
            if ( ufBin > 0. ) {
                (hist+h)->SetBinError(1, i, TMath::Sqrt(ufBin));
            }
            if ( ofBin > 0. ) {
                (hist+h)->SetBinError(nBinsx, i, TMath::Sqrt(ofBin));
            }
        }
    }
    
} 

// Set bin Content
void TopSVDFunctions::SVD_SetBinValErr1D(TH1D* histo, int bin, double value, double error, int numHist)
{
    // Existence of Objects
    if ( histo == NULL ) return;  
    int bins = histo->GetNbinsX();
    
    // Loop over all histograms in the array 
    for ( int h = 0 ; h < numHist ; h++ ) {
        (histo+h)->SetBinContent(bin, value);
        (histo+h)->SetBinError(bins, error);
    }
    
}

// This empties the dedicated side bins.
// Should be called for reconstructed quantities.
// Histogram side bins are also emptied.
// Can handle arrays of histograms
void TopSVDFunctions::SVD_EmptyRecSideBins1D(TH1D* histo, bool cutLowerRecSideBin, bool cutUpperRecSideBin, int numHist) 
{
    // Existence of Objects
    if ( histo == NULL ) return;  
    int bins = histo->GetNbinsX();
    
    // Loop over all histograms in the array 
    for ( int h = 0 ; h < numHist ; h++ ) {
        (histo+h)->SetBinContent(0, 0.);
        (histo+h)->SetBinContent(bins+1, 0.);
        (histo+h)->SetBinError(0, 0.);
        (histo+h)->SetBinError(bins+1, 0.);
        if ( cutLowerRecSideBin == true ) {
            (histo+h)->SetBinContent(1, 0.);
            (histo+h)->SetBinError(1, 0.);
        }
        if ( cutUpperRecSideBin == true ) {
            (histo+h)->SetBinContent(bins, 0.);    
            (histo+h)->SetBinError(bins, 0.);
        } 
    }
    
}


// This empties the dedicated side bins.
// Should be called for generated quantities.
// Histogram side bins are also emptied.
// Can handle arrays of histograms
void TopSVDFunctions::SVD_EmptyGenSideBins1D(TH1D* histo, bool cutLowerGenSideBin, bool cutUpperGenSideBin, int numHist) 
{
    // Existence of Objects
    if ( histo == NULL ) return;  
    int bins = histo->GetNbinsX();
    
    // Loop over all histograms in the array 
    for ( int h = 0 ; h < numHist ; h++ ) {
        (histo+h)->SetBinContent(0, 0.);
        (histo+h)->SetBinContent(bins+1, 0.);
        (histo+h)->SetBinError(0, 0.);
        (histo+h)->SetBinError(bins+1, 0.);
        if ( cutLowerGenSideBin == true ) {
            (histo+h)->SetBinContent(1, 0.);
            (histo+h)->SetBinError(1, 0.);
        }
        if ( cutUpperGenSideBin == true ) {
            (histo+h)->SetBinContent(bins, 0.);    
            (histo+h)->SetBinError(bins, 0.);
        } 
    }
    
}

// This empties the dedicated side bins.
// Should be called for reconstructed quantities,
// which are assumed to be on X axis.
// Histogram side bins are also emptied.
// Can handle arrays of histograms
void TopSVDFunctions::SVD_EmptyRecSideBins2D(TH2D* histo, bool cutLowerRecSideBin, bool cutUpperRecSideBin, int numHist)
{
    // Existence of Objects
    if ( histo == NULL ) return;   
    
    int binsx = histo->GetNbinsX();
    int binsy = histo->GetNbinsY();
    
    // Loop over all histograms in the array 
    for ( int h = 0 ; h < numHist ; h++ ) {

        // Loop over generated bins (on Y axis)
        for ( int j = 0 ; j < binsy+2 ; j++ ) { 
            (histo+h)->SetBinContent(0, j, 0.);
            (histo+h)->SetBinContent(binsx+1, j, 0.);
            (histo+h)->SetBinError(0, j, 0.);
            (histo+h)->SetBinError(binsx+1, j, 0.); 
            if ( cutLowerRecSideBin == true ) {
                (histo+h)->SetBinContent(1, j, 0.);
                (histo+h)->SetBinError(1, j, 0.);
            } 
            if ( cutUpperRecSideBin == true ) {
                (histo+h)->SetBinContent(binsx, j, 0.);
                (histo+h)->SetBinError(binsx, j, 0.);
            } 
        } 
    } 
}

// This empties the dedicated side bins.
// Should be called for genetared quantities,
// which are assumed to be on Y axis.
// Histogram side bins are also emptied.
// Can handle arrays of histograms
void TopSVDFunctions::SVD_EmptyGenSideBins2D(TH2D* histo, bool cutLowerGenSideBin, bool cutUpperGenSideBin, int numHist)
{
    // Existence of Objects
    if ( histo == NULL ) return;   
    
    int binsx = histo->GetNbinsX();
    int binsy = histo->GetNbinsY();
    
    // Loop over all histograms in the array 
    for ( int h = 0 ; h < numHist ; h++ ) {

        // Loop over reconstructed bins (on X axis)
        for ( int i = 0 ; i < binsx+2 ; i++ ) { 
            (histo+h)->SetBinContent(i, 0, 0.);
            (histo+h)->SetBinContent(i, binsy+1, 0.);
            (histo+h)->SetBinError(i, 0, 0.);
            (histo+h)->SetBinError(i, binsy+1, 0.);
            if ( cutLowerGenSideBin == true ) {
                (histo+h)->SetBinContent(i, 1, 0.);
                (histo+h)->SetBinError(i, 1, 0.); 
            } 
            if ( cutUpperGenSideBin == true ) {
                (histo+h)->SetBinContent(i, binsy, 0.);
                (histo+h)->SetBinError(i, binsy, 0.); 
            } 
        } 
    } 
}


// Add up histograms
// This function adds up the content of the right histo to the left histo.
// This means, the left histogram is CHANGED, the right one is not changed.
void TopSVDFunctions::SVD_AddRightToLeft(TH1D*& histoLeft, TH1D* histoRight, double factorLeft, double factorRight, int numHist)
{
     
    // Existence of Objects
    if ( histoLeft == NULL ) return;
    if ( histoRight == NULL ) return;
    
    
    // Create new objects 
    for ( int h = 0 ; h < numHist ; h++ ) { 
          
        int binsx = (histoLeft+h)->GetNbinsX();  
        for ( int i = 1 ; i <= binsx ; i++ ) {
            
            // Get Values
            double left = (histoLeft+h)->GetBinContent(i);
            double leftError = (histoLeft+h)->GetBinError(i);
            double right = (histoRight+h)->GetBinContent(i);
            double rightError = (histoRight+h)->GetBinError(i);
            
            
            // Calculate new Values
            double value = left*factorLeft+right*factorRight;
            double errorsq = factorLeft*factorLeft*leftError*leftError+factorRight*factorRight*rightError*rightError;
            double error = SVD_Sqrt(errorsq);
            
            
            // Save it
            (histoLeft+h)->SetBinContent(i, value);
            (histoLeft+h)->SetBinError(i, error); 
        } 
    } 
}

// Multiply histograms
// This function multiplies the content of the right histo to the left histo.
// This means, the left histogram is CHANGED, the right one is not changed.
void TopSVDFunctions::SVD_MultRightToLeft(TH1D*& histoLeft, TH1D* histoRight, double factorLeft, double factorRight, int numHist)
{
     
    // Existence of Objects
    if ( histoLeft == NULL ) return;
    if ( histoRight == NULL ) return;
    
    
    // Create new objects 
    for ( int h = 0 ; h < numHist ; h++ ) { 
          
        int binsx = (histoLeft+h)->GetNbinsX();  
        for ( int i = 1 ; i <= binsx ; i++ ) {
            
            // Get Values
            double left = (histoLeft+h)->GetBinContent(i);
            double leftError = (histoLeft+h)->GetBinError(i);
            double right = (histoRight+h)->GetBinContent(i);
            double rightError = (histoRight+h)->GetBinError(i); 
            
            
            // Calculate new Values
            double value = left*factorLeft*right*factorRight;  
            double errorsq =  factorLeft*factorRight*(right*right*leftError*leftError+left*left*rightError*rightError);
            double error = SVD_Sqrt(errorsq);
            
            
            // Save it
            (histoLeft+h)->SetBinContent(i, value);
            (histoLeft+h)->SetBinError(i, error); 
        } 
    } 
}

// Divide histograms
// This function multiplies the bin wise inverse of the right histo to the left histo.
// This means, the left histogram is CHANGED, the right one is not changed.
void TopSVDFunctions::SVD_MultInvRightToLeft(TH1D*& histoLeft, TH1D* histoRight, double factorLeft, double factorRight, double divZeroDefault,  int numHist)
{
     
    // Existence of Objects
    if ( histoLeft == NULL ) return;
    if ( histoRight == NULL ) return;
    
    
    // Create new objects 
    for ( int h = 0 ; h < numHist ; h++ ) { 
          
        int binsx = (histoLeft+h)->GetNbinsX();  
        for ( int i = 1 ; i <= binsx ; i++ ) {
            
            // Get Values
            double left = (histoLeft+h)->GetBinContent(i);
            double leftError = (histoLeft+h)->GetBinError(i);
            double right = (histoRight+h)->GetBinContent(i);
            double rightError = (histoRight+h)->GetBinError(i); 
            
            // Default Values
            double value = divZeroDefault;
            double error = 0.;
            
            // Divide carefully
            double denom = right*factorRight;
            double nom = left*factorLeft;
            if ( denom != 0.) { 
                value = nom / denom;
                double errorsq = ((factorLeft*factorLeft)/(denom*denom))*(leftError*leftError + ((left*left)/(right*right))*rightError*rightError);
                error = SVD_Sqrt(errorsq);
            }  
            
            // Save it
            (histoLeft+h)->SetBinContent(i, value);
            (histoLeft+h)->SetBinError(i, error); 
        } 
    } 
}

// Add up bins from histograms
// This function adds up the content of the right histo to the left histo.
// This means, the left histogram is CHANGED, the right one is not changed.
void TopSVDFunctions::SVD_AddBinRightToLeft(TH1D*& histoLeft, TH1D* histoRight, int bin, double factorLeft, double factorRight, int numHist)
{
     
    // Existence of Objects
    if ( histoLeft == NULL ) return;
    if ( histoRight == NULL ) return;
    
    
    // Create new objects 
    for ( int h = 0 ; h < numHist ; h++ ) {  

        // Get Values
        double left = (histoLeft+h)->GetBinContent(bin);
        double leftError = (histoLeft+h)->GetBinError(bin);
        double right = (histoRight+h)->GetBinContent(bin);
        double rightError = (histoRight+h)->GetBinError(bin);
        
        
        // Calculate new Values
        double value = left*factorLeft+right*factorRight;
        double errorsq = factorLeft*factorLeft*leftError*leftError+factorRight*factorRight*rightError*rightError;
        double error = SVD_Sqrt(errorsq);
        
        
        // Save it
        (histoLeft+h)->SetBinContent(bin, value);
        (histoLeft+h)->SetBinError(bin, error);   
    } 
}




// Get the diagonal elements from a matrix
TH1D* TopSVDFunctions::SVD_DiagElements(TH2D* histo, int numHist)
{ 
    
    // Existence of Objects
    if ( histo == NULL ) return NULL;   
     
    // Create new objects
    TH1D* outputhists = new TH1D[numHist];
    for ( int h = 0 ; h < numHist ; h++ ) { 
         
        // Bins
        bool isVariableBinSize = (histo+h)->GetXaxis()->IsVariableBinSize();
        int binsx = (histo+h)->GetNbinsX();
        if ( isVariableBinSize == true ) {
            const double* theBinsx = (histo+h)->GetXaxis()->GetXbins()->GetArray();
            (outputhists+h)->SetBins(binsx, theBinsx);
        } else { 
            double xmin = (histo+h)->GetXaxis()->GetXmin();
            double xmax = (histo+h)->GetXaxis()->GetXmax();
            (outputhists+h)->SetBins(binsx, xmin, xmax);
        } 
        
        // Loop over bins
        for ( int i = 1 ; i <= binsx ; i++ ) {
            
            // Get Values
            double value = histo->GetBinContent(i, i);
            double error = histo->GetBinError(i, i);
            
            // Save it
            (outputhists+h)->SetBinContent(i, value);
            (outputhists+h)->SetBinError(i, error); 
        }
            
             
        // Make sure the new histogram is not on a TDirectory
        // You will most likely run into garbage collector problems,
        (outputhists+h)->SetDirectory(NULL);
    }
    
    // Return
    return outputhists;
}


// Empties all bins and all errors in a histogram
// Can handle arrays of histograms
void TopSVDFunctions::SVD_EmptyHistogram1D(TH1D* histo, int numHist)
{
    // Existence of Objects
    if ( histo == NULL ) return;   
    
    
    int nbins = histo->GetNbinsX();
    
    // Loop over all histograms in the array 
    for ( int h = 0 ; h < numHist ; h++ ) {
        for ( int i = 0 ; i < nbins ; i++ ) {
            (histo+h)->SetBinContent(i, 0.);
            (histo+h)->SetBinError(i, 0.);
        }
    }
}



// Empties all bins and all errors in a histogram
// Can handle arrays of histograms
void TopSVDFunctions::SVD_EmptyHistogram2D(TH2D* histo, int numHist)
{ 
    // Existence of Objects
    if ( histo == NULL ) return;  
    
     
    int nbinsx = histo->GetNbinsX();
    int nbinsy = histo->GetNbinsY();
    int nbins = (nbinsx+2) * (nbinsy+2);
    
    for ( int h = 0 ; h < numHist ; h++ ) {
        for ( int i = 0 ; i < nbins ; i++ ) {
            (histo+h)->SetBinContent(i, 0.);
            (histo+h)->SetBinError(i, 0.);
        }
    }
}



// Empties all errors in a histogram
// Can handle arrays of histograms
void TopSVDFunctions::SVD_EmptyHistoErrors1D(TH1D* histo, int numHist)
{
    // Existence of Objects
    if ( histo == NULL ) return;   
    
    int nbins = histo->GetNbinsX();
    
    // Loop over all histograms in the array 
    for ( int h = 0 ; h < numHist ; h++ ) {
        for ( int i = 0 ; i < nbins ; i++ ) { 
            (histo+h)->SetBinError(i, 0.);
        }
    }
}

// Empties all errors in a histogram
// Can handle arrays of histograms
void TopSVDFunctions::SVD_EmptyHistoErrors2D(TH2D* histo, int numHist)
{ 
    // Existence of Objects
    if ( histo == NULL ) return;  
    
    
    int nbinsx = histo->GetNbinsX();
    int nbinsy = histo->GetNbinsY();
    int nbins = (nbinsx+2) * (nbinsy+2);
    
    // Loop over all histograms in the array 
    for ( int h = 0 ; h < numHist ; h++ ) {
        for ( int i = 0 ; i < nbins ; i++ ) { 
            (histo+h)->SetBinError(i, 0.);
        }
    }
}


// Creates a histogram with just the errors
// Attention:
// (1) Depicts relative uncertainties in PERCENT
// (2) TH1D-side bins are exluded
TH1D* TopSVDFunctions::SVD_HistoErrors1D(TH1D* histo, int numHist)
{
    // Existence of Objects
    if ( histo == NULL ) return NULL;       
    
    // Create new objects
    TH1D* outputhists = new TH1D[numHist];
    for ( int h = 0 ; h < numHist ; h++ ) {
         
        // Bins
        bool isVariableBinSize = (histo+h)->GetXaxis()->IsVariableBinSize();
        int binsx = (histo+h)->GetNbinsX();
        if ( isVariableBinSize == true ) {
            const double* theBinsx = (histo+h)->GetXaxis()->GetXbins()->GetArray();
            (outputhists+h)->SetBins(binsx, theBinsx);
        } else { 
            double xmin = (histo+h)->GetXaxis()->GetXmin();
            double xmax = (histo+h)->GetXaxis()->GetXmax();
            (outputhists+h)->SetBins(binsx, xmin, xmax);
        } 
        
        // Loop over bins
        for ( int i = 1 ; i <= binsx ; i++ ) {
             
            
            // Calculate Uncertainty
            double errorAbs = (histo+h)->GetBinError(i);
            double value = (histo+h)->GetBinContent(i);
            double error = 0.;
            if ( value > 0. ) {
                error = 100* errorAbs / value;
            } 
            
            // Save it
            (outputhists+h)->SetBinContent(i, error);
            (outputhists+h)->SetBinError(i, 0.); 
        }
            
             
        // Make sure the new histogram is not on a TDirectory
        // You will most likely run into garbage collector problems,
        (outputhists+h)->SetDirectory(NULL);
    }
    
    // Return array
    return outputhists;
}
 


// Creates a histogram with just the errors
// Attention:
// (1) Depicts relative uncertainties in PERCENT
// (2) TH2D-side bins are exluded
TH2D* TopSVDFunctions::SVD_HistoErrors2D(TH2D* histo, int numHist)
{
    
    
    // Existence of Objects
    if ( histo == NULL ) return NULL;   
    
    
    
    // Create new objects
    TH2D* outputhists = new TH2D[numHist];
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Bins
        bool isVariableBinSizeX = (histo+h)->GetXaxis()->IsVariableBinSize(); 
        bool isVariableBinSizeY = (histo+h)->GetXaxis()->IsVariableBinSize(); 
        int binsx = (histo+h)->GetNbinsX();
        int binsy = (histo+h)->GetNbinsY();
        if ( isVariableBinSizeX == true && isVariableBinSizeY == true ) { 
            const double* theBinsx = (histo+h)->GetXaxis()->GetXbins()->GetArray(); 
            const double* theBinsy = (histo+h)->GetYaxis()->GetXbins()->GetArray(); 
            (outputhists+h)->SetBins(binsx, theBinsx, binsy, theBinsy); 
        } else {  
            double xmin = (histo+h)->GetXaxis()->GetXmin();
            double xmax = (histo+h)->GetXaxis()->GetXmax();
            double ymin = (histo+h)->GetYaxis()->GetXmin();
            double ymax = (histo+h)->GetYaxis()->GetXmax(); 
            (outputhists+h)->SetBins(binsx, xmin, xmax, binsy, ymin, ymax);
        } 
        
        // Loop over bins
        for ( int i = 1 ; i <= binsx ; i++ ) {
            for ( int j = 1 ; j <= binsy ; j++ ) {
             
                
                // Calculate Uncertainty
                double errorAbs = (histo+h)->GetBinError(i, j);
                double value = (histo+h)->GetBinContent(i, j);
                double error = 0.;
                if ( value > 0. ) {
                    error = 100* errorAbs / value;
                } 
                
                // Save it
                (outputhists+h)->SetBinContent(i, j, error);
                (outputhists+h)->SetBinError(i, j, 0.);
        
                // Make sure the new histogram is not on a TDirectory
                // You will most likely run into garbage collector problems,
                (outputhists+h)->SetDirectory(NULL);
            }
        }
        
        // Make sure the new histogram is not on a TDirectory
        // You will most likely run into garbage collector problems,
        (outputhists+h)->SetDirectory(NULL);
    }
    
    // Return array
    return outputhists;
}


// Return new number of bins
int TopSVDFunctions::SVD_NumNewBins(int oldNumBins, bool addLowerSideBin, bool addUpperSideBin)
{ 
    int nbins = oldNumBins;
    if ( addLowerSideBin == true ) nbins++;
    if ( addUpperSideBin == true ) nbins++; 
    
    return nbins;
}


// New Binning 
void TopSVDFunctions::SVD_NewBins(double*& theNewBins, int newNumBins, const double* theOldBins, int oldNumBins, bool addLowerSideBin, bool addUpperSideBin)
{
    // Width of Side bins
    double rangeWidth = theOldBins[oldNumBins] - theOldBins[0];
    double ofWidth = rangeWidth / ((double) oldNumBins);
     
    // Potentially smaller UF Bin
    double ufWidth = ofWidth;
    if ( theOldBins[0] > 0. && theOldBins[0] < ufWidth ) ufWidth = theOldBins[0];

    // Calculate new bins 
    if ( addLowerSideBin == true) {
        theNewBins[0] = theOldBins[0] - ufWidth;
        for ( int i = 0 ; i <= oldNumBins ; i++ ) {
            theNewBins[i+1] = theOldBins[i]; 
        }
    } else {
        for ( int i = 0 ; i <= oldNumBins ; i++ ) {
            theNewBins[i] = theOldBins[i];
        } 
    } 
    if ( addUpperSideBin == true) { 
        theNewBins[newNumBins] = theOldBins[oldNumBins]+ofWidth;
    } 
    
}



// Rebinning of 1D Histograms. 
// ATTENTION! ATTENTION! ATTENTION! ATTENTION! ATTENTION!
// (1) Creates new histograms on the heap. 
//     So dont forget to delete them sometimes. 
// (2) In 'binsx' the bin boundaries are specified!
//     'nbinsx' gives the number of specified bins!
// (3) If wished, one can add lower and upper side bins
//     (This is highly recommended!) 
//     Note: If no side bins are added, all events  
//     outside the central bins will be discarded!!
// (5) Each bin boundary in 'binsx' MUST match one in   
//     the input-histograms, otherwise double counting 
//     will occur.
// (6) Can also handle Arrays of histograms, if numHist > 1.
//     All input histograms must have the same binning!
TH1D* TopSVDFunctions::SVD_Rebin1D(TH1D* input, int nbinsx, const double* binsx, bool addLowerSideBin, bool addUpperSideBin, int numHist) 
{ 
    
    bool debugThis = false;
    
    // Debug Output
    if ( debugThis == true ) {
        cout << endl;
        cout << "************************************************************" << endl;
        cout << "1D-Rebinning of " <<  input->GetName() << endl; 
    }
    
    
    addLowerSideBin = true;
    addUpperSideBin = true;
    
    // Existence of Objects
    if ( input == NULL ) return NULL;   
     
    // Binning of old histogram
    int nbinsGenRecX =input->GetNbinsX();
     
    // How many bins do you need
    int numnewbins = SVD_NumNewBins(nbinsx, addLowerSideBin, addUpperSideBin);
     
    // Calculate new bins
    double arrayOfNewBinsX[numnewbins+1];
    double* newbins = arrayOfNewBinsX;
    SVD_NewBins(newbins, numnewbins, binsx, nbinsx, addLowerSideBin, addUpperSideBin);    
    
    
    // Create new objects
    TH1D* outputhists = new TH1D[numHist];
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Debug Output
        if ( debugThis == true ) {
            cout << "    Histo Number " << h << endl;
        } 
        
        // Clean fresh histogram
        TString nameNewHist = (input+h)->GetName();
        nameNewHist.Append("Rebinned");
        (outputhists+h)->SetNameTitle(nameNewHist, nameNewHist);
        (outputhists+h)->SetBins(numnewbins, newbins);
        SVD_EmptyHistogram1D(outputhists+h, 1); 
    
    
    
        // Do NOT use the side bins of th TH1D-Class! 
        for(int i=1; i<=numnewbins; ++i){ 
        
            // Variables for the bin boundaries
            double xlow  = 0.;
            double xhigh = 0.;
            double newxlow = 0.;
            double newxhigh = 0.;
            double binwidth = 0.; 
            bool alarm = false;
            double alarmTresh = 0.0001;
            
            // Which treatment for this bin?
            bool isUFBin = false;
            bool isOFBin = false; 
            if ( i == 1 && addLowerSideBin == true ) isUFBin = true; 
            if ( i == numnewbins && addUpperSideBin == true ) isOFBin = true; 
                
            // Debug Output
            if ( debugThis == true ) {
                cout << "        Bin i/j =  " << i;
                if ( isUFBin == true ) cout << ", Underflow in X";
                if ( isOFBin == true ) cout << ", Overflow in X"; 
                cout << endl;
            }  
            
            
            // Find bin boundaries in X
            int nxlow = -1;
            int nxhigh = -1;
            if ( isUFBin == true) {    
                xlow = newbins[0]; 
                xhigh = newbins[1];  
                nxlow = 0;
                nxhigh = SVD_FindBin((input+h)->GetXaxis(), xhigh) - 1;   
                newxlow = xlow;
                newxhigh = (input+h)->GetXaxis()->GetBinLowEdge(nxhigh+1);
                binwidth = xhigh - xlow; 
                if ( TMath::Abs(newxhigh-xhigh)/binwidth >  alarmTresh ) alarm = true;
            }
            if ( isUFBin == false && isOFBin == false) {  
                xlow = newbins[i-1];
                xhigh = newbins[i];
                nxlow = SVD_FindBin((input+h)->GetXaxis(), xlow);
                nxhigh = SVD_FindBin((input+h)->GetXaxis(), xhigh) - 1;  
                newxlow = (input+h)->GetXaxis()->GetBinLowEdge(nxlow);
                newxhigh = (input+h)->GetXaxis()->GetBinLowEdge(nxhigh+1);
                binwidth = xhigh - xlow;
                if ( TMath::Abs(newxhigh-xhigh)/binwidth >  alarmTresh ) alarm = true;
                if ( TMath::Abs(newxlow-xlow)/binwidth >  alarmTresh ) alarm = true;
            }
            if ( isOFBin == true ) {   
                xlow = newbins[numnewbins-1]; 
                xhigh = newbins[numnewbins]; 
                nxlow = SVD_FindBin((input+h)->GetXaxis(), xlow);
                nxhigh = nbinsGenRecX+1;
                newxlow = (input+h)->GetXaxis()->GetBinLowEdge(nxlow);
                newxhigh = xhigh;
                binwidth = xhigh - xlow;
                if ( TMath::Abs(newxlow-xlow)/binwidth >  alarmTresh ) alarm = true;
            }
            
                
            // Debug Output
            if ( debugThis == true ) {
                cout << "            X / Range:     " << xlow << " ... " << xhigh << endl; 
                cout << "            X / Old Bins:  " << nxlow << "-" << nxhigh << endl;
                cout << "            X / Old Range: " << newxlow << " ... " << newxhigh << endl;
            }

    
            // Alarm, if the bin boundaries do not match
            if ( alarm == true ) {
                cout << "*****************************************************************" << endl;
                cout << "Error in TopSVDFunctions::SVD_Rebin1D" << endl;
                cout << "The new bin boundaries cannot be matched to the ones in the" << endl;
                cout << "input histogram." << endl;
                cout << "Old bins in X: " << endl;
                for ( int k = 1; k <= nbinsGenRecX+1 ; k++ ) {
                    cout << " " <<  (input+h)->GetXaxis()->GetBinLowEdge(k);
                } 
                cout << endl;
                cout << "New bins in X: " << endl;
                for ( int k = 0; k < numnewbins+1 ; k++ ) {
                    cout << " " << newbins[k];
                }
                cout << endl;
                cout << "Exiting ... " << endl;
                cout << "*****************************************************************" << endl;
                exit(1);
            }
                 
     
            // Integrate over all the original bins
            double sum = 0.;
            double sumErrSq = 0.;
            if ( nxlow <= nxhigh ) {
                for ( int b = nxlow ; b <= nxhigh ; b++ ) {
                    double binContent = ((TH1D*) (input+h))->GetBinContent(b);
                    double binError = ((TH1D*) (input+h))->GetBinError(b);
                    sum += binContent;
                    sumErrSq += binError*binError;
                }
            } else {
                sum = 0.;
                sumErrSq = 0.;
            }
                
                
            // Debug Output
            if ( debugThis == true ) {
                cout << "            Integral::     " << sum << " +/- " << TMath::Sqrt(sumErrSq) << endl;
            }
     

            // Save it
            (outputhists+h)->SetBinContent(i, sum);
            (outputhists+h)->SetBinError(i, TMath::Sqrt(sumErrSq));
        }
            
             
        // Make sure the new histogram is not on a TDirectory
        // You will most likely run into garbage collector problems,
        (outputhists+h)->SetDirectory(NULL);
        
    }

    
    
    if ( debugThis == true ) {
        cout << endl;
        cout << "Done. " << endl;
        cout << "************************************************************" << endl; 
        cout << endl;
    }
    
    // Return the array
    return outputhists;
}

    


// Rebinning of a 2D Histogram. 
// ATTENTION! ATTENTION! ATTENTION! ATTENTION! ATTENTION!
// (1) Creates a new histogram on the heap. 
//     So dont forget to delete it sometimes.
// (2) If needed, you can transpose the matrix, too. 
// (3) In 'binsx' and 'binsy' the bin boundaries are specified!
//     'nbinsx' and 'nbinsy' the number of specified bins!
// (4) If wished, one can add lower and upper side bins
//     (This is highly recommended!) 
//     Note: If no side bins are added, all events  
//     outside the central bins will be discarded!! 
// (5) Each (but the first and last) bin boundary in 'binsx' and 'binsy' 
//     MUST match one in the input-histogram, otherwise double counting 
//     will occur.
TH2D* TopSVDFunctions::SVD_Rebin2D(TH2D* input, int nbinsx, const double* binsx, int nbinsy, const double* binsy, bool addLowerSideBin, bool addUpperSideBin, int numHist, bool transpose)
{ 
    bool debugThis = false;
    
    // Debug Output
    if ( debugThis == true ) {
        cout << endl;
        cout << "************************************************************" << endl;
        cout << "2D-Rebinning of " <<  input->GetName() << endl; 
    }
    
    addLowerSideBin = true;
    addUpperSideBin = true;
    
    
    // Existence of Objects
    if ( input == NULL ) return NULL;  
    
    // Binning of old histogram
    int nbinsGenRecX =input->GetNbinsX();
    int nbinsGenRecY =input->GetNbinsY();
     
     
    // How many bins do you need
    int numnewbinsx = SVD_NumNewBins(nbinsx, addLowerSideBin, addUpperSideBin);
    int numnewbinsy = SVD_NumNewBins(nbinsy, addLowerSideBin, addUpperSideBin);
     
     
    // Calculate new bins
    double arrayOfNewBinsX[numnewbinsx+1];
    double* newbinsx = arrayOfNewBinsX;
    SVD_NewBins(newbinsx, numnewbinsx, binsx, nbinsx, addLowerSideBin, addUpperSideBin);   
    double arrayOfNewBinsY[numnewbinsy+1];
    double* newbinsy = arrayOfNewBinsY;
    SVD_NewBins(newbinsy, numnewbinsy, binsy, nbinsy, addLowerSideBin, addUpperSideBin); 
    
    
    
    // Create new objects
    TH2D* outputhists = new TH2D[numHist];
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Debug Output
        if ( debugThis == true ) {
            cout << "    Histo Number " << h << endl;
        } 
        
        // Clean fresh histogram
        TString nameNewHist = (input+h)->GetName();
        nameNewHist.Append("Rebinned");
        (outputhists+h)->SetNameTitle(nameNewHist, nameNewHist);
        (outputhists+h)->SetBins(numnewbinsx, newbinsx, numnewbinsy, newbinsy);
        SVD_EmptyHistogram2D(outputhists+h, 1);  
    
        // Do NOT use the side bins of th TH2D-Class! 
        for(int i=1; i<=numnewbinsx; ++i){ 
            for(int j=1; j<=numnewbinsy; ++j){ 
                
                // Variables for the bin boundaries
                double xlow  = 0.;
                double xhigh = 0.;
                double ylow  = 0.;
                double yhigh = 0.;
                double newxlow = 0.;
                double newxhigh = 0.;
                double newylow = 0.;
                double newyhigh = 0.;
                double binwidthx = 0.; 
                double binwidthy = 0.; 
                bool alarm = false;
                double alarmTresh = 0.0001;
            
            
                // Which treatment for this bin?
                bool isUFBinX = false;
                bool isOFBinX = false; 
                if ( i == 1 && addLowerSideBin == true ) isUFBinX = true; 
                if ( i == numnewbinsx && addUpperSideBin == true ) isOFBinX = true; 
                bool isUFBinY = false;
                bool isOFBinY = false; 
                if ( j == 1 && addLowerSideBin == true ) isUFBinY = true; 
                if ( j == numnewbinsy && addUpperSideBin == true ) isOFBinY = true; 
                
                // Debug Output
                if ( debugThis == true ) {
                    cout << "        Bin i/j =  " << i << "/" << j;
                    if ( isUFBinX == true ) cout << ", Underflow in X";
                    if ( isOFBinX == true ) cout << ", Overflow in X";
                    if ( isUFBinY == true ) cout << ", Underflow in Y";
                    if ( isOFBinY == true ) cout << ", Overflow in Y";
                    cout << endl;
                }  
                
                // Find bin boundaries in X
                int nxlow = -1;
                int nxhigh = -1;
                if ( isUFBinX == true ) {    
                    xlow = newbinsx[0]; 
                    xhigh = newbinsx[1];  
                    nxlow = 0;
                    nxhigh = SVD_FindBin((input+h)->GetXaxis(), xhigh) - 1;   
                    newxlow = xlow;
                    newxhigh = (input+h)->GetXaxis()->GetBinLowEdge(nxhigh+1);
                    binwidthx = xhigh - xlow; 
                    if ( TMath::Abs(newxhigh-xhigh)/binwidthx >  alarmTresh ) alarm = true; 
                }
                if ( isUFBinX == false && isOFBinX == false) {  
                    xlow = newbinsx[i-1];
                    xhigh = newbinsx[i];
                    nxlow = SVD_FindBin((input+h)->GetXaxis(), xlow);
                    nxhigh = SVD_FindBin((input+h)->GetXaxis(), xhigh) - 1;  
                    newxlow = (input+h)->GetXaxis()->GetBinLowEdge(nxlow);
                    newxhigh = (input+h)->GetXaxis()->GetBinLowEdge(nxhigh+1);
                    binwidthx = xhigh - xlow;
                    if ( TMath::Abs(newxhigh-xhigh)/binwidthx >  alarmTresh ) alarm = true;
                    if ( TMath::Abs(newxlow-xlow)/binwidthx >  alarmTresh ) alarm = true; 
                }
                if ( isOFBinX == true ) {   
                    xlow = newbinsx[numnewbinsx-1]; 
                    xhigh = newbinsx[numnewbinsx]; 
                    nxlow = SVD_FindBin((input+h)->GetXaxis(), xlow);
                    nxhigh = nbinsGenRecX+1;
                    newxlow = (input+h)->GetXaxis()->GetBinLowEdge(nxlow);
                    newxhigh = xhigh;
                    binwidthx = xhigh - xlow;
                    if ( TMath::Abs(newxlow-xlow)/binwidthx >  alarmTresh ) alarm = true; 
                }
                
                // Debug Output
                if ( debugThis == true ) {
                    cout << "            X / Range:     " << xlow << " ... " << xhigh << endl; 
                    cout << "            X / Old Bins:  " << nxlow << "-" << nxhigh << endl;
                    cout << "            X / Old Range: " << newxlow << " ... " << newxhigh << endl;
                }
                    
    
                // Find bin boundaries in Y
                int nylow = -1;
                int nyhigh = -1;
                if ( isUFBinY == true  ) {  // UF bin in Y
                    ylow = newbinsy[0]; 
                    yhigh = newbinsy[1];  
                    nylow = 0;
                    nyhigh = SVD_FindBin((input+h)->GetYaxis(), yhigh) - 1;   
                    newylow = ylow;
                    newyhigh = (input+h)->GetYaxis()->GetBinLowEdge(nyhigh+1);
                    binwidthy = yhigh - ylow; 
                    if ( TMath::Abs(newyhigh-yhigh)/binwidthy >  alarmTresh ) alarm = true;   
                }
                if ( isUFBinY == false && isOFBinY == false ) { // Center bin in Y
                    ylow = newbinsy[j-1];
                    yhigh = newbinsy[j];
                    nylow = SVD_FindBin((input+h)->GetYaxis(), ylow);
                    nyhigh = SVD_FindBin((input+h)->GetYaxis(), yhigh) - 1;  
                    newylow = (input+h)->GetYaxis()->GetBinLowEdge(nylow);
                    newyhigh = (input+h)->GetYaxis()->GetBinLowEdge(nyhigh+1);
                    binwidthy = yhigh - ylow;
                    if ( TMath::Abs(newyhigh-yhigh)/binwidthy >  alarmTresh ) alarm = true;
                    if ( TMath::Abs(newylow-ylow)/binwidthy >  alarmTresh ) alarm = true;  
                } 
                if ( isOFBinY == true) {   // UF bin in Y
                    ylow = newbinsy[numnewbinsy-1]; 
                    yhigh = newbinsy[numnewbinsy]; 
                    nylow = SVD_FindBin((input+h)->GetYaxis(), ylow);
                    nyhigh = nbinsGenRecY+1;
                    newylow = (input+h)->GetYaxis()->GetBinLowEdge(nylow);
                    newyhigh = yhigh;
                    binwidthy = yhigh - ylow;
                    if ( TMath::Abs(newylow-ylow)/binwidthy >  alarmTresh ) alarm = true;  
                }
                
                // Debug Output
                if ( debugThis == true ) {
                    cout << "            Y / Range:     " << ylow << " ... " << yhigh << endl; 
                    cout << "            Y / Old Bins:  " << nylow << "-" << nyhigh << endl;
                    cout << "            Y / Old Range: " << newylow << " ... " << newyhigh << endl;
                }
                
                
                
        
                // Alarm, if the bin boundaries do not match
                if ( alarm == true ) {
                    cout << "*****************************************************************" << endl;
                    cout << "Error in TopSVDFunctions::SVD_Rebin2D" << endl;
                    cout << "The new bin boundaries cannot be matched to the ones in the" << endl;
                    cout << "input histogram." << endl;
                    cout << "Old bins in X: " << endl;
                    for ( int k = 1; k <= nbinsGenRecX+1 ; k++ ) {
                        cout << " " <<  (input+h)->GetXaxis()->GetBinLowEdge(k);
                    } 
                    cout << endl;
                    cout << "New bins in X: " << endl;
                    for ( int k = 0; k < numnewbinsx+1 ; k++ ) {
                        cout << " " << newbinsx[k];
                    }
                    cout << endl;
                    cout << "Old bins in Y: " << endl;
                    for ( int k = 1; k <= nbinsGenRecY+1 ; k++ ) {
                        cout << " " <<  (input+h)->GetYaxis()->GetBinLowEdge(k);
                    } 
                    cout << endl;
                    cout << "New bins in Y: " << endl;
                    for ( int k = 0; k < numnewbinsy+1 ; k++ ) {
                        cout << " " << newbinsy[k];
                    }
                    cout << endl;
                    cout << "Exiting ... " << endl;
                    cout << "*****************************************************************" << endl;
                    exit(1);
                }
                
                
                // Integrate over all the original bins        
                double sum = 0.;
                double sumErrSq = 0.;
                if ( nxlow <= nxhigh && nylow <= nyhigh ) {
                    for ( int bx = nxlow ; bx <= nxhigh ; bx++ ) {
                        for ( int by = nylow ; by <= nyhigh ; by++ ) { 
                            double binContent = ((TH1D*) (input+h))->GetBinContent(bx, by);
                            double binError = ((TH1D*) (input+h))->GetBinError(bx, by);
                            sum += binContent;
                            sumErrSq += binError*binError;
                        }
                    }
                } else { 
                       sum = 0.;
                       sumErrSq = 0.;
                }
                
                
                // Debug Output
                if ( debugThis == true ) {
                    cout << "            Integral::     " << sum << " +/- " << TMath::Sqrt(sumErrSq) << endl;
                }
                
    
                // Save it, thereby transpose if needed
                if ( transpose == true ) {
                    (outputhists+h)->SetBinContent(j, i, sum);
                    (outputhists+h)->SetBinError(j, i, TMath::Sqrt(sumErrSq)); 
                } else { 
                    (outputhists+h)->SetBinContent(i, j, sum);
                    (outputhists+h)->SetBinError(i, j, TMath::Sqrt(sumErrSq)); 
                }
            }
        }
        
        // Make sure the new histogram is not on a TDirectory
        // You will most likely run into garbage collector problems,
        (outputhists+h)->SetDirectory(NULL);
    }
    
    
    if ( debugThis == true ) {
        cout << endl;
        cout << "Done. " << endl;
        cout << "************************************************************" << endl; 
        cout << endl;
    }
     
    // return it
    return outputhists;
                    
}




// Create Correlation Matrix from Covariance Matrix
// Creates a TH2D on the heap, so delete it eventually
// Attention: There are some cases, where one of the errors for the denominator
// is zero. This is okay, if the covariance is zero itself. In this case, the
// correlation should be zero. In all other cases, a warning message is thrown.   
TH2D* TopSVDFunctions::SVD_Cov2Corr(TH2D* input, TString name, TString axisName, TString type)
{
    int nbins = input->GetNbinsX(); // Same for X and Y

    // Clone input histogram
    TH2D* output = SVD_CloneHists2D(input);
    input->SetDirectory(NULL);
    output->SetName(name);
    output->SetTitle(name);
    TString corrStrZ(type + " Correlations in %");
    output->GetXaxis()->SetTitle(axisName);
    output->GetYaxis()->SetTitle(axisName);
    output->GetZaxis()->SetTitle(corrStrZ);
    for ( int i = 1; i <= nbins ; i++ ) {
        for ( int j = 1; j <= nbins ; j++ ) {
            double cov = input->GetBinContent(i,j);
            double errsq1 = input->GetBinContent(i,i);
            double errsq2 = input->GetBinContent(j,j);
            double entry = 0.;
            if ( i != j ) {
                if ( errsq1 > 0. && errsq2 > 0. ) {
                    entry = 100*cov/(TMath::Sqrt(errsq1*errsq2));
                } else {
                    if ( cov == 0. &&  errsq1 >= 0. && errsq2 >= 0. && (errsq1 == 0. || errsq2 == 0.)) {
                        entry = 0.;
                    } else {
                        cout << "*********************************************************" << endl;
                        cout << "WARNING! TopSVDFunctions::SVD_Cov2Corr()" << endl;
                        cout << "Couldn't calculate correlation for i=" << i << ", j=" << j << endl;
                        cout << "    in Histogram: " << name << endl;
                        cout << "Reason: Funny numbers.  " << endl;
                        cout << "    Covariance(" << i << "," <<  j << ") = " << cov << endl;
                        cout << "    ErrorSquared(" << i << ") = " << errsq1 << endl;
                        cout << "    ErrorSquared(" << j << ") = " << errsq2 << endl; 
                        cout << "Setting Correlation to Zero for now. Please fix this." << endl;
                        cout << "*********************************************************" << endl;
                    } 
               }
           }
           output->SetBinContent(i,j, entry);
       }
    }

    // Return pointer
    return output;

}




// Create Error Plot from Covariance Matrix
// Creates a TH1D on the heap, so delete it eventually
TH1D* TopSVDFunctions::SVD_Cov2Err(TH2D* input, TH1D* means, TString name, TString axisName, TString type)
{
    int nbins = input->GetNbinsX(); 

    // Clone input histogram
    TH1D* output = SVD_CloneHists1D(means); 
    output->SetName(name);
    output->SetTitle(name);
    TString corrStrY(type + " Error in %");
    output->GetXaxis()->SetTitle(axisName);
    output->GetYaxis()->SetTitle(corrStrY);
    for ( int i = 1 ; i <= nbins ; i++ ) {
        double variance = input->GetBinContent(i,i);
        double value = means->GetBinContent(i);
        double error = 0.;
        if ( value > 0. ) error = 100*TMath::Sqrt(variance) / value;
        output->SetBinContent(i, error);
        output->SetBinError(i, 0.);
    }

    // Return pointer
    return output;
}






// Formats a 2D histogram according to current style
// Corrects formatting in a manner specific to 2D histos  
void TopSVDFunctions::SVD_Draw2D(TH2D* hist, TString options, double rangeL, double rangeR, bool setNoExp)
{
    // Set the Style    
    hist->UseCurrentStyle();

    // Z axis title offset
    double zTitleOffset = 1.6;
    hist->GetZaxis()->SetTitleOffset(zTitleOffset);

    // Text Style depending on maximal value
    double max = hist->GetMaximum();
    double min = hist->GetMinimum();
    TString oldTextStule = gStyle->GetPaintTextFormat();
    if ( max > 1000. || min < 0.001) {
        hist->SetMarkerSize(0.7);
    }
    
    // If histogram is empty, set range manually
    int nbinsX = hist->GetNbinsX();
    int nbinsY = hist->GetNbinsY();
    bool isEmpty = true;
    for ( int i = 1 ; i <= nbinsX ; i++ ) {
        for ( int j = 1 ; j <= nbinsY ; j++ ) {
            double binContent = hist->GetBinContent(i,j);
            //double binError= hist->GetBinError(i,j);
            if ( binContent != 0. ) {
                isEmpty = false;
            }
            if ( isEmpty == false ) break; 
        }
        if ( isEmpty == false ) break; 
    }
    if ( isEmpty == true ) {
        hist->GetZaxis()->SetRangeUser(0., 1.);
    }

    // Draw it, thereby create pad primitives
    hist->Draw(options);
    
        // restrict range if desired
    if(rangeL!=-9999 && rangeR!=-9999) {
      hist->GetXaxis()->SetRangeUser(rangeL, rangeR);
      hist->GetYaxis()->SetRangeUser(rangeL, rangeR);
    }
    if(setNoExp) {
      hist->GetXaxis()->SetNoExponent(true);
      hist->GetYaxis()->SetNoExponent(true);
    }

    // Change the position of the palette
    gPad->Update();
    TPaletteAxis* palette = NULL;
    palette = (TPaletteAxis*) hist->GetListOfFunctions()->FindObject("palette");
    if ( palette != NULL ) {
        palette->SetY1NDC(gStyle->GetPadBottomMargin()+0.05);
    }
     
    
}


// Formats a 1D histogram according to current style
// Corrects formatting in a manner specific to 1D histos
// If color < 1, then color will be set autmatically 
void TopSVDFunctions::SVD_Draw1D(TH1D* hist, TString options,  int color)
{
    // Use current style
    hist->UseCurrentStyle();
    if (color > 0) hist->SetLineColor(color);
    if (color > 0) hist->SetMarkerColor(color);

    
    // If histogram is empty, set range manually
    int nbinsX = hist->GetNbinsX(); 
    bool isEmpty = true;
    for ( int i = 1 ; i <= nbinsX ; i++ ) { 
        double binContent = hist->GetBinContent(i);
        //double binError= hist->GetBinError(i);
        if ( binContent != 0. ) {
            isEmpty = false; 
        }
        if ( isEmpty == false ) break; 
    }
    if ( isEmpty == true ) {
        hist->GetYaxis()->SetRangeUser(0., 1.);
    }
    
    // Draw it
    hist->Draw(options);

    // Correct Y Range
    double ymin = hist->GetMinimum();
    double ymax = hist->GetMaximum();
    double newYRange = 1.3 * (ymax -ymin);
    hist->SetMaximum(ymin+newYRange);

}

// Get Range of a Graph
void TopSVDFunctions::SVD_RangeGraph(TGraph* graph, double& xmin, double& xmax, double& ymin, double& ymax)
{
    
    int nPoints = graph->GetN();
     
    // Get Range
    xmax = -FLT_MAX;
    ymax = -FLT_MAX;
    xmin = +FLT_MAX;
    ymin = +FLT_MAX;
    for ( int i = 0; i < nPoints ; i++ ) {
         double theX = 0.;
         double theY = 0.;
         graph->GetPoint(i, theX, theY);
         if ( theX > xmax ) xmax = theX;
         if ( theY > ymax ) ymax = theY;
         if ( theX < xmin ) xmin = theX;
         if ( theY < ymin ) ymin = theY;
    }
     
}

// Get Range of a Histogram
// Attention:
// (1)  Does not look at OF bins, if showOF set to false
// (2)  Excludes EMPTY bins from the range calculation
// (3)  Errors bands are considered in the range calculation
// (4)  Can also handle an array of histograms
void TopSVDFunctions::SVD_RangeHisto(TH1D* histo, double& ymin, double& ymax, bool showOF, int numHist)
{
    // Existence of Objects
    if ( histo == NULL ) return; 
    
    
    // Error handling
    bool withErrors = true;
    
    // Number of bins
    int nPoints = histo->GetNbinsX();
     
    // Get Range 
    ymax = -FLT_MAX;
    ymin = +FLT_MAX;
     
    // Loop over all histograms in the array
    for ( int h = 0 ; h < numHist ; h++ ) {
    
        // Loop over all bins except OF  bins
        // Exclude both, TH1D side bins and the
        // side bins intrinsic to your unfolding code
        for ( int i = 1; i <= nPoints ; i++ ) {
            if ( i == 1 && showOF == false) continue;
            if ( i == nPoints && showOF == false) continue;         
            
            double theY = (histo+h)->GetBinContent(i);
            if ( theY == 0. ) continue;
            
            double theYUp = theY;
            double theYDown = theY;
            if ( withErrors == true ) {
                theYUp += (histo+h)->GetBinError(i);
                theYDown -= (histo+h)->GetBinError(i);
            } 
            if ( theYUp > ymax ) ymax = theYUp;
            if ( theYDown < ymin ) ymin = theYDown; 
        }
    
    } 
}

// Get Range of a Histogram
// Can also handle an array of histograms
// (1)  Does not look at OF bins, if showOF set to false
// (2)  Excludes EMPTY bins from the range calculation
// (3)  Errors bands are considered in the range calculation
void TopSVDFunctions::SVD_RangeStack(THStack* stack, double& ymin, double& ymax, bool showOF)
{
    // Initialize Variables 
    ymax = -FLT_MAX;
    ymin = +FLT_MAX;
       
    
    // Loop over all histos 
    TH1D* histo = NULL; 
    TList* listHistos = stack->GetHists();
    TListIter iter(listHistos);
    while ( (histo = (TH1D*) iter.Next()) != NULL ) {
        
        // Read out a single histogram
        double currYmax = -FLT_MAX;
        double currYmin = +FLT_MAX;
        SVD_RangeHisto(histo, currYmin, currYmax, showOF, 1);
        
        // Adapt Range 
        if ( currYmax > ymax ) ymax = currYmax; 
        if ( currYmin < ymin ) ymin = currYmin;
        
    } 
}


// Makes a graph from two vectors and formats it
// Corrects formatting in a manner specific to 1D histos
// Attention:
// (1) If color < 1, then color will be set autmatically
// (2) First, an empty histogram is drawn and the pointer is
//     returned in bgrHisto.
// (3) If you already have a histogram, than provide it in 
//     'bgrHisto'  
void TopSVDFunctions::SVD_DrawGraphAutoRange(TGraph* graph, TH1D*& bgrHisto, TString options,  int color, double marginFact)
{ 
    // Draw Graph with automatically  
    SVD_DrawGraph(graph, bgrHisto, options, color, marginFact);
      
}


// Makes a graph from two vectors and formats it
// Corrects formatting in a manner specific to 1D histos
// Attention:
// (1) If color < 1, then color will be set autmatically
// (2) First, an empty histogram is drawn and the pointer is
//     returned in bgrHisto.
// (3) If you already have a histogram, than provide it in 
//     'bgrHisto'  
void TopSVDFunctions::SVD_DrawGraphRange(TGraph* graph, TH1D*& bgrHisto, double ymin, double ymax, TString options,  int color)
{
    
    // Draw Graph with automatically set ranges FIRST
    SVD_DrawGraph(graph, bgrHisto, options, color);
    
    
    // Now make the correction to the range
    bgrHisto->SetMaximum(ymax);
    bgrHisto->SetMinimum(ymin); 
     
}


// Makes a graph from two vectors and formats it
// Corrects formatting in a manner specific to 1D histos
// Attention:
// (1) If color < 1, then color will be set autmatically
// (2) First, an empty histogram is drawn and the pointer is
//     returned in bgrHisto.
// (3) If you already have a histogram, than provide it in 
//     'bgrHisto'  
void TopSVDFunctions::SVD_DrawGraphZero(TGraph* graph, TH1D*& bgrHisto, TString options,  int color)
{
    
    // Draw Graph with automatically set ranges FIRST
    SVD_DrawGraph(graph, bgrHisto, options, color);
     
    // Now make the correction to the range
    bgrHisto->SetMinimum(0.);  
    
}

// Makes a graph from two vectors and formats it
// Corrects formatting in a manner specific to 1D histos
// Attention:
// (1) If color < 1, then color will be set autmatically
// (2) First, an empty histogram is drawn and the pointer is
//     returned in bgrHisto.
// (3) If you already have a histogram, than provide it in 
//     'bgrHisto'  
void TopSVDFunctions::SVD_DrawGraph(TGraph* graph, TH1D*& bgrHisto, TString options,  int color, double marginFact)
{
     
    // Get Range of Graph in x automatically
    double autoxmax = 0.;
    double autoxmin = 0.;
    double autoymax = 0.;
    double autoymin = 0.;
    SVD_RangeGraph(graph, autoxmin, autoxmax, autoymin, autoymax);
    
    // Add some extra margin on both sides in x 
    double newXmin = autoxmin;
    double newXmax = autoxmax;
    SVD_NewRange(autoxmin, autoxmax, newXmin, newXmax, gPad->GetLogx()); 
    
    // Add some extra margin on both sides in y 
    double newYmin = autoymin;
    double newYmax = autoymax;
    SVD_NewRange(autoymin, autoymax, newYmin, newYmax, gPad->GetLogx(), marginFact); 
     
    
    
    // Make new Histo
    if ( bgrHisto == NULL) { 
        
        // Create the object
        bgrHisto = new TH1D("bgrHist", "bgrHist", 10, newXmin, newXmax);
        bgrHisto->UseCurrentStyle();
        
        // Draw it
        bgrHisto->Draw();
        
        // Correct Y Range
        bgrHisto->SetMinimum(newYmin);
        bgrHisto->SetMaximum(newYmax);
    }
    
    // Draw options
    TString drawOptions = options;
    drawOptions.Append(" SAME");
    
    // Use current style
    graph->UseCurrentStyle();
    if (color > 0) graph->SetLineColor(color);
    if (color > 0) graph->SetMarkerColor(color);

    // Draw it
    graph->Draw(drawOptions);
      

}
 

// A little shift on the position of a text string in X
double TopSVDFunctions::SVD_TextPosX(TH1D* hist, double xpos, bool log, int orientation)
{ 
    // Get Axis Range
    double min = hist->GetXaxis()->GetXmin();
    double max = hist->GetXaxis()->GetXmax(); 
    if ( min >= max ) return 0.;
     
     
    // Position
    double thePos = xpos; 
    if ( thePos < min ) thePos = min;
    if ( thePos > max ) thePos = max; 
    
    
    // Log
    if ( log == true ) {
        if ( min > 0. ) { 
            min = TMath::Log10(min);     

        } else {
            return 0.;
        }
        if ( max > 0. ) { 
            max = TMath::Log10(max); 
        } else {
            return 0.;
        }
        if ( thePos > 0. ) {
            thePos = TMath::Log10(thePos); 
        } else {
            return 0.;
        }
    } 
    
    // Orientation ... take first digit
    int theOrientation = orientation / 10;
    
    
    // Shift
    double shift = 0.05*(max - min); 
    if ( theOrientation == 3 ) thePos = thePos - shift; 
    if ( theOrientation == 2 ) thePos = thePos; 
    if ( theOrientation == 1 ) thePos = thePos + shift; 
    
     
    // Log
    if ( log == true ) {
        thePos = TMath::Power(10., thePos);
    }
     
    // Return
    return thePos; 
}
 

// A little shift on the position of a text string in Y
double TopSVDFunctions::SVD_TextPosY(TH1D* hist, double ypos, bool log, int orientation)
{ 
    // Get Axis Range
    double min = hist->GetMinimum();
    double max = hist->GetMaximum(); 
    if ( min >= max ) return 0.;
     
      
    // Position
    double thePos = ypos; 
    if ( thePos < min ) thePos = min;
    if ( thePos > max ) thePos = max; 
    
    // Log
    if ( log == true ) {
        if ( min > 0. ) {
            min = TMath::Log10(min);
        } else {
            return 0.;
        }
        if ( max > 0. ) {
            max = TMath::Log10(max);
        } else {
            return 0.;
        }
        if ( thePos > 0. ) {
            thePos = TMath::Log10(thePos);
        } else {
            return 0.;
        }
    } 
    
    // Orientation ... take second digit
    int theOrientation = orientation % 10;
    
    
    // Shift
    double shift = 0.05*(max - min);
    if ( theOrientation == 3 ) thePos = thePos - shift; 
    if ( theOrientation == 2 ) thePos = thePos; 
    if ( theOrientation == 1 ) thePos = thePos + shift; 
    
    
    // Log
    if ( log == true ) {
        thePos = TMath::Power(10., thePos);
    }
    
    // Return
    return thePos; 
    
} 
 

// Draws a String somewhere in your pad
// If a TLatex object is given, then it will be used. Otherwise, it will be created.
TLatex* TopSVDFunctions::SVD_DrawText(TString theText, double xpos, double ypos, int orientation, TLatex*& mylatex)
{
    // Create Latex Object
    if ( mylatex == NULL ) {
        mylatex = new  TLatex();
        mylatex->SetTextSize(0.05);  
    }
     
    // Set Orientation
    mylatex->SetTextAlign(orientation);
    
    
    // Draw it
    mylatex->DrawLatex(xpos, ypos, theText);
    
    
    // Make TPaveText
    //TPaveText* pt = new TPaveText();
    
    
    // Return it;
    return mylatex;
        
}


// Draws a vertical line from top to bottom of bgrHisto.
// If color < 1, then color will be set autmatically 
TLine* TopSVDFunctions::SVD_DrawVertLine(TH1D*& bgrHisto, double xpos, int color)
{
    // Get Range of Histogram
    double Ymin = bgrHisto->GetMinimum();
    double Ymax = bgrHisto->GetMaximum();
    
    // New Line
    TLine* line = new TLine(xpos, Ymin, xpos, Ymax);
    line->SetLineColor(color);
    line->SetLineWidth(1);
    line->SetLineStyle(1);
    
    // Test it, whether it is in the range of the histogram
    bool draw = true;
    int nbins = bgrHisto->GetNbinsX();
    double xmin = bgrHisto->GetXaxis()->GetBinLowEdge(1);
    double xmax = bgrHisto->GetXaxis()->GetBinLowEdge(nbins+1);
    if ( xpos < xmin ) draw = false;
    if ( xpos > xmax ) draw = false;
    
    // Draw it
    if ( draw == true ) {
        line->Draw("SAME");
    }
    
    // REturn it
    return line;
    
}

// Draws a horizontal line from left to right of bgrHisto.
// If color < 1, then color will be set autmatically 
TLine* TopSVDFunctions::SVD_DrawHorizLine(TH1D*& bgrHisto, double ypos, int color)
{
    // Get Range of Histogram
    int nbins = bgrHisto->GetNbinsX();
    double Xmin = bgrHisto->GetXaxis()->GetBinLowEdge(1);
    double Xmax = bgrHisto->GetXaxis()->GetBinLowEdge(nbins+1);
    
    // New Line
    TLine* line = new TLine(Xmin,  ypos, Xmax, ypos);
    line->SetLineColor(color);
    line->SetLineWidth(1);
    line->SetLineStyle(1);
    
    // Test it, whether it is in the range of the histogram
    bool draw = true; 
    double ymin = bgrHisto->GetMinimum();
    double ymax = bgrHisto->GetMaximum();
    if ( ypos < ymin ) draw = false;
    if ( ypos > ymax ) draw = false;
    
    // Draw it
    if ( draw == true ) {
        line->Draw("SAME");
    }
    
    // REturn it
    return line;
    
}


// Draws a vertical lines from top to bottom of bgrHisto. 
// If color < 1, then color will be set autmatically 
void TopSVDFunctions::SVD_DrawVertLines(TH1D*& bgrHisto, TVectorD* xpositions, int color)
{
    
    // Make Array of Pointers
    int nPoints = xpositions->GetNoElements();   
    
    
    // Loop
    for ( int i = 0 ; i < nPoints ; i++ ) {
        
        // Get Position
        double xpos = (*xpositions)[i];
        
        // Draw Line
        SVD_DrawVertLine(bgrHisto, xpos, color);
          
    }  
    
} 


// Makes a graph from two vectors and formats it 
// Attention: Creates a Graph on the heap. Dont forget to delete it!
TGraph* TopSVDFunctions::SVD_Vect2Graph(TVectorD* vectorX, TVectorD* vectorY)
{
    //Test, if vectorX and vectorY have same number of entries
    int nelements = vectorX->GetNoElements();
    int nelementsY = vectorY->GetNoElements();
    if ( nelements <= 0 || nelementsY <= 0 || nelements != nelementsY ) {
        cout << "Error in TopSVDFunctions::SVD_Vect2Graph()" << endl;
        cout << "   The two vectors have strange sizes." << endl;
        cout << "   vectorX has " << nelements  << " entries." << endl;
        cout << "   vectorY has " << nelementsY << " entries." << endl;
        cout << "Exiting ... " << endl;
        exit(1);
    }
    
    
    // New Graph
    TGraph* graph = new TGraph(*vectorX, *vectorY);
    
    // Return
    return graph;
    
}


// Makes a graph from a the coordinates of a point
// Attention: Creates a Graph on the heap. Dont forget to delete it!
TGraph* TopSVDFunctions::SVD_Point2Graph(double x, double y)
{
    // Create new vectors
    TVectorD vectorX(1);
    TVectorD vectorY(1);
    vectorX[0] = x;
    vectorY[0] = y;
    
    // New Graph
    TGraph* graph = new TGraph(vectorX, vectorY);
    
    // Return
    return graph;
    
}


// SVD Unfolding Helper FUnction
// Format a Ratio plot according to current style
void TopSVDFunctions::SVD_DrawRange(TH1D* histo, double ymin, double ymax, TString options,  int color)
{
    // Clone Histo
    TH1D* hist = SVD_CloneHists1D(histo);
    hist->SetDirectory(NULL);
 
    // Hide values you dont want
    for ( int i = 0; i < hist->GetNbinsX()+2; i++) {
         double val = hist->GetBinContent(i);
         if ( val > ymax || val < ymin ) {
             hist->SetBinContent(i, 2*ymin - ymax);
             hist->SetBinError(i, 0.);
         }
    }

 
    // Use current style
    hist->UseCurrentStyle();
    hist->SetLineColor(histo->GetLineColor());
    hist->SetMarkerColor(histo->GetMarkerColor());
    if (color > 0) hist->SetLineColor(color);
    if (color > 0) hist->SetMarkerColor(color);
    
    // Draw it
    hist->Draw(options);

    // Range
    double yMin = ymin;
    double yMax = ymax;
    double newYRange = yMax - yMin;
    hist->SetMaximum(ymin+newYRange);
    hist->SetMinimum(ymin);
   
}



// SVD Unfolding Helper FUnction
// Creates a Ratio plot of two histograms
// Errors: relative errors from the NUMERATOR are taken
// Boundaries: If ymax < ymin, then range will be set automatically
TH1D* TopSVDFunctions::SVD_MakeRatioNum(TH1D* numerator, TH1D* denominator, int numHist)
{
     
    // Existence of Objects
    if ( numerator == NULL ) return NULL;
    if ( denominator == NULL ) return NULL;  
      

    // Create new Histograms
    TH1D* hist = SVD_CloneHists1D(numerator, numHist);

    
    // Loop over all histograms in the array 
    for ( int h = 0 ; h < numHist ; h++ ) {

        // Fill Histogram
        int numbins = hist->GetNbinsX();
        for ( int i = 0 ; i < numbins ; i++) {
            double num = (numerator+h)->GetBinContent(i);
            double denom = (denominator+h)->GetBinContent(i);
            double abs = (numerator+h)->GetBinError(i);
            double val = 0.;
            double err = 0.;
            if ( denom != 0. && num != 0.) {
                val = num / denom;
                err = abs / num;
            }
            (hist+h)->SetBinContent(i, val);
            (hist+h)->SetBinError(i, err);
        }
    }
    
    // Return
    return hist;

}        


// SVD Unfolding Helper FUnction
// Creates a Ratio plot of two histograms
// Errors: relative errors from the DENOMINATOR are taken
// Boundaries: If ymax < ymin, then range will be set automaticallyf
TH1D* TopSVDFunctions::SVD_MakeRatioDenom(TH1D* numerator, TH1D* denominator, int numHist)
{
  
    // Existence of Objects
    if ( numerator == NULL ) return NULL;
    if ( denominator == NULL ) return NULL; 
     
     
    // Create new Histograms
    TH1D* hist = SVD_CloneHists1D(numerator, numHist);
 
    // Loop over all histograms in the array 
    for ( int h = 0 ; h < numHist ; h++ ) { 
 
        // Fill Histogram
        int numbins = hist->GetNbinsX();
        for ( int i = 0 ; i < numbins ; i++) {
            double num = (numerator+h)->GetBinContent(i);
            double denom = (denominator+h)->GetBinContent(i);
            double abs = (denominator+h)->GetBinError(i);
            double val = 0.;
            double err = 0.;
            if ( denom != 0. ) {
                val = num / denom;
                err = abs / denom;
            }
            (hist+h)->SetBinContent(i, val);
            (hist+h)->SetBinError(i, err);
        }
 
    }
    // Return
    return hist;

}
             

// SVD Unfolding Helper Function
// Creates a Ratio plot of two histograms
// Errors of ratio plot are all set to zero
// Boundaries: If ymax < ymin, then range will be set automatically
TH1D* TopSVDFunctions::SVD_MakeRatioZero(TH1D* numerator, TH1D* denominator, int numHist)
{
  
    // Existence of Objects
    if ( numerator == NULL ) return NULL;
    if ( denominator == NULL ) return NULL; 
    
    // Create new Histograms
    TH1D* hist = SVD_CloneHists1D(numerator, numHist);
 
    // Loop over all histograms in the array 
    for ( int h = 0 ; h < numHist ; h++ ) { 
   
        // Fill Histogram
        int numbins = hist->GetNbinsX();
        for ( int i = 0 ; i < numbins ; i++) {
            double num = (numerator+h)->GetBinContent(i);
            double denom = (denominator+h)->GetBinContent(i);
            double val = 0.;
            double err = 0.;
            if ( denom != 0. && num != 0.) {
                val = num / denom; 
            }
            (hist+h)->SetBinContent(i, val);
            (hist+h)->SetBinError(i, err);
        }
    }
    
    return hist;

}  



// Formats a Stack of Histos to the current style
// The last argument takes integers between 1 and 9 for the colors, i.e.
// cols=0   means: Do not touch the colors!
// cols=243 means: First Histogram red, second blue, and third and all following ones green 
// cols<0   means: Standard palette.
void TopSVDFunctions::SVD_DrawStackAutoRange(THStack* stack, TLegend* leg, TString xTitle, TString yTitle, TString options, int col, bool showOF, bool log, double rangeL, double rangeR, bool setNoExp)
{
    // Get Range
    double min = 0;
    double max = 0;
    SVD_RangeStack(stack, min, max, showOF);
    
    // Add margin
    double newmin = 0;
    double newmax = 0;
    SVD_NewRange(min, max, newmin, newmax, log);
    
    // Draw with new Range
    SVD_DrawStackRange(stack, leg, xTitle, yTitle, newmin, newmax, options, col, showOF, rangeL, rangeR, setNoExp);
}

// Formats a Stack of Histos to the current style
// The last argument takes integers between 1 and 9 for the colors, i.e.
// cols=0   means: Do not touch the colors!
// cols=243 means: First Histogram red, second blue, and third and all following ones green 
// cols<0   means: Standard palette.
void TopSVDFunctions::SVD_DrawStackZero(THStack* stack, TLegend* leg, TString xTitle, TString yTitle, TString options, int col, bool showOF, double rangeL, double rangeR, bool setNoExp)
{
    
    // Get Range of stack
    double min = FLT_MAX;
    double max = -FLT_MIN;
    SVD_RangeStack(stack, min, max, showOF);
    
    // Add UPPER margin
    // 'newmin' must be 0. 
    double newmax = -FLT_MIN;
    double newmin = 0;
    SVD_NewUpperRange(max, newmax);
    
    // Draw with new Range
    SVD_DrawStackRange(stack, leg, xTitle, yTitle, newmin, newmax, options, col, showOF, rangeL, rangeR, setNoExp);
}


// Formats a Stack of Histos to the current style
// The last argument takes integers between 1 and 9 for the colors, i.e.
// cols=0   means: Do not touch the colors!
// cols=243 means: First Histogram red, second blue, and third and all following ones green 
// cols<0   means: Standard palette.
// If ymax < ymin, the range will be calculated automatically  
void TopSVDFunctions::SVD_DrawStackRange(THStack* stack, TLegend* leg, TString xTitle, TString yTitle, double ymin, double ymax, TString options, int col, bool showOF, double rangeL, double rangeR, bool setNoExp)
{ 
    // Colors 
    bool useStandardPalette = false;
    bool useUserPalette = false;
    bool useOldColors = false;
    int* theColors = NULL; 
    int numdigits = 1;
    if ( col > 0 ) { 
        useUserPalette = true;
         
        // Create an array of colors 
        int colstmp1 = col;
        while ( (colstmp1=colstmp1/10)>0 ) { numdigits++; }
        theColors = new int[numdigits];
        
         
        // Save colors in that array
        int colstmp2 = col;
        for ( int i = 0 ; i < numdigits ; i++ ) { 
            int position = numdigits-1-i;
            int lastdigit = colstmp2 % 10;
            theColors[position] = lastdigit;
            colstmp2=colstmp2/10; 
        }
    } else if ( col < 0 ) { 
        useStandardPalette = true;
    } else {
        useOldColors = true;
    }
  
 
    // Range Issues
    bool doSetRange = false;
    if ( ymax > ymin ) doSetRange = true;
  
  
    // Loop over all histos
    // Set Colors and make cosmetic changes 
    TH1D* histo = NULL;
    int cnt = 0;
    TList* listHistos = stack->GetHists();
    TListIter iter(listHistos);
    while ( (histo = (TH1D*) iter.Next()) != NULL ) {
         
        
        // Hide Unwanted Values
        // These can be values outside the range or the OF bins 
        if ( doSetRange == true ) {    
            for ( int i = 0; i < histo->GetNbinsX()+2; i++) {
                 bool doHide = false;
                 double val = histo->GetBinContent(i);
                 if ( val > ymax || val < ymin ) doHide = true;
                 if ( showOF == false ) {
                     if ( i <=1 ) doHide = true;
                     if ( i >= histo->GetNbinsX() ) doHide = true;
                 }
                 if ( doHide == true ) {
                    histo->SetBinContent(i, 2*ymin - ymax);
                    histo->SetBinError(i, 0.);
                }
            }
        }
        
        // Set Color
        int currentcolor = 0;
        if ( useStandardPalette == true ) currentcolor = cnt+1;
        if ( useUserPalette == true ) currentcolor = theColors[cnt];
        if ( useOldColors == true ) currentcolor = histo->GetMarkerColor(); 
        double saveMarkerSize = histo->GetMarkerSize();
        double saveLineStyle = histo->GetLineStyle();
        histo->UseCurrentStyle();  
        histo->SetLineColor(currentcolor);
        histo->SetLineStyle(saveLineStyle);
        histo->SetMarkerColor(currentcolor);
        histo->SetMarkerSize(saveMarkerSize);  
        
        cnt++;
    } 
  

    // Draw it
    options.Append("nostack");
    stack->Draw(options);
    
    // restrict range if desired
    if(rangeL!=-9999 && rangeR!=-9999) stack->GetHistogram()->GetXaxis()->SetRangeUser(rangeL, rangeR);
    if(setNoExp) stack->GetHistogram()->GetXaxis()->SetNoExponent(true);

    // Style
    stack->GetHistogram()->UseCurrentStyle();

    // Titles
    stack->GetHistogram()->GetXaxis()->SetTitle(xTitle);
    stack->GetHistogram()->GetYaxis()->SetTitle(yTitle);

    // Y Range
    if ( doSetRange == false ) {
        double yMin = stack->GetHistogram()->GetMinimum();
        double yMax = stack->GetHistogram()->GetMaximum();
        double yRange = yMax-yMin;
        double yUpperMargin = 0.3*yRange;
        double yLowerMargin = 0.;
        if ( yMin < 0. ) {
            yLowerMargin = 0.05*yRange;
        } 
        stack->SetMaximum(yMax+yUpperMargin);
        stack->SetMinimum(yMin-yLowerMargin);
    } else { 
        double yMin = ymin;
        double yMax = ymax;
        //double yRange = yMax - yMin;
        double yUpperMargin = 0.;
        double yLowerMargin = 0.; 
        stack->SetMaximum(yMax+yUpperMargin);
        stack->SetMinimum(yMin-yLowerMargin); 
    } 
    
    // Legend
    if ( leg != NULL ) {
        leg->Draw("SAME");
    }

}



// Set Name, Title and Axis Titles
// note: The Title is set to the value of name
void TopSVDFunctions::SVD_SetTitles1D(TH1D* histo, TString name, TString xTitle, TString yTitle, int numHist)
{
    // Existence of Objects
    if ( histo == NULL ) return;    
    
    // Loop over all histograms in the array
    for ( int h = 0 ; h < numHist ; h++ ) {  
        (histo+h)->SetName(name);
        (histo+h)->SetTitle(name); 
        (histo+h)->GetXaxis()->SetTitle(xTitle);
        (histo+h)->GetYaxis()->SetTitle(yTitle);
    }
}




// Set Name, Title and Axis Titles
// note: The Title is set to the value of name
void TopSVDFunctions::SVD_SetTitles2D(TH2D* histo, TString name, TString xTitle, TString yTitle, TString zTitle, int numHist)
{
    // Existence of Objects
    if ( histo == NULL ) return;       
    
    // Loop over all histograms in the array
    for ( int h = 0 ; h < numHist ; h++ ) { 
        (histo+h)->SetName(name);
        (histo+h)->SetTitle(name);
        (histo+h)->GetXaxis()->SetTitle(xTitle);
        (histo+h)->GetYaxis()->SetTitle(yTitle);
        (histo+h)->GetZaxis()->SetTitle(zTitle);
    }
}
 


// Calculate ChiSquare between Data and Refolded Distribution
// Note: The side bins on rec level will NOT be considered!
// Empty bins will also NOT be considered!
double TopSVDFunctions::SVD_Chi2A(TH1D* dataHist, TH1D* refoldHist)
{
    double chi2a = 0.;
    for ( int i = 1 ; i <= dataHist->GetNbinsX(); i++) {
        double dataValue = dataHist->GetBinContent(i);
        double refoldedValue = refoldHist->GetBinContent(i);
        double dataErrorSquared = dataHist->GetBinContent(i);
        if ( dataValue <= 0. ) continue;
        double current = (dataValue - refoldedValue) * (dataValue - refoldedValue) / dataErrorSquared;
        chi2a += current;
    }
    return chi2a;
}

 


// Calculate the non-empty and non-side rec level bins
int TopSVDFunctions::SVD_GoodRecBins(TH1D* dataHist)
{
    int numbins = 0;
    for ( int i = 1 ; i <= dataHist->GetNbinsX(); i++) {
        double dataValue = dataHist->GetBinContent(i);
        if ( dataValue <= 0. ) continue;
        numbins++;
    }
    
    return numbins;
}




// Calculate the non-empty and gen level bins
int TopSVDFunctions::SVD_GoodGenBins(TH1D* xiniHist)
{
    int numbins = 0;
    for ( int i = 1 ; i <= xiniHist->GetNbinsX(); i++) {
        double xiniValue = xiniHist->GetBinContent(i);
        if ( xiniValue <= 0. ) continue;
        numbins++;
    }

    return numbins;
}



// Calculate NDOF
int TopSVDFunctions::SVD_Ndof(TH1D* dataHist, TH1D* xiniHist)
{
    int ndof = SVD_GoodRecBins(dataHist) - SVD_GoodGenBins(xiniHist);
    return ndof;
}



// Tau Scan
// Calculate Scan Points
TVectorD* TopSVDFunctions::SVD_CalcScanPoints(double firstTau, double lastTau, int nScanPoints)
{
    
    // Use 3 scan points minimum
    if ( nScanPoints < 3 ) nScanPoints = 3;
    
    // Setup Vector 
    TVectorD* vScanPoints = new TVectorD(nScanPoints);
    
    // Find the scan points
    // Use equidistant steps on a logarithmic scale
    double steplog = (TMath::Log10(lastTau) - TMath::Log10(firstTau)) / ((double) (nScanPoints-1));
    for ( int i = 0 ; i < nScanPoints ; i++ ) {
        double posScanPoint = TMath::Power(10., (TMath::Log10(firstTau) + ((double) i) * steplog));
        (*vScanPoints)[i] = posScanPoint;
    }
    
    return vScanPoints;
}

// Tau Scan
// Scan Global Correlation    
double TopSVDFunctions::SVD_ScanGlobalCorrelation(TH2D* statCovHist, TH1D* dataHist)
{  
  
    // Calculate Global Correlation in all bins
    TH1D* globCorrHist = SVD_CalcGlobCorr(statCovHist, dataHist);
        
    

    // Averaged Global Correlation
    double globcorr = SVD_AvgGlobCorr(globCorrHist); 
 

    // Delete old objects 
    //delete globCorrHist; 
    
    // return value 
    return globcorr;
    
}

// Tau Scan
// Scan Chi Squared
double TopSVDFunctions::SVD_ScanChiSquared(TH1D* weightHist, TH1D* dataHist, TH2D* mcHist)
{  
    // Refold it 
    TH1D* refoldHist = SVD_Refold(weightHist, NULL, mcHist, true);  
     
    // Calculate Chi2A
    double chisq = SVD_Chi2A(dataHist, refoldHist);
     
    // Delete old objects
    // delete refoldHist; 
      
    // Return
    return chisq;
    
}


// Tau Scan
// Scan Curvature
double TopSVDFunctions::SVD_ScanCurvature(double curv)
{  
 
    // return vector
    return curv; 
}

// This is a crude thing.
// It calculates the statistical uncertainty of the
// unfolding result, averaged over all bins except
// empty bins and side bins 
double TopSVDFunctions::SVD_ScanAvgSqErr(TH2D* covMatrix)
{
    // nbins
    int nbins = covMatrix->GetNbinsX(); 
    
    // Loop
    double sumSquaredErrors = 0.;
    int numRealBins = 0;
    
    // Recall: You have two sets of OF bins!
    for ( int i = 2 ; i < nbins ; i++ ) {
        
        // Get trace element of covariance Matrix
        double binContent = covMatrix->GetBinContent(i,i);
        
        // Do not count this bin if it is empty
        if ( binContent == 0. ) continue;
        
        sumSquaredErrors += binContent;
        numRealBins++;
    }
    
    // Average
    double avgSqErr = sumSquaredErrors / ((double) numRealBins );
    
    // Take Square Root
    avgSqErr = TMath::Sqrt(avgSqErr);
    
    // Return
    return avgSqErr;
} 


// This is a crude thing.
// It calculates the mean value of the
// unfolding result, averaged over all bins except
// empty bins and side bins
double TopSVDFunctions::SVD_ScanAvgMean(TH1D* unfHist)
{
    
    // nbins
    int nbins = unfHist->GetNbinsX(); 
    
    // Loop
    double sumElements = 0.; 
    int numRealBins = 0;
    
    // Recall: You have two sets of OF bins!
    for ( int i = 2 ; i < nbins ; i++ ) {
        
        // Get trace element of covariance Matrix
        double binContent = unfHist->GetBinContent(i);
        
        // Do not count this bin if it is empty
        if ( binContent == 0. ) continue;
        
        sumElements += binContent;
        numRealBins++;
    }
    
    // Average
    double avgMean =  sumElements / ((double) numRealBins );
    return avgMean; 
    
}

// Helper Function
// If you have a range you want to depict in a histogram,
// than this function gives you a good range 
void TopSVDFunctions::SVD_NewRange(double min, double max, double& newmin, double& newmax, bool log, double marginFact)
{
    // Set Min and Max Values
    double theMin = min;
    double theMax = max;
    if ( log == true ) {
        if ( min > 0 ) theMin = TMath::Log10(min);
        if ( max > 0 ) theMax = TMath::Log10(max);
    }
    
    // Get margin
    double margin = marginFact*(theMax - theMin);
    
    // New Min and Max
    double theNewMin = theMin - margin;
    double theNewMax = theMax + margin;
    
    
    // Return
    newmin = theNewMin;
    newmax = theNewMax;
    if ( log == true ) {
        newmin = TMath::Power(10., theNewMin);
        newmax = TMath::Power(10., theNewMax);
    } 
}
    
// Helper Function
// If you have a range you want to depict in a histogram,
// than this function gives you a good range 
// Here, Only the upper range will be changed
void TopSVDFunctions::SVD_NewUpperRange(double max, double& newmax)
{
    // Set Min and Max Values
    double theMin = 0.;
    double theMax = max; 
    
    // Get margin
    double margin = 0.3*(theMax - theMin);
    
    // New Min and Max 
    double theNewMax = theMax + margin;
    
     
    // Return 
    newmax = theNewMax;  
    
} 
    

// GLOBAL CORRELATION
// Attention: 
// (1) We do not use the side bins for the computation
//     of the global correlation !!!
// (2) We also exclude all empty rows/columns int he
//     covariance matrix.
// (3) We also exclude bins that are empty on REC level.
//     These are only determined by the periphery.
// (2) For excluded bins, the global correlation
//     will be set to zero.
TH1D* TopSVDFunctions::SVD_CalcGlobCorr(TH2D* statCovHist, TH1D* dataHist)
{
    // Number of bins
    int nbins = statCovHist->GetNbinsX();
     
    // Map Histo To Matrix
    TArrayI binMap(nbins);
    int numBinsToSkip = 0;
    
    // Search for bins that have to be skipped
    int bincounter = 0;
    for ( int i = 1; i <= nbins ; i++ ) {
        double data = dataHist->GetBinContent(i);
        if ( data <= 0. ) {
            // Through out bin with no data
            binMap[i-1] = -1;
            numBinsToSkip++;
        } else if ( i == 1 ) {
            // Through out underflow bin
            binMap[i-1] = -1;
            numBinsToSkip++;
        } else if ( i == nbins ) {
            // Through out overflow bin
            binMap[i-1] = -1;
            numBinsToSkip++;
        } else {
            // Search for  bins with empty rows/columns
            bool skipThisBin = true;
            for ( int j = 2; j <= nbins-1 ; j++ ) {
                double value = statCovHist->GetBinContent(i,j); 
                if ( value != 0. ) {
                    skipThisBin = false; 
                }  
            } 
            // Through out bins with empty rows/columns
            if ( skipThisBin == true ) {
                binMap[i-1] = -1;
                numBinsToSkip++;
            } else {
                binMap[i-1] = bincounter;
                bincounter++;
            } 
        }
    }     

    // Create new Matrix
    int matrixdim = nbins - numBinsToSkip;
    TMatrixDSym statCovMat(matrixdim);
    
    
    
    // New Matrix
    // Beware the side bins of the problem
    // AND the side bins of the TH2D object 
    for ( int i = 2; i <= nbins-1 ; i++ ) {
        for ( int j = 2; j <= nbins-1 ; j++ ) {
            
            // Is this bin to be skipped?
            bool skipThisBin = false;
            if (binMap[i-1] == -1 ) skipThisBin = true;
            if (binMap[j-1] == -1 ) skipThisBin = true;       
            // Set Element
            if ( skipThisBin == false ) {
                double value = statCovHist->GetBinContent(i,j);
                int binnrI = binMap[i-1];
                int binnrJ = binMap[j-1];
                statCovMat[binnrI][binnrJ] = value;  
            }
        }
    }  
    
    // Determinant
    double detStatCovMat[1]; 
    
    // Invert the whole thing
    TMatrixDSym statCovMatInv = statCovMat;
    statCovMatInv.Invert(detStatCovMat);
    
    // Check Invertibility
    bool isInvertible = (bool) ((*detStatCovMat) != 0.);
    if ( isInvertible == false ) {
        cout << "Error in TopSVDFunctions::SVD_CalcGlobCorr() " << endl;
        cout << "Covariance Matrix cannot be inverted." << endl;
        cout << "Check the reason for this now." << endl;
        exit(1);
    }
     
 
    // Create new Histo for global correlation 
    TH1D* glcHist = new TH1D[1];
    glcHist->SetNameTitle("glcHist", "glcHist");
    bool isVariableBinSize =  statCovHist->GetXaxis()->IsVariableBinSize();
    if ( isVariableBinSize == true ) {
        const TArrayD* xbins = statCovHist->GetXaxis()->GetXbins();
        glcHist->SetBins(nbins, xbins->GetArray());
    } else {
        double xmin = statCovHist->GetXaxis()->GetXmin();
        double xmax = statCovHist->GetXaxis()->GetXmax();
        glcHist->SetBins(nbins, xmin, xmax); 
    }
     
    
    // Fill the histo you just created
    for ( int i = 1 ; i <= nbins ; i++ ) { 
        double glc = 0.; 
        
        // Find out the "true" bin number
        int binnr = binMap[i-1];
        
        // Skip bad bins
        bool skipThis = false;
        if ( binnr == -1 ) skipThis = true;
        
        
        // Run over good bins 
        if ( skipThis == false ) {
            double cov = statCovMat[binnr][binnr]; 
            double covinv = statCovMatInv[binnr][binnr]; 
            // The product cov*covinv should be greater than zero            
            double hhh = SVD_Divide(1., cov * covinv);
            double glcsq = 0.;
            if ( hhh > 0.) glcsq = 1. - hhh;
            glc = 100.*SVD_Sqrt(glcsq); 
        } else {
            glc = 0.;
        }
        
        // Set the value
        glcHist->SetBinContent(i, glc); 
        glcHist->SetBinError(i, 0.); 
    }
    
    // return histo
    return glcHist;
}  
 
// Calculate Averaged Global Correlation
double TopSVDFunctions::SVD_AvgGlobCorr(TH1D* globCorrHist)
{
    // Number of bins
    int nbins = globCorrHist->GetNbinsX();
    
    double sum = 0.;
    double avg = 0.;
    int bincounter = 0;
    
    // Do not use the side bins!
    // Neither of the TH1D class
    // nor the extra side bins
    for ( int i = 1; i <= nbins; i++ ) {
        if ( i == 1 ) continue;
        if ( i == nbins ) continue;
        double globcorr = globCorrHist->GetBinContent(i);
        sum += globcorr;
        bincounter++;
    }
    
    // Averaging
    if ( bincounter > 0 ) {
        avg = sum / ((double) bincounter);
    }
     
    return avg;
    
}


// Find minimum of a graph, given as two vectors.
// The coordinates of the minimum are returned in xBest and yBest
// The position in the vectors is returned as an integer.
int TopSVDFunctions::SVD_FindMinimum(TVectorD* xVec, TVectorD* yVec, double& xBest, double& yBest)
{
    //Test, if vectorX and vectorY have same number of entries
    int nelements = xVec->GetNoElements();
    int nelementsY = yVec->GetNoElements();
    if ( nelements <= 0 || nelementsY <= 0 || nelements != nelementsY ) {
        cout << "Error in TopSVDFunctions::SVD_FindMinimum()" << endl;
        cout << "   The two vectors have strange sizes." << endl;
        cout << "   xVec has " << nelements  << " entries." << endl;
        cout << "   yVec has " << nelementsY << " entries." << endl;
        cout << "Exiting ... " << endl;
        exit(1);
    }
    
     
    // Find smallest value
    double minYValue = FLT_MAX;
    double minXValue = 0.;
    int minPosition = -1; 
    for ( int i = 0 ; i < nelements ; i++ ) {
        double yvalue = (*yVec)[i];
        double xvalue = (*xVec)[i];
        if ( yvalue < minYValue ) {
            minYValue = yvalue;
            minXValue = xvalue;
            minPosition = i;
        }
    }
    
    
    // Return it all
    xBest = minXValue;
    yBest = minYValue;
    return minPosition;
    
}
 
 
// Helper Function
// Adds histos from an array to a stack
// Color Handling:
//      col=0   means: The colors will not be touched
//      col=2   means: All histos get the color 2
//      col=243 means: First Histogram red, second blue, third and all following ones green 
void TopSVDFunctions::SVD_Array2Stack(THStack* stack, TLegend* leg, TH1D* histo, TString label, TString drawOptions, TString legOptions, int col, int numHist, int lineSty)
{  
    // Existence of Objects
    if ( histo == NULL ) return; 
    
    // Calculate Colors
    bool doColors = false;
    int* theColors = NULL; 
    int numdigits = 1;
    if ( col > 0 ) {
        doColors = true; 
        int colstmp1 = col;
        while ( (colstmp1=colstmp1/10)>0 ) { numdigits++; }
        theColors = new int[numdigits]; 
        int colstmp2 = col;
        for ( int i = 0 ; i < numdigits ; i++ ) {  
            int position = numdigits-1-i;
            int lastdigit = colstmp2 % 10;
            theColors[position] = lastdigit;
            colstmp2=colstmp2/10; 
        }
    }
    
     
    // Loop over all histograms in the array
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        
        // Make a copy of the histogram before
        // adding it to the stack. This way you can
        // make changes to the histogram in the stack
        // without affecting the original histogram
        TH1D* stackHisto = (TH1D*) (histo+h)->Clone("stackHisto");
        stackHisto->SetDirectory(NULL);
        
        
        // Set Colors
        if ( doColors == true ) {
            int position = h;
            if ( position > numdigits-1 ) position = numdigits-1;
            int currentcolor = theColors[position];
            
            stackHisto->SetLineColor(currentcolor);
            stackHisto->SetMarkerColor(currentcolor);  
            
        }
        
        // Draw Options
        // For Systematic shifts, draw Options should ALWAYS
        // be just HIST
        TString theDrawOptions = drawOptions;
        if ( h > 0 ) {
            theDrawOptions = "HIST";
        }
        
        // Marker Style
        if ( drawOptions.Contains("P") == false ) {
            stackHisto->SetMarkerStyle(1);
            stackHisto->SetMarkerSize(0.);
        }  
        
        // Line Style
    if(lineSty<0) stackHisto->SetLineStyle(h+1);
    else          stackHisto->SetLineStyle(lineSty);  
         
        // Draw Options 
        stack->Add(stackHisto, theDrawOptions);


        
        // Legend entry
        // Only for first element
        if ( leg != NULL ) {
            if ( h == 0 ) {  
                TString theString = label;
                TString theLegOption = legOptions;
                if ( theLegOption.CompareTo("") == 0 ) theLegOption = "l"; //"lpf";
                leg->AddEntry(stackHisto, theString, theLegOption);
            }
        }
    }  
    
    
}


// Helper Function
// Adds histos from an array to a stack
void TopSVDFunctions::SVD_ArrayScale(TH1D* histo, double scale, int numHist)
{ 
    if ( histo == NULL ) return; 
    
    // Loop over all histograms in the array 
    for ( int h = 0 ; h < numHist ; h++ ) { 
        (histo+h)->Scale(scale);
    } 
    
}
 


// REFOLDING - without Background
// Two possibilities:
// If useWeights==false (Standard)
//     Divide "unfHist" by "xiniHist", then refold with event matrix "mcHist". 
//     ATTENTION: If there are no generated events in a generator bin, then
//     the division by "xini" doesnt work. But there is a way out:
//     In such a case, the whole column of the event matrix will be empty! 
//     Multiplication of x/xini with the event matrix will give zero.
//     This justifies, that those bins can just be skipped in the sum!
//     However: If the column of the event matrix is not empty for some
//     reason, than something is odd and you have to check this.
//     Overflow bins: Remember, all side bins are in regular binning.
//     Indexing: i=Rec=(xAxis in mcHist), j=Gen=(yAxis in mcHist)
// If useWeights==true
//     "unfHist" is taken to be the ratio between unfolded data and generated distribution,
//     which means "unfHist" is directly refolded with event matrix "mcHist". 
//     "xiniHist" will not be used, can be set to NULL if you wish.
// Attention: A new histogram ist created on the heap
TH1D* TopSVDFunctions::SVD_Refold(TH1D* unfHist, TH1D* xiniHist, TH2D* mcHist, bool useWeights, int numHist) 
{
    
    // Existence of Objects    
    if ( unfHist == NULL ) return NULL;
    if ( mcHist == NULL ) return NULL;
    if ( xiniHist == NULL && useWeights == false ) return NULL;
    
    // Create new Histograms
    TH1D* refoldHist = SVD_CloneHists1D(unfHist, numHist);
    int nbins = unfHist->GetNbinsX(); 
 
    // Loop over all histograms in the array 
    for ( int h = 0 ; h < numHist ; h++ ) {  
        
        // Important Check:
        // Make sure that generator level bins
        // with empty xiniHist-bins
        // have also empty mcHist-rows
        // GenLevelIndex=j, DetLevelIndex=i
        if ( useWeights == false ) {
            for ( int j = 1 ; j <= nbins ; j++ ) {
                double numerator = (xiniHist+h)->GetBinContent(j);
                if ( numerator > 0. ) continue; 
                bool doWarning = true;
                for ( int i = 1 ; i <= nbins ; i++ ) {  
                    double matrixelement = (mcHist+h)->GetBinContent(i,j); 
                    if ( matrixelement == 0. ) continue; 
                    if ( doWarning == true ) {
                        cout << endl;
                        cout << "******************************************************************" << endl;
                        cout << "WARNING! WARNING! WARNING! " << endl;
                        cout << "Severe problem in TopSVDFunctions::SVD_Refold() !!!" << endl;
                        cout << "For the generator level bin " << j << " ... " << endl;
                        cout << "    ... xiniHist is empty ... " << endl;
                        cout << "    ... but in the response matrix there are non-empty elements." << endl;
                        cout << "Unfolding might work, but refolding won't !!!" << endl;
                        cout << "This means, your chi square will be busted and no one will believe " <<  endl;
                        cout << "your results !!!" << endl;
                        cout << "******************************************************************" << endl;
                        cout << endl; 
                        doWarning = false;
                    }
                }
            }
        }      
        
        // Run over all bins
        // GenLevelIndex=j, DetLevelIndex=i 
        for ( int i = 1 ; i <= nbins ; i++ ) { 
            double refoldedvalue = 0.; 
            for ( int j = 1 ; j <= nbins ; j++ ) { 
                double xini = 1.;
                if ( useWeights == false ) {
                    xini = (xiniHist+h)->GetBinContent(j);
                }
                double unf = (unfHist+h)->GetBinContent(j); 
                double matrixelement = (mcHist+h)->GetBinContent(i,j); 
                double current = 0.; 
                if ( xini > 0. ) {  
                    current = matrixelement * ( unf / xini );
                } else {
                    current = 0.;
                }
                refoldedvalue += current;
            }
            (refoldHist+h)->SetBinContent(i, refoldedvalue);
            (refoldHist+h)->SetBinError(i, 0.);
        }
    }
     
    // Return the histogram
    return refoldHist;

}

// Calculate global lumi scale factor from reconstructed
// distributions. Exclude overflow bins! Use background free
// data distribution!
double TopSVDFunctions::SVD_LumiScaleFactor(TH1D* dataHist, TH1D* recHist)
{
    int nbins = dataHist->GetNbinsX();
    double sumData = 0.;
    double sumMC = 0.;
    
    for ( int i = 1 ; i <= nbins ; i++ ) {
        sumData = sumData + dataHist->GetBinContent(i);
        sumMC = sumMC + recHist->GetBinContent(i);
    }
    
    if ( sumMC <= 0. ) return 1.;
    if ( sumData <= 0. ) return 1.;
    
    double scalefactor = sumData / sumMC ;
    return scalefactor;
    
}

// Count the number of weighted events in a histo
double TopSVDFunctions::SVD_Integral1D(TH1D* hist, int syst, bool doOF)
{ 
    if ( hist == NULL ) return 0.; 
    int nbins = hist->GetNbinsX();
    double sumData = 0.;
    
    // Sum it up
    for ( int i = 0 ; i <= nbins+1 ; i++ ) {
        if ( i == 0       && doOF == false ) continue;
        if ( i == nbins+1 && doOF == false ) continue;
        double value = (hist+syst)->GetBinContent(i);
        bool NaN_found = false;
        if ( std::isnan(value) == true ) NaN_found = true;
        if ( NaN_found == true ) {
        	cout << endl;
        	cout << "Error in TopSVDFunctions::SVD_Integral()" << endl;
        	cout << "NaN Found in Bin " << i << endl;
        	cout << "Exiting ... " << endl;
        	cout << endl;
        	exit(1);
        }
         
       
       
        sumData += value;
    }
        
    return sumData;    
}

// Count the number of weighted events in a histo
double TopSVDFunctions::SVD_Integral2D(TH2D* hist, int syst, bool doOF)
{
    if ( hist == NULL ) return 0.;
    int nbinsx = hist->GetNbinsX();
    int nbinsy = hist->GetNbinsY();
    double sumData = 0.;
    
    // Sum it up
    for ( int i = 0 ; i <= nbinsx+1 ; i++ ) {
        if ( i == 0        && doOF == false ) continue;
        if ( i == nbinsx+1 && doOF == false ) continue;
        for ( int j = 0 ; j <= nbinsy+1 ; j++ ) {
            if ( j == 0        && doOF == false ) continue;
            if ( j == nbinsy+1 && doOF == false ) continue;
            double value = (hist+syst)->GetBinContent(i,j);
	        bool NaN_found = false;
	        if ( std::isnan(value) == true ) NaN_found = true;
	        if ( NaN_found == true ) {
	        	cout << endl;
	        	cout << "Error in TopSVDFunctions::SVD_Integral()" << endl;
	        	cout << "NaN Found in Bin " << i << "/" << j << endl;
	        	cout << "Exiting ... " << endl;
	        	cout << endl;
	        	exit(1);
	        }
            sumData += value; 
            
        }
    }
        
    return sumData;    
}

// Integrates histos in an array
// Creates a new array of doubles on the heap.
// Don't forget to delete it sometimes.  
double* TopSVDFunctions::SVD_ArrayIntegral1D(TH1D* hist, bool doOF, int numHist)
{
    // Existence of Objects 
    if ( numHist < 1 ) return NULL; 
 
    // Craeate Object
    double* array = new double[numHist];
     
     
    // Full Array 
    for ( int h = 0 ; h < numHist ; h++ ) { 
        double integral = SVD_Integral1D(hist, h, doOF); 
        *(array+h) = integral; 
    } 
 
     // return
     return array;
 
}

// Find bin number on an axis
// You give a bin boundary 'bnd' and this function gives
// you the number of the bin that STARTS at 'bnd'
int TopSVDFunctions::SVD_FindBin(TAxis* axis, double bnd)
{
    // First attempt to find bin number
    int binnr = axis->FindBin(bnd);
    
    // Possibly, the binnr is by 1 too low.
    // This can happen due to rounding errors
    // busting the comparison between 'bnd' and the
    // bin boundaries of the axis.
    // To make sure that we found the right bin,
    // we check which bin boundary is actually closer.
    double lowxFirstTry  = axis->GetBinLowEdge(binnr);
    double lowxStepUp    = axis->GetBinLowEdge(binnr+1);
    double distanceFirstTry = TMath::Abs(bnd - lowxFirstTry);
    double distanceStepUp = TMath::Abs(bnd - lowxStepUp);
    if ( distanceStepUp < distanceFirstTry ) binnr++;
    
    // Return
    return binnr;
     
    
}


// Setup Unfolding Tool
// This creates an array of TopSVDUnfold objects on the heap.
// Each of them is set up with a different systematic
// All other settings will be the same.
TopSVDUnfold* TopSVDFunctions::SVD_SetupTool(TH1D* dataHist, TH1D* biniHist, TH1D* xiniHist, TH2D* mcHist, double theTau, int numHist)
{ 
     
    
    TopSVDUnfold* mySVDUnfold = new TopSVDUnfold[numHist];
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) { 
        
        // Create Unfolding Tool
        TopSVDUnfold tmpTool(dataHist+h, biniHist+h, xiniHist+h, mcHist+h); 
        *(mySVDUnfold+h) = tmpTool; 
        
        // Setup
        (mySVDUnfold+h)->SetNormalize(false); 
        (mySVDUnfold+h)->SetTau(theTau);
        
    }
    
    // Return
    return mySVDUnfold;
}


// Call TopSVDUnfold sequentially
// Creates new objects on the heap!
TH1D* TopSVDFunctions::SVD_CallSVDUnfolding(TopSVDUnfold* mySVDUnfold, int theKReg, int numHist)
{  
    // Existence of Objects    
    if ( mySVDUnfold == NULL ) return NULL;  
    
    // Create new histograms
    TH1D* unfHist = new TH1D[numHist]; 
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) { 
        TH1D* tmpHist = (mySVDUnfold+h)->Unfold(theKReg); 
        tmpHist->Copy(*(unfHist+h)); 
        (unfHist+h)->SetName("unfHist"); 
        
        // Make sure the new histogram is not on a TDirectory
        // You will most likely run into garbage collector problems,
        (unfHist+h)->SetDirectory(NULL);
    } 
    
    return unfHist;
}


// Retrieve Weights
// Creates new objects on the heap! 
TH1D* TopSVDFunctions::SVD_RetrieveWeights(TopSVDUnfold* mySVDUnfold, int numHist)
{  
    // Existence of Objects    
    if ( mySVDUnfold == NULL ) return NULL;
    
    // Create new histograms
    TH1D* weightHist = new TH1D[numHist];
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) { 
        TH1D* tmpHist = (mySVDUnfold+h)->GetWeights();
        tmpHist->Copy(*(weightHist+h));
        (weightHist+h)->SetName("weightHist");
        
        // Make sure the new histogram is not on a TDirectory
        // You will most likely run into garbage collector problems,
        (weightHist+h)->SetDirectory(NULL);
    } 
    
    return weightHist;
}


// Retrieve D Value Plot
// Creates new objects on the heap! 
TH1D* TopSVDFunctions::SVD_RetrieveDD(TopSVDUnfold* mySVDUnfold, int numHist)
{  
    // Existence of Objects    
    if ( mySVDUnfold == NULL ) return NULL;
    
    // Create new histograms
    TH1D* ddHist = new TH1D[numHist];
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) { 
        TH1D* tmpHist = (mySVDUnfold+h)->GetD();
        tmpHist->Copy(*(ddHist+h));
        (ddHist+h)->SetName("ddHist");
        
        // Make sure the new histogram is not on a TDirectory
        // You will most likely run into garbage collector problems,
        (ddHist+h)->SetDirectory(NULL);
    } 
    
    return ddHist;
}


// Retrieve Singular Values
// Creates new objects on the heap! 
TH1D* TopSVDFunctions::SVD_RetrieveSV(TopSVDUnfold* mySVDUnfold, int numHist)
{ 
    // Existence of Objects    
    if ( mySVDUnfold == NULL ) return NULL;
    
    
    // Create new histograms
    TH1D* svHist = new TH1D[numHist];
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) { 
        TH1D* tmpHist = (mySVDUnfold+h)->GetSV();
        tmpHist->Copy(*(svHist+h));
        (svHist+h)->SetName("svHist");
        
        // Make sure the new histogram is not on a TDirectory
        // You will most likely run into garbage collector problems,
        (svHist+h)->SetDirectory(NULL);
    } 
    
    return svHist;
}


// Take an array of systematics
// and make relative shifts in % out of it
// NOTE: The output array will have one histogram LESS!
TH1D* TopSVDFunctions::SVD_ArrayToShifts(TH1D* array, int numHist)
{
    // Existence of Objects    
    if ( array == NULL ) return NULL;
    if ( numHist < 2 ) return NULL;
    
    // Create new histograms
    // ... one less 
    TH1D* output = new TH1D[numHist-1];
    
    // Loop over all histograms
    // ... exclude the first
    for ( int h = 1 ; h < numHist ; h++ ) { 
        
        // Substract the nominal values
        TH1D* absDiv = (TH1D*) (array+h)->Clone("tmp");
        absDiv->SetDirectory(NULL);
        absDiv->Add(array, -1.);
        
        
        // Division by nominal values
        TH1D* relDiv = SVD_MakeRatioZero(absDiv, array, 1);
        
        // Scale to percent values
        relDiv->Scale(100.);
        
        // Save it in new array
        *(output-1+h) = *relDiv;
        
        
        
        // Make sure the new histogram is not on a TDirectory
        // You will most likely run into garbage collector problems,
        (output-1+h)->SetDirectory(NULL);
        
    }
    
    return output;
    
}

// BACKGROUND HANDLING
// The raw data are given plus the simulated
// backgrounds from MC. The result is returned in
// the placeholder parameter 'dataHist'.
// Notice:
// (1) Errors from background will be neglected! 
// (2) If 'ttbgrInputHist' is given, it will be removed by
//     a multiplicative approach that is xSec independent!   
// (3) If 'bgrInputHist' is given, it will be removed from the
//     data by substraction.
// (4) If both 'bgrInputHist' and 'ttbgrInputHist' are given,
//     the difference between 'bgrInputHist' and 'ttbgrInputHist'
//     will be substracted from the data and the background given
//     int 'ttbgrInputHist' will be treated as described in (2).
// (5) If neither 'bgrInputHist' nor 'ttbgrInputHist' are given,
//     no background will be considered.
void TopSVDFunctions::SVD_BackgrHandling(TH1D*& dataHist, TH1D* bgrHist, TH1D* ttbgrHist, TH1D* biniHist, TH1D* rawHist, bool beForgiving, int numHist)
{ 
    
    int flag_verbose = 1; 
    

    // Steer background handling
    bool doBgr = false;
    if ( bgrHist != NULL ) doBgr = true;
    bool doTtBgr = false;
    if ( ttbgrHist != NULL ) doTtBgr = true;
    
    // Nbins
    int nbins = rawHist->GetNbinsX();
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
         
        // loop bins
        for ( int i = 1 ; i <= nbins ; i++) { 

            // get bin and errors for data
            // Notice: Only ONE histogram will be used here!
            double value_data = rawHist->GetBinContent(i);
            double err_data = rawHist->GetBinError(i);  
            // Get background value to be substracted! 
            double value_bgr = 0.; 
            if ( doBgr == true   ) value_bgr = (bgrHist+h)->GetBinContent(i);  
            // if ttbar signal fraction method is used: consider only non ttbar BG in check
            double relevant_bgr = value_bgr;
            if ( doTtBgr == true ) relevant_bgr-=(ttbgrHist+h)->GetBinContent(i);
            // check for data to be smaller as BG that will be subtracted
            if(value_data<relevant_bgr){
                if ( beForgiving == false ) {
                    std::cout << "ERROR in TopSVDFunctions::SVD_BackgrHandling: " << std::endl;
                    std::cout << "N_MC BG > N_data in bin " << i << std::endl;
                    for ( int k = 1 ; k <= nbins ; k++ ) {
                        std::cout << "   Bin " << k << ":  Data=" << rawHist->GetBinContent(k) << ",    Bgr=";
                        double relevant_bgr_bin =(bgrHist+h)->GetBinContent(k);
                        if ( doTtBgr == true ) relevant_bgr_bin-=(ttbgrHist+h)->GetBinContent(k);
                        std::cout << relevant_bgr_bin;
                        if ( relevant_bgr > value_data ) {
                            cout << " !@#$%^&*! " << endl;
                        } 
                        else {
                            cout << endl;
                        }
                    }  
                    exit(1);
                } else {
                    relevant_bgr = value_data;
                } 
            }
                 
            // Debug output
            if(flag_verbose>=3){
                std::cout << std::endl << "bin " << i << std::endl;
                std::cout << "value_data: " << value_data << std::endl;
                std::cout << "value_bgr: "  << value_bgr << std::endl;
                std::cout << "err_data: "   << err_data << std::endl;
            }
                 
            
            // Get background scale factor!
            double sigFrac = 1.;
            if ( doTtBgr == true ) {
                double value_ttSig=(biniHist+h)->GetBinContent(i);
                double value_ttBgr=(ttbgrHist+h)->GetBinContent(i);  
                if(value_ttBgr>value_bgr){
                    if ( beForgiving == false ) {
                        std::cout << "ERROR in TopSVDFunctions::SVD_BackgrHandling: " << std::endl;
                        std::cout << "The TtBar-Background is larger than the complete Background!" << std::endl;
                        std::cout << "This happens in bin " << i << " (range " << rawHist->GetBinLowEdge(i) << ",";
                        std::cout << rawHist->GetBinLowEdge(i+1) << " )" << std::endl;
                        std::cout << "    TTBar-Bgr: " << value_ttBgr << std::endl;
                        std::cout << "    Bgr:       " << value_bgr << std::endl;
                        exit(0);
                    } else {
                        value_ttBgr=value_bgr;
                    }
                }
                
                // don't subtract ttbar BG !!!!!!!
                value_bgr-=value_ttBgr;
                
                // calculate signal fraction
                if ( value_ttSig+value_ttBgr > 0. ) {
                    sigFrac=value_ttSig/(value_ttSig+value_ttBgr);
                }
                  
                // debug output
                if(flag_verbose>=3){
                    std::cout << "value_ttBgr: " << value_ttBgr << std::endl;
                    std::cout << "value_bgr(no ttbar): " << value_bgr << std::endl;
                    std::cout << "value_ttSig: " << value_ttSig << std::endl;
                    std::cout << "sigFrac: " << sigFrac << std::endl;
                } 
            
            }
                
            // New Values
            double value_new = (value_data - value_bgr)*sigFrac;
            double err_new = err_data*sigFrac; 
             
              
            // debug output
            if(flag_verbose>=3){
                std::cout << "data_value_new: " << value_new << std::endl;
                std::cout << "data_err_new: "   << err_new   << std::endl;
            }
              
            // Save new values
            (dataHist+h)->SetBinContent(i, value_new);
            (dataHist+h)->SetBinError(i, err_new); 
              
        }
    }   
}   
    
    
// Calculate Purity
// Will be stored in 'purHist'
// Note:
// (1) Purity is set to zero, if number of reconstructed
//     events in a bin is zero
// (2) The purity is given in PERCENT
void TopSVDFunctions::SVD_Pur(TH1D*& purHist, TH2D* mcHist, int numHist)
{   
     
    // Existence of Objects 
    if ( mcHist == NULL ) return;   
    
    // Nbins
    int nbins = mcHist->GetNbinsX();
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Loop over all bins
        for ( int i = 1 ; i <= nbins ; i++ ) {
            double numerator = (mcHist+h)->GetBinContent(i, i);
    
            // Sum over all Gen Level bins INCLUDING OF bins
            // NOTE: The generator level quantities are
            // ... on the Y Axis
            // ... that is the second index
            // ... which is index j here.
            double denominator = 0.;
            for ( int j = 1 ; j < nbins+1 ; j++ ) {
                denominator += (mcHist+h)->GetBinContent(i, j);
            }
    
            //double denominator = biniHist->GetBinContent(i);
            double purity = SVD_Divide(numerator, denominator);
            purity *= 100.; 
            
            // Save it
            (purHist+h)->SetBinContent(i,purity);
            (purHist+h)->SetBinError(i, 0.); 
        }
    }
}
    
    
// Sum over Elements in Response Matrix    
// INCLUDING OF Bins
void TopSVDFunctions::SVD_SumOverGenBins(TH1D*& sumHist, TH2D* mcHist, int numHist)
{   
     
    // Existence of Objects 
    if ( mcHist == NULL ) return;   
    
    // Nbins
    int nbins = mcHist->GetNbinsX();
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Loop over all bins
        for ( int i = 1 ; i <= nbins ; i++ ) {
    
            // Sum over all Gen Level bins INCLUDING OF bins
            // NOTE: The generator level quantities are
            // ... on the Y Axis
            // ... that is the second index
            // ... which is index j here.
            double sum = 0.;
            for ( int j = 1 ; j < nbins+1 ; j++ ) {
                sum += (mcHist+h)->GetBinContent(i, j);
            }
            
            // Save it
            (sumHist+h)->SetBinContent(i,sum);
            (sumHist+h)->SetBinError(i, 0.); 
        }
    }
}
    
    
// Sum over Elements in Response Matrix
// INCLUDING OF Bins    
void TopSVDFunctions::SVD_SumOverRecBins(TH1D*& sumHist, TH2D* mcHist, int numHist)
{   
     
    // Existence of Objects 
    if ( mcHist == NULL ) return;   
    
    // Nbins
    int nbins = mcHist->GetNbinsX();
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Loop over all bins
        for ( int i = 1 ; i <= nbins ; i++ ) {
    
            // Sum over all Rec Level bins INCLUDING OF bins
            // NOTE: The rec level quantities are
            // ... on the X Axis
            // ... that is the first index
            // ... which is index j here.
            double sum = 0.;
            for ( int j = 1 ; j < nbins+1 ; j++ ) { 
                sum += (mcHist+h)->GetBinContent(j, i);
            } 
            
            // Save it
            (sumHist+h)->SetBinContent(i,sum);
            (sumHist+h)->SetBinError(i, 0.); 
        }
    }
}
 

// Calculate Smearin
// Will be stored in 'smearinHist'
void TopSVDFunctions::SVD_LowerSmearin(TH1D*& smearinHist, TH2D* mcHist, int numHist)
{   
     
    // Existence of Objects 
    if ( mcHist == NULL ) return;   
    
    // Nbins
    int nbins = mcHist->GetNbinsX();
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Loop over all bins
        for ( int i = 1 ; i <= nbins ; i++ ) {
            
            // Calculate Smearin
            // NOTE: The generator level quantities are
            // ... on the Y Axis
            // ... that is the second index
            double smearin = 0.;
            smearin += (mcHist+h)->GetBinContent(i, 1); 
                    
            // Save it
            (smearinHist+h)->SetBinContent(i,smearin);
            (smearinHist+h)->SetBinError(i, 0.); 
        }
    }
}



// Calculate Smearin
// Will be stored in 'smearinHist'
void TopSVDFunctions::SVD_UpperSmearin(TH1D*& smearinHist, TH2D* mcHist, int numHist)
{   
     
    // Existence of Objects 
    if ( mcHist == NULL ) return;   
    
    // Nbins
    int nbins = mcHist->GetNbinsX();
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Loop over all bins
        for ( int i = 1 ; i <= nbins ; i++ ) {
            
            // Calculate Smearin
            // NOTE: The generator level quantities are
            // ... on the Y Axis
            // ... that is the second index
            double smearin = 0.; 
            smearin += (mcHist+h)->GetBinContent(i, nbins);
                    
            // Save it
            (smearinHist+h)->SetBinContent(i,smearin);
            (smearinHist+h)->SetBinError(i, 0.); 
        }
    }
}   
       
// Calculate Smearout
// Will be stored in 'smearoutHist'
void TopSVDFunctions::SVD_LowerSmearout(TH1D*& smearoutHist, TH2D* mcHist, int numHist)
{   
     
    // Existence of Objects 
    if ( mcHist == NULL ) return;   
    
    // Nbins
    int nbins = mcHist->GetNbinsX();
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Loop over all bins
        for ( int i = 1 ; i <= nbins ; i++ ) {
            
            // Calculate smearout
            // NOTE: The generator level quantities are
            // ... on the Y Axis
            // ... that is the second index
            double smearout = 0.;
            smearout += (mcHist+h)->GetBinContent(1, i);    
                    
            // Save it
            (smearoutHist+h)->SetBinContent(i,smearout);
            (smearoutHist+h)->SetBinError(i, 0.); 
        }
    }
} 
       
// Calculate Smearout
// Will be stored in 'smearoutHist'
void TopSVDFunctions::SVD_UpperSmearout(TH1D*& smearoutHist, TH2D* mcHist, int numHist)
{   
     
    // Existence of Objects 
    if ( mcHist == NULL ) return;   
    
    // Nbins
    int nbins = mcHist->GetNbinsX();
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Loop over all bins
        for ( int i = 1 ; i <= nbins ; i++ ) {
            
            // Calculate smearout
            // NOTE: The generator level quantities are
            // ... on the Y Axis
            // ... that is the second index
            double smearout = 0.; 
            smearout += (mcHist+h)->GetBinContent(nbins, i); 
                    
            // Save it
            (smearoutHist+h)->SetBinContent(i,smearout);
            (smearoutHist+h)->SetBinError(i, 0.); 
        }
    }
}      
    
    
// Calculate Stability
// Will be stored in 'stabHist'
// Note: 
// (1) Stability is set to zero, if there are no
//     gerenated events in a bin
// (2) If there is no reconstructed event in this bin,
//     then the stability is naturally zero, too
// (3) The stability is given in PERCENT
void TopSVDFunctions::SVD_Stab(TH1D*& stabHist, TH2D* mcHist, int numHist)
{ 
    // Existence of Objects 
    if ( mcHist == NULL ) return;   
      
    // Nbins
    int nbins = mcHist->GetNbinsX();
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Loop over all bins
        for ( int i = 1 ; i <= nbins ; i++ ) {
            double numerator = (mcHist+h)->GetBinContent(i, i);
    
            // Sum over all Rec Level bins EXCEPT OF bins
            // NOTE: The rec level quantities are
            // ... on the X Axis
            // ... that is the first index
            // ... which is index j here.
            double denominator = 0.;
            for ( int j = 1 ; j < nbins+1 ; j++ ) {
                if ( j == 1 ) continue;
                if ( j == nbins ) continue;
                denominator += (mcHist+h)->GetBinContent(j, i);
            }
                 
            double stability = SVD_Divide(numerator, denominator);
            stability *= 100.; 
            
            // Save it
            (stabHist+h)->SetBinContent(i,stability);
            (stabHist+h)->SetBinError(i, 0.);
        }
    }
}


// Calculate Real Efficiency    
// Will be stored in 'effHist'
// Attention:
// (1) If no generated events are available in a bin,
//     then the efficiency is set to zero.
// (2) Statistical Errors are calculated in a binomial manner.
//     This is different from bbb Efficiency.
// (3) The efficiency is given in PERCENT
void TopSVDFunctions::SVD_Eff(TH1D*& effHist, TH2D* mcHist, TH1D* xiniHist, int numHist)
{

    // Existence of Objects 
    if ( mcHist == NULL ) return;  
    if ( xiniHist == NULL ) return;  
    
    
    // Nbins
    int nbins = mcHist->GetNbinsX(); 
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Loop over all bins
        for ( int i = 1 ; i <= nbins ; i++ ) {
            
            double denominator = (xiniHist+h)->GetBinContent(i);
            double numerator = 0.;
            for ( int j = 1; j <= nbins ; j++ ) {
                // Sum over all Rec Level Bins EXEPT OF bins!
                // So Please, please, please do not forget to exclude
                // the side bins
                // Note: Rec Level bins are on x axis of mcHist
                double elem = (mcHist+h)->GetBinContent(j, i);
                if ( j == 1 ) elem = 0.;
                if ( j == nbins ) elem = 0.;
                numerator += elem;
            }
            double efficiency = SVD_Divide(numerator,denominator);
            double effErr = SVD_BinomialError(numerator, denominator);
            
            // scale to PERCENT
            efficiency *= 100.;
            effErr *= 100.; 
            
            // Save it
            (effHist+h)->SetBinContent(i,efficiency);
            (effHist+h)->SetBinError(i, effErr);
        } 
    }
}


// Calculate BBB Efficiency    
// Will be stored in 'beffHist'
// Attention:
//     For the statistical uncertainty, relative errors from 
//     reconstructed and generated events counts are added in quadrature.
//     This is different from the errors on the real efficiency,
void TopSVDFunctions::SVD_BBBEff(TH1D*& beffHist, TH1D* biniHist, TH1D* xiniHist, int numHist)
{

    // Existence of Objects 
    if ( biniHist == NULL ) return;  
    if ( xiniHist == NULL ) return;  
    
    
    // Nbins
    int nbins = biniHist->GetNbinsX(); 
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Loop over all bins
        for ( int i = 1 ; i < nbins+1 ; i++ ) {
              
            // Sum over all Gen Level bins Including OF bins
            double numerator = (biniHist+h)->GetBinContent(i);  
            double numeratorErr = (biniHist+h)->GetBinError(i);
            double denominator = (xiniHist+h)->GetBinContent(i); 
            double denominatorErr = (xiniHist+h)->GetBinError(i);
            double bbbeff = 0.;
            double bbbeffErr = 0.;
            if ( denominator > 0. && numerator > 0.) {
                bbbeff = numerator / denominator;
                double bbbeffErrRelSq = 0.;
                bbbeffErrRelSq += TMath::Power(numeratorErr/numerator, 2.);
                bbbeffErrRelSq += TMath::Power(denominatorErr/denominator, 2.);
                bbbeffErr = bbbeff * TMath::Sqrt(bbbeffErrRelSq);
                
                // Scale to PERCENT
                bbbeff *= 100.;
                bbbeffErr *= 100.;
            }
            
            // Save it
            (beffHist+h)->SetBinContent(i,bbbeff);
            (beffHist+h)->SetBinError(i,bbbeffErr);
             
        }
    }
}


// Do BBB Unfolding 
// Output stored in 'bbbHist'
// ATTENTION:
// (1) Use BBB Efficiency.
// (2) Although this might be tempting,
//     do not include the error from the efficiency.
//     This would bust your BBB / SVD comparison.
//     Instead, add the error on the efficiency
//     later on as a systematic MC uncertainty! 
void TopSVDFunctions::SVD_BBBUnf(TH1D*& bbbHist, TH1D* dataHist, TH1D* beffHist, int numHist)
{
    
    // Existence of Objects 
    if ( dataHist == NULL ) return;  
    if ( beffHist == NULL ) return; 
    
    // Nbins
    int nbins = beffHist->GetNbinsX(); 
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Loop over all bins  
        for ( int i = 1 ; i <= nbins ; i++ ) {
            
            // Get Input Data
            // Notice: BBB Efficiency is given in PERCENT
            double data = (dataHist+h)->GetBinContent(i);
            double bbbeff  = 0.01 * (beffHist+h)->GetBinContent(i);
            double dataErr = (dataHist+h)->GetBinError(i);
            double bbbunf = 0.;
            double bbbunfErr = 0.;
            if ( data > 0. && bbbeff > 0.) {
                bbbunf = data / bbbeff;
                bbbunfErr = dataErr / bbbeff;
            }
            (bbbHist+h)->SetBinContent(i,bbbunf);
            (bbbHist+h)->SetBinError(i,bbbunfErr);
        } 
    }
}


// Converts the event matrix into the probability matrix
// Attention:
// (1) Here we assume, that unfolding side bins have already bin added 
// (2) We assume, that the reconstructed quantities are on the x axis
//     (first index) and the generated quantities on the y axis (second
//     index).  
// (3) The Probabilities will be given in PERCENT!
TH2D* TopSVDFunctions::SVD_ProbMatrix(TH2D* mcHist, TH1D* genHist, int numHist)  
{  
    // Existence of Objects
    if ( mcHist == NULL ) return NULL;
    if ( genHist == NULL ) return NULL;   
    
    // Create new objects
    TH2D* outputhists = new TH2D[numHist];
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Bins
        bool isVariableBinSizeX = (mcHist+h)->GetXaxis()->IsVariableBinSize(); 
        bool isVariableBinSizeY = (mcHist+h)->GetXaxis()->IsVariableBinSize(); 
        int binsx = (mcHist+h)->GetNbinsX();
        int binsy = (mcHist+h)->GetNbinsY();
        if ( isVariableBinSizeX == true && isVariableBinSizeY == true ) { 
            const double* theBinsx = (mcHist+h)->GetXaxis()->GetXbins()->GetArray(); 
            const double* theBinsy = (mcHist+h)->GetYaxis()->GetXbins()->GetArray(); 
            (outputhists+h)->SetBins(binsx, theBinsx, binsy, theBinsy); 
        } else {  
            double xmin = (mcHist+h)->GetXaxis()->GetXmin();
            double xmax = (mcHist+h)->GetXaxis()->GetXmax();
            double ymin = (mcHist+h)->GetYaxis()->GetXmin();
            double ymax = (mcHist+h)->GetYaxis()->GetXmax(); 
            (outputhists+h)->SetBins(binsx, xmin, xmax, binsy, ymin, ymax);
        } 
        
        // Loop over all bins
        for ( int i = 1 ; i <= binsx ; i++ ) {
            for ( int j = 1 ; j <= binsy ; j++ ) {
                  
                // Attention: 
                // Detector Level == xAxis == First Index == Index i
                // Generator Level == yAxis == Second Index == Index j
                double mcHistElement = (mcHist+h)->GetBinContent(i,j);
                double generatedEvents = (genHist+h)->GetBinContent(j);
                
                // Divide
                double probHistElement = SVD_Divide(mcHistElement, generatedEvents);
                
                // Scale to PERCENTs
                probHistElement *= 100.;
                
                // Fill it
                (outputhists+h)->SetBinContent(i, j, probHistElement);
                (outputhists+h)->SetBinError(i, j, 0.); 
            }
        } 
        
        // Make sure the new histogram is not on a TDirectory
        // You will most likely run into garbage collector problems,
        (outputhists+h)->SetDirectory(NULL);
    }
    
    // REturn
    return outputhists;
    
    
} 


// Gives the diagonal elements of the Probability Matrix 
// See documentation there.
TH1D* TopSVDFunctions::SVD_ProbHist(TH2D* mcHist, TH1D* genHist, int numHist)  
{ 

    // Existence of Objects
    if ( mcHist == NULL ) return NULL;
    if ( genHist == NULL ) return NULL;    
    
    // Prob Matrix
    TH2D* probMatrix = SVD_ProbMatrix(mcHist, genHist, numHist);
       
    // Get diagonal elements
    TH1D* outputhists = SVD_DiagElements(probMatrix, numHist);
  
      // Remove obsolete Objects
    SVD_DeleteHists2D(probMatrix, numHist);
   
    // REturn
    return outputhists; 
} 


// Write Array of 1D histograms to Root file
void TopSVDFunctions::SVD_WriteHists1D(TH1D* hists, int numHist)
{
    // Existence of Objects
    if ( hists == NULL ) return; 
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Root File
        (hists+h)->Write();
        
    } 
}


// Write Array of 2D histograms to Root file
void TopSVDFunctions::SVD_WriteHists2D(TH2D* hists, int numHist)
{
    // Existence of Objects
    if ( hists == NULL ) return;
    
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        (hists+h)->Write();
    } 
}  


// Write Array of 1D histograms to ASCII file
void TopSVDFunctions::SVD_Hists1DToASCII(TH1D* hists, TString textOutputFolderName, TString numberFormat, TString errSep, TString lineSep, int numHist)
{
    // Existence of Objects
    if ( hists == NULL ) return; 
    
    // Get Histogram name
    TString histoname = hists->GetName();
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Create File     
        TString filename = "";
        filename.Append(histoname);
        if ( h > 0 ) {
            filename.Append("_");
            filename.Append(h);
        } 
        filename.Append(".txt");       
        SVD_MakeFolder(textOutputFolderName);
        filename.Prepend("/");
        filename.Prepend(textOutputFolderName);
//         SVD_RemoveFile(filename);
//         SVD_TouchFile(filename);
        
        // Size of Histogram
        int nbins = hists->GetNbinsX();
        
        
        // Fill that new file with numbers
        FILE* outputFile; 
        outputFile = fopen(filename, "w"); 
        for ( int i = 1 ; i <= nbins ; i++ ) {
            double value = (hists+h)->GetBinContent(i);
            double error = (hists+h)->GetBinError(i);
            
            TString FormatString = "";
            FormatString.Append(numberFormat);
            FormatString.Append(errSep);
            FormatString.Append(numberFormat);
            FormatString.Append(lineSep);
              
            fprintf(outputFile, FormatString, value, error); 
        
        }
        fclose(outputFile);
    } 
}

// Write Array of 1D histograms to ASCII file
void TopSVDFunctions::SVD_Hists2DToASCII(TH2D* hists, TString textOutputFolderName, TString numberFormat, TString errSep, TString lineSep, TString colSep, int numHist)
{
    // Existence of Objects
    if ( hists == NULL ) return; 
    
    // Get Histogram name
    TString histoname = hists->GetName();
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Create File
        TString filename = "";
        filename.Append(histoname);
        if ( h > 0 ) {
            filename.Append("_");
            filename.Append(h);
        } 
        filename.Append(".txt");        
        SVD_MakeFolder(textOutputFolderName);
        filename.Prepend("/");
        filename.Prepend(textOutputFolderName);
        
//         SVD_RemoveFile(filename);
//         SVD_TouchFile(filename);
        
        // Size of Histogram
        int nbins = hists->GetNbinsX();
        
        
        // Fill that new file with numbers
        FILE* outputFile; 
        outputFile = fopen(filename, "w"); 
        for ( int i = 1 ; i <= nbins ; i++ ) {
            for ( int j = 1 ; j <= nbins ; j++ ) {
                double value = (hists+h)->GetBinContent(i, j);
                double error = (hists+h)->GetBinError(i, j);
                
                TString FormatString = "";
                FormatString.Append(numberFormat);
                FormatString.Append(errSep);
                FormatString.Append(numberFormat);
                if (j < nbins )  FormatString.Append(colSep);
                if ( j == nbins) FormatString.Append(lineSep);
                fprintf(outputFile, FormatString, value, error); 
            }
        }
        fclose(outputFile);
    } 
}


// Remove from TDirectory
void TopSVDFunctions::SVD_RmDir1D(TH1D* hists, int numHist)
{
    // Existence of Objects
    if ( hists == NULL ) return;
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        (hists+h)->SetDirectory(NULL);
    } 
}  

// Get Double From Array
double TopSVDFunctions::SVD_DoubleFromArray(double* arr, int pos)
{
    if ( arr == NULL ) return 0.;
    if ( pos < 0 ) return 0.;
    double returnvalue = *(arr+pos);
    return returnvalue;
}

// Remove from TDirectory
void TopSVDFunctions::SVD_RmDir2D(TH2D* hists, int numHist)
{
    // Existence of Objects
    if ( hists == NULL ) return;
    
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        (hists+h)->SetDirectory(NULL);
    } 
}  


// Find plain file name without path
TString TopSVDFunctions::SVD_FindFile(TString filepath)
{    
    
    
    // Find the folder the file belongs to 
    // and DO NOT Chop the last trailing slash
    TString outputpath = filepath;
    bool foundFolder = false; 
    bool foundSlash = false;
    while ( foundSlash == false ) {
        if ( outputpath.EndsWith("/") == true ) {
            foundSlash = true;
            foundFolder = true; 
        } else {
            foundSlash = false;
            outputpath.Chop(); 
        }
    }
    
    // if no path is given, return
    if ( foundFolder == false ) {
        return filepath;
    }
    
    // Otherwise, remove the directory string
    int lengthOfFolder = (int) outputpath.Length();
    TString filename = filepath;
    filename.Remove(0, lengthOfFolder);
    
    // return the truncated filename
    return filename;        
}


// Find parent folder of a file
// If just a filename is given, "." is returned
TString TopSVDFunctions::SVD_FindFolder(TString filepath)
{    
    
    // Find the folder the file belongs to 
    TString outputpath = filepath;
    bool foundFolder = false; 
    bool foundSlash = false;
    while ( foundSlash == false ) {
        if ( outputpath.EndsWith("/") == true ) {
            foundSlash = true;
            foundFolder = true;
            outputpath.Chop(); 
        } else {
            foundSlash = false;
            outputpath.Chop(); 
        }
    }
    

    // Return folder
    TString dot = ".";
    if ( foundFolder == true ) {
        return outputpath;
    } else {
        return dot;
    }
        
}


// Make Folder
void TopSVDFunctions::SVD_MakeFolder(TString outputpath)
{
    // Use the TSystem class for this 
    // Do it RECURSIVELY!!
    gSystem->mkdir(outputpath, true); 
//    bool success = gSystem->mkdir(outputpath, true) == 0; 
//    if (!success) {
//        std::cerr << "Cannot create directory " << outputpath << std::endl;
//        std::exit(1);
//    }
}


// Make Folder and a file
void TopSVDFunctions::SVD_TouchFile(TString filepath)
{
    
    // Create folder
    TString outputpath = SVD_FindFolder(filepath);
    SVD_MakeFolder(outputpath);
    
    // Find file
    TString filename = SVD_FindFile(filepath);
     
    // Touch the File
//     TString touchCommand = "cd ";
//     touchCommand.Append(outputpath);
//     touchCommand.Append(" ; touch ");
//     touchCommand.Append(filename);
//     gSystem->Exec(touchCommand); 
    
    //WTF?! Why would one call external programs?!
    std::ofstream touch((outputpath + filename).Data(), std::ios::app);
}


// Remove file
void TopSVDFunctions::SVD_RemoveFile(TString filepath)
{
    
    // Get the folder 
    TString outputpath = SVD_FindFolder(filepath);
    
    // Find file
    TString filename = SVD_FindFile(filepath);
    
    // Remove command
//     TString touchCommand = "cd ";
//     touchCommand.Append(outputpath);
//     touchCommand.Append(" ; rm ");
//     touchCommand.Append(filename);
//     touchCommand.Append(" 2> /dev/null");
//     gSystem->Exec(touchCommand); 
    std::remove((outputpath+filename).Data());
    
}


// Get a digit from an int
int TopSVDFunctions::SVD_GetDigit(TString steering, int digit, int standard)
{
    
    int neededLength = 17;
    TString oldsteering = steering;
    
    // Length
    int length = steering.Length();
    
    // First Plausibility Check
    if ( length < 0 || digit < 1 || digit > neededLength || standard < 1 || standard > 9) {
        cout << "**********************************************************************" << endl;
        cout << "Error in TopSVDFunctions::SVD_GetDigit()" << endl;
        cout << "These parameters do not make sense: " << endl;
        cout << "steering = " << oldsteering << endl;
        cout << "digit = " << digit << endl;
        cout << "standard = " << standard << endl;
        cout << "**********************************************************************" << endl;
        exit(1);
    }
    
    
     
    // Fill up Zeros
    while ( steering.Length() < neededLength ) {
        steering.Prepend("0");
    }
    length = steering.Length();
     
    
    // Remove digits from right
    while ( steering.Length() > neededLength+1-digit ) {
        steering.Chop();
    }
    length = steering.Length();
     
    
    // Remove digits from left
    while ( steering.Length() > 1 ) {
        steering.Remove(0,1);
    }
    length = steering.Length(); 
     
    // Get Number
    int number = steering.Atoi(); 

    
    // Second Plausibility Check
    if ( number < 0 || digit > neededLength || length < 1 || length  > 1 ) {
        cout << "**********************************************************************" << endl;
        cout << "Error in TopSVDFunctions::SVD_GetDigit()" << endl;
        cout << "These parameters do not make sense: " << endl;
        cout << "steering = " << oldsteering << endl;
        cout << "digit = " << digit << endl;
        cout << "standard = " << standard << endl;
        cout << "**********************************************************************" << endl;
        exit(1);
    }
 
    // Standard Value
    if ( number == 0 ) number = standard;
 
    // REturn it
    return number;
}     

// String to file
// Options: "a" is append, "w" is write
void TopSVDFunctions::SVD_LineToFile(TString string, TString filepath, TString option)
{ 
    
    // Find and make folder
    TString outputpath = SVD_FindFolder(filepath);
    SVD_MakeFolder(outputpath);
    
    // Touch the file ... this won't hurt at all 
    SVD_TouchFile(filepath);
     
    // Open File
    FILE* theFile;
    const char* c_option = option.Data();
    const char* c_path = filepath.Data();
    theFile = fopen(c_path,c_option);
    
    // Warning
    if ( theFile == NULL ) {
        cout << "***********************************************************************" << endl;
        cout << "Error in  TopSVDFunctions::SVD_LineToFile " << endl;
        cout << "Not able to open file " << filepath << endl;
        cout << "Check this now!" << endl;
        cout << "***********************************************************************" << endl;
        exit(1);
    }
    
    // Append newline character to string
    string.Append("\n");
    
    // Print to File 
    const char* c_string = string.Data();
    fputs(c_string, theFile);
     
    
    
    // Close file
    fclose(theFile);  
        
}

// line from file
// Search for a line that starts with startstring and return it
TString TopSVDFunctions::SVD_LineFromFile(TString key, TString filepath)
{ 
    // Read from file
    FILE* theFile;
    char linebuffer[256];
    theFile = fopen(filepath,"r");
    
    // Startstring
    // Append a whitespace to avoid being sensitive
    // to keys that start other keys
    TString startstring = key;
    startstring.Append(" ");
    
    // Did the file opening succeed?
    if ( theFile == NULL ) {
        cout << "***********************************************************************" << endl;
        cout << "Error in  TopSVDFunctions::SVD_LineFromFile() " << endl;
        cout << "Not able to open file " << filepath << endl;
        cout << "Check this now!" << endl;
        cout << "***********************************************************************" << endl;
        exit(1);
    }

    // Search for full line
    bool foundKey = false; 
    TString theFullLine = "";
     
    // Loop over all lines
    while ( fgets(linebuffer, 256, theFile) ) {

        // Get the line
        TString theLine(linebuffer);
        
        // Check the beginning
        if ( theLine.BeginsWith(startstring) == true ) {
            foundKey = true;
            theFullLine = theLine;
            break;
        }
    }
    
    // Return what was found or give an error message
    if ( foundKey == true ) {
        return theFullLine;
    } else {  
        cout << "***********************************************************************" << endl;
        cout << "Error in  TopSVDFunctions::SVD_LineFromFile() " << endl; 
        cout << "We were searching for the key" << endl;
        cout << "'" << startstring << "'" << endl;
        cout << "but didn't find it in the file " << endl;
        cout << filepath << " " << endl;
        cout << "Check this now!" << endl;
        cout << "***********************************************************************" << endl;
        exit(1); 
        return theFullLine;
    }
}


// Print Page
// Note: Empty filename means: Do not print to this file 
// Note: Options can be OPEN and CLOSE
void TopSVDFunctions::SVD_PrintPage(TCanvas*& canvas, TString outputfilenamePs, TString outputfilenameEps, TString description, TString option)
{
    
    // What to print?
    bool doEps = false;
    if ( outputfilenameEps.CompareTo("") != 0 ) doEps = true;
    bool doPs = false;
    if ( outputfilenamePs.CompareTo("") != 0 ) doPs = true;
    if ( doEps == false && doPs == false ) return;
    
    // Options
    bool doOpen = false;
    if ( option.CompareTo("OPEN") == 0 ) doOpen = true;
    bool doClose = false;
    if ( option.CompareTo("CLOSE") == 0 ) doClose = true;
     
    
    // Filename EPS
    TString fulloutputfilenameEps = outputfilenameEps;
    fulloutputfilenameEps.Chop();
    fulloutputfilenameEps.Chop();
    fulloutputfilenameEps.Chop();
    fulloutputfilenameEps.Chop(); 
    fulloutputfilenameEps.Append("_");
    fulloutputfilenameEps.Append(description);
    fulloutputfilenameEps.Append(".eps"); 
    
    
    // Filename C 
    TString fulloutputfilenameC = outputfilenameEps;
    fulloutputfilenameC.Chop();
    fulloutputfilenameC.Chop();
    fulloutputfilenameC.Chop();
    fulloutputfilenameC.Chop(); 
    fulloutputfilenameC.Append("_");
    fulloutputfilenameC.Append(description);
    fulloutputfilenameC.Append(".C"); 
     
     
    // Filename PS 
    TString fulloutputfilenamePs = outputfilenamePs;
    fulloutputfilenamePs.Chop();
    fulloutputfilenamePs.Chop();
    fulloutputfilenamePs.Chop();
    fulloutputfilenamePs.Append(".ps");
    if ( doOpen == true  ) fulloutputfilenamePs.Append("(");
    if ( doClose == true ) fulloutputfilenamePs.Append(")");
    
    
    // Print plot   
    if ( doEps == true ) {
        canvas->Print(fulloutputfilenameEps);
    //    canvas->Print(fulloutputfilenameC);
    }
    if ( doPs == true ) { 
        canvas->Print(fulloutputfilenamePs);   
    }
    
    
    
}


// Delete obsolete TopSVDUnfold objects
void TopSVDFunctions::SVD_DeleteSVD(TopSVDUnfold* SVDs, int numSVDs)
{
    if ( SVDs == NULL ) return;
    if ( numSVDs <= 0 ) return;
    
    delete[] SVDs;  
    SVDs = NULL;

}



// Normalization
// Calculate Global Cross Section
void TopSVDFunctions::SVD_GlobalCrossSection(double*& globCrossSection, double*& globCrossSectionErr, double* globEventYield, double* globEventYieldErr, double* globEfficiency, int numHist)
{
     
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Get Input
        double globalEvYield    = SVD_DoubleFromArray(globEventYield, h);
        double globalEvYieldErr = SVD_DoubleFromArray(globEventYieldErr, h);
        double globalEfficiency = SVD_DoubleFromArray(globEfficiency, h);
        
        
        // Calculate Cross Section and Uncertainty
        double globalXSec    = SVD_Divide(globalEvYield,    globalEfficiency);
        double globalXSecErr = SVD_Divide(globalEvYieldErr, globalEfficiency);
     
        // Save it
        *(globCrossSection+h) = globalXSec;
        *(globCrossSectionErr+h) = globalXSecErr;
    }
    
}
    
// Normalization
// Calculate Global Event Yield
void TopSVDFunctions::SVD_GlobalEventYield(double*& globEvYield, double*& globEvYieldErr, double* totalDataEvents, double* totalBgrEvents, double* totalTtBgrEvents, double* totalRecEvents, int numHist)
{
     
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Get Input
        // For data, use only first entry in array
        double data = SVD_DoubleFromArray(totalDataEvents, 0); 
        double bgr = SVD_DoubleFromArray(totalBgrEvents, h);
        double ttbgr = SVD_DoubleFromArray(totalTtBgrEvents, h);
        double signal =  SVD_DoubleFromArray(totalRecEvents, h);
        
        // Calculate Yield
        double signalFraction = signal / (signal + ttbgr);
        double yield = (data - (bgr-ttbgr)) * signalFraction;
        
        // Input Error
        double dataErrSq = data;
        double bgrErrSq = 0.;
        double ttbgrErrSq = 0.;
        
        // Calculate Yield Error
        double yieldErrSq = dataErrSq + bgrErrSq + ttbgrErrSq;
        double yieldErr = TMath::Sqrt(yieldErrSq);
         
        // Save it
        *(globEvYield+h) = yield;
        *(globEvYieldErr+h) = yieldErr;
        
    } 
}

// Normalization
// Normalize SVD unfolded Event counts, using a global Efficiency and a global event yield
// This creates new histograms on the heap. Do not forget to delete them sometimes.
void TopSVDFunctions::SVD_ExtNormalizeSVDDistribution(TH1D* outputResultHist, TH2D*& outputCovHist, TH1D* inputHist, TH2D* probMatrixHist, TH2D* statCovMatrix, double* globalEfficiency, double* globalEventYield, double* globalEventYieldErr, int numHist, int verbose)
{ 
    // Thomas
    if (0) std::cout << probMatrixHist << statCovMatrix << verbose << std::endl; // Construct to avoid compilation warning
     
    // Existence of Objects
    if ( inputHist == NULL ) return;  
    if ( outputResultHist == NULL ) return;  
    if ( outputCovHist == NULL ) return;  
    
    // Number of Bins
    int nbins = inputHist->GetNbinsX();
    
    // Loop over distributions
    for ( int h = 0 ; h < numHist ; h++ ) { 
        
        double globalYield = *(globalEventYield+h);
        double globalYieldErr = *(globalEventYieldErr+h);
        double globalEff = *(globalEfficiency+h);  
        
        // Calculate the Covariance Matrix
        // only for nominal sample
        if ( h == 0 ) {
              
            // Loop over bins, including OF
            for ( int i = 1 ; i <= nbins ; i++ ) { 
                for ( int j = 1 ; j<= nbins ; j++ ) {
                
                    // First term
                    double cov_first = (statCovMatrix)->GetBinContent(i,j);
                    
                    
                    // Second term
                    double cov_second = 0.;
                    for ( int a = 1 ; a <= nbins ; a++ ) { 
                        for ( int b = 1 ; b<= nbins ; b++ ) {
                            double response_term = 0.01*probMatrixHist->GetBinContent(a,b);
                            double cov_term = (statCovMatrix)->GetBinContent(i,b); 
                            cov_second += response_term*cov_term;  
                        }
                    }
                    cov_second = cov_second * (inputHist)->GetBinContent(j);
                    cov_second = cov_second * (-1.) / globalYield;
                    
                    
                    // Third term
                    double cov_third = 0.;
                    for ( int a = 1 ; a <= nbins ; a++ ) { 
                        for ( int b = 1 ; b<= nbins ; b++ ) {
                            double response_term = 0.01*probMatrixHist->GetBinContent(a,b);
                            double cov_term = (statCovMatrix)->GetBinContent(j,b); 
                            cov_third += response_term*cov_term; 
                        }
                    }
                    cov_third = cov_third * (inputHist)->GetBinContent(i);
                    cov_third = cov_third * (-1.) / globalYield;
                    
                    // Fourth term
                    double cov_fourth = globalYieldErr * globalYieldErr / globalYield / globalYield;
                    cov_fourth = cov_fourth * (inputHist)->GetBinContent(i);
                    cov_fourth = cov_fourth * (inputHist)->GetBinContent(j);  

                    
                    // Sum up and normalize
                    double cov_new = cov_first + cov_second + cov_third + cov_fourth;
                    cov_new = cov_new * (globalEff / globalYield) * (globalEff / globalYield); 
                    (outputCovHist)->SetBinContent(i, j, cov_new);
                    (outputCovHist)->SetBinError(i, j, 0.);
                }
            }
        }
        
        
        // Loop over bins, including OF
        for ( int i = 1 ; i <= nbins ; i++ ) {
               
            double value_old = (inputHist+h)->GetBinContent(i); 
                
            double value_new = value_old * globalEff / globalYield; 
            
            // Get Error from Covariance Matrix
            double error_new = TMath::Sqrt(outputCovHist->GetBinContent(i, i)); 
           
            (outputResultHist+h)->SetBinContent(i, value_new);
            (outputResultHist+h)->SetBinError(i, error_new);
        }
        
    }  
}

// Normalization
// Normalize BBB unfolded Event counts, using a global Efficiency and a global event yield
// This creates new histograms on the heap. Do not forget to delete them sometimes.
void TopSVDFunctions::SVD_ExtNormalizeBBBDistribution(TH1D* outputResultHist, TH1D* inputHist, double* globalEfficiency, double* globalEventYield, double* globalEventYieldErr, int numHist, int verbose)
{
    
    // Thomas  
    if (0) std::cout << globalEventYieldErr << verbose << std::endl; // Construct to avoid compilation warning
     
    // Existence of Objects
    if ( inputHist == NULL ) return;   
    if ( outputResultHist == NULL ) return;   
    
    // Number of Bins
    int nbins = inputHist->GetNbinsX();
    
    // Loop over distributions
    for ( int h = 0 ; h < numHist ; h++ ) { 
        
        double globalYield = *(globalEventYield+h);
        double globalEff = *(globalEfficiency+h);  
        
        // Loop over bins, including OF
        for ( int i = 1 ; i <= nbins ; i++ ) {
               
            double value_old = (inputHist+h)->GetBinContent(i);
            double error_old = (inputHist+h)->GetBinError(i);
                
            double value_new = value_old * globalEff / globalYield; 
            double error_new = error_old * globalEff / globalYield;
           
            (outputResultHist+h)->SetBinContent(i, value_new);
            (outputResultHist+h)->SetBinError(i, error_new);
        }
    }   
    
}

// Normalization
// Normalize Generator Level distribution
// This creates new histograms on the heap. Do not forget to delete them sometimes.
void TopSVDFunctions::SVD_ExtNormalizeGenDistribution(TH1D* outputResultHist, TH1D* inputHist, double* totalGenEvents, int numHist, int verbose)
{
    
    // Thomas  
    if (0) std::cout <<  verbose << std::endl; // Construct to avoid compilation warning
     
    // Existence of Objects
    if ( inputHist == NULL ) return;    
    if ( outputResultHist == NULL ) return;   
    
    // Number of Bins
    int nbins = inputHist->GetNbinsX();
    
    // Loop over distributions
    for ( int h = 0 ; h < numHist ; h++ ) {  
        
        double totalEvents = *(totalGenEvents+h);
        
        // Loop over bins, including OF
        for ( int i = 1 ; i <= nbins ; i++ ) {
               
            double value_old = (inputHist+h)->GetBinContent(i); 
               
            double value_new = value_old / totalEvents;
            double error_new = 0.;
           
            (outputResultHist+h)->SetBinContent(i, value_new);
            (outputResultHist+h)->SetBinError(i, error_new);
        }
    }   
    
}

// Normalization
// Normalize SVD unfolded Event counts, using a global Efficiency and a global event yield
// This creates new histograms on the heap. Do not forget to delete them sometimes.
void TopSVDFunctions::SVD_IntNormalizeSVDDistribution(TH1D* outputResultHist, TH2D*& outputCovHist, TH1D* inputHist, TH2D* statCovMatrix, int numHist, int verbose)
{ 
    
    // Thomas  
    if (0) std::cout << statCovMatrix << verbose << std::endl; // Construct to avoid compilation warning
 
    // Existence of Objects
    if ( inputHist == NULL ) return;   
    if ( outputResultHist == NULL ) return;  
    if ( outputCovHist == NULL ) return;  
    
    // Number of Bins
    int nbins = inputHist->GetNbinsX();
    
    // Loop over distributions
    for ( int h = 0 ; h < numHist ; h++ ) { 
        
        // Get Integral 
        bool doOF = true;
        double integral = SVD_Integral1D(inputHist, h, doOF);  
        
        // Calculate the Covariance Matrix
        // only for nominal sample
        if ( h == 0 ) {
              
            // Loop over bins, including OF
            for ( int i = 1 ; i <= nbins ; i++ ) { 
                for ( int j = 1 ; j<= nbins ; j++ ) {
                
                    double cov_old = (statCovMatrix)->GetBinContent(i,j);
                    
                    // First Term
                    double cov_first = integral * integral * cov_old;
                    
                    // Second Term
                    double cov_second = 0.;
                    for ( int a = 1 ; a <= nbins; a++ ) {
                         cov_second += (statCovMatrix)->GetBinContent(a,i);
                    }
                    cov_second = cov_second * (inputHist)->GetBinContent(j);
                    cov_second = (-1.) * cov_second * integral ; 
                    
                    // Third Term
                    double cov_third = 0.;
                    for ( int a = 1 ; a <= nbins; a++ ) {
                         cov_third += (statCovMatrix)->GetBinContent(a,j);
                    }
                    cov_third = cov_third * (inputHist)->GetBinContent(i);
                    cov_third = (-1.) * cov_third * integral ; 
                    
                    // Fourth Term
                    double cov_fourth = 0.;
                    for ( int a = 1 ; a <= nbins; a++ ) { 
                        for ( int b = 1 ; b <= nbins; b++ ) {
                             cov_fourth += (statCovMatrix)->GetBinContent(a,b);
                        }
                    }
                    cov_fourth = cov_fourth * (inputHist)->GetBinContent(i);
                    cov_fourth = cov_fourth * (inputHist)->GetBinContent(j);
                    
                    // Normalization
                    double cov_new = cov_first + cov_second + cov_third + cov_fourth; 
                    if ( integral > 0. ) cov_new = cov_new / integral / integral / integral / integral;
                    
                    (outputCovHist)->SetBinContent(i, j, cov_new);
                    (outputCovHist)->SetBinError(i, j, 0.);
                }
            }
        }
        
        
        // Loop over bins, including OF
        for ( int i = 1 ; i <= nbins ; i++ ) {
                
            double value_old = (inputHist+h)->GetBinContent(i); 
                
            double value_new = value_old;
            if ( integral > 0. ) value_new = value_new / integral ;  
            
            // Get Error from Covariance Matrix
            double error_new = TMath::Sqrt(outputCovHist->GetBinContent(i, i)); 
           
            (outputResultHist+h)->SetBinContent(i, value_new);             
            (outputResultHist+h)->SetBinError(i, error_new);
        }
        
        
    }   
}

// Normalization
// Normalize BBB unfolded Event counts, using a global Efficiency and a global event yield
// This creates new histograms on the heap. Do not forget to delete them sometimes.
void TopSVDFunctions::SVD_IntNormalizeBBBDistribution(TH1D* outputResultHist, TH1D* inputHist, int numHist, int verbose)
{
    
    
    // Thomas  
    if (0) std::cout <<  verbose << std::endl; // Construct to avoid compilation warning
     
    // Existence of Objects
    if ( inputHist == NULL ) return;   
    if ( outputResultHist == NULL ) return;   
    
    // Number of Bins
    int nbins = inputHist->GetNbinsX();
    
    // Loop over distributions
    for ( int h = 0 ; h < numHist ; h++ ) { 
        
        // Get Integral 
        bool doOF = true;
        double integral = SVD_Integral1D(inputHist, h, doOF);  
        
        // Loop over bins, including OF
        for ( int i = 1 ; i <= nbins ; i++ ) {
               
            double value_old = (inputHist+h)->GetBinContent(i);
            double error_old = (inputHist+h)->GetBinError(i);
                
            double value_new = value_old;
            if ( integral > 0. ) value_new = value_new / integral ;  
            double error_new = error_old;
            if ( integral > 0. ) error_new = error_new / integral ;  
           
            (outputResultHist+h)->SetBinContent(i, value_new);
            (outputResultHist+h)->SetBinError(i, error_new);
        }
    }  
    
}

// Normalization
// Normalize Generator Level distribution
// This creates new histograms on the heap. Do not forget to delete them sometimes.
void TopSVDFunctions::SVD_IntNormalizeGenDistribution(TH1D* outputResultHist, TH1D* inputHist, int numHist, int verbose)
{
     
    
    // Thomas  
    if (0) std::cout <<  verbose << std::endl; // Construct to avoid compilation warning
    
    
    // Existence of Objects
    if ( inputHist == NULL ) return;  
    if ( outputResultHist == NULL ) return;   
    
    // Number of Bins
    int nbins = inputHist->GetNbinsX();
    
    // Loop over distributions
    for ( int h = 0 ; h < numHist ; h++ ) { 
        
        // Get Integral 
        bool doOF = true;
        double integral = SVD_Integral1D(inputHist, h, doOF);  
        
        // Loop over bins, including OF
        for ( int i = 1 ; i <= nbins ; i++ ) {
               
            double value_old = (inputHist+h)->GetBinContent(i);
            double error_old = (inputHist+h)->GetBinError(i);
                
            double value_new = value_old;
            if ( integral > 0. ) value_new = value_new / integral ;  
            double error_new = error_old;
            if ( integral > 0. ) error_new = error_new / integral ;  
           
            (outputResultHist+h)->SetBinContent(i, value_new);
            (outputResultHist+h)->SetBinError(i, error_new);
        }
    }  
    
}


// Normalization
// Normalize SVD unfolded Event counts, using a global Efficiency and a global event yield
// This creates new histograms on the heap. Do not forget to delete them sometimes.
void TopSVDFunctions::SVD_PseudatNormalizeSVDDistribution(TH1D* outputResultHist, TH2D*& outputCovHist, TH1D* inputHist, TH2D* statCovMatrixNorm, int numHist, int verbose)
{ 
    // Thomas  
    if (0) std::cout << statCovMatrixNorm << verbose << std::endl; // Construct to avoid compilation warning
 
    // Existence of Objects
    if ( inputHist == NULL ) return;   
    if ( outputResultHist == NULL ) return;  
    if ( outputCovHist == NULL ) return;  
    
    // Number of Bins
    int nbins = inputHist->GetNbinsX(); 
    
    // Loop over distributions
    for ( int h = 0 ; h < numHist ; h++ ) { 
        
        // Get Integral 
        bool doOF = true;
        double integral = SVD_Integral1D(inputHist, h, doOF);  
    
        //---
        // parameters for mathematical calculation of error by error propagation
        double errorIntegral2 = 0.;
        double errorIntegralNaive2 = 0.;
        for ( int i = 1 ; i <= nbins ; i++ ) {
            for ( int j = 1 ; j <= nbins ; j++ ) {
                errorIntegral2 += statCovMatrixNorm->GetBinContent(i,j);
            }
            errorIntegralNaive2 += statCovMatrixNorm->GetBinContent(i,i);
        }
        if(verbose>=2) cout << "errorIntegral2 = " << errorIntegral2 << "; errorIntegralNaive2= " << errorIntegralNaive2 <<endl;
        //---
        
        // Loop over bins, including OF
        for ( int i = 1 ; i <= nbins ; i++ ) {
               
            double value_old = (inputHist+h)->GetBinContent(i);
            double error_old = (inputHist+h)->GetBinError(i);
                    
            double value_new = value_old;
            if ( integral > 0. ) value_new = value_new / integral ;  
            double errorFromDiv = error_old;
            double errorFromCov = 0.;
            if ( integral > 0. )   errorFromDiv = errorFromDiv / integral ;
            if (statCovMatrixNorm) errorFromCov = TMath::Sqrt(statCovMatrixNorm->GetBinContent(i,i));
              
            //---
            // mathematical calculation of error by error propagation
            double covIntegralSigmaI = 0.;
            for ( int j = 1 ; j <= nbins ; j++ ) {
                covIntegralSigmaI += statCovMatrixNorm->GetBinContent(j,i);
            }
            // derived equation
            double errorFromCovCalculated =  value_new*value_new * TMath::Sqrt(errorIntegral2/(value_old*value_old) - 2*integral*covIntegralSigmaI/TMath::Power(value_old,3) + integral*integral *error_old*error_old/TMath::Power(value_old,4)  );
            //---
            
            (outputResultHist)->SetBinContent(i, value_new);
            (outputResultHist)->SetBinError(i, errorFromCov);
            
            
            // Copy covariance matrix
            if ( h == 0 ) {
                for ( int k = 1 ; k <= nbins ; k++ ) { 
                    for ( int l = 1 ; l<= nbins ; l++ ) {  
                        outputCovHist->SetBinContent(k,l, statCovMatrixNorm->GetBinContent(k,l));
                    }
                }
            }
            
            // comparison of different errors in verbose case
            if(verbose>=2){
                double value = value_new;
                if(value==0) value=1e10;
                double relDifference = 0.;
                if(errorFromCov!=0) relDifference = (errorFromCov-errorFromDiv)/errorFromCov;
                cout << "rel. err. from division;            bin " << i << "; err/value= " << errorFromDiv/value <<endl;
                cout << "rel. err. from covariance matrix;   bin " << i << "; err/value= " << errorFromCov/value <<endl;
                cout << "rel. err. from covariance m. calc.; bin " << i << "; err/value= " << errorFromCovCalculated/value <<endl;
                cout << "difference = " << relDifference <<endl;
            }  
    }
    }  
}

// Normalization
// Normalize BBB unfolded Event counts, using a global Efficiency and a global event yield
// This creates new histograms on the heap. Do not forget to delete them sometimes.
void TopSVDFunctions::SVD_PseudatNormalizeBBBDistribution(TH1D* outputResultHist, TH1D* inputHist, int numHist, int verbose)
{
    
    
    // Thomas  
    if (0) std::cout <<  verbose << std::endl; // Construct to avoid compilation warning
     
    // Existence of Objects
    if ( inputHist == NULL ) return;   
    if ( outputResultHist == NULL ) return;   
    
    // Number of Bins
    int nbins = inputHist->GetNbinsX();
    
    // Loop over distributions
    for ( int h = 0 ; h < numHist ; h++ ) { 
        
        // Get Integral 
        bool doOF = true;
        double integral = SVD_Integral1D(inputHist, h, doOF);  
        
        // Loop over bins, including OF
        for ( int i = 1 ; i <= nbins ; i++ ) {
               
            double value_old = (inputHist+h)->GetBinContent(i);
            double error_old = (inputHist+h)->GetBinError(i);
                
            double value_new = value_old;
            if ( integral > 0. ) value_new = value_new / integral ;  
            double error_new = error_old;
            if ( integral > 0. ) error_new = error_new / integral ;  
           
            (outputResultHist+h)->SetBinContent(i, value_new);
            (outputResultHist+h)->SetBinError(i, error_new);
        }
    }  
    
}

// Normalization
// Normalize Generator Level distribution
// This creates new histograms on the heap. Do not forget to delete them sometimes.
void TopSVDFunctions::SVD_PseudatNormalizeGenDistribution(TH1D* outputResultHist, TH1D* inputHist, int numHist, int verbose)
{
     
    
    // Thomas  
    if (0) std::cout <<  verbose << std::endl; // Construct to avoid compilation warning
    
    
    // Existence of Objects
    if ( inputHist == NULL ) return;  
    if ( outputResultHist == NULL ) return;   
    
    // Number of Bins
    int nbins = inputHist->GetNbinsX();
    
    // Loop over distributions
    for ( int h = 0 ; h < numHist ; h++ ) { 
        
        // Get Integral 
        bool doOF = true;
        double integral = SVD_Integral1D(inputHist, h, doOF);  
        
        // Loop over bins, including OF
        for ( int i = 1 ; i <= nbins ; i++ ) {
               
            double value_old = (inputHist+h)->GetBinContent(i);
            double error_old = (inputHist+h)->GetBinError(i);
                
            double value_new = value_old;
            if ( integral > 0. ) value_new = value_new / integral ;  
            double error_new = error_old;
            if ( integral > 0. ) error_new = error_new / integral ;  
           
            (outputResultHist+h)->SetBinContent(i, value_new);
            (outputResultHist+h)->SetBinError(i, error_new);
        }
    }  
    
}
// Normalization
// Calculate Globel Efficiency
void TopSVDFunctions::SVD_GlobalEfficiency(double*& globalEff, double* totalRecEvents, double* totalGenEvents, int numHist)
{
     
    // Loop over all histograms
    for ( int h = 0 ; h < numHist ; h++ ) {
        
        // Get Input
        double rec = SVD_DoubleFromArray(totalRecEvents, h);
        double gen = SVD_DoubleFromArray(totalGenEvents, h); 
        
        
        // Zero Gen?
        if ( gen <= 0.) {
            cout << "Error in TopSVDFUnctions::SVD_GlobalEfficiency" << endl;
            cout << "No generator level events found for systematic " << h << endl;
            cout << "Exiting ... " << endl;
            exit(1);
        }
        
        // Calculate Efficiency
        double eff = rec / gen ;
         
        // Save it
        *(globalEff+h) = eff;
        
    }
    
}

// Count Bins of a histogram 
int TopSVDFunctions::SVD_NumBins1D(TH1D* histo, bool doCntOF)
{
    if ( histo == NULL ) return 0;
    int numbins = histo->GetXaxis()->GetNbins();
    if ( doCntOF == true ) numbins += 2;
    return numbins;
}

// Count Bins of a histogram 
int TopSVDFunctions::SVD_NumBinsX2D(TH2D* histo, bool doCntOF)
{
    if ( histo == NULL ) return 0;
    int numbins = histo->GetXaxis()->GetNbins();
    if ( doCntOF == true ) numbins += 2;
    return numbins;
}

// Count Bins of a histogram 
int TopSVDFunctions::SVD_NumBinsY2D(TH2D* histo, bool doCntOF)
{
    if ( histo == NULL ) return 0;
    int numbins = histo->GetYaxis()->GetNbins();
    if ( doCntOF == true ) numbins += 2;
    return numbins;
}

// Save Bin Boundaries in Array
void TopSVDFunctions::SVD_GetBinBoundaries1D(double* boundaries, TH1D* histo)
{
    if ( boundaries == NULL ) return;
    if ( histo == NULL ) return;
    
    int numBnd = SVD_NumBins1D(histo, false) + 1;
    for ( int i = 0 ; i < numBnd ; i++ ) {
        boundaries[i] = histo->GetXaxis()->GetBinLowEdge(i+1);
    }
}

// Save Bin Boundaries in Array
void TopSVDFunctions::SVD_GetBinBoundariesX2D(double* boundaries, TH2D* histo)
{
    if ( boundaries == NULL ) return;
    if ( histo == NULL ) return;
    
    int numBnd = SVD_NumBinsX2D(histo, false) + 1;
    for ( int i = 0 ; i < numBnd ; i++ ) {
        boundaries[i] = histo->GetXaxis()->GetBinLowEdge(i+1);
    }
}

// Save Bin Boundaries in Array
void TopSVDFunctions::SVD_GetBinBoundariesY2D(double* boundaries, TH2D* histo)
{
    if ( boundaries == NULL ) return;
    if ( histo == NULL ) return;
    
    int numBnd = SVD_NumBinsY2D(histo, false) + 1;
    for ( int i = 0 ; i < numBnd ; i++ ) {
        boundaries[i] = histo->GetYaxis()->GetBinLowEdge(i+1);
    }
}

// PERFORM UNFOLDING 
// Steering options in parameter 'steering' 
//    (1)  REGMODE (1. digit from right)
//         0 means: Default setting. Same as 2
//         1 means: Bin by Bin Unfolding
//         2 means: SVD Unfolding, minimal regularization.
//              Regularization will be done with a k value
//              and the k value will be set to the number of
//              bins (not counting OF). Parameter 'regPar'
//              is ignored.
//         3 means: SVD unfolding, regularization according
//              to the parameter 'regPar'. Whether this one
//              is taken as k value or tau value, can be 
//              specified with the second digit of 'steering'
//              (default)
//         4 means: SVD unfolding according to the file
//              given in 'regParFile'. Whether this one
//              is taken as k value or tau value, can be 
//              specified with the second digit of this
//              parameter.Parameter 'regPar'
//              is ignored.
//    (2)  REGULARIZATION PARAMETER (2. digit from right)
//         0 means: Default setting, same as 1
//         1 means: interpret 'regPar' as k value (default)
//         2 means: interpret 'regPar' as tau value
//    (3)  SCAN ( 3. digit from right)
//         0 means: Default setting, same as 2
//         1 means: no scan for optimal tau is perforemd,
//         2 means: tau scan  (default).
//              A scan for the optimal tau parameter is
//              performed. The scan is performed around
//              a "center value" to be specified in 'regpar'.
//              Whether this one is taken
//              as k value or tau value, can be 
//              specified with the second digit of this
//              parameter.
//              Note: The scan may take a while!    
//    (4)  PS/EPS PLOTTING  (4. digit from right)
//         0 means: Default value, same as 4
//         1 means: no plotting at all
//         2 means: standard plots
//         3 means: standard plots + k scan plots
//         4 means: standard plots + k scan plots + tau scan plots (default)
//         Note: EPS are PS printed if their corresponding fileNames are specified.
//    (5)  ROOT FILE ( 5. digit from right)
//         0 means: Default value, same as 1
//         1 means: no root file will be written (default)
//         2 means: standard plots to root file  
//    (6)  TEXT FILE ( 6. digit from right)
//         0 means: Default value, same as 2
//         1 means: no text file is written
//         2 means: text file with histogram entries is written (default)
//         3 means: text file is written if (syst==sysNo)
//    (7)  VERBOSITY (7. digit from right)
//         0 means: Default value, same as 2
//         1 means: no output at all
//         2 means: standard output (default)
//         3 means: debugging output
//    (8)  SCANPOINTS (8. digit from right)
//         0 means: Default value, same as 3
//         1 means: 5 scan points
//         2 means: 25 scan points
//         3 means: 125 scan points (default)
//         4 means: 625 scan points
//    (9)  SCANRANGE (9. digit from right)
//         0 means: Default value, same as 2
//         1 means: Tau+ / Tau- = 100  
//         2 means: Tau+ / Tau- = 10000
//         3 means: Tau+ / Tau- = 1000000
//    (10) LOWER SIDE BIN (10. digit from right)
//         0 means: Default value, same as 3
//         1 means: Regard as regular bin (not encouraged!)
//         2 means: Regard as regular bin, keep bin content fixed to MC (not encouraged!)
//         3 means: Cut away on Rec Level, unfold to Gen Level (default)
//         4 means: Cut away on Rec Level, unfold to Gen Level, keep bin content fixed to MC 
//         5 means: Cut away on Rec Level, ignore on Gen Level (not encouraged!) 
//    (11) UPPER SIDE BIN (11. digit from right)
//         0 means: Default value, same as 3
//         1 means: Regard as regular bin (not encouraged!)
//         2 means: Regard as regular bin, keep bin content fixed to MC (not encouraged!)
//         3 means: Cut away on Rec Level, unfold to Gen Level (default)
//         4 means: Cut away on Rec Level, unfold to Gen Level, keep bin content fixed to MC 
//         5 means: Cut away on Rec Level, ignore on Gen Level (not encouraged!) 
//    (12) ORIENTATION OF RESPONSE MATRIX (12. digit from right)
//         0 means: Default value, same as 2
//         1 means: Do not transpose input response matrix during rebinning
//         2 means: Do transpose input response matrix during rebinning (default)
//    (13) NORMALIZATION (13. digit from right)
//         0 means: Default value, same as 1
//         1 means: Extrinsic Normalization (default)
//             This means, the normalization is done with a global inclusive cross
//             section which is calculated from the parameters totalDataEvents, 
//             totalBgrEvents, totalTtBgrEvents, totalRecEvents and totalGenEvents.
//             Analytic Error Propagation.
//         2 means: Intrinsic Normalization.
//             Each unfolded distribution is normalized with its integral. 
//             Analytic Error Propagation.
//         3 means: Intrinsic Normalization.
//             Each unfolded distribution is normalized with its integral. 
//             Error Propagation on the basis of pseudo experiments
//    (14) CLOSURE TEST (14. digit from right)
//         0 means: Default value, same as 1
//         1 means: No closure test, unfold real data (default)
//         2 means: Closure test, unfold the reconstructed MC 
//    (15) PRE-WEIGHTING (15. digit from right)
//         0 means: Default value, same as 1
//         1 means: No preweighting of MC is performed (default)
//         2 means: MC is reweighted to unfolded data (1 iteration)
//         3 means: MC is reweighted to unfolded data (2 iteration)
//         4 means: MC is reweighted to unfolded data (3 iteration)
//         5 means: MC is reweighted to unfolded data (4 iteration)
//         6 means: MC is reweighted to unfolded data (5 iteration)
//         7 means: MC is reweighted to unfolded data (6 iteration)
//         8 means: MC is reweighted to unfolded data (7 iteration)
//         9 means: MC is reweighted to unfolded data (8 iteration)
//    (16) BACKGROUND REDUCTION (16. digit from right)
//         0 means: Default value, same as 1
//         1 means: Be non-forgiving
//         2 means: If Background>Data, set Data to zero.
//    (17) CONTROL PLOT STYLE (17. digit from right)
//         0 means: Default value, same as 1
//         1 means: default style (e.g. show all bins incl. UF/OF, show vertical lines for tau scan etc.)
//         2 means: Uni HH style (e.g. show only bins of measurement, no vertical lines etc.)
//
// Return value: 
//        Best value of tau if scan is performed, -1. otherwise
// Systematics Handling: 
//        If numSys > 0, the TH1-parameters 
//        'respInputHist', 'recInputHist', 'bgrInputHist',  
//        'ttbgrInputHist' and 'genInputHist' will be interpreted as 
//        arrays with numSys+1 entries, with nominal sample at position
//        0 followed by the shifted samples.
//        Exception: Since you do not want to shift the data distribution,
//        this one will NOT be considered as an array!
//        The placeholder parameter 'unfolded' must also be an
//        array with numSys+1 elements. 
double TopSVDFunctions::SVD_Unfold(
        // Data Input (RAW Data including the background)
        TH1D* dataInputHist,   
        // Background (reco level, all BG, will be substracted from data)
        TH1D* bgrInputHist,
        // ttbar background only (used to calculate a ttbar signal 
        // fraction instead of subtracting the yield which depends 
        // on the inclusive ttbar cross section.) 
        // Note: if 0 pointer is handed over 
        TH1D* ttbgrInputHist, 
        // Generated signal MC
        TH1D* genInputHist,
        // Reconstructed signal MC
        TH1D* recInputHist,        
        // Response Matrix 
        TH2D* respInputHist,    
        // Returned: Unfolded Distribution.    
        TH1D*& unfolded,    
        // Returned: Normalized unfolded Distribution.    
        TH1D*& unfoldedNorm,    
        // For Normalization: Total Data Events
        // If set to zero, will be taken from the integral of dataInputHist
        double* totalDataEvents,
        // For Normalization: Total Reconstructed MC Events
        // If set to NULL, will be taken from the integral of bgrInputHist
        double* totalBgrEvents,
        // For Normalization: Total Reconstructed MC Events
        // If set to NULL, will be taken from the integral of ttbgrInputHist
        double* totalTtBgrEvents,
        // For Normalization: Total Reconstructed MC Events
        // If set to NULL, will be taken from the integral of recInputHist
        double* totalRecEvents,
        // For Normalization: Total Generated MC Events
        // If set to NULL, will be taken from the integral of genInputHist
        double* totalGenEvents,
        // Binning for the unfolding (bin sizes !!!)
        const double thebins[],     
        // Number of bins for unfolding (not counting OF bins !!!)
        const int numbins,  
        // Regularization parameter
        double regPar,                          
        // Regularization Modus
        TString steering,  
        // Specify the number of systematic samples to unfold 
        const int numSys, 
        // Specify Name for the Channel ("mumu", "emu", "ee" ...)
        TString channel,  
        // Specify Name for the Physics Object ("Top", "Jets", "Leptons")
        TString particle, 
        // Specify Name for the Quantity ("Eta", "Pt", or "Mass");
        TString quantity,
        // Specify Name for special run of unfolding
        TString special,
        // Specify Name of systematic
        // If you run over the nominal sample, provide NULL
        TString syst,
        // Nicely formatted name for the channel
        TString channelTex,
        // Nicely formatted name for the physics object
        TString particleTex,
        // Nicely formatted name for ithe quantity
        TString quantityTex,
        // Nicely formatted name indicating some special condition 
        TString specialTex,
        // Array of Names for the different systematics
        // if you run only over the nominal sample, provide NULL
        TString systTex,
        // If specified, plots will be saved in ROOT File
        // Please give filename WITHOUT extension, i.e.
        // "/asdf/asdf/file" instead of "/asdf/asdf/file.root" 
        TString rootFile,
        // If specified, plots will be saved in PS File 
        // Note: If EPS files are produces, the filename will be 
        // extended by the number of the plot. 
        TString psFile,
        // If specified, plots will be saved in EPS Files 
        // The filename will be extended by the number of the plot
        // and the file extension. 
        TString epsFile,
        // The optimal Reg Parameters will be written to this file.
        // The file will NOT be overwritten, rather the result will be appended.
        // The following data will be saved in this order: 
        // (1) the optimal tau, (2) the two nearest k values,
        // (3) the k value from the d value method
        // (4) the number of bins including side bins 
        TString txtFile,
        // In the following file, regularization parameters can be specified
        // for the unfolding. It may contain the optimal results from a previous run.
        // It should have the same file structure as txtFile, but it will be READ
        // rather than written. 
        // Note: It will only be read if the steering option 'REGMODE' is set to 4.
        TString regParFile
)
{  

 
    // Number of Pseudo Experiments for the error calculation    
    int nExperiments = 1000;
     
   
      
    ///////////////////////////////////////////////////////////////////  
    /////////  S E T   F L A G S   //////////////////////////////////// 
    /////////////////////////////////////////////////////////////////// 

     
    // REG MODE
    int flag_regmode = SVD_GetDigit(steering, 1, 2); 
      
      
    // REG PAR
    int flag_regpar = SVD_GetDigit(steering, 2, 1);  
     
     
    // SCAN
    int flag_scan = SVD_GetDigit(steering, 3, 2); 
     
     
    // PS FILE
    int flag_ps = SVD_GetDigit(steering, 4, 4);  
    if ( flag_scan == 1 && flag_ps == 4 ) flag_ps = 3;
    if ( psFile.CompareTo("")==0 ) flag_ps = 1;
     
     
    // ROOT FILE
    int flag_root = SVD_GetDigit(steering, 5, 1); 
    if ( rootFile.CompareTo("")==0 ) flag_root = 1;
      
     
    // TEXT FILE
    int flag_text = SVD_GetDigit(steering, 6, 2);  
    bool writeHistToText=false;
    if(flag_text==2) writeHistToText=true;
    if(flag_text==3 && syst=="sysNo") writeHistToText=true;


    // VERBOSITY
    int flag_verbose = SVD_GetDigit(steering, 7, 2);    
 
 
    // SCAN POINTS
    int flag_scanpoints = SVD_GetDigit(steering, 8, 3); 
    int nScanPoints = 1;
    for ( int i = 1 ; i <= flag_scanpoints ; i++ ) {
        nScanPoints = nScanPoints * 5;
    } 
    
    
    // SCAN RANGE
    int flag_range = SVD_GetDigit(steering, 9, 2); 
    double rangefactor = 1.;
    for ( int i = 1 ; i <= flag_range ; i++ ) {
        rangefactor = rangefactor * 10.;
    } 
     
     
     
    // LOWER SIDE BIN
    int flag_lowerOF = SVD_GetDigit(steering, 10, 3); 
    bool cutLowerRecSideBin = (flag_lowerOF >= 3); 
    bool fixLowerSideBin = (flag_lowerOF == 2 || flag_lowerOF == 4); 
    bool cutLowerGenSideBin = (flag_lowerOF == 5); 
    
     
     
     
    // GEN LEVEL OF BINS
    int flag_upperOF = SVD_GetDigit(steering, 11, 3); 
    bool cutUpperRecSideBin = (flag_upperOF >= 3);
    bool fixUpperSideBin = (flag_upperOF == 2 || flag_upperOF == 4); 
    bool cutUpperGenSideBin = (flag_upperOF == 5);
     
     
     
     // RESP ORIENTATION
    int flag_respOr = SVD_GetDigit(steering, 12, 2); 
     
     
     
    // NORMALIZATION
    int flag_norm = SVD_GetDigit(steering, 13, 1);  
     
     
     
    // CLOSURE TEST 
    int flag_closure = SVD_GetDigit(steering, 14, 1);  
     
     
     
    // PREWEIGHTING
    int flag_preweighting = SVD_GetDigit(steering, 15, 1);  
     
     
     
    // BACKGROUND HANDLING
    int flag_bgrhandling = SVD_GetDigit(steering, 16, 1);  
    
    
    
    // CONTROL PLOT STYLE
    int flag_plotStyle = SVD_GetDigit(steering, 17, 1);
    if(flag_verbose>=3) std::cout << "flag_plotStyle="<<flag_plotStyle<<std::endl;
    // plot range; per default (-9999) not restricted
    double rangeL = -9999.;
    double rangeR = -9999.;
    if(flag_plotStyle==2) {
      // restrict plot range to relevant bins
      rangeL = thebins[0]+0.001;
      rangeR = thebins[numbins]-0.001;
      if(flag_verbose>=3) std::cout << "Range of axis restricted to: rangeL="<<rangeL<< "; rangeR="<<rangeR<<std::endl;
    }
    
      
 
    // Systematics Flags
    int numberSyst = numSys;
    if ( numberSyst < 0 ) numberSyst = 0;
    bool doSystematics = false;
    if ( numberSyst > 0 ) doSystematics = true;
    
         
    ///////////////////////////////////////////////////////////////////  
    /////////  S E T U P   R O O T //////////////////////////////////// 
    /////////////////////////////////////////////////////////////////// 



    // Automatically create structure for errors
    // in every histogram.
    bool save_sumw2 = TH1::GetDefaultSumw2();
    TH1::SetDefaultSumw2(true);
    
    
    // Turn off graphics mode
    bool save_isbatch = gROOT->IsBatch(); 
    gROOT->SetBatch(true);
    
    
    // Turn of adding histograms to directory
    bool save_addDir = TH1::AddDirectoryStatus();
    TH1::AddDirectory(false);
    
    
    // Adapt gErrorIgnoreLevel to reduce output
    int save_gErrorLv=gErrorIgnoreLevel;
    if (flag_verbose==2) gErrorIgnoreLevel=kWarning;
    if (flag_verbose==1) gErrorIgnoreLevel=kFatal; 

   
    // Output Plot File Name and Path 
    TString outputfilenamePs = psFile;
    TString outputfilenameEps = epsFile;
    TString outputfilenameTxt = txtFile; 
    TString fullRegParFile = regParFile;   
    TString outputfilenameRoot = rootFile;   
      
      
    /////////////////////////////////////////////////////////////////// 
    /////////  R E A D   R E G P A R   //////////////////////////////// 
    /////////////////////////////////////////////////////////////////// 
  
    // Strings to search for
    TString cpqss = SVD_CPQSS(channel, particle, quantity, special, syst);
    TString thekey = SVD_CPQSS(channel, particle, quantity, special, syst);
    TString thekeynom = SVD_CPQSS(channel, particle, quantity, special, "");
     
    
    // Get regularization parameter
    // for minimal regularization
    if ( flag_regmode == 2 ) {
        regPar = (double) numbins;
    }

    // Get regularization paramter from file     
    if ( flag_regmode == 4 ) { 
    
        // Get the Line from File by searching for a key
        TString theLine = SVD_LineFromFile(thekeynom, fullRegParFile);
         
            
        // Read Tau
        TString theTauStr = theLine;
        theTauStr.Remove(113, 100000);
        theTauStr.Remove(0,100);
        double theTau = (double) theTauStr.Atof();
        
        
        // Read K Value
        TString theKStr = theLine;
        theKStr.Remove(125,100000);
        theKStr.Remove(0,121);
        int theK = (int) theKStr.Atoi();
            
            
        // Save Reg Par
        if ( flag_regpar == 2 ) {
            regPar = theTau;
        } else {
            regPar = (double) theK ;
        }
    
        // Output
        if ( flag_verbose > 1 ) { 
            cout << "********************************************************************************************************************" << endl;
            cout << "Reading Regularization Parameter from File: " << endl; 
            cout << "    File:               " << fullRegParFile << endl;
            cout << "    Key:                " << thekey << endl;
            if ( flag_regpar == 2 ) { 
                cout << "    RegPar (tau-Value): " << regPar << endl;
            } else {
                cout << "    RegPar (k-Value):   " << regPar << endl;
            }
            cout << "********************************************************************************************************************" << endl;
        }  
    }
    
     

    // Regularization
    int theKReg = -1;
    double theTau = -1.;
    if ( flag_regpar == 2 ) {
        theKReg = -1;
        theTau = regPar;
    } else {
        theKReg = (int) regPar; 
        theTau = 0.;
    }  
         
    
    /////////////////////////////////////////////////////////////////// 
    /////////  P R E P A R E   R E B I N N I N G  ///////////////////// 
    ///////////////////////////////////////////////////////////////////
 

    // Save old binnings for further reference
    int oldNumBinsData = SVD_NumBins1D(dataInputHist, false);  
    double oldBinBndsData[oldNumBinsData+1];
    SVD_GetBinBoundaries1D(oldBinBndsData, dataInputHist);
    int oldNumBinsBgr = SVD_NumBins1D(bgrInputHist, false);  
    double oldBinBndsBgr[oldNumBinsBgr+1];
    SVD_GetBinBoundaries1D(oldBinBndsBgr, bgrInputHist);
    int oldNumBinsTtbgr = SVD_NumBins1D(ttbgrInputHist, false);  
    double oldBinBndsTtbgr[oldNumBinsTtbgr+1];
    SVD_GetBinBoundaries1D(oldBinBndsTtbgr, ttbgrInputHist);
    int oldNumBinsGen = SVD_NumBins1D(genInputHist, false);  
    double oldBinBndsGen[oldNumBinsGen+1];
    SVD_GetBinBoundaries1D(oldBinBndsGen, genInputHist);
    int oldNumBinsRec = SVD_NumBins1D(recInputHist, false);  
    double oldBinBndsRec[oldNumBinsRec+1];
    SVD_GetBinBoundaries1D(oldBinBndsRec, recInputHist);
    int oldNumBinsRespX = SVD_NumBinsX2D(respInputHist, false);  
    double oldBinBndsRespX[oldNumBinsRespX+1];
    SVD_GetBinBoundariesX2D(oldBinBndsRespX, respInputHist);
    int oldNumBinsRespY = SVD_NumBinsY2D(respInputHist, false);  
    double oldBinBndsRespY[oldNumBinsRespY+1];
    SVD_GetBinBoundariesY2D(oldBinBndsRespY, respInputHist);
    
    
    // How many bins do you need
    //int nbins = SVD_NumNewBins(numbins, !cutLowerGenSideBin, !cutUpperGenSideBin); 
    int nbins = SVD_NumNewBins(numbins, true, true);
     
    // Calculate new bins
    double arrayOfNewBins[nbins+1];
    double* bins = arrayOfNewBins; 
    SVD_NewBins(bins, nbins, thebins, numbins, true, true);   
    
    
    // Create a string of the binnings
    TString oldBinningStr = "";
    for ( int i = 0 ; i <= numbins ; i++ ) { 
        oldBinningStr.Append(TString::Format("%.5f  ", thebins[i]));
    }
    TString newBinningStr = "";
    for ( int i = 0 ; i <= nbins ; i++ ) { 
        newBinningStr.Append(TString::Format("%.5f  ", bins[i]));
    }
        
 
    
    /////////////////////////////////////////////////////////////////// 
    /////////  R E B I N N I N G   //////////////////////////////////// 
    ///////////////////////////////////////////////////////////////////
 
     
    
    // Remove the input histograms from the TDirectory
    SVD_RmDir1D(dataInputHist, 1);
    SVD_RmDir1D(bgrInputHist, 1+numberSyst);
    SVD_RmDir1D(ttbgrInputHist, 1+numberSyst);
    SVD_RmDir1D(genInputHist, 1+numberSyst);
    SVD_RmDir1D(recInputHist, 1+numberSyst);
    SVD_RmDir2D(respInputHist, 1+numberSyst);
    

    // data (including background)
    // ATTENTION: Hier we assume to have only ONE histo in the array1
    TH1D* rawHist = SVD_Rebin1D(dataInputHist, numbins, thebins, !cutLowerGenSideBin, !cutUpperGenSideBin, 1);
    SVD_EmptyRecSideBins1D(rawHist, cutLowerRecSideBin, cutUpperRecSideBin); 
    
    
    // All Background
    // Set errors to zero!
    TH1D* bgrHist = NULL;
    if ( bgrInputHist != NULL ) {
        bgrHist = SVD_Rebin1D((TH1D*) bgrInputHist, numbins, thebins, !cutLowerGenSideBin, !cutUpperGenSideBin, numberSyst+1); 
        SVD_EmptyHistoErrors1D(bgrHist, numberSyst+1);
        SVD_EmptyRecSideBins1D(bgrHist, cutLowerRecSideBin, cutUpperRecSideBin, numberSyst+1);
    } 
  
    // ttbar background only
    TH1D* ttbgrHist = NULL;
    if ( ttbgrInputHist != NULL ) {
        ttbgrHist= SVD_Rebin1D((TH1D*) ttbgrInputHist, numbins, thebins, !cutLowerGenSideBin, !cutUpperGenSideBin, numberSyst+1); 
        SVD_EmptyHistoErrors1D(ttbgrHist, numberSyst+1);
        SVD_EmptyRecSideBins1D(ttbgrHist, cutLowerRecSideBin, cutUpperRecSideBin, numberSyst+1);
    } 
     
    // Response Matrix
    // ... thereby transposing it, if needed
    bool doTranspose = false;
    if ( flag_respOr == 2 ) doTranspose = true; 
    TH2D* mcHist = SVD_Rebin2D((TH2D*) respInputHist, numbins, thebins, numbins, thebins, !cutLowerGenSideBin, !cutUpperGenSideBin, numberSyst+1, doTranspose);  
    SVD_EmptyRecSideBins2D(mcHist, cutLowerRecSideBin, cutUpperRecSideBin, numberSyst+1); 
    SVD_EmptyGenSideBins2D(mcHist, cutLowerGenSideBin, cutUpperGenSideBin, numberSyst+1);  

 
    
    // MC truth signal
    // ATTENTION: Do not empty the side bins here!
    TH1D* xiniHist = SVD_Rebin1D((TH1D*) genInputHist, numbins, thebins, !cutLowerGenSideBin, !cutUpperGenSideBin, numberSyst+1);  
    SVD_EmptyGenSideBins1D(xiniHist, cutLowerGenSideBin, cutUpperGenSideBin, numberSyst+1);  
     
    // Reconstructed MC signal
    TH1D* biniHist = SVD_Rebin1D((TH1D*) recInputHist, numbins, thebins, !cutLowerGenSideBin, !cutUpperGenSideBin, numberSyst+1);
       
 
     
    /////////////////////////////////////////////////////////////////// 
    /////////  B A C K G R O U N D   H A N D L I N G  /////////////////
    ///////////////////////////////////////////////////////////////////
    
    // This is done for all systematics!
     
    // ATTENTION! 
    // (1) Errors from background will be neglected! 
    // (2) If 'ttbgrInputHist' is given, it will be removed by
    //     a multiplicative approach that is xSec independent!   
    // (3) If 'bgrInputHist' is given, it will be removed from the
    //     data by substraction.
    // (4) If both 'bgrInputHist' and 'ttbgrInputHist' are given,
    //     the difference between 'bgrInputHist' and 'ttbgrInputHist'
    //     will be substracted from the data and the background given
    //     int 'ttbgrInputHist' will be treated as described in (2).
    // (5) If neither 'bgrInputHist' nor 'ttbgrInputHist' are given,
    //     no background will be considered.
 
    // ATTENTION!
    // If this is supposed to be a closure test, then the
    // data histogram should just be the biniHist
 
 
    // Placeholder for background free data
    TH1D* dataHist = SVD_CloneHists1D(biniHist,numberSyst+1);
 
    // Background reduction 
    if ( flag_closure != 2 ) {
        SVD_BackgrHandling(dataHist, bgrHist, ttbgrHist, biniHist, rawHist, (flag_bgrhandling == 2), numberSyst+1);
    }
    
    
    // Empty Rec Side Bins 
    SVD_EmptyRecSideBins1D(dataHist, cutLowerRecSideBin, cutUpperRecSideBin, numberSyst+1);
     
     
    ///////////////////////////////////////////////////////////////////
    ////////////   P R E W E I G H T I N G  ///////////////////////////
    /////////////////////////////////////////////////////////////////// 
    
    // Here, the MC model gets reweighted to the data.
    // For this purpose, preliminary unfoldings are performed
    // and then the MC model is reweighted to fit the data.
    // This is an iterative procedure, and the number of iterations
    // is set in the 'steering' parameter.
    
    // Attention:  
    // (1) The unfolding for all iterations is done with the same
    //     regularization parameter.
    // (2) Background Reduction will be done on the basis of the
    //     unreweighted MC!
    // (3) Side bins will never be fixed to MC for preliminary 
    //     unfoldings.
    
      
    // Determine the number of iterations prior to regular unfolding
    int numberPreweights = flag_preweighting - 1; 
    
    // Preweighting factors to be applied to
    // the total number of events
    double* preweightingFactorGen   = new double[1+numberSyst];
    double* preweightingFactorRec   = new double[1+numberSyst];
    for ( int i = 0 ; i < 1+numberSyst ; i++) {
        *(preweightingFactorGen+i) = 1.;
        *(preweightingFactorRec+i) = 1.;
    }
    
    
    // Make Iterations
    TopSVDUnfold*  preliminarySVDUnfold = NULL;
    TH1D* prelUnfHist = NULL;
    TH1D* prelBbbHist = NULL;
    for ( int iter = 0 ; iter < numberPreweights ; iter++) {
         
        // Don't do all this for BBB 
        if ( flag_regmode == 1 ) break;
        
        // Output
        if (flag_verbose>1 && iter ==  0) { 
            cout << endl;
            cout << "*******************************************************************************************************************" << endl;
            cout << "P R E W E I G H T I N G" << endl;
            cout << endl;
        }  
        
        
        // Output
        if (flag_verbose>1 ) {
            cout << "" << endl;
            cout << "    Preweighting, iteration " << iter+1 << " of " << numberPreweights << endl;
            cout << "        Integrals before Reweighting:" << endl;
            cout << "            Background free data (dataHist):                     " <<  SVD_Integral1D(dataHist, 0, true) << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << SVD_Integral1D(dataHist, i, true) << endl;
            }   
            cout << "            Reconstructed MC (biniHist):                         " <<  SVD_Integral1D(biniHist, 0, true) << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << SVD_Integral1D(biniHist, i, true) << endl;
            }  
            cout << "            Generated MC (xiniHist):                             " <<  SVD_Integral1D(xiniHist, 0, true) << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << SVD_Integral1D(xiniHist, i, true) << endl;
            }  
            cout << "            Response Matrix (mcHist):                            " <<  SVD_Integral2D(mcHist, 0, true) << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << SVD_Integral2D(mcHist, i, true) << endl;
            }  
        }
        
        // Save integrals of bini and xini
        double* intBiniBefore   = new double[1+numberSyst];
        double* intXiniBefore   = new double[1+numberSyst];
        for ( int i = 0 ; i < 1+numberSyst ; i++) {
            *(intBiniBefore+i) = SVD_Integral1D(biniHist, i, true);
            *(intXiniBefore+i) = SVD_Integral1D(xiniHist, i, true); 
        }
        
        // BBB Unfolding for comparison
        TH1D* preliminaryBBBBeffHist = SVD_CloneHists1D(xiniHist,numberSyst+1);
        SVD_BBBEff(preliminaryBBBBeffHist, biniHist, xiniHist, numberSyst+1); 
        prelBbbHist = SVD_CloneHists1D(xiniHist,numberSyst+1);
        SVD_BBBUnf(prelBbbHist, dataHist, preliminaryBBBBeffHist, numberSyst+1);  
        
        // Setup tool 
        preliminarySVDUnfold = SVD_SetupTool(dataHist, biniHist, xiniHist, mcHist, theTau, numberSyst+1);  
        
        // Unfolding 
        prelUnfHist = SVD_CallSVDUnfolding(preliminarySVDUnfold, theKReg, numberSyst+1); 
        
        // Calculate Gen Weights 
        TH1D* reweightsGen =  SVD_CloneHists1D(prelUnfHist, numberSyst+1); 
        SVD_MultInvRightToLeft(reweightsGen, xiniHist, 1., 1., 1., numberSyst+1);
        
        // Reweight xiniHist   
        SVD_MultRightToLeft(xiniHist, reweightsGen, 1., 1., numberSyst+1); 
        
        
        // Sum up  events in mcHist before Reweighting
        TH1D* sumOverGenBeforeRW = SVD_CloneHists1D(prelUnfHist, numberSyst+1); 
        SVD_SumOverGenBins(sumOverGenBeforeRW, mcHist, numberSyst+1); 
        
        // Apply weights to mcHist 
        // NOTE: The generator level quantities are
        // ... on the Y Axis
        // ... that is the second index
        // ... which is index j here.
        for ( int h = 0; h < 1+numberSyst ; h++ ) {
            for ( int i = 1; i <= nbins ; i++ ) {
                for ( int j = 1; j <= nbins ; j++ ) {
                    double value = (mcHist+h)->GetBinContent(i,j);
                    double error = (mcHist+h)->GetBinError(i,j);
                    double weight = reweightsGen->GetBinContent(j);
                    value = value * weight;
                    error = error * weight;
                    (mcHist+h)->SetBinContent(i,j, value);
                    (mcHist+h)->SetBinError(i,j, error);
                }
            }
        }
        
        // Sum up  events in mcHist after Reweighting
        TH1D* sumOverGenAfterRW = SVD_CloneHists1D(prelUnfHist, numberSyst+1); 
        SVD_SumOverGenBins(sumOverGenAfterRW, mcHist, numberSyst+1);
        
         
        // Calculate rec level weights
        TH1D* reweightsRec =  SVD_CloneHists1D(sumOverGenAfterRW, numberSyst+1); 
        SVD_MultInvRightToLeft(reweightsRec, sumOverGenBeforeRW, 1., 1., 1., numberSyst+1);
        
        
        // Apply weights to biniHist
        SVD_MultRightToLeft(biniHist, reweightsRec, 1., 1., numberSyst+1); 
        
        // Output
        if (flag_verbose>1 ) {
            cout << "        Integrals after Reweighting:" << endl;
            cout << "            Background free data (dataHist):                     " <<  SVD_Integral1D(dataHist, 0, true) << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << SVD_Integral1D(dataHist, i, true) << endl;
            }   
            cout << "            Reconstructed MC (biniHist):                         " <<  SVD_Integral1D(biniHist, 0, true) << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << SVD_Integral1D(biniHist, i, true) << endl;
            }  
            cout << "            Generated MC (xiniHist):                             " <<  SVD_Integral1D(xiniHist, 0, true) << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << SVD_Integral1D(xiniHist, i, true) << endl;
            }  
            cout << "            Response Matrix (mcHist):                            " <<  SVD_Integral2D(mcHist, 0, true) << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << SVD_Integral2D(mcHist, i, true) << endl;
            }  
            
            cout << "        Helper Distributions:       " << endl;
            cout << "            Unfolding Result (prelUnfHist):                      " <<  SVD_Integral1D(prelUnfHist, 0, true) << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << SVD_Integral1D(prelUnfHist, i, true) << endl;
            }  
            cout << "            Bin by BinResult (prelBbbHist):                      " <<  SVD_Integral1D(prelBbbHist, 0, true) << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << SVD_Integral1D(prelBbbHist, i, true) << endl;
            }  
            cout << "            Generator Level Weights (reweightsGen):              " <<  SVD_Integral1D(reweightsGen, 0, true) << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << SVD_Integral1D(reweightsGen, i, true) << endl;
            }  
            cout << "            Detector Level Weights (reweightsRec):               " <<  SVD_Integral1D(reweightsRec, 0, true) << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << SVD_Integral1D(reweightsRec, i, true) << endl;
            }   
            cout << "            Sum over Gen Bins before RW (sumOverGenBeforeRW):    " <<  SVD_Integral1D(sumOverGenBeforeRW, 0, true) << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << SVD_Integral1D(sumOverGenBeforeRW, i, true) << endl;
            }     
            cout << "            Sum over Gen Bins after RW (sumOverGenAfterRW):      " <<  SVD_Integral1D(sumOverGenAfterRW, 0, true) << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << SVD_Integral1D(sumOverGenAfterRW, i, true) << endl;
            }   
        }
         
        // Save integrals of bini and xini after reweighting
        double* intBiniAfter    = new double[1+numberSyst];
        double* intXiniAfter    = new double[1+numberSyst]; 
        for ( int i = 0 ; i < 1+numberSyst; i++ ) {
            *(intBiniAfter+i) = SVD_Integral1D(biniHist, i, true);
            *(intXiniAfter+i) = SVD_Integral1D(xiniHist, i, true); 
        }
        
        // Global Preweighting Factors for this Iteration          
        double* preweightingFactorGenThisIteration    = new double[1+numberSyst];
        double* preweightingFactorRecThisIteration    = new double[1+numberSyst]; 
        for ( int i = 0 ; i < 1+numberSyst ; i++) {
            *(preweightingFactorGenThisIteration+i) = SVD_Divide(*(intXiniAfter+i),*(intXiniBefore+i));
            *(preweightingFactorRecThisIteration+i) = SVD_Divide(*(intBiniAfter+i),*(intBiniBefore+i));
        }
    
        // Global Preweighting Factors 
        for ( int i = 0 ; i < 1+numberSyst ; i++) { 
            *(preweightingFactorGen+i) = *(preweightingFactorGen+i) * *(preweightingFactorGenThisIteration+i);
            *(preweightingFactorRec+i) = *(preweightingFactorRec+i) * *(preweightingFactorRecThisIteration+i);
        } 
        
        // Output
        if (flag_verbose>1 ) {
            cout << "        Global Factors for Total Event Counts:       " << endl;
            cout << "            Factor for this Iteration (Generator Level):         " <<  *preweightingFactorGenThisIteration << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << *(preweightingFactorGenThisIteration+i) << endl;
            } 
            cout << "            Factor for this Iteration (Reconstruction Level):    " <<  *preweightingFactorRecThisIteration << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << *(preweightingFactorRecThisIteration+i) << endl;
            }
            cout << "            Product of Factors (Generator Level):                " <<  *preweightingFactorGen << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << *(preweightingFactorGen+i) << endl;
            } 
            cout << "            Product of Factors (Reconstruction Level):           " <<  *preweightingFactorRec << endl;
            for ( int i = 1 ; i <= numberSyst ; i++ ) {
                cout << "                Syst. Sample " << i << "                                   " << *(preweightingFactorRec+i) << endl;
            }
            cout << "        Number of Bins:                                          " << nbins << endl; 
            cout << "" << endl;
        }
        
        
        // Delete obsolete Objects
        SVD_DeleteSVD(preliminarySVDUnfold, numberSyst+1);
        SVD_DeleteHists1D(sumOverGenAfterRW, numberSyst+1);
        SVD_DeleteHists1D(sumOverGenBeforeRW, numberSyst+1);
        SVD_DeleteHists1D(reweightsGen, numberSyst+1); 
        SVD_DeleteHists1D(reweightsRec, numberSyst+1);
        SVD_DeleteHists1D(preliminaryBBBBeffHist, numberSyst+1);
        delete[] intXiniBefore;
        delete[] intBiniBefore;
        delete[] intXiniAfter;
        delete[] intBiniAfter;
        delete[] preweightingFactorGenThisIteration;
        delete[] preweightingFactorRecThisIteration;
        
         
        // Output
        if (flag_verbose>1 && iter ==  numberPreweights - 1) { 
            cout << endl;
            cout << "*******************************************************************************************************************" << endl; 
            cout << endl;
        }  
    }
    
    
 
 
 
    ///////////////////////////////////////////////////////////////////
    ////////////   L U M I   S C A L I N G    /////////////////////////
    ///////////////////////////////////////////////////////////////////
 

    // Scaling data and simulation to same statistics to enable comparisons
     
    // Get Lumi Scale Factor
    // This is done only for the nominal sample
    // MUST BE DONE AFTER THE BACKGROUND HANDLING
    double lumiScaleFactor =  SVD_LumiScaleFactor(dataHist, biniHist);  
     
 
    // xini Hist, scaled
    TH1D* xiniScaledHist = SVD_CloneHists1D(xiniHist,numberSyst+1);
    SVD_ArrayScale(xiniScaledHist, lumiScaleFactor, numberSyst+1);
      

    // Scaled Data
    TH1D* dataScaledHist = SVD_CloneHists1D(dataHist,numberSyst+1);
    SVD_ArrayScale(dataScaledHist, 1./lumiScaleFactor, numberSyst+1);
    
    
    
        
    /////////////////////////////////////////////////////////////////// 
    /////////  P R E P A R E   N O R M A L I Z A T I O N  ///////////// 
    ///////////////////////////////////////////////////////////////////
          
          
          
    // Data Normalization 
    if ( totalDataEvents == NULL ) {
          totalDataEvents = SVD_ArrayIntegral1D(dataInputHist, false, 1);
    }  
          
    // Bgr Normalization 
    if ( totalBgrEvents == NULL ) {
          totalBgrEvents = SVD_ArrayIntegral1D(bgrInputHist, false, 1+numberSyst);
    }  
          
    // TtBgr Normalization 
    if ( totalTtBgrEvents == NULL ) {
          totalTtBgrEvents = SVD_ArrayIntegral1D(ttbgrInputHist, false, 1+numberSyst);
    }  
          
    // Rec Normalization 
    if ( totalRecEvents == NULL ) {
          totalRecEvents = SVD_ArrayIntegral1D(recInputHist, false, 1+numberSyst);
    }  
          
    // Gen Normalization 
    if ( totalGenEvents == NULL ) {
          totalGenEvents = SVD_ArrayIntegral1D(genInputHist, false, 1+numberSyst);
    }  
    
    // Apply EXTRA Factors to account for the preweighting 
    for ( int i = 1 ; i <= numberSyst ; i++ ) {
        *(totalGenEvents+i) = *(totalGenEvents+i) * *(preweightingFactorGen+i);
        *(totalRecEvents+i) = *(totalRecEvents+i) * *(preweightingFactorRec+i);
    }  
     
    // Global Event Count
    // (means: Data after background removal)
    double* globalEventYield    = new double[1+numberSyst];
    double* globalEventYieldErr = new double[1+numberSyst];
    SVD_GlobalEventYield(globalEventYield, globalEventYieldErr, totalDataEvents, totalBgrEvents, totalTtBgrEvents, totalRecEvents, 1+numberSyst);
     
     
    // Global Efficiency
    double* globalEfficiency    = new double[1+numberSyst];
    SVD_GlobalEfficiency(globalEfficiency, totalRecEvents, totalGenEvents, 1+numberSyst);
    
    
    // Global Cross Section
    double* globalCrossSection    = new double[1+numberSyst];
    double* globalCrossSectionErr = new double[1+numberSyst];
    SVD_GlobalCrossSection(globalCrossSection, globalCrossSectionErr, globalEventYield, globalEventYieldErr, globalEfficiency, 1+numberSyst);
     
    
      
     
    /////////////////////////////////////////////////////////////////// 
    /////////  C O N T R O L   O U T P U T   ////////////////////////// 
    ///////////////////////////////////////////////////////////////////
    
             
  
    // Screen Output
    if ( flag_verbose > 1 ) {
        cout << endl;
        cout << "********************************************************************************************************************" << endl;
        cout << endl;
        cout << "XX   XX  XX   XX  XXXXXX   XXXXX   XX      XXXXX   XXXX  XX   XX   XXXXX    " << endl;
        cout << "XX   XX  XXX  XX  XX      XX   XX  XX      XX  XX   XX   XXX  XX  XX        " << endl;
        cout << "XX   XX  XX X XX  XXXXX   XX   XX  XX      XX  XX   XX   XX X XX  XX XXXX   " << endl;
        cout << "XX   XX  XX  XXX  XX      XX   XX  XX      XX  XX   XX   XX  XXX  XX   XX   " << endl;
        cout << " XXXXX   XX   XX  XX       XXXXX   XXXXXX  XXXXX   XXXX  XX   XX   XXXXX    " << endl;
        cout << " " << endl; 
        cout << " " << endl;
        cout << "Distribution:" << endl;
        cout << "    Channel                                    " << channel << endl;
        cout << "    Particle                                   " << particle << endl;
        cout << "    Quantity:                                  " << quantity << endl;
        cout << "    Special condition:                         " << special << endl;
        cout << "    Systematic:                                " << syst << endl;
        cout << "        Number Syst Samples:                   " << numberSyst << endl;
        cout << "    " << endl;
        cout << "Files:" << endl;
        cout << "    Output ROOT File:                          " << outputfilenameRoot << endl;
        cout << "        Write to ROOT:                         " << (flag_root == 2 ) << endl;
        cout << "    Output PS File:                            " << outputfilenamePs << endl;
        cout << "    Output EPS File:                           " << outputfilenameEps << endl;
        cout << "        Create PS/EPS Files:                   " << (flag_ps >= 2) << endl; 
        cout << "        Write K Scan Plots to PS/EPS:          " << (flag_ps >= 3 ) << endl;
        cout << "        Write Tau Scan Plots to PS/EPS:        " << (flag_ps >= 4 ) << endl;
        cout << "    Output TXT File (with scan output):        " << outputfilenameTxt << endl;
	cout << "        Write to Text File:                    " << (writeHistToText ) << endl;
        cout << "    File with regularization parameters:       " << fullRegParFile << endl;
        cout << "        Key to search for:                     " << thekey << endl; 
        cout << "    " << endl;
        cout << "Regularization and Scan:" << endl;
        cout << "    Regularization Parameter:                  " << regPar << endl;
        cout << "    Use SVD:                                   " << (flag_regmode >= 2) << endl;
        cout << "    Use K:                                     " << (flag_regpar == 1) << endl;
        cout << "    Use Tau:                                   " << (flag_regpar == 2) << endl;
        cout << "    Minimal Regul.:                            " << (flag_regmode == 2) << endl;
        cout << "    Manual Reg. Par.:                          " << (flag_regmode == 3) << endl;
        cout << "    Reg. Par. from File:                       " << (flag_regmode == 4) << endl;
        cout << "    Scan for best tau:                         " << (flag_scan == 2) << endl;
        cout << "        Number of scan points:                 " << nScanPoints << endl;
        cout << "        Range of scan:                         " << 1/rangefactor << " ... " << rangefactor << endl;
        cout << "    " << endl;
        cout << "Rebinning: " << endl;
        cout << "    Binning without OF:" << endl;
        cout << "        Number of Bins:                        " << numbins << endl;
        cout << "        Binning:                               " << oldBinningStr << endl;
        cout << "    Binning with OF:" << endl;
        cout << "        Number of Bins:                        " << nbins << endl;
        cout << "        Binning:                               " << newBinningStr << endl;
        cout << "    Lower Side Bin:" << endl;
        cout << "        Cut out on rec level:                  " << cutLowerRecSideBin << endl;
        cout << "        Cut out on gen level:                  " << cutLowerGenSideBin << endl;
        cout << "        Fix to MC:                             " << fixLowerSideBin << endl;
        cout << "    Upper Side Bin:" << endl;
        cout << "        Cut out on rec level:                  " << cutUpperRecSideBin << endl;
        cout << "        Cut out on gen level:                  " << cutUpperGenSideBin << endl;
        cout << "        Fix to MC:                             " << fixUpperSideBin << endl; 
        cout << " " << endl;
        cout << "More Steering Options: " << endl;
        cout << "    Verbosity:                                 " << endl;
        cout << "        Standard Output:                       " << (flag_verbose>=2) << endl;
        cout << "        Debugging Output:                      " << (flag_verbose==3) << endl;
        cout << "    Transpose Response Matrix:                 " << (flag_respOr == 2) << endl; 
        cout << "    Plain Steering Parameter:                  " << steering << endl;
        cout << "    Closure Test (Unfold MC):                  " << (flag_closure == 2) << endl;
        cout << "    Number of preweighting iterations:         " << flag_preweighting - 1 << endl;
        cout << " " << endl; 
        cout << "Integrals of Input Distributions:" << endl;
        cout << "    Data:                                      " << SVD_Integral1D(dataInputHist, 0, true) << endl;
        cout << "    Background:                                " << SVD_Integral1D(bgrInputHist, 0, true) << endl;
        for ( int i = 1 ; i <= numberSyst ; i++ ) {
            cout << "        Syst. Sample " << i << "                         " << SVD_Integral1D(bgrInputHist, i, true) << endl;
        }           
        cout << "    Tt-Background:                             " << SVD_Integral1D(ttbgrInputHist, 0, true)  << endl;
        for ( int i = 1 ; i <= numberSyst ; i++ ) {
            cout << "        Syst. Sample " << i << "                         " << SVD_Integral1D(ttbgrInputHist, i, true) << endl;
        }          
        cout << "    Generated MC:                              " << SVD_Integral1D(genInputHist, 0, true) << endl;
        for ( int i = 1 ; i <= numberSyst ; i++ ) {
            cout << "        Syst. Sample " << i << "                         " << SVD_Integral1D(genInputHist, i, true) << endl;
        }   
        cout << "    Reconstructed MC:                          " << SVD_Integral1D(recInputHist, 0, true) << endl;
        for ( int i = 1 ; i <= numberSyst ; i++ ) {
            cout << "        Syst. Sample " << i << "                         " << SVD_Integral1D(recInputHist, i, true) << endl;
        }   
        cout << "    Response Matrix:                           " << SVD_Integral2D(respInputHist, 0, true) << endl;
        for ( int i = 1 ; i <= numberSyst ; i++ ) {
            cout << "        Syst. Sample " << i << "                         " << SVD_Integral2D(respInputHist, i, true) << endl;
        }  
        cout << " " << endl; 
        cout << "Normalization:                                " << endl; 
        cout << "    Extrinsic Normalization:                  " << (flag_norm == 1 ) << endl;
        cout << "    Intrinsic Normalization:                  " << (flag_norm == 2 ) << endl;
        cout << "    Intrinsic Normalization (Pseudo-Exp.):    " << (flag_norm == 3 ) << endl;
        cout << "    Total number of data events:              " << SVD_DoubleFromArray(totalDataEvents, 0)  << endl; 
        cout << "    Total number of background events:        " << SVD_DoubleFromArray(totalBgrEvents, 0) << endl;
        for ( int i = 1 ; i <= numberSyst ; i++ ) {
            cout << "        Syst. Sample " << i << "                        " << SVD_DoubleFromArray(totalBgrEvents, i) << endl;
        } 
        cout << "    Total number of tt background events:     " << SVD_DoubleFromArray(totalTtBgrEvents, 0) << endl;
        for ( int i = 1 ; i <= numberSyst ; i++ ) {
            cout << "        Syst. Sample " << i << "                        " << SVD_DoubleFromArray(totalTtBgrEvents, i) << endl;
        } 
        cout << "    Total number of reconstructed events:     " << SVD_DoubleFromArray(totalRecEvents, 0) << endl;
        for ( int i = 1 ; i <= numberSyst ; i++ ) {
            cout << "        Syst. Sample " << i << "                        " << SVD_DoubleFromArray(totalRecEvents, i) << endl;
        } 
        cout << "    Total number of generated events:         " << SVD_DoubleFromArray(totalGenEvents, 0) << endl;
        for ( int i = 1 ; i <= numberSyst ; i++ ) {
            cout << "        Syst. Sample " << i << "                        " << SVD_DoubleFromArray(totalGenEvents, i) << endl;
        }  
        cout << "    Global Event Yield:                       " << SVD_DoubleFromArray(globalEventYield, 0) << " +/- " << SVD_DoubleFromArray(globalEventYieldErr, 0) << endl;
        for ( int i = 1 ; i <= numberSyst ; i++ ) {
            cout << "        Syst. Sample " << i << "                        " << SVD_DoubleFromArray(globalEventYield, i) << " +/- " << SVD_DoubleFromArray(globalEventYieldErr, i) << endl;
        } 
        cout << "    Global Unfolded Event Yield:              " << SVD_DoubleFromArray(globalCrossSection, 0) << " +/- " << SVD_DoubleFromArray(globalCrossSectionErr, 0) << endl;
        for ( int i = 1 ; i <= numberSyst ; i++ ) {
            cout << "        Syst. Sample " << i << "                        " << SVD_DoubleFromArray(globalCrossSection, i) << " +/- " << SVD_DoubleFromArray(globalCrossSectionErr, i) << endl;
        }   
        cout << "    Global Efficiency:                        " << SVD_DoubleFromArray(globalEfficiency, 0) << endl;
        for ( int i = 1 ; i <= numberSyst ; i++ ) {
            cout << "        Syst. Sample " << i << "                        " << SVD_DoubleFromArray(globalEfficiency, i) << endl;
        }  
        cout << endl;
        if ( flag_verbose > 2 ) {
            cout <<"Input Binnings:" << endl;
            cout <<"     Data Histogram (" << oldNumBinsData << "+2 Bins):"; 
            for ( int i = 0; i < oldNumBinsData+1 ; i++ ) {
                if ( i % 20 == 0 ) { cout << endl; cout << "        ";}
                cout << oldBinBndsData[i] << " " ;
            }
            cout << endl; 
            cout <<"     Bgr Histogram (" << oldNumBinsBgr << "+2 Bins):"; 
            for ( int i = 0; i < oldNumBinsBgr+1 ; i++ ) {
                if ( i % 20 == 0 ) { cout << endl; cout << "        ";}
                cout << oldBinBndsBgr[i] << " " ;
            }
            cout << endl; 
            cout <<"     Ttbgr Histogram (" << oldNumBinsTtbgr << "+2 Bins):"; 
            for ( int i = 0; i < oldNumBinsTtbgr+1 ; i++ ) {
                if ( i % 20 == 0 ) { cout << endl; cout << "        ";}
                cout << oldBinBndsTtbgr[i] << " " ;
            }
            cout << endl; 
            cout <<"     Gen Histogram (" << oldNumBinsGen << "+2 Bins):"; 
            for ( int i = 0; i < oldNumBinsGen+1 ; i++ ) {
                if ( i % 20 == 0 ) { cout << endl; cout << "        ";}
                cout << oldBinBndsGen[i] << " " ;
            }
            cout << endl; 
            cout <<"     Rec Histogram (" << oldNumBinsRec << "+2 Bins):"; 
            for ( int i = 0; i < oldNumBinsRec+1 ; i++ ) {
                if ( i % 20 == 0 ) { cout << endl; cout << "        ";}
                cout << oldBinBndsRec[i] << " " ;
            }
            cout << endl; 
            cout <<"     Resp Histogram in X (" << oldNumBinsRespX << "+2 Bins):"; 
            for ( int i = 0; i < oldNumBinsRespX+1 ; i++ ) {
                if ( i % 20 == 0 ) { cout << endl; cout << "        ";}
                cout << oldBinBndsRespX[i] << " " ;
            }
            cout << endl; 
            cout <<"     Resp Histogram in Y (" << oldNumBinsRespY << "+2 Bins):"; 
            for ( int i = 0; i < oldNumBinsRespY+1 ; i++ ) {
                if ( i % 20 == 0 ) { cout << endl; cout << "        ";}
                cout << oldBinBndsRespY[i] << " " ;
            } 
        }
        cout << endl; 
        cout << "********************************************************************************************************************" << endl; 
        
    }        
    
 
    ///////////////////////////////////////////////////////////////////
    ////////////   S M E A R - I N - O U T  ///////////////////////////
    ///////////////////////////////////////////////////////////////////
     
    // Smaer out: Fraction of events generated inside but 
    // reconstructed outside the chosen binning.

    // Smear in: Fraction of generated events outside but
    // reconstructed inside the chosen binning.

    
     
    // Make Copies of input histograms 
    TH1D* dataHist_forTSVD = SVD_CloneHists1D(dataHist,numberSyst+1); 
    TH1D* biniHist_forTSVD = SVD_CloneHists1D(biniHist,numberSyst+1);
    TH1D* xiniHist_forTSVD = SVD_CloneHists1D(xiniHist,numberSyst+1);
    TH2D* mcHist_forTSVD = SVD_CloneHists2D(mcHist,numberSyst+1); 
    
    
    // LOWER SMEAR IN
    TH1D* lowerSmearinHist = SVD_CloneHists1D(biniHist, numberSyst+1);
    SVD_LowerSmearin(lowerSmearinHist, mcHist, numberSyst+1); 
    TString lowerSmearinStr = SVD_PlotName(channel, particle, quantity, special, syst, "LSI");
    SVD_SetTitles1D(lowerSmearinHist, lowerSmearinStr, quantityTex, "Lower Smear-In", numberSyst+1);
    
    
    // UPPER SMEAR IN
    TH1D* upperSmearinHist = SVD_CloneHists1D(biniHist, numberSyst+1);
    SVD_UpperSmearin(upperSmearinHist, mcHist, numberSyst+1); 
    TString upperSmearinStr = SVD_PlotName(channel, particle, quantity, special, syst, "USI");
    SVD_SetTitles1D(upperSmearinHist, lowerSmearinStr, quantityTex, "Upper Smear-In", numberSyst+1);
    
    
    // LOWER SMEAR OUT
    TH1D* lowerSmearoutHist = SVD_CloneHists1D(biniHist, numberSyst+1);
    SVD_LowerSmearout(lowerSmearoutHist, mcHist, numberSyst+1); 
    TString lowerSmearoutStr = SVD_PlotName(channel, particle, quantity, special, syst, "LSO");
    SVD_SetTitles1D(lowerSmearoutHist, lowerSmearoutStr, quantityTex, "Lower Smear-In", numberSyst+1);
    
    
    // UPPER SMEAR OUT
    TH1D* upperSmearoutHist = SVD_CloneHists1D(biniHist, numberSyst+1);
    SVD_UpperSmearout(upperSmearoutHist, mcHist, numberSyst+1); 
    TString upperSmearoutStr = SVD_PlotName(channel, particle, quantity, special, syst, "USO");
    SVD_SetTitles1D(upperSmearoutHist, lowerSmearoutStr, quantityTex, "Upper Smear-In", numberSyst+1);
    
    
    // Remove Smearin from lower side bin 
    if ( fixLowerSideBin == true ) {
        SVD_AddRightToLeft(dataHist_forTSVD, lowerSmearinHist, 1., -1.*lumiScaleFactor, numberSyst+1); 
        SVD_AddRightToLeft(biniHist_forTSVD, lowerSmearinHist, 1., -1., numberSyst+1);
        SVD_EmptyGenSideBins1D(xiniHist_forTSVD, true, false, numberSyst+1);
        SVD_EmptyGenSideBins2D(mcHist_forTSVD, true, false, numberSyst+1);  
    }
    
    
    // Remove Smearin from upper side bin 
    if ( fixUpperSideBin == true ) {
        SVD_AddRightToLeft(dataHist_forTSVD, upperSmearinHist, 1., -1.*lumiScaleFactor, numberSyst+1); 
        SVD_AddRightToLeft(biniHist_forTSVD, upperSmearinHist, 1., -1., numberSyst+1);
        SVD_EmptyGenSideBins1D(xiniHist_forTSVD, false, true, numberSyst+1);
        SVD_EmptyGenSideBins2D(mcHist_forTSVD, false, true, numberSyst+1);  
    }
        
        
    
       
 
    ///////////////////////////////////////////////////////////////////
    ////////////   B I N   B Y   B I N   U N F O L D I N G   //////////
    ///////////////////////////////////////////////////////////////////
 

    // This is done for all systematics!
 
    // BBB-Efficiency
    TH1D* beffHist = SVD_CloneHists1D(xiniHist,numberSyst+1);
    SVD_BBBEff(beffHist, biniHist, xiniHist, numberSyst+1); 
  
  
    // Do BBB Unfolding 
    // Output stored in 'bbbHist'
    // ATTENTION:
    // (1) Use BBB Efficiency.
    // (2) Although this might be tempting,
    //     do not include the error from the efficiency.
    //     This would bust your BBB / SVD comparison.
    //     Instead, add the error on the efficiency
    //     later on as a systematic MC uncertainty! 
    TH1D* bbbHist = SVD_CloneHists1D(xiniHist,numberSyst+1);
    SVD_BBBUnf(bbbHist, dataHist, beffHist, numberSyst+1); 

     
    ///////////////////////////////////////////////////////////////////
    ////////////   U N F O L D I N G //////////////////////////////////
    ///////////////////////////////////////////////////////////////////
    
    
    // This is done for all systematics! 
    
    

    // DATA COVARIANCE
    // ATTENTION! This establishes the statistical uncertainty!
    // I think ...
    // ... it should not be shifted for systematics!
    TH2D* dataCovHist = SVD_CloneHists2D(mcHist); 
    for ( int i=1; i<=nbins; i++ ) {
        for ( int j = 1; j <= nbins ; j++ ) {
            dataCovHist->SetBinContent(i,j,0.);
            dataCovHist->SetBinError(i,j,0.);
        } 
    } 
    for (int i=1; i<=nbins; i++ ) {
       double variance = TMath::Power(dataHist->GetBinError(i), 2.);
       dataCovHist->SetBinContent(i,i,variance);
    } 
    

    // Setup Unfolding Tools
    // One for each Systematic Sample
    // This way, they are all independent from each other 
    // UNFOLDINGPROGRAM
    TopSVDUnfold*  mySVDUnfold = SVD_SetupTool(dataHist_forTSVD, biniHist_forTSVD, xiniHist_forTSVD, mcHist_forTSVD, theTau, numberSyst+1); 
   
    // Unfolding 
    TH1D* unfHist = NULL;
    if ( flag_regmode > 1 ) {
        // UNFOLDINGPROGRAM
        unfHist = SVD_CallSVDUnfolding(mySVDUnfold, theKReg, numberSyst+1);
    } else {
        // The BBB Only Case
        unfHist = SVD_CloneHists1D(bbbHist,numberSyst+1); 
    } 
    
    // Get more output
    TH1D* weightHist = NULL;
    TH1D* ddHist = NULL;
    TH1D* svHist = NULL;
    if ( flag_regmode > 1 ) {
        // UNFOLDINGPROGRAM
        weightHist = SVD_RetrieveWeights(mySVDUnfold, numberSyst+1); // Weights = x^unf / x^gen, TSVDUnfold-specific
        ddHist = SVD_RetrieveDD(mySVDUnfold, numberSyst+1); // D-Values
        svHist = SVD_RetrieveSV(mySVDUnfold, numberSyst+1); // SV-Values
    } else { 
        // The BBB Only Case 
        weightHist = SVD_CloneHists1D(unfHist,numberSyst+1);
        SVD_MultInvRightToLeft(weightHist, xiniHist_forTSVD, 1., 1., 1., numberSyst+1); 
        ddHist = SVD_CloneHists1D(unfHist,numberSyst+1); 
        SVD_EmptyHistogram1D(ddHist, numberSyst+1);
        svHist = SVD_CloneHists1D(unfHist,numberSyst+1);
        SVD_EmptyHistogram1D(svHist, numberSyst+1);  
    }
     

    // Add Smearin from lower side bin
    if ( fixLowerSideBin == true ) { 
        double theWeightSideBin = lumiScaleFactor;  
        SVD_AddBinRightToLeft(unfHist, xiniHist, 1, 0., lumiScaleFactor, numberSyst+1);  
        SVD_SetBinValErr1D(weightHist, 1, theWeightSideBin, 0.,  numberSyst+1);
    }
         
    // Add Smearin from upper side bin
    if ( fixUpperSideBin == true ) { 
        double theWeightSideBin = lumiScaleFactor;  
        SVD_AddBinRightToLeft(unfHist, xiniHist, nbins, 0., lumiScaleFactor, numberSyst+1);  
        SVD_SetBinValErr1D(weightHist, nbins, theWeightSideBin, 0.,  numberSyst+1);
    } 
    
     
    ///////////////////////////////////////////////////////////////////
    ////////////   C O V A R I A N C E    /////////////////////////////
    ///////////////////////////////////////////////////////////////////
    
    // This is done for the nominal sample only!
    
     
    // Calculate TOTAL COVARIANCE
    TH2D* statCovHist = NULL;
    TH2D* statCovHistNorm = NULL;
    TH2D* mcCovHist = NULL;
    if ( flag_regmode > 1 ) {
        // UNFOLDINGPROGRAM
        statCovHist = SVD_CloneHists2D(mySVDUnfold->GetUnfoldCovMatrix(dataCovHist, nExperiments), 1);
        statCovHistNorm = SVD_CloneHists2D(mySVDUnfold->GetUnfoldCovMatrixNorm(dataCovHist, nExperiments, 1, flag_norm, flag_verbose), 1);
        mcCovHist = SVD_CloneHists2D(mySVDUnfold->GetAdetCovMatrix(nExperiments), 1);
    } else {
        // The BBB only case 
        // Note: MC Cov hist is empty in this case.
        statCovHist = SVD_CloneHists2D(mcHist,1);   
        SVD_EmptyHistogram2D(statCovHist, 1); 
        for (int i=1; i<=nbins; i++ ) {
           double variance = TMath::Power(bbbHist->GetBinError(i), 2.);
           statCovHist->SetBinContent(i,i,variance);
        } 
        mcCovHist = SVD_CloneHists2D(mcHist,1);
        SVD_EmptyHistogram2D(mcCovHist, 1); 
    }    
     
    TH2D* totCovHist = SVD_CloneHists2D(statCovHist); 
    totCovHist->Add(mcCovHist);
 
    // STAT ERROR ON RESULT
    // ATTENTION:
    // Use only the statistical uncertainty!
    // Not the MC uncertainty or anything else!
    // Otherwise, comparisons to standard BBB unfolding
    // will be very difficut. 
    for ( int i = 1 ; i <= nbins ; i++ ) {
        double errorsq = statCovHist->GetBinContent(i,i);
        double error = 0.;
        if ( errorsq > 0. ) {
            error = TMath::Sqrt(errorsq);
        } 
        unfHist->SetBinError(i, error); 
    }
    
 
    
    ///////////////////////////////////////////////////////////////////
    ////////////   K   S C A N     ////////////////////////////////////
    /////////////////////////////////////////////////////////////////// 
    
    // This is done for the nominal sample only!
    
    
    // Count the Singular Values
    int nSingularValues = 0;
    int nSignSingularValues = 0;
    int nScanSingularValues = 0;  
    for ( int i = 0 ; i < nbins ; i++ ) { 
        double sv = svHist->GetBinContent(i+1);
        double dd = ddHist->GetBinContent(i+1);
        if ( sv > 0. ) {
            nSingularValues++ ;
            if ( dd > 1. ) nSignSingularValues++ ;
            if ( dd > 0.01 ) nScanSingularValues++ ; 
        }
    }  
      
     
    // Save Singular Values in TVector  
    TVectorD vSingularValues(nSingularValues);
    for ( int i = 0 ; i < svHist->GetNbinsX() ; i++ ) {
        double sv = svHist->GetBinContent(i+1);  
        if ( sv > 0. ) vSingularValues[i] = sv; 
    }  
    
    
    // Setup Arrays to hold the histograms
    TH1D* arrK_StatError = NULL;
    TH1D* arrK_Weight = NULL;
    TH1D* arrK_Mean = NULL;
    TH1D* arrK_GlC = NULL;
      
    
    if ( flag_ps >= 3 && nScanSingularValues > 1) { 
        
        // Create new Objects
        arrK_StatError = new TH1D[nScanSingularValues-1];
        arrK_Weight = new TH1D[nScanSingularValues-1];
        arrK_GlC = new TH1D[nScanSingularValues-1];
        arrK_Mean = new TH1D[nScanSingularValues-1];


        // Do Scan
        int kscancounter = 0;
        for ( int tmpK = 2 ; tmpK <= nScanSingularValues ; tmpK++ ) {  
            
            // Don't do this for the BBB only case
            if ( flag_regmode == 1 ) break;
            
            // Unfold with temporary k
            TH1D* tmpUnfResult = SVD_CloneHists1D(mySVDUnfold->Unfold(tmpK));  
            SVD_AddBinRightToLeft(tmpUnfResult, xiniHist, 1, 0., lumiScaleFactor, 1);  
            TString unfStr = SVD_PlotName(channel, particle, quantity, special, syst, "UNF");
            SVD_SetTitles1D(tmpUnfResult, unfStr, quantityTex, "Events");
                
            // Get Weights
            TH1D* tmpWeights = SVD_CloneHists1D(mySVDUnfold->GetWeights()); 
            TString weightStr = SVD_PlotName(channel, particle, quantity, special, syst, "WGT");
            SVD_SetTitles1D(tmpWeights, weightStr, quantityTex, "Weights");  
              
            // Add Smearin from lower side bin
            if ( fixLowerSideBin == true ) { 
                double theWeightSideBin = lumiScaleFactor;   
                SVD_AddBinRightToLeft(tmpUnfResult, xiniHist, 1, 0., lumiScaleFactor, 1);  
                SVD_SetBinValErr1D(tmpWeights, 1, theWeightSideBin, 0., 1);
            }
                 
            // Add Smearin from upper side bin
            if ( fixUpperSideBin == true ) { 
                double theWeightSideBin = lumiScaleFactor;  
                SVD_AddBinRightToLeft(tmpUnfResult, xiniHist, nbins, 0., lumiScaleFactor, 1);  
                SVD_SetBinValErr1D(tmpWeights, nbins, theWeightSideBin, 0., 1);
            } 
            
            // Scale Weights 
            tmpWeights->Scale(1./lumiScaleFactor);
    
            // Calculate the Error matrix (internally, more unfoldings ... )
            TH2D* tmpCovHist = SVD_CloneHists2D(mySVDUnfold->GetUnfoldCovMatrix(dataCovHist, nExperiments));  
            TH1D* tmpErr = SVD_Cov2Err(tmpCovHist, tmpUnfResult, TString::Format("arrK_Mean_%i", tmpK), quantityTex, "Stat.");
            TString statErrStr = SVD_PlotName(channel, particle, quantity, special, syst, "STATERR"); 
            SVD_SetTitles1D(tmpErr, statErrStr, quantityTex, "Stat. Err. in %");  
              
            
            // Global Correlation
            TH1D* tmpGlC = SVD_CalcGlobCorr(tmpCovHist, dataHist); 
            TString glcStr = SVD_PlotName(channel, particle, quantity, special, syst, "GLOBC");
            SVD_SetTitles1D(tmpGlC, glcStr, quantityTex, "Glob. Corr. in %");
                   
 
            // Save Objects in Array
            arrK_StatError[kscancounter] = *tmpErr;
            arrK_Weight[kscancounter] = *tmpWeights;
            arrK_GlC[kscancounter] = *tmpGlC;
            arrK_Mean[kscancounter] = *tmpUnfResult;
 
 
            // Remove them from a potential TDirectory
            (arrK_StatError+kscancounter)->SetDirectory(NULL);
            (arrK_Weight+kscancounter)->SetDirectory(NULL);
            (arrK_GlC+kscancounter)->SetDirectory(NULL); 
            (arrK_Mean+kscancounter)->SetDirectory(NULL);
             
                  
            // Delete obsolete objects
            SVD_DeleteHists1D(tmpUnfResult);
            SVD_DeleteHists1D(tmpWeights);
            SVD_DeleteHists2D(tmpCovHist);
            SVD_DeleteHists1D(tmpErr); 
            SVD_DeleteHists1D(tmpGlC); 
     
            // Increment Counter
            kscancounter++;
        } 
         
    }
   
    ///////////////////////////////////////////////////////////////////
    ////////////   T A U    S C A N     ///////////////////////////////
    ///////////////////////////////////////////////////////////////////  

 
    // This is done for the nominal sample only!
 
    // Create Scan Vectors
    TVectorD* vScanPoints = NULL; 
    TVectorD* vGlobCorr = NULL; 
    TVectorD* vChiSq = NULL; 
    TVectorD* vCurv = NULL;    
    TVectorD* vAvgSqErr = NULL;
    TVectorD* vAvgMean = NULL;    
    double optimalTauX     = -1.;
    double optimalTauY     = -1.;
    double optimalChiSq    = -1.;
    double optimalCurv     = -1.;
    double optimalAvgSqErr = -1.;
    double optimalAvgMean  = -1.;
    int optimalTauLowK = -1;
    int optimalTauHighK = -1;

    // Doubles for BBB values
    double bbbAvgSqErr = 0.;
    double bbbAvgMean = 0.;  
 
  
    if ( flag_scan == 2 ) {
            
             
        // Range for Scan
        double lowTau = 0.;
        double highTau = 0.;
        if ( flag_regpar == 2 ) {
            lowTau = (1./rangefactor) * theTau;
            highTau = rangefactor * theTau;
        } else{
            // Find interesting singular value 
            // Limit the kReg to the Modes which are actually there! 
            int effectiveKReg = theKReg;
            if ( effectiveKReg >= nSingularValues ) effectiveKReg = nSingularValues-1;
            double singularValueK =  vSingularValues[effectiveKReg];
    
            lowTau = (1./rangefactor) * singularValueK;
            highTau = rangefactor * singularValueK;
        }  

        // Slow scan for plots
        if ( flag_ps >= 4 ) {       
              
            // Get the Scan Points
            vScanPoints = SVD_CalcScanPoints(lowTau, highTau, nScanPoints);
    
            // Create the vectors
            vGlobCorr = new TVectorD(nScanPoints);
            vChiSq = new TVectorD(nScanPoints);
            vCurv = new TVectorD(nScanPoints);
            vAvgSqErr = new TVectorD(nScanPoints);
            vAvgMean = new TVectorD(nScanPoints);
    
            // Output
            if (flag_verbose>1) { 
                cout << "*******************************************************************************************************************" << endl;
                cout << "Perform Tau Scan from " << lowTau << " to " << highTau << " for plotting purpose" << endl;
            }  
    
            // Perform the Scan
            for ( int i = 0 ; i < nScanPoints ; i++ ) {
       
                // Do the unfolding with k = -1
                double tau = (*vScanPoints)[i];
                mySVDUnfold->SetTau(tau);  // UNFOLDINGPROGRAM
                TH1D* tmpUnfResult = SVD_CloneHists1D(mySVDUnfold->Unfold(-1)); 
                TH1D* tmpWeights = SVD_CloneHists1D(mySVDUnfold->GetWeights()); 
                double tmpCurv = mySVDUnfold->GetCurv();
                 
                // Add Smearin from lower side bin
                if ( fixLowerSideBin == true ) { 
                    double theWeightSideBin = lumiScaleFactor;  
                    SVD_AddBinRightToLeft(tmpUnfResult, xiniHist, 1, 0., lumiScaleFactor, 1);  
                    SVD_SetBinValErr1D(tmpWeights, 1, theWeightSideBin, 0., 1);
                }
                     
                // Add Smearin from upper side bin
                if ( fixUpperSideBin == true ) { 
                    double theWeightSideBin = lumiScaleFactor;  
                    SVD_AddBinRightToLeft(tmpUnfResult, xiniHist, nbins, 0., lumiScaleFactor, 1);  
                    SVD_SetBinValErr1D(tmpWeights, nbins, theWeightSideBin, 0., 1);
                } 
                   
                // Calculate the Error matrix (internally, more unfoldings ... )
                TH2D* tmpCovHist = SVD_CloneHists2D(mySVDUnfold->GetUnfoldCovMatrix(dataCovHist, nExperiments)); 
                   
                // Global Correlation 
                (*vGlobCorr)[i] = SVD_ScanGlobalCorrelation(tmpCovHist, dataHist);
                   
                // Chi Squared
                (*vChiSq)[i] = SVD_ScanChiSquared(tmpWeights, dataHist, mcHist);
                   
                // Curvature
                (*vCurv)[i] = SVD_ScanCurvature(tmpCurv);
                   
                // Avg Err
                (*vAvgSqErr)[i] = SVD_ScanAvgSqErr(tmpCovHist);
                   
                // Avg Mean
                (*vAvgMean)[i] = SVD_ScanAvgMean(tmpUnfResult);
                   
                // Output
                if(flag_verbose>1){
                    cout << "    Scanpoint " << i << " of " << nScanPoints << " ->  Tau = " << tau << " -> GlobCorr: " << (*vGlobCorr)[i] << endl;
                }
                  
                // delete things
                SVD_DeleteHists1D(tmpUnfResult);
                SVD_DeleteHists1D(tmpWeights);
                SVD_DeleteHists2D(tmpCovHist); 
            }
    
            // Output
            if (flag_verbose>1) { 
                cout << "*******************************************************************************************************************" << endl; 
            }  
             
        }

        // Before starting the golden section search, we have to reset the unfolding class
        mySVDUnfold->Unfold(-1);

        // Determine optimal tau value with golden section     
        double goldSec = (3.-sqrt(5.))/2.;
        double globCorrLowTau = 1000.;
        double globCorrHighTau = 1000.;
        double newLowTau  = TMath::Power(10.,TMath::Log10(lowTau) + goldSec*(TMath::Log10(highTau)-TMath::Log10(lowTau)));
        double newHighTau = TMath::Power(10.,TMath::Log10(highTau) - goldSec*(TMath::Log10(highTau)-TMath::Log10(lowTau)));
    
    
        // Output
        if (flag_verbose>1) { 
            cout << "*******************************************************************************************************************" << endl;
            cout << "Perform Tau Scan from " << lowTau << " to " << highTau << " with Golden Section Search" << endl;
        }  
            
        for ( int i = 0 ; i <= 100 ; i++ ) {
            if( optimalTauX!=newLowTau ){
                // Do the unfolding for lowTau
                mySVDUnfold->SetTau(newLowTau);
                // Calculate the Error matrix (internally, more unfoldings ... )
                TH2D* tmpCovHist = SVD_CloneHists2D(mySVDUnfold->GetUnfoldCovMatrix(dataCovHist, nExperiments)); 
                // Global Correlation 
                globCorrLowTau = SVD_ScanGlobalCorrelation(tmpCovHist, dataHist);
                // delete Error matrix
                SVD_DeleteHists2D(tmpCovHist);
                // print global correlation for new tau
                if(flag_verbose>1)cout << "    newLowTau: " << newLowTau << " -> globCorr: " << globCorrLowTau << endl;
            }
            if(optimalTauX!=newHighTau){
                // Do the unfolding for highTau
                mySVDUnfold->SetTau(newHighTau);
                // Calculate the Error matrix (internally, more unfoldings ... ) 
                TH2D* tmpCovHist = SVD_CloneHists2D(mySVDUnfold->GetUnfoldCovMatrix(dataCovHist, nExperiments));  
                // Global Correlation 
                globCorrHighTau = SVD_ScanGlobalCorrelation(tmpCovHist, dataHist);
                // delete Error matrix
                SVD_DeleteHists2D(tmpCovHist); 
                // print global correlation for new tau
                if(flag_verbose>1)cout << "    newHighTau: " << newHighTau << " -> globCorr: " << globCorrHighTau << endl;
            }
            if(globCorrLowTau<globCorrHighTau){
                highTau=newHighTau;
                optimalTauX=newLowTau;
                optimalTauY=globCorrLowTau;
                newHighTau=newLowTau;
                if(TMath::Abs(globCorrLowTau-globCorrHighTau)/(globCorrLowTau+globCorrHighTau)<0.000001)break;
                globCorrHighTau=globCorrLowTau;
                newLowTau=TMath::Power(10.,TMath::Log10(lowTau) + goldSec*(TMath::Log10(highTau)-TMath::Log10(lowTau)));
            } else {
                lowTau=newLowTau;
                optimalTauX=newHighTau;
                optimalTauY=globCorrHighTau;
                newLowTau=newHighTau;
                if(TMath::Abs(globCorrLowTau-globCorrHighTau)/(globCorrLowTau+globCorrHighTau)<0.000001)break;
                globCorrLowTau=globCorrHighTau;
                newHighTau=TMath::Power(10.,TMath::Log10(highTau) - goldSec*(TMath::Log10(highTau)-TMath::Log10(lowTau)));
            }
            if(i==100)cout << "Warning: Tau scan ended before reaching break condition" << endl;
        }
        if(flag_verbose>1)cout << "Optimal Tau = " << optimalTauX << endl;
     
    
        // Output
        if (flag_verbose>1) { 
            cout << "*******************************************************************************************************************" << endl; 
        }  
     
        // Find closest K values
        optimalTauLowK = 0;
        for ( int i = 0 ; i < nSingularValues ; i++ ) {
            double sv = vSingularValues[i];
            if ( sv > optimalTauX ) optimalTauLowK++;
        }
        optimalTauHighK=optimalTauLowK+1;
     
        // Evaluation at optimal tau
        mySVDUnfold->SetTau(optimalTauX); 
        TH1D* tmpUnfResult = SVD_CloneHists1D(mySVDUnfold->Unfold(-1));  
        TH1D* tmpWeights = SVD_CloneHists1D(mySVDUnfold->GetWeights());  
        double tmpCurv = mySVDUnfold->GetCurv();
        
        // Calculate the Error matrix (internally, more unfoldings ... )
        TH2D* tmpCovHist = SVD_CloneHists2D(mySVDUnfold->GetUnfoldCovMatrix(dataCovHist, nExperiments)); 
        
        // Chi Squared
        optimalChiSq = SVD_ScanChiSquared(tmpWeights, dataHist, mcHist);
        
        // Curvature
        optimalCurv = SVD_ScanCurvature(tmpCurv);
        
        // Avg Err
        optimalAvgSqErr = SVD_ScanAvgSqErr(tmpCovHist);
        
        // Avg Mean
        optimalAvgMean = SVD_ScanAvgMean(tmpUnfResult);
    
        // delete things
        SVD_DeleteHists1D(tmpUnfResult);
        SVD_DeleteHists1D(tmpWeights);
        SVD_DeleteHists2D(tmpCovHist);  
    
        // Calculations for BBB
        // Do not count the OF bins!
        // Recall: You have 2 sets of OF bins
        bbbAvgSqErr = 0.;
        bbbAvgMean = 0.;
        int bincounter = 0;
        for ( int i = 2 ; i < nbins ; i++ ) {
            double bbb_mean = bbbHist->GetBinContent(i);
            double bbb_err = bbbHist->GetBinError(i);
            if ( bbb_mean <= 0.) continue;
                  
            bincounter++;
            bbbAvgSqErr += bbb_err*bbb_err;
            bbbAvgMean += bbb_mean;
        }
        bbbAvgSqErr = bbbAvgSqErr / ((double) bincounter );
        bbbAvgSqErr = bbbAvgMean / ((double) bincounter );
        bbbAvgSqErr = TMath::Sqrt(bbbAvgSqErr); 
        
    }
     
     
    ///////////////////////////////////////////////////////////////////
    ////////////   S A V E   B E S T   V A L U E S   //////////////////
    ///////////////////////////////////////////////////////////////////    
 
    if ( flag_scan == 2) {
         
        // Create String to Print
        TString lineStrList =  thekey;
        lineStrList.Append(" ");
        for ( int i = lineStrList.Length() ; i < 100 ; i++ ) lineStrList.Append(" ");   
        lineStrList.Append(TString::Format("%.5f", optimalTauX)); 
        for ( int i = lineStrList.Length() ; i < 112 ; i++ ) lineStrList.Append(" ");    
        lineStrList.Append(TString::Format("%i", optimalTauLowK)); 
        for ( int i = lineStrList.Length() ; i < 116 ; i++ ) lineStrList.Append(" ");  
        lineStrList.Append(TString::Format("%i", optimalTauHighK)); 
        for ( int i = lineStrList.Length() ; i < 120 ; i++ ) lineStrList.Append(" ");  
        lineStrList.Append(TString::Format("%i", nSignSingularValues)); 
        for ( int i = lineStrList.Length() ; i < 124 ; i++ ) lineStrList.Append(" ");  
        lineStrList.Append(TString::Format("%i", nbins));
        
        // Print to File, with option 'append'
        SVD_LineToFile(lineStrList, outputfilenameTxt, "a"); 
        
    }
   
    ///////////////////////////////////////////////////////////////////
    ////////////   N O R M A L I Z A T I O N  /////////////////////////
    /////////////////////////////////////////////////////////////////// 
      
     
    // Get Probability Matrix
    TH2D* probMatrixHist = SVD_ProbMatrix(mcHist, xiniHist, numberSyst+1); 
    
     
    // Normalize BBB Unfolded histogram
    TH1D* normBBBHist = SVD_CloneHists1D(bbbHist, numberSyst+1);
    if ( flag_norm == 1 ) SVD_ExtNormalizeBBBDistribution(normBBBHist, bbbHist, globalEfficiency, globalEventYield, globalEventYieldErr, numberSyst+1);
    if ( flag_norm == 2 ) SVD_IntNormalizeBBBDistribution(normBBBHist, bbbHist, numberSyst+1);
    if ( flag_norm == 3 ) SVD_PseudatNormalizeBBBDistribution(normBBBHist, bbbHist, numberSyst+1);
    
    // Normalize SVD Unfolded histogram
    TH1D* normUnfHist = SVD_CloneHists1D(bbbHist, numberSyst+1);
    TH2D* normUnfCovHist = SVD_CloneHists2D(mcHist, 1); 
    if ( flag_regmode > 1 ) {
        if ( flag_norm == 1 ) SVD_ExtNormalizeSVDDistribution(normUnfHist, normUnfCovHist, unfHist, probMatrixHist, statCovHist, globalEfficiency, globalEventYield, globalEventYieldErr, numberSyst+1);
        if ( flag_norm == 2 ) SVD_IntNormalizeSVDDistribution(normUnfHist, normUnfCovHist, unfHist, statCovHist, numberSyst+1); 
        if ( flag_norm == 3 ) SVD_PseudatNormalizeSVDDistribution(normUnfHist, normUnfCovHist, unfHist, statCovHistNorm, numberSyst+1); 
    } else {
        delete[] normUnfHist;
        normUnfHist =  SVD_CloneHists1D(normBBBHist, numberSyst+1); 
        for ( int i=1; i<=nbins; i++ ) {
            for ( int j = 1; j <= nbins ; j++ ) {
                normUnfCovHist->SetBinContent(i,j,0.);
                normUnfCovHist->SetBinError(i,j,0.);
            } 
        } 
        for (int i=1; i<=nbins; i++ ) {
           double variance = TMath::Power(normBBBHist->GetBinError(i), 2.);
           normUnfCovHist->SetBinContent(i,i,variance);
        }   
    }
     
    // Normalize the SCALED generator level histogram
    TH1D* normGenHist = SVD_CloneHists1D(bbbHist, numberSyst+1);
    if ( flag_norm == 1 ) SVD_ExtNormalizeGenDistribution(normGenHist, xiniScaledHist, totalGenEvents, numberSyst+1);
    if ( flag_norm == 2 ) SVD_IntNormalizeGenDistribution(normGenHist, xiniScaledHist, numberSyst+1);
    if ( flag_norm == 3 ) SVD_PseudatNormalizeGenDistribution(normGenHist, xiniScaledHist, numberSyst+1);
      
      
     
    ///////////////////////////////////////////////////////////////////
    ////////////   R E T U R N   O B J E C T S  ///////////////////////
    ///////////////////////////////////////////////////////////////////    
 
    // Format everything that will be returned in all cases 
 
    // Strings
    if ( channel.CompareTo("") == 0 ) channel.Append("channel");
    if ( particle.CompareTo("") == 0 ) particle.Append("particle");
    if ( quantity.CompareTo("") == 0 ) quantity.Append("quantity");
    if ( special.CompareTo("") == 0 ) special.Append("");
    if ( syst.CompareTo("") == 0 ) syst.Append("");
    if ( channelTex.CompareTo("") == 0 ) channelTex.Append("channelTex");
    if ( particleTex.CompareTo("") == 0 ) particleTex.Append("particleTex");
    if ( quantityTex.CompareTo("") == 0 ) quantityTex.Append("quantityTex");
    if ( specialTex.CompareTo("") == 0 ) specialTex.Append("");
    if ( systTex.CompareTo("") == 0 ) systTex.Append("");
 
      
    // Format BBB Plot
    TString bbbStr = SVD_PlotName(channel, particle, quantity, special, syst, "BBB");
    SVD_SetTitles1D(bbbHist, bbbStr, quantityTex, "Events", numberSyst+1);

         
    // Format SVD Plot (Weights)
    TString weightStr = SVD_PlotName(channel, particle, quantity, special, syst, "WGT");
    SVD_SetTitles1D(weightHist, weightStr, quantityTex, "Weights", numberSyst+1);
        
         
    // Format SVD Plot
    TString unfStr = SVD_PlotName(channel, particle, quantity, special, syst, "UNF");
    SVD_SetTitles1D(unfHist, unfStr, quantityTex, "Events", numberSyst+1);
  
     
    // Format Normalized BBB Plot
    TString normBBBStr = SVD_PlotName(channel, particle, quantity, special, syst, "BBBnorm");
    SVD_SetTitles1D(normBBBHist, normBBBStr, quantityTex, "{#sigma}_{Bin} / {#sigma}_{tot}", numberSyst+1);

         
    // Format Normalized SVD Plot
    TString normUnfStr = SVD_PlotName(channel, particle, quantity, special, syst, "UNFnorm");
    SVD_SetTitles1D(normUnfHist, normUnfStr, quantityTex, "{#sigma}_{Bin} / {#sigma}_{tot}", numberSyst+1);
 
  
 
    ///////////////////////////////////////////////////////////////////
    ////////////   S A V E   T H E   O U T P U T    ///////////////////
    ///////////////////////////////////////////////////////////////////  
 

    // Prepare Output 
    unfolded = SVD_CloneHists1D(unfHist, numberSyst+1); 
    unfoldedNorm = SVD_CloneHists1D(normUnfHist, numberSyst+1);

   
    ///////////////////////////////////////////////////////////////////
    ////////////   F O R M A T   H I S T O G R A M S   ////////////////
    ///////////////////////////////////////////////////////////////////  
    
    // Format all histos we have already created!
     
    // Format Covariance Matrix
    TString dataCovStr = SVD_PlotName(channel, particle, quantity, special, syst, "DATACOV");
    SVD_SetTitles2D(dataCovHist, dataCovStr, quantityTex, quantityTex, "Covariance on Measurements", 1);  
        

    // Format Raw Data
    // Note: Here we have only 1 histogram in the array!
    TString rawStr = SVD_PlotName(channel, particle, quantity, special, syst, "YRAW");
    SVD_SetTitles1D(rawHist, rawStr, quantityTex, "Events", 1);

    // Format Bgr    
    TString bgrStr = SVD_PlotName(channel, particle, quantity, special, syst, "YBGR");
    SVD_SetTitles1D(bgrHist, bgrStr, quantityTex, "Events", numberSyst+1);


    // Format Data
    TString dataStr = SVD_PlotName(channel, particle, quantity, special, syst, "YDAT");
    SVD_SetTitles1D(dataHist, dataStr, quantityTex, "Events", numberSyst+1);


    // Format Scaled Data
    TString dataScaledStr = SVD_PlotName(channel, particle, quantity, special, syst, "YDATsc");
    SVD_SetTitles1D(dataScaledHist, dataScaledStr, quantityTex, "Events", numberSyst+1);


    // Format Response Matrix
    TString mcStr = SVD_PlotName(channel, particle, quantity, special, syst, "AA");
    TString mcStrX = quantityTex;
    mcStrX.Append(" (Reconstructed)");
    TString mcStrY = quantityTex;
    mcStrY.Append(" (Generated)");
    TString mcStrZ("Entries");
    SVD_SetTitles2D(mcHist, mcStr, mcStrX, mcStrY, mcStrZ, numberSyst+1);
 
  
    // Format Probability Matrix 
    TString probMatrixStr = SVD_PlotName(channel, particle, quantity, special, syst, "PRM"); 
    SVD_SetTitles2D(probMatrixHist, probMatrixStr, mcStrX, mcStrY, "Transition Probabilities in %", numberSyst+1);
    
    
    // Format Gen MC
    TString xiniStr = SVD_PlotName(channel, particle, quantity, special, syst, "XGEN");
    SVD_SetTitles1D(xiniHist, xiniStr, quantityTex, "Events", numberSyst+1);

 
    // Format Gen MC
    TString xiniScaledStr = SVD_PlotName(channel, particle, quantity, special, syst, "XGEN");
    SVD_SetTitles1D(xiniScaledHist, xiniScaledStr, quantityTex, "Events", numberSyst+1);

 
    // Format Normalized Gen MC
    TString normGenStr = SVD_PlotName(channel, particle, quantity, special, syst, "XGENnorm");
    SVD_SetTitles1D(normGenHist, normGenStr, quantityTex, "{#sigma}_{Bin} / {#sigma}_{tot}", numberSyst+1);
  

    // Format Rec MC
    TString biniStr = SVD_PlotName(channel, particle, quantity, special, syst, "XREC");
    SVD_SetTitles1D(biniHist, biniStr, quantityTex, "Events", numberSyst+1);


    // Format BBB Eff
    TString beffStr = SVD_PlotName(channel, particle, quantity, special, syst, "BBBEFF");        
    SVD_SetTitles1D(beffHist, beffStr, quantityTex, "Bin by Bin Efficiency", numberSyst+1);

        
    // Format Stat COV
    // Note: Here we have only 1 histogram in the array!
    TString statCovStr = SVD_PlotName(channel, particle, quantity, special, syst, "STATCOV");
    SVD_SetTitles2D(statCovHist, statCovStr, quantityTex, quantityTex, "Statistical Covariance (Absolute)", 1);


    // Format MC COV 
    // Note: Here we have only 1 histogram in the array!
    TString mcCovStr = SVD_PlotName(channel, particle, quantity, special, syst, "MCCOV");
    SVD_SetTitles2D(mcCovHist, mcCovStr, quantityTex, quantityTex, "Monte Carlo Covariance", 1);


    // Format Tot COV
    // Note: Here we have only 1 histogram in the array!
    TString totCovStr = SVD_PlotName(channel, particle, quantity, special, syst, "TOTCOV");
    SVD_SetTitles2D(totCovHist, totCovStr, quantityTex, quantityTex, "Total Covariance", 1);

        
    // Format Stat COV NORM 
    // Note: Here we have only 1 histogram in the array!
    TString normUnfCovStr = SVD_PlotName(channel, particle, quantity, special, syst, "STATCOVNORM");
    SVD_SetTitles2D(normUnfCovHist, normUnfCovStr, quantityTex, quantityTex, "Normalized Statistical Covariance", 1);
      

    // Format D-Distribution 
    TString ddStr = SVD_PlotName(channel, particle, quantity, special, syst, "DD");
    SVD_SetTitles1D(ddHist, ddStr, "Index i", "d-Value d_{i}", numberSyst+1);


    // Format Singular Values
    TString svStr = SVD_PlotName(channel, particle, quantity, special, syst, "SV");
    SVD_SetTitles1D(svHist, svStr, "Index i", "Singular Value s_{i}", numberSyst+1);
   
        
    /////////////////////////////////////////////////////////////////// 
    ////////////   M I G R A T I O N   P L O T S    /////////////////// 
    ///////////////////////////////////////////////////////////////////
    
    
    // SUM OVER REC BINS
    TH1D* sumOverRecHist = SVD_CloneHists1D(biniHist, numberSyst+1);
    SVD_SumOverRecBins(sumOverRecHist, mcHist, numberSyst+1);
    TString sumOverRecStr = SVD_PlotName(channel, particle, quantity, special, syst, "SUMREC");
    SVD_SetTitles1D(sumOverRecHist, sumOverRecStr, quantityTex, "#Sigma _{Rec}", numberSyst+1);
    
    
    // SUM OVER GEN BINS
    TH1D* sumOverGenHist = SVD_CloneHists1D(biniHist, numberSyst+1);
    SVD_SumOverGenBins(sumOverGenHist, mcHist, numberSyst+1);
    TString sumOverGenStr = SVD_PlotName(channel, particle, quantity, special, syst, "SUMGEN");
    SVD_SetTitles1D(sumOverGenHist, sumOverGenStr, quantityTex, "#Sigma _{Gen}", numberSyst+1);
    

    
    // TRANSITION PROBABILITY
    TH1D* probHist = SVD_ProbHist(mcHist, xiniHist, numberSyst+1);
    TString probStr = SVD_PlotName(channel, particle, quantity, special, syst, "PROB"); 
    SVD_SetTitles1D(probHist, probStr, quantityTex, "Trans. Prob.", numberSyst+1);
    
  
    // PURITY
    // Note: Purity is set to zero, if number of reconstructed
    // events in a bin is zero
    TH1D* purHist = SVD_CloneHists1D(biniHist, numberSyst+1);
    SVD_Pur(purHist, mcHist, numberSyst+1); 
    TString purStr = SVD_PlotName(channel, particle, quantity, special, syst, "PUR");
    SVD_SetTitles1D(purHist, purStr, quantityTex, "Purity", numberSyst+1);
 


    // STABILITY
    // Note: Stability is set to zero, if there are no
    // gerenated events in a bin
    // Note: If there is no reconstructed event in this bin,
    // then the stability is naturally zero, too
    TH1D* stabHist = SVD_CloneHists1D(biniHist, numberSyst+1); 
    SVD_Stab(stabHist, mcHist, numberSyst+1); 
    TString stabStr = SVD_PlotName(channel, particle, quantity, special, syst, "STAB");
    SVD_SetTitles1D(stabHist, stabStr, quantityTex, "Stability", numberSyst+1);
 
 

    // EFFICIENCY
    // Note: If no generated events are available in a bin,
    // then the efficiency is set to zero
    TH1D* effHist = SVD_CloneHists1D(biniHist, numberSyst+1); 
    SVD_Eff(effHist, mcHist, xiniHist, numberSyst+1);  
    TString effStr = SVD_PlotName(channel, particle, quantity, special, syst, "EFF"); 
    SVD_SetTitles1D(effHist, effStr, quantityTex, "Efficiency", numberSyst+1);
    
 

    // RATIO: Unfolded versus BBB
    TH1D* histRatioUnfBBB = SVD_MakeRatioNum(unfHist, bbbHist, numberSyst+1);
    TString ratioUnfBBBStr = SVD_PlotName(channel, particle, quantity, special, syst, "RATIOUNFBBB"); 
    SVD_SetTitles1D(histRatioUnfBBB, ratioUnfBBBStr, quantityTex, "Ratio", numberSyst+1);

 
 
         
    /////////////////////////////////////////////////////////////////// 
    ////////////   R E F O L D I N G     ////////////////////////////// 
    ///////////////////////////////////////////////////////////////////

 

    // REFOLDING - without Background 
    // Note, here we need to use the weightHist, since the
    // xiniHist contains empty gen level bins 
    TH1D* refoldHist = SVD_Refold(weightHist, NULL, mcHist, true, numberSyst+1); 
    TString refoldStr = SVD_PlotName(channel, particle, quantity, special, syst, "RFLD");
    SVD_SetTitles1D(refoldHist, refoldStr, "Bin i", "Refolded Value d_{i}", numberSyst+1);


    // RATIO: Refolded versus Data
    TH1D* histRatioRefDat = SVD_MakeRatioDenom(refoldHist, dataHist, numberSyst+1);
    TString ratioRefDatStr = SVD_PlotName(channel, particle, quantity, special, syst, "RATIOREFDAT"); 
    SVD_SetTitles1D(histRatioRefDat, ratioRefDatStr, quantityTex, "Ratio", numberSyst+1);



    // Calculate Chi2A 
    // For the nominal histo!
    double chi2A = SVD_Chi2A(dataHist, refoldHist);   
    int goodRecBins = SVD_GoodRecBins(dataHist);
    TString chi2AStr = TString::Format("#chi^{2} = %.3f", chi2A); 
    TString nDOFStr = TString::Format("nDOF=%i", goodRecBins);  
    
    
         
    /////////////////////////////////////////////////////////////////// 
    ////////////   E R R O R   P L O T S ////////////////////////////// 
    ///////////////////////////////////////////////////////////////////
 
 
    // These are done only for the nominal sample!

    // ERROR PLOT (Statistics)
    TString statErrStr = SVD_PlotName(channel, particle, quantity, special, syst, "STATERR");
    TH1D* statErrHist = SVD_Cov2Err(statCovHist, unfHist, statErrStr, quantityTex, "Statistical");

   
    // ERROR PLOT (MC)
    TString mcErrStr = SVD_PlotName(channel, particle, quantity, special, syst, "MCERR");
    TH1D* mcErrHist = SVD_Cov2Err(mcCovHist, unfHist, mcErrStr, quantityTex, "MC");



    // ERROR PLOT (Total)
    TString totErrStr = SVD_PlotName(channel, particle, quantity, special, syst, "TOTERR");
    TH1D* totErrHist = SVD_Cov2Err(totCovHist, unfHist, totErrStr, quantityTex, "Total");
 
 

    // ERROR PLOT (BBB)
    // Only for nominal sample!
    TH1D* bbbErrHist = SVD_HistoErrors1D(bbbHist, 1); 
    TString bbbErrStr = SVD_PlotName(channel, particle, quantity, special, syst, "BBBERR"); 
    SVD_SetTitles1D(bbbErrHist, bbbErrStr, quantityTex, "Statistical Error in %");

    
    // ERROR PLOT (Data Statistics, Normalized)
    // Only for nominal sample! 
    TH1D* normStatErrHist = SVD_HistoErrors1D(normUnfHist, 1);
    TString normStatErrStr = SVD_PlotName(channel, particle, quantity, special, syst, "STATERRnorm"); 
    SVD_SetTitles1D(normUnfHist, normStatErrStr, quantityTex, "Statistical Error on {#sigma}_{bin} / {#sigma}_{tot} in % (SVD)"); 
    
    
    // ERROR PLOT (BBB Data Statistics, Normalized)
    // Only for nominal sample!
    TH1D* normBBBErrHist = SVD_HistoErrors1D(normBBBHist, 1);
    TString normBBBErrStr = SVD_PlotName(channel, particle, quantity, special, syst, "BBBERRnorm"); 
    SVD_SetTitles1D(normBBBErrHist, normBBBErrStr, quantityTex, "Statistical Error on {#sigma}_{bin} / {#sigma}_{tot} in % (BBB)"); 
    

    // RATIO: UnfErrors versus BBBErrors
    TH1D* histRatioErrors = SVD_MakeRatioZero(statErrHist, bbbErrHist);
    TString ratioErrorsStr = SVD_PlotName(channel, particle, quantity, special, syst, "RATIOERROR"); 
    SVD_SetTitles1D(histRatioErrors, ratioErrorsStr, quantityTex, "Ratio");
    

    // RATIO: UnfErrors versus BBBErrors (Normalized)
    TH1D* normRatioErrors = SVD_MakeRatioZero(normStatErrHist, normBBBErrHist);
    TString rormRatioErrorsStr = SVD_PlotName(channel, particle, quantity, special, syst, "RATIONORMERROR"); 
    SVD_SetTitles1D(normRatioErrors, rormRatioErrorsStr, quantityTex, "Ratio");


    // CORRELATION STAT
    TString statCorrStr = SVD_PlotName(channel, particle, quantity, special, syst, "STATCORR");
    TH2D* statCorrHist = SVD_Cov2Corr(statCovHist, statCorrStr, quantityTex, "Statistical");


    // CORRELATION MC    
    TString mcCorrStr = SVD_PlotName(channel, particle, quantity, special, syst, "MCCORR");
    TH2D* mcCorrHist = SVD_Cov2Corr(mcCovHist, mcCorrStr, quantityTex, "MC");


    // CORRELATION TOT
    TString totCorrStr = SVD_PlotName(channel, particle, quantity, special, syst, "TOTCORR");
    TH2D* totCorrHist = SVD_Cov2Corr(totCovHist, totCorrStr, quantityTex, "Total");


    // CORRELATION STAT NORMALIZED
    TString normUnfCorrStr = SVD_PlotName(channel, particle, quantity, special, syst, "STATCORRORM");
    TH2D* normUnfCorrHist = SVD_Cov2Corr(normUnfCovHist, normUnfCorrStr, quantityTex, "Normalized Statistical");

 
    // GLOBAL CORRELATION
    // Attention: Do not use the side bins for the computation
    // of the global correlation !!!    
    // The global correlation in the side bins will be set 
    // to zero!
    TH1D* glcHist = SVD_CalcGlobCorr(statCovHist, dataHist);
    TString glcStr = SVD_PlotName(channel, particle, quantity, special, syst, "GLOBC");
    glcHist->SetName(glcStr);
    glcHist->SetTitle(glcStr);
    SVD_SetTitles1D(glcHist, glcStr, quantityTex, "Glob. Corr. in %");



    // Systematic Shifts with Unfolding
    TH1D* unfShiftHist = SVD_ArrayToShifts(unfHist, numberSyst+1);
    TString unfShiftHistStr  = SVD_PlotName(channel, particle, quantity, special, syst, "UNFSHIFT");
    SVD_SetTitles1D(unfShiftHist, unfShiftHistStr, quantityTex, "Syst. Shift in %", numberSyst);
    
    
    // Systematic Shifts with BBB
    TH1D* bbbShiftHist = SVD_ArrayToShifts(bbbHist, numberSyst+1);
    TString bbbShiftHistStr  = SVD_PlotName(channel, particle, quantity, special, syst, "BBBSHIFT");
    SVD_SetTitles1D(bbbShiftHist, bbbShiftHistStr, quantityTex, "Syst. Shift in %", numberSyst);
    
    
    // Comparison of systematic shifts
    TH1D* ratioShiftHist = SVD_MakeRatioZero(unfShiftHist, bbbShiftHist, numberSyst);
    TString ratioShiftHistStr = SVD_PlotName(channel, particle, quantity, special, syst, "RATIOSHIFT"); 
    SVD_SetTitles1D(ratioShiftHist, ratioShiftHistStr, quantityTex, "Ratio of Syst. Shifts", numberSyst); 
    
     
    // Systematic Shifts with Unfolding, Normalized
    TH1D* normUnfShiftHist = SVD_ArrayToShifts(normUnfHist, numberSyst+1);
    TString normUnfShiftHistStr  = SVD_PlotName(channel, particle, quantity, special, syst, "NORMUNFSHIFT");
    SVD_SetTitles1D(normUnfShiftHist, normUnfShiftHistStr, quantityTex, "Syst. Shift on {#sigma}_{bin} / {#sigma}_{tot} in %", numberSyst);
    
     
    // Systematic Shifts with BBB, Normalized
    TH1D* normBBBShiftHist = SVD_ArrayToShifts(normBBBHist, numberSyst+1);
    TString normBbbShiftHistStr  = SVD_PlotName(channel, particle, quantity, special, syst, "NORMBBBSHIFT");
    SVD_SetTitles1D(normBBBShiftHist, normBbbShiftHistStr, quantityTex, "Syst. Shift on {#sigma}_{bin} / {#sigma}_{tot} in %", numberSyst);
    
    
    


    ///////////////////////////////////////////////////////////
    ///////////     P S    F I L E   //////////////////////////
    /////////////////////////////////////////////////////////// 

 
 
    // DRAW IT ALL
    if ( flag_ps >= 2 ) { 
        
        // Setup Style
        #ifdef STANDALONE_COMPILE // needed for setHHStyle in dilepton analysis, where it is in a namespace
        common::setHHStyle(*gStyle);
        #else
        setHHStyle(*gStyle);
        #endif
        gStyle->SetPadTopMargin(0.15);
        gStyle->SetPadRightMargin(0.2);
        gStyle->SetPaintTextFormat(".0f");

        // Setup Canvas
        TCanvas* canvas = new TCanvas("Unfolding", "Unfolding");
        gROOT->ForceStyle();
        TGaxis::SetMaxDigits(3);
         
 
        // Make Output Path
        TString outputpath = "";
        bool makeFolder = false;
        if ( outputpath.CompareTo("") == 0 ) outputpath.Append(psFile);
        if ( outputpath.CompareTo("") == 0 ) outputpath.Append(rootFile);
        bool foundSlash = false;
        while ( foundSlash == false ) {
            if ( outputpath.EndsWith("/") == true ) {
                foundSlash = true;
                makeFolder = true;
            }
            if ( foundSlash == false ) outputpath.Chop();
        }
        if ( makeFolder == true) { 
            TString mkdirCommand = "mkdir -p ";
            mkdirCommand.Append(outputpath);
            mkdirCommand.Append(" 2> /dev/null");
            gSystem->Exec(mkdirCommand);
        } 

        // Tex / Channel, Particle and Quantity
        TString CPQTex = ""; 
        CPQTex.Append(channelTex); 
        CPQTex.Append(", "); 
        CPQTex.Append(particleTex);
        CPQTex.Append(", ");
        CPQTex.Append(quantityTex);     
        if(flag_plotStyle==2) CPQTex=channelTex;

        // Tex / Systematics
        TString SystTex = syst;
        SystTex = systTex;
        if(flag_plotStyle==2) SystTex="";

        // Tex / Regularization
        TString RegTex = ""; 
        if ( flag_regpar == 2 ) { 
            RegTex.Append("#tau = ");
            RegTex.Append(TString::Format("%.3f", regPar));  
        } else { 
            RegTex.Append("k = ");
            RegTex.Append(TString::Format("%.0f", regPar)); 
        } 
        if(flag_plotStyle==2) RegTex="";
    
        // Stacks and Legends
        THStack* theRegStack = NULL;
        THStack* theDValStack = NULL;
        TLegend* theLegend = NULL;

 
        // Draw Response Matrix
        SVD_Draw2D(mcHist, "COLZ TEXT", rangeL, rangeR, true);
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "InputEventMatrix", "OPEN");  


        // Draw Probability Matrix
        SVD_Draw2D(probMatrixHist, "COLZ TEXT", rangeL, rangeR, true);
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "InputResponseMatrix");
        


        // Draw Input Distributions
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, "");
        SVD_Array2Stack(theRegStack, theLegend, xiniHist, "Gen", "HIST E", "", 4, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, biniHist, "Rec", "HIST E", "", 2, numberSyst+1); 
        SVD_Array2Stack(theRegStack, theLegend, dataScaledHist, "Data, scaled", "HIST E", "", 1, numberSyst+1); 
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Entries", "", 0, true, rangeL, rangeR, true);   
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "InputDataGenRec");         
 
 

        // Draw Background related Distributions 
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, "");
        SVD_Array2Stack(theRegStack, theLegend, dataHist, "Data", "HIST E", "", 1, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, rawHist, "Raw", "HIST E", "", 2, 1); // Only 1 Histogram
        SVD_Array2Stack(theRegStack, theLegend, bgrHist, "Bgr", "HIST", "", 4, numberSyst+1); 
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Entries", "", 0, true, rangeL, rangeR, true);   
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "InputDataRawBgr");
 

        // Draw Weight distributions          
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, "");
        SVD_ArrayScale(weightHist, 1./lumiScaleFactor, numberSyst+1);
        TString weightHistLegendEntry = "Weights w_{i}";  
        SVD_Array2Stack(theRegStack, theLegend, weightHist, "Weights", "HIST", "", 1, numberSyst+1); 
        SVD_DrawStackAutoRange(theRegStack, theLegend, quantityTex, "Weights", "", 0, true, rangeL, rangeR, true);   
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultWeights");

  
        // Draw Unfolded distributions 
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex); 
        SVD_Array2Stack(theRegStack, theLegend, unfHist, "Unfolded", "HIST E", "", 1, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, dataHist, "Measured", "HIST E", "", 3, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, bbbHist, "BBB", "HIST E", "", 4, numberSyst+1); 
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Entries", "", 0, true, rangeL, rangeR, true);   
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultUnfBbbData");



        // Draw Unfolded versus MC prediction  
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex);  
        SVD_Array2Stack(theRegStack, theLegend, unfHist, "Unfolded", "HIST E", "", 1, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, xiniScaledHist, "Gen, sc.", "HIST E", "", 3, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, bbbHist, "BBB", "HIST E", "", 4, numberSyst+1); 
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Entries", "", 0, true, rangeL, rangeR, true);   
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultUnfBbbGen");
     
     
        // RATIO: Unfolded versus BBB 
        // Exclude OF bins here!
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex);  
        SVD_Array2Stack(theRegStack, theLegend, histRatioUnfBBB, "Unf/BBB", "HIST E", "", 1, numberSyst+1);  
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Ratio: SVD over BBB, Error from SVD", "", 0, false, rangeL, rangeR, true);  
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "RatioUnfBBB");
     
    
        // Draw Refolded Distribution
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex);
    if(flag_plotStyle==2){
        SVD_Array2Stack(theRegStack, theLegend, dataHist, "Measured", "HIST E", "", 1, numberSyst+1, 1);
        SVD_Array2Stack(theRegStack, theLegend, refoldHist, "Refolded", "HIST", "", 2, numberSyst+1, 2);
    }
    else{
        SVD_Array2Stack(theRegStack, theLegend, refoldHist, "Refolded", "HIST", "", 3, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, dataHist, "Measured", "HIST E", "", 1, numberSyst+1); 
    } 
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Entries", "", 0, true, rangeL, rangeR, true);
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultRefoldedData"); 
           
     
    
        // RATIO: Refolded versus Data 
        // Exclude OF bins here!
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex);  
        SVD_Array2Stack(theRegStack, theLegend, histRatioRefDat, "Ref/Data", "HIST E", "", 1, numberSyst+1);  
        SVD_DrawStackAutoRange(theRegStack, theLegend, quantityTex, "Ratio: Refolded over Data, Error from Data", "", 0, false, rangeL, rangeR, true);   
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "RatioRefoldedData");
    
      
        // Draw D-Plot  
        SVD_ClearStackLeg(theDValStack, theLegend, CPQTex, SystTex, "");  
        if ( flag_regmode > 1) gPad->SetLogy(true); 
        SVD_Array2Stack(theDValStack, theLegend, ddHist, "d-Values", "HIST", "", 1, numberSyst+1); 
        SVD_DrawStackAutoRange(theDValStack, theLegend, quantityTex, "d-Values", "", 0, true, gPad->GetLogy(), rangeL, rangeR, true);   
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "SvdD");
        gPad->SetLogy(false); 
   
   
        // Draw SV Plot  
        SVD_ClearStackLeg(theDValStack, theLegend, CPQTex, SystTex, "");  
        if ( flag_regmode > 1) gPad->SetLogy(true);
        SVD_Array2Stack(theDValStack, theLegend, svHist, "Sing. Values", "HIST", "", 1, numberSyst+1); 
        SVD_DrawStackAutoRange(theDValStack, theLegend, quantityTex, "Sing. Values", "", 0, true, gPad->GetLogy(), rangeL, rangeR, true);   
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "SvdSV");
        gPad->SetLogy(false); 
     
    
    
        // Draw Stab, Eff, Pur 
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, "");    
        SVD_Array2Stack(theRegStack, theLegend, purHist, "Pur.", "HIST", "", 3, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, stabHist, "Stab.", "HIST", "", 2, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, effHist, "Eff.", "HIST", "", 1, numberSyst+1); 
        SVD_Array2Stack(theRegStack, theLegend, probHist, "Prob. i #rightarrow i", "HIST", "", 4, numberSyst+1);
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Pur., Stab., Eff. in %", "", 0, true, rangeL, rangeR, true); 
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "MigrationEffStabPurr"); 
     
    
    
        // Smear-In and Smear-Out
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, ""); 
        SVD_Array2Stack(theRegStack, theLegend, lowerSmearinHist, "Smear-In Low", "HIST", "", 1, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, upperSmearinHist, "Smear-In Up", "HIST", "", 2, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, lowerSmearoutHist, "Smear-Out Low", "HIST", "", 3, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, upperSmearoutHist, "Smear-Out Up", "HIST", "", 4, numberSyst+1);
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Smear-In and Smear-Out", "", 0, true, rangeL, rangeR, true); 
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "MigrationSmear"); 
        
    
    
        // Draw BBB Efficiency 
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, "");   
        SVD_Array2Stack(theRegStack, theLegend, beffHist, "Gen. Eff.", "HIST", "", 2, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, effHist, "Efficiency", "HIST", "", 1, numberSyst+1); 
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Efficiency %", "", 0, true, rangeL, rangeR, true);  
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "MigrationEfficiency");
    
    
        // Covariance Matrix (DATA)
        // Only for the nominal sample
        SVD_Draw2D(dataCovHist, "COLZ TEXT", rangeL, rangeR, true);
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "InputCov");
    
    
        // Covariance Matrix (STAT)
        // Only for the nominal sample
        SVD_Draw2D(statCovHist, "COLZ TEXT", rangeL, rangeR, true);
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultCovStat");
    
    
        // Covariance Matrix (MC)
        // Only for the nominal sample
        SVD_Draw2D(mcCovHist, "COLZ TEXT", rangeL, rangeR, true);
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultCovMC");
    
    
        // Covariance Matrix (Total)
        // Only for the nominal sample
        SVD_Draw2D(totCovHist, "COLZ TEXT", rangeL, rangeR, true);
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultCovTot");
    
    
        // Covariance Matrix Normalized (STAT)
        // Only for the nominal sample
        SVD_Draw2D(normUnfCovHist, "COLZ TEXT", rangeL, rangeR, true);
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultNormCovStat");
     
    
        // Correlation Matrix (STAT)
        // Only for the nominal sample
        SVD_Draw2D(statCorrHist, "COLZ TEXT", rangeL, rangeR, true);
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultCorrStat");
    
    
        // Correlation Matrix (MC)
        // Only for the nominal sample
        SVD_Draw2D(mcCorrHist, "COLZ TEXT", rangeL, rangeR, true);
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultCorrMC");
    
    
        // Correlation Matrix (Total)
        // Only for the nominal sample
        SVD_Draw2D(totCorrHist, "COLZ TEXT", rangeL, rangeR, true);
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultCorrTot");
    
    
        // Correlation Matrix Normalized (STAT)
        // Only for the nominal sample
        SVD_Draw2D(normUnfCorrHist, "COLZ TEXT", rangeL, rangeR, true);
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultNormCorrStat"); 
         
    
     
        // GLOBAL CORRELATION STAT
        // Only for the nominal sample
        glcHist->UseCurrentStyle(); 
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex);   
        SVD_Array2Stack(theRegStack, theLegend, glcHist, "#rho_{i} in %", "HIST", "", 4, 1);  
        SVD_DrawStackAutoRange(theRegStack, theLegend, quantityTex, "Glob. Corr. #rho_{i} in %", "HIST ", 0, true, rangeL, rangeR, true);
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultGlobalCorr");
 
 
    
        // Errors (STAT, BBB)
        // Only for the nominal sample 
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex);  
        SVD_Array2Stack(theRegStack, theLegend, statErrHist, "Unf. Unc.", "HIST", "", 3, 1);  
        SVD_Array2Stack(theRegStack, theLegend, bbbErrHist, "BBB Unc.", "HIST", "", 4, 1);  
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Stat. Error in %", "HIST ", 0, true, rangeL, rangeR, true); 
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultErrorUnfBbb");
            
    
    
        // Errors (STAT, MC, TOTAL)
        // Only for the nominal sample
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex);   
        SVD_Array2Stack(theRegStack, theLegend, totErrHist, "Tot. Unc.", "HIST", "", 1, 1);  
        SVD_Array2Stack(theRegStack, theLegend, statErrHist, "Unf. Unc.", "HIST", "", 3, 1);  
        SVD_Array2Stack(theRegStack, theLegend, mcErrHist, "MC Unc.", "HIST", "", 4, 1);   
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "MC Stat. Error in %", "HIST ", 0, true, rangeL, rangeR, true); 
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultErrorStatMcTotal");
   

    
        // RATIO: UnfErrors versus BBBErrors
        // Only for the nominal sample
        // Do not show OF bins here!
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex);   
        SVD_Array2Stack(theRegStack, theLegend, histRatioErrors, "Err. Ratio", "HIST", "", 1, 1);  
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Err. Ratio", "", 0, false, rangeL, rangeR, true); 
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "RatioErrorUnfBbb");
    
          
        // Draw Unfolded versus MC prediction  
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex);  
        SVD_Array2Stack(theRegStack, theLegend, normUnfHist, "Unf for #sigma_{bin}/#sigma_{tot}", "HIST E", "", 1, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, normGenHist, "Gen for #sigma_{bin}/#sigma_{tot}", "HIST E", "", 3, numberSyst+1);
        SVD_Array2Stack(theRegStack, theLegend, normBBBHist, "BBB for #sigma_{bin}/#sigma_{tot}", "HIST E", "", 4, numberSyst+1); 
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Entries", "", 0, true, rangeL, rangeR, true);   
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultNormUnfBbbGen");
        
        
        // Errors (STAT, BBB)
        // Only for the nominal sample 
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex);  
        SVD_Array2Stack(theRegStack, theLegend, normStatErrHist, "Unf. Err. #sigma_{bin}/#sigma_{tot}", "HIST", "", 3, 1);  
        SVD_Array2Stack(theRegStack, theLegend, normBBBErrHist, "BBB Err. #sigma_{bin}/#sigma_{tot}", "HIST", "", 4, 1);  
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Error in %", "HIST ", 0, true, rangeL, rangeR, true); 
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultNormErrorUnfBbb");
   
 
    
        // RATIO: UnfErrors versus BBBErrors
        // Only for the nominal sample
        // Do not show OF bins here!
        SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex);   
        SVD_Array2Stack(theRegStack, theLegend, normRatioErrors, "Err. Rat. #sigma_{bin}/#sigma_{tot}", "HIST", "", 1, 1);  
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Err. Ratio", "", 0, false, rangeL, rangeR, true); 
        SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "RatioNormErrorUnfBbb");
            
              
        
        // Systematics Plots 
        if ( doSystematics == true ) {
         
            // Syst Shifts in Comparison
            SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex);   
            SVD_Array2Stack(theRegStack, theLegend, unfShiftHist, "Unfolding",  "HIST", "", 1, numberSyst);
            SVD_Array2Stack(theRegStack, theLegend, bbbShiftHist, "BBB", "HIST", "", 2, numberSyst);  
        SVD_DrawStackAutoRange(theRegStack, theLegend, quantityTex, "Syst. Shift in %", "", 0, false, rangeL, rangeR, true);  
            SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultShifts");
        
            // Shift Hist as Ratio
            SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex);   
            SVD_Array2Stack(theRegStack, theLegend, ratioShiftHist, "Shift Unf/BBB", "HIST", "", 1, numberSyst); 
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Syst. Shift Ratio", "", 0, false, rangeL, rangeR, true);  
            SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "RatioShifts");  
                
            // Syst Shifts in Comparison
            SVD_ClearStackLeg(theRegStack, theLegend, CPQTex, SystTex, RegTex);   
            SVD_Array2Stack(theRegStack, theLegend, normUnfShiftHist, "Unf.",  "HIST", "", 1, numberSyst);
            SVD_Array2Stack(theRegStack, theLegend, normBBBShiftHist, "BBB", "HIST", "", 2, numberSyst);  
        SVD_DrawStackAutoRange(theRegStack, theLegend, quantityTex, "Syst. Shift on #sigma_{bin} / #sigma_{tot} in %", "", 0, false, rangeL, rangeR, true);  
            SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "ResultNormShifts"); 
                 
        }  
            
    
        // K Scan Plots 
        // Only for the nominal sample
        if ( flag_ps >= 3 && flag_regmode > 1 ) {
         
             
            // Draw Unfolded distributions    
            SVD_ClearStackLeg(theRegStack, theLegend, "", "", "");
            SVD_FormatLegendKScan(theLegend);  
            for ( int tmpK = 2 ; tmpK <= nScanSingularValues ; tmpK++ ) {  
                TH1D* tmpUnfResult = SVD_CloneHists1D(arrK_Mean + tmpK - 2);
                SVD_Array2Stack(theRegStack, theLegend, tmpUnfResult, TString::Format("k=%i", tmpK), "HIST", "", tmpK, 1);  
                SVD_DeleteHists1D(tmpUnfResult);
            }
            SVD_Array2Stack(theRegStack, theLegend, unfHist, "Opt.", "P", "p", 1, 1);
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Unfolding Result x_{i}", "HIST ", 0, true, rangeL, rangeR, true);  
            SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "KScanResults"); 
            
             
            // Draw Error distributions   
            SVD_ClearStackLeg(theRegStack, theLegend, "", "", "");
            SVD_FormatLegendKScan(theLegend);   
            for ( int tmpK = 2 ; tmpK <= nScanSingularValues ; tmpK++ ) {  
                TH1D* tmpErr = SVD_CloneHists1D(arrK_StatError + tmpK - 2); 
                SVD_Array2Stack(theRegStack, theLegend, tmpErr, TString::Format("k=%i", tmpK), "HIST", "", tmpK, 1);  
                SVD_DeleteHists1D(tmpErr); 
            }
            SVD_Array2Stack(theRegStack, theLegend, statErrHist, "Opt.", "P", "p", 1, 1);
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Stat. Unc. #delta x_{i} in %", "HIST ", 0, true, rangeL, rangeR, true);  
            SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "KScanErrors"); 
         
             
            // Draw Glob Corr distributions   
            SVD_ClearStackLeg(theRegStack, theLegend, "", "", "");
            SVD_FormatLegendKScan(theLegend);    
            for ( int tmpK = 2 ; tmpK <= nScanSingularValues ; tmpK++ ) {  
                TH1D* tmpGlC = SVD_CloneHists1D(arrK_GlC + tmpK - 2); 
                SVD_Array2Stack(theRegStack, theLegend, tmpGlC, TString::Format("k=%i", tmpK), "HIST", "", tmpK, 1);     
                SVD_DeleteHists1D(tmpGlC); 
            }
            SVD_Array2Stack(theRegStack, theLegend, glcHist, "Opt.", "P", "p", 1, 1);
        SVD_DrawStackAutoRange(theRegStack, theLegend, quantityTex, "Global Correlation #rho_{i} in %", "HIST ", 0, true, rangeL, rangeR, true);  
            SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "KScanGlobCorr"); 
            
             
            // Draw Weight distributions   
            SVD_ClearStackLeg(theRegStack, theLegend, "", "", "");
            SVD_FormatLegendKScan(theLegend);     
            for ( int tmpK = 2 ; tmpK <= nScanSingularValues ; tmpK++ ) {  
                TH1D* tmpWeights = SVD_CloneHists1D(arrK_Weight + tmpK - 2); 
                SVD_Array2Stack(theRegStack, theLegend, tmpWeights, TString::Format("k=%i", tmpK), "HIST", "", tmpK, 1);    
                SVD_DeleteHists1D(tmpWeights); 
            } 
            SVD_Array2Stack(theRegStack, theLegend, weightHist, "Opt.", "P", "p", 1, 1);
        SVD_DrawStackZero(theRegStack, theLegend, quantityTex, "Weights w_{i}", "HIST ", 0, true, rangeL, rangeR, true);  
            SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "KScanWeights"); 
            
            
            // Delete all the obsolete objects
            SVD_DeleteHists1D(arrK_StatError, nScanSingularValues-1);
            SVD_DeleteHists1D(arrK_Weight, nScanSingularValues-1);
            SVD_DeleteHists1D(arrK_GlC, nScanSingularValues-1);
            SVD_DeleteHists1D(arrK_Mean, nScanSingularValues-1); 

        }
     
        // Scan Plots
        // Only for the nominal sample
        if ( flag_ps >= 4 && flag_regmode > 1) {
                  
            // Logscale
            gPad->SetLogx(true); 
                  
                  
            // Placeholders for Ranges       
            double xRangeMin = 0.;       
            double xRangeMax = 0.;
            double yRangeMin = 0.;
            double yRangeMax = 0.;     
            double newyRangeMin = 0.;
            double newyRangeMax = 0.;
            
            
            // Global objects for scan plots
            TH1D* theBgrHisto = NULL;
            TLatex* latexBestPoint = new TLatex();
            
                    
            // Global Correlation Plot
            canvas->Clear();
            delete theBgrHisto; theBgrHisto = NULL; 
            TGraph* gGlobCorr = SVD_Vect2Graph(vScanPoints, vGlobCorr);
            if(flag_plotStyle==2) SVD_DrawGraphAutoRange(gGlobCorr,  theBgrHisto, "P", 1, 0.2); 
            else                  SVD_DrawGraphAutoRange(gGlobCorr,  theBgrHisto, "P", 1); 
            TString gGlobCorrStr = SVD_PlotName(channel, particle, quantity, special, syst, "scanGLOBC");
            if(flag_plotStyle==2) SVD_SetTitles1D(theBgrHisto, gGlobCorrStr, "Parameter #tau", "RMS Global Correlation  #bar{#rho}  [%]");
            else                  SVD_SetTitles1D(theBgrHisto, gGlobCorrStr, "Parameter #tau", "Averaged Global Correlation  #bar{#rho}  in %");  
            TGraph* bestPoint = SVD_Point2Graph(optimalTauX, optimalTauY); 
            SVD_DrawGraph(bestPoint,  theBgrHisto, "P", 2); 
            bestPoint->SetMarkerSize(2.5); 
            bestPoint->SetMarkerStyle(3); 
            if(flag_plotStyle!=2) SVD_DrawVertLines(theBgrHisto, &vSingularValues, 2);
            TString textBestPoint = TString::Format("#tau = %.3f", optimalTauX);
            int textOrientationBestPoint = 23; 
            latexBestPoint->SetTextSize(0.03); 
            latexBestPoint->SetTextColor(2); 
            double textPosBestPointX =  SVD_TextPosX(theBgrHisto, optimalTauX, gPad->GetLogx(), textOrientationBestPoint);
            double textPosBestPointY =  SVD_TextPosY(theBgrHisto, optimalTauY, gPad->GetLogy(), textOrientationBestPoint);
            SVD_DrawText(textBestPoint, textPosBestPointX, textPosBestPointY, textOrientationBestPoint, latexBestPoint); 
            SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "TauScanGlobCorr");  
            delete gGlobCorr;
            delete bestPoint;  
                   
        
            // Chi Squared
            canvas->Clear();
            delete theBgrHisto; theBgrHisto = NULL; 
            TGraph* gChiSq = SVD_Vect2Graph(vScanPoints, vChiSq);
            SVD_DrawGraphZero(gChiSq,  theBgrHisto, "P", 4); 
            theBgrHisto->SetMaximum(2.*optimalChiSq); // Range Correction
            TString gChiSqStr = SVD_PlotName(channel, particle, quantity, special, syst, "scanCHI2A");
            SVD_SetTitles1D(theBgrHisto, gChiSqStr, "Parameter #tau", "#chi^{2}");
            bestPoint = SVD_Point2Graph(optimalTauX, optimalChiSq);
            SVD_DrawGraph(bestPoint,  theBgrHisto, "P", 2); 
            bestPoint->SetMarkerSize(2.5); 
            bestPoint->SetMarkerStyle(3);
            SVD_DrawVertLines(theBgrHisto, &vSingularValues, 2);
            textBestPoint = TString::Format("#tau = %.3f", optimalTauX);
            textOrientationBestPoint = 31; 
            latexBestPoint->SetTextSize(0.03);
            latexBestPoint->SetTextColor(2); 
            textPosBestPointX =  SVD_TextPosX(theBgrHisto, optimalTauX, gPad->GetLogx(), textOrientationBestPoint);
            textPosBestPointY =  SVD_TextPosY(theBgrHisto, optimalChiSq, gPad->GetLogy(), textOrientationBestPoint); 
            SVD_DrawText(textBestPoint, textPosBestPointX, textPosBestPointY, textOrientationBestPoint, latexBestPoint);  
            SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "TauScanChiSquare");
            delete gChiSq;
            delete bestPoint;  
                  
                  
            // Curvature
            canvas->Clear();
            delete theBgrHisto; theBgrHisto = NULL; 
            TGraph* gCurv = SVD_Vect2Graph(vScanPoints, vCurv);
            SVD_DrawGraphZero(gCurv,  theBgrHisto, "P", 7);  
            TString gCurvStr = SVD_PlotName(channel, particle, quantity, special, syst, "scanCURV");
            SVD_SetTitles1D(theBgrHisto, gCurvStr, "Parameter #tau", "Curvature x^{T}Cx");
            bestPoint = SVD_Point2Graph(optimalTauX, optimalCurv);
            SVD_DrawGraph(bestPoint,  theBgrHisto, "P", 2); 
            bestPoint->SetMarkerSize(2.5); 
            bestPoint->SetMarkerStyle(3);
            SVD_DrawVertLines(theBgrHisto, &vSingularValues, 2);
            textBestPoint = TString::Format("#tau = %.3f", optimalTauX);
            textOrientationBestPoint = 11; 
            latexBestPoint->SetTextSize(0.03);
            latexBestPoint->SetTextColor(2); 
            textPosBestPointX =  SVD_TextPosX(theBgrHisto, optimalTauX, gPad->GetLogx(), textOrientationBestPoint);
            textPosBestPointY =  SVD_TextPosY(theBgrHisto, optimalCurv, gPad->GetLogy(), textOrientationBestPoint); 
            SVD_DrawText(textBestPoint, textPosBestPointX, textPosBestPointY, textOrientationBestPoint, latexBestPoint);  
            SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "TauScanCurvature");
            delete gCurv;
            delete bestPoint;  
                  
                  
                  
                  
            // Avg Squared Error
            canvas->Clear();
            delete theBgrHisto; theBgrHisto = NULL; 
            TGraph* gAvgSqErr = SVD_Vect2Graph(vScanPoints, vAvgSqErr);
            SVD_DrawGraph(gAvgSqErr,  theBgrHisto, "P", 3); 
            SVD_RangeGraph(gAvgSqErr, xRangeMin, xRangeMax, yRangeMin, yRangeMax); 
            TString gAvgSqErrStr = SVD_PlotName(channel, particle, quantity, special, syst, "scanAVGERR");
            SVD_SetTitles1D(theBgrHisto, gAvgSqErrStr, "Parameter #tau", "Avg. Sq. Stat. Error");
            bestPoint = SVD_Point2Graph(optimalTauX, optimalAvgSqErr);
            SVD_DrawGraph(bestPoint,  theBgrHisto, "P", 2); 
            bestPoint->SetMarkerSize(2.5); 
            bestPoint->SetMarkerStyle(3);
            SVD_DrawVertLines(theBgrHisto, &vSingularValues, 2);
            SVD_DrawHorizLine(theBgrHisto, bbbAvgSqErr, 3);
            textBestPoint = TString::Format("#tau = %.3f", optimalTauX);
            textOrientationBestPoint = 11; 
            latexBestPoint->SetTextSize(0.03);
            latexBestPoint->SetTextColor(2); 
            textPosBestPointX =  SVD_TextPosX(theBgrHisto, optimalTauX, gPad->GetLogx(), textOrientationBestPoint);
            textPosBestPointY =  SVD_TextPosY(theBgrHisto, optimalAvgSqErr, gPad->GetLogy(), textOrientationBestPoint); 
            SVD_DrawText(textBestPoint, textPosBestPointX, textPosBestPointY, textOrientationBestPoint, latexBestPoint);  
            SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "TauScanAvgError");
            delete gAvgSqErr;
            delete bestPoint;  
                  
                  
                  
                  
            // Avg Mean
            canvas->Clear();
            delete theBgrHisto; theBgrHisto = NULL; 
            TGraph* gAvgMean = SVD_Vect2Graph(vScanPoints, vAvgMean);
            SVD_DrawGraph(gAvgMean,  theBgrHisto, "P", 6);
            SVD_RangeGraph(gAvgMean, xRangeMin, xRangeMax, yRangeMin, yRangeMax);
            SVD_NewRange(yRangeMin, yRangeMax, newyRangeMin, newyRangeMax, gPad->GetLogy()); 
            theBgrHisto->SetMaximum(newyRangeMax); 
            theBgrHisto->SetMinimum(newyRangeMin); 
            TString gAvgMeanStr = SVD_PlotName(channel, particle, quantity, special, syst, "scanAVGMEAN");
            SVD_SetTitles1D(theBgrHisto, gAvgMeanStr, "Parameter #tau", "Average Bin Content");
            bestPoint = SVD_Point2Graph(optimalTauX, optimalAvgMean);
            SVD_DrawGraph(bestPoint,  theBgrHisto, "P", 2); 
            theBgrHisto->SetMaximum(1.3*(*vAvgMean).Max()); // Range 
            theBgrHisto->SetMinimum(0.8*(*vAvgMean).Min()); // Range 
            bestPoint->SetMarkerSize(2.5); 
            bestPoint->SetMarkerStyle(3);
            SVD_DrawVertLines(theBgrHisto, &vSingularValues, 2);
            SVD_DrawHorizLine(theBgrHisto, bbbAvgMean, 3);
            textBestPoint = TString::Format("#tau = %.3f", optimalTauX);
            textOrientationBestPoint = 23; 
            latexBestPoint->SetTextSize(0.03);
            latexBestPoint->SetTextColor(2); 
            textPosBestPointX =  SVD_TextPosX(theBgrHisto, optimalTauX, gPad->GetLogx(), textOrientationBestPoint);
            textPosBestPointY =  SVD_TextPosY(theBgrHisto, optimalAvgMean, gPad->GetLogy(), textOrientationBestPoint); 
            SVD_DrawText(textBestPoint, textPosBestPointX, textPosBestPointY, textOrientationBestPoint, latexBestPoint);  
            SVD_PrintPage(canvas, outputfilenamePs, outputfilenameEps, "TauScanAvgMean");
            delete gAvgMean;
            delete bestPoint;   
                  
                   
             
            // Delete obsolete Objects
            delete latexBestPoint;     
            delete theBgrHisto;      
                  
            // Logscale
            gPad->SetLogx(false);
        } 
    
   
        // Last page empty
        canvas->Clear();
        SVD_PrintPage(canvas, outputfilenamePs, "", "", "CLOSE");
    
    
        // Reset Style
        #ifdef STANDALONE_COMPILE // needed for setHHStyle in dilepton analysis, where it is in a namespace
        common::setHHStyle(*gStyle);
        #else
        setHHStyle(*gStyle);
        #endif
        TGaxis::SetMaxDigits(2);
    
    
        // Delete Global Objects
        delete theRegStack;
        delete theDValStack;
        delete theLegend; 
            
    
            
        // Delete Canvas
        canvas->Clear();
        delete canvas;
        
    
    }
     
     
  
    //////////////////////////////////////////////////////////
    ///////////   SAVE TO TXT FILE  //////////////////////////
    //////////////////////////////////////////////////////////
    
    
    // Save all relevant Plots in ROOT File
    if ( writeHistToText ) {
    
     
        TString textOutputFolderName = "./SVD/Text";
        TString errSep = "+/-";
        TString lineSep = "\n";
        TString colSep = ", ";
    
        // Write histograms
        SVD_Hists1DToASCII(rawHist, textOutputFolderName, "%2.5f", errSep, lineSep, 1); 
        SVD_Hists1DToASCII(bgrHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(ttbgrHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(dataHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(dataScaledHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists2DToASCII(mcHist, textOutputFolderName, "%2.5f", errSep, lineSep, colSep, numberSyst+1); 
        SVD_Hists1DToASCII(xiniHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(xiniScaledHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(biniHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(beffHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1);  
        SVD_Hists1DToASCII(bbbHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists2DToASCII(dataCovHist, textOutputFolderName, "%2.5f", errSep, lineSep, colSep, 1); 
        SVD_Hists1DToASCII(weightHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(ddHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(svHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists2DToASCII(statCovHist, textOutputFolderName, "%2.5f", errSep, lineSep, colSep, 1);
        SVD_Hists2DToASCII(statCovHistNorm, textOutputFolderName, "%2.5f", errSep, lineSep, colSep, 1);
        SVD_Hists2DToASCII(mcCovHist, textOutputFolderName, "%2.5f", errSep, lineSep, colSep, 1); 
        SVD_Hists2DToASCII(totCovHist, textOutputFolderName, "%2.5f", errSep, lineSep, colSep, 1); 
        SVD_Hists1DToASCII(probHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(purHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(stabHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(effHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(histRatioUnfBBB, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(refoldHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(histRatioRefDat, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst+1); 
        SVD_Hists1DToASCII(statErrHist, textOutputFolderName, "%2.5f", errSep, lineSep, 1); 
        SVD_Hists1DToASCII(mcErrHist, textOutputFolderName, "%2.5f", errSep, lineSep, 1); 
        SVD_Hists1DToASCII(totErrHist, textOutputFolderName, "%2.5f", errSep, lineSep, 1); 
        SVD_Hists1DToASCII(bbbErrHist, textOutputFolderName, "%2.5f", errSep, lineSep, 1); 
        SVD_Hists1DToASCII(histRatioErrors, textOutputFolderName, "%2.5f", errSep, lineSep, 1); 
        SVD_Hists2DToASCII(statCorrHist, textOutputFolderName, "%2.5f", errSep, lineSep, colSep, 1); 
        SVD_Hists2DToASCII(normUnfCorrHist, textOutputFolderName, "%2.5f", errSep, lineSep, colSep, 1);  
        SVD_Hists2DToASCII(mcCorrHist, textOutputFolderName, "%2.5f", errSep, lineSep, colSep, 1); 
        SVD_Hists2DToASCII(totCorrHist, textOutputFolderName, "%2.5f", errSep, lineSep, colSep, 1); 
        SVD_Hists1DToASCII(glcHist, textOutputFolderName, "%2.5f", errSep, lineSep, 1); 
        SVD_Hists1DToASCII(ratioShiftHist, textOutputFolderName, "%2.5f", errSep, lineSep, numberSyst);  
 
    }  
  
    //////////////////////////////////////////////////////////
    ///////////   SAVE TO ROOT FILE //////////////////////////
    //////////////////////////////////////////////////////////


    // Save all relevant Plots in ROOT File
    if ( flag_root == 2 ) {
    
       
        // Open a ROOT file        
        TFile* file = new TFile(outputfilenameRoot, "RECREATE");
    
        // Write histograms
        SVD_WriteHists1D(rawHist, 1); 
        SVD_WriteHists1D(bgrHist, numberSyst+1); 
        SVD_WriteHists1D(ttbgrHist, numberSyst+1); 
        SVD_WriteHists1D(dataHist, numberSyst+1); 
        SVD_WriteHists1D(dataScaledHist, numberSyst+1); 
        SVD_WriteHists2D(mcHist, numberSyst+1); 
        SVD_WriteHists1D(xiniHist, numberSyst+1); 
        SVD_WriteHists1D(xiniScaledHist, numberSyst+1); 
        SVD_WriteHists1D(biniHist, numberSyst+1); 
        SVD_WriteHists1D(beffHist, numberSyst+1);  
        SVD_WriteHists1D(bbbHist, numberSyst+1); 
        SVD_WriteHists2D(dataCovHist, 1); 
        SVD_WriteHists1D(weightHist, numberSyst+1); 
        SVD_WriteHists1D(ddHist, numberSyst+1); 
        SVD_WriteHists1D(svHist, numberSyst+1); 
        SVD_WriteHists2D(statCovHist, 1); 
        SVD_WriteHists2D(mcCovHist, 1); 
        SVD_WriteHists2D(totCovHist, 1);
        SVD_WriteHists2D(normUnfCovHist, 1);  
        SVD_WriteHists1D(probHist, numberSyst+1); 
        SVD_WriteHists1D(purHist, numberSyst+1); 
        SVD_WriteHists1D(stabHist, numberSyst+1); 
        SVD_WriteHists1D(effHist, numberSyst+1); 
        SVD_WriteHists1D(histRatioUnfBBB); 
        SVD_WriteHists1D(refoldHist, numberSyst+1); 
        SVD_WriteHists1D(histRatioRefDat, numberSyst+1); 
        SVD_WriteHists1D(statErrHist, 1); 
        SVD_WriteHists1D(mcErrHist, 1); 
        SVD_WriteHists1D(totErrHist, 1); 
        SVD_WriteHists1D(bbbErrHist, 1); 
        SVD_WriteHists1D(histRatioErrors, 1); 
        SVD_WriteHists2D(statCorrHist, 1); 
        SVD_WriteHists2D(mcCorrHist, 1); 
        SVD_WriteHists2D(totCorrHist, 1); 
        SVD_WriteHists2D(normUnfCorrHist, 1);  
        SVD_WriteHists1D(glcHist, 1); 
        SVD_WriteHists1D(bbbShiftHist, numberSyst); 
        SVD_WriteHists1D(ratioShiftHist, numberSyst); 
    
        // Close file
        file->Close();
 
 
    }      
    
    //////////////////////////////////////////////////////////
    ///////////   DELETE OLD HISTOGRAMS  /////////////////////
    //////////////////////////////////////////////////////////
               
                 
    // Delete all old histograms
    SVD_DeleteHists1D(rawHist, 1); 
    SVD_DeleteHists1D(bgrHist, numberSyst+1); 
    SVD_DeleteHists1D(ttbgrHist, numberSyst+1); 
    SVD_DeleteHists1D(dataHist, numberSyst+1); 
    SVD_DeleteHists1D(dataHist_forTSVD, numberSyst+1); 
    SVD_DeleteHists1D(dataScaledHist, numberSyst+1); 
    SVD_DeleteHists2D(mcHist, numberSyst+1);  
    SVD_DeleteHists2D(mcHist_forTSVD, numberSyst+1);  
    SVD_DeleteHists1D(xiniHist, numberSyst+1); 
    SVD_DeleteHists1D(xiniHist_forTSVD, numberSyst+1); 
    SVD_DeleteHists1D(xiniScaledHist, numberSyst+1); 
    SVD_DeleteHists1D(biniHist, numberSyst+1); 
    SVD_DeleteHists1D(biniHist_forTSVD, numberSyst+1); 
    SVD_DeleteHists1D(beffHist, numberSyst+1);  
    SVD_DeleteHists1D(bbbHist, numberSyst+1); 
    SVD_DeleteHists2D(dataCovHist, 1); 
    SVD_DeleteHists1D(weightHist, numberSyst+1);  
    SVD_DeleteHists1D(ddHist, numberSyst+1); 
    SVD_DeleteHists1D(svHist, numberSyst+1); 
    SVD_DeleteHists2D(statCovHist, 1); 
    SVD_DeleteHists2D(mcCovHist, 1); 
    SVD_DeleteHists2D(totCovHist, 1); 
    SVD_DeleteHists2D(statCovHist, 1);  
    SVD_DeleteHists1D(probHist, numberSyst+1);  
    SVD_DeleteHists1D(purHist, numberSyst+1); 
    SVD_DeleteHists1D(stabHist, numberSyst+1); 
    SVD_DeleteHists1D(effHist, numberSyst+1); 
    SVD_DeleteHists1D(histRatioUnfBBB); 
    SVD_DeleteHists1D(refoldHist, numberSyst+1); 
    SVD_DeleteHists1D(histRatioRefDat, numberSyst+1); 
    SVD_DeleteHists1D(statErrHist, 1);  
    SVD_DeleteHists1D(mcErrHist, 1); 
    SVD_DeleteHists1D(totErrHist, 1); 
    SVD_DeleteHists1D(bbbErrHist, 1); 
    SVD_DeleteHists1D(histRatioErrors, 1); 
    SVD_DeleteHists2D(statCorrHist, 1); 
    SVD_DeleteHists2D(mcCorrHist, 1); 
    SVD_DeleteHists2D(totCorrHist, 1); 
    SVD_DeleteHists2D(statCorrHist, 1);   
    SVD_DeleteHists1D(glcHist, 1);  
    SVD_DeleteHists1D(bbbShiftHist, numberSyst);  
    SVD_DeleteHists1D(ratioShiftHist, numberSyst);  


    // Delete TopSVDUnfold-Objekt
    SVD_DeleteSVD(mySVDUnfold, numberSyst);
    
    
  
    //////////////////////////////////////////////////////////
    ///////////   RESTORE GLOBAL SETTINGS ////////////////////
    //////////////////////////////////////////////////////////
    

    // Restore gErrorIgnoreLv
    gErrorIgnoreLevel=save_gErrorLv;
    
    // Restore Graphics Modee 
    gROOT->SetBatch(save_isbatch);
     
    // Restore adding histos to TDirectory
    TH1::AddDirectory(save_addDir);
    
    // Restore handling of error structures
    TH1::SetDefaultSumw2(save_sumw2);


    // Screen Output
    if ( flag_verbose > 1 ) {
        cout << endl;
        cout << "********************************************************************************************************************" << endl;
        cout << "Unfolding finished." << endl;
        cout << "********************************************************************************************************************" << endl;
    }

    // return 
    return optimalTauX;
    
     
    
}
 

 
