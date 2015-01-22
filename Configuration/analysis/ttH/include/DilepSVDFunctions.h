#ifndef DILEPSVDFUNCTIONS_H
#define DILEPSVDFUNCTIONS_H

// Include parent class
#include "../../unfolding/TopSVDFunctions.h"





class DilepSVDFunctions : public TopSVDFunctions
{

public:
    DilepSVDFunctions();
    
    // Setup
    void SetOutputPath(TString path) { outpath = path;}


    // Helpers 
    TString SVD_GetSteering(TString channel, TString particle, TString quantity, TString special, TString syst);
    void SVD_Tex(TString channel, TString particle, TString quantity, TString special, TString& channelTex, TString& particleTex, TString& quantityTex, TString& specialTex);
    TString SVD_GetOutputPath();
    TString SVD_GetOutputFileNamePs(TString channel, TString particle, TString quantity, TString syst, TString special);
    TString SVD_GetOutputFileNameEps(TString channel, TString particle, TString quantity, TString syst, TString special);
    TString SVD_GetOutputFileNameRoot(TString channel, TString particle, TString quantity, TString syst, TString special); 
    TString SVD_GetOutputFileNameTxt(); 
    
    // Unfolding
    double SVD_DoUnfold(
        TH1D* dataInputHist, 
        TH1D* bgrInputHist, 
        TH1D* ttbgrInputHist, 
        TH1D* genInputHist, 
        TH1D* recInputHist, 
        TH2D* respInputHist,
        double* totalDataEvents,
        double* totalBgrEvents, 
        double* totalTtBgrEvents, 
        double* totalRecEvents,
        double* totalGenEvents,
        const double thebins[], 
        const int numbins, 
        TH1D*& unfolded, 
        TH1D*& unfoldedNorm, 
        const int numSys, 
        TString channel,
        TString particle,
        TString quantity,
        TString special,
        TString syst);
        
    // Unfolding with Systematics
    double SVD_DoUnfoldSys(
        TH1D* dataInputHist,        
        TH1D* bgrNom, TH1D* bgrUp, TH1D* bgrDown,  
        TH1D* ttbgrNom, TH1D* ttbgrUp, TH1D* ttbgrDown,       
        TH1D* genNom, TH1D* genUp, TH1D* genDown,         
        TH1D* recNom, TH1D* recUp, TH1D* recDown,      
        TH2D* respNom, TH2D* respUp, TH2D* respDown,  
        double totalDataEventsNom,  
        double totalBgrEventsNom,    double totalBgrEventsUp,    double totalBgrEventsDown, 
        double totalTtBgrEventsNom,  double totalTtBgrEventsUp,  double totalTtBgrEventsDown, 
        double totalRecEventsNom,    double totalRecEventsUp,    double totalRecEventsDown, 
        double totalGenEventsNom,    double totalGenEventsUp,    double totalGenEventsDown,  
        const double thebins[],                 
        const int numbins,                       
        TH1D*& shifts,                        
        TString channel,                       
        TString particle,                      
        TString quantity,                      
        TString special,
        TString syst);
            
 
 
private: 
    TString outpath;
     
     
};

#endif

