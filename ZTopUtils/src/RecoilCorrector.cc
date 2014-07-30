#include "../interface/RecoilCorrector.h"

#include <fstream>

#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"


namespace ztop{

Double_t FitFuncGauss(Double_t* x, Double_t* par)
{
    Double_t a = (x[0]-par[1])/par[2];
    Double_t b = (x[0]-par[1])/par[3];
    Double_t c = (x[0]-par[1])/par[4];
    
    Double_t aLeft = (x[0]-par[1])/par[5];
    Double_t bLeft = (x[0]-par[1])/par[6];
    Double_t cLeft = (x[0]-par[1])/par[7];
    
    Double_t aGauss = par[8]*TMath::Exp(-0.5*a*a);
    Double_t bGauss = (1-par[8])*par[9]*TMath::Exp(-0.5*b*b);
    Double_t cGauss = (1-par[8])*(1-par[9])*TMath::Exp(-0.5*c*c);
    
    if (x[0]<par[1]) {
        aGauss = par[8]*TMath::Exp(-0.5*aLeft*aLeft);
        bGauss = (1-par[8])*par[9]*TMath::Exp(-0.5*bLeft*bLeft);
        cGauss = (1-par[8])*(1-par[9])*TMath::Exp(-0.5*cLeft*cLeft);
    }
    
    return par[0]*(aGauss+bGauss+cGauss);
}



RecoilCorrector::RecoilCorrector():
basedir_(""),
cannotCorrectMet_(0),
method_(-1),
dataset_(-1),
debug_(0)
{}



RecoilCorrector::RecoilCorrector(const std::string& basedir, const int method, const int dataset, const int debug):
basedir_(""),
cannotCorrectMet_(0),
method_(-1),
dataset_(-1),
debug_(0)
{
    basedir_ = basedir;
    method_ = method;
    dataset_ = dataset;
    debug_ = debug;
}



RecoilCorrector::~RecoilCorrector()
{
    delete metZParalData_;
    delete metZPerpData_;
    delete metZParalMC_;
    delete metZPerpMC_;
}



void RecoilCorrector::setFiles(const std::string& fileData, const std::string& fileMC)
{
    
    std::cout<<"\n--- Beginning preparation of Mva Met recoil correction"<<std::endl;
    
    std::string fileData_ = basedir_+fileData;
    std::string fileMC_ = basedir_+fileMC;
    
    std::ifstream file;
    file.open(fileData_.c_str(), std::ifstream::in);
    if (!file.is_open()) {
        std::cout<<"RecoilCorrector::SetFiles"<<std::endl;
        std::cout<<"The file '"<<fileData_<<"'doesn't exit"<<std::endl;
        std::cout<<"The MvaMet will __NOT BE__ corrected"<<std::endl;
        cannotCorrectMet_ = 1;
    }
    file.close();
    
    file.open(fileMC_.c_str(), std::ifstream::in);
    if (!file.is_open()) {
        std::cout<<"RecoilCorrector::SetFiles"<<std::endl;
        std::cout<<"The file '"<<fileMC_<<"'doesn't exit"<<std::endl;
        std::cout<<"The MvaMet will __NOT BE__ corrected"<<std::endl;
        cannotCorrectMet_ = 1;
    }
    file.close();
    
    if (cannotCorrectMet_) return;
    
    std::cout<<"Found all ROOT files needed for recoil correction"<<std::endl;
    
    data_ = new TFile(fileData_.c_str());
    mc_ = new TFile(fileMC_.c_str());
    
    std::cout<<"Calculating integrals of the different fits on Data and MC"<<std::endl;
    
    // hardcoded Z pt and njets bins to match the histograms in the ROOT files
    const int zptbins[] = {0, 10, 20, 30, 50, 1500};
    ZPtBins_ = std::vector<int>(zptbins, zptbins+sizeof(zptbins)/sizeof(zptbins[0]));
    const int njets[]= {0,1,2};
    nJets_ = std::vector<int>(njets, njets+sizeof(njets)/sizeof(njets[0]));
    
    metZParal_boundaries_vect_ = std::vector<std::vector<std::vector<double> > >(nZPtBins_ ,std::vector<std::vector<double> >( nJetsBins_, std::vector<double>(nbins_for_integral_, 0.)));
    metZPerp_boundaries_vect_ = std::vector<std::vector<std::vector<double> > >(nZPtBins_ ,std::vector<std::vector<double> >( nJetsBins_, std::vector<double>(nbins_for_integral_, 0.)));
    
    metZParalMC_integral_vect_ = std::vector<std::vector<std::vector<double> > >(nZPtBins_ ,std::vector<std::vector<double> >( nJetsBins_, std::vector<double>(nbins_for_integral_, 0.)));
    metZPerpMC_integral_vect_ = std::vector<std::vector<std::vector<double> > >(nZPtBins_ ,std::vector<std::vector<double> >( nJetsBins_, std::vector<double>(nbins_for_integral_, 0.)));
    metZParalData_integral_vect_ = std::vector<std::vector<std::vector<double> > >(nZPtBins_ ,std::vector<std::vector<double> >( nJetsBins_, std::vector<double>(nbins_for_integral_, 0.)));
    metZPerpData_integral_vect_ = std::vector<std::vector<std::vector<double> > >(nZPtBins_ ,std::vector<std::vector<double> >( nJetsBins_, std::vector<double>(nbins_for_integral_, 0.)));
    for (int ZPtBin=0; ZPtBin<nZPtBins_; ++ZPtBin) {
        for (int jetBin=0; jetBin<nJetsBins_; ++jetBin) {
            // this is needed to match the proper histogram names
            std::string lastBin = std::to_string(ZPtBins_.at(ZPtBin+1));
            if (ZPtBin == nZPtBins_-1) lastBin = std::string("Inf");
            
            std::string parallel = "metZParal"+std::to_string(ZPtBins_.at(ZPtBin))+"to"+lastBin+"_"+std::to_string(nJets_.at(jetBin))+"jetsH";
            metZParalData_[ZPtBin][jetBin] = (TF1*)data_->Get(parallel.c_str());
            metZParalMC_[ZPtBin][jetBin]   = (TF1*)mc_->Get(parallel.c_str());
            
            std::string perpendicular = "metZPerp"+std::to_string(ZPtBins_.at(ZPtBin))+"to"+lastBin+"_"+std::to_string(nJets_.at(jetBin))+"jetsH";
            metZPerpData_[ZPtBin][jetBin]  = (TF1*)data_->Get(perpendicular.c_str());
            metZPerpMC_[ZPtBin][jetBin]    = (TF1*)mc_->Get(perpendicular.c_str());
            
            // Get different constant parameters as mean and RMS if different parameters
            double xminD,xmaxD;
            // Met Paral Data
            metZParalData_[ZPtBin][jetBin]->GetRange(xminD,xmaxD);
            xminMetZParalData_[ZPtBin][jetBin] = float(xminD);
            xmaxMetZParalData_[ZPtBin][jetBin] = float(xmaxD);
            
            TF1 * func = getFuncRecoil(metZParalData_[ZPtBin][jetBin],true);
            rmsLeftMetZParalData_[ZPtBin][jetBin] = std::sqrt(func->CentralMoment(2,xminMetZParalData_[ZPtBin][jetBin],xmaxMetZParalData_[ZPtBin][jetBin]));
            delete func;
            
            func = getFuncRecoil(metZParalData_[ZPtBin][jetBin],false);
            rmsRightMetZParalData_[ZPtBin][jetBin] = std::sqrt(func->CentralMoment(2,xminMetZParalData_[ZPtBin][jetBin],xmaxMetZParalData_[ZPtBin][jetBin]));
            delete func;
            
            meanMetZParalData_[ZPtBin][jetBin] = metZParalData_[ZPtBin][jetBin]->GetParameter(1);
            rmsMetZParalData_[ZPtBin][jetBin] = std::sqrt(metZParalData_[ZPtBin][jetBin]->CentralMoment(2,xminMetZParalData_[ZPtBin][jetBin],xmaxMetZParalData_[ZPtBin][jetBin]));
            
            // Met Perp Data
            metZPerpData_[ZPtBin][jetBin]->GetRange(xminD,xmaxD);
            xminMetZPerpData_[ZPtBin][jetBin] = float(xminD);
            xmaxMetZPerpData_[ZPtBin][jetBin] = float(xmaxD);
            
            func = getFuncRecoil(metZPerpData_[ZPtBin][jetBin],true);
            rmsLeftMetZPerpData_[ZPtBin][jetBin] = std::sqrt(func->CentralMoment(2,xminMetZPerpData_[ZPtBin][jetBin],xmaxMetZPerpData_[ZPtBin][jetBin]));
            delete func;
            
            func = getFuncRecoil(metZPerpData_[ZPtBin][jetBin],false);
            rmsRightMetZPerpData_[ZPtBin][jetBin] = std::sqrt(func->CentralMoment(2,xminMetZPerpData_[ZPtBin][jetBin],xmaxMetZPerpData_[ZPtBin][jetBin]));
            delete func;
            
            meanMetZPerpData_[ZPtBin][jetBin] = metZPerpData_[ZPtBin][jetBin]->GetParameter(1);
            rmsMetZPerpData_[ZPtBin][jetBin] = std::sqrt(metZPerpData_[ZPtBin][jetBin]->CentralMoment(2,xminMetZPerpData_[ZPtBin][jetBin],xmaxMetZPerpData_[ZPtBin][jetBin]));
            
            // Met Paral MC
            metZParalMC_[ZPtBin][jetBin]->GetRange(xminD,xmaxD);
            xminMetZParalMC_[ZPtBin][jetBin] = float(xminD);
            xmaxMetZParalMC_[ZPtBin][jetBin] = float(xmaxD);
            
            func = getFuncRecoil(metZParalMC_[ZPtBin][jetBin],true);
            rmsLeftMetZParalMC_[ZPtBin][jetBin] = std::sqrt(func->CentralMoment(2,xminMetZParalMC_[ZPtBin][jetBin],xmaxMetZParalMC_[ZPtBin][jetBin]));
            delete func;
            
            func = getFuncRecoil(metZParalMC_[ZPtBin][jetBin],false);
            rmsRightMetZParalMC_[ZPtBin][jetBin] = std::sqrt(func->CentralMoment(2,xminMetZParalMC_[ZPtBin][jetBin],xmaxMetZParalMC_[ZPtBin][jetBin]));
            delete func;
            
            meanMetZParalMC_[ZPtBin][jetBin] = metZParalMC_[ZPtBin][jetBin]->GetParameter(1);
            rmsMetZParalMC_[ZPtBin][jetBin] = std::sqrt(metZParalMC_[ZPtBin][jetBin]->CentralMoment(2,xminMetZParalMC_[ZPtBin][jetBin],xmaxMetZParalMC_[ZPtBin][jetBin]));
            
            // Met Perp MC
            metZPerpMC_[ZPtBin][jetBin]->GetRange(xminD,xmaxD);
            xminMetZPerpMC_[ZPtBin][jetBin] = float(xminD);
            xmaxMetZPerpMC_[ZPtBin][jetBin] = float(xmaxD);
            
            func = getFuncRecoil(metZPerpMC_[ZPtBin][jetBin],true);
            rmsLeftMetZPerpMC_[ZPtBin][jetBin] = std::sqrt(func->CentralMoment(2,xminMetZPerpMC_[ZPtBin][jetBin],xmaxMetZPerpMC_[ZPtBin][jetBin]));
            delete func;
            
            func = getFuncRecoil(metZPerpMC_[ZPtBin][jetBin],false);
            rmsRightMetZPerpMC_[ZPtBin][jetBin] = std::sqrt(func->CentralMoment(2,xminMetZPerpMC_[ZPtBin][jetBin],xmaxMetZPerpMC_[ZPtBin][jetBin]));
            delete func;
            
            meanMetZPerpMC_[ZPtBin][jetBin] = metZPerpMC_[ZPtBin][jetBin]->GetParameter(1);
            rmsMetZPerpMC_[ZPtBin][jetBin] = std::sqrt(metZPerpMC_[ZPtBin][jetBin]->CentralMoment(2,xminMetZPerpMC_[ZPtBin][jetBin],xmaxMetZPerpMC_[ZPtBin][jetBin]));
            
            xminMetZParal_[ZPtBin][jetBin] = std::max(xminMetZParalData_[ZPtBin][jetBin],xminMetZParalMC_[ZPtBin][jetBin]);
            xmaxMetZParal_[ZPtBin][jetBin] = std::min(xmaxMetZParalData_[ZPtBin][jetBin],xmaxMetZParalMC_[ZPtBin][jetBin]);
            
            xminMetZPerp_[ZPtBin][jetBin] = std::max(xminMetZPerpData_[ZPtBin][jetBin],xminMetZPerpMC_[ZPtBin][jetBin]);
            xmaxMetZPerp_[ZPtBin][jetBin] = std::min(xmaxMetZPerpData_[ZPtBin][jetBin],xmaxMetZPerpMC_[ZPtBin][jetBin]);
            
            for (int iter = 0; iter<nbins_for_integral_; iter++) {
                
                // calculate the boundaries of the parallel component
                metZParal_boundaries_vect_[ZPtBin][jetBin][iter] = xminMetZParal_[ZPtBin][jetBin] + 1.* iter * (xmaxMetZParal_[ZPtBin][jetBin] - xminMetZParal_[ZPtBin][jetBin]) / nbins_for_integral_;
                
                // calculate the integral of the parallel component on MC
                metZParalMC_integral_vect_[ZPtBin][jetBin][iter] = metZParalMC_[ZPtBin][jetBin]->Integral(xminMetZParalMC_[ZPtBin][jetBin],metZParal_boundaries_vect_[ZPtBin][jetBin][iter], (const Double_t *)0, 1.e-5);
                
                // calculate the integral of the parallel component on Data
                metZParalData_integral_vect_[ZPtBin][jetBin][iter] = metZParalData_[ZPtBin][jetBin]->Integral(xminMetZParalData_[ZPtBin][jetBin],metZParal_boundaries_vect_[ZPtBin][jetBin][iter], (const Double_t *)0, 1.e-5);
                
                // calculate the boundaries of the perpendicular components
                metZPerp_boundaries_vect_[ZPtBin][jetBin][iter] = xminMetZPerp_[ZPtBin][jetBin] + 1.* iter * (xmaxMetZPerp_[ZPtBin][jetBin] - xminMetZPerp_[ZPtBin][jetBin]) / nbins_for_integral_;
                
                // calculate the integral of the perpendicular component on MC
                metZPerpMC_integral_vect_[ZPtBin][jetBin][iter] = metZPerpMC_[ZPtBin][jetBin]->Integral(xminMetZPerpMC_[ZPtBin][jetBin],metZPerp_boundaries_vect_[ZPtBin][jetBin][iter], (const Double_t *)0, 1.e-5);
                
                // calculate the integral of the perpendicular component on Data
                metZPerpData_integral_vect_[ZPtBin][jetBin][iter] = metZPerpData_[ZPtBin][jetBin]->Integral(xminMetZPerpData_[ZPtBin][jetBin],metZPerp_boundaries_vect_[ZPtBin][jetBin][iter], (const Double_t *)0, 1.e-5);
            }
        }
    }
    
    data_->Close(); delete data_;
    mc_->Close();   delete mc_;
    std::cout<<"=== Finish preparation of Mva Met recoil correction\n"<<std::endl;
}



float RecoilCorrector::correctMet(float& MetPx, float& MetPy,
                                  const float genZPx, const float genZPy,
                                  const float diLepPx, const float diLepPy,
                                  int njets)const
{
    // input parameters
    // MetPx, MetPy - missing transverse momentum 
    //                ( corrected replaces uncorrected )
    // genZPx, genZPy - generated transverse momentum of Z
    // diLepPx, diLepPy - dilepton transverse momentum 
    // njets - number of jets 
    // method : 2 - corrections by width w(MC)=w(Data)/w(MC) w(Process)
    // method : 3 - corrections by sampling 
    //              ( calculations of quantiles )
    
    if (cannotCorrectMet_) return 1;
    
    if (debug_) {
        std::cout<<"*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"<<std::endl;
        std::cout<<"Correcting the MET"<<std::endl;
    }
    
    float Zpt = std::sqrt(genZPx*genZPx + genZPy*genZPy);

    if (debug_>1) {
        std::cout<<"Reco Z pt = "<<Zpt<<std::endl;
    }
    
    float U1 = 0, U2 = 0;
    float metU1 = 0, metU2 = 0;
    
    this->CalculateU1U2FromMet(MetPx, MetPy, genZPx, genZPy, diLepPx, diLepPy, U1, U2, metU1, metU2);
    
    if (Zpt>1000) { Zpt = 999; }
    if (njets>=nJetsBins_) { njets = nJetsBins_-1;}
    
    int ZptBin = this->getBinNumber(Zpt, ZPtBins_);
    
    if (method_==3) {
        if (debug_) { std::cout<<"Correcting the MET by using the quantiles"<<std::endl; }
        
        if (U1>xminMetZParal_[ZptBin][njets] && U1<xmaxMetZParal_[ZptBin][njets]) {
            size_t u1bin = this->getBinFromVector(metZParal_boundaries_vect_.at(ZptBin).at(njets), U1);
            double sumProb[1] = {metZParalMC_integral_vect_.at(ZptBin).at(njets).at(u1bin)};
            
            if (sumProb[0]<0) { sumProb[0] = 1e-4; }
            if (sumProb[0]>1) { sumProb[0] = 0.9999; }
            
            const size_t bin = this->getBinFromVector(metZParalData_integral_vect_.at(ZptBin).at(njets), sumProb[0]);
            const float U1reco = metZParal_boundaries_vect_.at(ZptBin).at(njets).at(bin);
            
            if (U1reco>xminMetZParal_[ZptBin][njets] && U1reco<xmaxMetZParal_[ZptBin][njets]) {
                U1 = U1reco;
            }
        }
        
        if (U2>xminMetZPerp_[ZptBin][njets] && U2<xmaxMetZPerp_[ZptBin][njets]) {
            size_t u2bin = this->getBinFromVector(metZPerp_boundaries_vect_.at(ZptBin).at(njets), U2);
            double sumProb[1] = {metZPerpMC_integral_vect_.at(ZptBin).at(njets).at(u2bin)};
            
            if (sumProb[0]<0) { sumProb[0] = 1e-4; }
            if (sumProb[0]>1) { sumProb[0] = 0.9999; }
            
            const size_t bin = this->getBinFromVector(metZPerpData_integral_vect_.at(ZptBin).at(njets), sumProb[0]);
            const float U2reco = metZPerp_boundaries_vect_.at(ZptBin).at(njets).at(bin);
            
            if (U2reco>xminMetZPerp_[ZptBin][njets] && U2reco<xmaxMetZPerp_[ZptBin][njets]) {
                U2 = U2reco;
            }
        }
    }
    else {
        if (debug_) { std::cout<<"Correcting the MET by using the width of the mean values of the MET in data and MC"<<std::endl; }
        this->U1U2CorrectionsByWidth(U1, U2, ZptBin, njets);
    }

    this->CalculateMetFromU1U2(U1, U2, genZPx, genZPy, diLepPx, diLepPy, MetPx, MetPy);

    float weight = 1;
    if (method_==3) {
        const float MetX = std::sqrt(MetPx*MetPx+MetPy*MetPy);
        if (dataset_==0) {
        // Run2011A SingleMu
            if (njets==0 && MetX>52 && MetX<=60) {
                weight *= 0.634;
            }
            if (njets==2 && MetX>40 && MetX<=60) {
            weight *= 0.752;
            }
        }

        if (dataset_==1) {
        // Run2011A DoubleMu
            if (njets==0 && MetX>52 && MetX<=64) {
                weight *= 0.706;
            }
        }

        if (dataset_==2) {
        // Run2011B
            if (njets==0) {
                if (MetX>56&&MetX<=72) {
                    weight *= 1.271;
                }
                if (MetX>72&&MetX<=80) {
                    weight *= 6.236;
                }
            }
            if (njets==1 && MetX>56&&MetX<=72) {
                weight *= 1.596;
            }
            if (njets==2 && MetX>68&&MetX<=72) {
                weight *= 0.1;
            }
        }
    }
    return weight;
}



size_t RecoilCorrector::getBinFromVector(const std::vector<double>& vector_of_boundaries, const float value)const
{
    size_t bin = std::lower_bound(vector_of_boundaries.begin(), vector_of_boundaries.end(), value) - vector_of_boundaries.begin();
    if (bin == vector_of_boundaries.size()) return --bin;
    
    return bin;
}



void RecoilCorrector::U1U2CorrectionsByWidth(float& U1, float& U2, const int ZptBin, int njets)const
{
    if (njets>=nJetsBins_) {njets = nJetsBins_-1;}
    
    float width = U1 - meanMetZParalMC_[ZptBin][njets];
    
    if (width<0) {
        width *= rmsLeftMetZParalData_[ZptBin][njets]/rmsLeftMetZParalMC_[ZptBin][njets];
    } else {
        width *= rmsRightMetZParalData_[ZptBin][njets]/rmsRightMetZParalMC_[ZptBin][njets];
    }
    U1 = meanMetZParalData_[ZptBin][njets] + width;

    // ********* U2 *************
    width = U2 - meanMetZPerpMC_[ZptBin][njets];
    
    if (width<0) {
        width *= rmsLeftMetZPerpData_[ZptBin][njets]/rmsLeftMetZPerpMC_[ZptBin][njets];
    } else {
        width *= rmsRightMetZPerpData_[ZptBin][njets]/rmsRightMetZPerpMC_[ZptBin][njets];
    }
    U2 = meanMetZPerpData_[ZptBin][njets] + width;
}



void RecoilCorrector::CalculateU1U2FromMet(const float metPx, const float metPy,
                                           const float genZPx, const float genZPy,
                                           const float diLepPx, const float diLepPy,
                                           float& U1, float& U2,
                                           float& metU1, float& metU2)const
{
    if (debug_) { std::cout<<"Calculating the U1 and U2 components from the MET"<<std::endl; }
    
    const float hadRecX = metPx + diLepPx - genZPx;
    const float hadRecY = metPy + diLepPy - genZPy;
    
    const float hadRecPt = std::sqrt(hadRecX*hadRecX+hadRecY*hadRecY);
    const float phiHadRec = std::atan2(hadRecY,hadRecX);
    
    const float phiDiMuon = std::atan2(diLepPy,diLepPx);
    const float phiMEt = std::atan2(metPy,metPx);
    
    const float metPt = std::sqrt(metPx*metPx+metPy*metPy);
    const float phiZ = std::atan2(genZPy,genZPx);
    
    U1 = hadRecPt * std::cos(phiHadRec - phiZ);
    U2 = hadRecPt * std::sin(phiHadRec - phiZ);
    
    if (debug_>1) {
        std::cout<<"Parallel component of Zpt, U1 = "<<U1<<std::endl;
        std::cout<<"Orthogonal component of Zpt, U2 = "<<U2<<std::endl;
    }
    
    metU1 = metPt * std::cos(phiMEt - phiDiMuon);
    metU2 = metPt * std::sin(phiMEt - phiDiMuon);
    if (debug_>1) {
        std::cout<<"Parallel component of MET, U1 = "<<metU1<<std::endl;
        std::cout<<"Orthogonal component of MET, U2 = "<<metU2<<std::endl;
    }
    
}



void RecoilCorrector::CalculateMetFromU1U2(const float U1, const float U2,
                                           const float genZPx, const float genZPy,
                                           const float diLepPx, const float diLepPy,
                                           float& metPx, float& metPy)const
{
    if (debug_) { std::cout<<"Calculating the U1 and U2 components from the MET"<<std::endl; }
    
    const float hadRecPt = std::sqrt(U1*U1+U2*U2);
    
    const float deltaPhiZHadRec = std::atan2(U2,U1);
    const float phiZ = std::atan2(genZPy,genZPx);
    
    const float phiHadRec = phiZ + deltaPhiZHadRec;
    
    const float hadRecX = hadRecPt*std::cos(phiHadRec);
    const float hadRecY = hadRecPt*std::sin(phiHadRec);
    
    metPx = hadRecX + genZPx - diLepPx;
    metPy = hadRecY + genZPy - diLepPy;
    if (debug_>1) {
        std::cout<<"x component of MET, px = "<<metPx<<std::endl;
        std::cout<<"y component of MET, py = "<<metPy<<std::endl;
    }
}



int RecoilCorrector::getBinNumber(const float x, const std::vector<int>& bins)const
{
    for (size_t iB=0; iB<bins.size(); ++iB) {
        if (x>=1.*bins.at(iB) && x<1.*bins.at(iB+1)) return iB;
    }
    return 0;
}



TF1* RecoilCorrector::getFuncRecoil(const TF1* const initFunc, const bool left)const
{
    double xminD;
    double xmaxD;

    initFunc->GetRange(xminD,xmaxD);

    float xmin = float(xminD);
    float xmax = float(xmaxD);

    TF1 * func = new TF1("func",FitFuncGauss,xmin,xmax,10);
    func->SetParameter(0,initFunc->GetParameter(0));
    func->SetParameter(1,initFunc->GetParameter(1));
    if (left) {
        func->SetParameter(2,initFunc->GetParameter(5));
        func->SetParameter(3,initFunc->GetParameter(6));
        func->SetParameter(4,initFunc->GetParameter(7));
        func->SetParameter(5,initFunc->GetParameter(5));
        func->SetParameter(6,initFunc->GetParameter(6));
        func->SetParameter(7,initFunc->GetParameter(7));
    } else {
        func->SetParameter(2,initFunc->GetParameter(2));
        func->SetParameter(3,initFunc->GetParameter(3));
        func->SetParameter(4,initFunc->GetParameter(4));
        func->SetParameter(5,initFunc->GetParameter(2));
        func->SetParameter(6,initFunc->GetParameter(3));
        func->SetParameter(7,initFunc->GetParameter(4));
    }
    func->SetParameter(8,initFunc->GetParameter(8));
    func->SetParameter(9,initFunc->GetParameter(9));

    return func;
}



}
