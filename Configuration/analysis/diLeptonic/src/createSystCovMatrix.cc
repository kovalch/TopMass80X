#include <TSystem.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>

#include "fstream"
#include "iostream"
#include "iomanip"

const TString basedir = TString("Plots/");
const TString baseOutDir = TString("correlationMatrix/");

std::vector<TString> createVecFiles(TString channel, const std::vector<TString>systematics)
{
    std::vector<TString> result;

    for(auto syst:systematics){
        TString tmp_file;
        if(!syst.Contains("HAD_")){ tmp_file = basedir.Copy().Append(syst+"/"+channel+"/DiffXS_");}
        else{
            if(syst.Contains("HAD_UP")){ tmp_file = basedir.Copy().Append("MCATNLO/"+channel+"/DiffXS_");}
            else if(syst.Contains("HAD_DOWN")){ tmp_file = basedir.Copy().Append("POWHEG/"+channel+"/DiffXS_");}
        }
        result.push_back(tmp_file);
    }
    return result;
}


TH1* getHistoFromFile(TString filename, TString variable)
{
    TString histogram = TString(variable);
    TFile *f = TFile::Open(filename);
    if(!filename) {
        std::cout<<"File '"<<filename<<"' not valid"<<std::endl;
        return nullptr;
    };

    TH1* histo = dynamic_cast<TH1*>(f->Get(histogram));
    histo->SetDirectory(0);
    if(!histo){
        std::cout<<"Histogram '"<<histogram<<"' not valid"<<std::endl;
        return nullptr;
    };
    f->Close();
    delete f;
    return histo;
}

TGraphAsymmErrors* getGraphFromFile(TString filename, TString variable)
{
    TString histogram = TString(variable);
    TFile *f = TFile::Open(filename);
    if(!filename) {
        std::cout<<"File '"<<filename<<"' not valid"<<std::endl;
        return nullptr;
    };
    TGraphAsymmErrors* gr = dynamic_cast<TGraphAsymmErrors*>(f->Get(histogram));
    if(!gr){
        std::cout<<"Histogram '"<<gr<<"' not valid"<<std::endl;
        return nullptr;
    };
    f->Close();
    delete f;
    return gr;
}




TH1* convertGraphToHisto(TH1 *histo , TGraphAsymmErrors* graph)
{
    TH1 *htmp = dynamic_cast<TH1*>(histo->Clone("htmp"));
    for (Int_t iter=0; iter<(Int_t)graph->GetN(); iter++)
    {
        int bin = histo->FindBin(graph->GetX()[iter]);
        htmp->SetBinContent(bin, graph->GetY()[iter]);
    }
    htmp->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    return htmp;
}

/// Calculate the binwise difference between 2 histograms: h_final (i) = h0(i)- h1(i)
/// Please note that if 'absValue == true' the absolute value of the difference is returned
TH1* calculateDifference(TH1* h0, TH1* h1, bool absValue)
{
    TH1* htmp = dynamic_cast<TH1*> (h0->Clone("htmp"));
    for(auto iter=0; iter<= h1->GetNbinsX()+1; iter++)
    {
        double value = h0->GetBinContent(iter) - h1->GetBinContent(iter);
        if(absValue){htmp->SetBinContent(iter, std::fabs(value));}
        else {htmp->SetBinContent(iter, value);};
    }
    return htmp;
}

/// This function returns the value of: hCent + factor *symmError
/// where symmError = 0.5*(abs(hUp - hCent) + abs(hDo - hCent)) * sign(hUp - hCent)

TH1* getSymmError(TH1* hCent, TH1* hUp, TH1* hDo, int sign)
{
    if(std::abs(sign) != 1) {std::cout<<"getSymmError: factor not '1' or '-1'. Not valid"<<std::endl; exit(14);};

    TH1* tmp_Cent = dynamic_cast<TH1*>(hCent->Clone("tmp_Cent"));
    TH1* tmp_Up = dynamic_cast<TH1*>(hUp->Clone("tmp_Up"));
    TH1* tmp_Do = dynamic_cast<TH1*>(hDo->Clone("tmp_Do"));

    for (auto it2 =0; it2<= 1+tmp_Cent->GetNbinsX(); it2++){
        /// get values from histograms
        double valCent = tmp_Cent->GetBinContent(it2);
        double valUp = tmp_Up->GetBinContent(it2);
        double valDo = tmp_Do->GetBinContent(it2);

        /// calculate variations
        double diffUp = valUp - valCent;
        double diffDo = valDo - valCent;
        double factor = (diffUp>=0) ? 1 : -1;
        double totalDiff = factor * 0.5 * (std::fabs(diffUp) + std::fabs(diffDo));
        /// fill new variations into histogram
        tmp_Up->SetBinContent(it2, valCent + sign * totalDiff);
        tmp_Do->SetBinContent(it2, valCent - sign * totalDiff);
    }
    delete tmp_Cent;
    delete tmp_Do;
    return  tmp_Up;
}

const TH2D combineStatCovMatrices(const TH2D ee_StatCov, const double ee_statError[], const TH2D emu_StatCov, const double emu_statError[], const TH2D mumu_StatCov, const double mumu_statError[])
{
    TH2D result_h = ee_StatCov;
//     TH2D *result_h = (TH2D*) ee_StatCov.Clone("result");
    const int nbins = ee_StatCov.GetNbinsX();
    
    for (int iterX = 1; iterX<nbins+1; iterX++){
        double ee_Stat = ee_statError[iterX-1];
        double emu_Stat = emu_statError[iterX-1];
        double mumu_Stat = mumu_statError[iterX-1];
        
        for (int iterY = 1; iterY<nbins+1; iterY++){
            const double ee_Value = ee_StatCov.GetBinContent(iterX, iterY);
            const double emu_Value = emu_StatCov.GetBinContent(iterX, iterY);
            const double mumu_Value = mumu_StatCov.GetBinContent(iterX, iterY);
        
            const double numerator = ee_Value/(ee_Stat*ee_Stat) + emu_Value/(emu_Stat*emu_Stat) + mumu_Value/(mumu_Stat*mumu_Stat);
            const double denominator = 1./(ee_Stat*ee_Stat) + 1./(emu_Stat*emu_Stat) + 1./(mumu_Stat*mumu_Stat);
            
            result_h.SetBinContent(iterX, iterY, numerator/denominator);
        }
    }
    return result_h;
}


const TH2D statCovValues(TString variable, TString channel)
{
    TString particle_ = "", variable_ = "";

    if (variable.Contains("Eta"))           { variable_ = "Eta";}
    else if (variable.Contains("Rapidity")) { variable_ = "Rapidity";}
    else if (variable.Contains("Mass"))     { variable_ = "Mass";}
    else if (variable.Contains("pT"))       { variable_ = "Pt";}
    else if (variable.Contains("Mult"))     { variable_ = "Mult";}

    if (variable.Contains("HypLepton"))         { particle_ = "Leptons";}
    else if (variable.Contains("HypLeptonBjet")){ particle_ = "Leptons";}
    else if (variable.Contains("HypLLBar"))     { particle_ = "LepPair";}
    else if (variable.Contains("HypBJet"))      { particle_ = "BJets";}
    else if (variable.Contains("HypBBBar"))     { particle_ = "BBbar";}
    else if (variable.Contains("HypTop"))       { particle_ = "TopQuarks";}
    else if (variable.Contains("HypTTBar"))     { particle_ = "TtBar";}
    else if (variable.Contains("HypJet"))       { particle_ = "Jets";}
    
    if(channel == "combined")
    {
        // Statistically combine the individual channel __NORMALIZED__ covariance matrices
        TH2D ee_StatCov = statCovValues(variable, "ee");
        TH2D emu_StatCov = statCovValues(variable, "emu");
        TH2D mumu_StatCov = statCovValues(variable, "mumu");
        
        ifstream eeResult ("UnfoldingResults/Nominal/ee/"+variable+"Results.txt");
        ifstream emuResult ("UnfoldingResults/Nominal/emu/"+variable+"Results.txt");
        ifstream mumuResult ("UnfoldingResults/Nominal/mumu/"+variable+"Results.txt");
        
        const int nbins = ee_StatCov.GetNbinsX();
        TString Dummy_str="";
        double dummy=0;
        double ee_statError[nbins], emu_statError[nbins], mumu_statError[nbins];
        
        for(int iter=0;iter<nbins; iter++){
            eeResult>>Dummy_str>>dummy>>Dummy_str>>dummy>>Dummy_str>>dummy>>Dummy_str>>dummy>>Dummy_str>>ee_statError[iter]>>Dummy_str>>dummy;
            emuResult>>Dummy_str>>dummy>>Dummy_str>>dummy>>Dummy_str>>dummy>>Dummy_str>>dummy>>Dummy_str>>emu_statError[iter]>>Dummy_str>>dummy;
            mumuResult>>Dummy_str>>dummy>>Dummy_str>>dummy>>Dummy_str>>dummy>>Dummy_str>>dummy>>Dummy_str>>mumu_statError[iter]>>Dummy_str>>dummy;
        }
        return combineStatCovMatrices(ee_StatCov, ee_statError, emu_StatCov, emu_statError, mumu_StatCov, mumu_statError);
    }
    
    TString file_ = "SVD/Unfolding_"+channel+"_"+particle_+"_"+variable_+"_"+variable+".root";
    ifstream file(file_);
    if(!file.is_open()){
        std::cout<<"File '"<<file_<<"'\n cannot be opened"<<std::endl;
        std::cout<<"Exiting!!"<<std::endl;
        exit(12);
    }
    file.close();
    
    TFile f(file_);
    TH2D *histo = (TH2D*)f.Get("SVD_"+channel+"_"+particle_+"_"+variable_+"_"+variable+"_STATCOVNORM");
    const int nbins = histo->GetNbinsX() - 2;
    Double_t bins[nbins-1];
    for (int iter=0; iter<=nbins; iter++) {bins[iter] = histo->GetXaxis()->GetBinLowEdge(iter+2);}
    TH2D newStatCov = TH2D("newHisto", "newHisto", nbins, bins, nbins, bins);
    for (int iterX=2; iterX<=nbins+1; iterX++) {
        for (int iterY=2; iterY<=nbins+1; iterY++) {
            newStatCov.SetBinContent(iterX-1, iterY-1, histo->GetBinContent(iterX, iterY));
        }
    }
    return newStatCov;
}

void fillMatrix(TString channel, TString variable, std::vector<TString> files, const TString outDir, const TH2D statCovariance)
{
    std::vector<TH1*> vec_histos;

    // Get results from ROOT file and store them as histograms
    for (auto file:files){
        file.Append(variable+"_source.root");
        TGraphAsymmErrors* tmp_graph = getGraphFromFile(file, "data_staterror_only");
        TH1*tmp_Histo = getHistoFromFile(file, "mc");
        vec_histos.push_back(convertGraphToHisto(tmp_Histo, tmp_graph));
    }

    std::cout<<"Creating full covariance matrix for  "<<std::setw(22)<<variable;
    std::cout<<" in channel "<<std::setw(6)<<channel<<" saving results in '"<<outDir<<"'"<<std::endl;

    
    // Calcualte symmetrized difference
    // skip 0th element of vector because is the central result which is treated differently
    for(std::size_t iter = 1; iter <=(vec_histos.size()-1)/2; iter++){
        int upHistoNr = 2*iter -1;
        int doHistoNr = 2*iter;
        vec_histos.at(upHistoNr) = getSymmError(vec_histos.at(0), vec_histos.at(upHistoNr), vec_histos.at(doHistoNr), 1);
        vec_histos.at(doHistoNr) = getSymmError(vec_histos.at(0), vec_histos.at(upHistoNr), vec_histos.at(doHistoNr), -1);
        // Normalize the varied histograms histograms
        vec_histos.at(upHistoNr)->Scale(1./vec_histos.at(upHistoNr)->Integral("width"));
        vec_histos.at(doHistoNr)->Scale(1./vec_histos.at(doHistoNr)->Integral("width"));
        // obtain the symmetrized error from the renormalized variations
        vec_histos.at(upHistoNr) = getSymmError(vec_histos.at(0), vec_histos.at(upHistoNr), vec_histos.at(doHistoNr), 1);
        vec_histos.at(doHistoNr) = getSymmError(vec_histos.at(0), vec_histos.at(upHistoNr), vec_histos.at(doHistoNr), -1);
    }

    // Prepare 2D correlation matrix
    const Int_t nbins = vec_histos.at(0)->GetNbinsX();
    double binRanges[nbins+1];
    for (Int_t iter=0; iter<=nbins; iter++){binRanges[iter] = vec_histos.at(0)->GetXaxis()->GetBinLowEdge(iter+1);};
    TH2D* correlationMatrix = new TH2D("matrix", "Correlation Matrix", nbins, binRanges, nbins, binRanges);
    correlationMatrix->SetDirectory(0);
    correlationMatrix->SetTitle("Total (Stat.+Syst.) Correlation Matrix [%]");
    correlationMatrix->GetXaxis()->SetTitle(vec_histos.at(0)->GetXaxis()->GetTitle());
    correlationMatrix->GetYaxis()->SetTitle(vec_histos.at(0)->GetXaxis()->GetTitle());
    
    if(nbins!=statCovariance.GetNbinsX() || nbins!=statCovariance.GetNbinsY())
    {
        std::cout<<"Stat. Covariance Matrix has different number of bins"<<std::endl;
        std::cout<<"EXIT"<<std::endl;
        exit(13);
    }
    
    // Fill 2D correlation matrix
    for (Int_t iterx = 1; iterx<=1+correlationMatrix->GetNbinsX(); iterx++){
        for (Int_t itery = 1; itery<=1+correlationMatrix->GetNbinsY(); itery++){
            double value = 0;
            for(size_t iterSyst = 1; iterSyst< vec_histos.size(); iterSyst +=2){
                /// get only the difference 'Up-Nom'
                double dX = (vec_histos.at(iterSyst)->GetBinContent(iterx) - vec_histos.at(0)->GetBinContent(iterx));
                double dY = (vec_histos.at(iterSyst)->GetBinContent(itery) - vec_histos.at(0)->GetBinContent(itery));
                value += dX * dY;
            }
            // Add the statistical covariance value
            value += statCovariance.GetBinContent(iterx, itery);
            // the correlation matrix is filled with the absolute errors
            correlationMatrix->SetBinContent(iterx, itery, value);
        }
    }
    TH2D *finalMatrix = (TH2D*)correlationMatrix->Clone("finalMatrix");
    TH2D *fakeMatrix = (TH2D*)correlationMatrix->Clone("fakeMatrix");
    // Normalize the covariance matrix
    for (Int_t iterx = 1; iterx<=1+correlationMatrix->GetNbinsX(); iterx++){
        for (Int_t itery = 1; itery<=1+correlationMatrix->GetNbinsY(); itery++){
            double value = correlationMatrix->GetBinContent(iterx, itery);
            double norm = std::sqrt(correlationMatrix->GetBinContent(iterx, iterx) * correlationMatrix->GetBinContent(itery, itery));
            if (norm==0) norm = 1;
            fakeMatrix->SetBinContent(iterx, itery, 100);
            finalMatrix->SetBinContent(iterx, itery, 100*value/norm);
        }
    }
    vec_histos.clear();

    TCanvas *c = new TCanvas();
    gStyle->SetPaintTextFormat("3.1f");
    gStyle->SetOptStat(0);
    fakeMatrix->Draw("box");
    finalMatrix->Draw("text,same");
    c->Print(outDir+variable+".eps");
    c->Print(outDir+variable+".C");

    TFile out_root(outDir+variable+"_source.root", "RECREATE");
    finalMatrix->Write("correlationMatrix_"+variable);
    c->Write("correlationMatrix_"+variable+ "_canvas");
    out_root.Close();

    c->Clear();

    delete c;
    delete correlationMatrix;
    delete finalMatrix;
}



int main()
{
    
    // Details about how to estimate the full correlatio matrices are provided in:
    // https://indico.desy.de/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=8972   (passw: bottom)
    // AN-13-266  and   AN-13-267
    
    gErrorIgnoreLevel = 1001;
    const std::vector<TString> channels = {"mumu", "emu", "ee", "combined"};
    const std::vector<TString> systematics = {"Nominal", "DY_UP", "DY_DOWN", "BG_UP", "BG_DOWN", "KIN_UP", "KIN_DOWN",
                                              "PU_UP", "PU_DOWN", "TRIG_UP", "TRIG_DOWN", "LEPT_UP", "LEPT_DOWN",
                                              "JER_UP", "JER_DOWN", "JES_UP", "JES_DOWN", 
                                              "BTAG_UP", "BTAG_DOWN", "BTAG_LJET_UP", "BTAG_LJET_DOWN",
                                              "BTAG_PT_UP", "BTAG_PT_DOWN", "BTAG_ETA_UP", "BTAG_ETA_DOWN",
                                              "BTAG_LJET_PT_UP", "BTAG_LJET_PT_DOWN", "BTAG_LJET_ETA_UP", "BTAG_LJET_ETA_DOWN",
                                              "MASS_UP", "MASS_DOWN", "MATCH_UP", "MATCH_DOWN", "SCALE_UP", "SCALE_DOWN","HAD_UP", "HAD_DOWN"};

    const std::vector<TString> variables = {"HypLeptonpT", //"HypLeptonpTLead","HypLeptonpTNLead",
                                            "HypLeptonEta",//"HypLeptonEtaLead","HypLeptonEtaNLead",
                                            "HypLLBarpT","HypLLBarMass",//"HypLLBarDPhi",
                                            "HypBJetpT",//"HypBJetpTLead","HypBJetpTNLead",
                                            "HypBJetEta",//"HypBJetEtaLead","HypBJetEtaNLead",
                                            "HypTopRapidity","HypTopRapidityLead","HypTopRapidityNLead",
                                            "HypToppT","HypToppTLead","HypToppTNLead",
                                            "HypTTBarRapidity","HypTTBarpT","HypTTBarMass","HypTTBarDeltaPhi","HypTTBarDeltaRapidity",
                                            "HypToppTTTRestFrame",
                                            "HypLeptonBjetMass",
                                            "HypBBBarpT","HypBBBarMass"};
    
    for(auto chan: channels){
        std::vector<TString> files = createVecFiles(chan, systematics);
        TString outdir = baseOutDir.Copy().Append(chan+"/");
        gSystem->mkdir(outdir, kTRUE);
        for(auto variable: variables){
            fillMatrix(chan, variable, files, outdir, statCovValues(variable, chan));
        };
    };

    return 0;
}
