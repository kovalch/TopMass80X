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

TH1* calculateDifference(TH1* h0, TH1* h1, bool absValue)
{
    /// Calculate the binwise difference between 2 histograms: h_final (i) = h0(i)- h1(i)
    /// Please note that if 'absValue == true' the absolute value of the difference is returned
    TH1* htmp = dynamic_cast<TH1*> (h0->Clone("htmp"));
    for(auto iter=0; iter<= h1->GetNbinsX()+1; iter++)
    {
        double value = h0->GetBinContent(iter) - h1->GetBinContent(iter);
        if(absValue){htmp->SetBinContent(iter, std::fabs(value));}
        else {htmp->SetBinContent(iter, value);};
    }
    return htmp;
}

TH1* getSymmError(TH1* hCent, TH1* hUp, TH1* hDo, int sign)
{
    /// This function returns the value of: hCent + factor *symmError
    /// where symmError = 0.5*(abs(hUp - hCent) + abs(hDo - hCent)) * sign(hUp - hCent)

    if(std::abs(sign) != 1) {std::cout<<"getSymmError: factor not '1' or '-1'. Not valid"<<std::endl; exit(89);};

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


void fillMatrix(TString variable, std::vector<TString> files, const TString outDir)
{
    std::vector<TH1*> vec_histos;

    /// Get results from ROOT file and store them as histograms
    for (auto file:files){
        file.Append(variable+"_source.root");
        TGraphAsymmErrors* tmp_graph = getGraphFromFile(file, "data_staterror_only");
        TH1*tmp_Histo = getHistoFromFile(file, "mc");
        vec_histos.push_back(convertGraphToHisto(tmp_Histo, tmp_graph));
    }

    /// Calcualte symmetrized difference
    /// skip 0th element of vector because is the cetral result
    for(std::size_t iter = 1; iter <=(vec_histos.size()-1)/2; iter++){
        int upHistoNr = 2*iter -1;
        int doHistoNr = 2*iter;
        vec_histos.at(upHistoNr) = getSymmError(vec_histos.at(0), vec_histos.at(upHistoNr), vec_histos.at(doHistoNr), 1);
        vec_histos.at(doHistoNr) = getSymmError(vec_histos.at(0), vec_histos.at(upHistoNr), vec_histos.at(doHistoNr), -1);
        /// Normalize the varied histograms histograms
        vec_histos.at(upHistoNr)->Scale(1./vec_histos.at(upHistoNr)->Integral("width"));
        vec_histos.at(doHistoNr)->Scale(1./vec_histos.at(doHistoNr)->Integral("width"));
        /// obtain the symmetrized error from the renormalized variations
        vec_histos.at(upHistoNr) = getSymmError(vec_histos.at(0), vec_histos.at(upHistoNr), vec_histos.at(doHistoNr), 1);
        vec_histos.at(doHistoNr) = getSymmError(vec_histos.at(0), vec_histos.at(upHistoNr), vec_histos.at(doHistoNr), -1);
    }

    /// Prepare 2D correlation matrix
    const Int_t nbins = vec_histos.at(0)->GetNbinsX();
    double binRanges[nbins+1];
    for (Int_t iter=0; iter<=nbins; iter++){binRanges[iter] = vec_histos.at(0)->GetXaxis()->GetBinLowEdge(iter+1);};
    TH2D* correlationMatrix = new TH2D("matrix", "Correlation Matrix", nbins, binRanges, nbins, binRanges);
    correlationMatrix->SetDirectory(0);
    correlationMatrix->SetTitle("Systematic Correlation Matrix");
    correlationMatrix->GetXaxis()->SetTitle(vec_histos.at(0)->GetXaxis()->GetTitle());
    correlationMatrix->GetYaxis()->SetTitle(vec_histos.at(0)->GetXaxis()->GetTitle());

    /// Fill 2D correlation matrix
    for (Int_t iterx = 1; iterx<=1+correlationMatrix->GetNbinsX(); iterx++){
        for (Int_t itery = 1; itery<=1+correlationMatrix->GetNbinsY(); itery++){
            double value = 0;
            for(size_t iterSyst = 1; iterSyst< vec_histos.size(); iterSyst +=2){
                /// get only the difference 'Up-Nom'
                double dX = (vec_histos.at(iterSyst)->GetBinContent(iterx) - vec_histos.at(0)->GetBinContent(iterx));
                double dY = (vec_histos.at(iterSyst)->GetBinContent(itery) - vec_histos.at(0)->GetBinContent(itery));
                value += dX * dY;
            }
            /// the correlation matrix is filled with the absolute errors
            correlationMatrix->SetBinContent(iterx, itery, value);
        }
    }
    vec_histos.clear();

    TCanvas *c = new TCanvas();
    gStyle->SetOptStat(0);
    correlationMatrix->Draw("text");
    c->Print(outDir+variable+".eps");
    c->Print(outDir+variable+".C");

    TFile out_root(outDir+variable+"_source.root", "RECREATE");
    correlationMatrix->Write("correlationMatrix_"+variable);
    c->Write("correlationMatrix_"+variable+ "_canvas");
    out_root.Close();

    c->Clear();

    delete c;
    delete correlationMatrix;
}



int main()
{
    gErrorIgnoreLevel = 1001;
    const std::vector<TString> channels = {"mumu", "emu", "ee"};
    const std::vector<TString> systematics = {"Nominal", "DY_UP", "DY_DOWN", "BG_UP", "BG_DOWN", "KIN_UP", "KIN_DOWN",
                                              "PU_UP", "PU_DOWN", "TRIG_UP", "TRIG_DOWN", "LEPT_UP", "LEPT_DOWN",
                                              "JER_UP", "JER_DOWN", "JES_UP", "JES_DOWN", 
                                              "BTAG_UP", "BTAG_DOWN", "BTAG_LJET_UP", "BTAG_LJET_DOWN",
                                              "BTAG_PT_UP", "BTAG_PT_DOWN", "BTAG_ETA_UP", "BTAG_ETA_DOWN",
                                              "BTAG_LJET_PT_UP", "BTAG_LJET_PT_DOWN", "BTAG_LJET_ETA_UP", "BTAG_LJET_ETA_DOWN",
                                              "MASS_UP", "MASS_DOWN", "MATCH_UP", "MATCH_DOWN", "SCALE_UP", "SCALE_DOWN","HAD_UP", "HAD_DOWN"};

    const std::vector<TString> variables = {"HypLeptonpT", "HypLeptonpTLead","HypLeptonpTNLead",
                                            "HypLeptonEta","HypLeptonEtaLead","HypLeptonEtaNLead",
                                            "HypLLBarpT","HypLLBarMass",//"HypLLBarDPhi",
                                            "HypBJetpT","HypBJetpTLead","HypBJetpTNLead",
                                            "HypBJetEta","HypBJetEtaLead","HypBJetEtaNLead",
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
            std::cout<<"Creating correlation matrix for  "<<std::setw(22)<<variable<<"   in channel "<<std::setw(6)<<chan<<"   saving results in '"<<outdir<<"'"<<std::endl;
            fillMatrix(variable, files, outdir);
        };
    };

    return 0;
}
