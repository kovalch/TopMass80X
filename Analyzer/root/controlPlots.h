#ifndef TOPMASSCONTROLPLOTS_H
#define TOPMASSCONTROLPLOTS_H

#include "TChain.h"
#include "TTreeFormula.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "Helper.h"

#include <vector>
#include <boost/algorithm/string/replace.hpp>

enum sampleType {kData, kSig, kBkg, kSigVar};

class TopMassControlPlots{
public:
  TopMassControlPlots();
  ~TopMassControlPlots() {}

private:

  class MySample{
  public:
    MySample(std::string name_, std::string file_, int type_, int color_, int line_ = 1, double scale_ = 1.) :
      name(name_), file(file_), type(type_), color(color_), line(line_), scale(scale_) {}
    
    std::string name, file;
    int type, color, line;
    double scale;
  };

  class MyHistogram{
  public:
    MyHistogram(std::string name_, std::string formulax_, std::string selection_, std::string title, int nBins, double min, double max) :
      varx(0), vary(0),sel(0),
      name(name_), formulax(formulax_), formulay("1."), selection(selection_),
      data(new TH1F((std::string("hD")+name).c_str(), (std::string("Data")+title).c_str(), nBins, min, max)),
      histogramDimension(1)
    {
      data->SetLineWidth(1);
      data->SetLineColor(kBlack);
      data->SetMarkerStyle(20);
      data->SetMarkerColor(kBlack);
    }

    MyHistogram(std::string name_, std::string formulax_, std::string formulay_, std::string selection_, std::string title, int x_nBins, double x_min, double x_max, int y_nBins, double y_min, double y_max/*, double y_plotmin =0., double y_plotmax=-1.*/) :
      varx(0), vary(0),sel(0),
      name(name_), formulax(formulax_), formulay(formulay_), selection(selection_),
      data(new TH2F((std::string("h2D")+name).c_str(), (std::string("Data")+title).c_str(), x_nBins, x_min, x_max, y_nBins, y_min, y_max)),
      histogramDimension(2)
    {
      data->SetLineWidth(1);
      data->SetLineColor(kBlack);
      data->SetMarkerStyle(20);
      data->SetMarkerColor(kBlack);
    }

    void Init(TChain* chain, std::string topBranchName)
    {
      boost::replace_all(formulax,  "top.", topBranchName);
      boost::replace_all(formulay,  "top.", topBranchName);
      boost::replace_all(selection, "top.", topBranchName);
      varx = new TTreeFormula((std::string("fx")+name).c_str(), formulax.c_str() , chain);
      vary = new TTreeFormula((std::string("fy")+name).c_str(), formulay.c_str() , chain);
      if(varx->GetNdim() == 0 || vary->GetNdim() == 0){
        SetInvalid();
      }
      if (selection.size() > 0) sel = new TTreeFormula((std::string("s")+name).c_str(), selection.c_str(), chain);
    }
    void AddSignal(MySample* sample)
    {
      int colorShift[] =  {-11, -8, 0};
      std::string combinationType[] = {" unmatched", " wrong", " correct"};
      for(int i = 0; i < 3; ++i) {
        sig.push_back((TH1*)data->Clone());
        sig.back()->Reset();
        sig.back()->SetLineWidth(1);
        sig.back()->SetFillColor(sample->color+colorShift[i]);
        sig.back()->SetName(HelperFunctions::cleanedName(data->GetName()+sample->name + combinationType[i]).c_str());
        sig.back()->SetTitle((sample->name + combinationType[i]).c_str());
      }
    }
    void AddBackground(MySample* sample)
    {
      bkg.push_back((TH1*)data->Clone());
      bkg.back()->Reset();
      bkg.back()->SetLineWidth(1);
      bkg.back()->SetFillColor(sample->color);
      bkg.back()->SetName(HelperFunctions::cleanedName(data->GetName()+sample->name).c_str());
      bkg.back()->SetTitle(sample->name.c_str());
    }
    void AddSignalVariation(MySample* sample)
    {
      sigvar.push_back((TH1*)data->Clone());
      sigvar.back()->Reset();
      sigvar.back()->SetLineWidth(2);
      sigvar.back()->SetLineColor(sample->color);
      sigvar.back()->SetLineStyle(sample->line);
      sigvar.back()->SetName(HelperFunctions::cleanedName(data->GetName()+sample->name).c_str());
      sigvar.back()->SetTitle(sample->name.c_str());
    }
    
    std::vector<TH1F*> Sigvar1D(){
      std::vector<TH1F*> sigvar1D;
      if(histogramDimension == 1){
        for (size_t i=0; i<sigvar.size(); i++)sigvar1D.push_back((TH1F*) sigvar.at(i));
      }
      return sigvar1D;
    }
    std::vector<TH1F*> Sig1D(){
      std::vector<TH1F*> sig1D;
      if(histogramDimension == 1){
        for (size_t i=0; i<sig.size(); i++)sig1D.push_back((TH1F*) sig.at(i));
      }
      return sig1D;
    }
    std::vector<TH1F*> Bkg1D(){
      std::vector<TH1F*> bkg1D;
      if(histogramDimension == 1){
        for (size_t i=0; i<bkg.size(); i++)bkg1D.push_back((TH1F*) bkg.at(i));
      }
      return bkg1D;
    }
    TH1F* Data1D(){
      return (TH1F*) data;
    }

    std::vector<TH2F*> Sigvar2D(){
      std::vector<TH2F*> sigvar2D;
      if(histogramDimension == 2){
        for (size_t i=0; i<sigvar.size(); i++)sigvar2D.push_back((TH2F*) sigvar.at(i));
      }
      return sigvar2D;
    }
    std::vector<TH2F*> Sig2D(){
      std::vector<TH2F*> sig2D;
      if(histogramDimension == 2){
        for (size_t i=0; i<sig.size(); i++)sig2D.push_back((TH2F*) sig.at(i));
      }
      return sig2D;
    }
    std::vector<TH2F*> Bkg2D(){
      std::vector<TH2F*> bkg2D;
      if(histogramDimension == 2){
        for (size_t i=0; i<bkg.size(); i++)bkg2D.push_back((TH2F*) bkg.at(i));
      }
      return bkg2D;
    }
    TH2F* Data2D(){
      return (TH2F*) data;
    }

    const short Dimension(){
      return histogramDimension;
    }

    void SetInvalid(){
      histogramDimension = -1;
    }
    
    TTreeFormula* varx;
    TTreeFormula* vary;
    TTreeFormula* sel;
    std::string name, formulax, formulay, selection;

  private:
    TH1 *data;
    std::vector<TH1*> sig, bkg, sigvar;
    short histogramDimension;
  };

  class MyBRegVarInfo{
  public:
    MyBRegVarInfo(){
      init();
    }
    void addVariable(std::string varName, std::string varForm, float xMin, float xMax){
      varNames.push_back(varName);
      varForms.push_back(varForm);
      xMins.push_back(xMin);
      xMaxs.push_back(xMax);
    }
    void init(){
      varNames.clear();
      varForms.clear();
      xMins.clear();
      xMaxs.clear();


      addVariable("corr. p_{T}"                             , "BRegJet.jetPtCorr"                 , 0   , 500   );                                 
      addVariable("jet #eta"                                , "BRegJet.jetEta"                    , -3  , 3     );                              
      addVariable("#rho_{25}"                               , "BRegJet.Rho25"                     , 0   , 50    );                             
      addVariable("jet area"                                , "BRegJet.jetArea"                   , 0.5 , 1     );                               
      addVariable("width in #varphi"                        , "BRegJet.EtWeightedSigmaPhi"        , 0   , 0.3   );  //5 width and PU                                          
      addVariable("leading chargedConstPt"                  , "BRegJet.leadingChargedConstPt"     , 0   , 200   );                                             
      addVariable("CHF"                                     , "jet.fChargedHadron"                , 0   , 1     );                                  
      addVariable("electron fraction"                       , "jet.fElectron"                     , 0   , 1     );                             
      addVariable("muon fraction"                           , "jet.fMuon"                         , 0   , 1     );  //9 charged constituents 
      addVariable("electron fraction_2"                     , "jet.fElectron"                     , 0.01, 1     );                             
      addVariable("muon fraction_2"                         , "jet.fMuon"                         , 0.01, 1     );  //9 charged constituents 
      addVariable("SV flight length"                        , "BRegJet.SV3DLength"                , 0   , 5     );                                  
      addVariable("SV flight length unc."                   , "BRegJet.SV3DLengthError"           , 0   , 1     );                                       
      addVariable("SV mass"                                 , "BRegJet.SVMass"                    , 0   , 20    );                              
      addVariable("SV p_{T}"                                , "BRegJet.SVPt"                      , 0   , 200   );                            
      addVariable("CSV b-tag"                               , "jet.bTagCSV"                       , 0   , 1     );  //14 SVtx and btag                         
      addVariable("p_{T}^{SoftMuon}"                        , "BRegJet.SoftMuonPt"                , 0   , 200   );                                  
      addVariable("relative p_{T}^{SoftMuon} fraction"      , "BRegJet.SoftMuonRatioRel"          , 0   , 2     );                                        
      addVariable("#Delta R(SoftMuon)"                      , "BRegJet.SoftMuonDeltaR"            , 0   , 0.5   );                                      
      addVariable("p_{T}^{SoftElectron}"                    , "BRegJet.SoftElectronPt"            , 0   , 200   );                                      
      addVariable("relative p_{T}^{SoftElectron} fraction"  , "BRegJet.SoftElectronRatioRel"      , 0   , 2     );                                            
      addVariable("#Delta R(SoftElectron)"                  , "BRegJet.SoftElectronDeltaR"        , 0   , 0.5   );  //20 soft lepton variables 
      addVariable("raw p_{T}"                               , "BRegJet.jetPtRaw"                  , 0   , 500   );                                
      addVariable("m_{T}"                                   , "BRegJet.jetMt"                     , 0   , 500   );                             
      addVariable("JES uncertainty"                         , "BRegJet.jesTotUnc"                 , 0   , 0.05  );         
      addVariable("number of constituents"                  , "jet.nConstituents"                 , 0   , 100   );         
      addVariable("number of charged hadrons"               , "jet.nChargedHadrons"               , 0   , 50    );           
      addVariable("number of charged PF candidates"         , "BRegJet.nChargedPFConstituents"    , 0   , 50    );  //26 more jet info, multiplicities
      addVariable("RlbReco"                                 , "BRegJet.RlbReco"                   , 0   , 2     );  //27 add ATLAS 3D variable into training
      //add QGL and PU Id variables
      addVariable("BRegJet.QGaxis1"      ,"BRegJet.QGaxis1"      , 0    ,   0.3);
      addVariable("BRegJet.QGaxis2"      ,"BRegJet.QGaxis2"      , 0    ,   0.3);
      addVariable("BRegJet.QGMult"       ,"BRegJet.QGMult"       , 0    ,   30 );
      addVariable("BRegJet.QGPtD"        ,"BRegJet.QGPtD"        , 0    ,   1  );
      addVariable("BRegJet.QGMLP"        ,"BRegJet.QGMLP"        , 0    ,   1  );

      addVariable("BRegJet.PUIddZ"       ,"BRegJet.PUIddZ"       , 0    ,   2  );
      addVariable("BRegJet.PUIddRMean"   ,"BRegJet.PUIddRMean"   , 0    ,   0.3  );
      addVariable("BRegJet.PUIddr2Mean"  ,"BRegJet.PUIddr2Mean"  , 0    ,   0.1  );
      addVariable("BRegJet.PUIdfrac01"   ,"BRegJet.PUIdfrac01"   , 0    ,   1  );
      addVariable("BRegJet.PUIdfrac02"   ,"BRegJet.PUIdfrac02"   , 0    ,   1  );
      addVariable("BRegJet.PUIdfrac03"   ,"BRegJet.PUIdfrac03"   , 0    ,   1  );
      addVariable("BRegJet.PUIdfrac04"   ,"BRegJet.PUIdfrac04"   , 0    ,   1  );
      addVariable("BRegJet.PUIdfrac05"   ,"BRegJet.PUIdfrac05"   , 0    ,   1  );
      addVariable("BRegJet.PUIdbeta"     ,"BRegJet.PUIdbeta"     , 0    ,   1  );
      addVariable("BRegJet.PUIdbetaStar" ,"BRegJet.PUIdbetaStar" , 0    ,   1  );
      addVariable("BRegJet.PUIdptD"      ,"BRegJet.PUIdptD"      , 0    ,   1  );

      addVariable("GBR factor"                              , "BRegJet.BRegGBRTrainResult"        , 0.5 , 1.5   );
      addVariable("TMVA factor"                             , "BRegJet.BRegResult"                , 0.5 , 1.5   );
    }
    
    std::vector<std::string> varNames;
    std::vector<std::string> varForms;
    std::vector<float> xMins;
    std::vector<float> xMaxs;
  };


  void doPlots();

  std::vector<MyHistogram> hists;
  std::vector<MySample>    samples;
};

#endif /* TOPMASSCONTROLPLOTS_H */
