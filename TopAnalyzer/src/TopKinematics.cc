#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLepEvtPartons.h"
#include "TopAnalysis/TopAnalyzer/interface/TopKinematics.h"
#include "TopAnalysis/TopAnalyzer/interface/TopAngles.h"
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"


/// default constructor for generator level analysis in fw lite
TopKinematics::TopKinematics() : hypoKey_(""), lepton_("muon"), matchForStabilityAndPurity_(false) 
{
  tree = 0;
  useTree_=false;
  ttbarInsteadOfLepHadTop_ = false;

  /// check if input is correct
  if (lepton_.compare("muon")!=0 && lepton_.compare("electron")!=0) throw edm::Exception( edm::errors::Configuration, "lepton specified incorrectly; has to be either 'muon' or 'electron'" ) ;
}

/// default constructor for reco level analysis in fw lite
TopKinematics::TopKinematics(const std::string& hypoKey, const std::string& lepton, const bool& matchForStabilityAndPurity) : hypoKey_(hypoKey), lepton_(lepton), matchForStabilityAndPurity_(matchForStabilityAndPurity)
{
  tree = 0;
  useTree_=false;
  ttbarInsteadOfLepHadTop_ = false;

  /// check if input is correct
  if (lepton_.compare("muon")!=0 && lepton_.compare("electron")!=0) throw edm::Exception( edm::errors::Configuration, "lepton specified incorrectly; has to be either 'muon' or 'electron'" ) ;
}

/// default constructor for full fw
TopKinematics::TopKinematics(const edm::ParameterSet& cfg) :
  hypoKey_( cfg.getParameter<std::string>("hypoKey") ),
  lepton_ ( cfg.getParameter<std::string>("lepton")  ),
  useTree_( cfg.getParameter<bool>       ("useTree") ),
  matchForStabilityAndPurity_( cfg.getParameter<bool>("matchForStabilityAndPurity") ),
  ttbarInsteadOfLepHadTop_   ( cfg.getParameter<bool>("ttbarInsteadOfLepHadTop"   ) ),
  maxNJets( cfg.getParameter<int> ("maxNJets") ),
  ndof( cfg.getParameter<int> ("ndof") )
{
  tree = 0;
  /// check if input is correct
  if (lepton_.compare("muon")!=0 && lepton_.compare("electron")!=0) throw edm::Exception( edm::errors::Configuration, "lepton specified incorrectly; has to be either 'muon' or 'electron'" ) ;

}

/// histogramm booking for fwlite 
void TopKinematics::book()
{

  /** 
      ttbar decay setup
  **/
  // decay channel
  hists_["decayChannel"] = new TH1F( "decayChannel", "decayChannel",  25,   -4.,   40.);

  /** 
      KinFit Monitoring Variables
  **/
  // fit probability of the best fit hypothesis
  hists_["prob"       ] = new TH1F( "prob"       , "prob"       ,   1000,   0.,     1.);
  // chi2 of the best fit hypothesis						      
  hists_["chi2"       ] = new TH1F( "chi2"       , "chi2"       ,   1000,   0.,   100.);
  // delta chi2 between best and second best fit hyothesis 
  hists_["delChi2"    ] = new TH1F( "delChi2"    , "delChi2"    ,   1000,   0.,   100.);

  /** 
      Top Variables for Cross Section Measurement
  **/
  // top pt (at the moment both top candidates are filled in one histogram)
  hists_["topPt"      ] = new TH1F( "topPt"      , "topPt"      ,  800,  0. ,  800.);
  // top y (at the moment both top candidates are filled in one histogram)
  hists_["topY"       ] = new TH1F( "topY"       , "topY"       ,  100, -5. ,    5.);
  // top phi (at the moment both top candidates are filled in one histogram)
  hists_["topPhi"     ] = new TH1F( "topPhi"     , "topPhi"     ,  800, -4  ,    4.);
  // ttbar pair pt
  hists_["ttbarPt"    ] = new TH1F( "ttbarPt"    , "ttbarPt"    ,  300,  0. ,  300.);
  // ttbar pair 
  hists_["ttbarY"     ] = new TH1F( "ttbarY"     , "ttbarY"     ,  100, -5. ,    5.);
  // ttbar pair phi
  hists_["ttbarPhi"   ] = new TH1F( "ttbarPhi"   , "ttbarPhi"   ,  800, -4  ,    4.);
  // ttbar pair invariant mass
  hists_["ttbarMass"  ] = new TH1F( "ttbarMass"  , "ttbarMass"  , 2500,  0. , 2500 );
  // deltaPhi between both top quarks
  hists_["ttbarDelPhi"] = new TH1F( "ttbarDelPhi", "ttbarDelPhi",  400,  0. ,    4.);
  // deltaY between both top quarks
  hists_["ttbarDelY"  ] = new TH1F( "ttbarDelY"  , "ttbarDelY"  , 1000, -5. ,    5.);
  // sum of y of both top quarks
  hists_["ttbarSumY"  ] = new TH1F( "ttbarSumY"  , "ttbarSumY"  ,  100, -5. ,    5.);
  // HT of the 4 jets assigned to the ttbar decay
  hists_["ttbarHT"    ] = new TH1F( "ttbarHT"    , "ttbarHT"    , 1500,  0. , 1500.);
  // ttbar phiStar
  // angular based quantity with same sensitivity as pT(ttbar)
  // (inspired by DY->ll: http://authors.library.caltech.edu/38127/1/1-s2.0-S0370269313000956-main.pdf)
  hists_["ttbarPhiStar"] = new TH1F( "ttbarPhiStar" , "ttbarPhiStar" , 200, 0., 2.);  

  /**
     Top Variables for Cross Checks
  **/
  // pt of the leading top candidate
  hists_["topPtLead"   ] = new TH1F( "topPtLead"     , "topPtLead"    , 800, 0., 800.);  
  // pt of the subleading top candidate
  hists_["topPtSubLead"] = new TH1F( "topPtSubLead"  , "topPtSubLead" , 800, 0., 800.);
  // pt of the top candidates in the ttbar restframe
  hists_["topPtTtbarSys"]= new TH1F( "topPtTtbarSys" , "topPtTtbarSys", 800, 0., 800.);
  // y of the leading top candidate
  hists_["topYLead"   ] = new TH1F( "topYLead"       , "topYLead"   ,   100, -5. ,    5.);
  // y of the subleading top candidate 
  hists_["topYSubLead"] = new TH1F( "topYSubLead"    , "topYSubLead",   100, -5. ,    5.);
  // pt of the hadronically decaying top candidate
  hists_["topPtHad"  ] = new TH1F( "topPtHad"  , "topPtHad"  ,  800,  0. ,  800.);
  // y  of the hadronically decaying top candidate
  hists_["topYHad"   ] = new TH1F( "topYHad"   , "topYHad"   ,  100, -5. ,    5.);
  // y  of the hadronically decaying top candidate
  hists_["topPhiHad" ] = new TH1F( "topPhiHad" , "topPhiHad" ,  800, -4  ,    4.);
  // pt of the leptonically decaying top candidate
  hists_["topPtLep"  ] = new TH1F( "topPtLep"  , "topPtLep"  ,  800,  0. ,  800.);
  // y  of the leptonically decaying top candidate
  hists_["topYLep"   ] = new TH1F( "topYLep"   , "topYLep"   ,  100, -5. ,    5.);
  // y  of the leptonically decaying top candidate
  hists_["topPhiLep" ] = new TH1F( "topPhiLep" , "topPhiLep" ,  800, -4. ,    4.);
  // hadronic Top mass
  hists_["hadTopMass"] = new TH1F( "hadTopMass", "hadTopMass", 1000,  0.,  1000.);
  // leptonic Top mass								
  hists_["lepTopMass"] = new TH1F( "lepTopMass", "lepTopMass", 1000,  0.,  1000.);
  // Top mass									
  hists_["topMass"   ] = new TH1F( "topMass"   , "topMass"   , 1000,  0.,  1000.);
  
  /**
     bquark system
  **/							
  hists_["bbbarPt"  ] = new TH1F( "bbbarPt"   , "bbbarPt"   , 1000,  0.,  1000.);
  hists_["bbbarY"   ] = new TH1F( "bbbarY"    , "bbbarY"    , 1000, -5.,  5.   );
  hists_["bbbarMass"] = new TH1F( "bbbarMass" , "bbbarMass" , 1500,  0.,  1500.);
  // leptonic b + lepton invariant mass
  hists_["lbMass"   ] = new TH1F( "lbMass"    , "lbMass"    , 500,   0.,  500. );
  // b-top hadronization variable
  // (see e.g.: http://link.springer.com/content/pdf/10.1140%2Fepjc%2Fs10052-009-1170-4.pdf)
  //hists_["xB"       ] = new TH1F( "xB"        , "xB"        , 200,  -0.5,  1.5 );

  /**
     asymmetry variables
  **/
  hists_["topPtPlus"  ] = new TH1F( "topPtPlus"  , "topPtPlus"  ,  800,  0. ,  800.);
  hists_["topPtMinus" ] = new TH1F( "topPtMinus" , "topPtMinus" ,  800,  0. ,  800.);
  hists_["lepPtPlus"  ] = new TH1F( "lepPtPlus"  , "lepPtPlus"  , 1200,  0. , 1200.);
  hists_["lepPtMinus" ] = new TH1F( "lepPtMinus" , "lepPtMinus" , 1200,  0. , 1200.);
  hists_["lepEtaPlus"     ] = new TH1F( "lepEtaPlus"   , "lepEtaPlus"   ,  100,  -5. ,  5. );
  hists_["lepEtaMinus"    ] = new TH1F( "lepEtaMinus"  , "lepEtaMinus"  ,  100,  -5. ,  5. );
  hists_["lepYPlus"       ] = new TH1F( "lepYPlus"     , "lepYPlus"     ,  100,  -5. ,  5. );
  hists_["lepYMinus"      ] = new TH1F( "lepYMinus"    , "lepYMinus"    ,  100,  -5. ,  5. );
  hists_["topEtaPlus"     ] = new TH1F( "topEtaPlus"   , "topEtaPlus"   ,  100,  -5. ,  5. );
  hists_["topEtaMinus"    ] = new TH1F( "topEtaMinus"  , "topEtaMinus"  ,  100,  -5. ,  5. );
  hists_["topYPlus"       ] = new TH1F( "topYPlus"     , "topYPlus"     ,  100,  -5. ,  5. );
  hists_["topYMinus"      ] = new TH1F( "topYMinus"    , "topYMinus"    ,  100,  -5. ,  5. );

  /**
     Angular distributions
  **/
  // angle between t candidate (detector rest frame)
  hists_["ttbarAngle"       ] = new TH1F( "ttbarAngle"       , "ttbarAngle"          , 315,  0.  ,  pi );
  // angle between b candidates (detector rest frame)
  hists_["bbbarAngle"       ] = new TH1F( "bbbarAngle"       , "bbbarAngle"          , 315,  0.  ,  pi );
  // angle between b candidates (ttbar rest frame)
  hists_["bbbarAngleTtRF"   ] = new TH1F( "bbbarAngleTtRF"   , "bbbarAngleTtR"       , 315,  0.  ,  pi );
  // angle between W bosons (ttbar rest frame)
  hists_["WWAngle"          ] = new TH1F( "WWAngle"          , "WWAngle"             , 315,  0.  ,  pi );
  // angle between the leptonically decaying top candidate and the corresponding W 
  // (Top in ttbar rest frame and W in top rest frame)
  hists_["topWAngleLep"     ] = new TH1F( "topWAngleLep"     , "topWAngleLep"        , 315,  0.  ,  pi );
  // angle between the hadronically decaying top candidate and the corresponding W 
  // (Top in ttbar rest frame and W in top rest frame)   
  hists_["topWAngleHad"     ] = new TH1F( "topWAngleHad"     , "topWAngleHad"        , 315,  0.  ,  pi );
  // angle between the leptonically decaying top candidate and the corresponding b (ttbar rest frame)
  hists_["topBAngleLep"     ] = new TH1F( "topBAngleLep"     , "topBAngleLep"        , 315,  0.  ,  pi );
  // angle between the hadronically decaying top candidate and the corresponding b (ttbar rest frame)
  hists_["topBAngleHad"     ] = new TH1F( "topBAngleHad"     , "topBAngleHad"        , 315,  0.  ,  pi );
  // angle between the leptonically decaying W candidate and the corresponding b (ttbar rest frame)
  hists_["bWAngleLep"       ] = new TH1F( "bWAngleLep"       , "bWAngleLep"          , 315,  0.  ,  pi );
  // angle between the hadronically decaying W candidate and the corresponding b (ttbar rest frame)
  hists_["bWAngleHad"       ] = new TH1F( "bWAngleHad"       , "bWAngleHad"          , 315,  0.  ,  pi );
  // angle between light quarks (ttbar rest frame)
  hists_["qqbarAngle"       ] = new TH1F( "qqbarAngle"       , "qqbarAngle"          , 315,  0.  ,  pi );
  // angle between light quarks and leptonic b (ttbar rest frame)
  hists_["qBlepAngle"       ] = new TH1F( "qBlepAngle"       , "qBlepAngle"          , 315,  0.  ,  pi );
  // angle between light quarks and hadronic b (corresponding top candidate rest frame)
  hists_["qBhadAngle"       ] = new TH1F( "qBhadAngle"       , "qBhadAngle"          , 315,  0.  ,  pi );
  // angle between lepton and leptonic b (corresponding top candidate rest frame)
  hists_["lepBlepAngle"     ] = new TH1F( "lepBlepAngle"     , "lepBlepAngle"        , 315,  0.  ,  pi );
  // angle between lepton and leptonic b ttbar rest frame)
  hists_["lepBlepAngleTtRF" ] = new TH1F( "lepBlepAngleTtRF" , "lepBlepAngleTtRF"    , 315,  0.  ,  pi );
  // angle between lepton and hadronic b (ttbar rest frame)
  hists_["lepBhadAngle"     ] = new TH1F( "lepBhadAngle"     , "lepBhadAngle"        , 315,  0.  ,  pi );
  // angle between lepton and light quarks (ttbar rest frame)
  hists_["lepQAngle"        ] = new TH1F( "lepQAngle"        , "lepQAngle"           , 315,  0.  ,  pi );
  // angle between muon and neutrino (in corresponding top rest frame)
  hists_["MuonNeutrinoAngle"] = new TH1F( "MuonNeutrinoAngle", "MuonNeutrinoAngle"   , 315,  0.  ,  pi );
  // angle between leptonic b and neutrino (in corresponding top rest frame)
  hists_["lepBNeutrinoAngle"] = new TH1F( "lepBNeutrinoAngle", "lepBNeutrinoAngle"   , 315,  0.  ,  pi );
  // angle between hadronic b and neutrino (in ttbar rest frame)
  hists_["hadBNeutrinoAngle"] = new TH1F( "hadBNeutrinoAngle", "hadBNeutrinoAngle"   , 315,  0.  ,  pi );
  // angle between light quarks and neutrino (in ttbar rest frame)
  hists_["qNeutrinoAngle"   ] = new TH1F( "qNeutrinoAngle"   , "qNeutrinoAngle"      , 315,  0.  ,  pi );
  // angle between lepton (in W rest frame) and corresponding W (in detector rest frame)
  hists_["lepWDir"          ] = new TH1F( "lepWDir"          , "lepWDir"             , 315,  0.  ,  pi );
  // angle between neutrino (in W rest frame) and corresponding W (in detector rest frame)
  hists_["nuWDir"           ] = new TH1F( "nuWDir"           , "nuWDir"              , 315,  0.  ,  pi );
  // angle between light quarks (in W rest frame) and corresponding W (in detector rest frame)
  hists_["qWDir"            ] = new TH1F( "qWDir"            , "qWDir"               , 315,  0.  ,  pi );

  /** 
      event shape variables
  **/
  // notation:
  // 0<=q1<=q2<=q3 are the eigenvalues of the momemtum tensor:
  // T_ab = sum{p_j[a]*p_j[b]}/sum{p_j**2} normalized to 1.

  // aplanarity of event: 
  // 1.5*q1 
  // Return values are 0.5 for spherical and 0 for plane and linear events
  hists_["aplanarity" ] = new TH1F("aplanarity" , "aplanarity" , 50 ,  0. , 0.5);
  // sphericity of event:
  // 1.5*(q1+q2) 
  // Return values are 1 for spherical, 3/4 for plane and 0 for linear events
  hists_["sphericity" ] = new TH1F("sphericity" , "sphericity" , 100,  0. , 1. );
  // C (three jet structure) of event:
  // 3.*(q1*q2+q1*q3+q2*q3)
  // Return value is between 0 and 1, C vanishes for a "perfect" 2-jet event
  hists_["C"          ] = new TH1F("C"          , "C"          , 100,  0. , 1. ); 
  // D (four jet structure) of event:
  // 27.*(q1*q2*q3)    
  // Return value is between 0 and 1, D vanishes for a planar event
  hists_["D"          ] = new TH1F("D"          , "D"          , 100,  0. , 1. );
  // circularity of event
  // the return value is 1 for spherical and 0 linear events in r-phi
  hists_["circularity"] = new TH1F("circularity", "circularity", 100,  0. , 1. ); 
  // isotropy of event:
  // the return value is 1 for spherical events and 0 for events linear in r-phi
  hists_["isotropy"   ] = new TH1F("isotropy"   , "isotropy"   , 100,  0. , 1. ); 


  /** 
      ttbar final state object distributions
  **/
  // leptonPt
  hists_["lepPt"      ] = new TH1F( "lepPt"      , "lepPt"      ,  1200,  0. ,  1200.);
  // leptonEta
  hists_["lepEta"     ] = new TH1F( "lepEta"     , "lepEta"     ,  100,  -5. ,  5.   );
  // neutrinoPt
  hists_["neutrinoPt" ] = new TH1F( "neutrinoPt" , "neutrinoPt" ,  1200,  0. ,  1200.);
  // neutrinoEta
  hists_["neutrinoEta"] = new TH1F( "neutrinoEta", "neutrinoEta",  100,  -5. ,  5.   );
  // light-quarks Pt
  hists_["lightqPt"   ] = new TH1F( "lightqPt"   , "lightqPt"   ,  1200,  0. ,  1200.);
  // light-quarks Eta
  hists_["lightqEta"  ] = new TH1F( "lightqEta"  , "lightqEta"  ,  100,  -5. ,  5.   );
  // b-quarks Pt
  hists_["bqPt"       ] = new TH1F( "bqPt"       , "bqPt"       ,  1200,  0. ,  1200.);
  // b-quarks Eta
  hists_["bqEta"      ] = new TH1F( "bqEta"      , "bqEta"      ,  100,  -5. ,  5.   );
  // leading quark Pt
  hists_["leadqPt"    ] = new TH1F( "leadqPt"    , "leadqPt"    ,  1200,  0. ,  1200.);
  // leading quark Eta
  hists_["leadqEta"   ] = new TH1F( "leadqEta"   , "leadqEta"   ,  100,  -5. ,  5.   );
 
  /** 
      wrong reconstructed quarks
  **/
  hists_["qAssignment"] = new TH1F( "qAssignment", "qAssignment",     10, -0.5,   9.5 );
  // set labels (assignment=0 corresponds to bin 1)
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(1, "ok"      );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(2, "bb"      );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(3, "blepq"   );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(4, "bhadq"   );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(5, "bbqlep"  );    
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(6, "bbqhad"  );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(7, "bbqq"    );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(8, "jmis"    );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(9, "wrongj"  );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(10,"nomatch" );

  /** 
      Correlation Plots
  **/
  if(!hypoKey_.compare("None")==0){
    // gen-rec level correlation top pt
    corrs_["topPt_"      ] = new TH2F( "topPt_"      , "topPt_"       , 800 ,    0.,  800.,     800,   0.,  800.);
    // gen-rec level correlation top y
    corrs_["topY_"       ] = new TH2F( "topY_"       , "topY_"        ,1000 ,   -5.,    5.,    1000,  -5.,    5.);
    // gen-rec level correlation top phi
    //corrs_["topPhi_"     ] = new TH2F( "topPhi_"     , "topPhi_"      , 628 , -3.14,  3.14,     628,-3.14,  3.14);
    // gen-rec level correlation ttbar pt
    corrs_["ttbarPt_"    ] = new TH2F( "ttbarPt_"    , "ttbarPt_"     , 300 ,    0.,  300.,     300,   0.,  300.);
    // gen-rec level correlation ttbar y
    corrs_["ttbarY_"     ] = new TH2F( "ttbarY_"     , "ttbarY_"      , 1000,   -5.,    5.,    1000,  -5.,    5.);
    // gen-rec level correlation ttbar mass
    corrs_["ttbarMass_"  ] = new TH2F( "ttbarMass_"  , "ttbarMass_"   , 2500 ,   0., 2500.,     2500,  0., 2500.);
    // gen-rec level correlation HT of the 4 jets assigned to the ttbar decay
    //corrs_["ttbarHT_"    ] = new TH2F( "ttbarHT_"    , "ttbarHT_"     , 1500,    0., 1500.,    1500,   0., 1500.);
    // gen-rec level correlation ttbar deltaPhi
    corrs_["ttbarDelPhi_"] = new TH2F( "ttbarDelPhi_", "ttbarDelPhi_" , 400 ,    0.,    4.,     400,   0.,    4.);
    // gen-rec level correlation ttbar deltaY
    corrs_["ttbarDelY_"  ] = new TH2F( "ttbarDelY_"  , "ttbarDelY_"   , 1000,   -5.,    5.,    1000,  -5.,    5.);
    // gen-rec level correlation ttbar sumY
    //corrs_["ttbarSumY_"  ] = new TH2F( "ttbarSumY_"  , "ttbarSumY_"   , 1000,   -5.,    5.,    1000,  -5.,    5.);
    // gen-rec level correlation angle between b jets
    //corrs_["bbbarAngle"  ] = new TH2F( "bbbarAngle"  , "bbbarAngle"   ,  315,    0.,  3.15,     315,   0.,  3.15);
    // gen-rec level correlation for lepton charge
    //corrs_["lepCharge_"  ] = new TH2F( "lepCharge_"  , "lepCharge_"   ,    3,  -1.5,   1.5,       3, -1.5,   1.5);
    // gen-rec level correlation for angle between the leptonically decaying top candidate and the neutrino
    //corrs_["MuonNeutrinoAngle_"] = new TH2F( "MuonNeutrinoAngle_", "MuonNeutrinoAngle_", 315,  0.  ,  3.15,   315,  0.  ,  3.15);
    // gen-rec level correlation for leptonPt
    //corrs_["lepPt_"      ] = new TH2F( "lepPt_"      , "lepPt_"      ,  1200,  0. ,  1200. ,  1200,  0. ,  1200.);
    // gen-rec level correlation for leptonEta
    //corrs_["lepEta_"     ] = new TH2F( "lepEta_"     , "lepEta_"     ,  100,  -5. ,  5.    ,  100,  -5. ,  5.   );
    // gen-rec level correlation for neutrinoPt
    //corrs_["neutrinoPt_" ] = new TH2F( "neutrinoPt_" , "neutrinoPt_" ,  1200,  0. ,  1200. ,  1200,  0. ,  1200.);
    // gen-rec level correlation for neutrinoEta
    //corrs_["neutrinoEta_"] = new TH2F( "neutrinoEta_", "neutrinoEta_",  100,  -5. ,  5.    ,  100,  -5. ,  5.   );
    // gen-rec level correlation for light-quarks Pt
    //corrs_["lightqPt_"   ] = new TH2F( "lightqPt_"   , "lightqPt_"   ,  1200,  0. ,  1200. ,  1200,  0. ,  1200.);
    // gen-rec level correlation for light-quarks Eta
    //corrs_["lightqEta_"  ] = new TH2F( "lightqEta_"  , "lightqEta_"  ,  100,  -5. ,  5.    ,  100,  -5. ,  5.   );
    // gen-rec level correlation for b-quarks Pt
    //corrs_["bqPt_"       ] = new TH2F( "bqPt_"       , "bqPt_"       ,  1200,  0. ,  1200. ,  1200,  0. ,  1200.);
    // gen-rec level correlation for b-quarks Eta
    //corrs_["bqEta_"      ] = new TH2F( "bqEta_"      , "bqEta_"      ,  100,  -5. ,  5.    ,  100,  -5. ,  5.   );
    // gen-rec level correlation for leading quark Pt
    //corrs_["leadqPt_"    ] = new TH2F( "leadqPt_"    , "leadqPt_"    ,  1200,  0. ,  1200. ,  1200,  0. ,  1200.);
    // gen-rec level correlation for leading quark Eta
    //corrs_["leadqEta_"   ] = new TH2F( "leadqEta_"   , "leadqEta_"   ,  100,  -5. ,  5.    ,  100,  -5. ,  5.   );
    // ttbar phiStar
    corrs_["ttbarPhiStar_"] = new TH2F( "ttbarPhiStar_" , "ttbarPhiStar_" , 200, 0., 2., 200, 0., 2.); 

    // splitting into leading and subleading top quarks
    corrs_["topPtLead_"      ] = new TH2F( "topPtLead_"    , "topPtLead_"    , 800 ,   0. , 800.,  800,   0. , 800.);
    corrs_["topPtSubLead_"   ] = new TH2F( "topPtSubLead_" , "topPtSubLead_" , 800 ,   0. , 800.,  800,   0. , 800.);
    corrs_["lepYLead_"       ] = new TH2F( "lepYLead_"     , "lepYLead_"     ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    corrs_["lepYSubLead_"    ] = new TH2F( "lepYSubLead_"  , "lepYSubLead_"  ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    // pt of the top candidates in the ttbar restframe
    corrs_["topPtTtbarSys_"  ] = new TH2F( "topPtTtbarSys_"    , "topPtTtbarSys_"    , 800 ,   0. , 800.,  800,   0. , 800.);
    // b-top hadronization variable
    //corrs_["xB_"             ] = new TH2F( "xB_"           , "xB_"           , 200 ,  -0.5 ,  1.5 ,  200, -0.5 , 1.5 );
    // asymmetry variables
    //corrs_["topPtPlus_"      ] = new TH2F( "topPtPlus_"    , "topPtPlus_"    , 800 ,   0. , 800.,  800,   0. , 800.);
    //corrs_["topPtMinus_"     ] = new TH2F( "topPtMinus_"   , "topPtMinus_"   , 800 ,   0. , 800.,  800,   0. , 800.);
    //corrs_["lepPtPlus_"      ] = new TH2F( "lepPtPlus_"    , "lepPtPlus_"    , 1200 ,  0. , 1200., 1200,  0. ,1200.);
    //corrs_["lepPtMinus_"     ] = new TH2F( "lepPtMinus_"   , "lepPtMinus_"   , 1200 ,  0. , 1200., 1200,  0. ,1200.);
    //corrs_["lepEtaPlus_"     ] = new TH2F( "lepEtaPlus_"   , "lepEtaPlus_"   ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    //corrs_["lepEtaMinus_"    ] = new TH2F( "lepEtaMinus_"  , "lepEtaMinus_"  ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    //corrs_["lepYPlus_"       ] = new TH2F( "lepYPlus_"     , "lepYPlus_"     ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    //corrs_["lepYMinus_"      ] = new TH2F( "lepYMinus_"    , "lepYMinus_"    ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    //corrs_["topEtaPlus_"     ] = new TH2F( "topEtaPlus_"   , "topEtaPlus_"   ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    //corrs_["topEtaMinus_"    ] = new TH2F( "topEtaMinus_"  , "topEtaMinus_"  ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    //corrs_["topYPlus_"       ] = new TH2F( "topYPlus_"     , "topYPlus_"     ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    //corrs_["topYMinus_"      ] = new TH2F( "topYMinus_"    , "topYMinus_"    ,  100,  -5. ,  5. ,  100,  -5. ,  5. );

  }
}

/// histogramm booking for fw
void TopKinematics::book(edm::Service<TFileService>& fs)
{

  /** 
      ttbar decay setup
  **/
  // decay channel
  hists_["decayChannel"] = fs->make<TH1F>( "decayChannel", "decayChannel",  25,   -4.5,   20.5);
  // set labels (decayChannel=-4 corresponds to bin 1)
  hists_.find("decayChannel")->second->GetXaxis()->SetBinLabel(1 , "??"      );
  hists_.find("decayChannel")->second->GetXaxis()->SetBinLabel(2 , "ntt"     );
  hists_.find("decayChannel")->second->GetXaxis()->SetBinLabel(3 , "ll?"     );
  hists_.find("decayChannel")->second->GetXaxis()->SetBinLabel(4 , "l?j"     );
  hists_.find("decayChannel")->second->GetXaxis()->SetBinLabel(6 , "ej"      );
  hists_.find("decayChannel")->second->GetXaxis()->SetBinLabel(7 , "#muj"    );
  hists_.find("decayChannel")->second->GetXaxis()->SetBinLabel(8 , "#tauj"   );
  hists_.find("decayChannel")->second->GetXaxis()->SetBinLabel(16, "ee"      );
  hists_.find("decayChannel")->second->GetXaxis()->SetBinLabel(17, "e#mu"    );
  hists_.find("decayChannel")->second->GetXaxis()->SetBinLabel(18, "e#tau"   );
  hists_.find("decayChannel")->second->GetXaxis()->SetBinLabel(19, "#mu#mu"  );
  hists_.find("decayChannel")->second->GetXaxis()->SetBinLabel(20, "#mu#tau" );
  hists_.find("decayChannel")->second->GetXaxis()->SetBinLabel(21, "#tau#tau");
  hists_.find("decayChannel")->second->GetXaxis()->SetBinLabel(25, "jj"      );
  /** 
      KinFit Monitoring Variables
  **/
  // fit probability of the best fit hypothesis
  hists_["prob"       ] = fs->make<TH1F>( "prob"        , "prob"        , 1000,   0.,    1. );
  // chi2 of the best fit hypothesis
  hists_["chi2"       ] = fs->make<TH1F>( "chi2"        , "chi2"        , 1000,   0.,   100.);
  // delta chi2 between best and second best fit hyothesis
  hists_["delChi2"    ] = fs->make<TH1F>( "delChi2"     , "delChi2"     , 1000,   0.,   100.);

 /** 
      Top Variables for Cross Section Measurement
  **/
  // top pt (at the moment both top candidates are filled in one histogram)
  hists_["topPt"      ] = fs->make<TH1F>( "topPt"      , "topPt"      ,  800,  0. ,  800.);
  // top y (at the moment both top candidates are filled in one histogram)
  hists_["topY"       ] = fs->make<TH1F>( "topY"       , "topY"       ,  100, -5. ,    5.);
  // top phi (at the moment both top candidates are filled in one histogram)
  hists_["topPhi"     ] = fs->make<TH1F>( "topPhi"     , "topPhi"     ,  628, -pi ,  pi  );
  // ttbar pair pt
  hists_["ttbarPt"    ] = fs->make<TH1F>( "ttbarPt"    , "ttbarPt"    ,  300,  0. ,  300.);
  // ttbar pair 
  hists_["ttbarY"     ] = fs->make<TH1F>( "ttbarY"     , "ttbarY"     ,  100, -5. ,    5.);
  // ttbar pair phi
  hists_["ttbarPhi"   ] = fs->make<TH1F>( "ttbarPhi"   , "ttbarPhi"   ,  628, -pi ,  pi  );
  // ttbar pair invariant mass
  hists_["ttbarMass"  ] = fs->make<TH1F>( "ttbarMass"  , "ttbarMass"  , 2500,  0. , 2500 );
  // deltaPhi between both top quarks
  hists_["ttbarDelPhi"] = fs->make<TH1F>( "ttbarDelPhi", "ttbarDelPhi",  400, 0.  ,    4.);
  // deltaY between both top quarks
  hists_["ttbarDelY"  ] = fs->make<TH1F>( "ttbarDelY"  , "ttbarDelY"  , 1000, -5. ,    5.);
  // sum of y of both top quarks
  hists_["ttbarSumY"  ] = fs->make<TH1F>( "ttbarSumY"  , "ttbarSumY"  ,  100, -5. ,    5.);
  // HT of the 4 jets assigned to the ttbar decay
  hists_["ttbarHT"    ] = fs->make<TH1F>( "ttbarHT"    , "ttbarHT"    , 1500,  0. , 1500.);
  // ttbar phiStar
  hists_["ttbarPhiStar"] = fs->make<TH1F>( "ttbarPhiStar" , "ttbarPhiStar" , 200, 0., 2.);
  
  /**
     Top Variables for Cross Checks
  **/
  // pt of the leading top candidate
  hists_["topPtLead"   ] =fs->make<TH1F>( "topPtLead"   , "topPtLead"   , 800, 0. , 800.);  
  // pt of the subleading top candidate
  hists_["topPtSubLead"] =fs->make<TH1F>( "topPtSubLead", "topPtSubLead", 800, 0. , 800.);
  // pt of the top candidates in the ttbar restframe
  hists_["topPtTtbarSys"] =fs->make<TH1F>( "topPtTtbarSys", "topPtTtbarSys", 800, 0., 800.);
  // y of the leading top candidate
  hists_["topYLead"   ] = fs->make<TH1F>( "topYLead"   , "topYLead"   ,  100, -5. ,   5.);
  // y of the subleading top candidate
  hists_["topYSubLead"] = fs->make<TH1F>( "topYSubLead", "topYSubLead",  100, -5. ,   5.);
  // pt of the hadronically decaying top candidate
  hists_["topPtHad"  ] = fs->make<TH1F>( "topPtHad"  , "topPtHad"  ,  800,  0. ,  800.);
  // y  of the hadronically decaying top candidate
  hists_["topYHad"   ] = fs->make<TH1F>( "topYHad"   , "topYHad"   ,  100, -5. ,    5.);
  // y  of the hadronically decaying top candidate
  hists_["topPhiHad" ] = fs->make<TH1F>( "topPhiHad" , "topPhiHad" ,  628, -pi ,  pi  );
  // pt of the leptonically decaying top candidate
  hists_["topPtLep"  ] = fs->make<TH1F>( "topPtLep"  , "topPtLep"  ,  800,  0. ,  800.);
  // y  of the leptonically decaying top candidate
  hists_["topYLep"   ] = fs->make<TH1F>( "topYLep"   , "topYLep"   ,  100, -5. ,    5.);
  // y  of the leptonically decaying top candidate
  hists_["topPhiLep" ] = fs->make<TH1F>( "topPhiLep" , "topPhiLep" ,  628, -pi ,  pi  );
  // hadronic Top mass
  hists_["hadTopMass"] = fs->make<TH1F>( "hadTopMass", "hadTopMass", 1000,  0.,  1000.);
  // leptonic Top mass								
  hists_["lepTopMass"] = fs->make<TH1F>( "lepTopMass", "lepTopMass", 1000,  0.,  1000.);
  // Top mass									
  hists_["topMass"   ] = fs->make<TH1F>( "topMass"   , "topMass"   , 1000,  0.,  1000.);

  /**
     bquark system
  **/							
  hists_["bbbarPt"   ] = fs->make<TH1F>( "bbbarPt"   , "bbbarPt"   , 1000,  0.,  1000.);
  hists_["bbbarY"    ] = fs->make<TH1F>( "bbbarY"    , "bbbarY"    , 100,  -5.,  5.   );
  hists_["bbbarMass" ] = fs->make<TH1F>( "bbbarMass" , "bbbarMass" , 1500,  0.,  1500.);
  // leptonic b + lepton invariant mass
  hists_["lbMass"    ] = fs->make<TH1F>( "lbMass"    , "lbMass"    , 500,   0.,  500. );
  // b-top hadronization variable
  //hists_["xB"        ] = fs->make<TH1F>( "xB"        , "xB"        , 200, -0.5,  1.5  );

  /**
     asymmetry variables
  **/
  hists_["topPtPlus"      ] =fs->make<TH1F>( "topPtPlus"    , "topPtPlus"    ,  800,   0. ,  800.);
  hists_["topPtMinus"     ] =fs->make<TH1F>( "topPtMinus"   , "topPtMinus"   ,  800,   0. ,  800.);
  hists_["lepPtPlus"      ] =fs->make<TH1F>( "lepPtPlus"    , "lepPtPlus"    , 1200,   0. , 1200.);
  hists_["lepPtMinus"     ] =fs->make<TH1F>( "lepPtMinus"   , "lepPtMinus"   , 1200,   0. , 1200.);
  hists_["lepEtaPlus"     ] =fs->make<TH1F>( "lepEtaPlus"   , "lepEtaPlus"   ,  100,  -5. ,  5.  );
  hists_["lepEtaMinus"    ] =fs->make<TH1F>( "lepEtaMinus"  , "lepEtaMinus"  ,  100,  -5. ,  5.  );
  hists_["lepYPlus"       ] =fs->make<TH1F>( "lepYPlus"     , "lepYPlus"     ,  100,  -5. ,  5.  );
  hists_["lepYMinus"      ] =fs->make<TH1F>( "lepYMinus"    , "lepYMinus"    ,  100,  -5. ,  5.  );
  hists_["topEtaPlus"     ] =fs->make<TH1F>( "topEtaPlus"   , "topEtaPlus"   ,  100,  -5. ,  5.  );
  hists_["topEtaMinus"    ] =fs->make<TH1F>( "topEtaMinus"  , "topEtaMinus"  ,  100,  -5. ,  5.  );
  hists_["topYPlus"       ] =fs->make<TH1F>( "topYPlus"     , "topYPlus"     ,  100,  -5. ,  5.  );
  hists_["topYMinus"      ] =fs->make<TH1F>( "topYMinus"    , "topYMinus"    ,  100,  -5. ,  5.  );

  /**
     Angular distributions
  **/
  // angle between t candidate (detector rest frame)
  hists_["ttbarAngle"       ] = fs->make<TH1F>( "ttbarAngle"       , "ttbarAngle"          , 315,  0.  ,  pi );
  // angle between b candidates (detector rest frame)
  hists_["bbbarAngle"       ] = fs->make<TH1F>( "bbbarAngle"       , "bbbarAngle"          , 315,  0.  ,  pi );
  // angle between b candidates (ttbar rest frame)
  hists_["bbbarAngleTtRF"   ] = fs->make<TH1F>( "bbbarAngleTtRF"   , "bbbarAngleTtR"       , 315,  0.  ,  pi );
  // angle between W bosons (ttbar rest frame)
  hists_["WWAngle"          ] = fs->make<TH1F>( "WWAngle"          , "WWAngle"             , 315,  0.  ,  pi );
  // angle between the leptonically decaying top candidate and the corresponding W 
  // (Top in ttbar rest frame and W in top rest frame)
  hists_["topWAngleLep"     ] = fs->make<TH1F>( "topWAngleLep"     , "topWAngleLep"        , 315,  0.  ,  pi );
  // angle between the hadronically decaying top candidate and the corresponding W 
  // (Top in ttbar rest frame and W in top rest frame)   
  hists_["topWAngleHad"     ] = fs->make<TH1F>( "topWAngleHad"     , "topWAngleHad"        , 315,  0.  ,  pi );
  // angle between the leptonically decaying top candidate and the corresponding b (ttbar rest frame)
  hists_["topBAngleLep"     ] = fs->make<TH1F>( "topBAngleLep"     , "topBAngleLep"        , 315,  0.  ,  pi );
  // angle between the hadronically decaying top candidate and the corresponding b (ttbar rest frame)
  hists_["topBAngleHad"     ] = fs->make<TH1F>( "topBAngleHad"     , "topBAngleHad"        , 315,  0.  ,  pi );
  // angle between the leptonically decaying W candidate and the corresponding b (ttbar rest frame)
  hists_["bWAngleLep"       ] = fs->make<TH1F>( "bWAngleLep"       , "bWAngleLep"          , 315,  0.  ,  pi );
  // angle between the hadronically decaying W candidate and the corresponding b (ttbar rest frame)
  hists_["bWAngleHad"       ] = fs->make<TH1F>( "bWAngleHad"       , "bWAngleHad"          , 315,  0.  ,  pi );
  // angle between light quarks (ttbar rest frame)
  hists_["qqbarAngle"       ] = fs->make<TH1F>( "qqbarAngle"       , "qqbarAngle"          , 315,  0.  ,  pi );
  // angle between light quarks and leptonic b (ttbar rest frame)
  hists_["qBlepAngle"       ] = fs->make<TH1F>( "qBlepAngle"       , "qBlepAngle"          , 315,  0.  ,  pi );
  // angle between light quarks and hadronic b (corresponding top candidate rest frame)
  hists_["qBhadAngle"       ] = fs->make<TH1F>( "qBhadAngle"       , "qBhadAngle"          , 315,  0.  ,  pi );
  // angle between lepton and leptonic b (corresponding top candidate rest frame)
  hists_["lepBlepAngle"     ] = fs->make<TH1F>( "lepBlepAngle"     , "lepBlepAngle"        , 315,  0.  ,  pi );
  // angle between lepton and leptonic b ttbar rest frame)
  hists_["lepBlepAngleTtRF" ] = fs->make<TH1F>( "lepBlepAngleTtRF" , "lepBlepAngleTtRF"    , 315,  0.  ,  pi );
  // angle between lepton and hadronic b (ttbar rest frame)
  hists_["lepBhadAngle"     ] = fs->make<TH1F>( "lepBhadAngle"     , "lepBhadAngle"        , 315,  0.  ,  pi );
  // angle between lepton and light quarks (ttbar rest frame)
  hists_["lepQAngle"        ] = fs->make<TH1F>( "lepQAngle"        , "lepQAngle"           , 315,  0.  ,  pi );
  // angle between muon and neutrino (in corresponding top rest frame)
  hists_["MuonNeutrinoAngle"] = fs->make<TH1F>( "MuonNeutrinoAngle", "MuonNeutrinoAngle"   , 315,  0.  ,  pi );
  // angle between leptonic b and neutrino (in corresponding top rest frame)
  hists_["lepBNeutrinoAngle"] = fs->make<TH1F>( "lepBNeutrinoAngle", "lepBNeutrinoAngle"   , 315,  0.  ,  pi );
  // angle between hadronic b and neutrino (in ttbar rest frame)
  hists_["hadBNeutrinoAngle"] = fs->make<TH1F>( "hadBNeutrinoAngle", "hadBNeutrinoAngle"   , 315,  0.  ,  pi );
  // angle between light quarks and neutrino (in ttbar rest frame)
  hists_["qNeutrinoAngle"   ] = fs->make<TH1F>( "qNeutrinoAngle"   , "qNeutrinoAngle"      , 315,  0.  ,  pi );
  // angle between lepton (in W rest frame) and corresponding W (in detector rest frame)
  hists_["lepWDir"          ] = fs->make<TH1F>( "lepWDir"          , "lepWDir"             , 315,  0.  ,  pi );
  // angle between neutrino (in W rest frame) and corresponding W (in detector rest frame)
  hists_["nuWDir"           ] = fs->make<TH1F>( "nuWDir"           , "nuWDir"              , 315,  0.  ,  pi );
  // angle between light quarks (in W rest frame) and corresponding W (in detector rest frame)
  hists_["qWDir"            ] = fs->make<TH1F>( "qWDir"            , "qWDir"               , 315,  0.  ,  pi );

  /** 
      event shape variables
  **/
  // notation:
  // 0<=q1<=q2<=q3 are the eigenvalues of the momemtum tensor:
  // T_ab = sum{p_j[a]*p_j[b]}/sum{p_j**2} normalized to 1.

  // aplanarity of event: 
  // 1.5*q1 
  // Return values are 0.5 for spherical and 0 for plane and linear events
  hists_["aplanarity"  ] = fs->make<TH1F> ("aplanarity" , "aplanarity" , 50 ,  0. , 0.5);
  // sphericity of event:
  // 1.5*(q1+q2) 
  // Return values are 1 for spherical, 3/4 for plane and 0 for linear events
  hists_["sphericity"  ] = fs->make<TH1F> ("sphericity" , "sphericity" , 100,  0. , 1. );
  // C (three jet structure) of event:
  // 3.*(q1*q2+q1*q3+q2*q3)
  // Return value is between 0 and 1, C vanishes for a "perfect" 2-jet event
  hists_["C"           ] = fs->make<TH1F> ("C"          , "C"          , 100,  0. , 1. ); 
  // D (four jet structure) of event:
  // 27.*(q1*q2*q3)    
  // Return value is between 0 and 1, D vanishes for a planar event
  hists_["D"           ] = fs->make<TH1F> ("D"          , "D"          , 100,  0. , 1. );
  // circularity of event
  // the return value is 1 for spherical and 0 linear events in r-phi
  hists_["circularity" ] = fs->make<TH1F> ("circularity", "circularity", 100,  0. , 1. ); 
  // isotropy of event:
  // the return value is 1 for spherical events and 0 for events linear in r-phi
  hists_["isotropy"    ] = fs->make<TH1F> ("isotropy"   , "isotropy"   , 100,  0. , 1. ); 

  /** 
      ttbar final state object distributions
  **/
  // leptonPt
  hists_["lepPt"      ] = fs->make<TH1F> ( "lepPt"      , "lepPt"      ,  1200,  0. ,  1200.);
  // leptonEta		     
  hists_["lepEta"     ] = fs->make<TH1F> ( "lepEta"     , "lepEta"     ,  100,  -5. ,  5.   );
  // neutrinoPt		
  hists_["neutrinoPt" ] = fs->make<TH1F> ( "neutrinoPt" , "neutrinoPt" ,  1200,  0. ,  1200.);
  // neutrinoEta	
  hists_["neutrinoEta"] = fs->make<TH1F> ( "neutrinoEta", "neutrinoEta",  100,  -5. ,  5.   );
  // light-quarks Pt	
  hists_["lightqPt"   ] = fs->make<TH1F> ( "lightqPt"   , "lightqPt"   ,  1200,  0. ,  1200.);
  // light-quarks Eta	
  hists_["lightqEta"  ] = fs->make<TH1F> ( "lightqEta"  , "lightqEta"  ,  100,  -5. ,  5.   );
  // b-quarks Pt	
  hists_["bqPt"       ] = fs->make<TH1F> ( "bqPt"       , "bqPt"       ,  1200,  0. ,  1200.);
  // b-quarks Eta	
  hists_["bqEta"      ] = fs->make<TH1F> ( "bqEta"      , "bqEta"      ,  100,  -5. ,  5.   );
  // leading quark Pt	
  hists_["leadqPt"    ] = fs->make<TH1F> ( "leadqPt"    , "leadqPt"    ,  1200,  0. ,  1200.);
  // leading quark Eta	
  hists_["leadqEta"   ] = fs->make<TH1F> ( "leadqEta"   , "leadqEta"   ,  100,  -5. ,  5.   );

  /** 
      wrong reconstructed quarks
  **/
  hists_["qAssignment"] = fs->make<TH1F>( "qAssignment", "qAssignment",     10, -0.5,   9.5 );
  // set labels (assignment=0 corresponds to bin 1)
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(1, "ok"      );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(2, "bb"      );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(3, "blepq"   );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(4, "bhadq"   );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(5, "bbqlep"  );    
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(6, "bbqhad"  );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(7, "bbqq"    );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(8, "jmis"    );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(9, "wrongj"  );
  hists_.find("qAssignment")->second->GetXaxis()->SetBinLabel(10,"nomatch" );

  /** 
      Correlation Plots
  **/
  if(!hypoKey_.compare("None")==0){
    // gen-rec level correlation top pt
    corrs_["topPt_"     ] = fs->make<TH2F>( "topPt_"      , "topPt_"     ,  800,    0.,  800.,     800,   0.,  800.);
    // gen-rec level correlation top y
    corrs_["topY_"      ] = fs->make<TH2F>( "topY_"       , "topY_"      ,  1000,  -5.,    5.,    1000,  -5.,    5.);
    // gen-rec level correlation top phi
    //corrs_["topPhi_"    ] = fs->make<TH2F>( "topPhi_"     , "topPhi_"    ,  628,   -pi,   pi ,     628, -pi ,   pi );
    // gen-rec level correlation ttbar pt
    corrs_["ttbarPt_"   ] = fs->make<TH2F>( "ttbarPt_"    , "ttbarPt_"   ,  300,    0.,  300.,     300,   0.,  300.);
    // gen-rec level correlation ttbar y
    corrs_["ttbarY_"    ] = fs->make<TH2F>( "ttbarY_"     , "ttbarY_"    , 1000,   -5.,    5.,    1000,  -5.,    5.);
    // gen-rec level correlation ttbar mass
    corrs_["ttbarMass_" ] = fs->make<TH2F>( "ttbarMass_"  , "ttbarMass_" , 2500,    0., 2500.,    2500,   0., 2500.);
    // gen-rec level correlation HT of the 4 jets assigned to the ttbar decay
    //corrs_["ttbarHT_"   ] = fs->make<TH2F>( "ttbarHT_"    , "ttbarHT_"   ,  150,    0., 1500.,     150,   0., 1500.);
    // gen-rec level correlation ttbar deltaPhi
    corrs_["ttbarDelPhi_"]= fs->make<TH2F>( "ttbarDelPhi_", "ttbarDelPhi_",  400,   0.,     4.,     400,  0.,   4.);
    // gen-rec level correlation ttbar deltaY
    corrs_["ttbarDelY_" ] = fs->make<TH2F>( "ttbarDelY_"  , "ttbarDelY_"  , 1000,   -5.,    5.,    1000,  -5.,   5.);
    // gen-rec level correlation ttbar sumY
    //corrs_["ttbarSumY_" ] = fs->make<TH2F>( "ttbarSumY_"  , "ttbarSumY_" ,   50,   -5.,    5.,      50,  -5.,    5.);
    // gen-rec level correlation ttbar phiStar
    corrs_["ttbarPhiStar_"] = fs->make<TH2F>( "ttbarPhiStar_" , "ttbarPhiStar_" , 200, 0., 2., 200, 0., 2.);  
    // gen-rec level correlation angle between b jets
    //corrs_["bbbarAngle_"] = fs->make<TH2F>( "bbbarAngle_" , "bbbarAngle_",  315,    0.,    pi,     315,   0.,   pi );
    // gen-rec level correlation for lepton charge
    //corrs_["lepCharge_" ] = fs->make<TH2F>( "lepCharge_"  , "lepCharge_" ,    3,  -1.5,   1.5,       3, -1.5,   1.5);
    // gen-rec level correlation for angle between the leptonically decaying top candidate and the neutrino
    //corrs_["MuonNeutrinoAngle_"] = fs->make<TH2F>( "MuonNeutrinoAngle_", "MuonNeutrinoAngle_", 315,  0. ,  pi,   315,  0. ,  pi);
    // gen-rec level correlation for leptonPt
    //corrs_["lepPt_"      ] = fs->make<TH2F>( "lepPt_"      , "lepPt_"      ,  1200,  0. ,  1200. ,  1200,  0. ,  1200.);
    // gen-rec level correlation for leptonEta
    //corrs_["lepEta_"     ] = fs->make<TH2F>( "lepEta_"     , "lepEta_"     ,  100,  -5. ,  5.    ,  100,  -5. ,  5.   );
    // gen-rec level correlation for neutrinoPt
    //corrs_["neutrinoPt_" ] = fs->make<TH2F>( "neutrinoPt_" , "neutrinoPt_" ,  120,  0. ,  1200.  ,  120,  0.  ,  1200.);
    // gen-rec level correlation for neutrinoEta
    //corrs_["neutrinoEta_"] = fs->make<TH2F>( "neutrinoEta_", "neutrinoEta_",  100,  -5. ,  5.    ,  100,  -5. ,  5.   );
    // gen-rec level correlation for light-quarks Pt
    //corrs_["lightqPt_"   ] = fs->make<TH2F>( "lightqPt_"   , "lightqPt_"   ,  1200,  0. ,  1200. ,  1200,  0. ,  1200.);
    // gen-rec level correlation for light-quarks Eta
    //corrs_["lightqEta_"  ] = fs->make<TH2F>( "lightqEta_"  , "lightqEta_"  ,  100,  -5. ,  5.    ,  100,  -5. ,  5.   );
    // gen-rec level correlation for b-quarks Pt
    //corrs_["bqPt_"       ] = fs->make<TH2F>( "bqPt_"       , "bqPt_"       ,  1200,  0. ,  1200. ,  1200,  0. ,  1200.);
    // gen-rec level correlation for b-quarks Eta
    //corrs_["bqEta_"      ] = fs->make<TH2F>( "bqEta_"      , "bqEta_"      ,  100,  -5. ,  5.    ,  100,  -5. ,  5.   );
    // gen-rec level correlation for leading quark Pt
    //corrs_["leadqPt_"    ] = fs->make<TH2F>( "leadqPt_"    , "leadqPt_"    ,  1200,  0. ,  1200. ,  1200,  0. ,  1200.);
    // gen-rec level correlation for leading quark Eta
    //corrs_["leadqEta_"   ] = fs->make<TH2F>( "leadqEta_"   , "leadqEta_"   ,  100,  -5. ,  5.    ,  100,  -5. ,  5.   );

    // splitting into leading and subleading top quarks
    corrs_["topPtLead_"   ] = fs->make<TH2F>("topPtLead_"    , "topPtLead_"    ,  800,   0. , 800.,  800,   0. , 800.);
    corrs_["topPtSubLead_"] = fs->make<TH2F>("topPtSubLead_" , "topPtSubLead_" ,  800,   0. , 800.,  800,   0. , 800.);
    corrs_["topYLead_"    ] = fs->make<TH2F>("topYLead_"     , "topYLead_"     ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    corrs_["topYSubLead_" ] = fs->make<TH2F>("topYSubLead_"  , "topYSubLead_"  ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    // pt of the top candidates in the ttbar restframe
    corrs_["topPtTtbarSys_" ] = fs->make<TH2F>( "topPtTtbarSys_", "topPtTtbarSys_", 800, 0., 800. , 800, 0., 800.);
    // b-top hadronization variable xB=(1/(1-mW^1/mt^2)) * 2pBpT/mt^2
    //corrs_["xB_"            ] = fs->make<TH2F>( "xB_"           , "xB_"           , 200, -0.5, 1.5, 200, -0.5, 1.5);
    // asymmetry variables
    //corrs_["topPtPlus_"  ] = fs->make<TH2F>( "topPtPlus_"    , "topPtPlus_"    ,  800,   0. , 800.,  800,   0. , 800.);
    //corrs_["topPtMinus_" ] = fs->make<TH2F>( "topPtMinus_"   , "topPtMinus_"   ,  800,   0. , 800.,  800,   0. , 800.);
    //corrs_["lepPtPlus_"  ] = fs->make<TH2F>( "lepPtPlus_"    , "lepPtPlus_"    , 1200 ,  0. , 1200., 1200,  0. ,1200.);
    //corrs_["lepPtMinus_" ] = fs->make<TH2F>( "lepPtMinus_"   , "lepPtMinus_"   , 1200 ,  0. , 1200., 1200,  0. ,1200.);
    //corrs_["lepEtaPlus_" ] = fs->make<TH2F>( "lepEtaPlus_"   , "lepEtaPlus_"   ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    //corrs_["lepEtaMinus_"] = fs->make<TH2F>( "lepEtaMinus_"  , "lepEtaMinus_"  ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    //corrs_["lepYPlus_"   ] = fs->make<TH2F>( "lepYPlus_"     , "lepYPlus_"     ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    //corrs_["lepYMinus_"  ] = fs->make<TH2F>( "lepYMinus_"    , "lepYMinus_"    ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    //corrs_["topEtaPlus_" ] = fs->make<TH2F>( "topEtaPlus_"   , "topEtaPlus_"   ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    //corrs_["topEtaMinus_"] = fs->make<TH2F>( "topEtaMinus_"  , "topEtaMinus_"  ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    //corrs_["topYPlus_"   ] = fs->make<TH2F>( "topYPlus_"     , "topYPlus_"     ,  100,  -5. ,  5. ,  100,  -5. ,  5. );
    //corrs_["topYMinus_"  ] = fs->make<TH2F>( "topYMinus_"    , "topYMinus_"    ,  100,  -5. ,  5. ,  100,  -5. ,  5. );

  }

  // book ttree entries
  if(useTree_){
    // ttbar decay channel
    bookVariable(fs, "decayChannel");
    // Kinfit performance
    bookVariable(fs, "prob"   );
    bookVariable(fs, "chi2"   );
    bookVariable(fs, "delChi2");
    // ttbar quantities
    bookVariable(fs, "topPtTtbarSys");
    bookVariable(fs, "ttbarPt"    );
    bookVariable(fs, "ttbarY"     );
    bookVariable(fs, "ttbarPhi"   );
    bookVariable(fs, "ttbarMass"  );
    bookVariable(fs, "ttbarDelPhi");
    bookVariable(fs, "ttbarDelY"  );
    bookVariable(fs, "ttbarSumY"  );
    bookVariable(fs, "ttbarHT"    );
    bookVariable(fs, "ttbarPhiStar");
    // hadronic top quantities
    bookVariable(fs, "topPtHad"  );
    bookVariable(fs, "topYHad"   );
    bookVariable(fs, "topPhiHad" );
    bookVariable(fs, "hadTopMass");
    //bookVariable(fs, "xBHad"     );
    // leptonic top quantities
    bookVariable(fs, "topPtLep"  );
    bookVariable(fs, "topYLep"   );  
    bookVariable(fs, "topPhiLep" );
    bookVariable(fs, "lepTopMass");
    //bookVariable(fs, "xBLep"     );
    // top quantities (= top plus)
    bookVariable(fs, "topPtPlus" );
    bookVariable(fs, "topEtaPlus");
    bookVariable(fs, "topYPlus"  );
    // antitop quantities (= top minus)
    bookVariable(fs, "topPtMinus" );
    bookVariable(fs, "topEtaMinus");
    bookVariable(fs, "topYMinus"  );
    // charge
    bookVariable(fs, "lepCharge"); 
    // angles
    bookVariable(fs, "topAngle"       );
    bookVariable(fs, "ttbarAngle"       );
    bookVariable(fs, "bbbarAngle"       );
    bookVariable(fs, "bbbarAngleTtRF"   );
    bookVariable(fs, "WWAngle"          );
    bookVariable(fs, "topWAngleLep"     );
    bookVariable(fs, "topWAngleHad"     );
    bookVariable(fs, "topBAngleLep"     );
    bookVariable(fs, "topBAngleHad"     );
    bookVariable(fs, "bWAngleLep"       );
    bookVariable(fs, "bWAngleHad"       );
    bookVariable(fs, "qqbarAngle"       );
    bookVariable(fs, "qBlepAngle"       );
    bookVariable(fs, "qBhadAngle"       );
    bookVariable(fs, "lepBlepAngle"     );
    bookVariable(fs, "lepBhadAngle"     );
    bookVariable(fs, "lepQAngle"        );
    bookVariable(fs, "MuonNeutrinoAngle");
    bookVariable(fs, "lepBNeutrinoAngle");
    bookVariable(fs, "hadBNeutrinoAngle");
    bookVariable(fs, "qNeutrinoAngle"   );
    bookVariable(fs, "lepWDir"          );
    bookVariable(fs, "nuWDir"           );
    bookVariable(fs, "qWDir"            );
    bookVariable(fs, "qAssignment"      );
    bookVariable(fs, "aplanarity"       );
    bookVariable(fs, "sphericity"       );
    bookVariable(fs, "circularity"      );
    bookVariable(fs, "isotropy"         ); 
    bookVariable(fs, "C"                );
    bookVariable(fs, "D"                ); 
    // final state objects
    bookVariable(fs, "lepPt"       );
    bookVariable(fs, "lepEta"      );
    bookVariable(fs, "neutrinoPt"  );
    bookVariable(fs, "neutrinoEta" );
    bookVariable(fs, "lightQPt"    );
    bookVariable(fs, "lightQbarPt" );
    bookVariable(fs, "bqPtLep"     );
    bookVariable(fs, "bqPtHad"     );  
    bookVariable(fs, "lightQEta"   );
    bookVariable(fs, "lightQbarEta");
    bookVariable(fs, "bqEtaLep"    );
    bookVariable(fs, "bqEtaHad"    );  
    // bbbar quantities    							
    bookVariable(fs, "bbbarPt"     );
    bookVariable(fs, "bbbarY"      );
    bookVariable(fs, "bbbarMass"   );
    bookVariable(fs, "lbMass"      );
    // parton truth value
    // ttbar quantities    
    bookVariable(fs, "topPtTtbarSysPartonTruth");
    bookVariable(fs, "ttbarPtPartonTruth"    );
    bookVariable(fs, "ttbarYPartonTruth"     );
    bookVariable(fs, "ttbarPhiPartonTruth"   );
    bookVariable(fs, "ttbarMassPartonTruth"  );
    bookVariable(fs, "ttbarDelPhiPartonTruth");
    bookVariable(fs, "ttbarDelYPartonTruth"  );
    bookVariable(fs, "ttbarSumYPartonTruth"  );
    bookVariable(fs, "ttbarHTPartonTruth"    );
    bookVariable(fs, "ttbarPhiStarPartonTruth"); 
    // hadronic top quantities
    bookVariable(fs, "topPtHadPartonTruth"  );
    bookVariable(fs, "topYHadPartonTruth"   );
    bookVariable(fs, "topPhiHadPartonTruth" );
    bookVariable(fs, "hadTopMassPartonTruth");
    //bookVariable(fs, "xBHadPartonTruth"     );
    // leptonic top quantities
    bookVariable(fs, "topPtLepPartonTruth"  );
    bookVariable(fs, "topYLepPartonTruth"   );  
    bookVariable(fs, "topPhiLepPartonTruth" );
    bookVariable(fs, "lepTopMassPartonTruth");
    //bookVariable(fs, "xBLepPartonTruth"     );
    // lepton quantities
    bookVariable(fs, "lepPtPartonTruth"     );
    bookVariable(fs, "lepEtaPartonTruth"    );
    // neutrino quantities
    bookVariable(fs, "neutrinoPtPartonTruth"  );
    bookVariable(fs, "neutrinoEtaPartonTruth" );
    // quark quantities
    bookVariable(fs, "lightQPtPartonTruth"    );
    bookVariable(fs, "lightQbarPtPartonTruth" );
    bookVariable(fs, "bqPtLepPartonTruth"      );
    bookVariable(fs, "bqPtHadPartonTruth"      );  
    bookVariable(fs, "lightQEtaPartonTruth"   );
    bookVariable(fs, "lightQbarEtaPartonTruth");
    bookVariable(fs, "bqEtaLepPartonTruth"     );
    bookVariable(fs, "bqEtaHadPartonTruth"     );   
    // bbbar quantities							
    bookVariable(fs, "bbbarPtPartonTruth"     );
    bookVariable(fs, "bbbarYPartonTruth"      );
    bookVariable(fs, "bbbarMassPartonTruth"   );
    bookVariable(fs, "lbMassPartonTruth"      );
    // top quantities (= top plus)
    bookVariable(fs, "topEtaPlusPartonTruth");
    bookVariable(fs, "topPtPlusPartonTruth" );
    bookVariable(fs, "topYPlusPartonTruth"  );
    // antitop quantities (= top minus)
    bookVariable(fs, "topEtaMinusPartonTruth");
    bookVariable(fs, "topPtMinusPartonTruth" );
    bookVariable(fs, "topYMinusPartonTruth"  );
    // top/antitop association with lep/had top
    bookVariable(fs, "lepTopIsTopPlus");
    // (sub)leading top quarks
    bookVariable(fs, "topPtLead"   );
    bookVariable(fs, "topPtSubLead");
    bookVariable(fs, "topYLead"    );
    bookVariable(fs, "topYSubLead" );
    bookVariable(fs, "topPtLeadPartonTruth"   );
    bookVariable(fs, "topPtSubLeadPartonTruth");
    bookVariable(fs, "topYLeadPartonTruth"    );
    bookVariable(fs, "topYSubLeadPartonTruth" );
  }
}


/// histogram filling interface for generator level for access with fwlite or full framework
void
TopKinematics::fill(const TtGenEvent& tops, const double& weight)
{
  if(useTree_) initializeTrees(-9999, weight);
  // make sure to have a ttbar pair belonging to the semi-leptonic decay channel with 
  // a muon in the final state and neglect events where top decay is not via Vtb
  if( ( (tops.isSemiLeptonic(WDecay::kMuon)&&lepton_.compare("muon")    ==0) || 
        (tops.isSemiLeptonic(WDecay::kElec)&&lepton_.compare("electron")==0) || 
         tops.isSemiLeptonic(WDecay::kTau)                                     ) 
      && tops.leptonicDecayB() 
      && tops.hadronicDecayB() ){
    // define leptonic/hadronic or positive/negative charged objects (B,W,t)
    bool switchLepAndHadTop = false;
    // if ttbarInsteadOfLepHadTop_ == true:
    // lepTop = Top     (positive charge)
    // hadTop = AntiTop (negative charge)
    if(((reco::LeafCandidate*)(tops.singleLepton()))->charge()<0) switchLepAndHadTop=true;
    const reco::GenParticle *lepTop= (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.hadronicDecayTop() : tops.leptonicDecayTop();
    const reco::GenParticle *hadTop= (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.leptonicDecayTop() : tops.hadronicDecayTop();
    const reco::GenParticle *lepW  = (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.hadronicDecayW  () : tops.leptonicDecayW  ();
    const reco::GenParticle *hadW  = (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.leptonicDecayW  () : tops.hadronicDecayW  ();
    const reco::GenParticle *lepB  = (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.hadronicDecayB  () : tops.leptonicDecayB  (); 
    const reco::GenParticle *hadB  = (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.leptonicDecayB  () : tops.hadronicDecayB  ();
    const reco::GenParticle *topPlus = switchLepAndHadTop ? tops.hadronicDecayTop() : tops.leptonicDecayTop();
    const reco::GenParticle *topMinus= switchLepAndHadTop ? tops.leptonicDecayTop() : tops.hadronicDecayTop();
    //std::cout << "hadronic b-quark: " << hadB->status() << std::endl;
    //std::cout << "leptonic b-quark:" << lepB->status() << std::endl;
    //std::cout << "light quark1:" << hadB->status() << std::endl;
    //std::cout << "light quark2:" << hadB->status() << std::endl;

    // define generated scalar sum of all jet pts
    double HT = lepB->pt() + hadB->pt()+tops.hadronicDecayQuark()->pt()+tops.hadronicDecayQuarkBar()->pt();
    // ---
    //    fill all 1D histos
    // ---
    fill(lepTop, hadTop, topPlus, topMinus, lepW, hadW, lepB, hadB, HT, tops.singleLepton()->charge(), weight);
    // ---
    //    fill 1D angle histos
    // ---
    const reco::GenParticle *lep   = tops.singleLepton();
    if(tops.isSemiLeptonic(WDecay::kTau)&&getFinalStateLepton(*tops.singleLepton())) lep=(const reco::GenParticle *) getFinalStateLepton(*tops.singleLepton());
    fillAngles(tops.hadronicDecayB()->p4(), tops.hadronicDecayQuark()->p4(), tops.hadronicDecayQuarkBar()->p4(),
	       tops.leptonicDecayB()->p4(), lep                      ->p4(), tops.singleNeutrino()       ->p4(),
	       weight);
    // ---
    //    fill  final state objects histo
    // ---
    fillFinalStateObjects(tops.hadronicDecayB()->p4(), tops.hadronicDecayQuark()->p4(), tops.hadronicDecayQuarkBar()->p4(),
			  tops.leptonicDecayB()->p4(), lep                      ->p4(), tops.singleNeutrino()       ->p4(),
			  lep->charge(),
			  weight);
    // save lepton charge
    fillValue( "lepCharge", ((reco::LeafCandidate*)(tops.singleLepton()))->charge(), weight );
    // decay channel
    fillValue( "decayChannel", getDecayChannel(tops), weight );
    // fill the tree, if any variable exists to put in
    if(treeVars_.size()) tree->Fill();
  }
}

/// histogram filling interface for reconstruction level for access with fwlite or full framework
void
TopKinematics::fill(const TtSemiLeptonicEvent& tops, const double& weight)
{

  double charge           =-9999;
  double assignment       =-9999;
  double genTtbarPt       =-9999;
  double genTtbarRapidity =-9999;
  double genTtbarPhi      =-9999;
  double genTtbarMass     =-9999;
  double genDeltaPhi      =-9999;
  double genDeltaRapidity =-9999;
  double genLepPt         =-9999;
  double genLepEta        =-9999;
  double genBbbarPt       =-9999;
  double genBbbarY        =-9999;
  double genBbbarMass     =-9999; 
  double genlbMass        =-9999;
  double genNuPt          =-9999;
  double genNuEta         =-9999;
  double genQPt           =-9999;
  double genQEta          =-9999;
  double genQbarPt        =-9999;
  double genQbarEta       =-9999;
  double genLepBPt        =-9999;
  double genLepBEta       =-9999;
  double genHadBPt        =-9999;
  double genHadBEta       =-9999;
  double sumRapidity      =-9999;
  double HTgen            =-9999;
  double hadTopGenPt      =-9999;
  double hadTopGenRapidity=-9999;
  double hadTopGenPhi     =-9999;
  double hadTopGenMass    =-9999;
  //double hadxBGen         =-9999;
  double lepTopGenPt      =-9999;
  double lepTopGenRapidity=-9999;  
  double lepTopGenPhi     =-9999;
  double lepTopGenMass    =-9999;
  //double lepxBGen         =-9999;
  double HTrec            =-9999;
  double prob             =-9999;
  double chi2             =-9999;
  double delChi2          =-9999;
  double genTopPlusY      =-9999;
  double genTopMinusY     =-9999;
  double genTopPlusEta    =-9999;
  double genTopMinusEta   =-9999;
  double genTopPlusPt     =-9999;
  double genTopMinusPt    =-9999;
  double genptLead        =-9999;
  double genptSubLead     =-9999;
  double genyLead         =-9999;
  double genySubLead      =-9999;
  double gentopPtTtbarSystem =-9999;
  double genTtbarPhiStar=-9999;
  double ttbarPhiStarRec=-9999;

  if(useTree_) initializeTrees(-9999, weight);
  // ttbar decay channel
  if(tops.genEvent().isAvailable()) fillValue( "decayChannel", getDecayChannel(*tops.genEvent()), weight );
  else fillValue( "decayChannel", -4, weight );
  // define leptonic/hadronic or positive/negative charged objects (B,W,t)
  bool switchLepAndHadTop = false; // default: false, leptonic/hadronic top
  // make sure to have a valid hypothesis on reconstruction level.
  if( tops.isHypoValid(hypoKey_) ){
    // if ttbarInsteadOfLepHadTop_ == true:
    // lepTop = Top     (positive charge)
    // lepTop = AntiTop (negative charge)
    if(((reco::LeafCandidate*)(tops.singleLepton(hypoKey_)))->charge()<0) switchLepAndHadTop=true;
    const reco::Candidate *lepTopRec= (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.hadronicDecayTop(hypoKey_) : tops.leptonicDecayTop(hypoKey_);
    const reco::Candidate *hadTopRec= (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.leptonicDecayTop(hypoKey_) : tops.hadronicDecayTop(hypoKey_);
    const reco::Candidate *lepWRec  = (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.hadronicDecayW  (hypoKey_) : tops.leptonicDecayW  (hypoKey_);
    const reco::Candidate *hadWRec  = (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.leptonicDecayW  (hypoKey_) : tops.hadronicDecayW  (hypoKey_);
    const reco::Candidate *lepBRec  = (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.hadronicDecayB  (hypoKey_) : tops.leptonicDecayB  (hypoKey_); 
    const reco::Candidate *hadBRec  = (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.leptonicDecayB  (hypoKey_) : tops.hadronicDecayB  (hypoKey_);
    const reco::Candidate *topPlusRec  = switchLepAndHadTop ? tops.hadronicDecayTop(hypoKey_) : tops.leptonicDecayTop(hypoKey_);
    const reco::Candidate *topMinusRec = switchLepAndHadTop ? tops.leptonicDecayTop(hypoKey_) : tops.hadronicDecayTop(hypoKey_);

    // define reconstructed scalar sum of all jet pts
    HTrec = lepBRec->pt() + hadBRec->pt() + tops.hadronicDecayQuark(hypoKey_)->pt() + tops.hadronicDecayQuarkBar(hypoKey_)->pt();
    // if the kinFit hypothesis is valid, fill KinFit quantities 
    if(hypoKey_=="kKinFit"){
      // fit probability of the best fit hypothesis
      if(ndof<0)prob=tops.fitProb();
      else prob=TMath::Prob(tops.fitChi2(),ndof);
      // chi2 of the best fit hypothesis
      chi2= tops.fitChi2();
      // make sure that a second best fit hypothesis exists to be able to fill these plots
      if( tops.fitChi2(1) >= 0 ){
	// delta chi2 between best and second best fit hyothesis
	delChi2=tops.fitChi2(1)-tops.fitChi2(0);
      }
    }
    // create indicator if plots are already filled
    bool oneDplotsFilled=false;
    // check if a generated ttbar semileptonic #mu event exists
    // neglect events where top decay is not via Vtb
    if( tops.genEvent().isAvailable() && 
        ( (tops.genEvent()->isSemiLeptonic(WDecay::kMuon)&&lepton_.compare("muon")    ==0) || 
          (tops.genEvent()->isSemiLeptonic(WDecay::kElec)&&lepton_.compare("electron")==0) ||
	   tops.genEvent()->isSemiLeptonic(WDecay::kTau)                                    ) &&
        lepBRec && hadBRec ){
      const reco::GenParticle *lepTopGen= (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.hadronicDecayTop() : tops.leptonicDecayTop();
      const reco::GenParticle *hadTopGen= (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.leptonicDecayTop() : tops.hadronicDecayTop();
      const reco::GenParticle *lepWGen  = (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.hadronicDecayW  () : tops.leptonicDecayW  ();
      const reco::GenParticle *hadWGen  = (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.leptonicDecayW  () : tops.hadronicDecayW  ();
      const reco::GenParticle *lepBGen  = (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.hadronicDecayB  () : tops.leptonicDecayB  (); 
      const reco::GenParticle *hadBGen  = (ttbarInsteadOfLepHadTop_==true&&switchLepAndHadTop) ? tops.leptonicDecayB  () : tops.hadronicDecayB  ();
      const reco::GenParticle *topPlusGen  = switchLepAndHadTop ? tops.hadronicDecayTop() : tops.leptonicDecayTop();
      const reco::GenParticle *topMinusGen = switchLepAndHadTop ? tops.leptonicDecayTop() : tops.hadronicDecayTop();

      // define generated scalar sum of all jet pts
      HTgen = lepBGen->pt() + hadBGen->pt() + tops.hadronicDecayQuark()->pt() + tops.hadronicDecayQuarkBar()->pt();
      /**
         fill 1D histos and N(gen&&rec)-histo for the determination of stability and purity
	 if matchForStabilityAndPurity_ is true
      **/
      if( matchForStabilityAndPurity_ ){
	fillGenRec(lepTopRec, lepTopGen, hadTopRec, hadTopGen, lepWRec, lepWGen, hadWRec, hadWGen, HTrec, HTgen, weight);
	oneDplotsFilled=true;
      }
      /**
	 fill 2D histos (rec versus gen level correlation plots)
      **/
      // create combined ttbar lorentz vector
      reco::Particle::LorentzVector genTtbar = hadTopGen->p4()+lepTopGen->p4();
      reco::Particle::LorentzVector recTtbar = hadTopRec->p4()+lepTopRec->p4();
      reco::Particle::LorentzVector genBbbar = lepBGen->p4()+hadBGen->p4();
      //reco::Particle::LorentzVector recBbbar = lepBRec->p4()+hadBRec->p4();
 
      /**
         fill parton level truth tree
      **/
      genTtbarPt       = genTtbar.pt();
      genTtbarRapidity = genTtbar.Rapidity();
      genTtbarPhi      = genTtbar.phi();
      genTtbarMass     = genTtbar.mass();
      genDeltaPhi      = std::abs(deltaPhi(lepTopGen->phi(), hadTopGen->phi()));
      genDeltaRapidity = lepTopGen->rapidity()-hadTopGen->rapidity();
      genTtbarPhiStar  = TMath::Tan(0.5*(pi-std::abs(genDeltaPhi)))/TMath::CosH(0.5*(lepTopGen->eta()-hadTopGen->eta()));
      ttbarPhiStarRec  = TMath::Tan(0.5*(pi-std::abs(deltaPhi(lepTopRec->phi(),hadTopRec->phi()))))/TMath::CosH(0.5*(lepTopRec->eta()-hadTopRec->eta()));
      sumRapidity      = lepTopGen->rapidity()+hadTopGen->rapidity();
      hadTopGenPt      = hadTopGen->pt();
      hadTopGenRapidity= hadTopGen->rapidity();
      hadTopGenPhi     = hadTopGen->phi();
      hadTopGenMass    = hadTopGen->mass();
      lepTopGenPt      = lepTopGen->pt();
      lepTopGenRapidity= lepTopGen->rapidity();
      lepTopGenPhi     = lepTopGen->phi();
      lepTopGenMass    = lepTopGen->mass();
      genLepPt         = tops.singleLepton()->pt();
      genLepEta        = tops.singleLepton()->eta();
      if(tops.genEvent().isAvailable()&&((*tops.genEvent()).isSemiLeptonic(WDecay::kTau))&&getFinalStateLepton(*tops.singleLepton())){
	genLepPt         = getFinalStateLepton(*tops.singleLepton())->pt();
	genLepEta        = getFinalStateLepton(*tops.singleLepton())->eta();
      }
      genNuPt          = tops.singleNeutrino()->pt();
      genNuEta         = tops.singleNeutrino()->eta();
      genQPt           = tops.hadronicDecayQuark()->pt();
      genQEta          = tops.hadronicDecayQuark()->eta();
      genQbarPt        = tops.hadronicDecayQuarkBar()->pt();
      genQbarEta       = tops.hadronicDecayQuarkBar()->eta();
      genLepBPt        = tops.leptonicDecayB()->pt();
      genLepBEta       = tops.leptonicDecayB()->eta();
      genHadBPt        = tops.hadronicDecayB()->pt();
      genHadBEta       = tops.hadronicDecayB()->eta();
      genBbbarPt       = genBbbar.pt();
      genBbbarY        = genBbbar.Rapidity();
      genBbbarMass     = genBbbar.mass();
      genlbMass        = (tops.leptonicDecayB()->p4()+tops.singleLepton()->p4()).mass();
      genTopPlusY      = topPlusGen ->rapidity();
      genTopMinusY     = topMinusGen->rapidity();
      genTopPlusEta    = topPlusGen ->eta();
      genTopMinusEta   = topMinusGen->eta();
      genTopPlusPt     = topPlusGen ->pt();
      genTopMinusPt    = topMinusGen->pt();
      genptLead        = topPlusGen ->pt();
      genptSubLead     = topMinusGen->pt();
      genyLead         = topPlusGen ->rapidity();
      genySubLead      = topMinusGen->rapidity();

      // find leading top
      if(genptLead<genptSubLead){// swop them
	double gentemp=genptLead;
	double gentemp2=genyLead;
	genptLead=genptSubLead;
	genyLead=genySubLead;
	genptSubLead=gentemp;
	genySubLead=gentemp2;
      }
      // b-top hadronization variable
      // (1/(1-mW^1/mt^2)) * pBpT/mt^2
      // gen
      //lepxBGen =(1/(1-((lepWGen->p4()).mass2())/((lepTopGen->p4()).mass2()))) * 2*((lepBGen->p4()).Dot(lepTopGen->p4()))/((lepTopGen->p4()).mass2());
      //hadxBGen =(1/(1-((hadWGen->p4()).mass2())/((hadTopGen->p4()).mass2()))) * 2*((hadBGen->p4()).Dot(hadTopGen->p4()))/((hadTopGen->p4()).mass2());
      // rec
      //double lepxBRec =(1/(1-((lepWRec->p4()).mass2())/((lepTopRec->p4()).mass2()))) * 2*((lepBRec->p4()).Dot(lepTopRec->p4()))/((lepTopRec->p4()).mass2());
      //double hadxBRec =(1/(1-((hadWRec->p4()).mass2())/((hadTopRec->p4()).mass2()))) * 2*((hadBRec->p4()).Dot(hadTopRec->p4()))/((hadTopRec->p4()).mass2());
      if(!hypoKey_.compare("None")==0){
	// fill xB correlation plot
	//corrs_.find("xB_")->second->Fill( lepxBGen, lepxBRec, weight );
	//corrs_.find("xB_")->second->Fill( hadxBGen, hadxBRec, weight );
	// fill pt correlation plot for ttbar pair
	corrs_.find("ttbarPt_"  )->second->Fill( genTtbar.pt()      , recTtbar.pt()       , weight );
	// fill y correlation plot for ttbar pair
	corrs_.find("ttbarY_"   )->second->Fill( genTtbar.Rapidity(), recTtbar.Rapidity() , weight );
	// fill mass correlation plot for ttbar pair
	corrs_.find("ttbarMass_")->second->Fill( genTtbar.mass()    , recTtbar.mass()     , weight );    
	// fill HT correlation plot for the 4 jets assigned to the ttbar decay
	//corrs_.find("ttbarHT_"  )->second->Fill( HTgen              , HTrec               , weight ); 
	// fill ttbar phi star 
	corrs_.find("ttbarPhiStar_")->second->Fill( genTtbarPhiStar, ttbarPhiStarRec, weight );
	// fill pt correlation plot for hadronic top candidate
	corrs_.find("topPt_")->second->Fill( hadTopGen->pt(), hadTopRec->pt(), weight );
	// fill pt correlation plot for leptonic top candidate
	corrs_.find("topPt_")->second->Fill( lepTopGen->pt(), lepTopRec->pt(), weight );
	
	// fill y correlation plot for hadronic top candidate
	corrs_.find("topY_")->second->Fill( hadTopGen->rapidity(), hadTopRec->rapidity(), weight );
	// fill y correlation plot for leptonic top candidate
	corrs_.find("topY_")->second->Fill( lepTopGen->rapidity(), lepTopRec->rapidity(), weight );

	// fill phi correlation plot for hadronic top candidate
	//corrs_.find("topPhi_"   )->second->Fill( hadTopGen->phi(), hadTopRec->phi(), weight );
	// fill phi correlation plot for leptonic top candidate							     
	//corrs_.find("topPhi_"   )->second->Fill( lepTopGen->phi(), lepTopRec->phi(), weight );
	// fill asymmetry variables
	//if(tops.singleLepton(hypoKey_)->charge()>0){
	//  corrs_.find("lepEtaPlus_" )->second->Fill( tops.singleLepton()->eta()     , tops.singleLepton(hypoKey_)->eta()     , weight );
	//  corrs_.find("lepYPlus_"   )->second->Fill( tops.singleLepton()->rapidity(), tops.singleLepton(hypoKey_)->rapidity(), weight );
	//  corrs_.find("lepPtPlus_"  )->second->Fill( tops.singleLepton()->pt()      , tops.singleLepton(hypoKey_)->pt()      , weight );
	//}
	//else{			 
	//  corrs_.find("lepEtaMinus_")->second->Fill( tops.singleLepton()->eta()     , tops.singleLepton(hypoKey_)->eta()     , weight );
	//  corrs_.find("lepYMinus_"  )->second->Fill( tops.singleLepton()->rapidity(), tops.singleLepton(hypoKey_)->rapidity(), weight );
	//  corrs_.find("lepPtMinus_" )->second->Fill( tops.singleLepton()->pt()      , tops.singleLepton(hypoKey_)->pt()      , weight );
	//}
	//corrs_.find("topEtaPlus_" )->second->Fill( topPlusGen ->eta()     , topPlusRec ->eta()     , weight );
	//corrs_.find("topEtaMinus_")->second->Fill( topMinusGen->eta()     , topMinusRec->eta()     , weight );
	//corrs_.find("topPtPlus_"  )->second->Fill( topPlusGen ->pt()      , topPlusRec ->pt()      , weight );
	//corrs_.find("topPtMinus_" )->second->Fill( topMinusGen->pt()      , topMinusRec->pt()      , weight );
	//corrs_.find("topYPlus_"   )->second->Fill( topPlusGen ->rapidity(), topPlusRec ->rapidity(), weight );
	//corrs_.find("topYMinus_"  )->second->Fill( topMinusGen->rapidity(), topMinusRec->rapidity(), weight );

	// find leading top
	bool plusGen=(topPlusGen->pt() > topMinusGen->pt() ? true : false);
	bool plusRec=(topPlusRec->pt() > topMinusRec->pt() ? true : false);

	double ptLeadGEN   =( plusGen ? topPlusGen->pt() : topMinusGen->pt());
	double ptSubLeadGEN=(!plusGen ? topPlusGen->pt() : topMinusGen->pt());
	double ptLeadREC   =( plusRec ? topPlusRec->pt() : topMinusRec->pt());
	double ptSubLeadREC=(!plusRec ? topPlusRec->pt() : topMinusRec->pt());
	double yLeadGEN    =( plusGen ? topPlusGen->rapidity() : topMinusGen->rapidity());
	double ySubLeadGEN =(!plusGen ? topPlusGen->rapidity() : topMinusGen->rapidity());
	double yLeadREC    =( plusRec ? topPlusRec->rapidity() : topMinusRec->rapidity());
	double ySubLeadREC =(!plusRec ? topPlusRec->rapidity() : topMinusRec->rapidity());

	corrs_.find("topPtLead_"   )->second->Fill( ptLeadGEN   , ptLeadREC   , weight);
	corrs_.find("topPtSubLead_")->second->Fill( ptSubLeadGEN, ptSubLeadREC, weight);
	corrs_.find("topYLead_"    )->second->Fill( yLeadGEN    , yLeadREC    , weight);
	corrs_.find("topYSubLead_" )->second->Fill( ySubLeadGEN , ySubLeadREC , weight);

	// fill deltaPhi correlation plot for ttbar pair
	if(switchLepAndHadTop) corrs_.find("ttbarDelPhi_")->second->Fill(std::abs(deltaPhi(hadTopGen->phi(), lepTopGen->phi())), std::abs(deltaPhi(hadTopRec->phi(), lepTopRec->phi())), weight);
	else                   corrs_.find("ttbarDelPhi_")->second->Fill(std::abs(deltaPhi(lepTopGen->phi(), hadTopGen->phi())), std::abs(deltaPhi(lepTopRec->phi(), hadTopRec->phi())), weight);
	// fill deltaY correlation plot for ttbar
	if(switchLepAndHadTop) corrs_.find("ttbarDelY_"  )->second->Fill(hadTopGen->rapidity()-lepTopGen->rapidity(), hadTopRec->rapidity()-lepTopRec->rapidity(), weight);
	else                   corrs_.find("ttbarDelY_"  )->second->Fill(lepTopGen->rapidity()-hadTopGen->rapidity(), lepTopRec->rapidity()-hadTopRec->rapidity(), weight);
	// fill sumY correlation plot for ttbar pair
	//corrs_.find("ttbarSumY_"  )->second->Fill(lepTopGen->rapidity()+hadTopGen->rapidity(),
	//lepTopRec->rapidity()+hadTopRec->rapidity(),
	//weight );
	// fill charge correlation plot for lepton
	//corrs_.find("lepCharge_" )->second->Fill( ((reco::LeafCandidate*)(tops.singleLepton(        )))->charge(), 
	//((reco::LeafCandidate*)(tops.singleLepton(hypoKey_)))->charge(), weight);
      }

      // ---
      //    fill correlation for ttbar final state object distributions
      // ---  
      // gen-rec level correlation for leptonPt
      //double pt  = tops.singleLepton()->pt();
      //double eta = tops.singleLepton()->eta();
      //if(tops.genEvent().isAvailable()&&((*tops.genEvent()).isSemiLeptonic(WDecay::kTau))&&getFinalStateLepton(*tops.singleLepton())){
      //	pt  = getFinalStateLepton(*tops.singleLepton())->pt();
      //	eta = getFinalStateLepton(*tops.singleLepton())->eta();
      //}
      //if(!hypoKey_.compare("None")==0){
      //	corrs_.find("lepPt_"      )->second->Fill(pt   , tops.singleLepton(hypoKey_)->pt()   , weight); 
      //	// gen-rec level correlation for leptonEta
      //	corrs_.find("lepEta_"     )->second->Fill(eta  , tops.singleLepton(hypoKey_)->eta()  , weight); 
      //	// gen-rec level correlation for neutrinoPt
      //	corrs_.find("neutrinoPt_" )->second->Fill(tops.singleNeutrino()->pt() , tops.singleNeutrino(hypoKey_)->pt() , weight); 
      //	// gen-rec level correlation for neutrinoEta
      //	corrs_.find("neutrinoEta_")->second->Fill(tops.singleNeutrino()->eta(), tops.singleNeutrino(hypoKey_)->eta(), weight); 
      //}
      // find q/qbar gen-reco association (closest in pt)
      //bool switchLightQuarks=false;
      //if((std::abs(tops.hadronicDecayQuark()->pt()-tops.hadronicDecayQuark(hypoKey_)->pt())+std::abs(tops.hadronicDecayQuarkBar()->pt()-tops.hadronicDecayQuarkBar(hypoKey_)->pt()))>(std::abs(tops.hadronicDecayQuark()->pt()-tops.hadronicDecayQuarkBar(hypoKey_)->pt())+std::abs(tops.hadronicDecayQuarkBar()->pt()-tops.hadronicDecayQuark(hypoKey_)->pt()))){
      //	switchLightQuarks=true;
      //}
      //const reco::GenParticle *QGen   =tops.hadronicDecayQuark   ();
      //const reco::GenParticle *QbarGen=tops.hadronicDecayQuarkBar();
      //const reco::Candidate *QReco   =(switchLightQuarks ? tops.hadronicDecayQuarkBar(hypoKey_) : tops.hadronicDecayQuark   (hypoKey_));
      //const reco::Candidate *QbarReco=(switchLightQuarks ? tops.hadronicDecayQuark   (hypoKey_) : tops.hadronicDecayQuarkBar(hypoKey_));
      //if(!hypoKey_.compare("None")==0){
      //	// gen-rec level correlation for light-quarks Pt
      //	corrs_.find("lightqPt_" )->second->Fill(QGen   ->pt(), QReco   ->pt(), weight);
      //	corrs_.find("lightqPt_" )->second->Fill(QbarGen->pt(), QbarReco->pt(), weight);
      //	// gen-rec level correlation for light-quarks Eta
      //	corrs_.find("lightqEta_")->second->Fill(QGen   ->eta(), QReco   ->eta(), weight);
      //	corrs_.find("lightqEta_")->second->Fill(QbarGen->eta(), QbarReco->eta(), weight);
      //	// gen-rec level correlation for b-quarks Pt
      //	corrs_.find("bqPt_"     )->second->Fill(tops.hadronicDecayB()->pt() , tops.hadronicDecayB(hypoKey_)->pt() , weight);
      //	corrs_.find("bqPt_"     )->second->Fill(tops.leptonicDecayB()->pt() , tops.leptonicDecayB(hypoKey_)->pt() , weight);
      //	// gen-rec level correlation for b-quarks Eta
      //	corrs_.find("bqEta_"    )->second->Fill(tops.hadronicDecayB()->eta(), tops.hadronicDecayB(hypoKey_)->eta(), weight);
      //	corrs_.find("bqEta_"    )->second->Fill(tops.leptonicDecayB()->eta(), tops.leptonicDecayB(hypoKey_)->eta(), weight);
      //}
      // find leading jet
      //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > leadqGen=tops.hadronicDecayB(        )->p4();
      //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > leadqRec=tops.hadronicDecayB(hypoKey_)->p4();
      //if(leadqGen.Pt()<(tops.leptonicDecayB()->pt())){
      //	leadqGen=tops.leptonicDecayB(        )->p4();
      //	leadqRec=tops.leptonicDecayB(hypoKey_)->p4();
      //}
      //if(leadqGen.Pt()<(tops.hadronicDecayQuarkBar()->pt())){
      //	leadqGen=tops.hadronicDecayQuarkBar(        )->p4();
      //	leadqRec=tops.hadronicDecayQuarkBar(hypoKey_)->p4();
      //}
      //if(leadqGen.Pt()<(tops.hadronicDecayQuark()->pt())){
      //	leadqGen=tops.hadronicDecayQuark(        )->p4();
      //	leadqRec=tops.hadronicDecayQuark(hypoKey_)->p4();
      //}
      //if(!hypoKey_.compare("None")==0){
      //	// gen-rec level correlation for leading quark Pt
      //	corrs_.find("leadqPt_"    )->second->Fill(leadqGen.Pt() , leadqRec.Pt() , weight);
      //	// gen-rec level correlation for leading quark Eta
      //	corrs_.find("leadqEta_"   )->second->Fill(leadqGen.Eta(), leadqRec.Eta(), weight);
      //}
      /**
	 fill 2D histos (angle correlation plots)
      **/
      // create boost into reconstructed/generated top quarks/ttbar system
      ROOT::Math::Boost CoMBoostRecHadTop(hadTopRec->p4().BoostToCM());
      ROOT::Math::Boost CoMBoostGenHadTop(hadTopGen->p4().BoostToCM());
      ROOT::Math::Boost CoMBoostRecLepTop(lepTopRec->p4().BoostToCM());
      ROOT::Math::Boost CoMBoostGenLepTop(lepTopGen->p4().BoostToCM());
      ROOT::Math::Boost CoMBoostRecTtbar (recTtbar       .BoostToCM());
      ROOT::Math::Boost CoMBoostGenTtbar (genTtbar       .BoostToCM());
      // get lorentz vectors of W and Top without boost
      reco::Particle::LorentzVector recHadronicDecayWBoosted   = hadWRec  ->p4();
      reco::Particle::LorentzVector genHadronicDecayWBoosted   = hadWGen  ->p4();
      reco::Particle::LorentzVector recLeptonicDecayWBoosted   = lepWRec  ->p4();
      reco::Particle::LorentzVector genLeptonicDecayWBoosted   = lepWGen  ->p4();
      reco::Particle::LorentzVector recHadronicDecayTopBoosted = hadTopRec->p4();
      reco::Particle::LorentzVector genHadronicDecayTopBoosted = hadTopGen->p4();
      reco::Particle::LorentzVector recLeptonicDecayTopBoosted = lepTopRec->p4();
      reco::Particle::LorentzVector genLeptonicDecayTopBoosted = lepTopGen->p4();
      reco::Particle::LorentzVector genLeptonicDecayBBoosted   = lepBGen  ->p4();
      reco::Particle::LorentzVector recLeptonicDecayBBoosted   = lepBRec  ->p4();
      reco::Particle::LorentzVector genHadronicDecayBBoosted   = hadBGen  ->p4();
      reco::Particle::LorentzVector recHadronicDecayBBoosted   = hadBRec  ->p4();

      // recalculate Lorent vectors using boosts defined above
      // boost W into top rest frame and and top and b-quarks into the ttbar rest frame
      recHadronicDecayWBoosted   = CoMBoostRecHadTop(recHadronicDecayWBoosted  );
      genHadronicDecayWBoosted   = CoMBoostGenHadTop(genHadronicDecayWBoosted  );
      recLeptonicDecayWBoosted   = CoMBoostRecLepTop(recLeptonicDecayWBoosted  );
      genLeptonicDecayWBoosted   = CoMBoostGenLepTop(genLeptonicDecayWBoosted  );
      recHadronicDecayTopBoosted = CoMBoostRecTtbar (recHadronicDecayTopBoosted);
      genHadronicDecayTopBoosted = CoMBoostGenTtbar (genHadronicDecayTopBoosted);
      recLeptonicDecayTopBoosted = CoMBoostRecTtbar (recLeptonicDecayTopBoosted);
      genLeptonicDecayTopBoosted = CoMBoostGenTtbar (genLeptonicDecayTopBoosted);
      genLeptonicDecayBBoosted   = CoMBoostGenTtbar (genLeptonicDecayBBoosted  );
      recLeptonicDecayBBoosted   = CoMBoostRecTtbar (recLeptonicDecayBBoosted  );
      genHadronicDecayBBoosted   = CoMBoostGenTtbar (genHadronicDecayBBoosted  );
      recHadronicDecayBBoosted   = CoMBoostRecTtbar (recHadronicDecayBBoosted  );
      gentopPtTtbarSystem=genHadronicDecayTopBoosted.pt();

      // fill correlation plot for top pt in ttbar system
      corrs_.find("topPtTtbarSys_")->second->Fill( gentopPtTtbarSystem, recHadronicDecayTopBoosted.pt(), weight );

      // fill correlation plot for the angle between b jets
      //corrs_.find("bbbarAngle_")->second->Fill( ROOT::Math::VectorUtil::Angle(genLeptonicDecayBBoosted, genHadronicDecayBBoosted), 
      //                                          ROOT::Math::VectorUtil::Angle(recLeptonicDecayBBoosted, recHadronicDecayBBoosted),
      //                                          weight);
      // fill correlation plot for muon - neutrino angle plot
      const reco::GenParticle *lep   = tops.singleLepton();
      if(tops.genEvent().isAvailable()&&((*tops.genEvent()).isSemiLeptonic(WDecay::kTau))&&getFinalStateLepton(*tops.singleLepton())){
	lep=(const reco::GenParticle *) getFinalStateLepton(*tops.singleLepton());
      }
      reco::Particle::LorentzVector genMuonBoosted     = CoMBoostGenTtbar(lep->p4());
      reco::Particle::LorentzVector genNeutrinoBoosted = CoMBoostGenTtbar(tops.singleNeutrino()->p4());
      reco::Particle::LorentzVector recMuonBoosted     = CoMBoostRecTtbar(tops.singleLepton  (hypoKey_)->p4());
      reco::Particle::LorentzVector recNeutrinoBoosted = CoMBoostRecTtbar(tops.singleNeutrino(hypoKey_)->p4());
      //corrs_.find("MuonNeutrinoAngle_") ->second->Fill(ROOT::Math::VectorUtil::Angle(genMuonBoosted, genNeutrinoBoosted), 
      //					   ROOT::Math::VectorUtil::Angle(recMuonBoosted, recNeutrinoBoosted),
      //					   weight);
      // fill 1D plot for angle between b-jets for purity and stability calculation
      if( matchForStabilityAndPurity_ ) match("bbbarAngle", ROOT::Math::VectorUtil::Angle(recLeptonicDecayBBoosted, recHadronicDecayBBoosted),
                                              ROOT::Math::VectorUtil::Angle(genLeptonicDecayBBoosted, genHadronicDecayBBoosted), 
                                              weight);
      // fill 1D plot for muon - neutrino angle plot for purity and stability calculation
      if( matchForStabilityAndPurity_ ) match("MuonNeutrinoAngle", ROOT::Math::VectorUtil::Angle(recMuonBoosted, recNeutrinoBoosted), 
					      ROOT::Math::VectorUtil::Angle(genMuonBoosted, genNeutrinoBoosted),
					      weight);

    }
    // if matchForStabilityAndPurity_=false: 1D plots for all events with valid hypothesis
    if( !oneDplotsFilled ){  
      // fill kinematic 1D histos
      fill(lepTopRec, hadTopRec, topPlusRec, topMinusRec, lepWRec, hadWRec, lepBRec, hadBRec, HTrec, tops.singleLepton(hypoKey_)->charge(), weight);
      // fill angle histos
      fillAngles(tops.hadronicDecayB       (hypoKey_)->p4(), tops.hadronicDecayQuark(hypoKey_)->p4(), 
		 tops.hadronicDecayQuarkBar(hypoKey_)->p4(), tops.leptonicDecayB    (hypoKey_)->p4(), 
		 tops.singleLepton         (hypoKey_)->p4(), tops.singleNeutrino    (hypoKey_)->p4(),
		 weight);
      // fill Event shape histos
      fillEventShapes(tops.hadronicDecayB       (hypoKey_)->p4(), tops.hadronicDecayQuark(hypoKey_)->p4(), 
		      tops.hadronicDecayQuarkBar(hypoKey_)->p4(), tops.leptonicDecayB    (hypoKey_)->p4(), 
		      tops.singleLepton         (hypoKey_)->p4(), tops.singleNeutrino    (hypoKey_)->p4(),
		      weight, false, false);
      // fill final state objects
      fillFinalStateObjects(tops.hadronicDecayB       (hypoKey_)->p4(), tops.hadronicDecayQuark(hypoKey_)->p4(), 
			    tops.hadronicDecayQuarkBar(hypoKey_)->p4(), tops.leptonicDecayB    (hypoKey_)->p4(), 
			    tops.singleLepton         (hypoKey_)->p4(), tops.singleNeutrino    (hypoKey_)->p4(),
			    tops.singleLepton(hypoKey_)->charge(),
			    weight);
    }
    // ---
    //    check quark assignment
    // ---
    // no parton jet parton match exists
    if(tops.isHypoValid("kKinFit")&&(!tops.isHypoValid("kGenMatch"))){
      assignment= 9.;
    }
    // if jet parton match exists:
    if(tops.isHypoValid("kKinFit")&& tops.isHypoValid("kGenMatch")){
      // indices for all quarks from Kinfit Hypothesis and genmatch
      int lepBIndex         = tops.jetLeptonCombination("kKinFit"  )[TtSemiLepEvtPartons::LepB     ];
      int hadBIndex         = tops.jetLeptonCombination("kKinFit"  )[TtSemiLepEvtPartons::HadB     ];
      int lightQIndex       = tops.jetLeptonCombination("kKinFit"  )[TtSemiLepEvtPartons::LightQ   ];
      int lightQBarIndex    = tops.jetLeptonCombination("kKinFit"  )[TtSemiLepEvtPartons::LightQBar];
      int lepBIndexGen      = tops.jetLeptonCombination("kGenMatch")[TtSemiLepEvtPartons::LepB     ];
      int hadBIndexGen      = tops.jetLeptonCombination("kGenMatch")[TtSemiLepEvtPartons::HadB     ];
      int lightQIndexGen    = tops.jetLeptonCombination("kGenMatch")[TtSemiLepEvtPartons::LightQ   ];
      int lightQBarIndexGen = tops.jetLeptonCombination("kGenMatch")[TtSemiLepEvtPartons::LightQBar];
      // calculate permutation
      // 0: nothing wrong
      if((lepBIndex==lepBIndexGen)&&(hadBIndex==hadBIndexGen)&&
	 (((lightQIndex==lightQIndexGen   )&&(lightQBarIndex==lightQBarIndexGen))||
	  ((lightQIndex==lightQBarIndexGen)&&(lightQBarIndex==lightQIndexGen   )))) assignment=0;
      else{
	// 1: b quarks switched
	if((lepBIndex==hadBIndexGen)&&(hadBIndex==lepBIndexGen)) assignment=1;
	// 2: leptonic b and light quark switched
	if(((lepBIndex==lightQIndexGen)||(lepBIndex==lightQBarIndexGen))&&
	   (((lightQIndex==lepBIndexGen)||(lightQBarIndex==lepBIndexGen)))) assignment=2;
	// 3: hadronic b and light quark switched/
	if(((hadBIndex==lightQIndexGen)||(hadBIndex==lightQBarIndexGen))&&
	   (((lightQIndex==hadBIndexGen)||(lightQBarIndex==hadBIndexGen)))) assignment=3;
	// 4: light quark->leptonic b & leptonic b->hadronic b & hadronic b-> light quark
	if(((lepBIndex==lightQIndexGen)||(lepBIndex==lightQBarIndexGen))&&(hadBIndex==lepBIndexGen)&&
	   ((lightQIndex==hadBIndexGen)||(lightQBarIndex==hadBIndexGen))) assignment=4;
	// 5: light quark->hadronic b & hadronic b->leptonic b & leptonic b-> light quark
	if(((hadBIndex==lightQIndexGen)||(hadBIndex==lightQBarIndexGen))&&(lepBIndex==hadBIndexGen)&&
	   ((lightQIndex==lepBIndexGen)||(lightQBarIndex==lepBIndexGen))) assignment=5;
	// 6: hadronic/leptonic b-> light quarks &  light quarks->hadronic/leptonic b
	if(((hadBIndex     ==lightQIndexGen)||(hadBIndex     ==lightQBarIndexGen))&&
	   ((lepBIndex     ==lightQIndexGen)||(lepBIndex     ==lightQBarIndexGen))&&
	   ((lightQIndex   ==lepBIndexGen  )||(lightQIndex   ==hadBIndexGen     ))&&
	   ((lightQBarIndex==lepBIndexGen  )||(lightQBarIndex==hadBIndexGen     ))) assignment=6;
	// make sure that no relevant jet is missing
	std::vector<int> genJets_, recoJets_;
	// list of genJets
	genJets_ .push_back(lepBIndexGen     );
	genJets_ .push_back(hadBIndexGen     );
	genJets_ .push_back(lightQIndexGen   );
	genJets_ .push_back(lightQBarIndexGen);
	std::sort( genJets_.begin(), genJets_.end());
	// list of recoJets
	recoJets_.push_back(lepBIndex);
	recoJets_.push_back(hadBIndex);
	recoJets_.push_back(lightQIndex);
	recoJets_.push_back(lightQBarIndex);
	std::sort( recoJets_.begin(), recoJets_.end());
	// compare recoJets and genJets
	for(unsigned int i=0; i<recoJets_.size(); ++i){
	  if( recoJets_[i]!=genJets_[i] ){ 
	    if(maxNJets<4){
	      std::cout << "ERROR: number of conidered jets can not be smaller than 4" << std::endl;
	      exit(1);
	    }


	    // 7: jet is missing
	    if( genJets_.back()>maxNJets-1 ) assignment=7;
	    // 8: wrong jet chosen (only valid if kinFitTtSemiLepEventHypothesis.maxNJets>4)
	    // e.g. took the wrong 4 out of 5 jets 
	    else{
	      if(assignment<0) assignment=8;
	    }
	    break;
	  }
	}
      }
    }
    // save lepton charge
    charge= ((reco::LeafCandidate*)(tops.singleLepton(hypoKey_)))->charge();
  }
  // fill trees and histos
  fillValue( "prob"       , prob      , weight);
  fillValue( "chi2"       , chi2      , weight);
  fillValue( "delChi2"    , delChi2   , weight);
  fillValue( "qAssignment", assignment, weight);
  fillValue( "lepCharge"  , charge    , weight);
  fillValue( "lepTopIsTopPlus"       , switchLepAndHadTop ,weight);
  fillValue( "topPtTtbarSysPartonTruth", gentopPtTtbarSystem, weight );
  fillValue( "ttbarPtPartonTruth"    , genTtbarPt      , weight );
  fillValue( "ttbarYPartonTruth"     , genTtbarRapidity, weight );
  fillValue( "ttbarPhiPartonTruth"   , genTtbarPhi     , weight );
  fillValue( "ttbarMassPartonTruth"  , genTtbarMass    , weight );
  fillValue( "ttbarDelPhiPartonTruth", genDeltaPhi     , weight );
  fillValue( "ttbarDelYPartonTruth"  , genDeltaRapidity, weight );
  fillValue( "ttbarSumYPartonTruth"  , sumRapidity     , weight );
  fillValue( "ttbarHTPartonTruth"    , HTgen           , weight );
  fillValue( "ttbarPhiStarPartonTruth", genTtbarPhiStar , weight );
  fillValue( "topPtHadPartonTruth"   , hadTopGenPt      , weight );
  fillValue( "topYHadPartonTruth"    , hadTopGenRapidity, weight );
  fillValue( "topPhiHadPartonTruth"  , hadTopGenPhi     , weight );
  fillValue( "hadTopMassPartonTruth" , hadTopGenMass    , weight );
  fillValue( "topPtLepPartonTruth"   , lepTopGenPt      , weight );
  fillValue( "topYLepPartonTruth"    , lepTopGenRapidity, weight );  
  fillValue( "topPhiLepPartonTruth"  , lepTopGenPhi     , weight );
  fillValue( "lepTopMassPartonTruth" , lepTopGenMass    , weight );
  fillValue( "lepPtPartonTruth"      , genLepPt , weight );      
  fillValue( "lepEtaPartonTruth"     , genLepEta, weight );  
  fillValue( "neutrinoPtPartonTruth"   , genNuPt    , weight );
  fillValue( "neutrinoEtaPartonTruth"  , genNuEta   , weight );
  fillValue( "lightQPtPartonTruth"     , genQPt     , weight );
  fillValue( "lightQbarPtPartonTruth"  , genQbarPt  , weight );
  fillValue( "bqPtLepPartonTruth"      , genLepBPt  , weight );
  fillValue( "bqPtHadPartonTruth"      , genHadBPt  , weight );  
  fillValue( "lightQEtaPartonTruth"    , genQEta    , weight );
  fillValue( "lightQbarEtaPartonTruth" , genQbarEta , weight );
  fillValue( "bqEtaLepPartonTruth"     , genLepBEta , weight );
  fillValue( "bqEtaHadPartonTruth"     , genHadBEta , weight );
  fillValue( "bbbarPtPartonTruth"     , genBbbarPt  , weight );
  fillValue( "bbbarYPartonTruth"      , genBbbarY   , weight );
  fillValue( "bbbarMassPartonTruth"   , genBbbarMass, weight );
  fillValue( "lbMassPartonTruth"      , genlbMass   , weight );
  fillValue( "topYPlusPartonTruth"    , genTopPlusY   , weight );
  fillValue( "topYMinusPartonTruth"   , genTopMinusY  , weight );
  fillValue( "topEtaPlusPartonTruth"  , genTopPlusEta , weight );
  fillValue( "topEtaMinusPartonTruth" , genTopMinusEta, weight );
  fillValue( "topPtPlusPartonTruth"   , genTopPlusPt  , weight );
  fillValue( "topPtMinusPartonTruth"  , genTopMinusPt , weight );
  fillValue( "topPtLeadPartonTruth"   , genptLead     , weight );
  fillValue( "topYLeadPartonTruth"    , genyLead      , weight );
  fillValue( "topPtSubLeadPartonTruth", genptSubLead  , weight );
  fillValue( "topYSubLeadPartonTruth" , genySubLead   , weight );
  //fillValue( "xBLepPartonTruth"       , lepxBGen      , weight );
  //fillValue( "xBHadPartonTruth"       , hadxBGen      , weight );

  // fill the tree, if any variable should be put in
  if(treeVars_.size()) tree->Fill();
}

/// histogram filling for kinematic 1D histos
void
TopKinematics::fill(const reco::Candidate* leptonicTop, const reco::Candidate* hadronicTop, 
		    const reco::Candidate* topPlus    , const reco::Candidate* topMinus   , 
		    const reco::Candidate* leptonicW  , const reco::Candidate* hadronicW, 
		    const reco::Candidate* leptonicB  , const reco::Candidate* hadronicB,
		    double HT, const double& charge, const double& weight)
{
  /** 
      calculate boosted Lorentz vectors
  **/

  // define lorentz vectors for ttbar system
  reco::Particle::LorentzVector ttBar = leptonicTop->p4() + hadronicTop->p4();
  // define lorentz vectors for bbbar system
  //reco::Particle::LorentzVector bbBar = leptonicB->p4() + hadronicB->p4();
  // create boost into reconstructed top quarks/ttbar system
  ROOT::Math::Boost CoMBoostHadTop(hadronicTop->p4().BoostToCM());
  ROOT::Math::Boost CoMBoostLepTop(leptonicTop->p4().BoostToCM());
  ROOT::Math::Boost CoMBoostTtbar (ttBar            .BoostToCM());
  // get lorentz vectors of W and Top without boost
  reco::Particle::LorentzVector HadronicDecayWBoosted   = hadronicW  ->p4();
  reco::Particle::LorentzVector LeptonicDecayWBoosted   = leptonicW  ->p4();
  reco::Particle::LorentzVector HadronicDecayTopBoosted = hadronicTop->p4();
  reco::Particle::LorentzVector LeptonicDecayTopBoosted = leptonicTop->p4();
  // recalculate Lorent vectors using boosts defined above
  // boost W into top rest frame and
  // and top into ttbar rest frame
  HadronicDecayWBoosted   = CoMBoostHadTop(HadronicDecayWBoosted  );
  LeptonicDecayWBoosted   = CoMBoostLepTop(LeptonicDecayWBoosted  );
  HadronicDecayTopBoosted = CoMBoostTtbar (HadronicDecayTopBoosted);
  LeptonicDecayTopBoosted = CoMBoostTtbar (LeptonicDecayTopBoosted);
      
  /** 
      Fill Top Variables for Cross Section Measurement
  **/

  // ---
  //    top variables (as no branch is created, they exist only as plot)
  // ---
  // fill top pt in ttbar restframe 
  fillValue( "topPtTtbarSys", HadronicDecayTopBoosted.pt() , weight );
  // fill top pt for leptonicTop candidate in combined histogram
  fillValue( "topPt", leptonicTop->p4().pt()     , weight );
  // fill top pt for hadronicTop candidate in combined histogram
  fillValue( "topPt", hadronicTop->p4().pt()     , weight );
  // fill top y for leptonicTop candidate in combined histogram
  fillValue( "topY", leptonicTop->p4().Rapidity(), weight );
  // fill top y for hadronicTop candidate in combined histogram
  fillValue( "topY", hadronicTop->p4().Rapidity(), weight );
  // fill top phi for leptonicTop candidate in combined histogram
  fillValue( "topPhi", leptonicTop->p4().phi()   , weight );
  // fill top phi for hadronicTop candidate in combined histogram
  fillValue( "topPhi", hadronicTop->p4().phi()   , weight );
  // fill leptonic Top mass in combined histogram
  fillValue( "topMass", leptonicTop->mass()      , weight );
  // fill hadronic Top mass in combined histogram
  fillValue( "topMass", hadronicTop->mass()      , weight );

  // ---
  //    top variables for leptonic top
  // ---
  // fill top pt for leptonicTop candidate in separate histogram
  fillValue( "topPtLep", leptonicTop->p4().pt()     , weight );
  // fill top y for leptonicTop candidate in separate histogram
  fillValue( "topYLep", leptonicTop->p4().Rapidity(), weight );
  // fill top phi for leptonicTop candidate in separate histogram
  fillValue( "topPhiLep", leptonicTop->p4().phi()   , weight );
  // fill leptonic Top mass in separate histogram
  fillValue( "lepTopMass", leptonicTop->mass()      , weight );


  // ---
  //    top variables for hadronic top
  // ---
  // fill top pt for hadronicTop candidate in separate histogram
  fillValue( "topPtHad", hadronicTop->p4().pt(), weight );
  // fill top y for hadronicTop candidate in separate histogram
  fillValue( "topYHad", hadronicTop->p4().Rapidity(), weight );
  // fill top phi for hadronicTop candidate in separate histogram
  fillValue( "topPhiHad", hadronicTop->p4().phi(), weight );
  // fill hadronic Top mass in separate histogram
  fillValue( "hadTopMass", hadronicTop->mass(), weight );

  // ---
  //    ttbar variables
  // ---
  // fill ttbar pt
  fillValue( "ttbarPt", ttBar.pt(), weight );
  // fill ttbar y
  fillValue( "ttbarY", ttBar.Rapidity(), weight );
  // fill ttbar phi
  fillValue( "ttbarPhi", ttBar.phi(), weight );
  // fill ttbar invariant mass
  fillValue( "ttbarMass", ttBar.mass(), weight );
  // fill deltaPhi between both top quarks 
  bool switchLepAndHadTop= charge<0 ? true : false;
  if(switchLepAndHadTop)  fillValue( "ttbarDelPhi", std::abs(deltaPhi(hadronicTop->phi(), leptonicTop->phi())), weight );
  else                    fillValue( "ttbarDelPhi", std::abs(deltaPhi(leptonicTop->phi(), hadronicTop->phi())), weight );
  // fill deltaY between both top quarks 
  if(switchLepAndHadTop) fillValue( "ttbarDelY", hadronicTop->rapidity()-leptonicTop->rapidity(), weight );
  else                   fillValue( "ttbarDelY", leptonicTop->rapidity()-hadronicTop->rapidity(), weight );
  // fill sum of y of both top quarks 
  fillValue( "ttbarSumY", leptonicTop->rapidity()+hadronicTop->rapidity(), weight );
  // fill HT of the 4 jets assigned to the ttbar decay
  fillValue( "ttbarHT", HT, weight );
  // fill ttbar phi star
  double ttbarPhiStar  = TMath::Tan(0.5*(pi-std::abs(deltaPhi(leptonicTop->phi(),hadronicTop->phi()))))/TMath::CosH(0.5*(leptonicTop->eta()-hadronicTop->eta()));
  fillValue( "ttbarPhiStar", ttbarPhiStar, weight );

  // ---
  //    asymmetry variables
  // ---
  fillValue( "topYPlus"   , topPlus->p4().Rapidity() , weight );
  fillValue( "topEtaPlus" , topPlus->p4().eta()      , weight );
  fillValue( "topYMinus"  , topMinus->p4().Rapidity(), weight );
  fillValue( "topEtaMinus", topMinus->p4().eta()     , weight );
  fillValue( "topPtPlus"  , topPlus->p4().pt()       , weight );
  fillValue( "topPtMinus" , topMinus->p4().pt()      , weight );

  // find leading top
  double ptLead   =leptonicTop->p4().pt();
  double ptSubLead=hadronicTop->p4().pt();
  double yLead   =leptonicTop->p4().Rapidity();
  double ySubLead=hadronicTop->p4().Rapidity();
  if(ptLead<ptSubLead){// swop them
    double temp=ptLead;
    double temp2=yLead;
    ptLead=ptSubLead;
    yLead=ySubLead;
    ptSubLead=temp;
    ySubLead=temp2;
  }
  // fill top pt for leading Top candidate
  fillValue( "topPtLead"   , ptLead   , weight );
  // fill top pt for subleading Top candidate
  fillValue( "topPtSubLead", ptSubLead, weight );
  // fill top y for leading Top candidate
  fillValue( "topYLead"    , yLead    , weight );
  // fill top y for subleading Top candidate
  fillValue( "topYSubLead" , ySubLead , weight );

  // fill  bbbar variables
  //fillValue( "bbbarPt"  , bbBar.pt()      , weight );
  //fillValue( "bbbarY"   , bbBar.Rapidity(), weight );
  //fillValue( "bbbarMass", bbBar.mass()    , weight );

  // fill b-top hadronization variable
  // (1/(1-mW^1/mt^2)) * 2pBpT/mt^2
  //double lepxBRec =(1/(1-((leptonicW->p4()).mass2())/((leptonicTop->p4()).mass2()))) * 2*((leptonicB->p4()).Dot(leptonicTop->p4()))/((leptonicTop->p4()).mass2());
  //double hadxBRec =(1/(1-((hadronicW->p4()).mass2())/((hadronicTop->p4()).mass2()))) * 2*((hadronicB->p4()).Dot(hadronicTop->p4()))/((hadronicTop->p4()).mass2());
  //fillValue( "xB"   , lepxBRec, weight );
  //fillValue( "xB"   , hadxBRec, weight );
  //fillValue( "xBLep", lepxBRec, weight );
  //fillValue( "xBHad", hadxBRec, weight );
}

/// histogram filling for kinematic 1D histos using only events with Ngen&&Nreco 
/// in the same bin wrt the binning defined above
void 
TopKinematics::fillGenRec(const reco::Candidate* leptonicTopRec, const reco::Candidate* leptonicTopGen, 
			  const reco::Candidate* hadronicTopRec, const reco::Candidate* hadronicTopGen, 
			  const reco::Candidate* leptonicWRec  , const reco::Candidate* leptonicWGen, 
			  const reco::Candidate* hadronicWRec  , const reco::Candidate* hadronicWGen,
			  double HTrec, double HTgen, const double& weight)
{
  /** 
      calculate boosted Lorentz vectors
  **/
  // define lorentz vectors for ttbar system
  reco::Particle::LorentzVector genTtBar = leptonicTopGen->p4() + hadronicTopGen->p4();
  reco::Particle::LorentzVector recTtBar = leptonicTopRec->p4() + hadronicTopRec->p4();
  // create boost into reconstructed/generated top quarks/ttbar system
  ROOT::Math::Boost CoMBoostRecHadTop(hadronicTopRec->p4().BoostToCM());
  ROOT::Math::Boost CoMBoostGenHadTop(hadronicTopGen->p4().BoostToCM());
  ROOT::Math::Boost CoMBoostRecLepTop(leptonicTopRec->p4().BoostToCM());
  ROOT::Math::Boost CoMBoostGenLepTop(leptonicTopGen->p4().BoostToCM());
  ROOT::Math::Boost CoMBoostRecTtbar (recTtBar            .BoostToCM());
  ROOT::Math::Boost CoMBoostGenTtbar (genTtBar            .BoostToCM());
  // get lorentz vectors of W and Top without boost
  reco::Particle::LorentzVector recHadronicDecayWBoosted   = hadronicWRec  ->p4();
  reco::Particle::LorentzVector genHadronicDecayWBoosted   = hadronicWGen  ->p4();
  reco::Particle::LorentzVector recLeptonicDecayWBoosted   = leptonicWRec  ->p4();
  reco::Particle::LorentzVector genLeptonicDecayWBoosted   = leptonicWGen  ->p4();
  reco::Particle::LorentzVector recHadronicDecayTopBoosted = hadronicTopRec->p4();
  reco::Particle::LorentzVector genHadronicDecayTopBoosted = hadronicTopGen->p4();
  reco::Particle::LorentzVector recLeptonicDecayTopBoosted = leptonicTopRec->p4();
  reco::Particle::LorentzVector genLeptonicDecayTopBoosted = leptonicTopGen->p4();
  // recalculate Lorent vectors using boosts defined above
  // boost W into top rest frame and top into ttbar rest frame
  recHadronicDecayWBoosted   = CoMBoostRecHadTop(recHadronicDecayWBoosted  );
  genHadronicDecayWBoosted   = CoMBoostGenHadTop(genHadronicDecayWBoosted  );
  recLeptonicDecayWBoosted   = CoMBoostRecLepTop(recLeptonicDecayWBoosted  );
  genLeptonicDecayWBoosted   = CoMBoostGenLepTop(genLeptonicDecayWBoosted  );
  recHadronicDecayTopBoosted = CoMBoostRecTtbar (recHadronicDecayTopBoosted);
  genHadronicDecayTopBoosted = CoMBoostGenTtbar (genHadronicDecayTopBoosted);
  recLeptonicDecayTopBoosted = CoMBoostRecTtbar (recLeptonicDecayTopBoosted);
  genLeptonicDecayTopBoosted = CoMBoostGenTtbar (genLeptonicDecayTopBoosted);

  /** 
      Fill 1D histos for events gen&&rec in the same bin, using function match()
  **/
  //std::cout << "reclepTop : " << recLeptonicDecayTopBoosted.pt() << ", rechadTop : " << recHadronicDecayTopBoosted.pt() << std::endl;
  //std::cout << "genlepTop : " << genLeptonicDecayTopBoosted.pt() << ", genhadTop : " << genHadronicDecayTopBoosted.pt() << std::endl;
  match( "topPtTtbarSys", recHadronicDecayTopBoosted.pt() , genHadronicDecayTopBoosted.pt(), weight );
  // fill top pt for topA candidate in combined histogram
  match( "topPt"       , leptonicTopRec->p4().pt ()     , leptonicTopGen->p4().pt ()      , weight );
  // fill top pt for topA candidate in separate histogram
  match( "topPtLep"    , leptonicTopRec->p4().pt ()      , leptonicTopGen->p4().pt ()     , weight );
  // fill top y for topA candidate in combined histogram
  match( "topY"        , leptonicTopRec->p4().Rapidity(), leptonicTopGen->p4().Rapidity() , weight );
  // fill top y for topA candidate in separate histogram
  match( "topYLep"     , leptonicTopRec->p4().Rapidity() , leptonicTopGen->p4().Rapidity(), weight );
  // fill top phi for topA candidate in combined histogram
  match( "topPhi"      , leptonicTopRec->p4().phi()     , leptonicTopGen->p4().phi()      , weight );
  // fill top phi for topA candidate in separate histogram
  match( "topPhiLep"   , leptonicTopRec->p4().phi()      , leptonicTopGen->p4().phi()     , weight );
  // fill leptonic Top mass for topB candidate in seperate histogram
  match( "lepTopMass"  , leptonicTopRec->mass()         , leptonicTopGen->mass()          , weight );
  // fill leptonic Top mass for topB candidate in combined histogram
  match( "topMass"     , leptonicTopRec->mass()         , leptonicTopGen->mass()          , weight );
  // fill top pt for topB candidate in combined histogram
  match( "topPt"       , hadronicTopRec->p4().pt ()      , hadronicTopGen->p4().pt ()      , weight );
  // fill top pt for topB candidate in separate histogram
  match( "topPtHad"    , hadronicTopRec->p4().pt ()      , hadronicTopGen->p4().pt ()      , weight );
  // fill top y for topB candidate in combined histogram
  match( "topY"        , hadronicTopRec->p4().Rapidity() , hadronicTopGen->p4().Rapidity() , weight );
  // fill top y for topB candidate in separate histogram
  match( "topYHad"     , hadronicTopRec->p4().Rapidity() , hadronicTopGen->p4().Rapidity() , weight );
  // fill top phi for topB candidate in combined histogram
  match( "topPhi"      , hadronicTopRec->p4().phi()      , hadronicTopGen->p4().phi()      , weight );
  // fill top phi for topB candidate in separate histogram
  match( "topPhiHad"   , hadronicTopRec->p4().phi()      , hadronicTopGen->p4().phi()      , weight );
  // fill hadronic Top mass for topB candidate in seperate histogram
  match( "hadTopMass"  , hadronicTopRec->mass()          , hadronicTopGen->mass()          , weight );
  // fill hadronic Top mass for topB candidate in combined histogram
  match( "topMass"     , hadronicTopRec->mass()          , hadronicTopGen->mass()          , weight );
  // fill ttbar pt
  match( "ttbarPt"     , recTtBar.pt  ()          , genTtBar.pt  ()          , weight );
  // fill ttbar y
  match( "ttbarY"      , recTtBar.Rapidity ()     , genTtBar.Rapidity ()     , weight );
  // fill ttbar phi
  match( "ttbarPhi"    , recTtBar.phi ()          , genTtBar.phi ()          , weight );
  // fill ttbar invariant mass
  match( "ttbarMass"   , recTtBar.mass()          , genTtBar.mass()          , weight );

  // fill deltaPhi between both top quarks 
  match( "ttbarDelPhi" , std::abs(deltaPhi(leptonicTopRec->phi(), hadronicTopRec->phi())), 
	                 std::abs(deltaPhi(leptonicTopGen->phi(), hadronicTopGen->phi())), weight );
  // fill deltaY between both top quarks 
  match( "ttbarDelY"   , leptonicTopRec->rapidity()-hadronicTopRec->rapidity() , 
                         leptonicTopGen->rapidity()-hadronicTopGen->rapidity() , weight );

  // fill sum of y of both top quarks 
  match( "ttbarSumY"   , leptonicTopRec->rapidity()+hadronicTopRec->rapidity() , 
	                 leptonicTopGen->rapidity()+hadronicTopGen->rapidity() , weight );
  // fill HT of the 4 jets assigned to the ttbar decay
  match( "ttbarHT"     , HTrec , HTgen , weight );
}

/// helper function to provid gen&&rec-histo for the determination of stability and purity
void 
TopKinematics::match(const std::string& histo, const double& recValue, const double& genValue, const double& weight)
{
  TH1* hist = hists_.find(histo)->second;
  // loop the bins of histogram 'hist'; note that the first bin
  // starts with '1' and not with '0'
  for(int bin=1; bin<=hist->GetNbinsX(); ++bin){
    // determine lower and upper edge of the corresponding bin
    double lowerEdge = hist->GetBinLowEdge(bin);
    double upperEdge = hist->GetBinLowEdge(bin)+hist->GetBinWidth(bin);
    // only fill the histogram when the genValue & recValue fall into the same bin ('<=' for upper edge)
    if( (lowerEdge<genValue && genValue<=upperEdge) && (lowerEdge<recValue && recValue<=upperEdge) ){
      fillValue( histo, recValue, weight );
      break;
    }
  }

}

///  helper functions to fill 1D angle histos
void
TopKinematics::fillAngles(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& hadB    ,
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& q       ,
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& qbar    ,
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& lepB    ,
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& muon    ,
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& neutrino,
			  const double& weight)
{ 
  // calculate angles from the 4 momentum vectors
  TopAngles angles = TopAngles(hadB, q   , qbar    ,  // branch1: hadronic
			       lepB, muon, neutrino); // branch2: leptonic
  // angle between t candidate (detector rest frame)
  fillValue("ttbarAngle"       , angles.ttDetFrame()  , weight);
  // angle between b candidates (detector rest frame)
  fillValue("bbbarAngle"       , angles.bbDetFrame()  , weight);
  // angle between b candidates (ttbar rest frame)
  fillValue("bbbarAngleTtRF"   , angles.bbTtbarFrame(), weight);
  // angle between W bosons (ttbar rest frame)
  fillValue("WWAngle"          , angles.WWTtbarFrame(), weight);
  // angle between the hadronically decaying top candidate and the corresponding b (ttbar rest frame)
  fillValue("topBAngleHad"     , angles.tBBranch1TtbarFrame(), weight);
  // angle between the leptonically decaying top candidate and the corresponding b (ttbar rest frame)
  fillValue("topBAngleLep"     , angles.tBBranch2TtbarFrame(), weight);
  // angle between the leptonically decaying W candidate and the corresponding b (ttbar rest frame)
  fillValue("bWAngleLep"       , angles.bWBranch2TtbarFrame(), weight);
  // angle between the hadronically decaying W candidate and the corresponding b (ttbar rest frame)
  fillValue("bWAngleHad"       , angles.bWBranch1TtbarFrame(), weight);
  // angle between the leptonically decaying top candidate and the corresponding W 
  // (Top in ttbar rest frame and W in top rest frame)
  fillValue("topWAngleLep"     , angles.tWBranch2TopInTtbarFrameWInTopFrame(), weight);
  // angle between the hadronically decaying top candidate and the corresponding W 
  // (Top in ttbar rest frame and W in top rest frame)   
  fillValue("topWAngleHad"     , angles.tWBranch1TopInTtbarFrameWInTopFrame(), weight);
  // angle between light quarks (ttbar rest frame)
  fillValue("qqbarAngle"       , angles.qQbarTopFrame()      , weight);
  // angle between light quarks and leptonic b (ttbar rest frame)
  fillValue("qBlepAngle"       , angles.blepQTtbarFrame()    , weight);
  fillValue("qBlepAngle"       , angles.blepQbarTtbarFrame() , weight);
  // angle between light quarks and hadronic b (corresponding top candidate rest frame)
  fillValue("qBhadAngle"       , angles.bhadQTopFrame()      , weight);
  fillValue("qBhadAngle"       , angles.bhadQbarTopFrame()   , weight);
  // angle between lepton and leptonic b (corresponding top candidate rest frame)
  fillValue("lepBlepAngle"     , angles.lepBlepTopFrame()    , weight);
  // angle between lepton and leptonic b (corresponding top candidate rest frame)
  fillValue("lepBlepAngleTtRF" , angles.lepBlepTtbarFrame()  , weight);
  // angle between lepton and hadronic b (ttbar rest frame)
  fillValue("lepBhadAngle"     , angles.lepBhadTtbarFrame()  , weight);
  // angle between lepton and light quarks (ttbar rest frame)
  fillValue("lepQAngle"        , angles.lepQTtbarFrame()     , weight);
  fillValue("lepQAngle"        , angles.lepQbarTtbarFrame()  , weight);
  // angle between muon and neutrino (in corresponding top rest frame)
  fillValue("MuonNeutrinoAngle", angles.lepNuTopFrame()      , weight);
  // angle between leptonic b and neutrino (in corresponding top rest frame)
  fillValue("lepBNeutrinoAngle", angles.nuBlepTopFrame()     , weight);
  // angle between hadronic b and neutrino (in ttbar rest frame)
  fillValue("hadBNeutrinoAngle", angles.nuBhadTtbarFrame()   , weight);
  // angle between light quarks and neutrino (in ttbar rest frame)
  fillValue("qNeutrinoAngle"   , angles.nuQTtbarFrame()      , weight);
  fillValue("qNeutrinoAngle"   , angles.nuQBarTtbarFrame()   , weight);
  // angle between lepton (in W rest frame) and corresponding W (in detector rest frame)
  fillValue("lepWDir"          , angles.lepWlepLepInWFrameWInDetFrame(), weight);
  // angle between neutrino (in W rest frame) and corresponding W (in detector rest frame)
  fillValue("nuWDir"           , angles.nuWlepNuInWFrameWInDetFrame()  , weight);
  // angle between light quarks (in W rest frame) and corresponding W (in detector rest frame)
  fillValue("qWDir"            , angles.qWhadQInWFrameWInDetFrame()    , weight);
  fillValue("qWDir"            , angles.qbarWhadQInWFrameWInDetFrame() , weight);
}


///  helper functions to fill event shape variables
void
TopKinematics::fillEventShapes(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& hadB    ,
			       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& q       ,
			       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& qbar    ,
			       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& lepB    ,
			       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& muon    ,
			       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& neutrino,
			       const double& weight, const bool useMu, const bool useNu)
{ 
  // ---
  //    convert 4 momentum vector into 3 momentum vector
  // ---
  std::vector<math::XYZVector> objects;
  math::XYZVector hadBjet(hadB.px(), hadB.py(), hadB.pz());
  objects.push_back(hadBjet);
  math::XYZVector lepBjet(lepB.px(), lepB.py(), lepB.pz());
  objects.push_back(lepBjet);
  math::XYZVector lightquarkjet1(q.px(), q.py(), q.pz());
  objects.push_back(lightquarkjet1);
  math::XYZVector lightquarkjet2(qbar.px(), qbar.py(), qbar.pz());
  objects.push_back(lightquarkjet2);
  if(useMu){
    math::XYZVector mu(muon.px(), muon.py(), muon.pz());
    objects.push_back(mu);
  }
  if(useNu){
    math::XYZVector nu(neutrino.px(), neutrino.py(), neutrino.pz());
    objects.push_back(nu);
  }
  // ---
  //    Calculate Event Shape Variables
  // ---
  EventShapeVariables eventshape(objects);
  double aplanarity_  = eventshape.aplanarity();
  double sphericity_  = eventshape.sphericity();
  double circularity_ = eventshape.circularity();
  double isotropy_    = eventshape.isotropy();
  double C_           = eventshape.C();
  double D_           = eventshape.D();

  // ---
  //    fill plots
  // ---
  // notation:
  // 0<=q1<=q2<=q3 are the eigenvalues of the momemtum tensor:
  // T_ab = sum{p_j[a]*p_j[b]}/sum{p_j**2} normalized to 1.
  // aplanarity of event: 
  // 1.5*q1 
  // Return values are 0.5 for spherical and 0 for plane and linear events
  fillValue("aplanarity" , aplanarity_ , weight);
  // sphericity of event:
  // 1.5*(q1+q2) 
  // Return values are 1 for spherical, 3/4 for plane and 0 for linear events
  fillValue("sphericity" , sphericity_ , weight);
  // C (three jet structure) of event:
  // 3.*(q1*q2+q1*q3+q2*q3)
  // Return value is between 0 and 1, C vanishes for a "perfect" 2-jet event
  fillValue("C"          , C_          , weight);
  // D (four jet structure) of event:
  // 27.*(q1*q2*q3)
  // Return value is between 0 and 1, D vanishes for a planar event
  fillValue("D"          , D_          , weight);
  // circularity of event
  // the return value is 1 for spherical and 0 linear events in r-phi
  fillValue("circularity", circularity_, weight);
  // isotropy of event:
  // the return value is 1 for spherical events and 0 for events linear in r-phi
  fillValue("isotropy"   , isotropy_   , weight);
}

///  helper functions to fill final state objects
void
TopKinematics::fillFinalStateObjects(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& hadB    ,
				     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& q       ,
				     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& qbar    ,
				     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& lepB    ,
				     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& lepton  ,
				     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& neutrino,
				     const double& charge,
				     const double& weight)
{
  // ---
  //    fill ttbar final state object distributions
  // ---  

  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > bbbar=lepB+hadB;

  // a) fill trees+histos for lepton and neutrino
  // leptonPt
  fillValue("lepPt"      , lepton.Pt()    , weight);
  // leptonEta				  
  fillValue("lepEta"     , lepton.Eta()	  , weight);
  // neutrinoPt				  
  fillValue("neutrinoPt" , neutrino.Pt()  , weight);
  // neutrinoEta			  
  fillValue("neutrinoEta", neutrino.Eta() , weight);
  // b) fill histos for quarks: sepearte only b/light quarks
  // light-quarks Pt		  
  fillValue("lightqPt"   , q.Pt()	  , weight);
  fillValue("lightqPt"   , qbar.Pt()	  , weight);
  // light-quarks Eta			  
  fillValue("lightqEta"  , q.Eta()	  , weight);
  fillValue("lightqEta"  , qbar.Eta()	  , weight);
  // b-quarks Pt			  
  fillValue("bqPt"       , lepB.Pt()	  , weight);
  fillValue("bqPt"       , hadB.Pt()	  , weight);
  // b-quarks Eta			  
  fillValue("bqEta"      , lepB.Eta()	  , weight);
  fillValue("bqEta"      , hadB.Eta()     , weight);
      							
  fillValue("bbbarPt"    , bbbar.Pt()      , weight);
  fillValue("bbbarY"     , bbbar.Rapidity(), weight);
  fillValue("bbbarMass"  , bbbar.mass()    , weight);
  // fill trees for (light/b)x(quark/antiquark)
  fillValue("lightQPt"    , q.Pt()    , weight);
  fillValue("lightQbarPt" , qbar.Pt() , weight);
  fillValue("bqPtLep"     , lepB.Pt() , weight);
  fillValue("bqPtHad"     , hadB.Pt() , weight);  
  fillValue("lightQEta"   , q.Eta()   , weight);
  fillValue("lightQbarEta", qbar.Eta(), weight);
  fillValue("bqEtaLep"    , lepB.Eta(), weight);
  fillValue("bqEtaHad"    , hadB.Eta(), weight);

  // fill asymmetry variables
  if(charge>0){
    fillValue("lepPtPlus" , lepton.pt()      , weight );
    fillValue("lepEtaPlus", lepton.eta()     , weight );
    fillValue("lepYPlus"  , lepton.Rapidity(), weight );
  }				 
  else{			 
    fillValue("lepPtMinus" , lepton.pt()      , weight );
    fillValue("lepEtaMinus", lepton.eta()     , weight );
    fillValue("lepYMinus"  , lepton.Rapidity(), weight );
  }
  
  // find leading jet
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > leadq=hadB;
  if(lepB.Pt()>leadq.Pt()) leadq=lepB;
  if(q   .Pt()>leadq.Pt()) leadq=q   ;
  if(qbar.Pt()>leadq.Pt()) leadq=qbar;
  // leading quark Pt
  fillValue("leadqPt" , leadq.Pt()  , weight);
  // leading quark Eta
  fillValue("leadqEta", leadq.Eta() , weight);
  // lepton b-quark system
  fillValue( "lbMass", (lepton+lepB).mass()   , weight );
}

// get the decay mode
double 
TopKinematics::getDecayChannel(const TtGenEvent& tops){
  double eventType=-4;
  // identify event by the following number coding:
  // -4 : undefined
  // -3 : non ttbar MC
  // -2 : dileptonic unknown    ttbar MC
  // -1 : semileptonic unknown  ttbar MC
  // 01 : semileptonic electron ttbar MC
  // 02 : semileptonic muon     ttbar MC
  // 03 : semileptonic tau      ttbar MC
  // 11 : dileptonic e   e      ttbar MC
  // 12 : dileptonic e   mu     ttbar MC
  // 13 : dileptonic e   tau    ttbar MC
  // 14 : dileptonic mu  mu     ttbar MC
  // 15 : dileptonic mu  tau    ttbar MC
  // 16 : dileptonic tau tau    ttbar MC
  // 20 : fully hadronic        ttbar MC
  if( !tops.isTtBar()        ) eventType = -3;
  else if(  tops.isFullHadronic() ) eventType = 20;
  else if(  tops.isSemiLeptonic() ) {
    switch( tops.semiLeptonicChannel() ) {
    case WDecay::kElec : eventType = 1; break;
    case WDecay::kMuon : eventType = 2; break;
    case WDecay::kTau  : eventType = 3; break;
    default            : eventType = -1; break;
    }
  }
  else if( tops.isFullLeptonic() ) {
    if     (  tops.fullLeptonicChannel().first  == WDecay::kElec  && 
	      tops.fullLeptonicChannel().second == WDecay::kElec    )    eventType = 11;
    else if( ( tops.fullLeptonicChannel().first  == WDecay::kElec && 
	       tops.fullLeptonicChannel().second == WDecay::kMuon   ) ||
	     ( tops.fullLeptonicChannel().first  == WDecay::kMuon && 
	       tops.fullLeptonicChannel().second == WDecay::kElec   )  ) eventType = 12;
    else if( ( tops.fullLeptonicChannel().first  == WDecay::kElec && 
	       tops.fullLeptonicChannel().second == WDecay::kTau    ) ||
	     ( tops.fullLeptonicChannel().first  == WDecay::kTau && 
	       tops.fullLeptonicChannel().second == WDecay::kElec   )  ) eventType = 13;
    else if(  tops.fullLeptonicChannel().first  == WDecay::kMuon && 
	      tops.fullLeptonicChannel().second == WDecay::kMuon  )      eventType = 14;
    else if( ( tops.fullLeptonicChannel().first  == WDecay::kMuon && 
	       tops.fullLeptonicChannel().second == WDecay::kTau    ) ||
	     ( tops.fullLeptonicChannel().first  == WDecay::kTau  &&
	       tops.fullLeptonicChannel().second  == WDecay::kMuon) )    eventType = 15;
    else if(  tops.fullLeptonicChannel().first  == WDecay::kTau &&
	      tops.fullLeptonicChannel().second == WDecay::kTau   )      eventType = 16;
    else eventType = -2;
  }
  return eventType;
}

// get the final state lepton in tau->lepton events 
// decay chain is: tau(status 3)->tau(status 2)->e/mu(status 1)
reco::GenParticle* 
TopKinematics::getFinalStateLepton(const reco::GenParticle& particle){
  //std::cout << "getFinalStateLepton called!" << std::endl;
  //std::cout << "which has " << particle.numberOfDaughters() << " daughters" << std::endl; 
  // loop daughters
  for(unsigned int i=0; i<particle.numberOfDaughters(); ++i){
    //std::cout << "daughter #" << i << std::endl;
    reco::GenParticle *daughter = (reco::GenParticle *) particle.daughter(i);
    //std::cout << daughter << std::endl;
    if(daughter){
      int ID=daughter->pdgId();
      //std::cout << "ID: " << ID << ", status: " << daughter->status() << std::endl;
      // if daughter is tau
      if(ID==15||ID==-15){
	//std::cout << "tau daughter found!" << std::endl;
	// bubble up
	return getFinalStateLepton(*daughter);
      }
      // if daughter is muon or electron return
      if(ID==-13||ID==13||ID==-11||ID==11){
	//std::cout << "e/#mu found!" << std::endl;
	return daughter;
      }
    }
  }
  return 0;
}
