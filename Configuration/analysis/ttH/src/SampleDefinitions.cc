#include <algorithm>
#include <TString.h>
#include <TColorWheel.h>

#include "SampleDefinitions.h"
#include "Sample.h"





std::map<TString, Sample> SampleDefinitions::samples8TeV(const int mergeLevel)
{
    // Define all samples as differential as they are needed
    // Samples with same legend will appear as one single sample (requires also same colour)
    // They are written to a map to enable sorting and including/excluding of these samples independent of the definitions here
    // Cross sections from:
    // - ttbar from TOP-14-016 (ATLAS+CMS combination)
    // - ttH from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt8TeV
    // - most samples from https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
    // - exceptions defined inline
    std::map<TString, Sample> result;
    
    result["data"] = Sample(
        "Data",
        kBlack,
        1.,
        0., -1.,
        {"run2012A.root", "run2012B.root", "run2012C.root", "run2012D.root"},
        Sample::data
    );
    
    // Several pseudodata samples with different normalisations can be specified
    // Are meant to be modified (reweighted) versions of original samples
    // Should have different names, same legend entry and added to the selectAndOrderSamples8TeV()
    result["pseudodata"] = Sample(
        "Pseudodata",
        kBlack,
        241.5,
        8.5/241.5, -1.,
        // Place for specific reweighted ROOT files to be used instead of standard MC files defined below (one per sample)
        {},
        Sample::pseudodata
    );
    
    // EXAMPLE CONFIGURATION: activated by adding to the list in selectAndOrderSamples8TeV()
    // ttbb sample with 2x larger normalisation in pseudodata than in MC stack
    result["pseudo_ttbarDileptonPlusBbbar"] = Sample(
        "Pseudodata",
        kBlack,
        // 2x increased cross section
        2.*241.5,
        8.5/241.5, -1.,
        {"ttbarDileptonNotauBbbar.root", "ttbarDileptonOnlytauBbbar.root"},
        Sample::pseudodata
    );
    
    result["ttbarDileptonPlusBbbar"] = Sample(
        "t#bar{t}b#bar{b}",
        18,
        241.5,
        8.5/241.5, -1.,
        {"ttbarDileptonNotauBbbar.root", "ttbarDileptonOnlytauBbbar.root"},
        Sample::ttbb
    );
    
    result["ttbarDileptonPlusB"] = Sample(
        "t#bar{t}b",
        12,
        241.5,
        8.5/241.5, -1.,
        {"ttbarDileptonNotauB.root", "ttbarDileptonOnlytauB.root"},
        Sample::ttb
    );
    
    result["ttbarDileptonPlus2B"] = Sample(
        "t#bar{t}2b",
        28,
        // Corrected by Data/Madgraph from JHEP12(2013)039: Z+bb vs dRbb at 7 TeV
        241.5*1.74,
        // Uncertainty from envelope of uncertainties in JHEP12(2013)039: Z+bb vs dRbb at 7 TeV
        (2.43-1.74)/1.74, (1.74-1.)/1.74,
        {"ttbarDileptonNotau2b.root", "ttbarDileptonOnlytau2b.root"},
        Sample::tt2b
    );
    
    result["ttbarDileptonPlusCcbar"] = Sample(
        "t#bar{t} Other",
        23,
        241.5,
        // Plus 50% up/down uncertainty
        0.5, -1.,
        {"ttbarDileptonPlustauCcbar.root"},
        Sample::ttcc
    );
    
    result["ttbarDileptonPlusOther"] = Sample(
        "t#bar{t} Other",
        23,
        241.5,
        8.5/241.5, -1.,
        {"ttbarDileptonPlustauOther.root"},
        Sample::ttother
    );
    
    result["ttbarNotdilepton"] = Sample(
        "t#bar{t} Other",
        23,
        241.5,
        8.5/241.5, -1.,
        {"ttbarNotdilepton.root"},
        Sample::ttNoDilepton
    );
    
    result["singletop"] = Sample(
        "Single t",
        kAzure-9,
        11.1,
        0.76/11.1, -1.,
        {"singletop_tw.root", "singleantitop_tw.root"}
    );
    
    
    result["ww"] = Sample(
        "Diboson",
        10,
        56.0,
        3.05/56., -1.,
        {"wwtoall.root"}
    );
    
    result["wz"] = Sample(
        "Diboson",
        10,
        33.6,
        1.84/33.6, -1.,
        {"wztoall.root"}
    );
    
    result["zz"] = Sample(
        "Diboson",
        10,
        17.,
        0.7/17., -1.,
        {"zztoall.root"}
    );
    
    result["www"] = Sample(
        "Triboson",
        kYellow-7,
        8.058E-2,
        0.047, 0.039,
        {"www.root"}
    );
    
    result["wwz"] = Sample(
        "Triboson",
        kYellow-7,
        5.795E-2,
        0.056, 0.046,
        {"wwz.root"}
    );
    
    result["zzz"] = Sample(
        "Triboson",
        kYellow-7,
        5.527E-3,
        0.027, 0.024,
        {"zzz.root"}
    );
    
    result["dyee1050"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kPink+8,
        860.5,
        // Relative uncertainty as for Mll>20 GeV in the TWiki + 50% (FIXME: which TWiki???)
        39.34/860.5, 38.84/860.5,
        {"dyee1050.root"},
        Sample::dyee
    );
    
    result["dyee50inf"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kPink+8,
        3532.8,
        39.25/3532.8, 38.97/3532.8,
        {"dyee50inf.root"},
        Sample::dyee
    );
    
    result["dymumu1050"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kPink+8,
        860.5,
        // Relative uncertainty as for Mll>20 GeV in the TWiki + 50% (FIXME: which TWiki???)
        39.34/860.5, 38.84/860.5,
        {"dymumu1050.root"},
        Sample::dymumu
    );
    
    result["dymumu50inf"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kPink+8,
        3532.8,
        39.25/3532.8, 38.97/3532.8,
        {"dymumu50inf.root"},
        Sample::dymumu
    );
    
    result["dytautau1050"] = Sample(
        "Z / #gamma* #rightarrow #tau#tau",
        kPink+1,
        860.5,
        // Relative uncertainty as for Mll>20 GeV in the TWiki + 50% (FIXME: which TWiki???)
        39.34/860.5, 38.84/860.5,
        {"dytautau1050.root"},
        Sample::dytautau
    );
    
    result["dytautau50inf"] = Sample(
        "Z / #gamma* #rightarrow #tau#tau",
        kPink+1,
        3532.8,
        39.25/3532.8, 38.97/3532.8,
        {"dytautau50inf.root"},
        Sample::dytautau
    );
    
    result["wlnu"] = Sample(
        "W+Jets",
        kSpring+2,
        36257.2,
        422.3/36257.2, 416.6/36257.2,
        {"wtolnu.root"}
    );
    
    result["qcdmu15"] = Sample(
        "QCD Multijet",
        kOrange-2,
        3.640E8*3.7E-4,
        0.5, -1.,
        {"qcdmu15.root"}
    );
    
    result["qcdmu2030"] = Sample(
        "QCD Multijet",
        kOrange-2,
        2.870E8*6.500E-3,
        0.5, -1.,
        {}
    );
    
    result["qcdmu3050"] = Sample(
        "QCD Multijet",
        kOrange-2,
        6.609E7*12.20E-3,
        0.5, -1.,
        {}
    );
    
    result["qcdmu5080"] = Sample(
        "QCD Multijet",
        kOrange-2,
        8.802E6*21.80E-3,
        0.5, -1.,
        {}
    );
    
    result["qcdmu80120"] = Sample(
        "QCD Multijet",
        kOrange-2,
        1.024E6*39.50E-3,
        0.5, -1.,
        {}
    );
    
    result["qcdmu120170"] = Sample(
        "QCD Multijet",
        kOrange-2,
        1.578E5*47.30E-3,
        0.5, -1.,
        {}
    );
    
    result["qcdem2030"] = Sample(
        "QCD Multijet",
        kOrange-2,
        2.886E8*10.10E-3,
        0.5, -1.,
        {"qcdem2030.root"}
    );
    
    result["qcdem3080"] = Sample(
        "QCD Multijet",
        kOrange-2,
        7.433E7*62.10E-3,
        0.5, -1.,
        {"qcdem3080.root"}
    );
    
    result["qcdem80170"] = Sample(
        "QCD Multijet",
        kOrange-2,
        1.191E6*153.9E-3,
        0.5, -1.,
        {"qcdem80170.root"}
    );
    
    result["qcdbcem2030"] = Sample(
        "QCD Multijet",
        kOrange-2,
        2.886E8*5.800E-4,
        0.5, -1.,
        {"qcdbcem2030.root"}
    );
    
    result["qcdbcem3080"] = Sample(
        "QCD Multijet",
        kOrange-2,
        7.424E7*2.250E-3,
        0.5, -1.,
        {"qcdbcem3080.root"}
    );
    
    result["qcdbcem80170"] = Sample(
        "QCD Multijet",
        kOrange-2,
        1.191E6*10.90E-3,
        0.5, -1.,
        {"qcdbcem80170.root"}
    );
    
    result["ttbarW"] = Sample(
        "t#bar{t}W",
        kTeal+4,
        0.232,
        0.073/0.232, -1.,
        {"ttbarW.root"}
    );
    
    result["ttbarZ"] = Sample(
        "t#bar{t}Z",
        kTeal+1,
        0.2057,
        0.031/0.2057, -1.,
        {"ttbarZ.root"},
        Sample::ttZ
    );
    
    result["ttbargjets"] = Sample(
        "t#bar{t}#gamma",
        kTeal+3,
        // Taken from ../../diLeptonic/src/UsefulTools.cc::SampleXSection()
        // FIXME: Update the uncertainty: 25% from the 7 TeV ATLAS measurement: http://dx.doi.org/10.3204/DESY-PROC-2014-02/24
        1.8,
        0.25, -1.,
        {"ttgjets.root"}
    );
    
    result["ttbarH125inclusiveOther"] = Sample(
        "t#bar{t}H Other",
        kPink-6,
        0.1293,
        0.09, 0.12, 
        {"ttbarH125inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH125inclusiveBbbar"] = Sample(
        "t#bar{t}H (b#bar{b} via incl.)",
        kSpring+9,
        0.1293,
        0.09, 0.12, 
        {"ttbarH125inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH125tobbbar"] = Sample(
        "t#bar{t}H (b#bar{b})",
        2,
        0.1293*0.577,
        0.09*0.577, 0.12*0.577,
        {"ttbarH125tobbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH110inclusiveOther"] = Sample(
        "t#bar{t}H110 Other",
        kTeal+4,
        0.1871,
        0.09, 0.12,
        {"ttbarH110inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH110inclusiveBbbar"] = Sample(
        "t#bar{t}H110 (b#bar{b} via incl.)",
        kSpring+10,
        0.1871,
        0.09, 0.12,
        {"ttbarH110inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH110tobbbar"] = Sample(
        "t#bar{t}H110 (b#bar{b})",
        3,
        0.1871*0.744,
        0.09*0.744, 0.12*0.744,
        {"ttbarH110tobbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH115inclusiveOther"] = Sample(
        "t#bar{t}H115 Other",
        kTeal+5,
        0.1651,
        0.09, 0.12,
        {"ttbarH115inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH115inclusiveBbbar"] = Sample(
        "t#bar{t}H115 (b#bar{b} via incl.)",
        kSpring+11,
        0.1651,
        0.09, 0.12,
        {"ttbarH115inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH115tobbbar"] = Sample(
        "t#bar{t}H115 (b#bar{b})",
        4,
        0.1651*0.703,
        0.09*0.703, 0.12*0.703,
        {"ttbarH115tobbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH120inclusiveOther"] = Sample(
        "t#bar{t}H120 Other",
        kTeal+6,
        0.1459,
        0.09, 0.12,
        {"ttbarH120inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH120inclusiveBbbar"] = Sample(
        "t#bar{t}H120 (b#bar{b} via incl.)",
        kSpring+12,
        0.1459,
        0.09, 0.12,
        {"ttbarH120inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH120tobbbar"] = Sample(
        "t#bar{t}H120 (b#bar{b})",
        5,
        0.1459*0.648,
        0.09*0.648, 0.12*0.648,
        {"ttbarH120tobbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH1225inclusiveOther"] = Sample(
        "t#bar{t}H122.5 Other",
        kTeal+7,
        0.1373,
        0.09, 0.12,
        {"ttbarH1225inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH1225inclusiveBbbar"] = Sample(
        "t#bar{t}H122.5 (b#bar{b} via incl.)",
        kSpring+13,
        0.1373,
        0.09, 0.12,
        {"ttbarH1225inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH1275inclusiveOther"] = Sample(
        "t#bar{t}H127.5 Other",
        kTeal+8,
        0.1218,
        0.09, 0.12,
        {"ttbarH1275inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH1275inclusiveBbbar"] = Sample(
        "t#bar{t}H127.5 (b#bar{b} via incl.)",
        kSpring+14,
        0.1218,
        0.09, 0.12,
        {"ttbarH1275inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH130inclusiveOther"] = Sample(
        "t#bar{t}H130 Other",
        kTeal+9,
        0.1149,
        0.09, 0.12,
        {"ttbarH130inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH130inclusiveBbbar"] = Sample(
        "t#bar{t}H130 (b#bar{b} via incl.)",
        kSpring+15,
        0.1149,
        0.09, 0.12,
        {"ttbarH130inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH130tobbbar"] = Sample(
        "t#bar{t}H130 (b#bar{b})",
        6,
        0.1149*0.494,
        0.09*0.494, 0.12*0.494,
        {"ttbarH130tobbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH135inclusiveOther"] = Sample(
        "t#bar{t}H135 Other",
        kTeal+10,
        0.1024,
        0.09, 0.12,
        {"ttbarH135inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH135inclusiveBbbar"] = Sample(
        "t#bar{t}H135 (b#bar{b} via incl.)",
        kSpring+16,
        0.1024,
        0.09, 0.12,
        {"ttbarH135inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH135tobbbar"] = Sample(
        "t#bar{t}H135 (b#bar{b})",
        7,
        0.1024*0.404,
        0.09*0.404, 0.12*0.404,
        {"ttbarH135tobbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH140inclusiveOther"] = Sample(
        "t#bar{t}H140 Other",
        kTeal+11,
        0.09150,
        0.09, 0.12,
        {"ttbarH140inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH140inclusiveBbbar"] = Sample(
        "t#bar{t}H140 (b#bar{b} via incl.)",
        kSpring+17,
        0.09150,
        0.09, 0.12,
        {"ttbarH140inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    
    // Update legends of samples that should be merged
    updateMergedLegendEntries(result, mergedLegendEntries(mergeLevel));

    return result;
}



std::vector<TString> SampleDefinitions::selectAndOrderSamples8TeV(const int pseudodata)
{
    std::vector<TString> result = {
        "data",
        "qcdmu15",
        "qcdmu2030",
        "qcdmu3050",
        "qcdmu5080",
        "qcdmu80120",
        "qcdmu120170",
        "qcdem2030",
        "qcdem3080",
        "qcdem80170",
        "qcdbcem2030",
        "qcdbcem3080",
        "qcdbcem80170",
        "wlnu",
        "dyee1050",
        "dyee50inf",
        "dymumu1050",
        "dymumu50inf",
        "dytautau1050",
        "dytautau50inf",
        "singletop",
        "ww",
        "wz",
        "zz",
        "www",
        "wwz",
        "zzz",
        "ttbarNotdilepton",
        "ttbarDileptonPlusOther",
        "ttbarDileptonPlusCcbar",
        "ttbarDileptonPlus2B",
        "ttbarDileptonPlusB",
        "ttbarDileptonPlusBbbar",
        "ttbargjets",
        "ttbarW",
        "ttbarZ",
        "ttbarH125inclusiveOther",
        "ttbarH125inclusiveBbbar",
        "ttbarH125tobbbar",
        "ttbarH110inclusiveOther",
        "ttbarH110inclusiveBbbar",
        "ttbarH110tobbbar",
        "ttbarH115inclusiveOther",
        "ttbarH115inclusiveBbbar",
        "ttbarH115tobbbar",
        "ttbarH120inclusiveOther",
        "ttbarH120inclusiveBbbar",
        "ttbarH120tobbbar",
        "ttbarH1225inclusiveOther",
        "ttbarH1225inclusiveBbbar",
        "ttbarH1275inclusiveOther",
        "ttbarH1275inclusiveBbbar",
        "ttbarH130inclusiveOther",
        "ttbarH130inclusiveBbbar",
        "ttbarH130tobbbar",
        "ttbarH135inclusiveOther",
        "ttbarH135inclusiveBbbar",
        "ttbarH135tobbbar",
        "ttbarH140inclusiveOther",
        "ttbarH140inclusiveBbbar"
    };
    
    if(pseudodata > 0) result.at(0) = "pseudodata";
    
    return result;
}



std::map<TString, Sample> SampleDefinitions::samples13TeV(const int mergeLevel)
{
    // Define all samples as differential as they are needed
    // Samples with same legend will appear as one single sample (requires also same colour)
    // They are written to a map to enable sorting and including/excluding of these samples independent of the definitions here
    // Cross sections from:
    // - ttbar from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
    // - ttH from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV
    // - QCD, ttbarW, ttbarZ from https://cms-pdmv.cern.ch/mcm/requests?page=0
    // - most samples from https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
    // - exceptions defined inline
    // FIXME: Update everywhere below correct cross sections AND uncertainties
    std::map<TString, Sample> result;
    
    result["data"] = Sample(
        "Data",
        kBlack,
        1.,
        0., -1.,
        {"run2015B.root"},
        Sample::data
    );
    
    // Several pseudodata samples with different normalisations can be specified
    // Are meant to be modified (reweighted) versions of original samples
    // Should have different names, same legend entry and added to the selectAndOrderSamples13TeV()
    result["pseudodata"] = Sample(
        "Pseudodata",
        kBlack,
        831.76,
        8.5/831.76, -1.,  //FIXME: Update uncertainties
        // Place for specific reweighted ROOT files to be used instead of standard MC files defined below (one per sample)
        {},
        Sample::pseudodata
    );
    
    // EXAMPLE CONFIGURATION: activated by adding to the list in selectAndOrderSamples13TeV()
    // ttbb sample with 2X larger normalisation in pseudodata than in MC stack
    result["pseudo_ttbarDileptonPlusBbbar"] = Sample(
        "Pseudodata",
        kBlack,
        // Scale cross section by factor 2
        2.*831.76,
        8.5/831.76, -1.,  //FIXME: Update uncertainties
        {"ttbarDileptonNotauBbbar.root", "ttbarDileptonOnlytauBbbar.root"},
        Sample::pseudodata
    );
    
    
    result["ttbarPlusBbbar"] = Sample(
        "t#bar{t}b#bar{b}",
        18,
        831.76,
        8.5/831.76, -1., //FIXME: Update uncertainties
        {"ttbarDileptonNotauBbbar.root", "ttbarDileptonOnlytauBbbar.root", "ttbarNotdileptonBbbar.root"},
        Sample::ttbb
    );
    
    result["ttbarPlusB"] = Sample(
        "t#bar{t}b",
        12,
        831.76,
        8.5/831.76, -1., //FIXME: Update uncertainties
        {"ttbarDileptonNotauB.root", "ttbarDileptonOnlytauB.root", "ttbarNotdileptonB.root"},
        Sample::ttb
    );
    
    result["ttbarPlus2B"] = Sample(
        "t#bar{t}2b",
        28,
        831.76,
        8.5/831.76, -1., //FIXME: Update uncertainties
        {"ttbarDileptonNotau2b.root", "ttbarDileptonOnlytau2b.root", "ttbarNotdilepton2b.root"},
        Sample::tt2b
    );
    
    result["ttbarPlusCcbar"] = Sample(
        "t#bar{t}c#bar{c}",
        kMagenta-10,
        831.76,
        // Plus 50% up/down uncertainty
        0.5, -1., //FIXME: Update uncertainties
        {"ttbarDileptonPlustauCcbar.root", "ttbarNotdileptonCcbar.root"},
        Sample::ttcc
    );
    
    result["ttbarPlusOther"] = Sample(
        "t#bar{t}Other",
        23,
        831.76,
        8.5/831.76, -1., //FIXME: Update uncertainties
        {"ttbarDileptonPlustauOther.root", "ttbarNotdileptonOther.root"},
        Sample::ttother
    );
    
    result["singletop"] = Sample(
        "Single Top",
        kAzure-9,
        35.6,
        1.92/35.6, -1.,
        {"singletop_tw.root", "singleantitop_tw.root"}
    );
    
    result["ww"] = Sample(
        "Diboson",
        10,
        11.0,
        3.05/11., -1., //FIXME: Update uncertainties
        {"wwtoall.root"}
    );
    
    result["wz"] = Sample(
        "Diboson",
        10,
        66.1,
        1.84/66.1, -1., //FIXME: Update uncertainties
        {"wztoall.root"}
    );
    
    result["zz"] = Sample(
        "Diboson",
        10,
        31.8,     
        0.7/31.8, -1., //FIXME: Update uncertainties
        {"zztoall.root"}
    );
    
    result["www"] = Sample(
        "Triboson",
        kYellow-7,
        8.058E-2,      //FIXME: Update XS
        0.047, 0.039, //FIXME: Update uncertainties
        {"www.root"}
    );
    
    result["wwz"] = Sample(
        "Triboson",
        kYellow-7,
        5.795E-2,      //FIXME: Update XS
        0.056, 0.046, //FIXME: Update uncertainties
        {"wwz.root"}
    );
    
    result["zzz"] = Sample(
        "Triboson",
        kYellow-7,
        5.527E-3,  //FIXME: Update XS
        0.027, 0.024,  //FIXME: Update uncertainties
        {"zzz.root"}
    );
    
    result["dyee1050"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kPink+8,
        3.*860.5,
        39.34/860.5, 38.84/860.5,  //FIXME: Update uncertainties
        {"dyee1050.root"},
        Sample::dyee
    );
    
    result["dyee50inf"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kPink+8,
        3.*2008.4,
        39.25/(3.*2008.4), 38.97/(3.*2008.4),  //FIXME: Update uncertainties
        {"dyee50inf.root"},
        Sample::dyee
    );
    
    result["dymumu1050"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kPink+8,
        3*860.5,
        39.34/(3*860.5), 38.84/(3*860.5),  //FIXME: Update uncertainties
        {"dymumu1050.root"},
        Sample::dymumu
    );
    
    result["dymumu50inf"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kPink+8,
        3.*2008.4,
        39.25/(3.*2008.4), 38.97/(3.*2008.4),  //FIXME: Update uncertainties
        {"dymumu50inf.root"},
        Sample::dymumu
    );
    
    result["dytautau1050"] = Sample(
        "Z / #gamma* #rightarrow #tau#tau",
        kPink+1,
        3*860.5,
        39.34/(3*860.5), 38.84/(3*860.5),  //FIXME: Update uncertainties
        {"dytautau1050.root"},
        Sample::dytautau
    );
    
    result["dytautau50inf"] = Sample(
        "Z / #gamma* #rightarrow #tau#tau",
        kPink+1,
        3.*2008.4,
        39.25/(3.*2008.4), 38.97/(3.*2008.4),  //FIXME: Update uncertainties
        {"dytautau50inf.root"},
        Sample::dytautau
    );
    
    result["wlnu"] = Sample(
        "W+Jets",
        kSpring+2,
        61526.7,
        422.3/61526.7, 416.6/61526.7,  //FIXME: Update uncertainties
        {"wtolnu.root"}
    );
    
    result["qcdmu15"] = Sample(
        "QCD Multijet",
        kOrange-2,
        3.640E8*3.7E-4,
        0.5, -1.,  
        {"qcdmu15.root"}
    );
    
    result["qcdmu2030"] = Sample(
        "QCD Multijet",
        kOrange-2,
        2.870E8*6.500E-3,
        0.5, -1.,
        {}
    );
    
    result["qcdmu3050"] = Sample(
        "QCD Multijet",
        kOrange-2,
        6.609E7*12.20E-3,
        0.5, -1.,
        {}
    );
    
    result["qcdmu5080"] = Sample(
        "QCD Multijet",
        kOrange-2,
        8.802E6*21.80E-3,
        0.5, -1.,
        {}
    );
    
    result["qcdmu80120"] = Sample(
        "QCD Multijet",
        kOrange-2,
        1.024E6*39.50E-3,
        0.5, -1.,
        {}
    );
    
    result["qcdmu120170"] = Sample(
        "QCD Multijet",
        kOrange-2,
        1.578E5*47.30E-3,
        0.5, -1.,
        {}
    );
    
    result["qcdem2030"] = Sample(
        "QCD Multijet",
        kOrange-2,
        557600000*0.0096,
        0.5, -1.,
        {"qcdem2030.root"}
    );
    
    result["qcdem3050"] = Sample(
        "QCD Multijet",
        kOrange-2,
        136000000*0.073,
        0.5, -1.,
        {"qcdem3050.root"}
    );
    
    result["qcdem5080"] = Sample(
        "QCD Multijet",
        kOrange-2,
        19800000*0.146,
        0.5, -1.,
        {"qcdem5080.root"}
    );
    
    result["qcdem80120"] = Sample(
        "QCD Multijet",
        kOrange-2,
        2800000*0.125,
        0.5, -1.,
        {"qcdem80120.root"}
    );
    
    result["qcdem120170"] = Sample(
        "QCD Multijet",
        kOrange-2,
        477000*0.132,
        0.5, -1.,
        {"qcdem120170.root"}
    );
    
    result["qcdem170300"] = Sample(
        "QCD Multijet",
        kOrange-2,
        114000*0.165,
        0.5, -1.,
        {"qcdem170300.root"}
    );
    
    result["qcdem300inf"] = Sample(
        "QCD Multijet",
        kOrange-2,
        9000*0.15,
        0.5, -1.,
        {"qcdem300inf.root"}
    );
    
    result["qcdbcem2030"] = Sample(
        "QCD Multijet",
        kOrange-2,
        2.886E8*5.800E-4,
        0.5, -1.,
        {"qcdbcem2030.root"}
    );
    
    result["qcdbcem3080"] = Sample(
        "QCD Multijet",
        kOrange-2,
        7.424E7*2.250E-3,
        0.5, -1.,
        {"qcdbcem3080.root"}
    );
    
    result["qcdbcem80170"] = Sample(
        "QCD Multijet",
        kOrange-2,
        1.191E6*10.90E-3,
        0.5, -1.,
        {"qcdbcem80170.root"}
    );
    
    result["ttbarW"] = Sample(
        "t#bar{t}W",
        kTeal+4,
        1.152,
        0.073/1.152, -1., //FIXME: Update uncertainties
        {"ttbarW.root"}
    );
    
    result["ttbarZ"] = Sample(
        "t#bar{t}Z",
        kTeal+1,
        2.232,
        0.031/2.232, -1.,  //FIXME: Update uncertainties
        {"ttbarZ.root"},
        Sample::ttZ
    );
    
    result["ttbargjets"] = Sample(
        "t#bar{t}#gamma",
        kTeal+3,
        1.8,         //FIXME: Update XS
        0.25, -1.,   //FIXME: Update uncertainties
        {"ttgjets.root"}
    );
    
    result["ttbarH125inclusiveOther"] = Sample(
        "t#bar{t}H Other",
        kSpring+9,
        0.5085,
        0.09, 0.12,  //FIXME: Update uncertainties
        {"ttbarH125inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH125inclusiveBbbar"] = Sample(
        "t#bar{t}H (b#bar{b} via incl.)",
        kRed,
        0.5085,
        0.09, 0.12,  //FIXME: Update uncertainties
        {"ttbarH125inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH125tobbbar"] = Sample(
        "t#bar{t}H (b#bar{b})",
        2,
        0.5085*0.577,
        0.09*0.577, 0.12*0.577,  //FIXME: Update uncertainties
        {"ttbarH125tobbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH120inclusiveOther"] = Sample(
        "t#bar{t}H120 Other",
        kTeal+6,
        0.1459,
        0.09, 0.12,   //FIXME: Update uncertainties
        {"ttbarH120inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH120inclusiveBbbar"] = Sample(
        "t#bar{t}H120 (b#bar{b} via incl.)",
        kSpring+12,
        0.1459,
        0.09, 0.12,  //FIXME: Update uncertainties
        {"ttbarH120inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH120tobbbar"] = Sample(
        "t#bar{t}H120 (b#bar{b})",
        5,
        0.1459*0.648,
        0.09*0.648, 0.12*0.648,  //FIXME: Update uncertainties
        {"ttbarH120tobbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH1225inclusiveOther"] = Sample(
        "t#bar{t}H122.5 Other",
        kTeal+7,
        0.1373,
        0.09, 0.12,  //FIXME: Update uncertainties
        {"ttbarH1225inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH1225inclusiveBbbar"] = Sample(
        "t#bar{t}H122.5 (b#bar{b} via incl.)",
        kSpring+13,
        0.1373,
        0.09, 0.12,  //FIXME: Update uncertainties
        {"ttbarH1225inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH1275inclusiveOther"] = Sample(
        "t#bar{t}H127.5 Other",
        kTeal+8,
        0.1218,
        0.09, 0.12,  //FIXME: Update uncertainties
        {"ttbarH1275inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH1275inclusiveBbbar"] = Sample(
        "t#bar{t}H127.5 (b#bar{b} via incl.)",
        kSpring+14,
        0.1218,
        0.09, 0.12,  //FIXME: Update uncertainties
        {"ttbarH1275inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH130inclusiveOther"] = Sample(
        "t#bar{t}H130 Other",
        kTeal+9,
        0.1149,
        0.09, 0.12,  //FIXME: Update uncertainties
        {"ttbarH130inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH130inclusiveBbbar"] = Sample(
        "t#bar{t}H130 (b#bar{b} via incl.)",
        kSpring+15,
        0.1149,
        0.09, 0.12,  //FIXME: Update uncertainties
        {"ttbarH130inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH130tobbbar"] = Sample(
        "t#bar{t}H130 (b#bar{b})",
        6,
        0.1149*0.494,
        0.09*0.494, 0.12*0.494,  //FIXME: Update uncertainties
        {"ttbarH130tobbbar.root"},
        Sample::ttHbb
    );
    
    
    // Update legends of samples that should be merged
    updateMergedLegendEntries(result, mergedLegendEntries(mergeLevel));

    return result;
}



std::vector<TString> SampleDefinitions::selectAndOrderSamples13TeV(const int pseudodata)
{
    std::vector<TString> result = {
        "data",
        "qcdmu15",
        "qcdmu2030",
        "qcdmu3050",
        "qcdmu5080",
        "qcdmu80120",
        "qcdmu120170",
        "qcdem2030",
        "qcdem3050",
        "qcdem5080",
        "qcdem80120",
        "qcdem120170",
        "qcdem170300",
        "qcdem300inf",
        "qcdbcem2030",
        "qcdbcem3080",
        "qcdbcem80170",
        "wlnu",
        "dyee1050",
        "dyee50inf",
        "dymumu1050",
        "dymumu50inf",
        "dytautau1050",
        "dytautau50inf",
        "singletop",
        "ww",
        "wz",
        "zz",
        "www",
        "wwz",
        "zzz",
        "ttbarPlusOther",
        "ttbarPlusCcbar",
        "ttbarPlus2B",
        "ttbarPlusB",
        "ttbarPlusBbbar",
        "ttbargjets",
        "ttbarW",
        "ttbarZ",
        "ttbarH125inclusiveOther",
        "ttbarH125inclusiveBbbar",
        "ttbarH125tobbbar",
        "ttbarH120inclusiveOther",
        "ttbarH120inclusiveBbbar",
        "ttbarH120tobbbar",
        "ttbarH1225inclusiveOther",
        "ttbarH1225inclusiveBbbar",
        "ttbarH1275inclusiveOther",
        "ttbarH1275inclusiveBbbar",
        "ttbarH130inclusiveOther",
        "ttbarH130inclusiveBbbar",
        "ttbarH130tobbbar",
    };
    
    if(pseudodata > 0) result.at(0) = "pseudodata";
    
    return result;
}



std::vector<TString> SampleDefinitions::legendList(const std::map<TString, Sample>& samples, const std::vector<TString>& sampleIdentifiers)
{
    std::vector<TString> v_legend;
    
    for(TString sampleIdentifier : sampleIdentifiers) {
        TString legend = samples.at(sampleIdentifier).legendEntry();
        if(std::find(v_legend.begin(), v_legend.end(), legend) != v_legend.end()) continue;
        v_legend.push_back(legend);
    }
    
    return v_legend;
}



bool SampleDefinitions::usingPseudodata(const std::map<TString, Sample>& samples, const std::vector<TString>& sampleIdentifiers)
{
    for(auto nameSamplePair : samples) {
        if(nameSamplePair.second.sampleType() != Sample::pseudodata) continue;
        if(std::find(sampleIdentifiers.begin(), sampleIdentifiers.end(), nameSamplePair.first) != sampleIdentifiers.end()) return true;
    }
    
    return false;
}



std::map<TString, SampleDefinitions::ColorLegends> SampleDefinitions::mergedLegendEntries(const int mergeLevel)
{
    std::map<TString, ColorLegends> m_legends;
    
    if(mergeLevel < 1) return m_legends;
    
    // ############################################ Meging level: 1
    if(mergeLevel >= 1) {
        // Merging Electroweak processes
        m_legends["Electroweak"] = ColorLegends(kPink+6, std::set<TString>({
            "Diboson", "Triboson", "Z / #gamma* #rightarrow ee/#mu#mu", "Z / #gamma* #rightarrow #tau#tau", "W+Jets"
        }));
        // Merging ttW and ttgamma
        m_legends["t#bar{t}W/#gamma"] = ColorLegends(kTeal+4, std::set<TString>({
            "t#bar{t}W", "t#bar{t}#gamma"
        }));
    }
    if(mergeLevel == 1) return m_legends;
    
    // ############################################ Merging level: 2
    std::map<TString, ColorLegends > m_legends_2;
    if(mergeLevel >= 2) {
        // Merging Electroweak, QCD, Single t, ttW and ttgamma
        m_legends_2["Minor bkg."] = ColorLegends(kOrange-4, m_legends.at("Electroweak").second);
        m_legends_2.at("Minor bkg.").second.insert({
            "QCD Multijet", "t#bar{t}W", "t#bar{t}#gamma", "Single t"
        });
    }
    if(mergeLevel == 2) return m_legends_2;
    
    return m_legends;
}



void SampleDefinitions::updateMergedLegendEntries(std::map<TString, Sample>& m_nameSample, 
                                                  const std::map<TString, ColorLegends>& m_mergedLegends)
{
    // Looping over the list of merged legend entries
    for(const auto& legendLegends : m_mergedLegends) {
        const TString& legend_merged = legendLegends.first;
        const ColorLegends& colorLegends = legendLegends.second;
        // Looping over all samples to update their legends appropriately
        int groupColor = colorLegends.first;
        for(auto& nameSample : m_nameSample) {
            Sample& sample = nameSample.second;
            if(colorLegends.second.count(sample.legendEntry()) < 1) continue;
            // Setting the color of the group to all merged samples
            sample.setColor(groupColor);
            sample.setLegendEntry(legend_merged);
        }
    }
}





