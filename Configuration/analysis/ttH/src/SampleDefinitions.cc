#include <TString.h>
#include <TColorWheel.h>

#include "SampleDefinitions.h"
#include "Sample.h"





std::map<TString, Sample> SampleDefinitions::samples8TeV()
{
    // Define all samples as differential as they are needed
    // Samples with same legend will appear as one single sample (requires also same colour)
    // They are written to a map to enable sorting and including/excluding of these samples independent of the definitions here
    // Cross sections from https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV with exceptions defined inline
    std::map<TString, Sample> result;
    
    result["data"] = Sample(
        "Data",
        kBlack,
        1.,
        0., -1.,
        {"run2012A.root", "run2012B.root", "run2012C.root", "run2012D.root"},
        Sample::data
    );
    
    result["pseudodata"] = Sample(
        "pseudoData [nominal]",
        kBlack,
        241.5,
        8.5/241.5, -1.,
        {   // Place for specific reweighted ROOT files to be used instead of standard MC files defined below
            "ttbarsignalPlusBbbar_nominal.root"
        },
        Sample::pseudodata
    );
    
    result["ttbarsignalPlusBbbar"] = Sample(
        "t#bar{t}b#bar{b}",
        18,
        // TOP-14-016 (ATLAS+CMS combination)
        241.5,
        8.5/241.5, -1.,
        {
            "ttbarsignalPlusBbbar.root",
            "ttbarsignalPlusBbbar_massup.root",
            "ttbarsignalPlusBbbar_massdown.root",
            "ttbarsignalPlusBbbar_matchingup.root",
            "ttbarsignalPlusBbbar_matchingdown.root",
            "ttbarsignalPlusBbbar_scaleup.root",
            "ttbarsignalPlusBbbar_scaledown.root",
            "ttbarsignalPlusBbbar_powheg.root",
            "ttbarsignalPlusBbbar_powhegHerwig.root",
            "ttbarsignalPlusBbbar_mcatnlo.root",
            "ttbarsignalPlusBbbar_Perugia11.root",
            "ttbarsignalPlusBbbar_Perugia11NoCR.root",
        },
        Sample::ttbb
    );

    result["ttbarsignalPlusB"] = Sample(
        "t#bar{t}b",
        12,
        // TOP-14-016 (ATLAS+CMS combination)
        241.5,
        8.5/241.5, -1.,
        {
            "ttbarsignalPlusB.root",
            "ttbarsignalPlusB_massup.root",
            "ttbarsignalPlusB_massdown.root",
            "ttbarsignalPlusB_matchingup.root",
            "ttbarsignalPlusB_matchingdown.root",
            "ttbarsignalPlusB_scaleup.root",
            "ttbarsignalPlusB_scaledown.root",
            "ttbarsignalPlusB_powheg.root",
            "ttbarsignalPlusB_powhegHerwig.root",
            "ttbarsignalPlusB_mcatnlo.root",
            "ttbarsignalPlusB_Perugia11.root",
            "ttbarsignalPlusB_Perugia11NoCR.root",
        },
        Sample::ttb
    );

    result["ttbarsignalPlus2B"] = Sample(
        "t#bar{t}2b",
        28,
        // TOP-14-016 (ATLAS+CMS combination)
        // Corrected by Data/Madgraph from JHEP12(2013)039: Z+bb vs dRbb at 7 TeV
        241.5*1.74,
        // Uncertainty from envelope of uncertainties in JHEP12(2013)039: Z+bb vs dRbb at 7 TeV
        (2.43-1.74)/1.74, (1.74-1.23)/1.74,
        {
            "ttbarsignalPlus2B.root",
            "ttbarsignalPlus2B_massup.root",
            "ttbarsignalPlus2B_massdown.root",
            "ttbarsignalPlus2B_matchingup.root",
            "ttbarsignalPlus2B_matchingdown.root",
            "ttbarsignalPlus2B_scaleup.root",
            "ttbarsignalPlus2B_scaledown.root",
            "ttbarsignalPlus2B_powheg.root",
            "ttbarsignalPlus2B_powhegHerwig.root",
            "ttbarsignalPlus2B_mcatnlo.root",
            "ttbarsignalPlus2B_Perugia11.root",
            "ttbarsignalPlus2B_Perugia11NoCR.root",
        },
        Sample::tt2b
    );

    result["ttbarsignalPlusCcbar"] = Sample(
        "t#bar{t}Other",
        23,
        // TOP-14-016 (ATLAS+CMS combination)
        241.5,
        // Plus 50% up/down uncertainty
        0.5, -1.,
        {
            "ttbarsignalPlusCcbar.root",
            "ttbarsignalPlusCcbar_massup.root",
            "ttbarsignalPlusCcbar_massdown.root",
            "ttbarsignalPlusCcbar_matchingup.root",
            "ttbarsignalPlusCcbar_matchingdown.root",
            "ttbarsignalPlusCcbar_scaleup.root",
            "ttbarsignalPlusCcbar_scaledown.root",
            "ttbarsignalPlusCcbar_powheg.root",
            "ttbarsignalPlusCcbar_powhegHerwig.root",
            "ttbarsignalPlusCcbar_mcatnlo.root",
            "ttbarsignalPlusCcbar_Perugia11.root",
            "ttbarsignalPlusCcbar_Perugia11NoCR.root",
        },
        Sample::ttcc
    );

    result["ttbarsignalPlusOther"] = Sample(
        "t#bar{t}Other",
        23,
        // TOP-14-016 (ATLAS+CMS combination)
        241.5,
        8.5/241.5, -1.,
        {
            "ttbarsignalPlusOther.root",
            "ttbarsignalPlusOther_massup.root",
            "ttbarsignalPlusOther_massdown.root",
            "ttbarsignalPlusOther_matchingup.root",
            "ttbarsignalPlusOther_matchingdown.root",
            "ttbarsignalPlusOther_scaleup.root",
            "ttbarsignalPlusOther_scaledown.root",
            "ttbarsignalPlusOther_powheg.root",
            "ttbarsignalPlusOther_powhegHerwig.root",
            "ttbarsignalPlusOther_mcatnlo.root",
            "ttbarsignalPlusOther_Perugia11.root",
            "ttbarsignalPlusOther_Perugia11NoCR.root",
        },
        Sample::ttother
    );
    
    result["ttbarbkg"] = Sample(
        "t#bar{t}Other",
        23,
        // TOP-14-016 (ATLAS+CMS combination)
        241.5,
        8.5/241.5, -1.,
        {
            "ttbarbg.root",
            "ttbarbg_massup.root",
            "ttbarbg_massdown.root",
            "ttbarbg_matchingup.root",
            "ttbarbg_matchingdown.root",
            "ttbarbg_scaleup.root",
            "ttbarbg_scaledown.root",
            "ttbarbg_powheg.root",
            "ttbarbg_powhegHerwig.root",
            "ttbarbg_mcatnlo.root",
        },
        Sample::ttNoDilepton
    );
    
    result["singletop"] = Sample(
        "Single Top",
        kViolet-3,
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
    
    result["dyee1050"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kAzure+2,
        860.5,
        // Relative uncertainty as for Mll>20 GeV in the TWiki + 50%
        39.34/860.5, 38.84/860.5,
        {"dyee1050.root"},
        Sample::dyee
    );
    
    result["dyee50inf"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kAzure+2,
        3532.8,
        39.25/3532.8, 38.97/3532.8,
        {"dyee50inf.root"},
        Sample::dyee
    );
    
    result["dymumu1050"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kAzure+2,
        860.5,
        // Relative uncertainty as for Mll>20 GeV in the TWiki + 50%
        39.34/860.5, 38.84/860.5,
        {"dymumu1050.root"},
        Sample::dymumu
    );
    
    result["dymumu50inf"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kAzure+2,
        3532.8,
        39.25/3532.8, 38.97/3532.8,
        {"dymumu50inf.root"},
        Sample::dymumu
    );
    
    result["dytautau1050"] = Sample(
        "Z / #gamma* #rightarrow #tau#tau",
        kAzure+10,
        860.5,
        // Relative uncertainty as for Mll>20 GeV in the TWiki + 50%
        39.34/860.5, 38.84/860.5,
        {"dytautau1050.root"},
        Sample::dytautau
    );
    
    result["dytautau50inf"] = Sample(
        "Z / #gamma* #rightarrow #tau#tau",
        kAzure+10,
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
        kViolet-4,
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
    
    result["ttbarH125inclusiveOther"] = Sample(
        "t#bar{t}H Other",
        kTeal+3,
        // ttH XS at https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt8TeV#ttH_Process
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
    
    return result;
}



std::vector<TString> SampleDefinitions::selectAndOrderSamples8TeV()
{
    const std::vector<TString> result = {
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
        "ttbarbkg",
        "ttbarsignalPlusOther",
        "ttbarsignalPlusCcbar",
        "ttbarsignalPlus2B",
        "ttbarsignalPlusB",
        "ttbarsignalPlusBbbar",
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
    
    return result;
}








