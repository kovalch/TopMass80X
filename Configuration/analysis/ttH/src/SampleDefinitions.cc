#include <TString.h>
#include <TColorWheel.h>

#include "SampleDefinitions.h"
#include "Sample.h"





std::map<TString, Sample> SampleDefinitions::samples8TeV()
{
    // Define all samples as differential as they are needed
    // Samples with same legend will appear as one single sample (requires also same colour)
    // They are written to a map to enable sorting and including/excluding of these samples independent of the definitions here
    std::map<TString, Sample> result;
    
    result["data"] = Sample(
        "Data",
        kBlack,
        1.,
        {"run2012A.root", "run2012B.root", "run2012C.root", "run2012D.root"},
        Sample::data
    );
    
    result["ttbarsignalPlusBbbar"] = Sample(
        "t#bar{t}b#bar{b}",
        18,
        234.0,
        {"ttbarsignalPlusBbbar.root"},
        Sample::ttbb
    );

    result["ttbarsignalPlusB"] = Sample(
        "t#bar{t}b",
        12,
        234.0,
        {"ttbarsignalPlusB.root"},
        Sample::ttb
    );

    result["ttbarsignalPlus2B"] = Sample(
        "t#bar{t}2b",
        28,
        234.0,
        {"ttbarsignalPlus2B.root"},
        Sample::tt2b
    );

    result["ttbarsignalPlusCcbar"] = Sample(
        "t#bar{t}cc",
        29,
        234.0,
        {"ttbarsignalPlusCcbar.root"},
        Sample::ttcc
    );

    result["ttbarsignalPlusOther"] = Sample(
        "t#bar{t}Other",
        23,
        234.0,
        {"ttbarsignalPlusOther.root"},
        Sample::ttother
    );
    
    result["ttbarbkg"] = Sample(
        "t#bar{t}Other",
        23,
        234.0,
        {"ttbarbg.root"},
        Sample::ttbg
    );
    
    result["singletop"] = Sample(
        "Single Top",
        kViolet-3,
        11.1,
        {"singletop_tw.root", "singleantitop_tw.root"}
    );
    
    result["ww"] = Sample(
        "Diboson",
        10,
        54.838,
        {"wwtoall.root"}
    );
    
    result["wz"] = Sample(
        "Diboson",
        10,
        33.21,
        {"wztoall.root"}
    );
    
    result["zz"] = Sample(
        "Diboson",
        10,
        17.654,
        {"zztoall.root"}
    );
    
    result["dyee1050"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kAzure+2,
        860.5,
        {"dyee1050.root"},
        Sample::dyee
    );
    
    result["dyee50inf"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kAzure+2,
        3532.8,
        {"dyee50inf.root"},
        Sample::dyee
    );
    
    result["dymumu1050"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kAzure+2,
        860.5,
        {"dymumu1050.root"},
        Sample::dymumu
    );
    
    result["dymumu50inf"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kAzure+2,
        3532.8,
        {"dymumu50inf.root"},
        Sample::dymumu
    );
    
    result["dytautau1050"] = Sample(
        "Z / #gamma* #rightarrow #tau#tau",
        kAzure+10,
        860.5,
        {"dytautau1050.root"},
        Sample::dytautau
    );
    
    result["dytautau50inf"] = Sample(
        "Z / #gamma* #rightarrow #tau#tau",
        kAzure+10,
        3532.8,
        {"dytautau50inf.root"},
        Sample::dytautau
    );
    
    result["wlnu"] = Sample(
        "W+Jets",
        kSpring+2,
        36257.2,
        {"wtolnu.root"}
    );
    
    result["qcdmu15"] = Sample(
        "QCD Multijet",
        kOrange-2,
        3.640E8*3.7E-4,
        {"qcdmu15.root"}
    );
    
    result["qcdmu2030"] = Sample(
        "QCD Multijet",
        kOrange-2,
        2.870E8*6.500E-3,
        {}
    );
    
    result["qcdmu3050"] = Sample(
        "QCD Multijet",
        kOrange-2,
        6.609E7*12.20E-3,
        {}
    );
    
    result["qcdmu5080"] = Sample(
        "QCD Multijet",
        kOrange-2,
        8.802E6*21.80E-3,
        {}
    );
    
    result["qcdmu80120"] = Sample(
        "QCD Multijet",
        kOrange-2,
        1.024E6*39.50E-3,
        {}
    );
    
    result["qcdmu120170"] = Sample(
        "QCD Multijet",
        kOrange-2,
        1.578E5*47.30E-3,
        {}
    );
    
    result["qcdem2030"] = Sample(
        "QCD Multijet",
        kOrange-2,
        2.886E8*10.10E-3,
        {"qcdem2030.root"}
    );
    
    result["qcdem3080"] = Sample(
        "QCD Multijet",
        kOrange-2,
        7.433E7*62.10E-3,
        {"qcdem3080.root"}
    );
    
    result["qcdem80170"] = Sample(
        "QCD Multijet",
        kOrange-2,
        1.191E6*153.9E-3,
        {"qcdem80170.root"}
    );
    
    result["qcdbcem2030"] = Sample(
        "QCD Multijet",
        kOrange-2,
        2.886E8*5.800E-4,
        {"qcdbcem2030.root"}
    );
    
    result["qcdbcem3080"] = Sample(
        "QCD Multijet",
        kOrange-2,
        7.424E7*2.250E-3,
        {"qcdbcem3080.root"}
    );
    
    result["qcdbcem80170"] = Sample(
        "QCD Multijet",
        kOrange-2,
        1.191E6*10.90E-3,
        {"qcdbcem80170.root"}
    );
    
    result["ttbarW"] = Sample(
        "t#bar{t}W",
        kViolet-4,
        0.232,
        {"ttbarW.root"}
    );
    
    result["ttbarZ"] = Sample(
        "t#bar{t}Z",
        kTeal+1,
        0.2057,
        {"ttbarZ.root"},
        Sample::ttZ
    );
    
    result["ttbarH125inclusiveOther"] = Sample(
        "t#bar{t}H Other",
        kTeal+3,
        0.1293,
        {"ttbarH125inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH125inclusiveBbbar"] = Sample(
        "t#bar{t}H (b#bar{b} via incl.)",
        kSpring+9,
        0.1293,
        {"ttbarH125inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH125tobbbar"] = Sample(
        "t#bar{t}H (b#bar{b})",
        2,
        0.1293*0.577,
        {"ttbarH125tobbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH110inclusiveOther"] = Sample(
        "t#bar{t}H110 Other",
        kTeal+4,
        0.1871,
        {"ttbarH110inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH110inclusiveBbbar"] = Sample(
        "t#bar{t}H110 (b#bar{b} via incl.)",
        kSpring+10,
        0.1871,
        {"ttbarH110inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH110tobbbar"] = Sample(
        "t#bar{t}H110 (b#bar{b})",
        3,
        0.1871*0.744,
        {"ttbarH110tobbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH115inclusiveOther"] = Sample(
        "t#bar{t}H115 Other",
        kTeal+5,
        0.1651,
        {"ttbarH115inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH115inclusiveBbbar"] = Sample(
        "t#bar{t}H115 (b#bar{b} via incl.)",
        kSpring+11,
        0.1651,
        {"ttbarH115inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH115tobbbar"] = Sample(
        "t#bar{t}H115 (b#bar{b})",
        4,
        0.1651*0.703,
        {"ttbarH115tobbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH120inclusiveOther"] = Sample(
        "t#bar{t}H120 Other",
        kTeal+6,
        0.1459,
        {"ttbarH120inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH120inclusiveBbbar"] = Sample(
        "t#bar{t}H120 (b#bar{b} via incl.)",
        kSpring+12,
        0.1459,
        {"ttbarH120inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH120tobbbar"] = Sample(
        "t#bar{t}H120 (b#bar{b})",
        5,
        0.1459*0.648,
        {"ttbarH120tobbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH1225inclusiveOther"] = Sample(
        "t#bar{t}H122.5 Other",
        kTeal+7,
        0.1373,
        {"ttbarH1225inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH1225inclusiveBbbar"] = Sample(
        "t#bar{t}H122.5 (b#bar{b} via incl.)",
        kSpring+13,
        0.1373,
        {"ttbarH1225inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH1275inclusiveOther"] = Sample(
        "t#bar{t}H127.5 Other",
        kTeal+8,
        0.1218,
        {"ttbarH1275inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH1275inclusiveBbbar"] = Sample(
        "t#bar{t}H127.5 (b#bar{b} via incl.)",
        kSpring+14,
        0.1218,
        {"ttbarH1275inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH130inclusiveOther"] = Sample(
        "t#bar{t}H130 Other",
        kTeal+9,
        0.1149,
        {"ttbarH130inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH130inclusiveBbbar"] = Sample(
        "t#bar{t}H130 (b#bar{b} via incl.)",
        kSpring+15,
        0.1149,
        {"ttbarH130inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH130tobbbar"] = Sample(
        "t#bar{t}H130 (b#bar{b})",
        6,
        0.1149*0.494,
        {"ttbarH130tobbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH135inclusiveOther"] = Sample(
        "t#bar{t}H135 Other",
        kTeal+10,
        0.1024,
        {"ttbarH135inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH135inclusiveBbbar"] = Sample(
        "t#bar{t}H135 (b#bar{b} via incl.)",
        kSpring+16,
        0.1024,
        {"ttbarH135inclusiveBbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH135tobbbar"] = Sample(
        "t#bar{t}H135 (b#bar{b})",
        7,
        0.1024*0.404,
        {"ttbarH135tobbbar.root"},
        Sample::ttHbb
    );
    
    result["ttbarH140inclusiveOther"] = Sample(
        "t#bar{t}H140 Other",
        kTeal+11,
        0.09150,
        {"ttbarH140inclusiveOther.root"},
        Sample::ttHother
    );
    
    result["ttbarH140inclusiveBbbar"] = Sample(
        "t#bar{t}H140 (b#bar{b} via incl.)",
        kSpring+17,
        0.09150,
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








