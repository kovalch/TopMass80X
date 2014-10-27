#include <TString.h>
#include <Rtypes.h>

#include "SampleDefinitions.h"
#include "Sample.h"





std::map<TString, Sample> SampleDefinitions::samples8TeV()
{
    
     //MC cross sections taken from:
    //  UsefulTools::SampleXSection in diLeptonic/src
    //  AN-12/194    AN-12/228
    double topxsec = 245.102;
    
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
    

    result["ttbarsignal"] = Sample(
        "t#bar{t} Signal",
        kRed+1,
        topxsec,
        {"ttbarsignalplustau.root"},
        Sample::ttother
    );
    
    result["ttbarbkg"] = Sample(
        "t#bar{t}Other",
        kRed-7,
        topxsec,
        {"ttbarbg.root", "ttbarbgviatau.root"},
        Sample::ttother
    );
    
    result["singletop"] = Sample(
        "Single Top",
        kMagenta,
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
        kAzure-2,
        860.5,
        {"dyee1050.root"},
        Sample::dyee
    );
    
    result["dyee50inf"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kAzure-2,
        3532.8,
        {"dyee50inf.root"},
        Sample::dyee
    );
    
    result["dymumu1050"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kAzure-2,
        860.5,
        {"dymumu1050.root"},
        Sample::dymumu
    );
    
    result["dymumu50inf"] = Sample(
        "Z / #gamma* #rightarrow ee/#mu#mu",
        kAzure-2,
        3532.8,
        {"dymumu50inf.root"},
        Sample::dymumu
    );
    
    result["dytautau1050"] = Sample(
        "Z / #gamma* #rightarrow #tau#tau",
        kAzure+8,
        860.5,
        {"dytautau1050.root"},
        Sample::dytautau
    );
    
    result["dytautau50inf"] = Sample(
        "Z / #gamma* #rightarrow #tau#tau",
        kAzure+8,
        3532.8,
        {"dytautau50inf.root"},
        Sample::dytautau
    );
    
    result["wlnu"] = Sample(
        "W+Jets",
        kGreen-3,
        36257.2,
        {"wtolnu.root"}
    );
    
    result["qcdmu15"] = Sample(
        "QCD Multijet",
        kYellow,
        3.640E8*3.7E-4,
        {"qcdmu15.root"}
    );
    
    result["qcdmu2030"] = Sample(
        "QCD Multijet",
        kYellow,
        2.870E8*6.500E-3,
        {}
    );
    
    result["qcdmu3050"] = Sample(
        "QCD Multijet",
        kYellow,
        6.609E7*12.20E-3,
        {}
    );
    
    result["qcdmu5080"] = Sample(
        "QCD Multijet",
        kYellow,
        8.802E6*21.80E-3,
        {}
    );
    
    result["qcdmu80120"] = Sample(
        "QCD Multijet",
        kYellow,
        1.024E6*39.50E-3,
        {}
    );
    
    result["qcdmu120170"] = Sample(
        "QCD Multijet",
        kYellow,
        1.578E5*47.30E-3,
        {}
    );
    
    result["qcdem2030"] = Sample(
        "QCD Multijet",
        kYellow,
        2.886E8*10.10E-3,
        {"qcdem2030.root"}
    );
    
    result["qcdem3080"] = Sample(
        "QCD Multijet",
        kYellow,
        7.433E7*62.10E-3,
        {"qcdem3080.root"}
    );
    
    result["qcdem80170"] = Sample(
        "QCD Multijet",
        kYellow,
        1.191E6*153.9E-3,
        {"qcdem80170.root"}
    );
    
    result["qcdbcem2030"] = Sample(
        "QCD Multijet",
        kYellow,
        2.886E8*5.800E-4,
        {"qcdbcem2030.root"}
    );
    
    result["qcdbcem3080"] = Sample(
        "QCD Multijet",
        kYellow,
        7.424E7*2.250E-3,
        {"qcdbcem3080.root"}
    );
    
    result["qcdbcem80170"] = Sample(
        "QCD Multijet",
        kYellow,
        1.191E6*10.90E-3,
        {"qcdbcem80170.root"}
    );
    
    result["ttbarW"] = Sample(
        "t#bar{t}+Z/W/#gamma",
        kOrange-2,
        0.232,
        {"ttbarW.root"}
    );
    
    result["ttbarZ"] = Sample(
        "t#bar{t}+Z/W/#gamma",
        kOrange-2,
        0.2057,
        {"ttbarZ.root"},
        Sample::ttZ
    );
    
    result["ttgjets"] = Sample(
        "t#bar{t}+Z/W/#gamma",
        kOrange-2,
        1.8,
        {"ttgjets.root"}
    );
    
    return result;
}



std::vector<TString> SampleDefinitions::selectAndOrderSamples8TeV()
{
    const std::vector<TString> result = {
        "data",
        "ttbarW",
        "ttbarZ",
        "ttgjets",
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
        "ww",
        "wz",
        "zz",
        "singletop",
        "ttbarbkg",
        "ttbarsignal"
    };
    
    return result;
}





/// selectionRoot/Nominal/emu/emu_ttbarbgviatau.root
/// selectionRoot/Nominal/emu/emu_ttbarsignalplustau.root



