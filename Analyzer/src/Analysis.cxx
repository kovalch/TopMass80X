#include "Analysis.h"

#include <string>
#include <stdlib.h>

#include "TSystem.h"
#include "TEventList.h"

#include "XMLConfigReader.h"
#include "Helper.h"
#include "GenMatchAnalyzer.h"
#include "MVAAnalyzer.h"
#include "IdeogramAnalyzer.h"
#include "RooFitTemplateAnalyzer.h"

#include "RandomSubsetCreatorAllJets.h"

#include "LHAPDF/LHAPDF.h"

typedef XMLConfigReader xml;

Analysis::Analysis(po::variables_map vm) :
  _fIdentifier(vm["input" ].as<std::string>()),
  _fMethod    (vm["method"].as<std::string>()),
  _fBins      (vm["bins"  ].as<int>()),
  _fChannel("alljets"),
  _fTree(0)
{
  CreateHistos();

  _fChannel = xml::GetParameter("decayChannel");
}

Analysis::Analysis(TString identifier, TString file, TString method, int bins, double lumi) :
  _fIdentifier(identifier),
  _fMethod(method), _fBins(bins),
  _fTree(0)
{
  CreateHistos();

  _fChannel = xml::GetParameter("decayChannel");
}

Analysis::~Analysis()
{
  delete _fTree;
  delete _hEntries;
  delete _hMass;
  delete _hMassError;
  delete _hJES;
  delete _hJESError;
  delete _hMassConstJES;
  delete _hMassConstJESError;
  delete _hFSig;
  delete _hFSigError;
  delete _hMassfSig;
  delete _hMassfSigError;
  delete _hJESfSig;
  delete _hJESfSigError;
}

void Analysis::Analyze(po::variables_map vm) {

  std::cout << "Analyze " << _fIdentifier << " with method " << _fMethod << std::endl;

  // random subset creation
  RandomSubsetCreator* fCreator = 0;
  if (!strcmp(_fChannel, "AllJets")) {
    fCreator = new RandomSubsetCreatorAllJets(vm);
  }
  _fTree = fCreator->CreateRandomSubset();

  MassAnalyzer* fAnalyzer = 0;

  if (!strcmp(_fMethod, "GenMatch")) {
    fAnalyzer = new GenMatchAnalyzer(_fIdentifier, _fTree);
  }
  else if (!strcmp(_fMethod, "MVA")) {
    fAnalyzer = new MVAAnalyzer(_fIdentifier, _fTree);
  }
  else if (!strcmp(_fMethod, "Ideogram")) {
    fAnalyzer = new IdeogramAnalyzer(_fIdentifier, _fTree);
  }
  else if (!strcmp(_fMethod, "RooFit")) {
    fAnalyzer = new RooFitTemplateAnalyzer(_fIdentifier, _fTree);
  }
  else {
    return;
  }

  Helper* helper = new Helper(_fBins);
  helper->SetTDRStyle();
  delete helper;

  TCanvas* canvas = new TCanvas("canvas", "Top mass", 900, 600);
  canvas->cd();

  double smearBins = 1.;
  double rangeX = 10.;
  double rangeY = 10.;

  double smearX = smearBins/_fBins*rangeX;
  double smearY = smearBins/_fBins*rangeY;

  double minEntries = 25;

  TString observableX = "dRbb";
  TString observableY = "dRbb";

  for(int i = 0; i < _fBins; i++) {
    for(int j = 0; j < _fBins; j++) {
      // calculate cuts
      std::stringstream stream;
      stream << (rangeX/_fBins)*(i)-smearX << "<" << observableX << "&"
             << observableX << "<" << (rangeX/_fBins)*(i+1)+smearX << " & "
             << rangeY/_fBins*(j)-smearY << "<" << observableY << "&"
             << observableY << "<" << rangeY/_fBins*(j+1)+smearY;
      TString cuts(stream.str());

      if (!strcmp(_fMethod, "GenMatch")) {
        cuts += " & comboType == 1";
      }
      else if (!strcmp(_fMethod, "MVA")) {
        cuts += " & mvaDisc > 0";
      }
      else if (!strcmp(_fMethod, "Ideogram")) {
//        //cuts += " & target == 1";
//        //cuts += " & (target == 0 | target == 1)";
//        //cuts += " & Njet == 6";
//        cuts += " & probs[0] > 0.09";
//        cuts += " & dRbb > 1.5";
//        //cuts += " & nVertex >= 5";
//        //cuts += " & nVertex <= 5";
//        cuts += " & topMasses[0] > 100. & topMasses[0] < 550.";
      }
      else if (!strcmp(_fMethod, "RooFit")) {
      }

      int entries = _fTree->GetEntries(cuts);
      std::cout << cuts << std::endl;
      std::cout << entries << std::endl;

      _hEntries->SetCellContent(i+1, j+1, entries);

      if (entries > minEntries) {
        fAnalyzer->Analyze(cuts, i, j, vm);

        _hMass     ->SetCellContent(i+1, j+1, fAnalyzer->GetMass());
        _hMassError->SetCellContent(i+1, j+1, fAnalyzer->GetMassError());
        _hJES      ->SetCellContent(i+1, j+1, fAnalyzer->GetJES());
        _hJESError ->SetCellContent(i+1, j+1, fAnalyzer->GetJESError());
        _hMassConstJES     ->SetCellContent(i+1, j+1, fAnalyzer->GetMassConstJES());
        _hMassConstJESError->SetCellContent(i+1, j+1, fAnalyzer->GetMassConstJESError());
        _hFSig     ->SetCellContent(i+1, j+1, fAnalyzer->GetFSig());
        _hFSigError->SetCellContent(i+1, j+1, fAnalyzer->GetFSigError());
        _hMassfSig     ->SetCellContent(i+1, j+1, fAnalyzer->GetMassfSig());
        _hMassfSigError->SetCellContent(i+1, j+1, fAnalyzer->GetMassfSigError());
        _hJESfSig      ->SetCellContent(i+1, j+1, fAnalyzer->GetJESfSig());
        _hJESfSigError ->SetCellContent(i+1, j+1, fAnalyzer->GetJESfSigError());

        std::cout << std::endl;
        std::cout << "Measured 1D mass: " << fAnalyzer->GetMassConstJES() << " +/- "
                  << fAnalyzer->GetMassConstJESError() << " GeV" << std::endl;
        std::cout << std::endl;
        std::cout << "Measured 2D mass: " << fAnalyzer->GetMass() << " +/- "
                  << fAnalyzer->GetMassError() << " GeV" << std::endl;
        std::cout << "Measured 2D JES: " << fAnalyzer->GetJES() << " +/- "
                  << fAnalyzer->GetJESError() << " " << std::endl;
        std::cout << std::endl;
        std::cout << "Measured 3D mass: " << fAnalyzer->GetMassfSig() << " +/- "
                  << fAnalyzer->GetMassfSigError() << " GeV" << std::endl;
        std::cout << "Measured 3D JES: " << fAnalyzer->GetJESfSig() << " +/- "
                  << fAnalyzer->GetJESfSigError() << " " << std::endl;
        std::cout << "Measured 3D fSig: " << fAnalyzer->GetFSig() << " +/- "
                  << fAnalyzer->GetFSigError() << " GeV" << std::endl;
        std::cout << std::endl;
      }
    }
  }

  canvas->Clear();
  canvas->Divide(2,2);

  canvas->cd(1);
  _hEntries->Draw("COLZ");

  canvas->cd(2);
  _hMass->Draw("COLZ,TEXT");
  _hMass->SetAxisRange(_hMass->GetMinimum(0.05), _hMass->GetMaximum(), "Z");

  canvas->cd(3);
  _hMassError->Draw("COLZ,TEXT");
  _hMassError->SetAxisRange(0.05, 5, "Z");

  TString path("plot/"); path += _fMethod; path += "_"; path += _fIdentifier; path += ".eps";
  canvas->Print(path);

  delete canvas;
  delete fAnalyzer;
  delete fCreator;
}

void Analysis::CreateHistos() {
  Helper* helper = new Helper(_fBins);

  _hEntries = helper->GetH2("Entries");

  _hMass = helper->GetH2("Mass");
  _hMassError = helper->GetH2("MassError");
  _hJES = helper->GetH2("JES");
  _hJESError = helper->GetH2("JESError");
  
  _hMassConstJES = helper->GetH2("MassConstJES");
  _hMassConstJESError = helper->GetH2("MassConstJESError");
  
  _hFSig = helper->GetH2("fSig");
  _hFSigError = helper->GetH2("fSigError");
  _hMassfSig = helper->GetH2("MassfSig");
  _hMassfSigError = helper->GetH2("MassfSigError");
  _hJESfSig = helper->GetH2("JESfSig");
  _hJESfSigError = helper->GetH2("JESfSigError");

  delete helper;
}

TH2F* Analysis::GetH2Mass() {
  return _hMass;
}

TH2F* Analysis::GetH2MassError() {
  return _hMassError;
}

TH2F* Analysis::GetH2JES() {
  return _hJES;
}

TH2F* Analysis::GetH2JESError() {
  return _hJESError;
}

TH2F* Analysis::GetH2MassConstJES() {
  return _hMassConstJES;
}

TH2F* Analysis::GetH2MassConstJESError() {
  return _hMassConstJESError;
}

TH2F* Analysis::GetH2FSig() {
  return _hFSig;
}

TH2F* Analysis::GetH2FSigError() {
  return _hFSigError;
}

TH2F* Analysis::GetH2MassfSig() {
  return _hMassfSig;
}

TH2F* Analysis::GetH2MassfSigError() {
  return _hMassfSigError;
}

TH2F* Analysis::GetH2JESfSig() {
  return _hJESfSig;
}

TH2F* Analysis::GetH2JESfSigError() {
  return _hJESfSigError;
}

TString Analysis::GetIdentifier() {
  return _fIdentifier;
}
