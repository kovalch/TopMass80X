#include "MassAnalyzer.h"

MassAnalyzer::MassAnalyzer(const TString& identifier, TTree* tree) : 
  fIdentifier(identifier), fTree(tree), fEntries(-10.),
  fMass(-10.), fMassError(-10.), //fMassSigma(-10.),
  fJES(-10.), fJESError(-10.),
  fMassConstJES(-10.), fMassConstJESError(-10.), //fMassConstJESSigma(-10.),
  ffSig(-10.), ffSigError(-10.),
  fMassfSig(-10.), fMassfSigError(-10.), //fMassfSigSigma(-10.),
  fJESfSig(-10.), fJESfSigError(-10.)  
{
}

double MassAnalyzer::GetMass() const {
  return fMass;
}

double MassAnalyzer::GetMassError() const{
  return fMassError;
}

//double MassAnalyzer::GetMassSigma() const{
//  return fMassSigma;
//}

double MassAnalyzer::GetJES() const{
  return fJES;
}

double MassAnalyzer::GetJESError() const {
  return fJESError;
}


double MassAnalyzer::GetMassConstJES() const {
  return fMassConstJES;
}

double MassAnalyzer::GetMassConstJESError() const{
  return fMassConstJESError;
}

//double MassAnalyzer::GetMassConstJESSigma() const{
//  return fMassConstJESSigma;
//}


double MassAnalyzer::GetFSig() const {
  return ffSig;
}

double MassAnalyzer::GetFSigError() const {
  return ffSigError;
}

double MassAnalyzer::GetMassfSig() const {
  return fMassfSig;
}

double MassAnalyzer::GetMassfSigError() const {
  return fMassfSigError;
}

//double MassAnalyzer::GetMassfSigSigma() const {
//  return fMassfSigSigma;
//}

double MassAnalyzer::GetJESfSig() const {
  return fJESfSig;
}

double MassAnalyzer::GetJESfSigError() const {
  return fJESfSigError;
}
 
