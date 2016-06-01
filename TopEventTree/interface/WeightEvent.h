/*
 * WeightEvent.h
 *
 *  Created on: Feb 7, 2013
 *      Author: eschliec
 */

#ifndef WEIGHTEVENT_H_
#define WEIGHTEVENT_H_

#include "TObject.h"
#include "TLorentzVector.h"
#include <vector>

class WeightEvent : public TObject {
 public:
  WeightEvent();
  void init();

  ClassDef(WeightEvent, 7);

  // WeightEvent data

  // MC weight
  double mcWeight;

  double mcGenWeight;

  // Generator weights (scale and pdf variations), check
  // https://twiki.cern.ch/twiki/bin/view/CMS/LHEReaderCMSSW for explanation
  std::vector<double> lheWeight;
  std::vector<unsigned int> lheWeightID;

  // MG BR correction
  double brWeight;

  // PU weight and control variables
  double puWeight;
  double puWeightUp;
  double puWeightDown;
  int nVertex;
  std::vector<int> nPU;
  double nPUTrue;

  // bTag efficiency and mistag rate weights
  std::vector<double> bTagEffJetWeight;
  std::vector<double> lTagEffJetWeight;
  double bTagWeight;
  double bTagWeight_bTagSFUp;
  double bTagWeight_bTagSFDown;
  double bTagWeight_misTagSFUp;
  double bTagWeight_misTagSFDown;

  // bJES neutrino fraction and fragmentation weights
  double bJESWeight_fNuUp;
  double bJESWeight_fNuDown;
  double bJESWeight_frag;
  double bJESWeight_fragHard;
  double bJESWeight_fragSoft;

  // trigger weight
  double triggerWeight;

  // lepton weights
  double muWeight;
  double elWeight;

  // PDF weights and variables needed for a recalculation of PDF weights
  std::vector<double> pdfWeight;
  double x1, x2;
  float Q;
  int id1, id2;

  // ME level top quarks
  TLorentzVector meTop1;
  TLorentzVector meTop2;

  // default combined event weight
  double combinedWeight;

  // for top pt systematics, not normalized!!!
  double topPtSysWeight;

 private:
};

#endif /* WEIGHTEVENT_H_ */
