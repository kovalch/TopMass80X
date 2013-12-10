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

  ClassDef(WeightEvent,2);

  // WeightEvent data

  // MC weight
  // currently only used for MC@NLO with 1 or -1
  double mcWeight;
  
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

  // default combined event weight
  double combinedWeight;


private:
};

#endif /* WEIGHTEVENT_H_ */
