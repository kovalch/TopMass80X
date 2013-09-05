#include "DataSample.h"

void DataSample::Clear()
{
  nEvents = 0;
  maxWeight = 0;

  topMasses.clear();
  wMasses  .clear();
  fitProbs .clear();
  weights  .clear();
  indices  .clear();
}

void DataSample::Fill(double topMass, double wMass, double prob, double weight, unsigned char index)
{
  topMasses.push_back(topMass);
  wMasses  .push_back(wMass);
  fitProbs .push_back(prob);
  weights  .push_back(weight);
  indices  .push_back(index);
  if(weight > maxWeight) maxWeight = weight;
}

DataSample& DataSample::operator+=(const DataSample& sample)
{
  nEvents += sample.nEvents;
  if(sample.maxWeight > maxWeight) maxWeight = sample.maxWeight;

  topMasses.insert(topMasses.end(),sample.topMasses.begin(),sample.topMasses.end());
  wMasses  .insert(wMasses  .end(),sample.wMasses  .begin(),sample.wMasses  .end());
  fitProbs .insert(fitProbs .end(),sample.fitProbs .begin(),sample.fitProbs .end());
  weights  .insert(weights  .end(),sample.weights  .begin(),sample.weights  .end());
  indices  .insert(indices  .end(),sample.indices  .begin(),sample.indices  .end());

  return *this;
}
