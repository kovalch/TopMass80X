#include "DataSample.h"


void DataSample::Clear()
{
  nEvents = 0;
  maxWeight = 0;

  events.clear();
}

void DataSample::Fill(double topMass, double wMass, double prob, int leptonFlavour, double weight, int index, int bin)
{
  if(index == 0){
    events.push_back(SimpleEvent());
    events.back().leptonFlavour = leptonFlavour;
    events.back().weight = weight;
    if(weight > maxWeight) maxWeight = weight;
    if(weight < minWeight) minWeight = weight;
    ++nEvents;
  }
  events.back().permutations.push_back({topMass, wMass, prob, bin});
}

void DataSample::AddEvent(const SimpleEvent& event)
{
  events.push_back(event);
  if(event.weight > maxWeight) maxWeight = event.weight;
  if(event.weight < minWeight) minWeight = event.weight;
  ++nEvents;
}

DataSample& DataSample::operator+=(const DataSample& sample)
{
  nEvents += sample.nEvents;
  if(sample.maxWeight > maxWeight) maxWeight = sample.maxWeight;
  if(sample.minWeight < minWeight) minWeight = sample.minWeight;

  events.insert(events.end(),sample.events.begin(),sample.events.end());

  return *this;
}
