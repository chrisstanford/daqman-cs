#include "EvalRois.hh"
#include "BaselineFinder.hh"
#include "ConvertData.hh"
#include "SumChannels.hh"
#include "RootWriter.hh"
#include "intarray.hh"
#include "ChannelData.hh"
#include "Integrator.hh"
#include "Pulse.hh"
#include "EvalPulse.hh"
#include "EventHandler.hh"
#include <algorithm>
#include <numeric>

EvalRois::EvalRois() : 
  ChannelModule(GetDefaultName(),
		"Measure the max, min, and integral of samples over a set of regions of interest defined by start and end times in microseconds")
{
  AddDependency<ConvertData>();
  AddDependency<BaselineFinder>();
  RegisterParameter("regions", _regions, "Start/end time pairs to evaluate");
  RegisterParameter("max_peak_time", max_peak_time=-1, 
		    "Maximum time from the start where peak is expected");
}

EvalRois::~EvalRois() {}

int EvalRois::Initialize()
{
  if(_regions.size() > 0){
    Message(DEBUG)<<"Saving info for "<<_regions.size()<<" regions.\n";
  }
  else{
    Message(ERROR)<<"No ROIs defined for EvalRois!\n";
    return 1;
  }  
  return 0;
}

int EvalRois::Finalize() { return 0; }

int EvalRois::Process(ChannelData* chdata)
{
  
  if(!(chdata->baseline.found_baseline)) return 0;
  for(size_t index = 0; index < _regions.size(); index++){
    Pulse region;
    region.start_index = chdata->TimeToSample(_regions[index].first, true);
    region.end_index = chdata->TimeToSample(_regions[index].second, true);
    EvaluatePulse(region, chdata, max_peak_time);
    //    region.Print(chdata->channel_id, index);
    chdata->regions.push_back(region);
  }
  
  return 0;
}
     
