#include "SumChannels.hh"
#include "ConvertData.hh"
#include "BaselineFinder.hh"
#include "EventHandler.hh"
#include <algorithm>
#include "RootWriter.hh"

SumChannels::SumChannels() :
  BaseModule(GetDefaultName(),"Create a virtual channel whose waveform is the sum of all other channels in the event")
{
  AddDependency<ConvertData>();
  AddDependency<BaselineFinder>();
}

SumChannels::~SumChannels()
{
  Finalize();
}

int SumChannels::Initialize()
{
  return 0;
}

int SumChannels::Finalize()
{
  return 0;
}

class add_with_scaling{
  double _scale_factor;
public:
  add_with_scaling(double scale_factor) : _scale_factor(scale_factor) {}
  double operator()(double x, double y) { return x + y*_scale_factor; }
};

int SumChannels::Process(EventPtr event)
{
  ChannelData sumdata;
  sumdata.channel_id = ChannelData::CH_SUM;
  sumdata.label = "sum";
  sumdata.sample_bits = 32;
  int n_channels_summed = 0;
  EventDataPtr data = event->GetEventData();
  if (data->channels.size() < 2)
    //No point in summing channels
    return 0;

  double sum_bl_mean=0, sum_bl_sigma=0;  
  for(size_t i=0; i<data->channels.size(); i++){
    const ChannelData& chdata = data->channels[i];
    if( _skip_channels.find( chdata.channel_id) != _skip_channels.end() ||
	chdata.channel_id < 0 || !(chdata.baseline.found_baseline) ){
      continue;
    }
    sum_bl_mean += chdata.baseline.mean/chdata.spe_mean;
    sum_bl_sigma += chdata.baseline.sigma/chdata.spe_mean;
    if(n_channels_summed == 0){
      sumdata.sample_rate = chdata.sample_rate;
      sumdata.trigger_index = chdata.trigger_index;
      sumdata.nsamps = chdata.nsamps;
      sumdata.waveform.resize(chdata.nsamps);
    }
    
    //line up the waveforms of the channels
    const int presamps = std::min(sumdata.trigger_index, chdata.trigger_index);
    const int postsamps = std::min(sumdata.nsamps - sumdata.trigger_index,
				   chdata.nsamps - chdata.trigger_index);
    //load the scale factor from the calibration database
    double scale_factor = 1./chdata.spe_mean;;
    
    //only sum the subtracted_waveform 
    //sum the two channels into a temporary vector
    std::vector<double> tempsum(presamps+postsamps);
    std::vector<double>::const_iterator sumit = sumdata.subtracted_waveform.begin();
    std::vector<double>::const_iterator chit = chdata.subtracted_waveform.begin();
    //sum from trigger_index - presamps to trigger_index+postsamps
    std::transform(sumit+sumdata.trigger_index-presamps,
		   sumit+sumdata.trigger_index+postsamps, 
		   chit+chdata.trigger_index-presamps,
		   tempsum.begin(),
		   add_with_scaling(scale_factor));
    //copy the temporary into the sumdata waveform vector
    sumdata.subtracted_waveform.assign(tempsum.begin(), tempsum.end());
    sumdata.trigger_index = presamps;
    sumdata.nsamps = presamps+postsamps;

    data->channels_summed.insert(chdata.channel_id);
    n_channels_summed++;
    // the sum is "saturated" if any single channel is
    if(chdata.saturated) sumdata.saturated = true;
  } //end loop over channels

  if (n_channels_summed > 0){
    sumdata.waveform.assign(sumdata.subtracted_waveform.begin(), sumdata.subtracted_waveform.end());
    for(int samp=0; samp<sumdata.nsamps; samp++) sumdata.waveform[samp] += sum_bl_mean;
    sumdata.baseline.mean = sum_bl_mean;
    sumdata.baseline.sigma = sum_bl_sigma/n_channels_summed;
    //reset the historical channel_start and end pointers
    sumdata.channel_start = (char*)(sumdata.GetWaveform());
    sumdata.channel_end = (char*)( sumdata.GetWaveform() + sumdata.nsamps );

    //find the max and min of the channel
    double* wave = sumdata.GetWaveform();
    double* max_samp = std::max_element(wave, wave+sumdata.nsamps);
    double* min_samp = std::min_element(wave, wave+sumdata.nsamps);
    //data is saturated if it hit 0 or maximum range
    sumdata.maximum = *max_samp;
    sumdata.minimum = *min_samp;
    sumdata.max_time = sumdata.SampleToTime(max_samp - wave);
    sumdata.min_time = sumdata.SampleToTime(min_samp - wave);
    data->channels.push_back(sumdata);
  }
  return 0;
  
}

