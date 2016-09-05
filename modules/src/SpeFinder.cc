#include "SpeFinder.hh"
//#include "EventHandler.hh"
#include "ConvertData.hh"
#include "PulseUtility.hh"
#include "BaselineFinder.hh"
#include <vector>
#include <numeric>

SpeFinder::SpeFinder():
  ChannelModule(GetDefaultName(), "Search for single photoelectrons in the tails of pulses identified by PulseFinder")
{
  AddDependency<BaselineFinder>();
  RegisterParameter("search_start_time", search_start_time = 1.,
                    "Time in [us] to start searching for spe");
  RegisterParameter("search_end_time", search_end_time = 50.,
                    "Time in [us] to end searching for spe");
  RegisterParameter("ch_spe_thresholds", ch_spe_thresholds,
		    "discriminator threshold values for spe finding");
  RegisterParameter("relative_bls_threshold", relative_bls_threshold = false,
		    "Is the threshold relative to baseline sigma");
  RegisterParameter("relative_spe_threshold", relative_spe_threshold = false,
		    "Is the threshold relative to spe size");
  RegisterParameter("spe_edge_add_nsamps", spe_edge_add_nsamps  = 3,
                    "Number of samples to be added at the edge of possible pulses");
  RegisterParameter("min_baseline_length_us", min_baseline_length_us = 0.04,
                    "Time in us of good baseline before and after pulse");
}

SpeFinder::~SpeFinder()
{
  Finalize();
}

int SpeFinder::Initialize()
{
  //make sure we don't process the sum channel
  _skip_sum = true;
  if(search_start_time>=search_end_time || spe_edge_add_nsamps<1 ||
     (relative_bls_threshold && relative_spe_threshold) ||
     min_baseline_length_us<0 ) return -1;
  return 0;
}

int SpeFinder::Finalize(){
  return 0;
}

int SpeFinder::Process(ChannelData* chdata){

  if(chdata->channel_id==ChannelData::CH_SUM || 
     !chdata->baseline.found_baseline) return 0;
	
  double threshold = 0;
  id_map::iterator it = ch_spe_thresholds.find(chdata->channel_id);
  if(it != ch_spe_thresholds.end()){
    if(it->second!=0) threshold = it->second;
    else return 0;
  }//end if
  else{
    it = ch_spe_thresholds.find(ChannelData::CH_DEFAULT);
    if(it != ch_spe_thresholds.end() && it->second!=0)
      threshold =it->second;
    else return 0; //we simply don't process this channel
  }//end else
  
  if(relative_bls_threshold){
    if(chdata->baseline.sigma>0) threshold *= chdata->baseline.sigma;
    else{
      Message(ERROR)<<"Channel "<<chdata->channel_id<< " has invalid baseline!\n";
      return 0;
    }
  }

  if(relative_spe_threshold){
    if(chdata->spe_mean>0){
      threshold *= chdata->spe_mean;
      if(chdata->spe_mean==1){
	Message(WARNING)<<"Channel "<<chdata->channel_id<< " uses spe_mean unspecified!\n";
      }//end if ==1
    }//end if >0
    else{
      Message(ERROR)<<"Channel "<<chdata->channel_id<< " has incorrectly specified spe_mean!\n";
      return -1;
    }//end else
  }//end relative spe

  int search_start_index  = chdata->TimeToSample(search_start_time)+spe_edge_add_nsamps;
  int search_end_index    = chdata->TimeToSample(search_end_time)-spe_edge_add_nsamps;
  int min_baseline_nsamps = min_baseline_length_us*chdata->sample_rate+1;

  //use the subtracted waveform to search for spe start and end
  double* subtracted = chdata->GetBaselineSubtractedWaveform();
  std::vector<int> start_index;
  std::vector<int> end_index;
  
  int first_baseline_samp = -1, last_baseline_samp = -1;
  bool in_spe = false;
  for(int index = search_start_index; index<search_end_index; index++){
    if(FirstAmplitudeIsSmaller(threshold, subtracted[index], threshold)){
      //index-1 is safe because we required spe_edge_add_nsamps>=1
      if(!FirstAmplitudeIsSmaller(threshold,subtracted[index-1],threshold)){//if just come to a spe
        in_spe = true;
	start_index.push_back(index-spe_edge_add_nsamps );
      }
      if(!FirstAmplitudeIsSmaller(threshold,subtracted[index+1],threshold)){//if about to leave spe
        if(in_spe) end_index.push_back(index+spe_edge_add_nsamps );
        in_spe = false;
      }
    }//end if subtracted index <=
    else{//if not in a spe
      if(first_baseline_samp<0) first_baseline_samp = index;
      last_baseline_samp = index;
    }
  }
  if(in_spe) start_index.pop_back();
  if(start_index.size()!=end_index.size() || start_index.size()==0
     || first_baseline_samp<0 || last_baseline_samp<0){
    start_index.clear();
    end_index.clear();
    return 0;
  }

  //evaluate all the pulses, use raw waveform
  double* wave = chdata->GetWaveform();
  int pre_pulse_end = -1, post_pulse_start = -1, peak_index=-1;
  double bl_mean = 0, bl_sigma = 0, spe_integral=0, peak_amplitude=-1;
  for(size_t i=0; i<start_index.size(); i++){
    if(i>0) pre_pulse_end = end_index.at(i-1);
    else pre_pulse_end = first_baseline_samp;
    if(i<start_index.size()-1) post_pulse_start = start_index.at(i+1);
    else post_pulse_start = last_baseline_samp;
    //check the good baseline length
    if(start_index.at(i)-pre_pulse_end<min_baseline_nsamps 
       || post_pulse_start-end_index.at(i)<min_baseline_nsamps){
      continue;
    }
    //calculate the local baseline mean and sigma
    bl_mean = 0;
    bl_sigma = 0;
    spe_integral = 0;
    peak_amplitude=chdata->baseline.mean;
    for(int index=start_index.at(i)-min_baseline_nsamps;
	index<=end_index.at(i)+min_baseline_nsamps; index++){
      if(index>=start_index.at(i) && index<=end_index.at(i)){
	spe_integral += wave[index];
	if(FirstAmplitudeIsSmaller(peak_amplitude,wave[index], threshold)){
	  peak_amplitude = wave[index];
	  peak_index = index;
	}
      }
      else{
	bl_mean += wave[index];
	bl_sigma += wave[index]*wave[index];
      }//end else
    }//end for index
    bl_mean /= 2.*min_baseline_nsamps;
    bl_sigma /= 2.*min_baseline_nsamps;
    bl_sigma = std::sqrt(bl_sigma-bl_mean*bl_mean);
    spe_integral -= bl_mean*(end_index.at(i)-start_index.at(i)+1);
    peak_amplitude = std::fabs(bl_mean - peak_amplitude);

    //store spe
    Spe aspe;
    aspe.start_time = chdata->SampleToTime(start_index.at(i));
    aspe.peak_time = chdata->SampleToTime(peak_index);
    aspe.width = (end_index.at(i)-start_index.at(i)+1.0)/chdata->sample_rate;
    aspe.peak_amplitude = peak_amplitude;
    aspe.integral = spe_integral;
    aspe.baseline_mean = bl_mean;
    aspe.baseline_sigma = bl_sigma;
    chdata->single_pe.push_back(aspe);
    
  }//end for i

  return 0;
}
