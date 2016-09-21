#include "Event.hh"
#include "PulseEvaluator.hh"
#include "ConvertData.hh"
#include "BaselineFinder.hh"
#include "Integrator.hh"
#include "PulseUtility.hh"
#include "Spike.hh"
#include "intarray.hh"
#include "RootWriter.hh"
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
//#include <functional>
#include <algorithm>

PulseEvaluator::PulseEvaluator():
  ChannelModule(GetDefaultName(), "Module to evaluate pulse/ROI quantities")
{
  AddDependency<ConvertData>();
  AddDependency<BaselineFinder>();
  AddDependency<Integrator>();
  
  //Register all the config handler parameters
  RegisterParameter("evaluate_pulses", evaluate_pulses = true,
		    "whether pulse quantities should be evaluated");
  RegisterParameter("evaluate_rois", evaluate_rois = true,
		    "whether roi quantities should be evaluated");
  RegisterParameter("max_peak_time_us", max_peak_time_us = -1,
		    "maximum peak time in us from the beginning of the pulse");
//   RegisterParameter("",  = 10,
// 		    "");

}

PulseEvaluator::~PulseEvaluator()
{;}

int PulseEvaluator::Initialize()
{
  //no point coming to this module if we don't evaluate anything
  if(!evaluate_pulses && !evaluate_rois) return -1;
  return 0;
}

int PulseEvaluator::Finalize()
{   
  return 0;
}

int PulseEvaluator::Process(ChannelData* chdata){

  if(evaluate_rois){  
    for (size_t i=0; i<chdata->regions.size(); i++){
      Pulse & pulse = chdata->regions.at(i);
      EvaluatePulse(pulse, chdata, max_peak_time_us);
    } // end for loop over pulses
  }//end if evaluate rios

  if(evaluate_pulses){  
    for (size_t i=0; i<chdata->pulses.size(); i++){
      Pulse & pulse = chdata->pulses.at(i);
      EvaluatePulse(pulse, chdata, max_peak_time_us);
      //evaluate the pulse parameters depending on other pulses
      if(i>0){
	pulse.time_since_last=pulse.half_max_time-chdata->pulses.at(i-1).end_time;
	if(pulse.start_index<=chdata->pulses.at(i-1).end_index)
	  {pulse.is_clean = false; chdata->pulses.at(i-1).is_clean = false;}
      }
      else{
	pulse.time_since_last=pulse.half_max_time-chdata->SampleToTime(0);
	if(pulse.start_index==0) pulse.is_clean = false;
      }
      
    } // end for loop over pulses
  }//end if evaluate pulses
  /*
  //get rid of noise pulses
  //	double noise_start=0, noise_end=0;
  for(size_t j=0; j<chdata->pulses.size(); j++){
  Pulse &pulse = chdata->pulses.at(j);

  if(pulse.integral+pulse.peak_amplitude>=0 ||
  //	     pulse.overshoot>5*chdata->baseline.sigma ||
  pulse.overshoot>0.5*pulse.peak_amplitude){
  pulse.evaluated=false;
  / *	    noise_start = pulse.start_time-0.2;
  noise_end = pulse.end_time+0.2;
  //get rid of the previous noise pulses
  for(int k=j-1; k>=0; k--){
  if(chdata->pulses.at(k).start_time>=noise_start &&
  chdata->pulses.at(k).start_time<=noise_end) chdata->pulses.at(k).evaluated=false;
  else break;
  }* /
  }//end if this is a bad pulse
  / *	  if(noise_end>noise_start){//if this pulse comes close after a noise
  if(chdata->pulses.at(j).start_time>=noise_start &&
  chdata->pulses.at(j).start_time<=noise_end) chdata->pulses.at(j).evaluated=false;
  else {noise_start=0; noise_end=0;}//out of noise window
  }//end if noise end* /
  }//end for loop
  */

  return 0;
}

//function to calculate pulse parameters, if max_peak_time>0
//the peak of the pulse should be within (start time, start time + max_peak_time)
int PulseEvaluator::EvaluatePulse(Pulse& pulse, ChannelData* chdata, double max_peak_time){
  
  if(!chdata->baseline.found_baseline ||  chdata->integral.empty() ||
     pulse.start_index<0 || pulse.end_index>=chdata->nsamps)
    return 1;
  //start time and end time
  pulse.start_time = chdata->SampleToTime(pulse.start_index);
  pulse.end_time = chdata->SampleToTime(pulse.end_index);
//   if(pulse.peak_index>0)
//     pulse.peak_time_raw = chdata->SampleToTime(pulse.peak_index);
  pulse.integral = chdata->integral[pulse.end_index] - 
    chdata->integral[pulse.start_index];
  pulse.npe = pulse.integral/chdata->spe_mean;

  const double* subtracted = chdata->GetBaselineSubtractedWaveform();
  //find peak in the time window of (start_time, start_time+max_peak_time)
  int max_peak_index = pulse.end_index;
  if(max_peak_time>0) 
    max_peak_index = chdata->TimeToSample(pulse.start_time+max_peak_time);
  if(max_peak_index>pulse.end_index) max_peak_index = pulse.end_index;
  //find the peak
  int max_index = std::max_element(subtracted + pulse.start_index,
				   subtracted + max_peak_index) - subtracted;
  int min_index = std::min_element(subtracted + pulse.start_index,
				   subtracted + max_peak_index) - subtracted;
  //warning: this will not work for bipolar pulses!
  //use pulse integral to determine if this is a positive or negative pulse
  if(pulse.integral<0){
    pulse.peak_index = min_index;
    pulse.peak_amplitude = -subtracted[min_index]; //pulse are neg, amplitude pos
    pulse.overshoot = subtracted[max_index];
  }
  else{
    pulse.peak_index = max_index;
    pulse.peak_amplitude = subtracted[max_index];
    pulse.overshoot = subtracted[min_index];
  }
  pulse.peak_time = chdata->SampleToTime(pulse.peak_index);

  //to look for 50% threshold time
  int index=-1;
  double threshold = 0.5*subtracted[pulse.peak_index];
  if(std::fabs(threshold)<=chdata->baseline.sigma) return 0;
  //if this is a small pulse, search backwards to be safe
  if(std::fabs(threshold)<4.*chdata->baseline.sigma){
    for(index=pulse.peak_index-1; index>=pulse.start_index; index--){
      if(FirstAmplitudeIsSmaller(subtracted[index],threshold,threshold)) break;
    }//end for index
  }//end if threshold
  else{//the pulse is not too small, search forward
    for(index=pulse.start_index+1; index<=pulse.peak_index; index++){
      if(FirstAmplitudeIsSmaller(threshold,subtracted[index],threshold)){
	index --;
	break;
      }//end if
    }//end for index
  }//end else
  
  //  std::cout<<"Index: "<<index<<'\t'<<chdata->SampleToTime(index)<<std::endl;
  //linear interpolation to calculate the time
  if(index>=pulse.start_index && index<pulse.peak_index && 
     FirstAmplitudeIsSmaller(subtracted[index], subtracted[index+1], threshold)){
    pulse.half_max_time = chdata->SampleToTime(index)+
      (subtracted[index]-threshold)/(subtracted[index]-subtracted[index+1])/chdata->sample_rate;
    // std::cerr<<"Pulse arriving Time:  "<<pulse.half_max_time
    //  	     <<" of channel "<<chdata->channel_id<<" in event at time "<<chdata->timestamp<<std::endl;
  }
  else{
    pulse.half_max_time = pulse.start_time;
    // std::cerr<<"Pulse arriving Time can not be evaluated for Pulse starting at time "<<pulse.start_time
    //  	     <<" of channel "<<chdata->channel_id<<" in event at time "<<chdata->timestamp<<std::endl;
  }
  pulse.is_clean = (std::fabs(subtracted[pulse.start_index])<std::fabs(threshold) && 
		    std::fabs(subtracted[pulse.end_index])  <std::fabs(threshold));
  
  //Check to see if peak is saturated
  const double* wave = chdata->GetWaveform();
  max_index = std::max_element(wave + pulse.start_index, 
			       wave + pulse.end_index) - wave;
  min_index = std::min_element(wave + pulse.start_index, 
			       wave + pulse.end_index) - wave;
  pulse.max = wave[max_index];
  pulse.min = wave[min_index];

  if(wave[max_index] >= chdata->GetVerticalRange()){
    pulse.saturated = true;
    if(pulse.integral>0){
      int max_end_index = max_index + 1;
      while (wave[max_end_index] >= chdata->GetVerticalRange() && 
	     max_end_index < pulse.end_index){
	max_end_index++;
      }
      pulse.peak_index = (int)(max_index + max_end_index)/2;
    }
  }
  if(wave[min_index] <= 0){
    pulse.saturated = true;
    if(pulse.integral<0){
      int min_end_index = min_index + 1;
      while (wave[min_end_index] == 0 && min_end_index < pulse.end_index){
	min_end_index++;
      }
      pulse.peak_index = (int)(min_index + min_end_index)/2;
    }
  }

  //evaluate the spike counting
  if(chdata->channel_id == ChannelData::CH_SUM){
    //sum channel is processed after all other channels
    std::set<int> &channels_summed = _current_event->GetEventData()->channels_summed;
    int nspikes=0;
    for(std::set<int>::iterator it=channels_summed.begin(); it!=channels_summed.end(); it++){
      ChannelData * ch = _current_event->GetEventData()->GetChannelByID(*it);
      bool pulse_identified=false;
      //first try the pulses
      for(size_t ii=0; ii<ch->pulses.size(); ii++){
	if(ch->pulses.at(ii).start_index==pulse.start_index &&
	   ch->pulses.at(ii).end_index==pulse.end_index){
	  pulse_identified=true;
	  nspikes += ch->pulses.at(ii).nspikes;
	  break;
	}//end if
      }//end for ii
      if(pulse_identified) continue;
      //try the regions
      for(size_t ii=0; ii<ch->regions.size(); ii++){
	if(ch->regions.at(ii).start_index==pulse.start_index &&
	   ch->regions.at(ii).end_index==pulse.end_index){
	  pulse_identified=true;
	  nspikes += ch->regions.at(ii).nspikes;
	  break;
	}//end if
      }//end for ii
      //here we can try other pulse objects if any
    }//end for it
    pulse.nspikes = nspikes;
  }
  //here it is not the sum channel
  else if(chdata->spikes.size()>0){
    int nspikes = 0;
    std::vector<Spike> &spikes=chdata->spikes;
    for(size_t ii=0; ii<spikes.size(); ii++){
      if(spikes.at(ii).peak_time>=pulse.start_time && spikes.at(ii).peak_time<=pulse.end_time)
	nspikes ++;
    }
    pulse.nspikes = nspikes;
  }//end else if

  pulse.evaluated = true;
  //  pulse.Print(chdata->channel_id, 0);
  return 0;
}
