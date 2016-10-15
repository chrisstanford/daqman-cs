#include "ChannelData.hh"
#include <algorithm>
#include <numeric>
#include <cmath>

int EvaluatePulse(Pulse& pulse, ChannelData* chdata, double max_peak_time){
  
  if(!chdata->baseline.found_baseline || pulse.start_index<0 ||
     pulse.end_index>=chdata->nsamps)
    return 1;
  //start time and end time
  pulse.start_time = chdata->SampleToTime(pulse.start_index);
  pulse.end_time = chdata->SampleToTime(pulse.end_index);

  const double* subtracted = chdata->GetBaselineSubtractedWaveform();
  //find peak in the time window of (start_time, start_time+max_peak_time)
  int max_peak_index = pulse.end_index;
  if(max_peak_time>0) 
    max_peak_index = chdata->TimeToSample(pulse.start_time+max_peak_time);
  if(max_peak_index>pulse.end_index) max_peak_index = pulse.end_index;
  //find the peak
  int peak_index = std::min_element(subtracted + pulse.start_index, 
				    subtracted + max_peak_index) - subtracted;
  //warning this is to test the noise rejection algorithm
  pulse.overshoot = *std::max_element(subtracted + pulse.start_index, 
				      subtracted + max_peak_index);
  pulse.peak_index = peak_index;
  pulse.peak_time = chdata->SampleToTime(pulse.peak_index);
  pulse.peak_amplitude = -subtracted[peak_index]; //pulse are neg, amplitude pos
  
  //to look for 50% threshold time
  int index=-1;
  //  double threshold = 0.5*subtracted[peak_index];
  double threshold = 0.1*subtracted[peak_index];
  //if this is a small pulse, search backwards to be safe
  if(threshold>-4){
    for(index=peak_index-1; index>=pulse.start_index; index--){
      if(subtracted[index]>=threshold) break;
    }//end for index
  }//end if threshold
  else{//the pulse is not too small, search forward
    for(index=pulse.start_index+1; index<=peak_index; index++){
      if(subtracted[index]<=threshold){
	index --;
	break;
      }//end if
    }//end for index
  }//end else
  //  std::cout<<"Index: "<<index<<'\t'<<chdata->SampleToTime(index)<<std::endl;
  //linear interpolation to calculate the time
  if(index>=pulse.start_index && index<peak_index && subtracted[index]>subtracted[index+1]){
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
  pulse.is_clean = (subtracted[pulse.start_index]>threshold) && (subtracted[pulse.end_index]>threshold);
  //calculate the pulse integral  
  if(chdata->integral.empty())
    pulse.integral = std::accumulate(subtracted+pulse.start_index,
				     subtracted+pulse.end_index,
				     0.);
  else pulse.integral = chdata->integral[pulse.end_index] - 
	 chdata->integral[pulse.start_index];
  pulse.npe = -pulse.integral/chdata->spe_mean;

  //Check to see if peak is saturated
  const double* wave = chdata->GetWaveform();
  int max_index = std::max_element(wave + pulse.start_index, 
				   wave + pulse.end_index) - wave;
  int min_index = std::min_element(wave + pulse.start_index, 
				   wave + pulse.end_index) - wave;
  pulse.max = wave[max_index];
  pulse.min = wave[min_index];

  if(wave[max_index] >= chdata->GetVerticalRange()){
    pulse.saturated = true;
  }
  if(wave[min_index] <= 0){
    pulse.saturated = true;
    int min_end_index = min_index + 1;
    while (wave[min_end_index] == 0 && min_end_index < pulse.end_index){
      min_end_index++;
    }
    pulse.peak_index = (int)(min_index + min_end_index)/2;
  }

  pulse.evaluated = true;
  //  if(chdata->integral.empty()) return 0;
  
  //look for the time it takes to reach X% of total integral
  //remember, integral is negative
  // int samp = pulse.start_index;
  // while( samp<pulse.end_index && 
  // 	 std::abs(chdata->integral[samp]-chdata->integral[pulse.start_index])<
  // 	 std::abs(pulse.integral)*0.05 ) samp++;
  // pulse.t05 = chdata->SampleToTime(samp)-pulse.start_time;
  // while( samp<pulse.end_index && 
  // 	 std::abs(chdata->integral[samp]-chdata->integral[pulse.start_index])<
  // 	 std::abs(pulse.integral)*0.10 ) samp++;
  // pulse.t10 = chdata->SampleToTime(samp)-pulse.start_time;
  // samp = pulse.end_index-1;
  // while( samp>=pulse.start_index && 
  // 	 std::abs(chdata->integral[samp]-chdata->integral[pulse.start_index])>
  // 	 std::abs(pulse.integral)*0.95 ) samp--;
  // pulse.t95 = chdata->SampleToTime(++samp)-pulse.start_time;
  // while( samp>=pulse.start_index && 
  // 	 std::abs(chdata->integral[samp]-chdata->integral[pulse.start_index])>
  // 	 std::abs(pulse.integral)*0.90 ) samp--;
  // pulse.t90 = chdata->SampleToTime(++samp)-pulse.start_time;
  

  return 0;
}
