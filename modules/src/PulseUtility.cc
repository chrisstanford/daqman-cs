#ifndef PULSE_UTILITY_c
#define PULSE_UTILITY_c

#include "ChannelData.hh"
#include "PulseUtility.hh"
#include <algorithm>
#include <numeric>
#include <cmath>

//determine if x1->x2 cross the relative threshold, either positive or negative
bool RelativeThresholdCrossed(double x1, double x2, double threshold){
  if(threshold>0) return ( (x2-x1)>threshold );
  else if(threshold<0) return ( (x2-x1)<threshold );
  else {
    std::cerr<<"*** Error *** Comparing with 0 threshold *** Error ***"<<std::endl;
    return false;
  }
}

//determine either x1 or x2 has smaller amplitude in the direction given by sign
bool FirstAmplitudeIsSmaller(double x1, double x2, double dumb_sign){
  if(dumb_sign>0) return (x1<x2);
  else if(dumb_sign<0) return (x1>x2);
  else {
    std::cerr<<"*** Error *** Comparing amplitudes without a given sign *** Error ***"<<std::endl;
    return false;
  }
}

//this function can not be part of the pulse class because pulse belongs to chdata
//it can not be part of pulse finder because EvalRoi calls it too
//function to calculate pulse parameters, if max_peak_time>0
//the peak of the pulse should be within (start time, start time + max_peak_time)
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

  return 0;
}

int RelativeThresholdSearch(std::vector<double> wave, double start_threshold, double end_threshold,
			    std::vector<int> & start_index, std::vector<int> & end_index, 
			    int pulse_edge_add, int step_size, int search_start, int search_end){
  std::cerr<<"Begining of RelativeThresholdSearch"<<std::endl;
  if(wave.size()<=1 || step_size<1 || start_threshold==0 || end_threshold==0) return -1;
  std::cerr<<"RelativeThresholdSearch"<<std::endl;
  const int nsamps = wave.size();
  if(search_start<0) search_start = 0;
  if(search_end<0 || search_end>=nsamps) search_end = nsamps-1;
  if(start_index.size()) start_index.clear();
  if(end_index.size()) end_index.clear();

  std::vector<int> mask(nsamps,0);
  //let's first mark the pulse starts and ends
  for(int i=search_start; i<=search_end; i++){
    if(i<1) continue;
    //the continue commands let pulse end overrule pulse start, and let nearby points overrule far away ones
    if(RelativeThresholdCrossed(wave[i-1],wave[i], end_threshold)) {mask[i]=-1; continue;}
    else if(RelativeThresholdCrossed(wave[i-1],wave[i], start_threshold)) {mask[i]=1; continue;}
    if(i<2*step_size) continue;
    if(RelativeThresholdCrossed(wave[i-2*step_size],wave[i], end_threshold)) {mask[i]=-1; continue;}
    else if(RelativeThresholdCrossed(wave[i-2*step_size],wave[i], start_threshold)) {mask[i]=1; continue;}
    if(i<3*step_size) continue;
    if(RelativeThresholdCrossed(wave[i-3*step_size],wave[i], end_threshold)) {mask[i]=-1; continue;}
    else if(RelativeThresholdCrossed(wave[i-3*step_size],wave[i], start_threshold)) {mask[i]=1; continue;}
  }//end for

  //combine pulses and remove false pulses
  //we are more agreesive here than in baseline finder because we don't want to miss pulses
  for(int i=0; i<nsamps; i++){
    if(mask[i]==1 && start_index.size()==end_index.size()){
      start_index.push_back(i);
      continue;
    }
    if(mask[i]==-1){
      if(start_index.size()==end_index.size() && end_index.size()>0) end_index.back()=i;
      else if(start_index.size()==end_index.size()+1) end_index.push_back(i);
    }
  }//end for

  if(start_index.size()!=end_index.size() || start_index.size()==0){
    std::cerr<<"start and end index vector size: "<<start_index.size()<<'\t'<<end_index.size()<<std::endl;
    if(start_index.size()) start_index.clear();
    if(end_index.size()) end_index.clear();
    return -1;
  }

  //let's resolve possible overlaps
  for(size_t index=0; index<start_index.size(); index++){
    // std::cout<<chdata->channel_id<<'\t'<<index<<'\t'<<chdata->SampleToTime(start_index.at(index))
    // 	     <<'\t'<<chdata->SampleToTime(end_index.at(index))<<std::endl;
    if(start_index.at(index)>=search_start+pulse_edge_add) start_index.at(index)-=pulse_edge_add;
    else start_index.at(index)=search_start;
    if(end_index.at(index)+pulse_edge_add<=search_end) end_index.at(index) += pulse_edge_add;
    else end_index.at(index)=search_end;
    if(index==0) continue;

    if(start_index.at(index)<=end_index.at(index-1)){
      int middle_index = (start_index.at(index)+end_index.at(index-1))/2;
      double min_value = wave[middle_index];
      //      for(int j=start_index.at(index); j<=end_index.at(index-1); j++){
      for(int j=end_index.at(index-1)+1-pulse_edge_add; j<start_index.at(index)+pulse_edge_add; j++){
	if(j<0) continue;
	else if(j>=nsamps) break;
	if(FirstAmplitudeIsSmaller(min_value,wave[j], start_threshold)) continue;
        min_value = wave[j];
        middle_index = j;
      }//end for j loop
      start_index.at(index) = middle_index;
      end_index.at(index-1) = middle_index;
    }//end if statement
  }//end for index loop

  return 0;
}

#endif
