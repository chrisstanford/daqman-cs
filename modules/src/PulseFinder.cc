#include "PulseFinder.hh"
#include "BaselineFinder.hh"
//#include "SumChannels.hh"
#include "Integrator.hh"
#include "ConvertData.hh"
//#include "RootWriter.hh"
#include "intarray.hh"
#include "TMath.h"
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cmath>

PulseFinder::PulseFinder() : 
  BaseModule(GetDefaultName(), "Search for individual physical scintillation pulses within the hardware trigger")
{
  AddDependency<ConvertData>();
  AddDependency<BaselineFinder>();
  AddDependency<Integrator>();
  
  //  RegisterParameter("search_mode", mode=DISCRIMINATOR);
  RegisterParameter("search_start_time", search_start_time = -10000,
		    "time in us to start searching for pulses");
  RegisterParameter("search_end_time", search_end_time = 10000,
		    "time in us to end searching for pulses");
  RegisterParameter("ch_pulse_filter_us", ch_pulse_filter_us,
		    "estimate of the pulse width in us");
  RegisterParameter("ch_pulse_thresholds", ch_pulse_thresholds,
		    "discriminator threshold values for pulse finding");
  RegisterParameter("ch_spike_thresholds", ch_spike_thresholds,
		    "discriminator threshold values for spike finding");
  RegisterParameter("relative_bls_threshold", relative_bls_threshold = false,
		    "Is the threshold relative to baseline sigma");
  RegisterParameter("relative_spe_threshold", relative_spe_threshold = false,
		    "Is the threshold relative to spe size");
  RegisterParameter("pulse_start_add_us", pulse_start_add_us = 0.02,
		    "time in us to add before pulse start ");
  RegisterParameter("pulse_end_add_us", pulse_end_add_us = 0.02,
		    "time in us to add after pulse end");
//   RegisterParameter("",  = ,
// 		    "");

}
  
PulseFinder::~PulseFinder()
{
}

int PulseFinder::Initialize()
{
  if(search_end_time<=search_start_time || 
     (relative_bls_threshold && relative_spe_threshold) ||
     pulse_start_add_us<0 || pulse_end_add_us<0)
    return -1;

  //pulse filter length value has to be positive
  for(id_map::iterator it = ch_pulse_filter_us.begin();
      it != ch_pulse_filter_us.end(); it++){
    if(it->second<=0){
      Message(ERROR)<<"Channel "<<it->first<< " has pulse width value set to be negative or 0!\n";
      return 1;
    }//end if
  }//end for

  //threshold can not be 0
  for(id_map::iterator it = ch_pulse_thresholds.begin();
      it != ch_pulse_thresholds.end(); it++){
    if(it->second==0){
      Message(ERROR)<<"Channel "<<it->first<< " has pulse search threshold set to be 0!\n";
      return 1;
    }//end if
  }//end for
  for(id_map::iterator it = ch_spike_thresholds.begin();
      it != ch_spike_thresholds.end(); it++){
    if(it->second==0){
      Message(ERROR)<<"Channel "<<it->first<< " has spike search threshold set to be 0!\n";
      return 1;
    }//end if
  }//end for


  return 0;
}

int PulseFinder::Finalize()
{
  return 0;
}

int PulseFinder::Process(EventPtr evt)
{
    EventDataPtr event = evt->GetEventData();

    //Loop over all channels and evaluate pulse edges individually
    for (size_t ch = 0; ch < event->channels.size(); ch++){
      
      //note: if we search pulses with the sum channel, we need to skip
      //the channels that are summed, then propagate pulses there
      ChannelData * chdata = &(event->channels[ch]);
      if(_skip_channels.find(chdata->channel_id) != _skip_channels.end()) continue;
      if(! (chdata->baseline.found_baseline) ||
	 (int(chdata->subtracted_waveform.size()) != chdata->nsamps) ) continue;
      
      int result = FindChannelPulses(chdata);
      //      if(result) return result; //we may not want to return if one channel fails
      if(result) continue;

      //note: we will not search for spikes in the sum channel
      //only add up the spikes from summed channels
      //spike class needs to specify which channel it comes from

    } //end loop over channels
    return 0;
}

int PulseFinder::FindChannelPulses(ChannelData* chdata){
  
  std::vector<int> start_index;
  std::vector<int> end_index;
  
  int pulse_start_add_nsamps = pulse_start_add_us*chdata->sample_rate+1;
  int pulse_end_add_nsamps = pulse_end_add_us*chdata->sample_rate+1;
	
  double threshold = 0;
  id_map::iterator it = ch_pulse_thresholds.find(chdata->channel_id);
  if(it != ch_pulse_thresholds.end() && it->second!=0)
    threshold = it->second;
  else{
    it = ch_pulse_thresholds.find(ChannelData::CH_DEFAULT);
    if(it != ch_pulse_thresholds.end() && it->second!=0)
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

  double filter_time = 0;
  //for real pulse finding, we would want to do a convulation here
  it = ch_pulse_filter_us.find(chdata->channel_id);
  if(it != ch_pulse_filter_us.end() && it->second>0)
    filter_time = it->second;
  else{
    it = ch_pulse_filter_us.find(ChannelData::CH_DEFAULT);
    if(it != ch_pulse_filter_us.end() && it->second>0)
      filter_time =it->second;
    else filter_time = 0;
  }//end else
  int filter_nsamps = filter_time*chdata->sample_rate+1;
  const double * wave = chdata->GetBaselineSubtractedWaveform();
  std::vector<double> smoothed;
  if(filter_nsamps>1){
    smoothed.resize(chdata->nsamps);
    //    bool is_good = SmoothWaveform(smoothed, chdata->nsamps, wave, filter_nsamps);
    bool is_good = RunningSumWaveform(smoothed, chdata->nsamps, wave, filter_nsamps);
    if(!is_good) return -1;
    wave = &(smoothed[0]);
  }
  
  int search_result = 
    DiscriminatorSearch(chdata, wave, start_index, end_index, threshold, 
			pulse_start_add_nsamps, pulse_end_add_nsamps);
  if(search_result) return search_result;	

  std::vector<int> start_index_split;
  std::vector<int> end_index_split;
  
  for(size_t index=0; index<start_index.size(); index++){
//     std::cout<<chdata->channel_id<<'\t'<<index<<'\t'<<chdata->SampleToTime(start_index.at(index))
// 	     <<'\t'<<chdata->SampleToTime(end_index.at(index))<<std::endl;
    if(start_index.at(index)>=end_index.at(index)) continue;
    std::vector<int> start_index_tmp;
    std::vector<int> end_index_tmp;
    int split_result= ResolvePileUps(chdata, wave,start_index_tmp,end_index_tmp, threshold/2, 
				     filter_nsamps, start_index.at(index),end_index.at(index));
    
    if(split_result) continue;
    //    std::cout<<"*** Great *** RelativeThresholdSearch succeeded! "<<std::endl;
    start_index_split.insert(start_index_split.end(), start_index_tmp.begin(), start_index_tmp.end());
    end_index_split.insert(end_index_split.end(), end_index_tmp.begin(), end_index_tmp.end());

  }//end for index  
  
  //todo: we may want to move everything to a evaluatepulse module
  int search_start_index=chdata->TimeToSample(search_start_time, true);
  chdata->pulses.reserve(start_index_split.size());
  for (size_t i = 0; i < start_index_split.size();  i++)
    {
      if (start_index_split[i] >= end_index_split[i]) continue;

      Pulse pulse;
      pulse.start_index = start_index_split[i];
      pulse.end_index = end_index_split[i];
      EvaluatePulse(pulse, chdata, -1);
      //evaluate the pulse parameters depending on other pulses
      if(chdata->pulses.size()){
	pulse.time_since_last=pulse.half_max_time-chdata->pulses.back().end_time;
	if(pulse.start_index<=chdata->pulses.back().end_index)
	  {pulse.is_clean = false; chdata->pulses.back().is_clean = false;}
      }
      else{
	pulse.time_since_last=pulse.half_max_time-chdata->SampleToTime(search_start_index);
	if(pulse.start_index<=search_start_index) pulse.is_clean = false;
      }
      //  pulse.Print(chdata->channel_id, chdata->pulses.size());
      chdata->pulses.push_back(pulse);
	    
    } // end for loop over pulses
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
  chdata->npulses = chdata->pulses.size();

  return 0;

}


int PulseFinder::FindChannelSpikes(ChannelData* chdata){

  std::vector<int> start_index;
  std::vector<int> end_index;
  
  double threshold = 0;
  id_map::iterator it = ch_spike_thresholds.find(chdata->channel_id);
  if(it != ch_spike_thresholds.end() && it->second!=0)
    threshold = it->second;
  else{
    it = ch_spike_thresholds.find(ChannelData::CH_DEFAULT);
    if(it != ch_spike_thresholds.end() && it->second!=0) 
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

  int spike_edge_add_nsamps = 5;
  int spike_split_comp_step = 3;
  const double * wave = chdata->GetBaselineSubtractedWaveform();
  int search_result = 
    DiscriminatorSearch(chdata, wave, start_index, end_index, threshold,
			spike_edge_add_nsamps, spike_edge_add_nsamps);
  if(search_result) return search_result;	

  std::vector<int> start_index_split;
  std::vector<int> end_index_split;

  for(size_t index=0; index<start_index.size(); index++){
    //     std::cout<<chdata->channel_id<<'\t'<<index<<'\t'<<chdata->SampleToTime(start_index.at(index))
    // 	     <<'\t'<<chdata->SampleToTime(end_index.at(index))<<std::endl;
    if(start_index.at(index)>=end_index.at(index)) continue;
    std::vector<int> start_index_tmp;
    std::vector<int> end_index_tmp;
    int split_result= ResolvePileUps(chdata, wave,start_index_tmp,end_index_tmp, threshold, 
				     spike_split_comp_step, start_index.at(index),end_index.at(index));
    
    if(split_result) continue;
    //    std::cout<<"*** Great *** RelativeThresholdSearch succeeded! "<<std::endl;
    start_index_split.insert(start_index_split.end(), start_index_tmp.begin(), start_index_tmp.end());
    end_index_split.insert(end_index_split.end(), end_index_tmp.begin(), end_index_tmp.end());

  }//end for index  

  //todo: we may want to have an evaluate spike function
  for (size_t i = 0; i < start_index_split.size();  i++){
    if (start_index_split[i] >= end_index_split[i]) continue;
    
    //  pulse.Print(chdata->channel_id, chdata->pulses.size());
    //    chdata->spikes.push_back(spike);
    
  } // end for loop over spikes
  
  return 0;
}

/// resolve pileup pulses/spikes
/// only trust the result if  within a pulse/spike
/// may pickup baseline fluctuations otherwise
int PulseFinder::ResolvePileUps(ChannelData* chdata, const double * wave,
				std::vector<int>& start_index, 
				std::vector<int>& end_index,
				double threshold, int step, 
				int search_start, int search_end){

  if(!chdata || !wave || step<1 || threshold==0 || search_start<0 || 
     search_end>=chdata->nsamps) return -1;
  if(step<3) step=3;  
  if(start_index.size()) start_index.clear();
  if(end_index.size()) end_index.clear();
  
  std::vector<int> turnpoints;
  std::vector<int> peak;
  turnpoints.reserve(100);  
  peak.reserve(100);  

  //we could use zero threshold to find the peaks and then compare peak to valley with threshold
  //try the simple threshold comparison first
  for(int samp=search_start+1; samp<=search_end-1; samp++){
    //here we come to a peak
    if((RelativeThresholdCrossed(wave[samp-1],wave[samp],threshold) ||
	(samp>=search_start+step && RelativeThresholdCrossed(wave[samp-step],wave[samp],threshold)) ) &&
       (RelativeThresholdCrossed(wave[samp+1],wave[samp],threshold) ||
	(samp+step<=search_end && RelativeThresholdCrossed(wave[samp+step],wave[samp],threshold)) ) ){
      //we need to make sure the first/last peak is complete
      if(peak.size()==0){
	turnpoints.push_back(search_start);
	peak.push_back(-1);
      }
      //now this operation is safe
      if(peak.back()==1){ //there is a peak just before this one, and there was not a valley
	if(FirstAmplitudeIsSmaller(wave[turnpoints.back()], wave[samp], threshold)) turnpoints.back()=samp;
      }//end if a peak before this one
      else {//if there is a valley right before it or nothing
	turnpoints.push_back(samp);
	peak.push_back(1);
      }//end else
    }//end if coming to a peak
    //here we come to a valley
    else if((RelativeThresholdCrossed(wave[samp],wave[samp-1],threshold) ||
	     (samp>=search_start+step && RelativeThresholdCrossed(wave[samp],wave[samp-step],threshold)) ) &&
	    (RelativeThresholdCrossed(wave[samp],wave[samp+1],threshold) ||
	     (samp+step<=search_end && RelativeThresholdCrossed(wave[samp],wave[samp+step],threshold)) ) ){
      //if there is another valley right before it but no peak in between
      if(peak.size()>0 && peak.back()==-1){ //there is a peak just before this one, and there was not a valley
	if(FirstAmplitudeIsSmaller(wave[samp], wave[turnpoints.back()], threshold)) turnpoints.back()=samp;
      }//end if a valley before this one
      else {//if there is a peak right before it or nothing
	turnpoints.push_back(samp);
	peak.push_back(-1);
      }//end else
    }//end if coming to a valley

  }//end for samp
  //we need to make sure the first/last peak is complete
  if(peak.back()==1){
    turnpoints.push_back(search_end);
    peak.push_back(-1);
  }

  //check the peak array: -1 1 -1 1 ... 1 -1
  //  std::cout<<std::endl<<"Pulse finding before combining"<<std::endl;
  for(size_t ii=0; ii<peak.size(); ii++){
    //    std::cout<<ii<<'\t'<<turnpoints.at(ii)<<'\t'<<peak.at(ii)<<std::endl;
    if((ii%2 && peak.at(ii)!=1) || (!(ii%2) && peak.at(ii)!=-1)) std::cout<<"Error in peak finding"<<std::endl;
  }

  // std::cout<<std::endl<<"Pulse finding after combining"<<std::endl;
  // for(size_t ii=0; ii<peak.size(); ii++)
  //   std::cout<<turnpoints.at(ii)<<'\t'<<peak.at(ii)<<std::endl;

  for(size_t ii=0; ii<turnpoints.size(); ii++){
    if(peak.at(ii)!=-1) continue;
    if((ii+1)!=turnpoints.size()) start_index.push_back(turnpoints.at(ii));
    if(ii) end_index.push_back(turnpoints.at(ii));
  }
  if(start_index.size()==end_index.size()) return 0;
  else{
    if(start_index.size()) start_index.clear();
    if(end_index.size()) end_index.clear();
    return -1;
  }
}

/// Search for pulses using a simple discrimination threshold
int PulseFinder::DiscriminatorSearch(ChannelData* chdata, const double * wave,
				     std::vector<int>& start_index, 
				     std::vector<int>& end_index, double threshold,
				     int pulse_start_add_nsamps, int pulse_end_add_nsamps){

  if(!chdata || !wave || threshold==0 || pulse_start_add_nsamps<0 || pulse_end_add_nsamps<0)
    return -1;

  if(start_index.size()) start_index.clear();
  if(end_index.size()) end_index.clear();
  bool in_pulse = false;
  int search_start_index = chdata->TimeToSample(search_start_time, true);
  int search_end_index = chdata->TimeToSample(search_end_time, true);

  for(int index = search_start_index; index<=search_end_index; index++){
    //we are beyond threshold
    if(FirstAmplitudeIsSmaller(threshold,wave[index], threshold)){
      //if just come to a pulse
      //if(!in_pulse){
      if(index==0 || FirstAmplitudeIsSmaller(wave[index-1],threshold, threshold)){
	in_pulse = true;
	if(index>=search_start_index+pulse_start_add_nsamps) start_index.push_back(index-pulse_start_add_nsamps );
	else start_index.push_back(search_start_index);
      }//end if just come to a pulse
      //if we are about to exit a pulse
      if(index==chdata->nsamps-1 || FirstAmplitudeIsSmaller(wave[index+1], threshold,threshold)){
	if(in_pulse){
	  if(index+pulse_end_add_nsamps<=search_end_index) end_index.push_back(index+pulse_end_add_nsamps);
	  else end_index.push_back(search_end_index);
	}//end if in pulse
	in_pulse = false;
      }//end if exit a pulse
    }//end if wave index <=
  }//end for int index
  if(in_pulse) end_index.push_back(search_end_index);

  if(start_index.size()!=end_index.size()){
    std::cerr<<start_index.size()<<" =/= "<<end_index.size()<<std::endl;
    std::cerr<<"PulseFinder failed to find correct pulses, give up this event."<<std::endl;
    start_index.clear();
    end_index.clear();
    return -1;
  }
 
  //resoving overlapping pulses
  for(size_t index=1; index<start_index.size(); index++){
    // std::cout<<chdata->channel_id<<'\t'<<index<<'\t'<<chdata->SampleToTime(start_index.at(index))
    // 	     <<'\t'<<chdata->SampleToTime(end_index.at(index))<<std::endl;
    if(start_index.at(index)<=end_index.at(index-1)){
      int middle_index = (start_index.at(index)+end_index.at(index-1))/2;
      double min_value = wave[middle_index];
      //      for(int j=start_index.at(index); j<=end_index.at(index-1); j++){
      for(int j=end_index.at(index-1)+1-pulse_end_add_nsamps; j<start_index.at(index)+pulse_start_add_nsamps; j++){
	if(j<0) continue;
	else if(j>=chdata->nsamps) break;
	if(FirstAmplitudeIsSmaller(min_value,wave[j], threshold)) continue;
        min_value = wave[j];
        middle_index = j;
      }//end for j loop
      start_index.at(index) = middle_index;
      end_index.at(index-1) = middle_index;
    }//end if statement
  }//end for index loop
  

    //  for(size_t index=0; index<start_index_tmp.size(); index++)
    //    std::cout<<chdata->channel_id<<'\t'<<index<<'\t'<<chdata->SampleToTime(start_index_tmp.at(index))
    //      	     <<'\t'<<chdata->SampleToTime(end_index_tmp.at(index))<<std::endl;
  return 0;
}


// std::ostream& operator<<(std::ostream& out, const PulseFinder::SEARCH_MODE& m)
// {
//   switch(m){
//   case PulseFinder::VARIANCE:
//     out<<"VARIANCE";
//     break;
//   case PulseFinder::DISCRIMINATOR:
//     out<<"DISCRIMINATOR";
//     break;
//   case PulseFinder::INTEGRAL:
//     out<<"INTEGRAL";
//     break;
//   case PulseFinder::CURVATURE:
//     out<<"CURVATURE";
//     break;
//   }
//   return out;    
// }

// std::istream& operator>>(std::istream& in, PulseFinder::SEARCH_MODE& m)
// {
//   std::string dummy;
//   in>>dummy;
//   if(dummy == "VARIANCE" || dummy == "variance")
//     m = PulseFinder::VARIANCE;
//   else if (dummy == "DISCRIMINATOR" || dummy == "discriminator")
//     m = PulseFinder::DISCRIMINATOR;
//   else if (dummy == "INTEGRAL" || dummy == "integral")
//     m = PulseFinder::INTEGRAL;
//   else if (dummy == "CURVATURE" || dummy == "curvature")
//     m = PulseFinder::CURVATURE;
//   else{
//     throw std::invalid_argument(dummy+"is not a valid value for search_mode!");
//   }
//   return in;
// }

