#include "PulseFinder.hh"
#include "BaselineFinder.hh"
//#include "SumChannels.hh"
#include "Integrator.hh"
#include "ConvertData.hh"
//#include "RootWriter.hh"
#include "PulseUtility.hh"
#include "intarray.hh"
#include "TMath.h"
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <cstdio>

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

    bool process_sum_ch = false;
    if(_skip_channels.find(ChannelData::CH_SUM)==_skip_channels.end() &&
       (ch_pulse_thresholds.find(ChannelData::CH_SUM)!= ch_pulse_thresholds.end() ||
	ch_pulse_thresholds.find(ChannelData::CH_DEFAULT)!= ch_pulse_thresholds.end()))
      process_sum_ch = true;
    std::set<int> &channels_summed = event->channels_summed;

    //Loop over all channels and evaluate pulse edges individually
    for (size_t ch = 0; ch < event->channels.size(); ch++){
      
      ChannelData * chdata = &(event->channels[ch]);
      if(_skip_channels.find(chdata->channel_id) != _skip_channels.end()) continue;
      if(! (chdata->baseline.found_baseline) ||
	 (int(chdata->subtracted_waveform.size()) != chdata->nsamps) ) continue;

      //find the spikes if the users specify a threshold
      //note: we will not search for spikes in the sum channel
      //only add up the spikes from summed channels
      int result = 0;
      if(chdata->channel_id != ChannelData::CH_SUM) FindChannelSpikes(chdata);

      //if we search pulses with the sum channel, we need to skip
      //the channels that are summed, then propagate pulses there
      if(process_sum_ch && channels_summed.find(chdata->channel_id)!=channels_summed.end())
	continue;
      
      result = FindChannelPulses(chdata);
      //      if(result) return result; //we may not want to return if one channel fails
      //      if(result) continue;
    } //end loop over channels

    //propagate sum channel pulses to the channels summed
    if(channels_summed.size()>0){
      //get the sum channel pulses
      std::vector<Pulse> & sum_pulses = event->GetChannelByID(ChannelData::CH_SUM)->pulses;
      if(sum_pulses.size()>0){

	for(std::set<int>::iterator it=channels_summed.begin(); it!=channels_summed.end(); it++){
	  //	for(size_t ch=0; ch<channels_summed.size(); ch++){
	  ChannelData * chdata = event->GetChannelByID(*it);
	  for(size_t jj=0; jj<sum_pulses.size(); jj++){
	    Pulse pulse;
	    pulse.start_index = sum_pulses.at(jj).start_index;
	    pulse.end_index = sum_pulses.at(jj).end_index;
	    chdata->pulses.push_back(pulse);
	  }//end for jj
	  chdata->npulses = chdata->pulses.size();
	}//end for it
      }//end if sum_pulses.size()
    }//end if channels_summed.size()>0
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
  
  //we don't trust the sum channel baseline sigma
  if(relative_bls_threshold){
    if(chdata->channel_id!=ChannelData::CH_SUM && chdata->baseline.sigma>0) 
      threshold *= chdata->baseline.sigma;
    else{
      Message(ERROR)<<"Channel "<<chdata->channel_id<< " has invalid baseline!\n";
      return 0;
    }
  }

  if(relative_spe_threshold){
    if(chdata->spe_mean>0){
      threshold *= chdata->spe_mean;
      if(chdata->channel_id!=ChannelData::CH_SUM && chdata->spe_mean==1){
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
  std::vector<double> summed;
  std::vector<double> smoothed;
  if(filter_nsamps>1){
    summed.resize(chdata->nsamps);
    //    bool is_good = SmoothWaveform(summed, chdata->nsamps, wave, filter_nsamps);
    bool is_good = RunningSumWaveform(summed, chdata->nsamps, wave, filter_nsamps);
    if(!is_good) return -1;
    is_good = SmoothWaveform(smoothed, chdata->nsamps, &(summed[0]),filter_nsamps/4);
    if(!is_good) return -1;

    //let's save the sum waveform for checks
    std::remove("event.txt");
    std::ofstream ofs ("event.txt", std::ofstream::out);
    for(size_t ab=0; ab<summed.size(); ab++)
      ofs <<chdata->SampleToTime(ab)<<'\t'<<summed.at(ab)<<'\t'<<smoothed.at(ab)<<std::endl;
    ofs.close();
    
    //    wave = &(summed[0]);
    wave = &(smoothed[0]);
  }
  
  int search_result = 
    DiscriminatorSearch(chdata, wave, start_index, end_index, threshold, 
			pulse_start_add_nsamps, pulse_end_add_nsamps);
  if(search_result) return search_result;	
  
  std::vector<int> start_index_split;
  std::vector<int> end_index_split;
  
  for(size_t index=0; index<start_index.size(); index++){
    std::cout<<"Pulse to resolve pileup: "<<chdata->channel_id<<'\t'<<index<<'\t'<<chdata->SampleToTime(start_index.at(index))
	      <<'\t'<<chdata->SampleToTime(end_index.at(index))<<std::endl;
    if(start_index.at(index)>=end_index.at(index)) continue;
    std::vector<int> start_index_tmp;
    std::vector<int> end_index_tmp;
    int split_result= ResolvePileUps(chdata, wave,start_index_tmp,end_index_tmp, threshold/2, 
				     filter_nsamps/2, start_index.at(index),end_index.at(index));
    
    if(split_result || start_index_tmp.size()==0) continue;
    std::cout<<"*** Great *** resolving pileup succeeded! "<<std::endl;
    start_index_split.insert(start_index_split.end(), start_index_tmp.begin(), start_index_tmp.end());
    end_index_split.insert(end_index_split.end(), end_index_tmp.begin(), end_index_tmp.end());

  }//end for index  
  
  chdata->pulses.reserve(start_index_split.size());
  for (size_t i = 0; i < start_index_split.size();  i++)
    {
      std::cout<<chdata->channel_id<<'\t'<<i<<'\t'<<chdata->SampleToTime(start_index_split.at(i))
	       <<'\t'<<chdata->SampleToTime(end_index_split.at(i))<<std::endl;
      if (start_index_split[i] >= end_index_split[i]) continue;
      Pulse pulse;
      pulse.start_index = start_index_split[i];
      pulse.end_index = end_index_split[i];
      chdata->pulses.push_back(pulse);
    } // end for loop over pulses
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
     search_end>=chdata->nsamps || search_start>=search_end) return -1;
  if(start_index.size()) start_index.clear();
  if(end_index.size()) end_index.clear();
  const int min_step = 3; //this is the best value for spike counting
  if(step<min_step) step=min_step;  
  //no point searching pileups if it is a small pulse
  if(search_end-search_start<=min_step){
    start_index.push_back(search_start);
    end_index.push_back(search_end);
    return 0;
  }
  
  std::cout<<"Resolving overlaps for pulse: "<<chdata->SampleToTime(search_start)<<'\t'<<chdata->SampleToTime(search_end)<<std::endl;

  std::vector<int> turnpoints;
  std::vector<int> peak;
  turnpoints.reserve(100);  
  peak.reserve(100);  

  for(int samp=search_start+1; samp<=search_end-1; samp++){
    //here we come to a possible peak, w/ bigger amplitude than +/-1, or +/-step, or +/-step/2
    if( (RelativeThresholdCrossed(wave[samp-1],wave[samp],threshold) ||
	 (samp>=search_start+step && RelativeThresholdCrossed(wave[samp-step],wave[samp],threshold)) ||
	 (step/2>min_step && samp>=search_start+step/2 && RelativeThresholdCrossed(wave[samp-step/2],wave[samp],threshold)) ) &&
	(RelativeThresholdCrossed(wave[samp+1],wave[samp],threshold) ||
	 (samp+step<=search_end && RelativeThresholdCrossed(wave[samp+step],wave[samp],threshold)) ||
	 (step/2>min_step && samp+step/2<=search_end && RelativeThresholdCrossed(wave[samp+step/2],wave[samp],threshold))) ){
      std::cout<<"Possible peak at: "<<chdata->SampleToTime(samp)<<std::endl;
      //we need to make sure the first/last peak is complete
      if(peak.size()==0){
	turnpoints.push_back(search_start);
	peak.push_back(-1);
      }
      //there is a peak just before this one, see which one is bigger
      if(peak.back()==1){
	if(FirstAmplitudeIsSmaller(wave[turnpoints.back()], wave[samp], threshold)) turnpoints.back()=samp;
      }
      else {//if there is a valley right before it
	turnpoints.push_back(samp);
	peak.push_back(1);
      }//end else
    }//end if coming to a peak
    //here we come to a valley
    else if((RelativeThresholdCrossed(wave[samp],wave[samp-1],threshold) ||
	     (samp>=search_start+step && RelativeThresholdCrossed(wave[samp],wave[samp-step],threshold)) ||
	     (step>min_step*2 && samp>=search_start+step/2 && RelativeThresholdCrossed(wave[samp],wave[samp-step/2],threshold)) ) &&
	    (RelativeThresholdCrossed(wave[samp],wave[samp+1],threshold) ||
	     (samp+step<=search_end && RelativeThresholdCrossed(wave[samp],wave[samp+step],threshold)) ||
	     (step>min_step*2 && samp+step/2<=search_end && RelativeThresholdCrossed(wave[samp],wave[samp+step/2],threshold))) ){
      std::cout<<"Possible valley at: "<<chdata->SampleToTime(samp)<<std::endl;
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
    std::cout<<ii<<'\t'<<chdata->SampleToTime(turnpoints.at(ii))<<'\t'<<peak.at(ii)<<std::endl;
    if((ii%2 && peak.at(ii)!=1) || (!(ii%2) && peak.at(ii)!=-1)) std::cout<<"Error in peak finding"<<std::endl;
  }

  for(size_t ii=0; ii<turnpoints.size(); ii++){
    if(peak.at(ii)!=-1) continue;
    if((ii+1)!=turnpoints.size()) start_index.push_back(turnpoints.at(ii));
    if(ii) end_index.push_back(turnpoints.at(ii));
  }

  if(start_index.size()==0){
    start_index.push_back(search_start);
    end_index.push_back(search_end);
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

  std::cerr<<"Searching for pulses in channel: "<<chdata->channel_id<<std::endl;
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
	std::cerr<<"Find pulse start: "<<chdata->SampleToTime(index)<<'\t'
		 <<chdata->SampleToTime(start_index.back())<<std::endl;
      }//end if just come to a pulse
      //if we are about to exit a pulse
      if(index==chdata->nsamps-1 || FirstAmplitudeIsSmaller(wave[index+1], threshold,threshold)){
	if(in_pulse){
	  if(index+pulse_end_add_nsamps<=search_end_index) end_index.push_back(index+pulse_end_add_nsamps);
	  else end_index.push_back(search_end_index);
	  std::cerr<<"Find pulse end: "<<chdata->SampleToTime(index)<<'\t'
		   <<chdata->SampleToTime(end_index.back())<<std::endl;
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
  for(size_t index=0; index<start_index.size(); index++){
    std::cout<<"Overlap? "<<chdata->channel_id<<'\t'<<index<<'\t'<<chdata->SampleToTime(start_index.at(index))
     	     <<'\t'<<chdata->SampleToTime(end_index.at(index))<<std::endl;
    if(index<1) continue;
    if(start_index.at(index)<=end_index.at(index-1)){
      int middle_index = (start_index.at(index)+pulse_start_add_nsamps+end_index.at(index-1)-pulse_end_add_nsamps)/2;
      double min_value = wave[middle_index];
      //      for(int j=start_index.at(index); j<=end_index.at(index-1); j++){
      std::cout<<"checking range: "<<chdata->SampleToTime(end_index.at(index-1)+1-pulse_end_add_nsamps)
	       <<'\t'<<chdata->SampleToTime(start_index.at(index)+pulse_start_add_nsamps)<<std::endl;
      for(int j=end_index.at(index-1)+1-pulse_end_add_nsamps; j<start_index.at(index)+pulse_start_add_nsamps; j++){
	if(j<0) continue;
	else if(j>=chdata->nsamps) break;
	if(FirstAmplitudeIsSmaller(min_value, wave[j], threshold)) continue;
        min_value = wave[j];
        middle_index = j;
      }//end for j loop
      std::cout<<"found minimum: "<<chdata->SampleToTime(middle_index)<<std::endl;
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

