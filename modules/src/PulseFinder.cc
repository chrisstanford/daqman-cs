#include "PulseFinder.hh"
#include "ConvertData.hh"
#include "BaselineFinder.hh"
#include "Integrator.hh"
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
  
PulseFinder::~PulseFinder(){;}

int PulseFinder::Initialize(){

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
  //we don't need this parameter to be user-specified
  spike_edge_add_nsamps=5;
  return 0;
}

int PulseFinder::Finalize(){
  return 0;
}

int PulseFinder::Process(EventPtr evt){

  EventDataPtr event = evt->GetEventData();

  bool process_sum_ch = false;
  if(_skip_channels.find(ChannelData::CH_SUM)==_skip_channels.end() &&
     (ch_pulse_thresholds.find(ChannelData::CH_SUM)!= ch_pulse_thresholds.end() ||
      ch_pulse_thresholds.find(ChannelData::CH_DEFAULT)!= ch_pulse_thresholds.end()))
    process_sum_ch = true;
  std::set<int> &channels_summed = event->channels_summed;

  //Loop over all channels 
  for (size_t ch = 0; ch < event->channels.size(); ch++){
      
    ChannelData * chdata = &(event->channels[ch]);
    if(_skip_channels.find(chdata->channel_id) != _skip_channels.end()) continue;
    if(! (chdata->baseline.found_baseline) ||
       (int(chdata->subtracted_waveform.size()) != chdata->nsamps) ||
       int(chdata->integral.size()) != chdata->nsamps) continue;

    //find the spikes if the users specify a threshold
    //note: we will not search for spikes in the sum channel
    //only add up the spikes from summed channels
    int result = 0;
    if(chdata->channel_id != ChannelData::CH_SUM) result=FindChannelSpikes(chdata);

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

  //let's first find the threshold for this channel  
  double threshold = 0;
  id_map::iterator it = ch_pulse_thresholds.find(chdata->channel_id);
  if(it != ch_pulse_thresholds.end()){
    if(it->second!=0) threshold = it->second;
    else return 0;//skip this channel is a 0 threshold is specified
  }//end if a threshold is specified
  else{//look at the default channel is this channel is not specified
    it = ch_pulse_thresholds.find(ChannelData::CH_DEFAULT);
    if(it != ch_pulse_thresholds.end() && it->second!=0)
      threshold =it->second;
    else return 0; //we simply don't process this channel
  }//end else

  //lets see what is the unit of the threshold for channels other than the sum channel
  //sum channel unit is already number of spes
  if(chdata->channel_id!=ChannelData::CH_SUM){  
    if(relative_bls_threshold){
      if(chdata->baseline.sigma>0) threshold *= chdata->baseline.sigma;
      else{
	Message(ERROR)<<"Channel "<<chdata->channel_id<< " has invalid baseline!\n";
	return 0;
      }//end else
    }//end if bls
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
  }//end if refining threshold checks

  //To find pulses, we will do a convulation here
  //first produce a running sum, and then smooth it (running average)
  double filter_time = 0;
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
  const double * integral = chdata->GetIntegralWaveform();
  std::vector<double> summed;
  std::vector<double> smoothed;
  if(filter_nsamps>1){
    //    bool is_good = RunningSumWaveform(summed, chdata->nsamps, wave, filter_nsamps);
    bool is_good = RunningSumFromIntegral(summed, chdata->nsamps, integral, filter_nsamps);
    if(!is_good) return -1;
    is_good = SmoothWaveform(smoothed, chdata->nsamps, &(summed[0]),filter_nsamps/4);
    if(!is_good) return -1;
    wave = &(smoothed[0]);

    //let's save the sum waveform for checks
    if(MessageHandler::GetInstance()->GetDefaultMessageThreshold()<INFO){
      //      std::cout<<"DefaultMessageThreshold: "<<MessageHandler::GetInstance()->GetDefaultMessageThreshold()<<std::endl;
      char wfm_fname[64];
      sprintf(wfm_fname, "waveform_%s.txt",chdata->label.c_str());
      std::remove(wfm_fname);
      std::ofstream ofs (wfm_fname, std::ofstream::out);
      for(size_t ab=0; ab<summed.size(); ab++)
	ofs <<chdata->SampleToTime(ab)<<'\t'<<summed.at(ab)<<'\t'<<smoothed.at(ab)<<std::endl;
      ofs.close();
    }//end if messagehandler
  }//end if filter nsamps>1
  
  //here we do a simple threshold discriminator search
  //find pulses that cross the absolute threshold
  std::vector<peak_t> peaks;
  int pulse_start_add_nsamps = pulse_start_add_us*chdata->sample_rate+1;
  int pulse_end_add_nsamps = pulse_end_add_us*chdata->sample_rate+1;
  int search_result = DiscriminatorSearch(chdata, wave, peaks, threshold, 
					  pulse_start_add_nsamps, pulse_end_add_nsamps);
  if(search_result) return search_result;	
  
  //some pulses overlap and not cross the absolute threshold
  //let's find the pulses that cross the relative threshold
  std::vector<peak_t> peaks_split;
  std::vector<peak_t> peaks_tmp;
  for(size_t index=0; index<peaks.size(); index++){
    //     std::cout<<"Pulse to resolve pileup: "<<chdata->channel_id<<'\t'<<index<<'\t'<<chdata->SampleToTime(peaks.at(index).startIndex)
    // 	      <<'\t'<<chdata->SampleToTime(peaks.at(index).endIndex)<<std::endl;
    //here we use the same threshold but use the half filter length because of +/-
    int split_result= PileUpPulses(chdata, wave,peaks_tmp, peaks.at(index),
				   threshold, filter_nsamps/2);
    if(split_result || peaks_tmp.size()==0) continue;
    peaks_split.insert(peaks_split.end(), peaks_tmp.begin(), peaks_tmp.end());
  }//end for index  

  //due to the way we define the summed waveform, there may be an offset between raw waveform minima
  //and summed waveform minima
  double offset_nsamps = filter_nsamps/(-2);
  
  //here we push the pulses to the channel data
  chdata->pulses.reserve(peaks_split.size());
  for (size_t i = 0; i < peaks_split.size();  i++){
    //     std::cout<<chdata->channel_id<<'\t'<<i<<'\t'<<chdata->SampleToTime(peaks_split.at(i).startIndex)
    // 	     <<'\t'<<chdata->SampleToTime(peaks_split.at(i).endIndex)<<std::endl;

    //if this pulse shares boundary with the next pulse
    if(offset_nsamps && i+1<peaks_split.size() && //we have a next pulse
       peaks_split.at(i).endIndex+1>=peaks_split.at(i+1).startIndex && //there is an overlap between pulses
       peaks_split.at(i+1).peakIndex>peaks_split.at(i+1).startIndex && //next pulse has a valid peak index
       FirstAmplitudeIsSmaller(0,wave[peaks_split.at(i+1).startIndex],wave[peaks_split.at(i+1).peakIndex]) ){ //double sure overlap
      double offset_fraction = wave[peaks_split.at(i+1).startIndex]/wave[peaks_split.at(i+1).peakIndex];
      //this is an empirical function that seems to work for thr ArTPC data
      offset_fraction = sqrt(offset_fraction*(2.-offset_fraction));
      //      std::cout<<"Overlap time corrected: "<<-0.5*filter_time*offset_fraction<<std::endl;
      AddOffsetWithBounds<int>(peaks_split.at(i).endIndex, offset_nsamps*offset_fraction, 0, chdata->nsamps-1);
      AddOffsetWithBounds<int>(peaks_split.at(i+1).startIndex, offset_nsamps*offset_fraction, 0, chdata->nsamps-1);
      //    AddOffsetWithBounds<int>(peaks_split.at(i).peakIndex, offset_nsamps*offset_fraction, 0, chdata->nsamps-1);
    }//end if

    Pulse pulse;
    pulse.start_index = peaks_split.at(i).startIndex;
    //    pulse.peak_index = peaks_split.at(i).peakIndex;
    pulse.end_index = peaks_split.at(i).endIndex;
    chdata->pulses.push_back(pulse);
  } // end for loop over pulses
  chdata->npulses = chdata->pulses.size();

  return 0;
}


int PulseFinder::FindChannelSpikes(ChannelData* chdata){

  //we simple don't process the sum channel for spikes finding
  if(chdata->channel_id == ChannelData::CH_SUM) return 0;
  
  //let's find the threshold for spikes
  double threshold = 0;
  id_map::iterator it = ch_spike_thresholds.find(chdata->channel_id);
  if(it != ch_spike_thresholds.end()){
    if( it->second!=0) threshold = it->second;
    else return 0;
  }//end if ch_spike_thresholds is specified
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
    }//end else
  }//end if relative to bls
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

  std::vector<peak_t> peaks;
  const double * wave = chdata->GetBaselineSubtractedWaveform();
  int search_result = DiscriminatorSearch(chdata, wave, peaks, threshold,
					  spike_edge_add_nsamps, spike_edge_add_nsamps);
  if(search_result) return search_result;	

  //now let's split the spikes that don't cross the absolute threshold
  //but cross the relative threshold
  std::vector<peak_t> peaks_split;
  std::vector<peak_t> peaks_tmp;
  for(size_t index=0; index<peaks.size(); index++){
    //     std::cout<<chdata->channel_id<<'\t'<<index<<'\t'<<chdata->SampleToTime(start_index.at(index))
    // 	     <<'\t'<<chdata->SampleToTime(end_index.at(index))<<std::endl;
    int split_result= PileUpSpikes(chdata, wave, peaks_tmp, peaks.at(index), threshold);
    if(split_result) continue;
    peaks_split.insert(peaks_split.end(), peaks_tmp.begin(), peaks_tmp.end());
  }//end for index  
  
  //evaluate the spikes and push to vector
  const double * integral = chdata->GetIntegralWaveform();
  chdata->spikes.reserve(peaks_split.size());
  for (size_t i = 0; i < peaks_split.size();  i++){
    Spike spike;
    spike.start_time = chdata->SampleToTime(peaks_split.at(i).startIndex);
    spike.peak_time = chdata->SampleToTime(peaks_split.at(i).peakIndex);
    spike.width = (peaks_split.at(i).endIndex-peaks_split.at(i).startIndex)/chdata->sample_rate;
    spike.peak_amplitude = std::fabs(wave[peaks_split.at(i).peakIndex]);  
    spike.integral = integral[peaks_split.at(i).endIndex]-integral[peaks_split.at(i).startIndex];
    chdata->spikes.push_back(spike);
  } // end for loop over spikes
  
  return 0;
}

/// Search for pulses using a simple discrimination threshold
int PulseFinder::DiscriminatorSearch(ChannelData* chdata, const double * wave,
				     std::vector<peak_t>& peaks, double threshold,
				     int pulse_start_add_nsamps, int pulse_end_add_nsamps){
  
  if(!chdata || !wave || threshold==0 || pulse_start_add_nsamps<0 || pulse_end_add_nsamps<0)
    return -1;

  bool spike_finding=true;
  if(pulse_start_add_nsamps>spike_edge_add_nsamps || 
     pulse_end_add_nsamps>spike_edge_add_nsamps) spike_finding=false;
  
  //  std::cerr<<"Searching for pulses in channel: "<<chdata->channel_id<<std::endl;
  if(peaks.size()) peaks.clear();
  int search_start_index = chdata->TimeToSample(search_start_time, true);
  int search_end_index = chdata->TimeToSample(search_end_time, true);

  for(int index = search_start_index; index<=search_end_index; index++){
    //we are beyond threshold
    if(FirstAmplitudeIsSmaller(threshold,wave[index], threshold)){
      //if just come to a pulse
      if(index==0 || !(FirstAmplitudeIsSmaller(threshold, wave[index-1],threshold))){
	peak_t apeak;
	apeak.startIndex=index;
	apeak.peakIndex = index;
	apeak.endIndex = -1;
	peaks.push_back(apeak);
	//	std::cerr<<"Find pulse start: "<<chdata->SampleToTime(index)<<std::endl;
	//Do NOT put continue here because this could also be the existing sample 
      }//end if just come to a pulse
      if(peaks.size() && peaks.back().endIndex==-1){
	if(FirstAmplitudeIsSmaller(wave[peaks.back().peakIndex],wave[index],threshold))
	  peaks.back().peakIndex=index;
	//if we are about to exit a pulse
	if(index==chdata->nsamps-1 || !FirstAmplitudeIsSmaller(threshold,wave[index+1], threshold)){
	  peaks.back().endIndex=index;
	  //	  std::cerr<<"Find pulse end: "<<chdata->SampleToTime(index)<<std::endl;
	}//end if exit a pulse
      }//end if there is an open pulse
    }//end if inside a pulse
  }//end for int index
  if(peaks.size() && peaks.back().endIndex==-1) peaks.back().endIndex=search_end_index;

  //resoving pulses that overlap due to the added samples on both ends
  int total_nsamps_add = pulse_start_add_nsamps+pulse_end_add_nsamps;
  for(size_t index=0; index<peaks.size(); index++){
    if(peaks.at(index).startIndex>peaks.at(index).endIndex){
      Message(ERROR)<<"PulseFinder::DiscriminatorSearch failed for channel "<<chdata->channel_id<<"! \n"
		    <<peaks.at(index).startIndex<<'\t'<<peaks.at(index).endIndex<<". \n";
      peaks.erase(peaks.begin()+index);
      index--;
      continue;
    }
    //     std::cout<<"Overlap? "<<chdata->channel_id<<'\t'<<index<<'\t'<<chdata->SampleToTime(peaks.at(index).startIndex)
    //      	     <<'\t'<<chdata->SampleToTime(peaks.at(index).endIndex)<<std::endl;
    if(index==0){//deal with the start index for the first one
      if(peaks.at(index).startIndex-pulse_start_add_nsamps>=search_start_index)
	peaks.at(index).startIndex-=pulse_start_add_nsamps;
      else peaks.at(index).startIndex=search_start_index;
    }
    if(index+1==peaks.size()){//deal with the end idnex for the last one
      if(peaks.at(index).endIndex+pulse_end_add_nsamps<=search_end_index)
	peaks.at(index).endIndex+=pulse_end_add_nsamps;
      else peaks.at(index).endIndex=search_end_index;
    }
    if(index<1) continue;
    if(peaks.at(index).startIndex-peaks.at(index-1).endIndex>total_nsamps_add){
      peaks.at(index).startIndex-=pulse_start_add_nsamps;
      peaks.at(index-1).endIndex+=pulse_end_add_nsamps;
      continue;
    }
    //if we get to here we have a possible overlap
    //leave it to the pileup function to resolve it
    if(!spike_finding){
      peaks.at(index-1).endIndex=peaks.at(index).endIndex;
      if(FirstAmplitudeIsSmaller(wave[peaks.at(index-1).peakIndex],wave[peaks.at(index).peakIndex],threshold))
	peaks.at(index-1).peakIndex=peaks.at(index).peakIndex;
      peaks.erase(peaks.begin()+index);
      index--;
      continue;
    }
    else{//resolve it for spike counting
      int middle_index = peaks.at(index).startIndex;
      double min_value = wave[middle_index];
      for(int j=peaks.at(index-1).endIndex; j<peaks.at(index).startIndex; j++){
	if(FirstAmplitudeIsSmaller(wave[j], min_value, threshold)){
	  min_value = wave[j];
	  middle_index = j;
	}//end if
      }//end for j loop
      //      std::cout<<"found minimum: "<<chdata->SampleToTime(middle_index)<<std::endl;
      peaks.at(index).startIndex = middle_index;
      peaks.at(index-1).endIndex = middle_index;
    }//end else

  }//end for index loop

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

/// resolve pileup pulses/spikes
/// only trust the result if  within a pulse/spike
/// may pickup baseline fluctuations otherwise
int PulseFinder::PileUpPulses(ChannelData* chdata, const double * wave,
			      std::vector<peak_t>& peaks,  peak_t range, 
			      double threshold, int step){
  
  if(!chdata || !wave || step<1 || threshold==0 || range.startIndex<0 || 
     range.endIndex>=chdata->nsamps || range.startIndex>=range.endIndex) return -1;
  if(peaks.size()) peaks.clear();
  peaks.reserve(10);  
  //no point searching pileups if it is a small pulse
  if(range.endIndex-range.startIndex<=step){
    peaks.push_back(range);
    return 0;
  }
  //  std::cout<<"Resolving overlaps for pulse: "<<chdata->SampleToTime(range.startIndex)<<'\t'<<chdata->SampleToTime(range.endIndex)<<std::endl;

  const int nstep=6;
  int steps[nstep] = {-step, step*2/(-3), step/(-3), step/3, step*2/3, step};

  for(int samp=range.startIndex+1; samp<=range.endIndex-1; samp++){
    //here we come to a possible local peak -- not smaller than local neighbours
    if(!FirstAmplitudeIsSmaller(wave[samp],wave[samp-1],threshold) &&
       !FirstAmplitudeIsSmaller(wave[samp],wave[samp+1],threshold) ){
      bool is_real_peak = true, good_on_left=false, good_on_right=false;
      for(int ii=0; ii<nstep; ii++){
	if(samp+steps[ii]>=0 && samp+steps[ii]<chdata->nsamps){
	  if(FirstAmplitudeIsSmaller(wave[samp],wave[samp+steps[ii]],threshold))
	    is_real_peak = false;
	  else if(FirstAmplitudeIsSmaller(wave[samp+steps[ii]],wave[samp],threshold)){
	    if(ii<nstep/2) good_on_left=true;
	    else good_on_right=true;
	  }
	}//end if within reange
      }//end for ii
      if(is_real_peak && good_on_left && good_on_right){
	//	std::cout<<"Possible peak at: "<<chdata->SampleToTime(samp)<<std::endl;
	if(peaks.size() && peaks.back().endIndex==-1){//if there is a peak not close
	  if(FirstAmplitudeIsSmaller(wave[peaks.back().peakIndex], wave[samp], threshold))
	    peaks.back().peakIndex=samp;
	}//end there is a peak not closed
	else{
	  peak_t apeak;
	  if(peaks.size()) apeak.startIndex=peaks.back().endIndex;
	  else apeak.startIndex=range.startIndex;
	  apeak.peakIndex=samp;
	  apeak.endIndex=-1;
	  peaks.push_back(apeak);
	}//end else
      }//end if is_real_peak
    }//end if a local peak
    //don't put else here, here we come to a possible local valley
    if(!FirstAmplitudeIsSmaller(wave[samp-1],wave[samp],threshold) &&
       !FirstAmplitudeIsSmaller(wave[samp+1],wave[samp],threshold) ){
      bool is_real_valley = true, good_on_left=false, good_on_right=false;
      for(int ii=0; ii<nstep; ii++){
	if(samp+steps[ii]>=0 && samp+steps[ii]<chdata->nsamps){ 
	  if(FirstAmplitudeIsSmaller(wave[samp+steps[ii]],wave[samp],threshold))
	    is_real_valley = false;
	  else if(FirstAmplitudeIsSmaller(wave[samp],wave[samp+steps[ii]],threshold)){
            if(ii<nstep/2) good_on_left=true;
            else good_on_right=true;
          }
	}//end if within range
      }//end for ii
      if(is_real_valley && good_on_left && good_on_right){
	//	std::cout<<"Possible valley at: "<<chdata->SampleToTime(samp)<<std::endl;
	if(peaks.size() && 
	   (peaks.back().endIndex==-1 ||
	    FirstAmplitudeIsSmaller(wave[samp], wave[peaks.back().endIndex], threshold)) )
	  peaks.back().endIndex=samp;
      }//end if real valley
    }//end if a local valley
  }//end for loop

  //we need to make sure the last peak always ends on the range.endIndex
  if(peaks.size()) peaks.back().endIndex=range.endIndex;
  
  //only keep pulses cross the threshold
  int peaks_merged = 0;
  do{
    peaks_merged = 0;
    for(size_t jj=0; jj<peaks.size(); jj++){
      if(peaks.at(jj).startIndex>=peaks.at(jj).endIndex){//specifically covers endIndex==-1
	Message(ERROR)<<"PulseFinder::PileUpPulses failed for channel "<<chdata->channel_id<<"! \n";
	peaks.erase(peaks.begin()+jj);
	jj--;//this should be okay for size_t 0
	continue;
      }
      //this peak doesn't cross the threshold
      if(!(RelativeThresholdCrossed(wave[peaks.at(jj).startIndex],wave[peaks.at(jj).peakIndex],threshold)) ||
	 !(RelativeThresholdCrossed(wave[peaks.at(jj).endIndex],wave[peaks.at(jj).peakIndex],threshold))){
	Message(DEBUG)<<"Small pulse in channel "<<chdata->channel_id<<" will be merged: "
		     <<chdata->SampleToTime(peaks.at(jj).startIndex)<<", "
		     <<chdata->SampleToTime(peaks.at(jj).endIndex)<<". \n";
	bool merge_left=false, merge_right=false;
	//merge left in this case
	if(jj>0 && (!FirstAmplitudeIsSmaller(wave[peaks.at(jj).startIndex], wave[peaks.at(jj).endIndex], threshold) ||
		    jj+1 == peaks.size())){
	  peaks.at(jj-1).endIndex=peaks.at(jj).endIndex;
	  if(FirstAmplitudeIsSmaller(wave[peaks.at(jj-1).peakIndex],wave[peaks.at(jj).peakIndex],threshold))
	    peaks.at(jj-1).peakIndex = peaks.at(jj).peakIndex;
	  merge_left = true;
	}//end if merge left
	//merge right in this case
	if(!merge_left && jj+1<peaks.size()){
	  peaks.at(jj+1).startIndex=peaks.at(jj).startIndex;
	  if(FirstAmplitudeIsSmaller(wave[peaks.at(jj+1).peakIndex],wave[peaks.at(jj).peakIndex],threshold))
	    peaks.at(jj+1).peakIndex = peaks.at(jj).peakIndex;
	  merge_right = true;
	}//end if merge right
	if(merge_left || merge_right){
//  	  std::cout<<"Peak merged "<<chdata->SampleToTime(peaks.at(jj).startIndex)<<", "
// 		   << chdata->SampleToTime(peaks.at(jj).endIndex)<<" to the left? "<<merge_left<<std::endl;
	  peaks_merged ++;
	  peaks.erase(peaks.begin()+jj);
	  jj--;//this should be okay for size_t 0
	}
      }//end if peak is too small
    }//end for jj loop
  }while(peaks_merged>0);

  //if nothing cross the relative threshold we still want to keep the original pulse
  if(peaks.size()==0) peaks.push_back(range);
  
  return 0;
}

/// resolve pileup spikes only trust the result if  within a pulse/spike
/// may pickup baseline fluctuations otherwise
int PulseFinder::PileUpSpikes(ChannelData* chdata, const double * wave,
			      std::vector<peak_t>& peaks, peak_t range,  double threshold){
  
  if(!chdata || !wave || threshold==0 || range.startIndex<0 || 
     range.endIndex>=chdata->nsamps || range.startIndex>=range.endIndex) return -1;
  if(peaks.size()) peaks.clear();
  peaks.reserve(100);  
  const int step = 2; //this is the best value for spike counting
  //no point searching pileups if it is a small pulse
  if(range.endIndex-range.startIndex<=step){
    peaks.push_back(range);
    return 0;
  }
  //  std::cout<<"Resolving overlaps for spike: "<<chdata->SampleToTime(range.startIndex)<<'\t'<<chdata->SampleToTime(range.endIndex)<<std::endl;

  for(int samp=range.startIndex+1; samp<=range.endIndex-1; samp++){
    //here we come to a possible peak, w/ bigger amplitude than +/-1, or +/-step
    if( (RelativeThresholdCrossed(wave[samp-1],wave[samp],threshold) ||
	 (samp>=range.startIndex+step && RelativeThresholdCrossed(wave[samp-step],wave[samp],threshold)) ) &&
	(RelativeThresholdCrossed(wave[samp+1],wave[samp],threshold) ||
	 (samp+step<=range.endIndex && RelativeThresholdCrossed(wave[samp+step],wave[samp],threshold))) ){
      //      std::cout<<"Possible peak at: "<<chdata->SampleToTime(samp)<<std::endl;
      if(peaks.size() && peaks.back().endIndex==-1){//if there is a peak not closed
	if(FirstAmplitudeIsSmaller(wave[peaks.back().peakIndex], wave[samp], threshold))
	  peaks.back().peakIndex=samp;
      }//end there is a peak not closed
      else{
	peak_t apeak;
	if(peaks.size()) apeak.startIndex=peaks.back().endIndex;
	else apeak.startIndex=range.startIndex;
	apeak.peakIndex=samp;
	apeak.endIndex=-1;
	peaks.push_back(apeak);
      }//end else
    }//end if coming to a peak
    //here we come to a valley
    else if((RelativeThresholdCrossed(wave[samp],wave[samp-1],threshold) ||
	     (samp>=range.startIndex+step && RelativeThresholdCrossed(wave[samp],wave[samp-step],threshold)) ) &&
	    (RelativeThresholdCrossed(wave[samp],wave[samp+1],threshold) ||
	     (samp+step<=range.endIndex && RelativeThresholdCrossed(wave[samp],wave[samp+step],threshold))) ){
      //      std::cout<<"Possible valley at: "<<chdata->SampleToTime(samp)<<std::endl;
      if(peaks.size() && 
	 (peaks.back().endIndex==-1 ||
	  FirstAmplitudeIsSmaller(wave[samp], wave[peaks.back().endIndex], threshold)) )
	peaks.back().endIndex=samp;
    }//end if coming to a valley
  }//end for samp
  //we need to make sure the first/last peak is complete
  if(peaks.size()) peaks.back().endIndex = range.endIndex;

  //make sure all spikes are valid
  for(size_t jj=0; jj<peaks.size(); jj++){
    if(peaks.at(jj).startIndex>=peaks.at(jj).endIndex){
      Message(ERROR)<<"PulseFinder::PileUpSpikes failed for channel "<<chdata->channel_id<<"! \n";
      peaks.erase(peaks.begin()+jj);
      jj --;
    continue;
    }
  }//end for

  if(peaks.size()==0) peaks.push_back(range);
  return 0;  
}

