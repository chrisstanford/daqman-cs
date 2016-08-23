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

    //Start search for pulse edges (start and end) *
    std::vector<int> start_index;
    std::vector<int> end_index;
    
    //Loop over all channels and evaluate pulse edges individually
    for (size_t ch = 0; ch < event->channels.size(); ch++)
      {
	ChannelData& chdata = event->channels[ch];
	//skip channels we've been told to explicitly skip
	if(_skip_channels.find(chdata.channel_id) != _skip_channels.end())
	  continue;
	//make sure that we check baseline for each channel
	if(! chdata.baseline.found_baseline)
	  continue;
	start_index.clear();
	end_index.clear();

	int start_add_nsamps = pulse_start_add_us*chdata.sample_rate+1;
	int end_add_nsamps = pulse_end_add_us*chdata.sample_rate+1;
	
	//	if(mode == DISCRIMINATOR)
	DiscriminatorSearch(&chdata, chdata.GetBaselineSubtractedWaveform(), 
			    start_index, end_index, start_add_nsamps, end_add_nsamps);
	
	int search_start_index=chdata.TimeToSample(search_start_time, true);
	for (size_t i = 0; i < start_index.size();  i++)
	{
	  if (start_index[i] >= end_index[i]) {
	    //	    return -1;
	    continue;
	  }
	    Pulse pulse;
	    pulse.start_index = start_index[i];
	    pulse.end_index = end_index[i];
	    EvaluatePulse(pulse, &chdata, 0.2);
	    //evaluate the pulse parameters depending on other pulses
	    if(chdata.pulses.size()){
	      pulse.time_since_last=pulse.half_max_time-chdata.pulses.back().end_time;
	      if(pulse.start_index<=chdata.pulses.back().end_index)
		{pulse.is_clean = false; chdata.pulses.back().is_clean = false;}
	    }
	    else{
	      pulse.time_since_last=pulse.half_max_time-chdata.SampleToTime(search_start_index);
	      if(pulse.start_index<=search_start_index) pulse.is_clean = false;
	    }
	    //  pulse.Print(chdata.channel_id, chdata.pulses.size());
	    chdata.pulses.push_back(pulse);
	    
	} // end for loop over pulses
	/*
	//get rid of noise pulses
	//	double noise_start=0, noise_end=0;
	for(size_t j=0; j<chdata.pulses.size(); j++){
	  Pulse &pulse = chdata.pulses.at(j);

	  if(pulse.integral+pulse.peak_amplitude>=0 ||
	     //	     pulse.overshoot>5*chdata.baseline.sigma ||
	     pulse.overshoot>0.5*pulse.peak_amplitude){
	    pulse.evaluated=false;
	    / *	    noise_start = pulse.start_time-0.2;
	    noise_end = pulse.end_time+0.2;
	    //get rid of the previous noise pulses
	    for(int k=j-1; k>=0; k--){
	      if(chdata.pulses.at(k).start_time>=noise_start &&
		 chdata.pulses.at(k).start_time<=noise_end) chdata.pulses.at(k).evaluated=false;
	      else break;
	      }* /
	  }//end if this is a bad pulse
	  / *	  if(noise_end>noise_start){//if this pulse comes close after a noise
	    if(chdata.pulses.at(j).start_time>=noise_start &&
	       chdata.pulses.at(j).start_time<=noise_end) chdata.pulses.at(j).evaluated=false;
	    else {noise_start=0; noise_end=0;}//out of noise window
	    }//end if noise end* /
	}//end for loop
	*/
	chdata.npulses = chdata.pulses.size();
      } //end loop over channels
    //End evaluation of pulse variables for each pulse on each channel
    return 0;
}


/// resolve pileup pulses/spikes
/// only trust the result if  within a pulse/spike
/// may pickup baseline fluctuations otherwise
int PulseFinder::ResolvePileUps(ChannelData* chdata, const double * wave,
				std::vector<int>& start_index, 
				std::vector<int>& end_index,
				double pileup_threshold, int step, 
				int search_start, int search_end)
{return 0;}

/// Search for pulses using a simple discrimination threshold
int PulseFinder::DiscriminatorSearch(ChannelData* chdata, const double * wave,
				     std::vector<int>& start_index, 
				     std::vector<int>& end_index,
				     int start_add_nsamps, int end_add_nsamps){

  if(!wave) return -1;
  
  double check_val = 0;
  id_map::iterator it = ch_pulse_thresholds.find(chdata->channel_id);
  if(it != ch_pulse_thresholds.end())
    check_val = it->second;
  else{
    it = ch_pulse_thresholds.find(ChannelData::CH_DEFAULT);
    if(it != ch_pulse_thresholds.end()) check_val =it->second;
    else return 0; //we simply don't process this channel
  }//end else

  if(relative_bls_threshold){
    if(chdata->baseline.sigma>0) check_val *= chdata->baseline.sigma;
    else{
      Message(WARNING)<<"Channel "<<chdata->channel_id<< " has invalid baseline!\n";
      return 0;
    }
  }

  if(relative_spe_threshold){
    if(chdata->spe_mean>0){
      check_val *= chdata->spe_mean;
      if(chdata->spe_mean==1){
	Message(WARNING)<<"Channel "<<chdata->channel_id<< " uses spe_mean unspecified!\n";
      }//end if ==1
    }//end if >0
    else{
      Message(ERROR)<<"Channel "<<chdata->channel_id<< " has incorrectly specified spe_mean!\n";
      return -1;
    }//end else
  }//end relative spe

  if(start_index.size()) start_index.clear();
  if(end_index.size()) end_index.clear();
  bool in_pulse = false;
  int search_start_index = chdata->TimeToSample(search_start_time, true);
  int search_end_index = chdata->TimeToSample(search_end_time, true);

  for(int index = search_start_index; index<=search_end_index; index++){
    //we are beyond threshold
    if(FirstAmplitudeIsSmaller(check_val,wave[index], check_val)){
      //if just come to a pulse
      //if(!in_pulse){
      if(index==0 || FirstAmplitudeIsSmaller(wave[index-1],check_val, check_val)){
	in_pulse = true;
	if(index>=search_start_index+start_add_nsamps) start_index.push_back(index-start_add_nsamps );
	else start_index.push_back(search_start_index);
      }//end if just come to a pulse
      //if we are about to exit a pulse
      if(index==chdata->nsamps-1 || FirstAmplitudeIsSmaller(wave[index+1], check_val,check_val)){
	if(in_pulse){
	  if(index+end_add_nsamps<=search_end_index) end_index.push_back(index+end_add_nsamps);
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
      for(int j=end_index.at(index-1)+1-end_add_nsamps; j<start_index.at(index)+start_add_nsamps; j++){
	if(j<0) continue;
	else if(j>=chdata->nsamps) break;
	if(FirstAmplitudeIsSmaller(min_value,wave[j], check_val)) continue;
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
  /*  
    std::vector<double>& baseform = chdata->subtracted_waveform;
  
    for(size_t index=1; index<start_index_tmp.size(); index++){
      std::cout<<chdata->channel_id<<'\t'<<index<<'\t'<<chdata->SampleToTime(start_index_tmp.at(index))
	       <<'\t'<<chdata->SampleToTime(end_index_tmp.at(index))<<std::endl;
      if(start_index_tmp.at(index)>=end_index_tmp.at(index)) continue;
      std::vector<int> start_index_tmp2;
      std::vector<int> end_index_tmp2;
      int test_result= RelativeThresholdSearch(baseform, -1.*std::fabs(discriminator_nsigma)*chdata->baseline.sigma,
					       std::fabs(discriminator_nsigma)*chdata->baseline.sigma,
					       start_index_tmp2, end_index_tmp2, discriminator_end_add, 1,
					       start_index_tmp.at(index),end_index_tmp.at(index));
      if(!test_result){
	std::cout<<"*** Great *** RelativeThresholdSearch succeeded! "<<std::endl;
	start_index.insert(start_index.end(), start_index_tmp2.begin(), start_index_tmp2.end());
	end_index.insert(end_index.end(), end_index_tmp2.begin(), end_index_tmp2.end());
      }

    }  
  */
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

