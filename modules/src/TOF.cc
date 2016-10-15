#include "TOF.hh"
#include "BaselineFinder.hh"
#include "EvalPulse.hh"
#include "PulseFinder.hh"
#include "EventHandler.hh"
#include "Integrator.hh"
#include "ConvertData.hh"
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cmath>

TOF::TOF() : 
  ChannelModule(GetDefaultName(), "Search for & combine pulses for time of flight analysis")
{
  AddDependency<ConvertData>();
  AddDependency<BaselineFinder>();
  AddDependency<Integrator>();
  AddDependency<PulseFinder>();

  RegisterParameter("ref_ch", ref_ch = -100 , "reference channel number");
  RegisterParameter("signal_ch", signal_ch = -100 , "signal channel number");
  RegisterParameter("ref_ch_offset", ref_ch_offset = 0, 
		    "time offset in us for the reference channel");
  RegisterParameter("signal_start_time", signal_start_time = -0.1 , 
		    "start time in us to look for signals");
  RegisterParameter("neutron_start_time", neutron_start_time = -0.1 , 
		    "start time in us to look for neutron signals");
  RegisterParameter("search_end_time", search_end_time = 0.1 , 
		    "end time in us to search for pulse rises");
  RegisterParameter("pulse_length_us", pulse_length_us = 1.,
		    "time in us for pulses to be combined for analysis");
  //  RegisterParameter("", ,"");
  
}

TOF::~TOF()
{
}

int TOF::Initialize()
{
  return 0;
}

int TOF::Finalize()
{
  return 0;
}

int TOF::Process(ChannelData* chdata)
{
  if(chdata->channel_id==ref_ch) ProcessRefChannel(chdata);
  else ProcessDataChannel(chdata);
  return 0;
}
				
int TOF::ProcessRefChannel(ChannelData* chdata){

  std::vector<Pulse> &pulses = chdata->pulses;
  if(!pulses.size()) return 0;

  int index_before_signal=-1, index_after_signal=-1;
  double time_before_signal=chdata->SampleToTime(0)-1;
  double time_after_signal=chdata->SampleToTime(chdata->nsamps)+1;
  for(int pulse_index = 0; pulse_index < chdata->npulses-1; pulse_index++){
    if(pulses.at(pulse_index).peak_amplitude<500) continue;
    if((pulses.at(pulse_index).half_max_time < ref_ch_offset) && 
       (pulses.at(pulse_index).half_max_time>time_before_signal)){
      time_before_signal = pulses.at(pulse_index).half_max_time;
      index_before_signal = pulse_index;
    }
    else if((pulses.at(pulse_index).half_max_time>=ref_ch_offset) && 
	    (pulses.at(pulse_index).half_max_time<time_after_signal)){
      time_after_signal = pulses.at(pulse_index).half_max_time;
      index_after_signal = pulse_index;
      break;
    }
  }
  if(index_before_signal>=0){
    chdata->tof.push_back(pulses.at(index_before_signal));
    chdata->tof.back().time_since_last = time_after_signal - time_before_signal;
  }
  if(index_after_signal>=0){
    chdata->tof.push_back(pulses.at(index_after_signal));
    chdata->tof.back().time_since_last = time_after_signal - time_before_signal;
  }

  return 0;
}

int TOF::ProcessDataChannel(ChannelData* chdata){

  std::vector<Pulse> &pulses = chdata->pulses;
  if(!pulses.size()) return 0;

  double search_start_time = neutron_start_time;
  if(chdata->channel_id == signal_ch) search_start_time = signal_start_time;
  
  Pulse pulse;
  //use evaluated as a temporary flag
  pulse.evaluated = false;
  int pulse_nsamps = 0;
  //use the start time of the first pulse
  //use summed pulse integral, hope noise can be rejected as bad pulses
  double half_max_time = -1, integral=0;
  for(int pulse_index = 0; pulse_index < chdata->npulses; pulse_index++){
    //    pulses.at(pulse_index).Print(chdata->channel_id, pulse_index);
    if(pulses.at(pulse_index).half_max_time<search_start_time) continue;
    //    else if(pulses.at(pulse_index).half_max_time<search_end_time) break;
    //    std::cout<<"\tPulse added\t"<<chdata->channel_id<<'\t'<<pulse_index<<std::endl;
    if(!pulse.evaluated){
      //no valid pulse within the signal region give it up
      if(pulses.at(pulse_index).half_max_time>search_end_time) break;
      //only record the basic information
      // std::cout<<"First Pulse in TOF: "<<pulses.at(pulse_index).start_time
      // 	       <<'\t'<<pulses.at(pulse_index).half_max_time<<std::endl;
      half_max_time = pulses.at(pulse_index).half_max_time;
      integral = pulses.at(pulse_index).integral;
      pulse.evaluated = true;
      pulse.start_index = pulses.at(pulse_index).start_index;
      pulse.end_index = pulses.at(pulse_index).end_index;
      pulse_nsamps += pulses.at(pulse_index).end_index-pulses.at(pulse_index).start_index;
      //use time_since_last to store the clean baseline information
      if(pulse_index>0) pulse.time_since_last = pulses.at(pulse_index).half_max_time-pulses.at(pulse_index-1).end_time;
      else pulse.time_since_last = pulses.at(pulse_index).half_max_time - chdata->SampleToTime(0);
    }
    else if (pulses.at(pulse_index).half_max_time<half_max_time+pulse_length_us){
      pulse_nsamps += pulses.at(pulse_index).end_index-pulses.at(pulse_index).start_index;
      pulse.end_index = pulses.at(pulse_index).end_index;
      integral += pulses.at(pulse_index).integral;
    }
    else break;
  }//end for
  
  if(pulse.evaluated){//we found at least 1 valid pulse in this region
    EvaluatePulse(pulse,chdata, 0.2);
    pulse.half_max_time = half_max_time;
    //    if(integral<0) std::cout<<pulse.integral-integral<<"\t,percentage="<<pulse.integral/integral<<std::endl;
    pulse.integral = integral;
    pulse.test     = pulse_nsamps;

    // //calculate the peaks -- this is only useful for NaI PSD studies
    // if(!(chdata->integral.empty())){
    //   for(int pulse_index = 0; pulse_index < chdata->npulses; pulse_index++){
    // 	//    pulses.at(pulse_index).Print(chdata->channel_id, pulse_index);
    // 	if(pulses.at(pulse_index).half_max_time<pulse.start_time) continue;
    // 	else if(pulses.at(pulse_index).end_time>pulse.end_time) break;
    // 	if(pulses.at(pulse_index).peak_amplitude>50 //here we look at pulses larger than 2 p.e. !!!
    // 	   || pulses.at(pulse_index).start_index==pulse.start_index){//or this is the first pulse in this TOF
    // 	  Peak pk;
    // 	  pk.peak_amplitude = pulses.at(pulse_index).peak_amplitude;
    // 	  pk.peak_time = pulses.at(pulse_index).peak_time;
    // 	  int pk_start_index = chdata->TimeToSample(pk.peak_time-0.05,true);//this is 2xNaI scintillation lifetime.
    // 	  int pk_end_index = chdata->TimeToSample(pk.peak_time+0.25,true);//this is 2xNaI scintillation lifetime.
    // 	  pk.integral = chdata->integral[pk_end_index]-chdata->integral[pk_start_index];
    // 	  pulse.peaks.push_back(pk);
    // 	}//end if pulses
    //   }//end for pulse index
    // }//end if integral
    
    chdata->tof.push_back(pulse);
    //    pulse.Print(chdata->channel_id, -1);
  }
  
  return 0;

}
