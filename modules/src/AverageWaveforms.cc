#include "AverageWaveforms.hh"
#include "ConvertData.hh"
#include "BaselineFinder.hh"
#include "EvalRois.hh"
#include "TGraphErrors.h"
#include "TFile.h"
#include "RootWriter.hh"
#include "EventData.hh"

#include <algorithm>
#include <numeric>
#include <iostream>

AverageWaveforms::AverageWaveforms() : 
  ChannelModule(GetDefaultName(), "Average the waveform for each channel over the entire run and save to the output root file")
{
  AddDependency<ConvertData>();
  AddDependency<BaselineFinder>();
  //we need at least one Roi to evaluate the pulse threshold time
  AddDependency<EvalRois>();
  //Register all the config handler parameters
  RegisterParameter("ref_roi_index", ref_roi_index = 0,
		    "Index of the ROI referred for pulse arrving time and averaging end time");
  RegisterParameter("aver_time_hr", aver_time_hr = 1,
		    "Maximum time duration in hours for one averaged waveform");
  RegisterParameter("pre_trigger_time", pre_trigger_time = 1,
		    "time before pulse rise to start summing waveforms");
  RegisterParameter("min_pulse_height", min_pulse_height = 0,
		    "Minimum pulse height for one waveform to be counted in average waveform");
  RegisterParameter("max_pulse_height", max_pulse_height = 1e100,
		    "Maximum pulse height for one waveform to be counted in average waveform");
  //warning: npe=1 if it is not read from database!!!
  //be careful when you apply the npe cut!!!!!!
  RegisterParameter("min_npe", min_npe = 0,
		    "Minimum npe for one waveform to be counted in average waveform");
  RegisterParameter("max_npe", max_npe = 1e100,
		    "Maximum npe for one waveform to be counted in average waveform");
  //set the pulse shape parameter values to be unphysical
  RegisterParameter("ref_fp_index", ref_fp_index = 0,
		    "Index of the Fprompt used for averaging waveform cuts");
  RegisterParameter("min_fprompt", min_fprompt = -100,
		    "Minimum fprompt for one waveform to be counted in average waveform");
  RegisterParameter("max_fprompt", max_fprompt = -100,
		    "Maximum fprompt for one waveform to be counted in average waveform");
  RegisterParameter("min_tof",  min_tof= -100, "min_tof");
  RegisterParameter("max_tof",  max_tof= -100, "max_tof");
  //do not register this parameter, will calculate later
  post_trigger_nsamps = -1;
}

AverageWaveforms::~AverageWaveforms()
{
  Cleanup();
}

void AverageWaveforms::Cleanup()
{
    std::map<int,TGraphErrors*>::iterator mapit = _plots.begin();
    for( ; mapit != _plots.end() ; mapit++){
	delete mapit->second;
    }
    _plots.clear();
    _sum.clear();
    _sum2.clear();
}

int AverageWaveforms::Initialize()
{ 
  if(aver_time_hr<0.01) aver_time_hr = 0.01;
  if(pre_trigger_time<=0.05) pre_trigger_time=0.05;
  return 0; 
}

int AverageWaveforms::Finalize()
{

  if(gFile && gFile->IsOpen()){
    std::map<int,TGraphErrors*>::iterator mapit = _plots.begin();
    for( ; mapit != _plots.end(); mapit++){
      TGraphErrors* graph = mapit->second;
      double* y_array = graph->GetY();
      double* ey_array = graph->GetEY();
      double a_sum  = _sum[mapit->first];
      double a_sum2 = _sum2[mapit->first];
      if(a_sum<=0)  a_sum=1;
      if(a_sum2<=0) a_sum2=1;
      for(int i=0; i < graph->GetN(); i++){
        y_array[i] /= a_sum;
        ey_array[i] = sqrt((ey_array[i]/a_sum-y_array[i]*y_array[i])*a_sum2)/a_sum;
      }
      (mapit->second)->SetFillStyle(3002);
      (mapit->second)->SetFillColor(kRed);
      (mapit->second)->Write();
    }
  }
  Cleanup();
  return 0;
}

//this is not a critical module, return 0 for all cases
int AverageWaveforms::Process(ChannelData *chdata){
  //we need at least one roi to know the threshold time to align pulses
  if(!(chdata->baseline.found_baseline) || chdata->saturated ||
     chdata->baseline.sigma>1.2 || chdata->regions.size()-ref_roi_index<=0)
    return 0;
  //  Pulse &pulse = chdata->regions.at(ref_roi_index);  
  Pulse pulse = chdata->regions.at(ref_roi_index);  
  if(!(pulse.evaluated) || pulse.integral>=0) return 0; //integral is negative
  //specify the end of the average, hopefully only do it once
  if(post_trigger_nsamps<0){
    post_trigger_nsamps = (pulse.end_time-pulse.half_max_time)*chdata->sample_rate;
    if(post_trigger_nsamps<=1){
      post_trigger_nsamps = -1; 
      return 0;
    }
  }

  //warning: add the cuts for TOF -- for NaI PSD only
  //only apply the cut when tof cuts are specified by the users
  if(max_tof>-100 && min_tof>-100){
    ChannelData *refch = _current_event->GetEventData()->GetChannelByID(7);
    if(!refch) return 0;
    //  std::cout<<chdata->channel_id<<": "<<chdata->tof.size()<<", "<<refch->tof.size()<<std::endl;
    if(chdata->tof.size()<1 || refch->tof.size()<2) return 0;
    //  std::cout<<'\t'<<chdata->tof[0].half_max_time-refch->tof[0].half_max_time<<std::endl;
    if(chdata->tof[0].time_since_last<0.2 || //psd parameters? 
       chdata->tof[0].half_max_time-refch->tof[0].half_max_time<min_tof ||
       chdata->tof[0].half_max_time-refch->tof[0].half_max_time>max_tof) return 0;
    pulse = chdata->tof.at(0);
    if(!(pulse.evaluated) || pulse.integral>=0) return 0; //integral is negative
    //  std::cout<<"candidate event!"<<std::endl;
    // ChannelData *nch = _current_event->GetEventData()->GetChannelByID(2);
    // if(!nch || nch->tof.size()<1) return 0;
    // //  std::cout<<"candidate event come here!"<<std::endl;
    // //this cut is for channel 5 in Run3
    // if(nch->tof[0].time_since_last<0.2 || nch->tof[0].peak_amplitude<100 || //nch->tof[0].t90<0.3 || //this can be dangerous
    //    nch->tof[0].half_max_time-refch->tof[0].half_max_time<-0.01 ||
    //    nch->tof[0].half_max_time-refch->tof[0].half_max_time>0.05) return 0;
     // std::cout<<"candidate event accepted!"<<std::endl;
     // std::cout<<'\t'<<chdata->tof[0].half_max_time-refch->tof[0].half_max_time<<std::endl;
  }

  //  get rid of pileup events
  if(pulse.min<1 || pulse.peak_amplitude>-0.08*pulse.integral
     //|| pulse.peak_amplitude<chdata->baseline.mean-pulse.min-0.5
     ) return 0;
  
  //apply amplitude and npe cut
  if(pulse.peak_amplitude>max_pulse_height || 
     pulse.peak_amplitude<min_pulse_height ||
     pulse.npe>max_npe || pulse.npe<min_npe)
    return 0;

  //PSD may not be necessary, so only apply the fprompt cut
  //when the cut values are specified by the users to be >-1 
  if(max_fprompt>-100 && min_fprompt>-100 &&
     (pulse.fparameters.size()-ref_fp_index<=0 ||
      pulse.fparameters.at(ref_fp_index)>max_fprompt ||
      pulse.fparameters.at(ref_fp_index)<min_fprompt))
    return 0;

  std::cout<<"average event: "<<_current_event->GetEventData()->event_id<<std::endl;

  int trigger_index=pulse.start_index;
  //  int trigger_index=chdata->TimeToSample(pulse.half_max_time+0.5/chdata->sample_rate);
  int pre_trigger_nsamps  = pre_trigger_time *chdata->sample_rate;
  const double* wave = chdata->GetBaselineSubtractedWaveform();

  //65536 is mainly to scale down the pulse integral close to unit
  //it need to be squared and summed, don't want the value to go out of dynamic range
  //65536=1024*64, or 64 samples with the value of 1024 ADC counts
  //the pulse also needs to be scaled down by 65536
  double pulse_integral = pulse.integral/(-65536);

  //calculate the average waveform
  //convention is to set time 0 to separate the bin lower than threshold
  //and the bin higher than threshold. best guess of threshold time

  //the error calculation is derived in Jingke's thesis Section 4.3.2.

  const int n_bins = pre_trigger_nsamps + post_trigger_nsamps + 1;
  //give a unique id to the graph, including both channel number and the graph number
  int graph_id = chdata->channel_id*1000 + 
    _current_event->GetEventData()->event_time/3.6e12/aver_time_hr;
  std::map<int,TGraphErrors*>::iterator prev = _plots.find(graph_id);
  if(prev == _plots.end()){
    //make a new TGraph and give it a name
    char name[25];
    if(chdata->channel_id>=0) sprintf(name,"average_channel%02d_%03d",graph_id/1000, graph_id%1000);
    else if(chdata->channel_id==ChannelData::CH_SUM) sprintf(name,"average_sum_%03d",graph_id%1000);
    else return 0;

    double x_array [n_bins]; 
    double y_array [n_bins]; 
    double ey_array[n_bins]; 
    for(int i=0; i < n_bins; i++){
      //set trigger index (1st bin above threshold) to be the first bin after 0
      //i=pre_trigger_nsamps when bin = trigger_index, x=0.5/chdata->sample_rate
      x_array[i] = double(i-pre_trigger_nsamps+0.5)/chdata->sample_rate;
      int bin = trigger_index - pre_trigger_nsamps + i;
      if(bin>=0 && bin<chdata->nsamps){
	double value = wave[bin]/(-65536.); //scale down pulse accordingly and convert to positive
	y_array[i]  = value;
	ey_array[i] = value*value/pulse_integral;
      }
      else{
	y_array[i] = 0;
	ey_array[i] = 0;
      }
      //      std::cout<<i<<'\t'<<x_array[i]<<'\t'<<y_array[i]<<'\t'<<ey_array[i]<<std::endl;
    }
    TGraphErrors* avg = new TGraphErrors(n_bins, x_array, y_array, NULL, ey_array);
    avg->SetName(name);
    avg->SetTitle(name);
    _plots.insert(std::make_pair(graph_id,avg));
    _sum[graph_id]  = pulse_integral;
    _sum2[graph_id] = pulse_integral*pulse_integral;
    std::cout<<"New Graph Added: "<<name<<std::endl;
  }
  else{
    TGraphErrors* avg = prev->second;
    double* y_array  = avg->GetY();
    double* ey_array = avg->GetEY();
    //add this waveform
    for(int i=0; i < n_bins; i++){
      int bin = trigger_index - pre_trigger_nsamps + i;
      //      std::cout<<i<<'\t'<<y_array[i]<<'\t'<<ey_array[i]<<std::endl;
      if(bin>=0 && bin<chdata->nsamps){
	double value = wave[bin]/(-65536.);
	y_array[i]  += value;
	ey_array[i] += value*value/pulse_integral;
      }
      //      std::cout<<chdata->SampleToTime(bin)<<'\t'<<y_array[i]<<'\t'<<ey_array[i]<<std::endl;
    }//end for
    _sum[graph_id]  += pulse_integral;
    _sum2[graph_id] += pulse_integral*pulse_integral;
  }
  
  return 0;
  
}

