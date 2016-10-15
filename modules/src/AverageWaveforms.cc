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
  RegisterParameter("normalize_by_total_waveforms", normalize_by_total_waveforms = 0,
		    "Normalize to the total number of waveforms rather than normalize to 1");
  RegisterParameter("one_pulse_only", one_pulse_only = 0,
		    "Only average events with 1 pulse");
  post_trigger_nsamps = -1;
  RegisterParameter("bipo_mode",  bipo_mode = 0, "Average BiPo events");
  RegisterParameter("run_id",  run_id=-1, "Run ID for run-specific cuts");
  RegisterParameter("spe",  spe=-1., "Spe amplitude for amplitude cut");
  RegisterParameter("spe_int",  spe_int=-1., "Spe integral for energy cut");
  RegisterParameter("include_errors_from_afterpulse",  include_errors_from_afterpulse=0, "Include the errors (and not the values) from events that fail the afterpulse cut");

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
    _total_waveforms.clear();

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
      std::cout<<"Total waveforms: "<<_total_waveforms[mapit->first]<<std::endl;
      for(int i=0; i < graph->GetN(); i++){
        y_array[i] /= a_sum;
        ey_array[i] = sqrt((ey_array[i]/a_sum-y_array[i]*y_array[i])*a_sum2)/a_sum;

	if (normalize_by_total_waveforms) {
	  y_array[i] = y_array[i]*_total_waveforms[mapit->first];
	  ey_array[i] = ey_array[i]*_total_waveforms[mapit->first];
	}
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
  if (one_pulse_only && chdata->npulses!=1) return 0;

  Pulse pulse;
  int contains_afterpulse = 0;  
  int contains_otherpulse = 0;
  
  if (run_id>0 && (spe<0 || spe_int<0)) {
    //    std::cout<<"SPE not set."<<std::endl;
    return 0;
  }
  if (bipo_mode) {   // BiPos (average second pulse)
    // CUT: Saturated
    if (chdata->saturated) return 0;
    double big_pulse_threshhold = 10*spe;
    double pre_trigger_time = 200;
    double post_trigger_time = 2000;
    std::vector<int> big_pulse_indices;
    for (int i=0; i<chdata->npulses; i++) {
      if (fabs(chdata->pulses.at(i).peak_amplitude)>big_pulse_threshhold)
	big_pulse_indices.push_back(i);
    }
    // CUT: Too many pulses
    if (big_pulse_indices.size()!=2) return 0;
    int i1 = big_pulse_indices.front();
    int i2 = big_pulse_indices.back();
    if (!chdata->pulses.at(i1).evaluated || !chdata->pulses.at(i2).evaluated) return 0;
    std::cout<<"Event "<<_current_event->GetEventData()->event_id<<". Two big pulses found at: "<<chdata->pulses.at(i1).start_time<<" "<<chdata->pulses.at(i2).start_time<<std::endl;

    // CUT: Make sure sure DAQ was triggered on the first big pulse
    if ( (chdata->pulses.at(i1).start_time < -0.1) ||
	 (chdata->pulses.at(i1).start_time > 0.1) ) {
      std::cout<<"Failed Start Time"<<std::endl;
      return 0;
    }
    // CUT: Beta Fprompt
    double fp = chdata->regions.at(0).integral/chdata->regions.at(1).integral;
    std::cout<<"1st pulse f90: "<<fp<<std::endl;
    if (fp<0.2 || fp>0.4) {
      std::cout<<"Failed F90"<<std::endl;
      return 0;
    }
    
    // CUT: Decay Time
    double decay_time = chdata->pulses.at(i2).start_time-chdata->pulses.at(i1).start_time;
    if (decay_time<pre_trigger_time ||// Stay away from beta tail
	decay_time>3000-post_trigger_time-10) {// Leave enough room for 2000us average waveform 
      std::cout<<"Failed Delay"<<std::endl;
      return 0; 
    }


    // CUT: Get rid of events with smaller pulses (afterpulses or random coincidences) that occur near the alpha event
    for (int i=0; i<chdata->npulses; i++) {
      if (i==i1 || i==i2) continue;
      if ( (chdata->pulses.at(i).start_time > chdata->pulses.at(i2).start_time-pre_trigger_time) && 
	   (chdata->pulses.at(i).start_time < chdata->pulses.at(i2).start_time+post_trigger_time) ) {
	std::cout<<"Failed Afterpulse at "<<chdata->pulses.at(i).start_time<<std::endl;
	//	return 0;
      }
    }
    
    post_trigger_nsamps = post_trigger_time*chdata->sample_rate; // Average over 2000us
    pulse = chdata->pulses.at(i2);

  } else { // Regular average
    std::cout<<"----------------------------"<<std::endl;
    std::cout<<"This event: "<<_current_event->GetEventData()->event_id<<std::endl;
    // CUT: Data quality
    std::cout<<"Checking Baseline"<<std::endl;
    if(!chdata->baseline.found_baseline) return 0;
    std::cout<<"Checking Saturated"<<std::endl;
    if(chdata->saturated) return 0;
    std::cout<<"Checking Baseline Sigma "<<chdata->baseline.sigma<<std::endl;
    if (chdata->baseline.sigma>2.0) return 0;
    std::cout<<"Checking Number of Regions"<<std::endl;
    if (chdata->regions.size()-ref_roi_index<=0) return 0;
    pulse = chdata->regions.at(ref_roi_index);
    std::cout<<"Checking Pulse Evaluated"<<std::endl;  
    if(!pulse.evaluated) return 0;
    std::cout<<"Checking Integral Sign"<<std::endl;
    if (pulse.integral>=0) return 0; //integral is negative


    //specify the end of the average, hopefully only do it once
    if(post_trigger_nsamps<0){
      post_trigger_nsamps = (pulse.end_time-pulse.half_max_time)*chdata->sample_rate;
      if(post_trigger_nsamps<=1){
	post_trigger_nsamps = -1; 
	return 0;
      }
    }

    // CUT: Amplitude and Energy
    std::cout<<"Checking Amplitude "<<pulse.peak_amplitude<<std::endl;
    if(pulse.peak_amplitude>max_pulse_height || 
       pulse.peak_amplitude<min_pulse_height ||
       pulse.npe>max_npe || pulse.npe<min_npe)
      return 0;
    // CUT: Pulse shape
    //PSD may not be necessary, so only apply the fprompt cut
    //when the cut values are specified by the users to be >-1 
    //    std::cout<<"Checking Pulse Shape"<<std::endl;
    if(max_fprompt>-100 && min_fprompt>-100 &&
       (pulse.fparameters.size()-ref_fp_index<=0 ||
	pulse.fparameters.at(ref_fp_index)>max_fprompt ||
	pulse.fparameters.at(ref_fp_index)<min_fprompt))
      return 0;
    
    // Run specific region cuts
    if (run_id>0) {
      std::vector<Pulse> reg = chdata->regions;
      double i0_09  = -reg.at(0).integral; 
      double i7     = -reg.at(1).integral;
      double i2000  = -reg.at(3).integral;

      // Peak amplitude in each region 
      const int nReg = reg.size();
      double k[nReg];
      for (int i=0; i<nReg; i++) {
	k[i] = reg.at(i).peak_amplitude;
      }
      //    std::cout<<"a"<<std::endl;
      double F90_2000 = i0_09/i2000;
      double F7_2000 = i7/i2000;
      double F90_7 = i0_09/i7;                                                                                                                                           

      //    std::cout<<"b"<<std::endl;
      // CUT: Pulse shape
      std::cout<<"Checking Pulse Shape"<<std::endl;
      if (run_id==247 && (F7_2000 > 0.8 || F7_2000 < 0.2 || F90_2000 < 0.05 || F90_2000 > 0.35)) return 0;
      if (run_id==249 && F90_2000 > 0.5) return 0;
      if (run_id==253 && F90_7 > 0.95) return 0;
      if (run_id==254 && F90_7 > 0.95) return 0;

      //    std::cout<<"c"<<std::endl;
      // CUT: Energy
      std::cout<<"Checking Energy"<<std::endl;
      if (i7/spe_int<25) return 0;

      std::cout<<"Checking Pulses"<<std::endl;
      // CUT: Numper of pulses
      if (chdata->npulses != 1) return 0;

      std::cout<<"Checking Afterpulse"<<std::endl;
      // CUT: Afterpulses
      // Region StartTime EndTime
      // 10 -2000.000 - -0.100
      // 11 -0.100 - 0.100
      // 12 0.100 - 0.196
      // 13 0.196 - 0.396
      // 14 0.396 - 0.792
      // 15 0.792 - 1.584
      // 16 1.584 - 3.160
      // 17 3.160 - 6.308
      // 18 6.308 - 12.588
      // 19 12.588 - 25.116
      // 20 25.116 - 50.116
      // 21 50.116 - 100.000
      // 22 100.000 - 199.524
      // 23 199.524 - 398.104
      // 24 398.104 - 794.328
      // 25 794.328 - 1584.888
      // 26 1584.888 - 3162.280
      // 27 3162.280 - 6309.568
      // 28 6309.568 - 12589.300
      // 29 12589.300 - 25118.900
      // 30 25118.900 - 29999.996
      // 31 29999.996 - 29999.996

      if (run_id==247) {
	if (k[10] > 42.267)  contains_otherpulse = 10;
	if (k[12] > 55.0885) contains_otherpulse = 12;
	if (k[13] > 50.5217) contains_otherpulse = 13;
	if (k[14] > 47.2261) contains_afterpulse = 14;
	if (k[15] > 45.346)  contains_afterpulse = 15;
	if (k[16] > 43.8861) contains_afterpulse = 16;
	if (k[17] > 42.306)  contains_otherpulse = 17;
	if (k[18] > 41.0423) contains_otherpulse = 18;
	if (k[19] > 38.8841) contains_otherpulse = 19;
	if (k[20] > 39.3019) contains_otherpulse = 20;
	if (k[21] > 38.7607) contains_otherpulse = 21;
	if (k[22] > 38.3034) contains_otherpulse = 22;
	if (k[23] > 38.2968) contains_otherpulse = 23;
	if (k[24] > 38.2057) contains_otherpulse = 24;
	if (k[25] > 38.7132) contains_otherpulse = 25;
	if (k[26] > 40.1554) contains_otherpulse = 26;

      }
      if (run_id==249) {
	if (k[10] > 52.9914) contains_otherpulse = 10;
	if (k[12] > 105.488) contains_otherpulse = 12;
	if (k[13] > 85.2509) contains_otherpulse = 13;
	if (k[14] > 71.2988) contains_afterpulse = 14;
	if (k[15] > 58.3709) contains_afterpulse = 15;
	if (k[16] > 49.0045) contains_afterpulse = 16;
	if (k[17] > 42.7381) contains_otherpulse = 17;
	if (k[18] > 39.8743) contains_otherpulse = 18;
	if (k[19] > 38.5458) contains_otherpulse = 19;
	if (k[20] > 37.9472) contains_otherpulse = 20;
	if (k[21] > 37.6178) contains_otherpulse = 21;
	if (k[22] > 37.821)  contains_otherpulse = 22;
	if (k[23] > 37.7643) contains_otherpulse = 23;
	if (k[24] > 39.7261) contains_otherpulse = 24;
	if (k[25] > 44.7066) contains_otherpulse = 25;
	if (k[26] > 42.915)  contains_otherpulse = 26;
      }
      if (run_id==253) {
	if (k[10] > 39.3316) contains_otherpulse = 10;
	if (k[12] > 56.6115) contains_otherpulse = 12;
	if (k[13] > 46.9317) contains_otherpulse = 13;
	if (k[14] > 41.1792) contains_afterpulse = 14;
	if (k[15] > 38.7479) contains_afterpulse = 15;
	if (k[16] > 37.8252) contains_afterpulse = 16;
	if (k[17] > 36.4579) contains_otherpulse = 17;
	if (k[18] > 36.1478) contains_otherpulse = 18;
	if (k[19] > 36.3116) contains_otherpulse = 19;
	if (k[20] > 36.5952) contains_otherpulse = 20;
	if (k[21] > 35.32)   contains_otherpulse = 21;
	if (k[22] > 35.7159) contains_otherpulse = 22;
	if (k[23] > 37.8858) contains_otherpulse = 23;
	if (k[24] > 39.0954) contains_otherpulse = 24;
	if (k[25] > 40.1744) contains_otherpulse = 25;
	if (k[26] > 40.5322) contains_otherpulse = 26;
      }
      if (run_id==254) {
	if (k[10] > 42.7789) contains_otherpulse = 10;
	if (k[12] > 51.4615) contains_otherpulse = 12;
	if (k[13] > 46.4879) contains_otherpulse = 13;
	if (k[14] > 43.9157) contains_afterpulse = 14;
	if (k[15] > 43.2425) contains_afterpulse = 15;
	if (k[16] > 42.2806) contains_afterpulse = 16;
	if (k[17] > 41.3621) contains_otherpulse = 17;
	if (k[18] > 40.9589) contains_otherpulse = 18;
	if (k[19] > 40.4873) contains_otherpulse = 19;
	if (k[20] > 40.6747) contains_otherpulse = 20;
	if (k[21] > 40.2588) contains_otherpulse = 21;
	if (k[22] > 40.2063) contains_otherpulse = 22;
	if (k[23] > 40.2976) contains_otherpulse = 23;
	if (k[24] > 41.4719) contains_otherpulse = 24;
	if (k[25] > 41.1702) contains_otherpulse = 25;
	if (k[26] > 41.5278) contains_otherpulse = 26;
      }    

      if (run_id==276) {
	double lim = 3*spe;
	if (k[10]>lim) contains_otherpulse = 10;
	if (k[12]>lim) contains_otherpulse = 12;
	if (k[13]>lim) contains_otherpulse = 13;
	if (k[14]>lim) contains_afterpulse = 14;
	if (k[15]>lim) contains_afterpulse = 15;
	if (k[16]>lim) contains_afterpulse = 16;
	if (k[17]>lim) contains_otherpulse = 17;
	if (k[18]>lim) contains_otherpulse = 18;
	if (k[19]>lim) contains_otherpulse = 19;
	if (k[20]>lim) contains_otherpulse = 20;
	if (k[21]>lim) contains_otherpulse = 21;
	if (k[22]>lim) contains_otherpulse = 22;
	if (k[23]>lim) contains_otherpulse = 23;
	if (k[24]>lim) contains_otherpulse = 24;
	if (k[25]>lim) contains_otherpulse = 25;
	if (k[26]>lim) contains_otherpulse = 26;
	if (k[27]>lim) contains_otherpulse = 27;
	if (k[28]>lim) contains_otherpulse = 28;
	if (k[29]>lim) contains_otherpulse = 29;
	if (k[30]>lim) contains_otherpulse = 30;
      }
      std::cout<<"Other Pulse: "<<contains_otherpulse<<" After Pulse: "<<contains_afterpulse<<std::endl;
      //    if (!run_id==247) {
      if (contains_otherpulse) return 0;
      if (contains_afterpulse && !include_errors_from_afterpulse) return 0;
      // }
      // if (run_id==247 && contains_afterpulse) std::cout<<"Should have cut this event!!!!!!!!!!"<<std::endl;
      std::cout<<"Passed Afterpulse"<<std::endl;

    }// End run-specific cuts

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
    //    if(pulse.min<1 //|| pulse.peak_amplitude>-0.08*pulse.integral
       //|| pulse.peak_amplitude<chdata->baseline.mean-pulse.min-0.5
      // ) return 0;
    //    std::cout<<"k"<<std::endl;
  }
  
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
    if (contains_afterpulse) return 0; // Don't start a new graph with an afterpulse event
    //make a new TGraph and give it a name
    char name[25];
    if(chdata->channel_id>=0) sprintf(name,"average_channel%02d_%03d",graph_id/1000, graph_id%1000);
    else if(chdata->channel_id==ChannelData::CH_SUM) sprintf(name,"average_sum_%03d",graph_id%1000);
    else return 0;
    _total_waveforms[graph_id] = 1; // First waveform in new TGraph

    double *x_array = new double[n_bins]; 
    double *y_array = new double[n_bins]; 
    double *ey_array= new double[n_bins]; 

    for(int i=0; i < n_bins; i++){
      //set trigger index (1st bin above threshold) to be the first bin after 0
      //i=pre_trigger_nsamps when bin = trigger_index, x=0.5/chdata->sample_rate
      x_array[i] = double(i-pre_trigger_nsamps+0.5)/chdata->sample_rate;
      int bin = trigger_index - pre_trigger_nsamps + i;
      if(bin>=0 && bin<chdata->nsamps){
	double value = wave[bin]/(-65536.); //scale down pulse accordingly and convert to positive
	y_array[i]  = value;
	ey_array[i] = value*value/pulse_integral;
	if (contains_afterpulse && include_errors_from_afterpulse) y_array[i] = 0; // only add the error for the event
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
    double* x_array  = avg->GetX();
    double* y_array  = avg->GetY();
    double* ey_array = avg->GetEY();
    //add this waveform
    for(int i=0; i < n_bins; i++){
      int bin = trigger_index - pre_trigger_nsamps + i;
      double t = x_array[i];
      //      std::cout<<i<<'\t'<<y_array[i]<<'\t'<<ey_array[i]<<std::endl;
      if(bin>=0 && bin<chdata->nsamps){
	double value = wave[bin]/(-65536.);
	y_array[i]  += value;
	ey_array[i] += value*value/pulse_integral;
	if (include_errors_from_afterpulse) {
	  if (contains_afterpulse) {
	    y_array[i] -= value; // subtact the average, keeping the error
	    //if (i%100==0 && t<10) std::cout<<"Subtracting Value from t "<<t<<std::endl;
	  }
	  if ((contains_afterpulse==14 && (t<0.396 || t>0.792)) || (contains_afterpulse==15 && (t<0.792 || t>1.584)) || (contains_afterpulse==16 && (t<1.584 || t>3.160))) {
	    ey_array[i] -= value*value/pulse_integral; // only add the error for the event
	    //if (i%100==0 && t<10) std::cout<<"Subracting error from t "<<t<<std::endl;

	  }
	}
      }
      //      std::cout<<chdata->SampleToTime(bin)<<'\t'<<y_array[i]<<'\t'<<ey_array[i]<<std::endl;
    }//end for
    if (!contains_afterpulse) {
      _sum[graph_id]  += pulse_integral;
      _sum2[graph_id] += pulse_integral*pulse_integral;
      _total_waveforms[graph_id]++;
    }
  }
  
  return 0;
  
}

