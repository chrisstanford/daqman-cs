#include "AverageWaveforms.hh"
#include "TGraphFunctions.hh"
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
  RegisterParameter("post_trigger_time", post_trigger_time = 10,
		    "time after pulse rise to start summing waveforms");
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
  RegisterParameter("bipo_mode",  bipo_mode = 0, "Average BiPo events");
  RegisterParameter("bipo_pulse_id",  bipo_pulse_id = 1, "BiPo Pulse to Average (0=beta,1=alpha)");
  RegisterParameter("run_id",  run_id=-1, "Run ID for run-specific cuts");
  RegisterParameter("spe_height_mean",  spe_height_mean=-1., "Spe amplitude for amplitude cut");
  RegisterParameter("spe_int_mean",  spe_int_mean=-1., "Spe integral for energy cut");
  RegisterParameter("pre_spe_nsamps",  pre_spe_nsamps=-1., "Number of samples before spe peak");
  RegisterParameter("post_spe_nsamps",  post_spe_nsamps=-1., "Number of samples avter spe peak");
  RegisterParameter("verbose",  verbose=0., "Print out statements");
  RegisterParameter("afterpulse_cut",  afterpulse_cut = 1, "Perform Afterpulse Cut");
  RegisterParameter("afterpulse_height_spe",  afterpulse_height_spe = 2.5, "Afterpulse height threshold in spe");

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
  waveforms = new TTree("waveforms","Waveforms to be averaged");
  waveforms->Branch("waveform",&waveform);
  spewfs = new TTree("spewfs","SPE waveforms");
  spewfs->Branch("spewf",&spewf);
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
      int nWaveforms = _total_waveforms[mapit->first];
      if(a_sum<=0)  a_sum=1;
      if(a_sum2<=0) a_sum2=1;
      std::cout<<"Total waveforms: "<<nWaveforms<<std::endl;
      for(int i=0; i < graph->GetN(); i++){
	// For RMS
	//	ey_array[i] = sqrt(1/nWaveforms * ey_array[i]);
        // For Variance
	ey_array[i] = sqrt( ey_array[i]/a_sum - y_array[i]*y_array[i]/(a_sum*a_sum) );
	ey_array[i] /= sqrt(nWaveforms);
        y_array[i] /= a_sum;
	//        ey_array[i] = sqrt((ey_array[i]/a_sum-y_array[i]*y_array[i])*a_sum2)/a_sum;
	if (normalize_by_total_waveforms) {
	  y_array[i] = y_array[i]*nWaveforms;
	  ey_array[i] = ey_array[i]*nWaveforms;
	}
      }
      (mapit->second)->SetFillStyle(3002);
      (mapit->second)->SetFillColor(kRed);
      (mapit->second)->Write();
    }
    waveforms->Write();
    spewfs->Write();
  }
  Cleanup();
  return 0;
}
double IntegrateWindow(std::vector<Spe> spes, double start_time, double length, double spe_int_mean) {
  double sum(0.);
  for (size_t s=0; s<spes.size(); s++) {
    Spe spe = spes.at(s);
    double t = spe.start_time;
    if (t<start_time || t>start_time+length) continue;
    double integral = spe.integral;
    sum += integral/spe_int_mean;
  }
  return sum;
}
bool AverageWaveforms::PassesCuts(const int run_id, ChannelData* channel, const double spe_int_mean, const double spe_height, TString mode = "normal", TString option = "") {
  if (!channel) return 0;
  // CUT: Found baseline
  if (!channel->baseline.found_baseline) return 0;
  // CUT: Saturated
  if (channel->saturated) return 0;

  if (mode=="normal") {
    if (channel->npulses!=1) return 0;
    Pulse pulse = channel->pulses.at(0);
    if (run_id == 270) {
      if (fabs(pulse.peak_amplitude) < 150) return 0;
      if (fabs(pulse.peak_amplitude) > 1000) return 0;
    }
    if (run_id >= 310 && run_id<320) { // UV 128nm runs
      if (fabs(pulse.peak_amplitude) < 1500) return 0;
      if (fabs(pulse.peak_amplitude) > 3000) return 0;
    }
    std::vector<Pulse> reg = channel->regions;
    double i0_09  = -reg.at(0).integral; 
    double i7     = -reg.at(1).integral;
    double i2000  = -reg.at(5).integral;
    double F90_2000 = i0_09/i2000;
    double F7_2000 = i7/i2000;
    double F90_7 = i0_09/i7;                                                                                                                                           

    if (run_id>=300 && run_id<310) { // Vacuum Po210 runs at room and LAr temp
      if (i7<1000 || i7>5000 || F90_7>0.8) return 0;
    }
    // CUT: Pulse shape
    if (run_id==247 && (F7_2000 > 0.8 || F7_2000 < 0.2 || F90_2000 < 0.05 || F90_2000 > 0.35)) return 0;
    if (run_id==249 && F90_2000 > 0.5) return 0;
    if (run_id==253 && F90_7 > 0.95) return 0;
    if (run_id==254 && F90_7 > 0.95) return 0;
    // CUT: Energy
    if (run_id==249 && -pulse.integral/spe_int_mean<50) return 0;
    //    if ((run_id==273||run_id==275||run_id==276) && -pulse.integral/spe_int_mean>150) return 0;

    if (afterpulse_cut) {
      std::vector<Spe> spes = channel->single_pe;
      for (size_t s=0; s<spes.size(); s++) {
	Spe spe = spes.at(s);
	double t = spe.peak_time;
	double height = spe.amplitude;
	// Conservative
	if (height>afterpulse_height_spe*spe_height && (t>0.2||t<-0.1)) return 0; 
      }
    }
  } else if (mode == "bipo") {
    std::vector<Pulse> pulses = channel->pulses;
    if (channel->npulses!=2) return 0;
    Pulse P0 = pulses.at(0);
    Pulse P1 = pulses.at(1);
    //    cout<<P0.start_time<<" "<<P0.peak_amplitude<<endl;
    //    cout<<P1.start_time<<" "<<P1.peak_amplitude<<endl;
    //  fparameter_times [ 0.09 , 7 , 15 , 50 , 100 , 2000 , 3000 ]
    //    double P0_F90_7 = P0.fparameters.at(0)/P0.fparameters.at(1);
    double P0_F7_P = P0.fparameters.at(1);
    double P1_F90_P = P1.fparameters.at(0);
    double P1_F7_P = P1.fparameters.at(1);
    double P1_F15_P = P1.fparameters.at(2);
    //    double P1_F50_2000 = P1.fparameters.at(1);
    //    double P1_50_2000 = P1.fparameters.at(3);
    double P1_F90_7 = P1_F90_P/P1_F7_P;
    //    double P1_F90_50 = P1_F90_2000/P1_F50_2000;

    //    if (P1_F7_2000>0.8) return 0;

    // if (P0_F90_7<0.2) return 0;
    // if (P0_F90_7>0.4) return 0;
    if (run_id == 340 && P1.peak_amplitude>600) return 0;
    if (run_id == 340 && P1_F15_P>0.75) return 0;

    if (P0_F7_P<0.8) return 0;
    // if (P0.start_time<-0.5) return 0;
    // if (P0.start_time>0.5) return 0;
    // if (-P1.integral/spe_int_mean<25) return 0;

    if (P1_F90_7>0.9) return 0;

    double delay = P1.start_time-P0.start_time;
    if (delay<25 || delay>999) return 0;

    std::vector<Spe> spes = channel->single_pe;
    double alpha_start_time = P1.start_time;
    for (size_t s=0; s<spes.size(); s++) {
      Spe spe = spes.at(s);
      double t = spe.peak_time;
      double height = spe.amplitude;
      // Pre-pulse spe
      if (IntegrateWindow(spes, 20, alpha_start_time-20-0.1, spe_int_mean)>10) return 0;
      // Cherenkov Pulses / Afterpulses
      if (afterpulse_cut) {
	if (height>afterpulse_height_spe*spe_height && (t<-1 || (t>5&&t<alpha_start_time-1) || t>alpha_start_time+0.3)) return 0;
      }
    }
  } else if (mode == "bipoP0alpha") {
    std::vector<Pulse> pulses = channel->pulses;
    if (channel->npulses!=1) return 0;
    Pulse P0 = pulses.at(0);
    //  regions [ -0.05 : 0.09 , -0.05 : 7 , -0.05 : 15 , -0.05 : 40 , -0.05 : 100 , -50 : 2000 , -100 : 3000 , -2000 : 30000 , -2000 : -0.05 , -0.05 : 1 ]
    //  fparameter_times [ 0.09 , 7 , 15 , 50 , 100 , 2000 , 3000 ]
    double P0_F90_3000 = channel->regions.at(0).integral/channel->regions.at(6).integral;
    double P0_F7_3000 = channel->regions.at(1).integral/channel->regions.at(6).integral;
    if (P0_F90_3000<0.1) return 0;
    if (P0_F90_3000>0.15) return 0;
    if (P0_F7_3000<0.77) return 0;
    if (P0_F7_3000>0.81) return 0;
    // if (P0.start_time<-0.5) return 0;
    // if (P0.start_time>0.5) return 0;
    // if (-P0.integral/spe_int_mean<25) return 0;

    std::vector<Spe> spes = channel->single_pe;
    //    double alpha_start_time = P0.start_time;

    // for (size_t s=0; s<spes.size(); s++) {
    //   Spe spe = spes.at(s);
    //   double t = spe.peak_time;
    //   double height = spe.amplitude;
    //      if (height>3*spe_height && t>alpha_start_time+1) return 0;
    //    }  
  // } else if (mode == "laser") {
  //   if (channel->npulses!=1) return 0;
  }
  //  if (P1_F90_2000>0.7) cout<<"outlier "<<-P1.integral/spe_int_mean<<" "<<P1_F90_2000<<" "<<option<<endl;

  return 1;
} 


//this is not a critical module, return 0 for all cases
int AverageWaveforms::Process(ChannelData *chdata){
  Pulse pulse;
  //bool PassesCuts(const int run_id, ChannelData* channel, const double spe_int_mean, const double spe_height, TString mode = "normal"
  if (run_id>0 && spe_height_mean>0 && spe_int_mean>0) { 
    if (bipo_mode) {   // BiPos (average second pulse)
      if (!PassesCuts(run_id,chdata,spe_int_mean,spe_height_mean,"bipo")) return 0;
      pulse = chdata->pulses.at(bipo_pulse_id);
    } else {
      if (!PassesCuts(run_id,chdata,spe_int_mean,spe_height_mean)) return 0;
      pulse = chdata->pulses.at(0);    
    }
  } else { // Regular average
    if (verbose) std::cout<<"----------------------------"<<std::endl;
    if (verbose) std::cout<<"This event: "<<_current_event->GetEventData()->event_id<<std::endl;
    // CUT: Data quality
    if (verbose) std::cout<<"Checking Baseline"<<std::endl;
    if(!chdata->baseline.found_baseline) return 0;
    if (verbose) std::cout<<"Checking Saturated"<<std::endl;
    //    if(chdata->saturated) return 0;
    if (verbose) std::cout<<"Checking Baseline Sigma "<<chdata->baseline.sigma<<std::endl;
    if (chdata->baseline.sigma>2.0) return 0;
    if (verbose) std::cout<<"Checking Pulses"<<std::endl;
    //    if (chdata->npulses!=1) return 0;
    if (chdata->npulses<1) return 0;
    pulse = chdata->pulses.at(0);
    if (verbose) std::cout<<"Checking Pulse Evaluated"<<std::endl;  
    if(!pulse.evaluated) return 0;
    if (verbose) std::cout<<"Checking Pulse Start Time "<<pulse.start_time<<std::endl;
    if (pulse.start_time<-0.1||pulse.start_time>0.2) return 0;
    if (verbose) std::cout<<"Checking Integral Sign"<<std::endl;
    if (pulse.integral>=0) return 0; //integral is negative
    // CUT: Amplitude and Energy
    if (verbose) std::cout<<"Checking Amplitude "<<pulse.peak_amplitude<<" and npe "<<pulse.npe<<std::endl;
    if(pulse.peak_amplitude>max_pulse_height || 
       pulse.peak_amplitude<min_pulse_height ||
       pulse.npe>max_npe || pulse.npe<min_npe)
      return 0;
    // CUT: Pulse shape
    // PSD may not be necessary, so only apply the fprompt cut
    // when the cut values are specified by the users to be >-1 
    //    if (verbose) std::cout<<"Checking Pulse Shape"<<std::endl;
    double fp = chdata->regions[ref_fp_index].integral/chdata->regions[ref_roi_index].integral;
    if(max_fprompt>-100 && min_fprompt>-100 &&
       (fp>max_fprompt || fp<min_fprompt))
       // (pulse.fparameters.size()-ref_fp_index<=0 ||
       // 	pulse.fparameters.at(ref_fp_index)>max_fprompt ||
       // 	pulse.fparameters.at(ref_fp_index)<min_fprompt))
      return 0;
    // CUT: Afterpulses
    if (afterpulse_cut) {
      if (verbose) std::cout<<"Checking Afterpulses"<<std::endl;
      if (chdata->pulses[0].has_cherenkov) return 0;
    }
  }

  std::cout<<"average event: "<<_current_event->GetEventData()->event_id<<std::endl;
  
  
  //  int trigger_index=pulse.start_index;
  int trigger_index=pulse.peak_index;
  //  int trigger_index=chdata->TimeToSample(pulse.half_max_time+0.5/chdata->sample_rate);
  int pre_trigger_nsamps  = pre_trigger_time *chdata->sample_rate;
  int post_trigger_nsamps  = post_trigger_time *chdata->sample_rate;
  //  std::cout<<pre_trigger_nsamps<<" "<<post_trigger_nsamps<<std::endl;
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

  // Save full waveform
  double *wf_x = new double[n_bins]; 
  double *wf_y = new double[n_bins]; 
  double *wf_ex = new double[n_bins]; 
  double *wf_ey = new double[n_bins]; 
  for(int i=0; i<n_bins; i++){
    //set trigger index (1st bin above threshold) to be the first bin after 0
    //i=pre_trigger_nsamps when bin = trigger_index, x=0.5/chdata->sample_rate
    wf_x[i] = double(i-pre_trigger_nsamps+0.5)/chdata->sample_rate;
    wf_y[i] = 0;
    wf_ex[i] = 0;
    wf_ey[i] = 0;
    int bin = trigger_index - pre_trigger_nsamps + i;
    if(bin>=0 && bin<chdata->nsamps)
      wf_y[i] = -wave[bin];
    else
      wf_y[i] = 0;
  }
  TGraphErrors* wf_orig = new TGraphErrors(n_bins,wf_x,wf_y,wf_ex,wf_ey);
  //  TGraphErrors* wf_post = postTriggerGraph(wf_orig,0.05);
  TGraphErrors* wf_post = postTriggerGraph(wf_orig,0.025);
  waveform = logCombineGraph(wf_post,1000);
  
  waveform->SetName(Form("ch%d_ev%d",chdata->channel_id,_current_event->GetEventData()->event_id));
  // Fill waveform tree
  waveforms->Fill();
  delete[] wf_x;
  delete[] wf_y;
  delete[] wf_ex;
  delete[] wf_ey;  
  delete wf_orig;
  delete wf_post;
  delete waveform;
  
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
    _total_waveforms[graph_id] = 1; // First waveform in new TGraph

    double *x_array = new double[n_bins]; 
    double *y_array = new double[n_bins]; 
    double *ey_array= new double[n_bins]; 

    for(int i=0; i < n_bins; i++){
      //set trigger index (1st bin above threshold) to be the first bin after 0
      //i=pre_trigger_nsamps when bin = trigger_index, x=0.5/chdata->sample_rate
      x_array[i] = double(i-pre_trigger_nsamps)/chdata->sample_rate;
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
      if(bin>=0 && bin<chdata->nsamps){
	double value = wave[bin]/(-65536.);
	y_array[i]  += value;
	ey_array[i] += value*value/pulse_integral;
	// For RMS
	//ey_array[i] += value*value/(pulse_integral*pulse_integral);
      }
      //      std::cout<<chdata->SampleToTime(bin)<<'\t'<<y_array[i]<<'\t'<<ey_array[i]<<std::endl;
    }//end for
    _sum[graph_id]  += pulse_integral;
    _sum2[graph_id] += pulse_integral*pulse_integral;
    _total_waveforms[graph_id]++;
  }
  if (pre_spe_nsamps>0 && post_spe_nsamps>0) {
    // SPES
    const int n_spebins = pre_spe_nsamps + post_spe_nsamps + 1;

    // loop over spes
    std::vector<Spe> spes = chdata->single_pe;
    for (size_t s=0; s<spes.size(); s++) {
      Spe spe = spes.at(s);
      double peak_index = spe.peak_index;
      double spe_height = spe.amplitude;
      int prev_end_index = s>0 ? spes.at(s-1).peak_index+post_spe_nsamps : 0;
      int next_start_index = s<spes.size()-1 ? spes.at(s+1).peak_index-pre_spe_nsamps : chdata->nsamps;
      //if (height>3*spe_height && (t>0.3||t<-0.1)) return 0;
      // Conservative
      if (spe_height<0.5*spe_height_mean || spe_height>1.5*spe_height_mean) continue; // Check height is consistent with spe
      if (peak_index-pre_spe_nsamps<=prev_end_index || peak_index+post_spe_nsamps>=next_start_index) continue; // Check overlap

      // Save waveform
      double *spewf_x = new double[n_spebins]; 
      double *spewf_y = new double[n_spebins]; 
      double *spewf_ex = new double[n_spebins]; 
      double *spewf_ey = new double[n_spebins]; 
    
      for(int i=0; i<n_spebins; i++){
	//set trigger index (1st bin above threshold) to be the first bin after 0
	//i=pre_trigger_nsamps when bin = trigger_index, x=0.5/chdata->sample_rate
	spewf_x[i] = double(i-pre_spe_nsamps)/chdata->sample_rate;
	spewf_y[i] = 0;
	spewf_ex[i] = 0;
	spewf_ey[i] = 0;
	int bin = peak_index - pre_spe_nsamps + i;
	if(bin>=0 && bin<chdata->nsamps)
	  spewf_y[i] = -wave[bin];
	else
	  spewf_y[i] = 0;
      }
      spewf = new TGraphErrors(n_spebins,spewf_x,spewf_y,spewf_ex,spewf_ey);
      spewf->SetName(Form("ch%d_ev%d_spe%d",chdata->channel_id,_current_event->GetEventData()->event_id,int(s)));
      // Fill spe tree
      spewfs->Fill();
      delete[] spewf_x;
      delete[] spewf_y;
      delete[] spewf_ex;
      delete[] spewf_ey;  
      delete spewf;
    }    
  }
  
  
  
  return 0;

}

/*
int AverageWaveforms::Process(ChannelData *chdata){
  bool verbose = true;
  if (one_pulse_only && chdata->npulses!=1) return 0;
  Pulse pulse;
  
  if (run_id>0 && (spe<0 || spe_int_mean<0)) {
    if (verbose) std::cout<<"SPE not set."<<std::endl;
    return 0;
  }
  if (bipo_mode) {   // BiPos (average second pulse)
    if (verbose) std::cout<<"----------------------------"<<std::endl;
    if (verbose) std::cout<<"This event: "<<_current_event->GetEventData()->event_id<<std::endl;
    if (verbose) std::cout<<"Checking baseline"<<std::endl;
    // CUT: Found baseline
    if (!chdata->baseline.found_baseline) return 0;
    if (verbose) std::cout<<"Checking saturated"<<std::endl;
    // CUT: Saturated
    if (chdata->saturated) return 0;

    // New method:
    std::vector<Pulse> pulses = chdata->pulses;
    if (verbose) std::cout<<"Checking npulses "<<pulses.size()<<std::endl;
    // CUT: npulses
    if (pulses.size()!=2) return 0;
    // CUT: beta f90
    double p0_f90_f7 = pulses.at(0).fparameters.at(0)/pulses.at(0).fparameters.at(1);
    if (verbose) std::cout<<"Checking beta f90 "<<p0_f90_f7<<std::endl;
    if (p0_f90_f7<0.15 || p0_f90_f7>0.45) return 0;
    if (verbose) std::cout<<"Checking alpha afterpulse"<<std::endl;
    // CUT: alpha afterpulse/cherenkov
    if (pulses.at(1).has_cherenkov) return 0;
    if (verbose) std::cout<<"Checking decay time"<<std::endl;
    // CUT: Decay Time
    double decay_time = chdata->pulses.at(1).start_time-chdata->pulses.at(0).start_time;
    if (decay_time<101) return 0; 
    if (verbose) std::cout<<"Success!"<<std::endl;
	
    double post_trigger_time = 2000;    
    post_trigger_nsamps = post_trigger_time*chdata->sample_rate;
    pulse = chdata->pulses.at(1);

  } else { // Regular average
    if (verbose) std::cout<<"----------------------------"<<std::endl;
    if (verbose) std::cout<<"This event: "<<_current_event->GetEventData()->event_id<<std::endl;
    // CUT: Data quality
    if (verbose) std::cout<<"Checking Baseline"<<std::endl;
    if(!chdata->baseline.found_baseline) return 0;
    if (verbose) std::cout<<"Checking Saturated"<<std::endl;
    if(chdata->saturated) return 0;
    if (verbose) std::cout<<"Checking Baseline Sigma "<<chdata->baseline.sigma<<std::endl;
    if (chdata->baseline.sigma>2.0) return 0;
    if (verbose) std::cout<<"Checking Number of Regions"<<std::endl;
    if (chdata->regions.size()-ref_roi_index<=0) return 0;
    pulse = chdata->regions.at(ref_roi_index);
    if (verbose) std::cout<<"Checking Pulse Evaluated"<<std::endl;  
    if(!pulse.evaluated) return 0;
    if (verbose) std::cout<<"Checking Integral Sign"<<std::endl;
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
    if (verbose) std::cout<<"Checking Amplitude "<<pulse.peak_amplitude<<std::endl;
    if(pulse.peak_amplitude>max_pulse_height || 
       pulse.peak_amplitude<min_pulse_height ||
       pulse.npe>max_npe || pulse.npe<min_npe)
      return 0;
    // CUT: Pulse shape
    //PSD may not be necessary, so only apply the fprompt cut
    //when the cut values are specified by the users to be >-1 
    //    if (verbose) std::cout<<"Checking Pulse Shape"<<std::endl;
    if(max_fprompt>-100 && min_fprompt>-100 &&
       (pulse.fparameters.size()-ref_fp_index<=0 ||
	pulse.fparameters.at(ref_fp_index)>max_fprompt ||
	pulse.fparameters.at(ref_fp_index)<min_fprompt))
      return 0;

    if (verbose) std::cout<<"Checking Pulses"<<std::endl;
    // CUT: Numper of pulses
    if (chdata->npulses != 1) return 0;

    if (verbose) std::cout<<"Checking Afterpulses"<<std::endl;
    // CUT: Afterpulses
    if (chdata->pulses[0].has_cherenkov) return 0;
    

    // Run specific cuts
    if (run_id>0) {
      // ( -0.05 , 0.09 )   #prompt pulse
      // ( -0.05 , 7 )
      // ( -0.05 , 15 )
      // ( -0.05 , 40 )
      // ( -0.05 , 100 )
      // ( -50 ,  2000 )
      // ( -100 , 3000 )
      // ( -2000 , 30000 )
      // ( -2000 , -0.05 )

      std::vector<Pulse> reg = chdata->regions;
      double i0_09  = -reg.at(0).integral; 
      double i7     = -reg.at(1).integral;
      double i2000  = -reg.at(5).integral;

      // // Peak amplitude in each region 
      // const int nReg = reg.size();
      // double k[nReg];
      // for (int i=0; i<nReg; i++) {
      // 	k[i] = reg.at(i).peak_amplitude;
      // }
      //    if (verbose) std::cout<<"a"<<std::endl;
      double F90_2000 = i0_09/i2000;
      double F7_2000 = i7/i2000;
      double F90_7 = i0_09/i7;                                                                                                                                           

      //    if (verbose) std::cout<<"b"<<std::endl;
      // CUT: Pulse shape
      if (verbose) std::cout<<"Checking Pulse Shape"<<std::endl;
      if (run_id==247 && (F7_2000 > 0.8 || F7_2000 < 0.2 || F90_2000 < 0.05 || F90_2000 > 0.35)) return 0;
      if (run_id==249 && F90_2000 > 0.5) return 0;
      if (run_id==253 && F90_7 > 0.95) return 0;
      if (run_id==254 && F90_7 > 0.95) return 0;

      //    if (verbose) std::cout<<"c"<<std::endl;
      // CUT: Energy
      if (verbose) std::cout<<"Checking Energy"<<std::endl;
      if (i7/spe_int_mean<25) return 0;


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

      // if (run_id==247) {
      // 	if (k[10] > 42.267)  contains_otherpulse = 10;
      // 	if (k[12] > 55.0885) contains_otherpulse = 12;
      // 	if (k[13] > 50.5217) contains_otherpulse = 13;
      // 	if (k[14] > 47.2261) contains_afterpulse = 14;
      // 	if (k[15] > 45.346)  contains_afterpulse = 15;
      // 	if (k[16] > 43.8861) contains_afterpulse = 16;
      // 	if (k[17] > 42.306)  contains_otherpulse = 17;
      // 	if (k[18] > 41.0423) contains_otherpulse = 18;
      // 	if (k[19] > 38.8841) contains_otherpulse = 19;
      // 	if (k[20] > 39.3019) contains_otherpulse = 20;
      // 	if (k[21] > 38.7607) contains_otherpulse = 21;
      // 	if (k[22] > 38.3034) contains_otherpulse = 22;
      // 	if (k[23] > 38.2968) contains_otherpulse = 23;
      // 	if (k[24] > 38.2057) contains_otherpulse = 24;
      // 	if (k[25] > 38.7132) contains_otherpulse = 25;
      // 	if (k[26] > 40.1554) contains_otherpulse = 26;

      // }
      // if (run_id==249) {
      // 	if (k[10] > 52.9914) contains_otherpulse = 10;
      // 	if (k[12] > 105.488) contains_otherpulse = 12;
      // 	if (k[13] > 85.2509) contains_otherpulse = 13;
      // 	if (k[14] > 71.2988) contains_afterpulse = 14;
      // 	if (k[15] > 58.3709) contains_afterpulse = 15;
      // 	if (k[16] > 49.0045) contains_afterpulse = 16;
      // 	if (k[17] > 42.7381) contains_otherpulse = 17;
      // 	if (k[18] > 39.8743) contains_otherpulse = 18;
      // 	if (k[19] > 38.5458) contains_otherpulse = 19;
      // 	if (k[20] > 37.9472) contains_otherpulse = 20;
      // 	if (k[21] > 37.6178) contains_otherpulse = 21;
      // 	if (k[22] > 37.821)  contains_otherpulse = 22;
      // 	if (k[23] > 37.7643) contains_otherpulse = 23;
      // 	if (k[24] > 39.7261) contains_otherpulse = 24;
      // 	if (k[25] > 44.7066) contains_otherpulse = 25;
      // 	if (k[26] > 42.915)  contains_otherpulse = 26;
      // }
      // if (run_id==253) {
      // 	if (k[10] > 39.3316) contains_otherpulse = 10;
      // 	if (k[12] > 56.6115) contains_otherpulse = 12;
      // 	if (k[13] > 46.9317) contains_otherpulse = 13;
      // 	if (k[14] > 41.1792) contains_afterpulse = 14;
      // 	if (k[15] > 38.7479) contains_afterpulse = 15;
      // 	if (k[16] > 37.8252) contains_afterpulse = 16;
      // 	if (k[17] > 36.4579) contains_otherpulse = 17;
      // 	if (k[18] > 36.1478) contains_otherpulse = 18;
      // 	if (k[19] > 36.3116) contains_otherpulse = 19;
      // 	if (k[20] > 36.5952) contains_otherpulse = 20;
      // 	if (k[21] > 35.32)   contains_otherpulse = 21;
      // 	if (k[22] > 35.7159) contains_otherpulse = 22;
      // 	if (k[23] > 37.8858) contains_otherpulse = 23;
      // 	if (k[24] > 39.0954) contains_otherpulse = 24;
      // 	if (k[25] > 40.1744) contains_otherpulse = 25;
      // 	if (k[26] > 40.5322) contains_otherpulse = 26;
      // }
      // if (run_id==254) {
      // 	if (k[10] > 42.7789) contains_otherpulse = 10;
      // 	if (k[12] > 51.4615) contains_otherpulse = 12;
      // 	if (k[13] > 46.4879) contains_otherpulse = 13;
      // 	if (k[14] > 43.9157) contains_afterpulse = 14;
      // 	if (k[15] > 43.2425) contains_afterpulse = 15;
      // 	if (k[16] > 42.2806) contains_afterpulse = 16;
      // 	if (k[17] > 41.3621) contains_otherpulse = 17;
      // 	if (k[18] > 40.9589) contains_otherpulse = 18;
      // 	if (k[19] > 40.4873) contains_otherpulse = 19;
      // 	if (k[20] > 40.6747) contains_otherpulse = 20;
      // 	if (k[21] > 40.2588) contains_otherpulse = 21;
      // 	if (k[22] > 40.2063) contains_otherpulse = 22;
      // 	if (k[23] > 40.2976) contains_otherpulse = 23;
      // 	if (k[24] > 41.4719) contains_otherpulse = 24;
      // 	if (k[25] > 41.1702) contains_otherpulse = 25;
      // 	if (k[26] > 41.5278) contains_otherpulse = 26;
      // }    

      // if (run_id==276) {
      // 	double lim = 3*spe;
      // 	if (k[10]>lim) contains_otherpulse = 10;
      // 	if (k[12]>lim) contains_otherpulse = 12;
      // 	if (k[13]>lim) contains_otherpulse = 13;
      // 	if (k[14]>lim) contains_afterpulse = 14;
      // 	if (k[15]>lim) contains_afterpulse = 15;
      // 	if (k[16]>lim) contains_afterpulse = 16;
      // 	if (k[17]>lim) contains_otherpulse = 17;
      // 	if (k[18]>lim) contains_otherpulse = 18;
      // 	if (k[19]>lim) contains_otherpulse = 19;
      // 	if (k[20]>lim) contains_otherpulse = 20;
      // 	if (k[21]>lim) contains_otherpulse = 21;
      // 	if (k[22]>lim) contains_otherpulse = 22;
      // 	if (k[23]>lim) contains_otherpulse = 23;
      // 	if (k[24]>lim) contains_otherpulse = 24;
      // 	if (k[25]>lim) contains_otherpulse = 25;
      // 	if (k[26]>lim) contains_otherpulse = 26;
      // 	if (k[27]>lim) contains_otherpulse = 27;
      // 	if (k[28]>lim) contains_otherpulse = 28;
      // 	if (k[29]>lim) contains_otherpulse = 29;
      // 	if (k[30]>lim) contains_otherpulse = 30;
      // }
      // if (verbose) std::cout<<"Other Pulse: "<<contains_otherpulse<<" After Pulse: "<<contains_afterpulse<<std::endl;
      // //    if (!run_id==247) {
      // if (contains_otherpulse) return 0;
      // if (contains_afterpulse && !include_errors_from_afterpulse) return 0;
      // // }
      // // if (run_id==247 && contains_afterpulse) if (verbose) std::cout<<"Should have cut this event!!!!!!!!!!"<<std::endl;
      // if (verbose) std::cout<<"Passed Afterpulse"<<std::endl;

    }// End run-specific cuts

    //warning: add the cuts for TOF -- for NaI PSD only
    //only apply the cut when tof cuts are specified by the users
    if(max_tof>-100 && min_tof>-100){
      ChannelData *refch = _current_event->GetEventData()->GetChannelByID(7);
      if(!refch) return 0;
      //  if (verbose) std::cout<<chdata->channel_id<<": "<<chdata->tof.size()<<", "<<refch->tof.size()<<std::endl;
      if(chdata->tof.size()<1 || refch->tof.size()<2) return 0;
      //  if (verbose) std::cout<<'\t'<<chdata->tof[0].half_max_time-refch->tof[0].half_max_time<<std::endl;
      if(chdata->tof[0].time_since_last<0.2 || //psd parameters? 
	 chdata->tof[0].half_max_time-refch->tof[0].half_max_time<min_tof ||
	 chdata->tof[0].half_max_time-refch->tof[0].half_max_time>max_tof) return 0;
      pulse = chdata->tof.at(0);
      if(!(pulse.evaluated) || pulse.integral>=0) return 0; //integral is negative
      //  if (verbose) std::cout<<"candidate event!"<<std::endl;
      // ChannelData *nch = _current_event->GetEventData()->GetChannelByID(2);
      // if(!nch || nch->tof.size()<1) return 0;
      // //  if (verbose) std::cout<<"candidate event come here!"<<std::endl;
      // //this cut is for channel 5 in Run3
      // if(nch->tof[0].time_since_last<0.2 || nch->tof[0].peak_amplitude<100 || //nch->tof[0].t90<0.3 || //this can be dangerous
      //    nch->tof[0].half_max_time-refch->tof[0].half_max_time<-0.01 ||
      //    nch->tof[0].half_max_time-refch->tof[0].half_max_time>0.05) return 0;
      // if (verbose) std::cout<<"candidate event accepted!"<<std::endl;
      // if (verbose) std::cout<<'\t'<<chdata->tof[0].half_max_time-refch->tof[0].half_max_time<<std::endl;
    }
    //  get rid of pileup events
    //    if(pulse.min<1 //|| pulse.peak_amplitude>-0.08*pulse.integral
       //|| pulse.peak_amplitude<chdata->baseline.mean-pulse.min-0.5
      // ) return 0;
    //    if (verbose) std::cout<<"k"<<std::endl;
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
      if(bin>=0 && bin<chdata->nsamps){
	double value = wave[bin]/(-65536.);
	y_array[i]  += value;
	ey_array[i] += value*value/pulse_integral;
      }
      //      std::cout<<chdata->SampleToTime(bin)<<'\t'<<y_array[i]<<'\t'<<ey_array[i]<<std::endl;
    }//end for
    _sum[graph_id]  += pulse_integral;
    _sum2[graph_id] += pulse_integral*pulse_integral;
    _total_waveforms[graph_id]++;
  }
  
  return 0;
  
}

*/
