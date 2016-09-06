#include "PSD.hh"
//#include "ConvertData.hh"
//#include "BaselineFinder.hh"
//#include "PulseFinder.hh"
//#include "TOF.hh"
#include "Pulse.hh"
//#include "intarray.hh"
//#include "TGraphErrors.h"
#include "TFile.h"
#include "TH1F.h"
//#include "RootWriter.hh"
//#include "EventData.hh"

#include <algorithm>
#include <numeric>
#include <iostream>

PSD::PSD() : 
  ChannelModule(GetDefaultName(), "Average the waveform for each channel over the entire run and save to the output root file")
{
  //  AddDependency<ConvertData>();
  //  AddDependency<BaselineFinder>();
  //warning: the following dependencies are added for PSD of NaI only
  //  AddDependency<PulseFinder>();
  //  AddDependency<TOF>();
  //Register all the config handler parameters
  RegisterParameter("fparameter_times",  fparameter_times, "Time values to calculate the fparameters");
  RegisterParameter("tparameter_ratios", tparameter_ratios, "Ratio values to calculate the tparameters");
  //gatti parameter stuff
  RegisterParameter("calculate_gatti", calculate_gatti = false, "if we want to calculate the gatti parameter");
  RegisterParameter("gatti_channels", gatti_channels, "list of channels to process gatti calculation");
  RegisterParameter("start_time_us", start_time_us = 0, "time relative to pulse rise as start of waveform PSD");
  RegisterParameter("end_time_us", end_time_us = 3, "time relative to pulse rise as end of waveform PSD");
  RegisterParameter("bin_width_us", bin_width_us = 0.04, "the bin width for waveform PSD");
  RegisterParameter("wfm_file", wfm_file="cuts/PSD.root", "the root file containing the waveform tempaltes for gatti calculation");
  // RegisterParameter("",  = , "");
  weights = NULL;
  rebined_pulse = NULL;
}

PSD::~PSD()
{
  Cleanup();
}

void PSD::Cleanup()
{
  if(weights) delete weights; 
  if(rebined_pulse) delete rebined_pulse; 
}

int PSD::Initialize()
{ 

  //order the tparameter ratios, make sure they are in increasing order
  if(tparameter_ratios.size()>1)
    std::sort (tparameter_ratios.begin(), tparameter_ratios.end());

  //prepare for gatti parameter calculations
  //convention is to set time 0 to separate the bin lower than threshold
  //and the bin higher than threshold. best guess of threshold time
  if(!calculate_gatti) return 0;

  TFile fin(wfm_file.c_str());
  if(!(fin.IsOpen())||fin.IsZombie()){
    std::cerr<<wfm_file<<" (the Gatti templates file) is not set properly!"<<std::endl;
    calculate_gatti = false;
    return 1;
  }

  TH1F *ref_a = (TH1F*)fin.Get("ref_a");
  TH1F *ref_b = (TH1F*)fin.Get("ref_b");
  if(!ref_a || !ref_a){
    std::cerr<<"Reference histogram ref_a and ref_a can not be found !"<<std::endl;
    fin.Close();
    calculate_gatti = false;
    return 1;
  }
  
  if(ref_a->GetNbinsX()<10 || ref_a->GetNbinsX()!=ref_b->GetNbinsX() || ref_a->GetBinWidth(1)<=0 ||
     ref_a->GetBinWidth(1)!=ref_b->GetBinWidth(1) || ref_a->GetBinLowEdge(1)!=ref_b->GetBinLowEdge(1)){
    fin.Close();
    calculate_gatti = false;
    return 1;
  }
  //make sure bin_width is multiples of 4ns
  double bin_width=(int)(bin_width_us/ref_a->GetBinWidth(1));
  if(bin_width<1) bin_width=1;
  bin_width *= ref_a->GetBinWidth(1);
  //make sure start_time is multiples of 4ns
  double start_time = (int)(start_time_us/ref_a->GetBinWidth(1))*ref_a->GetBinWidth(1);
  if(start_time<ref_a->GetBinLowEdge(1)) start_time = ref_a->GetBinLowEdge(1);

  double end_time = end_time_us;
  if(end_time>ref_a->GetBinLowEdge(ref_a->GetNbinsX())+ref_a->GetBinWidth(1)) 
    end_time=ref_a->GetBinLowEdge(ref_a->GetNbinsX())+ref_a->GetBinWidth(1);
  end_time = start_time + (int)((end_time-start_time)/bin_width)*bin_width;
  if(end_time-start_time<bin_width){
    fin.Close();
    calculate_gatti = false;
    return 1;
  }

  // std::cout<<"Refined start time: "<<start_time<<", end time: "<<end_time
  // 	   <<", bin width: "<<bin_width<<std::endl;

  //copy the histogram content to new histos with required binning
  TH1F *rebin_a = new TH1F("rebin_a", "rebined reference a", 
			   int((end_time-start_time)/ref_a->GetBinWidth(1)), start_time, end_time);
  TH1F *rebin_b = new TH1F("rebin_b", "rebined reference b", 
			   int((end_time-start_time)/ref_b->GetBinWidth(1)), start_time, end_time);
  for(int bin=1; bin<=rebin_a->GetNbinsX(); bin++){
    int orig_bin = ref_a->FindBin(rebin_a->GetBinCenter(bin));
    rebin_a->SetBinContent(bin, ref_a->GetBinContent(orig_bin));
    rebin_b->SetBinContent(bin, ref_b->GetBinContent(orig_bin));
  }
  // std::cout<<"test rebin histo: "<<rebin_a->GetBinContent(rebin_a->FindBin(-0.002))
  // 	   <<", next: "<<rebin_a->GetBinContent(rebin_a->FindBin(0.002))
  // 	   <<", peak: "<<rebin_a->GetBinContent(rebin_a->GetMaximumBin())<<std::endl;
  //rebin the hitograms
  rebin_a->Rebin((int)(bin_width/ref_a->GetBinWidth(1)));
  rebin_b->Rebin((int)(bin_width/ref_b->GetBinWidth(1)));
  //scale to unit integral
  if(rebin_a->Integral()<=0 || rebin_b->Integral()<=0){
    delete rebin_a;
    delete rebin_b;
    calculate_gatti = false;
    return 1;
  }
  else{
    rebin_a->Scale(1./rebin_a->Integral());
    rebin_b->Scale(1./rebin_b->Integral());
  }
  
  weights = (TH1F*)rebin_a->Clone("PSD_weights");
  weights->SetDirectory(0);
  weights->SetTitle("PSD weights");
  rebined_pulse = (TH1F*)rebin_a->Clone("Rebined_Pulse");
  rebined_pulse->SetDirectory(0);
  rebined_pulse->SetTitle("Rebined Pulse");
  //calculate the gatti weights
  for(int bin=1; bin<=weights->GetNbinsX(); bin++){
    rebined_pulse->SetBinContent(bin, 0);
    double a = rebin_a->GetBinContent(bin);
    double b = rebin_b->GetBinContent(bin);
    if(a+b>0) weights->SetBinContent(bin, (a-b)/(a+b));
    else weights->SetBinContent(bin, 0);
  }
  weights->SaveAs("Gatti_Weights.root");

  delete rebin_a;
  delete rebin_b;
  fin.Close();
  return 0; 
}

int PSD::Finalize()
{return 0;}

//this is not a critical module, return 0 for all cases
int PSD::Process(ChannelData *chdata){

  if(!(chdata->baseline.found_baseline)) return 0;

  //process the regions of interest
  for(size_t i=0; i<chdata->regions.size(); i++){
    Pulse& region = chdata->regions.at(i);
    if(region.evaluated) ProcessPulse(chdata, region);
  }

  //do not process pulses -- not being saved

  //process the tofs
  for(size_t i=0; i<chdata->tof.size(); i++){
    Pulse & tof = chdata->tof.at(i);
    if(tof.evaluated) ProcessPulse(chdata, tof);
  }

  return 0;
}

//this is not a critical module, return 0 for all cases
void PSD::ProcessPulse(ChannelData *chdata, Pulse & pulse){

  if(chdata->integral.empty()) return;

  //todo: need to test all these new algorithms
  //calculate the fparameters 
  for(size_t count=0; count<fparameter_times.size(); count++){
    int fp_index = chdata->TimeToSample(pulse.half_max_time+fparameter_times.at(count), true);
    if(fp_index>pulse.end_index) fp_index = pulse.end_index;
    double fp = chdata->integral[fp_index] - chdata->integral[pulse.start_index];
    // std::cout<<"fp calculation: "<<fparameter_times.at(count)<<'\t'<<chdata->SampleToTime(fp_index)
    // 	     <<'\t'<<fp<<'\t'<<pulse.integral<<std::endl;
    if(pulse.integral!=0) pulse.fparameters.push_back(fp/pulse.integral);
    else pulse.fparameters.push_back(-1);
  }

  //calculate the tparameters
  //look for the time it takes to reach X% of total integral
  //initialize the vector with unphysical value -1
  for(size_t count=0; count<tparameter_ratios.size(); count++){
    pulse.tparameters.push_back(-1);
  }
  //start from the beginning
  size_t tp_count = 0;
  double tp_threshold = 0;
  if(tparameter_ratios.size()) 
    tp_threshold = pulse.integral*tparameter_ratios.at(tp_count); //negative
  //remember, integral is negative
  for(int samp = pulse.start_index;
      (samp<pulse.end_index) && (tp_count<pulse.tparameters.size()); samp++){
    if((chdata->integral[samp]-chdata->integral[pulse.start_index]<=tp_threshold)){
      //need a loop here in case there are ratios too close in value
      while((tp_count<pulse.tparameters.size()) && 
	    (chdata->integral[samp]-chdata->integral[pulse.start_index]<=tp_threshold)){
	pulse.tparameters.at(tp_count) = chdata->SampleToTime(samp)-pulse.start_time;
	tp_count ++;
	if(tp_count<pulse.tparameters.size()) tp_threshold = pulse.integral*tparameter_ratios.at(tp_count);
      }//end while
    }//end if
  }//end for
  while(tp_count<pulse.tparameters.size())pulse.tparameters.at(tp_count++) = pulse.end_time;

  const double* wave = chdata->GetBaselineSubtractedWaveform();

  //calculate the mean time for NaI only, use test parameter for now
  //  int half_max_index = chdata->TimeToSample(pulse.half_max_time, true);
  int prompt_end_index = chdata->TimeToSample(pulse.half_max_time+0.15, true);//make it 0.15 for now
  double mean_time=0, integral=0;
  double l_mean_time=0, l_integral=0;
  for(int samp = pulse.start_index+1; samp<pulse.end_index; samp++){
    //    if((samp-pulse.start_index)/1.5/chdata->sample_rate<1){
      mean_time += wave[samp]*(samp-pulse.start_index);
      integral += wave[samp];
      //    }
      //this is to calculate the late mean time
      if(samp > prompt_end_index){
	l_mean_time += wave[samp]*(samp-pulse.start_index);
	l_integral += wave[samp];
      }//end if
  }
  if(integral<0) pulse.mean_time = mean_time/chdata->sample_rate/integral;
  if(l_integral<0) pulse.test = l_mean_time/chdata->sample_rate/l_integral;
  // for(size_t i=0; i<chdata->pulses.size(); i++){
  //   Pulse& ps = chdata->pulses[i];
  //   if(ps.start_index<pulse.start_index) continue;
  //   if(ps.start_index>pulse.end_index) break;
  //   for(int samp = ps.start_index; samp<ps.end_index; samp++){
  //     mean_time += wave[samp]*(samp-half_max_index);
  //   }
  // }
  // if(pulse.integral<0) pulse.test = mean_time/chdata->sample_rate/pulse.integral;
  /*
  //here test is a parameter to reject pileup events
  int search_max_index = chdata->TimeToSample(pulse.start_time+0.5, true);
  //  std::cerr<<">>Tests: "<<pulse.start_index<<'\t'<<search_max_index<<'\t'<<pulse.end_index<<std::endl;
  if(search_max_index<pulse.end_index){
    int peak_index = std::min_element(wave + search_max_index, 
				      wave + pulse.end_index) - wave;
    //warning this is to test the noise rejection algorithm
    // pulse.test = *std::max_element(subtracted + pulse.start_index, 
    // 				 subtracted + max_peak_index);
//    pulse.tail_peak = -wave[peak_index]; //pulse are neg, amplitude pos
  }
  */
  //calculate the gatti parameters
  if(!calculate_gatti || !weights || !rebined_pulse) return;
  //only calculate gatti for certain channels
  if( gatti_channels.find( chdata->channel_id ) == gatti_channels.end())
    return;
  //  const double* wave = chdata->GetBaselineSubtractedWaveform();

  //start_index is the index before the trigger index, add 1 to get trigger index later
  //  int start_index = chdata->TimeToSample(pulse.half_max_time-0.5/chdata->sample_rate);
  int start_index = pulse.start_index -1;
  // std::cout<<"test pulse: "<<pulse.half_max_time<<", "<<start_index<<", "
  // 	   <<-wave[start_index]<<", next: "<<-wave[start_index+1]
  // 	   <<", peak: "<<-wave[pulse.peak_index]<<std::endl;
  start_index += (int)(start_time_us*chdata->sample_rate);
  int pulse_nsamps = (end_time_us-start_time_us)*chdata->sample_rate;
  int bin_width = (int)(bin_width_us*chdata->sample_rate);
  int bin=0;
  double sum=0;
  for(int samp=1; samp<=pulse_nsamps; samp++){
    bin = start_index + samp;
    // if(samp==1) std::cout<<"first sample time: "<<chdata->SampleToTime(bin)
    // 			 <<", rebined histo time: "<<rebined_pulse->GetBinLowEdge(1)<<std::endl;
    if(bin>=pulse.start_index && bin<=pulse.end_index) sum -= wave[bin];
    if(!(samp%bin_width)){
      rebined_pulse->SetBinContent(samp/bin_width, sum);
      sum = 0.;
    }
  }
  if(_current_event->GetEventData()->event_id==1058) rebined_pulse->SaveAs("rebined_pulse.root");
  //use sum to store the summed gatti parameter
  sum = 0;
  if(rebined_pulse->Integral()>0){
    rebined_pulse->Scale(1./rebined_pulse->Integral());
    for(bin=1; bin<=rebined_pulse->GetNbinsX(); bin++){
      sum += rebined_pulse->GetBinContent(bin)*weights->GetBinContent(bin);
    }
  }
  pulse.gatti = sum;

  rebined_pulse->Reset();  
  return;
  
}

