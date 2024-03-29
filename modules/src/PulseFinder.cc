#include "PulseFinder.hh"
#include "EvalPulse.hh"
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

std::ostream& operator<<(std::ostream& out, const PulseFinder::SEARCH_MODE& m)
{
  switch(m){
  case PulseFinder::VARIANCE:
    out<<"VARIANCE";
    break;
  case PulseFinder::DISCRIMINATOR:
    out<<"DISCRIMINATOR";
    break;
  case PulseFinder::INTEGRAL:
    out<<"INTEGRAL";
    break;
  case PulseFinder::CURVATURE:
    out<<"CURVATURE";
    break;
  }
  return out;    
}

std::istream& operator>>(std::istream& in, PulseFinder::SEARCH_MODE& m)
{
  std::string dummy;
  in>>dummy;
  if(dummy == "VARIANCE" || dummy == "variance")
    m = PulseFinder::VARIANCE;
  else if (dummy == "DISCRIMINATOR" || dummy == "discriminator")
    m = PulseFinder::DISCRIMINATOR;
  else if (dummy == "INTEGRAL" || dummy == "integral")
    m = PulseFinder::INTEGRAL;
  else if (dummy == "CURVATURE" || dummy == "curvature")
    m = PulseFinder::CURVATURE;
  else{
    throw std::invalid_argument(dummy+"is not a valid value for search_mode!");
  }
  return in;
}


PulseFinder::PulseFinder() : 
  BaseModule(GetDefaultName(), "Search for individual physical scintillation pulses within the hardware trigger")
{
  AddDependency<ConvertData>();
  AddDependency<BaselineFinder>();
  AddDependency<Integrator>();
  
  ///@todo Provide helptext for PulseFinder parameters
  RegisterParameter("search_mode", mode=DISCRIMINATOR);
  RegisterParameter("search_start_time", search_start_time = -10000,
		    "time in us to start searching for pulses");
  RegisterParameter("search_end_time", search_end_time = 10000,
		    "time in us to end searching for pulses");
  ///parameters for curvature search
  RegisterParameter("start_window", start_window = 5);
  RegisterParameter("min_start_variance", min_start_variance = 100);
  RegisterParameter("min_resolution", min_resolution = 300);
  ///parameters for discriminator search
  RegisterParameter("discriminator_relative", discriminator_relative = true,
		    "Is the discriminator value relative to the baseline?");
  RegisterParameter("use_baseline_sigma", use_baseline_sigma = false,
		    "use the baseline fluctuation as discriminator level");
  RegisterParameter("discriminator_value", discriminator_value = -5,
		    "Value sample must cross to mark the start of the pulse");
  RegisterParameter("discriminator_nsigma", discriminator_nsigma = 4,
		    "Value sample must cross to mark the start of the pulse");
  RegisterParameter("discriminator_start_add", discriminator_start_add = 2,
		    "Number of samples to add before the detected start of the pulse");
  RegisterParameter("discriminator_end_add", discriminator_end_add = 2,
		    "Number of samples to add after the detected end of the pulse");
  RegisterParameter("discriminator_combine_overlapping_pulses", discriminator_combine_overlapping_pulses = false,
		    "Combine overlapping pulses into a single pulse");
  ///parameters for integral search
  // RegisterParameter("normalize_to_npe", normalize_to_npe=true,
  // 		    "normalize amplitude to spe before searching?");
  // RegisterParameter("integral_start_time", integral_start_time = 3, 
  // 		    "time in us over which photons arrive");
  // RegisterParameter("integral_end_time", integral_end_time = 3,
  // 		    "time over which to evaluate end of pulse to return to 0");
  // RegisterParameter("integral_start_threshold", integral_start_threshold = 5,
  // 		    "amount in npe integral must fall");
  // RegisterParameter("integral_end_threshold", integral_end_threshold = 3,
  // 		    "amount in npe integral must fall");
  // RegisterParameter("min_sep_time", min_sep_time = 2, 
  // 		    "minimum time between starts of two pulses");
  // RegisterParameter("multipulse_thresh_value", multipulse_thresh_value = 2,
  // 		    "secondary must be > this*previous integral+start_thresh");
  RegisterParameter("amplitude_start_threshold",amplitude_start_threshold = 0.3,
    		    "Raw signal must go below this at actual pulse start");  
  // RegisterParameter("amplitude_end_threshold", amplitude_end_threshold = 0.3,
  // 		    "Amplitude must fall below this to end pulse");
  // RegisterParameter("min_pulse_time", min_pulse_time = 7,
  // 		    "Minimum pulse width before it can be considered over");
  // RegisterParameter("lookback_time", lookback_time = 0.5,
  // 		    "Time to compare against for pileup");
  RegisterParameter("integral_window_nsamps", integral_window_nsamps = 20,
                    "number of samples in the moving integral window");
  RegisterParameter("moving_step_nsamps", moving_step_nsamps = 5,
                    "step to increase for the moving window");
  RegisterParameter("pileup_factor_rel", pileup_factor_rel = 2,
                    "pileup factor in multiple of previous sum");
  RegisterParameter("pileup_factor_abs", pileup_factor_abs = 30,
                    "pileup factor in the unit of adc samples");
  RegisterParameter("integral_start_threshold", integral_start_threshold = 20,
                    "amount in adc samples in the moving window");
  RegisterParameter("integral_end_threshold", integral_end_threshold = 20,
                    "amount in adc samples in the moving window");

  // Curvature Search
  RegisterParameter("down_sample_factor", down_sample_factor = 250,
                    "Reduce the integral vector size by this factor");
  RegisterParameter("pulse_start_curvature", pulse_start_curvature = -4,
                    "Curvature threshold to start a new pulse");
  RegisterParameter("pile_up_curvature", pile_up_curvature = -20,
                    "Curvature threshold to start a pile up pulse");
  RegisterParameter("pulse_end_slope", pulse_end_slope = -0.25,
                    "Slope threshold to end a pulse");

  // Afterpulse search
  RegisterParameter("afterpulse_RMS_factor", afterpulse_RMS_factor = 8,
                    "RMS factor to tag an afterpulse");
  RegisterParameter("afterpulse_window_nsamps", afterpulse_window_nsamps = 50,
                    "Number of samples used to calculate the RMS for afterulses");
  RegisterParameter("afterpulse_search_start", afterpulse_search_start = 0.3,
                    "Start time for afterpulse search");
  RegisterParameter("afterpulse_search_end", afterpulse_search_end = 3.0,
                    "End time for afterpulse search");
  RegisterParameter("cherenkov_discriminator_value", cherenkov_discriminator_value = -150,
                    "Decriminator value for cherenkov search");
  RegisterParameter("cherenkov_search_start", cherenkov_search_start = 3.0,
                    "Start time for cherenkov search");
  RegisterParameter("cherenkov_search_end", cherenkov_search_end = 99999.,
                    "End time for cherenkov search");

}
  
PulseFinder::~PulseFinder()
{

}

int PulseFinder::Initialize()
{
  return 0;
}

int PulseFinder::Finalize()
{
  return 0;
}


double PulseFinder::windowMean(double* wave, int start_index, int lookback, int nsamps) {
  int n=0;
  double v=0.;
  for (int i=start_index-lookback+1; i<=start_index; i++) {
    if (i<0 || i>=nsamps) continue;
    v+=wave[i];
    n++;
  }
  return n==0 ? wave[0] : v/n;
}

double PulseFinder::windowRMS(double* wave, int start_index, int lookback, int nsamps) {
  int n=0;
  double v=0.;
  for (int i=start_index-lookback+1; i<=start_index; i++) {
    if (i<0 || i>=nsamps) continue;
    double mean = windowMean(wave, start_index, lookback, nsamps);
    v+=(wave[i]-mean)*(wave[i]-mean);
    n++;
  }
  return n==0 ? 0 : sqrt(v/n);
}
/*
void PulseFinder::smoothArray(double& array, int nsamps, int smooth_radius=10) {
  double smooth_array = new double[nsamps];
  for (int i=0; i<nsamps; i++) {
    smooth_array[i] = windowMean(array, i+smooth_radius, smooth_radius, nsamps);
  }
  for (int i=0; i<nsamps; i++) {
    array[i] = smooth_array[i];
  }
  delete smooth_array;
  return;
}
*/
int PulseFinder::Process(EventPtr evt)
{
    EventDataPtr event = evt->GetEventData();

    //Start search for pulse edges (start and end) ************************************************************************************
    std::vector<int> start_index;
    std::vector<int> end_index;
    
    //Loop over all channels and evaluate pulse edges individually
    for (size_t ch = 0; ch < event->channels.size(); ch++)
      {
	ChannelData& chdata = event->channels[ch];
	//skip channels we've been told to explicitly skip
	if(_skip_channels.find(chdata.channel_id) != _skip_channels.end())
	  continue;
	if(! chdata.baseline.found_baseline)
	  continue;
	start_index.clear();
	end_index.clear();
	
	if(mode == DISCRIMINATOR)
	  DiscriminatorSearch(&chdata, start_index, end_index);
	else if (mode == VARIANCE)
	  VarianceSearch(&chdata, start_index, end_index);
	else if (mode == INTEGRAL)
	  IntegralSearch(&chdata, start_index, end_index);
	else if (mode == CURVATURE)
	  CurvatureSearch(&chdata, start_index, end_index);
	
	int search_start_index=chdata.TimeToSample(search_start_time, true);
	for (size_t i = 0; i < start_index.size();  i++)
	{
	  if (start_index[i] >= end_index[i]) {
	    std::cout<<i<<"/"<<start_index.size()<<" "<<start_index[i]<<" "<<end_index[i]<<std::endl;
	    return -1;
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


	// Search for afterpulses
	double* wave = chdata.GetBaselineSubtractedWaveform();
	  //	Baseline& baseline = chdata->baseline;
	const int nsamps = chdata.nsamps;
	// for(size_t p=0; p<chdata.pulses.size(); p++){
	//   std::cout<<"Pulse: "<<chdata.pulses.at(p).start_time<<" "<<chdata.pulses.at(p).end_time<<std::endl;
	// }
	/*
	// New method 2
	for(size_t p=0; p<chdata.pulses.size(); p++){
	  Pulse pulse = chdata.pulses.at(p);
	  int search_start_index = chdata.TimeToSample(pulse.start_time+afterpulse_search_start);
	  int search_end_index = chdata.TimeToSample(pulse.start_time+afterpulse_search_end);
	  const int afterpulse_window_nsamps = search_end_index-search_start_index;
	  // Make a waveform of the afterpulse search region
	  double ap_wave[afterpulse_window_nsamps];
	  for (int i=0; i<afterpulse_window_nsamps; i++) {
	    int wave_index = i+search_start_index;
	    ap_wave[i] = wave[wave_index];
	  }
	  // Downsample
	  const int ds_factor = 10;
	  const int ds_nsamps = afterpulse_window_nsamps/ds_factor+1;
	  double ds_ap_wave[ds_nsamps];
	  for (int i=0; i<ds_nsamps; i++) {
	    ds_ap_wave[i] = meanWindow(ap_wave, (i+1)*ds_factor, ds_factor, afterpulse_window_nsamps);
	  }
	  
	}


	// New method
	for(size_t p=0; p<chdata.pulses.size(); p++){
	  Pulse pulse = chdata.pulses.at(p);
	  int search_start_index = chdata.TimeToSample(pulse.start_time+afterpulse_search_start);
	  int search_end_index = chdata.TimeToSample(pulse.start_time+afterpulse_search_end);
	  const int afterpulse_window_nsamps = search_end_index-search_start_index;
	  // Make a smoothed out waveform of the afterpulse search region
	  double smooth_wave[afterpulse_window_nsamps];
	  for (int i=0; i<afterpulse_window_nsamps; i++) {
	    int wave_index = i+search_start_index;
	    smooth_wave[i] = wave[wave_index];
	  }
	  // Smooth the array many times
	  smoothArray(smooth_wave,nsamps);
	  smoothArray(smooth_wave,nsamps);
	  smoothArray(smooth_wave,nsamps);
	  smoothArray(smooth_wave,nsamps);
	  smoothArray(smooth_wave,nsamps);
	  for (int i=0; i<afterpulse_window_nsamps; i++) {
	    int wave_index = i+search_start_index;
	    
	  }
	  // Subract off the pulse from the waveform, now we just get the fluctuations
	  double pulse_subtracted_wave[afterpulse_window_nsamps];
	  for (int i=0; i<afterpulse_window_nsamps; i++) {
            int wave_index = i+search_start_index;
            pulse_subtracted_wave[i] = wave[wave_index]-smooth_wave[i];
          }
	}

	for(size_t p=0; p<chdata.pulses.size(); p++){
	  Pulse pulse = chdata.pulses.at(p);
	  int search_start_index = chdata.TimeToSample(pulse.start_time+afterpulse_search_start);
	  int search_end_index = chdata.TimeToSample(pulse.start_time+afterpulse_search_end);
	  for (int i=search_start_index; i<search_end_index; i++) {
	    if (i>nsamps) continue;
	    double preMean = windowMean(wave,i-10,afterpulse_window_nsamps,nsamps); // RMS of the waveforms in the previous window_nsamps
	    double preRMS = windowRMS(wave,i-10,afterpulse_window_nsamps,nsamps); // RMS of the waveforms in the previous window_nsamps
	    if ((preMean-wave[i])>afterpulse_RMS_factor*preRMS) {
	      pulse.has_afterpulse = true;
	      std::cout<<"afterpulse at "<<chdata.SampleToTime(i)<<" preMean: "<<preMean<<" preRMS: "<<preRMS<<" thisValue: "<<wave[i]<<std::endl;
	    }
	  }
	}
	*/
	// Search for cherenkov
	for(size_t p=0; p<chdata.pulses.size(); p++){
	  Pulse& pulse = chdata.pulses.at(p);
	  int search_start_index = chdata.TimeToSample(pulse.start_time+cherenkov_search_start);
	  int search_end_index = chdata.TimeToSample(pulse.start_time+cherenkov_search_end);
	  if (search_end_index>pulse.end_index) search_end_index = pulse.end_index;
	  for (int i=search_start_index; i<search_end_index && i<nsamps && !pulse.has_cherenkov; i++) {
	    if (fabs(wave[i])>fabs(cherenkov_discriminator_value)) {
	      pulse.has_cherenkov = true;
	      //	      std::cout<<"pulse "<<p<<" cherenkov at "<<chdata.SampleToTime(i)<<" thisValue: "<<wave[i]<<" discrim: "<<cherenkov_discriminator_value<<std::endl;
	    }
	  }
	}
	//get rid of noise pulses
	//	double noise_start=0, noise_end=0;
	// for(size_t j=0; j<chdata.pulses.size(); j++){
	//   Pulse &pulse = chdata.pulses.at(j);

	//   if(pulse.integral+pulse.peak_amplitude>=0 ||
	//      //	     pulse.overshoot>5*chdata.baseline.sigma ||
	//      pulse.overshoot>0.5*pulse.peak_amplitude){
	//     pulse.evaluated=false;
	//     /*	    noise_start = pulse.start_time-0.2;
	//     noise_end = pulse.end_time+0.2;
	//     //get rid of the previous noise pulses
	//     for(int k=j-1; k>=0; k--){
	//       if(chdata.pulses.at(k).start_time>=noise_start &&
	// 	 chdata.pulses.at(k).start_time<=noise_end) chdata.pulses.at(k).evaluated=false;
	//       else break;
	//       }*/
	//   }//end if this is a bad pulse
	//   /*	  if(noise_end>noise_start){//if this pulse comes close after a noise
	//     if(chdata.pulses.at(j).start_time>=noise_start &&
	//        chdata.pulses.at(j).start_time<=noise_end) chdata.pulses.at(j).evaluated=false;
	//     else {noise_start=0; noise_end=0;}//out of noise window
	//     }//end if noise end*/
	// }//end for loop

	chdata.npulses = chdata.pulses.size();
    } //end loop over channels
    //End evaluation of pulse variables for each pulse on each channel*************************************************************
    //      std::cout<<"Ending pulsefinder"<<std::endl;

    return 0;
}

void PulseFinder::VarianceSearch(ChannelData* chdata, 
				 std::vector<int>& start_index,
				 std::vector<int>& end_index)
{
  Baseline& baseline = chdata->baseline;
  int index=start_window;
  double* wave =  chdata->GetWaveform();
  double start_baseline;
  bool found_start;
  for(index = start_window; index < chdata->nsamps; index++)
    {
      if(start_index.size() > 5)
	{
	  break;
	}
      //look for the starts of pulses 
      //pulse must decrease by more than baseline variance for start_window 
      //consecutive samples, and be less than start_baseline-minvariance*var at end
      
      found_start = true;
      for(int samp=index-start_window; samp < index; samp++)
	{
	  if( wave[samp]-wave[samp+1] < baseline.sigma )
	    {
	      found_start = false;
	      break;
	    }
	}
      if (start_index.size() == 0)
	start_baseline = baseline.mean;
      else
	start_baseline = wave[index-start_window];
      
      if(!found_start || 
	 wave[index] > start_baseline-min_start_variance * baseline.sigma )
	continue;
      //if we get here, we have found a start point
      start_index.push_back(index-start_window);
      index = index + min_resolution;
    }//end loop through index
  
  //Look for ends of pulses
  
  for (size_t i = 0; i < start_index.size(); i++)
    {
      int limit_index;
      if (i == start_index.size() - 1)
	limit_index = chdata->nsamps;
      else
	limit_index = start_index[i+1];
      
      int min_index = std::min_element(wave + start_index[i], wave + limit_index) - wave;
      if(wave[min_index] > baseline.mean)
	end_index.push_back(min_index);
      else
	end_index.push_back(limit_index - 1);
    }
  
}

void PulseFinder::DiscriminatorSearch( ChannelData* chdata,
				       std::vector<int>& start_index,
				       std::vector<int>& end_index){
  
  double check_val = discriminator_value;
  if(use_baseline_sigma){
    check_val = -1.*std::abs(discriminator_nsigma)*chdata->baseline.sigma;
    if(!discriminator_relative) check_val += chdata->baseline.mean;
    if(check_val>discriminator_value) check_val=discriminator_value;
  }
  const double* wave = chdata->GetBaselineSubtractedWaveform();
  if(!discriminator_relative){
    wave = chdata->GetWaveform();
  }
  
  bool in_pulse = false;
  int search_start_index = chdata->TimeToSample(search_start_time, true);
  int search_end_index = chdata->TimeToSample(search_end_time, true);
  for(int index = search_start_index+discriminator_start_add; 
      index<search_end_index// -discriminator_end_add
	; index++){
    if(wave[index]<=check_val){
      if(wave[index-1]>check_val){//if just come to a pulse
        in_pulse = true;
	int this_start_index = index-discriminator_start_add;
	if(this_start_index<0)
	  this_start_index = 0;
	start_index.push_back(this_start_index);
	   //	std::cout<<"New Pulse starting at "<<chdata->SampleToTime(index)<<" wave[index]: "<<wave[index]<<" check_val: "<<check_val<<std::endl;
      }
      if(wave[index+1]>check_val){//if just leave a pulse
	int this_end_index = index+discriminator_end_add;
	if (this_end_index > chdata->nsamps-1) this_end_index = chdata->nsamps-1;
        if(in_pulse) end_index.push_back(this_end_index);
        in_pulse = false;
	//	std::cout<<"Pulse ending at "<<index+discriminator_end_add<<"/"<<search_end_index<<" wave[index]: "<<wave[index]<<" check_val: "<<check_val<<std::endl;
      }//end if wave index+1
    }//end if wave index <=
  }//end for int index
  if(in_pulse) end_index.push_back(search_end_index);

  // if(start_index.size()!=end_index.size()){
  //   std::cerr<<start_index.size()<<" =/= "<<end_index.size()<<std::endl;
  //   std::cerr<<"PulseFinder failed to find correct pulses, give up this event."<<std::endl;
  // start_index.clear();
  // end_index.clear();
  //    return;
  //  }
  // Resoving overlapping pulses
  if (start_index.size()<=1) return;
  if (discriminator_combine_overlapping_pulses) {
    std::vector<int> combined_start_index;
    std::vector<int> combined_end_index;
    for(size_t index=0; index<start_index.size();){
      combined_start_index.push_back(start_index.at(index));
      //      std::cout<<"New Start: "<<start_index.at(index)<<std::endl;
      while(++index && index<start_index.size() && start_index.at(index)<end_index.at(index-1));
      combined_end_index.push_back(end_index.at(index-1));
      //      std::cout<<"New End: "<<end_index.at(index-1)<<std::endl;
    }
    start_index.clear();
    end_index.clear();
    start_index = combined_start_index;
    end_index = combined_end_index;
  } else {
    // for(size_t index=0; index<start_index.size(); index++){
    //   std::cout<<"Pulse "<<index<<": "<<start_index.at(index)<<"-"<<end_index.at(index)<<std::endl;
    // }
    // Drop afterpulses
    for(size_t index=1; index<start_index.size(); ){
      double prev_start_t = chdata->SampleToTime(start_index.at(index-1));
      double start_t = chdata->SampleToTime(start_index.at(index));
      //      std::cout<<prev_start_t<<" "<<start_t<<" "<<afterpulse_search_start<<" "<<afterpulse_search_end<<std::endl;
      if (start_t > prev_start_t+afterpulse_search_start && start_t < prev_start_t+afterpulse_search_end) { 
	start_index.erase(start_index.begin() + index);
	end_index.erase(end_index.begin() + index);
      } else index++;
    }
    

    // Simple method: If pulses are overlapping, set end of early pulse to beginning of later pulse
    for(size_t index=1; index<start_index.size(); index++){
      if(start_index.at(index)<end_index.at(index-1)){
	end_index.at(index-1)=start_index.at(index);
      }
    }
    // Remove pulses that end before they start
    for(size_t index=0; index<start_index.size(); ){
      if(end_index.at(index)<start_index.at(index)) {
	start_index.erase(start_index.begin() + index);
	end_index.erase(end_index.begin() + index);
      } else index++;
    }
    
    
    // for(size_t index=0; index<start_index.size(); index++){
    //   std::cout<<"New Pulse "<<index<<": "<<start_index.at(index)<<"-"<<end_index.at(index)<<std::endl;
    // }
    // Advanced method
    // for(size_t index=1; index<start_index.size(); index++){
    //   // std::cout<<chdata->channel_id<<'\t'<<index<<'\t'<<chdata->SampleToTime(start_index.at(index))
    //   // 	     <<'\t'<<chdata->SampleToTime(end_index.at(index))<<std::endl;
    //   if(start_index.at(index)<=end_index.at(index-1)){
    // 	int middle_index = (start_index.at(index)+end_index.at(index-1))/2;
    // 	double max_value = wave[middle_index];
    // 	//      for(int j=start_index.at(index); j<=end_index.at(index-1); j++){
    // 	for(int j=end_index.at(index-1)-discriminator_end_add+1; j<start_index.at(index)+discriminator_start_add; j++){
    // 	  if(wave[j]<=max_value) continue;
    // 	  max_value = wave[j];
    // 	  middle_index = j;
    // 	}//end for j loop
    // 	start_index.at(index) = middle_index;
    // 	end_index.at(index-1) = middle_index;
    //   }//end if statement
    // }//end for index loop
  }

  // for(size_t index=0; index<start_index.size(); index++)
  //   std::cout<<chdata->channel_id<<'\t'<<index<<'\t'<<chdata->SampleToTime(start_index.at(index))
  //    	     <<'\t'<<chdata->SampleToTime(end_index.at(index))<<std::endl;
  
  //get rid of really small/fake pulses
  // for(int index=0; index<start_index.size()+0.; index++){
  //   if(end_index.at(index)-start_index.at(index)>discriminator_start_add+discriminator_end_add)
  //     continue;
  //   if(index>0 && end_index.at(index-1)>=start_index.at(index)-discriminator_end_add) 
  //     end_index.at(index-1) = end_index.at(index);
  //   else if(index<start_index.size()-1. && start_index.at(index+1)<=end_index.at(index)+discriminator_start_add) 
  //     start_index.at(index+1) = start_index.at(index);
  //   start_index.erase(start_index.begin()+index);
  //   end_index.erase(end_index.begin()+index);
  //   index --;
  // }
  
  return;
}

//jingke's integral search algorithm for short pulses
void PulseFinder::IntegralSearch( ChannelData* chdata,
				  std::vector<int>& start_index,
				  std::vector<int>& end_index){
  
  const double* integral = chdata->GetIntegralWaveform();
  std::vector<bool> pulse_pileup;
  double prev_moving_sum = 0;
  bool in_pulse = false;
  int search_start_index = chdata->TimeToSample(search_start_time, true);
  int search_end_index = chdata->TimeToSample(search_end_time, true);
  for(int index = search_start_index; 
      index<=search_end_index-integral_window_nsamps; index++){

    //look for pulses with step of moving_step_nsamps
    if(index%moving_step_nsamps) continue;
    //make the sum positive, so it is less confusing
    double moving_sum = integral[index]-integral[index+integral_window_nsamps];
    if(!in_pulse && moving_sum>integral_start_threshold){
      in_pulse = true;
      start_index.push_back(index);
      pulse_pileup.push_back(false);
    }
    else if(in_pulse && moving_sum<integral_end_threshold){
      in_pulse = false;
      end_index.push_back(index+integral_window_nsamps);
    }
    else if(in_pulse && (moving_sum>prev_moving_sum*pileup_factor_rel+pileup_factor_abs)
            && (index-start_index.back()>integral_window_nsamps*2)){
      end_index.push_back(index);
      start_index.push_back(index);
      pulse_pileup.back()=true;
      pulse_pileup.push_back(true);
    }
    prev_moving_sum = moving_sum;
  }
  if(in_pulse) end_index.push_back(search_end_index);

  //check if there is a problem
  // if(start_index.size()!=end_index.size()){
  //   std::cerr<<"PulseFinder failed to find correct pulses, give up this even."<<std::endl;
  //   start_index.clear();
  //   end_index.clear();
  //   return;
  // }

  //remove fake pulses and separate overlapping pulses
  //index can go negative by index--
  const double* wave = chdata->GetBaselineSubtractedWaveform();
  for(int index=0; index<start_index.size()-0.; index++){
    //this is a fake pulse integral/length<0.1
    if(integral[start_index.at(index)]-integral[end_index.at(index)]
       < 0.1*(end_index.at(index)-start_index.at(index))){
      start_index.erase(start_index.begin()+index);
      pulse_pileup.erase(pulse_pileup.begin()+index);
      end_index.erase(end_index.begin()+index);
      index --;
      continue;
    }
    //this pulse overlaps with the previous one    
    if(index>0 && start_index.at(index)<end_index.at(index-1)){
      int middle_index = (start_index.at(index)+end_index.at(index-1))/2;
      double max_value = wave[middle_index];
      for(int j=end_index.at(index-1)-integral_window_nsamps+1; 
	  j<start_index.at(index); j++){                              
	if(wave[j]<max_value) continue;
	max_value = wave[j];
	middle_index = j;
      }//end for j loop
      start_index.at(index) = middle_index;
      end_index.at(index-1) = middle_index;
    }//end if statement

  }//end for index loop
  
  return;
}

/*
void PulseFinder::IntegralSearch( ChannelData* chdata,
				  std::vector<int>& start_index,
				  std::vector<int>& end_index)
{
  double scale_factor = normalize_to_npe ? chdata->spe_mean : 1;
  
  double* integral = chdata->GetIntegralWaveform();
  double* wave = chdata->GetBaselineSubtractedWaveform();
  int start_samps = (int)(integral_start_time * chdata->sample_rate); 
  int end_samps = (int)(integral_end_time * chdata->sample_rate);
  int min_pulse_samps = (int)(min_pulse_time * chdata->sample_rate);
  if(start_samps <= 0) start_samps = 1;
  int samp = -1;
  //threshold is updated continuously if we are in a pulse based on the 
  // expected arrival of photons
  bool in_pulse = false;
  while ( ++samp < chdata->nsamps){
    int lookback_samps = std::min(samp,
				  (int)(lookback_time*chdata->sample_rate));
    double int_thresh = integral_start_threshold;
    double amp_thresh = amplitude_start_threshold;
    //actions are different depending on whether we are already in a pulse
    if(in_pulse){
      //first, test to see if we've hit the end of the pulse
      if(samp - start_index.back() >= min_pulse_samps){ //is it long enough?
	int test_samp = std::min(chdata->nsamps-1, samp+end_samps);
	if( (integral[samp] - integral[test_samp] )/scale_factor < 
	    integral_end_threshold  && 
	    (*std::min_element(wave+samp,wave+test_samp))/scale_factor > 
	    -amplitude_end_threshold ){
	  //we are at the end of this pulse
	  in_pulse = false;
	  end_index.push_back(samp);
	  continue;
	}
      }
      
      //now check for pileup      
      //see if the amplitude is generally still increasing - if so skip
      if(integral[samp-lookback_samps] - integral[samp] > 
	 (integral[samp-2*lookback_samps] - integral[samp-lookback_samps]) * 
	 multipulse_thresh_value)
	continue;
      //look for an excursion greater than the maximum in lookback samps
      //and with integral > max expected over that region
      double prev_max = -(*std::min_element(wave+samp-lookback_samps,
					    wave+samp)) / scale_factor;
      //ratio between lookback and next must be ratio of sizes * factor
      double halfback =  (integral[samp-lookback_samps/2] - integral[samp] );
      double ratio = halfback / 
	(integral[samp-lookback_samps] - 
	 integral[samp-lookback_samps/2] ) ;
      
      double integral_est = halfback*ratio*
	(2*start_samps/lookback_samps ) / scale_factor;
      amp_thresh += prev_max*multipulse_thresh_value;
      int_thresh += integral_est*multipulse_thresh_value;
    }

    //first look for something that crosses the signal threshold in amplitude
    if(wave[samp]/scale_factor > -amp_thresh) continue;
    //if we get here, signal is past threshold.  Make sure it's not isolated
    int test_samp = std::min(chdata->nsamps-1, samp+start_samps);
    if((integral[samp]-integral[test_samp])/scale_factor > 
       int_thresh){
      //we have found a pulse!
      int good_samp = std::max(0,samp-2);
      if(in_pulse)
	end_index.push_back(good_samp);
      in_pulse = true;
      start_index.push_back(good_samp);
      samp += (int)(min_sep_time * chdata->sample_rate);
      continue;
    }
  
  
  }//end loop over samples
  if(in_pulse)
    end_index.push_back(chdata->nsamps-1);
  
}
*/
void PulseFinder::CurvatureSearch( ChannelData* chdata,
                                   std::vector<int>& start_index,
                                   std::vector<int>& end_index) {
  double* integral = chdata->GetIntegralWaveform();

  int df = down_sample_factor;
  int n = chdata->nsamps/df;
  std::vector<double> sm(n);
  for(int i=0; i<n; i++){
    sm[i] = integral[i*df];
  }

  std::vector<double> diff(n);
  diff[0]= sm[1];
  for (int i=1; i<n-1; i++){
    diff[i] = sm[i+1]-sm[i-1];
  }

  std::vector<double> curve(n);
  curve[0] = diff[1];
  for (int i=1; i<n-1; i++){
    curve[i] = diff[i+1]-diff[i-1];
  }

  bool in_pulse = false;
  bool before_peak = true;  // before peak of diff
  int last = -1; // index of last local maximum on diff
  for (int i=0; i<n-1; i++){
    if (!in_pulse){
      if (curve[i] < pulse_start_curvature){
	in_pulse = true;
	int start = i*df;
	//while (integral[start+1] >= integral[start]){
	//start++;
	//}
	int loopcount=0;
	int maxloop = df;
	if(i<n-2) maxloop+= df;
	double* sub = chdata->GetBaselineSubtractedWaveform();
	while( ++loopcount<maxloop && -sub[start] < amplitude_start_threshold)
	  start++;
	start_index.push_back(start-2 > 0 ? start-2 : 0);
      }
    }
    else{ // in pulse
      if (before_peak){ // look for peak of diff
        if (curve[i] > 0){
          before_peak = false;
        }
      }
      else { // after peak of diff
        if (curve[i] < 0 && last < 0){ // keep track of last diff maximum
          last = i;
        }
        if (curve[i] > 0 && last > 0){ // last diff maximum not start of pileup
          last = -1;
        }

        if (curve[i] < pile_up_curvature){  // found pile up
          end_index.push_back(last*df);
          start_index.push_back(last*df);
          before_peak = true;
          last = -1;
        }
        else if (diff[i] > pulse_end_slope){ // pulse gradually ends
          end_index.push_back(i*df);
          in_pulse = false;
          before_peak = true;
          last = -1;
        }
      }
    }
  }
  if (in_pulse){
    end_index.push_back(chdata->nsamps-1);
  }  
}

