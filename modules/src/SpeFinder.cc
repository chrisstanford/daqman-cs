#include "SpeFinder.hh"
//#include "EventHandler.hh"
#include "ConvertData.hh"
//#include "RootWriter.hh"
#include "BaselineFinder.hh"
#include <vector>
#include <numeric>

SpeFinder::SpeFinder():
  ChannelModule(GetDefaultName(), "Search for single photoelectrons in the tails of pulses identified by PulseFinder")
{
  AddDependency<BaselineFinder>();
  RegisterParameter("search_start_time", search_start_time = 1.,
                    "Time in [us] to start searching for spe");
  RegisterParameter("search_end_time", search_end_time = 5.,
                    "Time in [us] to end searching for spe");
  RegisterParameter("threshold_nsigmas", threshold_nsigmas = -2.5,
                    "Threshold in the unit of baseline sigmas");
  RegisterParameter("pulse_edge_add", pulse_edge_add  = 3,
                    "Number of samples to be added at the edge of possible pulses");
  RegisterParameter("min_baseline_length", min_baseline_length = 0.04,
                    "Time in us of good baseline before and after pulse");
}

SpeFinder::~SpeFinder()
{
  Finalize();
}

int SpeFinder::Initialize()
{
  //make sure we don't process the sum channel
  _skip_sum = true;
  if(pulse_edge_add<1) pulse_edge_add = 1;
  if(threshold_nsigmas>0) threshold_nsigmas *= -1;
  if(threshold_nsigmas>-2) threshold_nsigmas = -2;
  return 0;
}

int SpeFinder::Finalize()
{
  return 0;
}

int SpeFinder::Process(ChannelData* chdata)
{
  //need a good baseline
  if(!chdata->baseline.found_baseline) {
    return 0;
  }

  int search_start_index  = chdata->TimeToSample(search_start_time);
  int search_end_index    = chdata->TimeToSample(search_end_time);
  if(search_end_index<search_start_index + pulse_edge_add*2) 
    return 0;

  int min_baseline_nsamps = min_baseline_length*chdata->sample_rate;
  if(min_baseline_nsamps<1) min_baseline_nsamps = 1;
  //if(min_baseline_nsamps<0) min_baseline_nsamps = 0;

  //use the subtracted waveform to search for spe start and end
  double* subtracted = chdata->GetBaselineSubtractedWaveform();
  double check_val = chdata->baseline.sigma*threshold_nsigmas;
  std::vector<int> start_index;
  std::vector<int> end_index;

  int first_baseline_samp = -1, last_baseline_samp = -1;
  bool in_pulse = false;
  for(int index = search_start_index+pulse_edge_add; 
      index<search_end_index-pulse_edge_add; index++){
    if(subtracted[index]<=check_val){
      if(subtracted[index-1]>check_val){//if just come to a pulse
        in_pulse = true;
	start_index.push_back(index-pulse_edge_add );
      }
      if(subtracted[index+1]>check_val){//if just leave a pulse
        if(in_pulse) end_index.push_back(index+pulse_edge_add );
        in_pulse = false;
      }
    }//end if subtracted index <=
    else{//if not in a pulse
      if(first_baseline_samp<0) first_baseline_samp = index;
      last_baseline_samp = index;
    }
  }
  if(in_pulse) start_index.pop_back();
  if(start_index.size()!=end_index.size() || start_index.size()==0
     || first_baseline_samp<0 || last_baseline_samp<0){
    start_index.clear();
    end_index.clear();
    return 0;
  }

  // Combine nearby pulses
  // Combine pulses that are not separated by more than min_baseline_nsamps
  for (size_t i=1; i<start_index.size(); ) {
    if (start_index.at(i)-end_index.at(i-1) < min_baseline_nsamps) {
      end_index.at(i-1) = end_index.at(i);      // Set previous pulse end to current pulse end
      start_index.erase(start_index.begin()+i); // Delete current pulse
      end_index.erase(end_index.begin()+i);
    } else {
      i++;
    }
  }



  //evaluate all the pulses, use raw waveform
  double* wave = chdata->GetWaveform();
  int pre_pulse_end = -1, post_pulse_start = -1, peak_index=-1;
  double bl_mean = 0, bl_sigma = 0, spe_integral=0, peak_amplitude=-1;
  for(size_t i=0; i<start_index.size(); i++){
    //    std::cout<<"Possible SPE pulse at: "<<chdata->SampleToTime(start_index.at(i))
    //     	     <<'\t'<<chdata->SampleToTime(end_index.at(i))<<std::endl;
    if(i>0) pre_pulse_end = end_index.at(i-1);
    else pre_pulse_end = first_baseline_samp;
    if(i<start_index.size()-1) post_pulse_start = start_index.at(i+1);
    else post_pulse_start = last_baseline_samp;
    //check the good baseline length
    if(start_index.at(i)-pre_pulse_end<min_baseline_nsamps 
       || post_pulse_start-end_index.at(i)<min_baseline_nsamps){
      //      std::cout<<"\tThe baseline is too short, give up"<<std::endl;
      //      std::cout<<start_index.at(i)<<" "<<pre_pulse_end<<" "<<post_pulse_start<<" "<<end_index.at(i)<<" "<<min_baseline_nsamps<<std::endl;
      continue;
    }
    //calculate the local baseline mean and sigma
    bl_mean = 0;
    bl_sigma = 0;
    spe_integral = 0;
    peak_amplitude=chdata->GetVerticalRange();
    for(int index=start_index.at(i)-min_baseline_nsamps;
	index<=end_index.at(i)+min_baseline_nsamps; index++){
      if(index>=start_index.at(i) && index<=end_index.at(i)){
	spe_integral += wave[index];
	if(wave[index]<peak_amplitude){
	  peak_amplitude = wave[index];
	  peak_index = index;
	}
      }
      else{
	bl_mean += wave[index];
	bl_sigma += wave[index]*wave[index];
      }//end else
    }//end for index
    bl_mean /= 2.*min_baseline_nsamps;
    bl_sigma /= 2.*min_baseline_nsamps;
    bl_sigma = std::sqrt(bl_sigma-bl_mean*bl_mean);
    spe_integral = bl_mean*(end_index.at(i)-start_index.at(i)+1) - spe_integral; 
    peak_amplitude = bl_mean - peak_amplitude;
    //maybe leave this cut to the root file?
    // if(bl_sigma>=1.2*chdata->baseline.sigma 
    //    || std::abs(bl_mean)*(end_index.at(i)-start_index.at(i)+1)>spe_integral*0.05)
    //   continue;

    //store spe
    Spe a_spe;
    a_spe.integral = spe_integral;
    a_spe.start_time = chdata->SampleToTime(start_index.at(i));
    a_spe.peak_time = chdata->SampleToTime(peak_index);
    a_spe.end_time = chdata->SampleToTime(end_index.at(i));
    a_spe.width = end_index.at(i)-start_index.at(i)+1;
    a_spe.amplitude = peak_amplitude;
    a_spe.amplitude_sigma = peak_amplitude/chdata->baseline.sigma;
    a_spe.start_index = start_index.at(i);
    a_spe.peak_index = peak_index;
    a_spe.end_index = end_index.at(i);
    a_spe.baseline_mean = bl_mean;
    a_spe.baseline_sigma = bl_sigma;
    chdata->single_pe.push_back(a_spe);
    
    // std::cout<<"Found SPE: "<<a_spe.integral<<'\t'<<a_spe.start_time
    // 	     <<'\t'<<a_spe.width<<'\t'<<a_spe.amplitude<<'\t'<<a_spe.peak_time
    // 	     <<'\t'<<a_spe.baseline_mean<<'\t'<<a_spe.baseline_sigma<<std::endl;
    
  }//end for i

  return 0;
}
