#include "Event.hh"
#include "BaselineFinder.hh"
#include "ConvertData.hh"
#include "intarray.hh"
#include "RootWriter.hh"
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
//#include <functional>
#include <algorithm>

bool IsThresholdCrossed(double x1, double x2, double threshold){
  if(threshold>0) return ( (x2-x1)>threshold );
  else if(threshold<0) return ( (x2-x1)<threshold );
  else {
    std::cerr<<"*** Error *** Comparing with 0 threshold *** Error ***"<<std::endl;
    return false;
  }
}

BaselineFinder::BaselineFinder():
  ChannelModule(GetDefaultName(), "Find the baseline (zero) of the channel in the samples read before the trigger")
{
  AddDependency<ConvertData>();
  
  //Register all the config handler parameters
  RegisterParameter("flat_channels",  flat_channels, "channels that use fixed baseline algorithm");
  RegisterParameter("relative_spe_threshold", relative_spe_threshold = false,
		    "whether the threshold values are relative to channel spe mean");
  RegisterParameter("suppress_zeros", suppress_zeros = false,
		    "whether we want to suppress zeros in baseline subtracted waveform");
  //   RegisterParameter("flat_baseline", flat_baseline = false,
  // 		    "Search for a flat baseline in pre-trigger window, otherwise search for a drifting baseline");
  RegisterParameter("pulse_start_time", pulse_start_time = -0.1,
                    "Time when the pulses are expected to arrive, baseline should end ");
  RegisterParameter("min_valid_nsamps", min_valid_nsamps = 10,
                    "Minimum samples for a baseline to be valid");
  RegisterParameter("pulse_edge_add", pulse_edge_add  = 2,
                    "Number of samples to be added at the edge of possible pulses");
  RegisterParameter("min_baseline", min_baseline = 3000,
                    "The minimum ADC sample that can be treated as a valid baseline sample");
  RegisterParameter("max_baseline", max_baseline = 4095,
                    "The maximum ADC sample that can be treated as a valid baseline sample");
  // RegisterParameter("",  = 10,
  //                   "");

  //parameters for fixed baseline
  RegisterParameter("pulse_threshold",  pulse_threshold = 3,
                    "absolute pulse threshold value, rough estimate of the baseline spread");
  RegisterParameter("group_adc_window", group_adc_window = 3,
                    "Window Length to be grouped for finding the most frequent adc count");

  //parameters for drifting baseline
  RegisterParameter("moving_window_nsamps",  moving_window_nsamps = 50,
                    "Number of samples in the moving baseline window");
  RegisterParameter("pulse_start_threshold", pulse_start_threshold = -3,
                    "Minimum ADC counts change over 1-3 samples to start a pulse");
  RegisterParameter("pulse_end_threshold", pulse_end_threshold = 2,
                    "Minimum ADC counts change over 1-3 samples to end a pulse");
  RegisterParameter("max_flat_nsamps", max_flat_nsamps = 6,
                    "Maximum number of counts staying flat within a pulse");
  RegisterParameter("min_good_fraction", min_good_fraction = 0.5,
		    "Minimum fraction of valid samples in a moving window");
  RegisterParameter("attenuated_ch", attenuated_ch = 1,
                    "Channel ID with attentuated input, will skip certain baseline regions");

}

BaselineFinder::~BaselineFinder()
{
  Finalize();
}

int BaselineFinder::Initialize()
{
  //we need to let the program fail if the user set unphysical parameters
  if(min_valid_nsamps<1 ||pulse_edge_add<1 || moving_window_nsamps<1 || 
     group_adc_window<1 || max_flat_nsamps<1 || min_good_fraction<=0 ||
     pulse_start_threshold==0 || pulse_end_threshold==0 || min_baseline>=max_baseline) 
    return -1;
  //we only use the absolute value of pulse_threshold for flat baseline search
  //so we allow users to use negative values if the want
  if(pulse_threshold<0)      pulse_threshold *= -1;
  return 0;
}

int BaselineFinder::Finalize()
{   
  return 0;
}

int BaselineFinder::Process(ChannelData* chdata)
{
  int run_result=-1;
  if( flat_channels.find( chdata->channel_id ) == flat_channels.end())
    run_result = DriftingBaseline(chdata);
  //we try the flat baseline if the drift baseline fails
  if(run_result) run_result = FlatBaseline(chdata);
  //we don't want to be overwhelmed by error messages due to baseline fail
  //dependent modules need to check whether baseline finder succeed or not
  return 0;
}

int BaselineFinder::DriftingBaseline(ChannelData* chdata)
{

  Baseline & baseline = chdata->baseline;
  double * wave = chdata->GetWaveform();
  const int nsamps = chdata->nsamps;
  if(nsamps<1) return -1;

  //find the maximum sample value within the pre_trigger area
  int pulse_start_index = chdata->TimeToSample(pulse_start_time) - 1;
  if(pulse_start_index<1) return -1;
  double max_pre_trig = *std::max_element(wave,wave+pulse_start_index);
  double min_pre_trig = *std::min_element(wave,wave+pulse_start_index);
  if(std::fabs(chdata->GetVerticalRange()-max_pre_trig )<0.1 || min_pre_trig<0.1)
    baseline.saturated = true;
  baseline.found_baseline = false;

  std::vector<double>& baseform = chdata->subtracted_waveform;
  baseform.resize(nsamps);
  std::vector<int> mask (nsamps, 0);

  //look for possible pulse starts and ends
  double start_threshold_adc, end_threshold_adc;
  if(relative_spe_threshold){
    if(chdata->spe_mean==1){
      Message(WARNING)<<"Channel "<<chdata->channel_id 
		      <<" attempts to use threshold values relative to uncalibrated spe_mean!\n";
    }
    start_threshold_adc = pulse_start_threshold*chdata->spe_mean;
    end_threshold_adc = pulse_end_threshold*chdata->spe_mean;
  }
  else{
    start_threshold_adc = pulse_start_threshold;
    end_threshold_adc = pulse_end_threshold;
  }

  //pulse start=1, end=-1, saturation=-2, middle of pulse=2
  for(int i=0; i<nsamps; i++){
    if(wave[i]>max_baseline || wave[i]<min_baseline || //PMT/electronics saturation
       std::fabs(chdata->GetVerticalRange()-wave[i])<1 || wave[i]<1) //ADC saturation
      mask[i] = -2; 
    if(i<1) continue;
    if(IsThresholdCrossed(wave[i-1],wave[i], end_threshold_adc)) mask[i]=-1;
    else if(IsThresholdCrossed(wave[i-1],wave[i], start_threshold_adc)) mask[i]=1;
//     if(wave[i]-wave[i-1]>pulse_end_threshold) mask[i]=-1;
//     else if(wave[i-1]-wave[i]>pulse_start_threshold) mask[i]=1;
    if(i<2) continue;
    if(IsThresholdCrossed(wave[i-2],wave[i], end_threshold_adc)) mask[i]=-1;
    else if(IsThresholdCrossed(wave[i-2],wave[i], start_threshold_adc)) mask[i]=1;
    //     if(wave[i]-wave[i-2]>pulse_end_threshold) mask[i]=-1;
    //     else if(wave[i-2]-wave[i]>pulse_start_threshold) mask[i]=1;
    if(i<3) continue;
    if(IsThresholdCrossed(wave[i-3],wave[i], end_threshold_adc)) mask[i]=-1;
    else if(IsThresholdCrossed(wave[i-3],wave[i], start_threshold_adc)) mask[i]=1;
    //     if(wave[i]-wave[i-3]>pulse_end_threshold) mask[i]=-1;
    //     else if(wave[i-3]-wave[i]>pulse_start_threshold) mask[i]=1;
  }//end for

  //combine pulses and remove false pulses
  int last_nonflat_index = -1;
  for(int i=0; i<nsamps; i++){
    //we are coming into a pulse
    if(mask[i]==1){
      //if there was a previous pulse rising
      if(last_nonflat_index>=0 && mask[last_nonflat_index]==1){
	if(i-last_nonflat_index<max_flat_nsamps)//if the previous rising is not far away
	  {int j=i; while(--j>last_nonflat_index) mask[j]=2;}
	else //the previous rising is fake, remove the mark
	  {int j=last_nonflat_index+1; while(--j>=0 && mask[j]==1) mask[j]=0;}
      }//end if last pulse index
      while(++i<nsamps && mask[i]==1);//move to the next non-rising sample
      last_nonflat_index = i-1;
    }//end mask==1
    //we are leaving a pulse
    else if(mask[i]==-1){
      //if the previous rising/falling/saturating is not far away
      if(i-last_nonflat_index<max_flat_nsamps)//allow last pulse index = -1
	{int j=i; while(--j>last_nonflat_index) mask[j]=2;}
      else{//this is a fake falling edge
	if(last_nonflat_index>=0 && mask[last_nonflat_index]==1)//the previous rising is fake too
	  {int j=last_nonflat_index+1;while(--j>=0 && mask[j]==1) mask[j]=0;}
	i--; while(++i<nsamps && mask[i]==-1) mask[i]=0;//get rid of the current fake falling
	continue;//do not update las_nonflat_index
      }//end else
      while(++i<nsamps && mask[i]==-1);
      last_nonflat_index = i-1;
    }//end else if mask = -1
    //here it saturates
    else if(mask[i]==-2){
      //if there was a previous pulse rising or saturation
      if(last_nonflat_index>=0 && (mask[last_nonflat_index]==1 || mask[last_nonflat_index]==-2)){
	if(i-last_nonflat_index<max_flat_nsamps)//if the previous rising/saturation is not far away
	  {int j=i; while(--j>last_nonflat_index) mask[j]=2;}
	else if(mask[last_nonflat_index]==1)//the previous rising is fake, remove the mark
	  {int j=last_nonflat_index+1; while(--j>=0 && mask[j]==1) mask[j]=0;}
      }//end if last pulse index
      while(++i<nsamps && mask[i]==-2);//move to the next non-rising sample
      last_nonflat_index = i-1;
    }//end mask==-2
  }//end combining pulses

  //tag edge of pulses, remove them from the bseline calculation
  int edge = 0;
  for(int samp=0; samp<nsamps; samp++){
    if(mask[samp]==-1) edge=pulse_edge_add+1;//tail may drop slowly
    else if(edge){ if(mask[samp]==0) mask[samp]=-1; edge--;}
  }
  edge = 0;
  for(int samp=nsamps-1; samp>=0; samp--){
    if(mask[samp]==1) edge=pulse_edge_add;
    else if(edge) {if(mask[samp]==0) mask[samp]=1; edge--;}
  }

  //get rid of the sparse baseline samples
  for(int samp=1; samp<nsamps-2; samp++){
    if(!mask[samp] && mask[samp-1] && (mask[samp+1]||mask[samp+2]))
      mask[samp]=1;//if the previous sample and the next 1 or 2 samples are not baseline
  }

  //warning: this is temporary
  double baseline_est = 3996;//this is for runs before Aug 2015
  int mask_begin = -1, mask_end = -1;
  double total_sum=0, total_n=0;
  if(chdata->channel_id==attenuated_ch){
    //estimate the baseline here
    for(int i=0; i<pulse_start_index; i++){
      if(!mask[i]){
	total_sum+=wave[i]; total_n++; 
      }
    }
    if(total_n>0) baseline_est=total_sum/total_n;
    else return -1;

    mask_begin = chdata->TimeToSample(0.15);
    mask_end = chdata->TimeToSample(0.25);
    double aver = 0;
    for(int samp=mask_begin; samp<mask_end; samp++) aver += wave[samp];
    if(mask_end>mask_begin){
      aver /= mask_end-mask_begin;
      //      std::cout<<"The average at ~0.2us is: "<<aver<<std::endl;
    }
    mask_begin = chdata->TimeToSample(-0.1);
    if(aver<baseline_est-2.) mask_end = chdata->TimeToSample(5.+0.035*(baseline_est-aver));
    for(int samp=mask_begin; samp<mask_end; samp++)
      mask[samp]=1;
  }

  //calculate the baseline
  int sum_nsamps=0, sum_samp=0, center_samp=0;
  double sum=0;
  total_sum=0, total_n=0;
  for(int i=0; i<nsamps; i++){
    baseform[i]=0;
    if(!mask[i]){sum+=wave[i];sum_samp+=i;sum_nsamps++; 
      if(mask_end>pulse_start_index && i<mask_end+pulse_start_index){ total_sum+=wave[i]; total_n++;}
    }
    if(i<moving_window_nsamps-1) continue;
    if(sum_nsamps>=moving_window_nsamps*min_good_fraction){
      center_samp = sum_samp/sum_nsamps;
      if(baseform[center_samp]==0 && std::abs(i-moving_window_nsamps/2-center_samp)<moving_window_nsamps*0.2)
	{baseform[center_samp]=sum/sum_nsamps;}
      //	baseform[i-moving_window_nsamps/2]=sum/sum_nsamps;
    }
    if(!mask[i-moving_window_nsamps+1])
      {sum-=wave[i-moving_window_nsamps+1]; sum_samp-=i-moving_window_nsamps+1;sum_nsamps--;}
  }//end for i

  //warning: this is only to resolve the channel 4 problem for Radon runs July 2015
  //this is a stupid fix, and should not be used for other data
  if(mask_begin>0 && mask_end>mask_begin && total_n>0){
    total_sum /= total_n;
    for(int samp=mask_begin; samp<mask_end; samp++) baseform[samp] = total_sum;
  }

  {  //smooth the baseline twice
    std::vector<double> smooth_bl (nsamps, 0);
    sum_nsamps=0; sum=0; sum_samp=0;
    for(int i=0; i<nsamps; i++){
      smooth_bl[i]=0;
      if(baseform[i]>0){sum+=baseform[i];sum_samp+=i;sum_nsamps++;}
      if(i<moving_window_nsamps-1) continue;
      if(sum_nsamps>=moving_window_nsamps*min_good_fraction){
	center_samp = sum_samp/sum_nsamps;
	if(smooth_bl[center_samp]==0 && std::abs(i-moving_window_nsamps/2-center_samp)<moving_window_nsamps*0.2)
	  {smooth_bl[center_samp]=sum/sum_nsamps;}
	//simply assume the middle point of the window
	//        smooth_bl[i-moving_window_nsamps/2]=sum/sum_nsamps;
      }
      if(baseform[i-moving_window_nsamps+1]>0)
        {sum-=baseform[i-moving_window_nsamps+1]; sum_samp-=i-moving_window_nsamps+1; sum_nsamps--;}
    }//end for i
    sum_nsamps=0; sum=0; sum_samp=0;
    for(int i=0; i<nsamps; i++){
      baseform[i]=0;
      if(smooth_bl[i]>0){sum+=smooth_bl[i];sum_samp+=i;sum_nsamps++;}
      if(i<moving_window_nsamps-1) continue;
      if(sum_nsamps>=moving_window_nsamps*min_good_fraction){
	center_samp = sum_samp/sum_nsamps;
	if(baseform[center_samp]==0 && std::abs(i-moving_window_nsamps/2-center_samp)<moving_window_nsamps*0.2)
	  {baseform[center_samp]=sum/sum_nsamps;}
	//simply assume the middle point of the window
	//        baseform[i-moving_window_nsamps/2]=sum/sum_nsamps;
      }
      if(smooth_bl[i-moving_window_nsamps+1]>0)
        {sum-=smooth_bl[i-moving_window_nsamps+1]; sum_samp-=i-moving_window_nsamps+1; sum_nsamps--;}
    }//end for i
  }

  //fill up the empty baseline
  int last_baseline_index=-1;
  for(int i=0; i<nsamps; i++){
    if(baseform[i]<1e-10) continue;
    int j=i;
    while(--j>last_baseline_index){
      if(last_baseline_index<0) baseform[j]=baseform[i];
      else baseform[j]=(baseform[i]*(j-last_baseline_index)+baseform[last_baseline_index]*(i-j))/(i-last_baseline_index);
    }//end while
    last_baseline_index=i;
  }//end for int i
  if(last_baseline_index<0) return -1;
  int j=nsamps;
  while(--j>last_baseline_index) baseform[j]=baseform[last_baseline_index];
  
  double sum2=0;
  sum=0, sum_nsamps=0;
  for(int i=0; i<pulse_start_index; i++){
    if(!mask[i]){
      sum+=wave[i]; sum2+=wave[i]*wave[i];sum_nsamps++; 
    }
  }

  if(sum_nsamps<min_valid_nsamps) return -1;
  baseline.found_baseline = true;
  baseline.mean = sum/sum_nsamps;
  baseline.sigma = std::sqrt(sum2/sum_nsamps-baseline.mean*baseline.mean);
  baseline.length = sum_nsamps;

  //subtract off the baseline
  for(int i=0; i<nsamps; i++){
    baseform[i] = wave[i]-baseform[i];
    if(suppress_zeros && !mask[i]) baseform[i] = 0;
  }
  return 0;

}

//Search for a flat baseline in the pre-trigger window
int BaselineFinder::FlatBaseline(ChannelData* chdata){

  Baseline & baseline = chdata->baseline;
  double * wave = chdata->GetWaveform();
  const int nsamps = chdata->nsamps;
  if(nsamps<1) return -1;

  //find the maximum sample value within the pre_trigger area
  int pulse_start_index = chdata->TimeToSample(pulse_start_time) - 1;
  if(pulse_start_index<1 || pulse_start_index>=nsamps-1) return -1;
  double max_pre_trig = *std::max_element(wave,wave+pulse_start_index);
  double min_pre_trig = *std::min_element(wave,wave+pulse_start_index);
  if(std::fabs(chdata->GetVerticalRange()-max_pre_trig )<1 || min_pre_trig<1)
    baseline.saturated = true;
  baseline.found_baseline = false;

  std::vector<double>& baseform = chdata->subtracted_waveform;
  baseform.resize(nsamps);

  int max_bl = int(max_pre_trig);
  int min_bl = int(min_pre_trig);
  if(max_bl>max_baseline) max_bl = max_baseline;
  if(min_bl<min_baseline) min_bl = min_baseline;
  if(min_bl<1) min_bl=1;//can not be 0 to avoid saturation pulses
  if(max_bl-min_bl<1) return -1;
  std::vector<int> stats (max_bl-min_bl+1, 0);
  for(int samp=0; samp<pulse_start_index; samp++){
    if(wave[samp]<min_bl || wave[samp]>max_bl) continue;
    stats.at(wave[samp]-min_bl) ++;
  }

  int moving_window=group_adc_window;
  if(moving_window*2.>stats.size()) moving_window=stats.size()/2;
  //max: maximum, mx: next maximum
  int max_sum=-1, mx_sum=-1, max_bin=-1, mx_bin=-1, moving_sum=0;
  for(int i=0; i<stats.size()+0.; i++){
    moving_sum += stats.at(i);
    if(i<moving_window-1) continue;
    if(moving_sum>=max_sum){//probably rising part (>=) update max and mx
      mx_sum = max_sum;
      mx_bin = max_bin;
      max_sum = moving_sum;
      max_bin = i;
    }
    else if(moving_sum>mx_sum){//probably falling part > instead of >=
      mx_sum = moving_sum;
      mx_bin = i;
    }
    moving_sum -= stats.at(i-moving_window+1);
  }
  //a quick baseline estimate is the most frequent adc count number
  double rough_bl = min_bl+max_bin-(moving_window-1)/2.;

  double threshold_adc;
  if(relative_spe_threshold){
    if(chdata->spe_mean==1){
      Message(WARNING)<<"Channel "<<chdata->channel_id 
		      <<" attempts to use threshold values relative to uncalibrated spe_mean!\n";
    }
    threshold_adc = pulse_threshold*chdata->spe_mean;
  }
  else{
    threshold_adc = pulse_threshold;
  }

  stats.resize(pulse_start_index);
  for(int samp=0; samp<pulse_start_index; samp++){
    if(std::abs(wave[samp]-rough_bl)>threshold_adc)
      stats[samp]=1;//1 means might be in pulse
    else stats[samp]=0;//0 means baseline
  }

  int edge = 0;
  for(int samp=0; samp<pulse_start_index; samp++){
    if(stats[samp]) edge=pulse_edge_add+1;
    else if(edge){stats[samp]=1; edge--;}
  }

  edge = 0;
  for(int samp=pulse_start_index-1; samp>=0; samp--){
    if(stats[samp]) edge=pulse_edge_add;
    else if(edge) {stats[samp]=1; edge--;}
  }

  double sum=0, sum2=0;
  int sum_nsamps=0;
  for(int samp=0; samp<pulse_start_index; samp++){
    if(stats[samp]==0){
      sum += wave[samp];
      sum2 += wave[samp]*wave[samp];
      sum_nsamps ++;
    }
  }

  if(std::abs(max_bin-mx_bin)<=threshold_adc){
    baseline.found_baseline = true;
    baseline.length = sum_nsamps;
    if(sum_nsamps>min_valid_nsamps){
      baseline.mean = sum/sum_nsamps;
      baseline.sigma = std::sqrt(sum2/sum_nsamps-baseline.mean*baseline.mean);
    }
    else{
      baseline.mean = rough_bl;
      baseline.sigma = 0;
    }
  }//if abs
  else return -1;

  //if we want to suppress zeros, we need to find the pulses and the edges
  if(suppress_zeros){
    //mark the pulses
    stats.resize(nsamps);
    for(int samp=0; samp<nsamps; samp++){
      if(std::abs(wave[samp]-baseline.mean)>threshold_adc)
	stats[samp]=1;//1 means might be in pulse
      else stats[samp]=0;//0 means baseline
    }
    //mark the edge of the pulses
    int edge = 0;
    for(int samp=0; samp<nsamps; samp++){
      if(stats[samp]) edge=pulse_edge_add+1;
      else if(edge){stats[samp]=1; edge--;}
    }
    edge = 0;
    for(int samp=nsamps-1; samp>=0; samp--){
      if(stats[samp]) edge=pulse_edge_add;
      else if(edge) {stats[samp]=1; edge--;}
    }
  }//end if suppress zeros

  //subtract off the baseline
  for(int samp=0; samp<nsamps; samp++){
    baseform[samp] = wave[samp]-baseline.mean;
    if(suppress_zeros && !stats[samp]) baseform[samp] = 0;
  }

  return 0;

}
