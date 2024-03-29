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

BaselineFinder::BaselineFinder():
  ChannelModule(GetDefaultName(), "Find the baseline (zero) of the channel in the samples read before the trigger")
{
  AddDependency<ConvertData>();
  
  //Register all the config handler parameters
  RegisterParameter("flat_baseline", flat_baseline = false,
		    "Search for a flat baseline in pre-trigger window, otherwise search for a drifting baseline");
  RegisterParameter("stanford_baseline", stanford_baseline = false,
		    "Use Chris Stanford's Baseline finder algorithm");
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
  RegisterParameter("baseline_spread",  baseline_spread = 3,
                    "Rough estimate of baseline spread + uncertainty");
  RegisterParameter("group_adc_window", group_adc_window = 3,
                    "Window Length to be grouped for finding the most frequent adc count");

  //parameters for drifting baseline
  RegisterParameter("moving_window_nsamps",  moving_window_nsamps = 50,
                    "Number of samples in the moving baseline window");
  RegisterParameter("pulse_start_inc", pulse_start_inc = 3,
                    "Minimum ADC counts increase over 1-3 samples to start a pulse");
  RegisterParameter("pulse_end_inc", pulse_end_inc = 2,
                    "Minimum ADC counts decrease over 1-3 samples to end a pulse");
  RegisterParameter("max_flat_nsamps", max_flat_nsamps = 6,
                    "Maximum number of counts staying flat within a pulse");
  RegisterParameter("min_good_fraction", min_good_fraction = 0.5,
		    "Minimum fraction of valid samples in a moving window");
 

  //parameters for Chris Stanford's baseline
  RegisterParameter("baseline_start_nsamps",  baseline_start_nsamps = 100,
                    "Number of samples to base the initial baseline RMS off of");
  RegisterParameter("start_RMS_factor",  start_RMS_factor = 6,
                    "If a sample is start_RMS_factor*baseline_RMS away from the baseline_mean, then start the pulse");
  RegisterParameter("end_RMS_factor",  end_RMS_factor = 5,
                    "If end_nsamps_within_RMS consecutive samples are within maseline_mean +/- end_RMS_factor*baseline_RMS, then end the pulse");
  RegisterParameter("end_nsamps_within_RMS",  end_nsamps_within_RMS = 2*moving_window_nsamps,
                    "If end_nsamps_within_RMS consecutive samples are within maseline_mean +/- end_RMS_factor*baseline_RMS, then end the pulse");
  RegisterParameter("min_nsamps_between_pulses",  min_nsamps_between_pulses = 3*moving_window_nsamps,
                    "Minimum number of samples between two pulses to consider them as separate");

}

BaselineFinder::~BaselineFinder()
{
  Finalize();
}

int BaselineFinder::Initialize()
{
  // if(min_valid_nsamps<10)    min_valid_nsamps = 10;
  // if(pulse_edge_add<1)       pulse_edge_add = 1;
  //  if(baseline_spread<3)      baseline_spread = 3;
  //  if(group_adc_window<3)     group_adc_window = 3;
  if(moving_window_nsamps<1) moving_window_nsamps = 1;
  // if(pulse_start_inc<3)      pulse_start_inc = 3;
  // if(pulse_end_inc<3)        pulse_end_inc = 3;
  if(max_flat_nsamps<1)      max_flat_nsamps = 1;
  //  if(min_good_fraction<0.1)  min_good_fraction = 0.1;
  return 0;
}

int BaselineFinder::Finalize()
{   
  return 0;
}

int BaselineFinder::Process(ChannelData* chdata)
{
  int run_result=-1;
  if(stanford_baseline) run_result=StanfordBaseline(chdata);
  if(run_result && !flat_baseline) run_result=DriftingBaseline(chdata);
  if(run_result) run_result=FlatBaseline(chdata);
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
  //pulse start=1, end=-1, saturation=-2, middle of pulse=2
  for(int i=0; i<nsamps; i++){
    if(wave[i]>max_baseline || wave[i]<min_baseline || //PMT/electronics saturation
       std::fabs(chdata->GetVerticalRange()-wave[i])<1 || wave[i]<1) //ADC saturation
      mask[i] = -2; 
    if(i<1) continue;
    if(wave[i]-wave[i-1]>pulse_end_inc) mask[i]=-1;
    else if(wave[i-1]-wave[i]>pulse_start_inc) mask[i]=1;
    if(i<2) continue;
    if(wave[i]-wave[i-2]>pulse_end_inc) mask[i]=-1;
    else if(wave[i-2]-wave[i]>pulse_start_inc) mask[i]=1;
    if(i<3) continue;
    if(wave[i]-wave[i-3]>pulse_end_inc) mask[i]=-1;
    else if(wave[i-3]-wave[i]>pulse_start_inc) mask[i]=1;
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

  //warning: this is only to resolve the channel 4 problem for Radon runs July 2015
  //this is a stupid fix, and should not be used for other data, assuming baseline = 3959
  //  int baseline_est = 3958.5;//this is for runs before Aug 2015
  int baseline_est = 3996;//this is for runs before Aug 2015
  int mask_begin = -1, mask_end = -1;
  double total_sum=0, total_n=0;
  if(chdata->channel_id==5){
    mask_begin = chdata->TimeToSample(0.15);
    mask_end = chdata->TimeToSample(0.25);
    double aver = 0;
    for(int samp=mask_begin; samp<mask_end; samp++) aver += wave[samp];
    if(mask_end>mask_begin){
      aver /= mask_end-mask_begin;
      //      std::cout<<"The average at ~0.2us is: "<<aver<<std::endl;
    }
    mask_begin = chdata->TimeToSample(-0.1);
    if(aver<baseline_est-3.5) mask_end = chdata->TimeToSample(5.+0.035*(baseline_est-aver));
    for(int samp=mask_begin; samp<mask_end; samp++)
      mask[samp]=1;
  }

  //calculate the baseline
  int sum_nsamps=0, sum_samp=0, center_samp=0;
  double sum=0;
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

  stats.resize(pulse_start_index);
  for(int samp=0; samp<pulse_start_index; samp++){
    if(std::abs(wave[samp]-rough_bl)>baseline_spread)
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

  if(std::abs(max_bin-mx_bin)<=baseline_spread){
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
  else return 1;

  //subtract off the baseline
  for(int samp=0; samp<nsamps; samp++){
    baseform[samp] = wave[samp]-baseline.mean;
  }

  return 0;

}


double windowMean(double* wave, int start_index, int lookback, int nsamps) {
  int n=0;
  double v=0.;
  for (int i=start_index-lookback+1; i<=start_index; i++) {
    if (i<0 || i>=nsamps) continue;
    v+=wave[i];
    n++;
  }
  return n==0 ? wave[0] : v/n;
}

double windowMean(std::vector<double> wave, size_t start_index, int lookback) {
  int n=0;
  double v=0.;
  for (size_t i=start_index-lookback+1; i<=start_index; i++) {
    if (i<0 || i>=wave.size()) continue;
    v+=wave[i];
    n++;
  }
  return n==0 ? wave[0] : v/n;
}


double windowRMS(double* wave, int start_index, int lookback, int nsamps) {
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

double withinRange(double* wave, int start_index, int lookback, double mean, double range, int nsamps) {
  for (int i=start_index-lookback+1; i<=start_index; i++) {
    if (i<0 || i>=nsamps) continue;
    if (fabs(wave[i]-mean) > range) return false;
  }
  return true;  
}

int BaselineFinder::StanfordBaseline(ChannelData* chdata)
{
  bool debug=0;
  // Method: Identify pulses by getting a baseline RMS and a running baseline mean defining the start of the pulse to be when the pulse goes outside of mean +/- start_RMS_factor*RMS. Then the end of the pulse is when end_samples_within_RMS consecutive samples are within the mean +/- RMS from before the pulse began.
  if (debug) std::cout<<"1"<<std::endl;

  Baseline & baseline = chdata->baseline;
  double* wave = chdata->GetWaveform();
  const int nsamps = chdata->nsamps;
  if(nsamps<1) return -1;

  std::vector<double>& baseform = chdata->subtracted_waveform;
  baseform.resize(nsamps);

  // Initialize variables
  double total_baseline_sum = 0.;
  double total_baseline_rms = 0.;
  int total_baseline_nsamps = 0;
  bool in_pulse = false;
  double pre_pulse_window_mean = 0.; // The mean taken just before the current pulse
  //  double pre_pulse_window_RMS = 0.; // The RMS taken just before the current pulse
  double baseline_start_RMS = windowRMS(wave,baseline_start_nsamps,baseline_start_nsamps,nsamps);
  double baseline_start_mean = windowMean(wave,20,20,nsamps);
  for (int samp=0; samp<30; samp++) {
    // Is the a pulse at the start of the baseline where the RMS is calculated?
    if (fabs(wave[samp]-baseline_start_mean) > start_RMS_factor*baseline_start_RMS)
      return 0; // without setting found_baseline = true
  }

  if (debug) std::cout<<"2"<<std::endl;

  //  if (debug) std::cout<<"baseline_start_RMS: "<<baseline_start_RMS<<std::endl;
  std::vector<int> pulse_start;
  std::vector<int> pulse_end;
  // Loop over samples
  //  for(int samp=0;samp<moving_window_nsamps/**/;samp++) baseform[samp] = wave[samp];
  for(int samp=0; samp<nsamps; samp++){
    double preMean = windowMean(wave,samp-1,moving_window_nsamps,nsamps);
    double thisMean = windowMean(wave,samp,moving_window_nsamps,nsamps);
    // Identify start of pulse
    if (in_pulse == false && fabs(wave[samp]-preMean) > start_RMS_factor*baseline_start_RMS) {
      //      if (debug) std::cout<<"Start pulse at "<<chdata->SampleToTime(samp)<<std::endl;
      //      if (debug) std::cout<<"preMean "<<preMean<<std::endl;
      //      if (debug) std::cout<<"preRMS "<<preRMS<<std::endl;
      
      pulse_start.push_back(samp-pulse_edge_add);
      pre_pulse_window_mean = preMean;
      //      pre_pulse_window_RMS = preRMS;
      in_pulse = true;
    }
    // Identify end of pulse
    //    if (in_pulse == true && withinRange(wave,samp,end_nsamps_within_RMS,pre_pulse_window_mean,end_RMS_factor*baseline_start_RMS,nsamps)) {
    if (in_pulse == true && withinRange(wave,samp,end_nsamps_within_RMS,thisMean,end_RMS_factor*baseline_start_RMS,nsamps)) { // NEW Method
      // std::cout<<windowMean(wave,samp,moving_window_nsamps/2,nsamps)<<std::endl;
      // std::cout<<windowMean(wave,samp-moving_window_nsamps/2,moving_window_nsamps/2,nsamps)<<std::endl;
      // std::cout<<pulse_start.back()+pulse_edge_add<<" "<<chdata->SampleToTime(pulse_start.back()+pulse_edge_add)<<std::endl;
      // std::cout<<samp-end_nsamps_within_RMS<<" "<<chdata->SampleToTime(samp-end_nsamps_within_RMS)<<std::endl;
      // std::cout<<samp-end_nsamps_within_RMS+(samp-end_nsamps_within_RMS-(pulse_start.back()+pulse_edge_add))<<" "<<chdata->SampleToTime(samp-end_nsamps_within_RMS+(samp-end_nsamps_within_RMS-(pulse_start.back()+pulse_edge_add)))<<std::endl<<std::endl;
      //      pulse_end.push_back(samp-end_nsamps_within_RMS+3*(samp-end_nsamps_within_RMS-(pulse_start.back()+pulse_edge_add)));
      pulse_end.push_back(samp);
      in_pulse = false;
      //      if (debug) std::cout<<"Pulse end at "<<chdata->SampleToTime(samp)<<std::endl;
    }
    baseform[samp] = thisMean;
    if (in_pulse == false) {
      total_baseline_sum += wave[samp];
      total_baseline_rms += (wave[samp]-preMean)*(wave[samp]-preMean);
      total_baseline_nsamps += 1;
    }
    // Debug
    int start_print =  chdata->TimeToSample(11.9);
    int end_print =  chdata->TimeToSample(12.4);
    if (samp>start_print && samp<end_print) {
      //      if (debug) std::cout<<"Sample: "<<samp<<" Time: "<<chdata->SampleToTime(samp)<<" wave: "<<wave[samp]<<" windowMean: "<<windowMean(wave,samp,moving_window_nsamps/**/)<<" in_pulse: "<<in_pulse<<std::endl;
    }
  }
  if (in_pulse) {
    pulse_end.push_back(nsamps-1);
    in_pulse = false;
  }
  if (pulse_start.size() != pulse_end.size()) // Something went terribly wrong
    return 0; // without setting found_baseline = true

  if (debug) std::cout<<"3"<<std::endl;

  // Combine pulses that are not separated by more than min_valid_nsamps
  int n_pulses = pulse_start.size();
  for (int p=1; p<n_pulses; ) {
    if (pulse_start.at(p)-pulse_end.at(p-1) < min_nsamps_between_pulses) {
      n_pulses--;
      pulse_end.at(p-1) = pulse_end.at(p); // Set previous pulse end to current pulse end
      pulse_start.erase(pulse_start.begin()+p); // Delete current pulse
      pulse_end.erase(pulse_end.begin()+p);
    } else {
      p++;
    }
  }
  if (debug) std::cout<<"4"<<std::endl;
  // Loop through pulses, setting the baseline in the pulse region to be a straight line from the ADC value at the start of the pulse to the ADC value at the end of the pulse
  for (size_t p=0; p<pulse_start.size(); p++) {
    int start_index = pulse_start.at(p);
    int end_index = pulse_end.at(p);
    int n_index = end_index-start_index;
    double start_baseline_adc = baseform[start_index];
    double end_baseline_adc = baseform[end_index];
    for (int i=start_index; i<=end_index; i++) {
      if (i<0 || i>=nsamps) continue;
      // Linear interpolation
      double slope = (end_baseline_adc-start_baseline_adc)/n_index;
      baseform[i]=start_baseline_adc+slope*(i-start_index);
    }
  } 

  if (debug) std::cout<<"5"<<std::endl;
  // // Smooth the baseline twice
  // double* rough_bl = &baseform[0];
  // double* smooth_bl = new double[nsamps];
  // for(int i=0; i<nsamps; i++) {
  //   // Set a sample's value to the mean of 2*moving_window_nsamps samples with this sample as its center.
  //   smooth_bl[i] = 0;
  //   smooth_bl[i] = windowMean(rough_bl, i+moving_window_nsamps, 2*moving_window_nsamps, nsamps);
  // }
  // if (debug) std::cout<<"6"<<std::endl;
  // for(int i=0; i<nsamps; i++) {
  //   // Set a sample's value to the mean of moving_window_nsamps samples with this sample as its center.
  //   baseform[i] = windowMean(smooth_bl, i+moving_window_nsamps, 2*moving_window_nsamps, nsamps);
  // }
  // if (debug) std::cout<<"7"<<std::endl;
  // delete smooth_bl;

  //subtract off the baseline
  for(int samp=0; samp<nsamps; samp++){
    //    if (samp%1000==0) if (debug) std::cout<<wave[samp]<<" "<<baseform[samp]<<std::endl;
    baseform[samp] = wave[samp]-baseform[samp];
  }
  if (debug) std::cout<<"8"<<std::endl;

  // Set parameters
  baseline.found_baseline=true;
  baseline.mean = total_baseline_sum/total_baseline_nsamps;
  baseline.sigma = sqrt(total_baseline_rms/total_baseline_nsamps);
  baseline.length = total_baseline_nsamps;
  
  return 0;
}
