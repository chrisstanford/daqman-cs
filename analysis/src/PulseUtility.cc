#ifndef PULSE_UTILITY_c
#define PULSE_UTILITY_c

#include "ChannelData.hh"
#include "PulseUtility.hh"
#include <algorithm>
#include <numeric>
#include <cmath>

//determine if x1->x2 cross the relative threshold, either positive or negative
bool RelativeThresholdCrossed(double x1, double x2, double threshold){
  if(threshold>0) return ( (x2-x1)>threshold );
  else if(threshold<0) return ( (x2-x1)<threshold );
  else {
    std::cerr<<"*** Error *** Comparing with 0 threshold *** Error ***"<<std::endl;
    return false;
  }
}

//determine either x1 or x2 has smaller amplitude in the direction given by sign
bool FirstAmplitudeIsSmaller(double x1, double x2, double dumb_sign){
  if(dumb_sign>0) return (x1<x2);
  else if(dumb_sign<0) return (x1>x2);
  else {
    std::cerr<<"*** Error *** Comparing amplitudes without a given sign *** Error ***"<<std::endl;
    return false;
  }
}

//calculate the gaussian function value for convolutions, etc.
double gaussian(double x, double mean, double sigma){
  if(std::fabs(x-mean)>=10*sigma) return 0;
  double arg = (x - mean)/sigma;
  arg *= -0.5 * arg;
  return 0.3989422804/sigma*exp(arg);
}

//smooth a waveform with a running average
bool SmoothWaveform(std::vector<double> &smoothed, int nsamps, const double *wave,  int sigma){
  
  if(nsamps<1 || sigma <=1) return false;
  smoothed.resize(nsamps);

  double result=0;
  int index=0, samp=-1, nsummed=0;  
  
  //  for (samp=0; samp<nsamps; samp++){
  while (++samp>=0){//while true
    if(samp<nsamps) result += wave[samp];
    if(samp>=sigma){
      result -= wave[samp-sigma];
      index = samp-sigma/2;
      nsummed = sigma;
      if(samp>=nsamps) nsummed += nsamps-samp-1;
    }
    else{
      index = samp/2;
      nsummed = samp+1;
    }
    //    if(samp<=sigma) std::cout<<samp<<'\t'<<index<<std::endl;
    if(index>=nsamps) break;
    if(nsummed>0) smoothed.at(index) = result/nsummed;
    else smoothed.at(index) = 0;
  }//end for samp

  return true;
}

// //smooth a waveform with a Gaussian function
// bool SmoothWaveform(std::vector<double> &smoothed, int nsamps, const double *wave,  int sigma){
  
//   if(nsamps<1 || sigma <=1) return false;
//   smoothed.resize(nsamps);
//   int start, end;
//   double result=0;

//   for (int ii=0; ii<nsamps; ii++){
//     start = ii - 3*sigma;
//     end   = ii + 3*sigma;
//     if(start<0)       start= 0;
//     if(end>=nsamps) end  = nsamps-1;
//     result = 0.;
//     for(int jj=start; jj<=end; jj++){
//       result += wave[jj]*gaussian(jj, ii, sigma);
//     }//end for jj
//     smoothed.at(ii) = result;
//   }//end for ii

//   return true;
// }

//smooth a waveform with a simple running sum
bool RunningSumWaveform(std::vector<double> &smoothed, int nsamps, const double *wave,  int sigma){
  
  if(nsamps<1 || sigma <=1) return false;
  smoothed.resize(nsamps);
  double result=0;
  int index=0, samp=-1;  
  
  //  for (samp=0; samp<nsamps; samp++){
  while (++samp>=0){//while true
    if(samp<nsamps) result += wave[samp];
    if(samp>=sigma){
      result -= wave[samp-sigma];
      index = samp-sigma/2;
    }
    else index = samp/2;
    //    if(samp<=sigma) std::cout<<samp<<'\t'<<index<<std::endl;
    if(index>=nsamps) break;
    smoothed.at(index) = result;
  }//end for samp

  return true;
}
// bool RunningSumWaveform(std::vector<double> &smoothed, int nsamps, const double *wave,  int sigma){
  
//   if(nsamps<1 || sigma <=1) return false;
//   smoothed.resize(nsamps);
//   int start, end;
//   double result=0;
  
//   for (int ii=0; ii<nsamps; ii++){
//     start = ii - sigma/2;
//     end   = ii + sigma - sigma/2;
//     if(start<0)       start= 0;
//     if(end>=nsamps) end  = nsamps-1;
//     result = 0.;
//     for(int jj=start; jj<=end; jj++){
//       result += wave[jj];
//     }//end for jj
//     smoothed.at(ii) = result;
//   }//end for ii

//   return true;
// }

int RelativeThresholdSearch(std::vector<double> wave, double start_threshold, double end_threshold,
			    std::vector<int> & start_index, std::vector<int> & end_index, 
			    int pulse_edge_add, int step_size, int search_start, int search_end){
  std::cerr<<"Begining of RelativeThresholdSearch"<<std::endl;
  if(wave.size()<=1 || step_size<1 || start_threshold==0 || end_threshold==0) return -1;
  std::cerr<<"RelativeThresholdSearch"<<std::endl;
  const int nsamps = wave.size();
  if(search_start<0) search_start = 0;
  if(search_end<0 || search_end>=nsamps) search_end = nsamps-1;
  if(start_index.size()) start_index.clear();
  if(end_index.size()) end_index.clear();

  std::vector<int> mask(nsamps,0);
  //let's first mark the pulse starts and ends
  for(int i=search_start; i<=search_end; i++){
    if(i<1) continue;
    //the continue commands let pulse end overrule pulse start, and let nearby points overrule far away ones
    if(RelativeThresholdCrossed(wave[i-1],wave[i], end_threshold)) {mask[i]=-1; continue;}
    else if(RelativeThresholdCrossed(wave[i-1],wave[i], start_threshold)) {mask[i]=1; continue;}
    if(i<2*step_size) continue;
    if(RelativeThresholdCrossed(wave[i-2*step_size],wave[i], end_threshold)) {mask[i]=-1; continue;}
    else if(RelativeThresholdCrossed(wave[i-2*step_size],wave[i], start_threshold)) {mask[i]=1; continue;}
    if(i<3*step_size) continue;
    if(RelativeThresholdCrossed(wave[i-3*step_size],wave[i], end_threshold)) {mask[i]=-1; continue;}
    else if(RelativeThresholdCrossed(wave[i-3*step_size],wave[i], start_threshold)) {mask[i]=1; continue;}
  }//end for

  //combine pulses and remove false pulses
  //we are more agreesive here than in baseline finder because we don't want to miss pulses
  for(int i=0; i<nsamps; i++){
    if(mask[i]==1 && start_index.size()==end_index.size()){
      start_index.push_back(i);
      continue;
    }
    if(mask[i]==-1){
      if(start_index.size()==0){
	if(i>search_start){
	  start_index.push_back(search_start);
	  end_index.push_back(i);
	}
      }
      else if(start_index.size()==end_index.size()) end_index.back()=i;
      else if(start_index.size()==end_index.size()+1) end_index.push_back(i);
    }
  }//end for
  if(start_index.size()==end_index.size()+1){
    if(start_index.back()<search_end) end_index.push_back(search_end);
    else start_index.pop_back();
  }
  
  if(start_index.size()!=end_index.size() || start_index.size()==0){
    std::cerr<<"*** Error *** start and end index vector size: "<<start_index.size()<<'\t'<<end_index.size()<<std::endl;
    if(start_index.size()) start_index.clear();
    if(end_index.size()) end_index.clear();
    return -1;
  }

  //let's resolve possible overlaps
  for(size_t index=0; index<start_index.size(); index++){
    // std::cout<<chdata->channel_id<<'\t'<<index<<'\t'<<chdata->SampleToTime(start_index.at(index))
    // 	     <<'\t'<<chdata->SampleToTime(end_index.at(index))<<std::endl;
    if(start_index.at(index)>=search_start+pulse_edge_add) start_index.at(index)-=pulse_edge_add;
    else start_index.at(index)=search_start;
    if(end_index.at(index)+pulse_edge_add<=search_end) end_index.at(index) += pulse_edge_add;
    else end_index.at(index)=search_end;
    if(index==0) continue;

    if(start_index.at(index)<=end_index.at(index-1)){
      int middle_index = (start_index.at(index)+end_index.at(index-1))/2;
      double min_value = wave[middle_index];
      //      for(int j=start_index.at(index); j<=end_index.at(index-1); j++){
      for(int j=end_index.at(index-1)+1-pulse_edge_add; j<start_index.at(index)+pulse_edge_add; j++){
	if(j<0) continue;
	else if(j>=nsamps) break;
	if(FirstAmplitudeIsSmaller(min_value,wave[j], start_threshold)) continue;
        min_value = wave[j];
        middle_index = j;
      }//end for j loop
      start_index.at(index) = middle_index;
      end_index.at(index-1) = middle_index;
    }//end if statement
  }//end for index loop

  return 0;
}

#endif
