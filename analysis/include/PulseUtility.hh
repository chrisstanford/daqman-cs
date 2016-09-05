#ifndef PULSE_UTILITY_h
#define PULSE_UTILITY_h

#include <vector>

//determine if x1->x2 cross the relative threshold, either positive or negative
bool RelativeThresholdCrossed(double x1, double x2, double threshold);
//determine either x1 or x2 has smaller amplitude in the direction given by sign
bool FirstAmplitudeIsSmaller(double x1, double x2, double dumb_sign=1, double precision=1e-4);
//calculate the gaussian function value for convolutions, etc.
// double gaussian(double x, double mean, double sigma);
//smooth a waveform with a Gaussian function
bool SmoothWaveform(std::vector<double> &smoothed, int nsamps, const double *wave, int sigma);
//smooth a waveform with a simple running sum
bool RunningSumWaveform(std::vector<double> &summed, int nsamps, const double *wave,  int sigma);
//calculate the running sum waveform from waveform integral
bool RunningSumFromIntegral(std::vector<double> &summed, int nsamps, const double *integral, int sigma);
//add offset to a variable with upper and lower bounds
template <typename T> 
inline void AddOffsetWithBounds(T &val, T offset, T lower, T upper){
  val += offset;
  if(val>upper) val=upper;
  if(val<lower) val=lower;
}

#endif
