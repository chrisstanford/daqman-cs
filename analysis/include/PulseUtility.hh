#ifndef PULSE_UTILITY_h
#define PULSE_UTILITY_h

//this files provides a function to evaluate pulse objects
//it is used by Rois, Pulses, and TOFs -- Jingke Xu

#include <vector>

//class Pulse;
//class ChannelData;

//this function can not be part of the pulse class because pulse belongs to chdata
//it can not be part of pulse finder because EvalRoi calls it too
//function to calculate pulse parameters, if max_peak_time>0
//the peak of the pulse should be within (start time, start time + max_peak_time)
//int EvaluatePulse(Pulse& pulse, ChannelData* chdata, double max_peak_time=-1);
//determine if x1->x2 cross the relative threshold, either positive or negative
bool RelativeThresholdCrossed(double x1, double x2, double threshold);
//determine either x1 or x2 has smaller amplitude in the direction given by sign
bool FirstAmplitudeIsSmaller(double x1, double x2, double dumb_sign=1);
//calculate the gaussian function value for convolutions, etc.
double gaussian(double x, double mean, double sigma);
//smooth a waveform with a Gaussian function
bool SmoothWaveform(std::vector<double> &smoothed, int nsamps, const double *wave, int sigma);
//smooth a waveform with a simple running sum
bool RunningSumWaveform(std::vector<double> &smoothed, int nsamps, const double *wave,  int sigma);
//function to search for pulses/spikes with a relative threshold
//the waveform can be smoothed/downsampled (pulse finding) or not (spike finding)
int RelativeThresholdSearch(std::vector<double> wave, double start_threshold, double end_threshold,
			    std::vector<int> & start_index, std::vector<int> & end_index, 
			    int pulse_edge_add=3, int step_size=1, int search_start=0, int search_end=-1);

#endif
