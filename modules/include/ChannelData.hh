/** @file ChannelData.hh
    @brief Defines the ChannelData storage class
    @author bloer
    @ingroup modules
    @ingroup daqroot
*/

#ifndef CHANNELDATA_h
#define CHANNELDATA_h

#include "Rtypes.h" //has the classdef macro
#include <string>
#include <vector>
#include <utility>
#include <stdint.h>

//subclasses
#include "Pulse.hh"
#include "Baseline.hh"
#include "Spe.hh"

class TGraph;
//notice: members with comment starting with ! are not saved
/** @class ChannelData
    @brief processed information for each channel and pointers to raw data
    @ingroup modules
    @ingroup daqroot
*/
class ChannelData{
  //interface
public:
  ChannelData() { Clear(); }
  virtual ~ChannelData(){}
  void Clear(); ///< reset all members
  /// Get a TGraph for drawing.  Needs to be deleted
  TGraph* GetTGraph(bool baseline_subtracted = false, int downsample=1) const; 
  /// Draw the channel on the current canvas
  void Draw(bool baseline_subtracted = false, int downsample=1,
	    bool autoscalex=true, bool autoscaley=true, double xmin=0,
	    double xmax=0, double ymin=0, double ymax=0); 
  /// Print the channel information
  void Print (int verbosity);  
  /// Convert a time in us to sample index. If checkrange, make sure that it is a valid sample number
  int TimeToSample(double time, bool checkrange=true);
  /// Convert a sample index to time in us
  double SampleToTime(int sample){ return (sample - trigger_index)/sample_rate;}
  /// Get a pointer to the channel's waveform data, preserver const-ness
  const double* GetWaveform() const { return &(waveform[0]); }
  /// Get a pointer to the channel's waveform data non-constly
  double* GetWaveform() { return &(waveform[0]); }
  /// Get a pointer to the waveform data after baseline subtraction constly
  const double* GetBaselineSubtractedWaveform() const 
  { return &(subtracted_waveform[0]); }
  /// Get a pointer to the waveform data after baseline subtraction
  double* GetBaselineSubtractedWaveform()
  { return &(subtracted_waveform[0]); }
  /// Get a pointer to the waveform data after integration
  const double* GetIntegralWaveform() const 
  { return &(integral[0]); }
  /// Get a pointer to the waveform data after integration constly
  double* GetIntegralWaveform()
  { return &(integral[0]); }
  
  /// Get the maximum vertical range of the digitizer
  double GetVerticalRange() const { return ((uint64_t)1<<sample_bits) - 1; } 

public:
  /// Define some 'fake' channel types
  enum ID_TYPES { CH_INVALID=-1, CH_SUM=-2, CH_DEFAULT=-3 };
public:
  //data members
  //raw-ish information
  int board_id;  ///< id of the board the channel is on
  int board_num; //!< number board is on in this event 
  int channel_num; ///< number of the channel on the board
  int channel_id; ///< global unique channel id: 10*board_id + channel_num;
  std::string label; ///< helpful label for each channel
  uint64_t timestamp; ///< trigger time for this channel (ns since run start)
  int sample_bits; //!< digitizer resolution, needed for readout
  double sample_rate; ///< samples per microsecondf
  int trigger_index; ///< sample at which trigger occurred
  int nsamps; //!< number of samples in the waveform
  std::vector<double> generic; ///< vector for generic data for testing
  bool saturated; ///< Did the signal hit the max or min range of the digitizer?
  double maximum; ///< Maximum value obtained by channel over range
  double minimum; ///< Minimum value obtained by channel over range
  double max_time; ///< time in us at which signal first achieved max
  double min_time; ///< time in us at which signal first achieved min
  const char* channel_start; //!< pointer to start of waveform
  const char* channel_end;   //!< pointer to end of waveform
  double spe_mean;     ///< mean photoelectron response read from database
  //  double spe_sigma;
  //vector waveforms
  std::vector<double> waveform; //!< Raw waveform as a double array
  std::vector<std::pair<int,int> > unsuppressed_regions; //!< list of begin,end sample of non-zero-suppressed regions in the waveform
  
  //processed information
  //baseline finding
  Baseline baseline;  ///< information about the channel's baseline
  std::vector<double> subtracted_waveform; //!< Channel waveform after baseline subtraction
  //individual pulses
  int npulses; ///< number of pulses found
  std::vector<Pulse> pulses;  //!< vector of individual pulses found 
  std::vector<Pulse> regions; ///< vector of ROIs for simple analysis
  std::vector<Pulse> tof;     ///< vector of Time of flight
  std::vector<Spe> single_pe; ///< vector of single photoelectron responses

  //integrator
  std::vector<double> integral; //!< Integral of the channel's waveform
  double integral_max;          ///< maximum value of integral along entire trigger
  double integral_min;          ///< minimum value of integral along entire trigger
  int integral_max_index;       ///< sample index at which max integral occurs
  int integral_min_index;       ///< sample index at which min integral occurs
  double integral_max_time;     ///< time at which max integral occurs
  double integral_min_time;     ///< time at which min integral occurs
  
  ClassDef(ChannelData,17)
};

inline void ChannelData::Clear()
{
  board_id = -1;
  board_num = -1;
  channel_num = -1;
  channel_id = -1;
  label = "";
  timestamp = 0;
  sample_bits = -1;
  sample_rate = -1;
  trigger_index = -1;
  nsamps = -1;
  saturated = false;
  maximum=-1;
  minimum=-1;
  generic.clear();
  channel_start = NULL;
  channel_end = NULL;
  spe_mean = 1.;
  //  spe_sigma = 0.;
  baseline.Clear();
  subtracted_waveform.clear();
  npulses=0;
  pulses.clear();
  regions.clear();
  tof.clear();
  single_pe.clear();
  integral.clear();
  integral_max = -1;
  integral_min = -1;
  integral_max_index = -1;
  integral_min_index = -1;
  integral_max_time = -1;
  integral_min_time = -1;
}

inline int ChannelData::TimeToSample(double time, bool checkrange)
{ 
  int samp = (int)(time*sample_rate + trigger_index);
  if(checkrange){
    if(samp < 0) samp = 0;
    if(samp > nsamps-1) samp = nsamps-1;
  }
  return samp;
}

#endif
