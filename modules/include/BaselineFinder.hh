/** @file BaselineFinder.hh
    @brief Defines BaselineFinder module.
    @author rsaldanha, jingkexu
    @ingroup modules
*/

#ifndef BASELINEFINDER_h
#define BASELINEFINDER_h

#include "ChannelModule.hh"

#include <vector>
/** @class BaselineFinder
    @brief searches the beginning of a channel's waveform to determine baseline
    @ingroup modules
    */
class BaselineFinder : public ChannelModule
{
public:
  BaselineFinder();
  ~BaselineFinder();
  
  int Initialize();
  int Finalize();
  int Process(ChannelData* chdata);
  
  static const std::string GetDefaultName(){ return "BaselineFinder";}

  //parameters
  bool flat_baseline;         /// use fixed baseline algorithm
  bool stanford_baseline;     /// use Chris Stanford's baseline algorithm
  double pulse_start_time;    /// time when the pulses are expected to arrive, baseline should end
  int min_valid_nsamps;       /// minimum valid sample number in the pre trigger window
  int pulse_edge_add;         /// sample number to be added to the edge of a possible pulse
  int min_baseline;           /// The minimum ADC sample that can be treated as a valid baseline sample
  int max_baseline;           /// The maximum ADC sample that can be treated as a valid baseline sample

  int baseline_spread;        /// rought estimate of the baseline spread
  int group_adc_window;       /// number of adc counts to be grouped in flat baseline search

  int moving_window_nsamps;   /// sample number to be used in the drifting baseline search
  int pulse_start_inc;        /// Minimum ADC counts increase over 1-3 samples to start a pulse
  int pulse_end_inc;          /// Minimum ADC counts decrease over 1-3 samples to end a pulse
  int max_flat_nsamps;        /// maximum sample number to stay "flat" inside a pulse
  double min_good_fraction;   /// Minimum fraction of valid samples in a moving window

  double baseline_start_nsamps;  /// Number of samples to base the initial baseline RMS off of
  double start_RMS_factor;    /// If a sample is start_RMS_factor*baseline_RMS away from the baseline_mean, then start the pulse
  double end_RMS_factor;      /// If a end_nsamps_within_RMS consecutive samples are within maseline_mean +/- end_RMS_factor*baseline_RMS, then end the pulse
  double end_nsamps_within_RMS;// If a end_nsamps_within_RMS consecutive samples are within maseline_mean +/- end_RMS_factor*baseline_RMS, then end the pulse
  double min_nsamps_between_pulses;      /// Minimum number of samples between two pulses to consider them as separate

private:
  int FlatBaseline(ChannelData* chdata);
  int DriftingBaseline(ChannelData* chdata);
  int StanfordBaseline(ChannelData* chdata);
};

#endif
