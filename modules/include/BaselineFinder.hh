/** @file BaselineFinder.hh
    @brief Defines BaselineFinder module.
    @author jingkexu
    @ingroup modules
*/

#ifndef BASELINEFINDER_h
#define BASELINEFINDER_h

#include "ChannelModule.hh"

#include <vector>
#include <set>
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
  double baseline_end_time;       /// time when the pulses are expected to arrive, baseline should end
  int min_baseline_adc;           /// The minimum ADC sample that can be treated as a valid baseline sample
  int max_baseline_adc;           /// The maximum ADC sample that can be treated as a valid baseline sample
  int min_valid_nsamps;           /// minimum valid sample number in the pre trigger window
  bool relative_spe_threshold;    /// whether the threshold values are relative to channel spe mean
  int pulse_edge_add_nsamps;      /// sample number to be added to the edge of a possible pulse
  bool suppress_zeros;            /// whether we want to suppress zeros in baseline subtracted waveform
  double suppress_nsigmas;        /// how many baseline sigmas do we want to suppress

  std::set<int> flat_channels;    /// channels that use fixed baseline algorithm
  double pulse_threshold;         /// absolute pulse threshold value, rough estimate of the baseline spread
  int adc_group_window;           /// number of adc counts to be grouped in flat baseline search

  double moving_window_length_us; /// length of the moving baseline window in us
  double pulse_start_threshold;   /// Minimum ADC counts change over 1-3 samples to start a pulse
  double pulse_end_threshold;     /// Minimum ADC counts change over 1-3 samples to end a pulse
  double max_flat_pulse_time;     /// maximum time in us for waveform to stay "flat" inside a pulse
  double min_good_fraction;       /// Minimum fraction of valid samples in a moving window

private:
  int FlatBaseline(ChannelData* chdata);
  int DriftingBaseline(ChannelData* chdata);

};

#endif
