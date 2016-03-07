/** @file PulseFinder.hh
    @brief Defines the PulseFinder module
    @author rsaldanha, bloer
    @ingroup modules
*/
#ifndef PULSEFINDER_h
#define PULSEFINDER_h

#include "BaseModule.hh"
#include <vector>

/** @class PulseFinder
    @brief Searches for individual scintillation events within a trigger
    @ingroup modules
*/
class PulseFinder : public BaseModule{
public:
  PulseFinder();
  ~PulseFinder();
   
  int Initialize();
  int Finalize();
  int Process(EventPtr evt);
  
  /// Search for pulses using a simple discrimination threshold
  void DiscriminatorSearch(ChannelData* chdata,
			   std::vector<int>& start_index, 
			   std::vector<int>& end_index);
  
  /// Search for pulses based on the measured variance of the baseline
  void VarianceSearch(ChannelData* chdata,
		      std::vector<int>& start_index, 
		      std::vector<int>& end_index);
  
  void IntegralSearch(ChannelData* chdata,
		      std::vector<int>& start_index,
		      std::vector<int>& end_index);
  
  /// Search for pulses based on the curvature of the integral
  void CurvatureSearch(ChannelData* chdata,
		       std::vector<int>& start_index,
		       std::vector<int>& end_index);
  
  static const std::string GetDefaultName(){ return "PulseFinder"; }
  
  enum SEARCH_MODE { VARIANCE , DISCRIMINATOR , INTEGRAL , CURVATURE };

private:
  ///parameters
  SEARCH_MODE mode;                ///< Which search function to use
  ///currently only used in discriminator and integral search
  double search_start_time;        ///< Time in us to start search
  double search_end_time;          ///< Time in us to start search
  ///parameters for variance search
  int start_window;                ///< sample window to start looking in
  double min_start_variance;       ///< min variance to define a start
  int min_resolution;              ///< What is this?
  ///parameters for discriminator search
  bool discriminator_relative;     ///< is disc value relative to baseline?
  double discriminator_value;      ///< discriminator treshold value in counts
  bool use_baseline_sigma;         ///< if we want to use the baseline fluctuation
  double discriminator_nsigma;     ///< discriminator treshold value in baseline sigmas
  int discriminator_start_add;     ///< n samples to add before start
  int discriminator_end_add;       ///< n samples to add after end
  //parameters for integral search
  // bool normalize_to_npe;           ///< Scale all searches by spemean?
  // double integral_start_time;      ///< time in us over which photons arrive
  // double integral_end_time;        ///< time at end of pulse to be below thresh
  // double integral_start_threshold; ///< amount in npe integral must fall
  // double integral_end_threshold;   ///< end when npe integral below thresh
  // double min_sep_time;             ///< minimum time two pulses must be apart
  // double multipulse_thresh_value;  ///< secondary must be > this*prev integral
  double amplitude_start_threshold;   ///< before broad integral search, look for signal above this threshold
  // double amplitude_end_threshold;  ///< signal must fall below this threshold to end pulse
  // double min_pulse_time;           ///< minimum length of pulse before ends naturally
  // double lookback_time;            ///< samples to look back for pileup
  ///parameters for integral search by jingke
  int integral_window_nsamps;      ///< number of samples in the integral window
  int moving_step_nsamps;          ///< number of samples to move after each integral step
  double integral_start_threshold; ///< integral threshold to start a pulse
  double integral_end_threshold;   ///< integral threshold to end a pulse
  double pileup_factor_rel;        ///< for a pileup if the later integral has to be X times larger
  double pileup_factor_abs;        ///< for a pileup if the later integral has to be X value larger
  //parameters for curvature search
  int down_sample_factor;          ///< reduce the integral vector size by this factor 
  double pulse_start_curvature;    ///< curvature threshold to start a new pulse  
  int pile_up_curvature;           ///< curvature threshold to start a pile up pulse
  double pulse_end_slope;          ///< slope threshold to end a pulse
};

//override stream ops for SearchMode
std::istream& operator>>(std::istream& in, PulseFinder::SEARCH_MODE& m);
std::ostream& operator<<(std::ostream& out, const PulseFinder::SEARCH_MODE& m);

#endif
