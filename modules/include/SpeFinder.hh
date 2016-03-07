/** @file SpeFinder.hh
    @brief Defines the SpeFinder module
    @author rsaldanha, jingkexu
    @ingroup modules
*/

#ifndef SPEFINDER_h
#define SPEFINDER_h

#include "ChannelModule.hh"

/** @class SpeFinder
    @brief Search for single-photoelectron events in the tails of scintillation
    @ingroup modules
*/
class SpeFinder : public ChannelModule
{
public:
  SpeFinder();
  ~SpeFinder();
  
  int Initialize();
  int Finalize();
  int Process(ChannelData* chdata);
  
  static const std::string GetDefaultName(){ return "SpeFinder";}
  
  //parameters
  double search_start_time; ///< Time from start of pulse to begin search [us]
  double search_end_time;   ///< Time from start of pulse to end search [us]
  double threshold_nsigmas; ///< number of baseline standard deviation for a pulse
  int pulse_edge_add;       ///< samples to be added to the edge of a pulse
  double min_baseline_length; ///< minimum length of baseline before and after pulses in us

};
#endif
