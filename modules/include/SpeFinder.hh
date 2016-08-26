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

  typedef std::map<int,double> id_map;
  
  //parameters
  double search_start_time;                ///< Time from start of pulse to begin search [us]
  double search_end_time;                  ///< Time from start of pulse to end search [us]
  //  double threshold_nsigmas;            ///< number of baseline standard deviation for a pulse
  std::map<int,double> ch_spe_thresholds;  ///< discriminator threshold values for spike finding
  bool relative_bls_threshold;             ///< Is the threshold relative to baseline sigma
  bool relative_spe_threshold;             ///< Is the threshold relative to spe size
  int spe_edge_add_nsamps;                 ///< samples to be added to the edge of a pulse
  double min_baseline_length_us;           ///< minimum length of baseline before and after pulses in us

};
#endif
