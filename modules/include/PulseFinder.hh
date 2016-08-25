/** @file PulseFinder.hh
    @brief Defines the PulseFinder module
    @author jingke
    @ingroup modules
*/
#ifndef PULSEFINDER_h
#define PULSEFINDER_h

#include "BaseModule.hh"
#include <vector>

/** @class PulseFinder
    @brief Searches for pulses defined by the users
    @ingroup modules
*/
class PulseFinder : public BaseModule{
public:
  PulseFinder();
  ~PulseFinder();
   
  int Initialize();
  int Finalize();
  int Process(EventPtr evt);
  static const std::string GetDefaultName(){ return "PulseFinder"; }
  
  typedef std::map<int,double> id_map;
  //for now I just want to use the discriminator search -- Jingke
  //we could smoothen the waveform and downsample to achieve complex algorithms
  //keep these in case we want to go back
  //  enum SEARCH_MODE { VARIANCE , DISCRIMINATOR , INTEGRAL , CURVATURE };

private:
  int FindChannelPulses(ChannelData* chdata);
  int FindChannelSpikes(ChannelData* chdata);
  /// Search for pulses using a simple discrimination threshold
  int DiscriminatorSearch(ChannelData* chdata, const double * wave,
			  std::vector<int>& start_index, 
			  std::vector<int>& end_index, double threshold,
			  int pulse_start_add_nsamps, int pulse_end_add_nsamps);
  /// resolve pileup pulses/spikes
  /// only trust the result if  within a pulse/spike
  /// may pickup baseline fluctuations otherwise
  int ResolvePileUps(ChannelData* chdata, const double * wave,
		      std::vector<int>& start_index, 
		      std::vector<int>& end_index,
		      double pileup_threshold, int step=3, 
		      int search_start=0, int search_end=-1);
  ///parameters
  //  SEARCH_MODE mode;                ///< Which search function to use
  double search_start_time;                  ///< Time in us to start search
  double search_end_time;                    ///< Time in us to end search
  ///parameters for discriminator search
  std::map<int,double> ch_pulse_filter_us;   ///< estimate of the pulse width in us
  std::map<int,double> ch_pulse_thresholds;  ///< discriminator threshold values for pulse finding
  std::map<int,double> ch_spike_thresholds;  ///< discriminator threshold values for spike finding
  bool relative_bls_threshold;               ///< Is the threshold relative to baseline sigma
  bool relative_spe_threshold;               ///< Is the threshold relative to spe size
  //if both relative to bls and spe are set to be true, initilize will fail
  int pulse_start_add_us;                    ///< time in us to add before pulse start
  int pulse_end_add_us;                      ///< time in us to add after pulse end
  //  int spike_edge_add_nsamps;                 ///< n samples to add on spike edge
};

//override stream ops for SearchMode
// std::istream& operator>>(std::istream& in, PulseFinder::SEARCH_MODE& m);
// std::ostream& operator<<(std::ostream& out, const PulseFinder::SEARCH_MODE& m);

#endif
