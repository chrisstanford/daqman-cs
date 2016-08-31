/** @file Pulse.hh
    @brief Defines the Pulse data class
    @author bloer, rsaldanha, jingkexu
    @ingroup daqroot
    @ingroup modules
*/

#ifndef PULSE_h
#define PULSE_h

#include "Message.hh"
#include "Spike.hh"
#include <vector>
#include <iostream>

/** @class Pulse
    @brief Stores information relevant to a single scintillation event
    @ingroup modules
    @ingroup daqroot
*/
class Pulse{
public:
  Pulse() { Clear(); }
  virtual ~Pulse() {  }
  void Clear();
  void Print(int channel, int pulse);
public:
  int start_index;       ///< index marking the start of the pulse
  int peak_index;        ///< index for the pulse peak
  int end_index;         ///< index for the end of the pulse
  double start_time;     ///< time since trigger at start of pulse
  double peak_time;      ///< time since trigger at pulse peak;
  double end_time;       ///< time since trigger at end of pulse
  double half_max_time;  ///< time when pulse rise to 50% of peak
  double max;            ///< Maximum value obtained over the region
  double min;            ///< Minimum value obtained over the region
  double overshoot; ///< overshoot of the peak
  double peak_amplitude; ///< amplitude of the peak
  double tail_peak;      ///< peak amplitude of the pulse tail
  double integral;       ///< integral of the peak
  double npe;            ///< integral scaled for single pe amplitude
  double time_since_last; ///< time from the end of the last pulse
  std::vector<double> fparameters; ///< fparameters for given prompt time
  std::vector<double> tparameters; ///< time for integral to reach certain fractions
  std::vector<Spike>  spikes; ///< spikes in this pulse
  double gatti;    ///< gatti parameter for this pulse
  double mean_time; ///< mean time of this pulse for PSD studies
  double test;     ///< test variable
  bool saturated;  ///< was the peak saturated?
  bool is_clean;   ///< clean = this pulse is not back to back with another, not cut off by the end of the scan, not an s2 trigger 
  bool evaluated;  ///< if this pulse has been evaluated
  
  ClassDef(Pulse,11);
};

inline void Pulse::Clear()
{
  start_index = -1;
  end_index = -1;
  peak_index = -1;
  start_time = 0;
  half_max_time = -10000;
  peak_time = 0;
  end_time = 0;
  max = -1;
  min = -1;
  overshoot = -1;
  peak_amplitude = -1;
  tail_peak = -1;
  integral=0;
  npe=0;
  time_since_last = 0;
  if(fparameters.size()) fparameters.clear();
  if(tparameters.size()) tparameters.clear();
  if(spikes.size()) spikes.clear();
  gatti = 0;
  mean_time = -1;
  test = 0;
  saturated = false;
  is_clean  = false;
  evaluated = false;
};

inline void Pulse::Print(int channel, int pulse)
{
    Message m(INFO);
    m<<std::endl;
    m<<"************************************************************************"<<std::endl;
    m<<"******************* CHANNEL "<<channel<<" PULSE "<<pulse<<"  INFORMATION *******************"<<std::endl;
    m<<"************************************************************************"<<std::endl;
    m<<"Start Index: "<<start_index<<std::endl;
    m<<"End Index: "<<end_index<<std::endl;
    m<<"Peak Index: "<<peak_index<<std::endl;
    m<<"Start Time: "<<start_time<<std::endl;
    m<<"Half Maximum Time: "<<half_max_time<<std::endl;
    m<<"Peak Time: "<<peak_time<<std::endl;
    m<<"End Time: "<<end_time<<std::endl;
    m<<"Maximum Value: "<<max<<std::endl;
    m<<"Minimum Value: "<<min<<std::endl;
    m<<"Overshoot Value: "<<overshoot<<std::endl;
    m<<"Peak Amplitude: "<<peak_amplitude<<std::endl;
    m<<"Tail Peak Amplitude: "<<tail_peak<<std::endl;
    m<<"Integral: "<<integral<<std::endl;
    m<<"Npe: "<<npe<<std::endl;
    m<<"Time since last pulse: "<<time_since_last<<std::endl;
    m<<"Size of fparameters: "<<fparameters.size()<<std::endl;
    m<<"Size of tparameters: "<<tparameters.size()<<std::endl;
    m<<"Gatti Variable: "<<test<<std::endl;
    m<<"Mean Time Variable: "<<test<<std::endl;
    m<<"Test Variable: "<<test<<std::endl;
    m<<"Pulse Saturated: "<<saturated<<std::endl;
    m<<"Is Clean: "<<is_clean<<std::endl;
    m<<"Pulse Evaluated: "<<evaluated<<std::endl;
    m<<"************************************************************************"<<std::endl;
}
#endif
