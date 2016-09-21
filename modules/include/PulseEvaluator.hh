/** @file PulseEvaluator.hh
    @brief Defines PulseEvaluator module that evaluate pulse/ROI quantities 
    @author jingkexu
    @ingroup modules
*/

#ifndef PULSEEVALUATOR_h
#define PULSEEVALUATOR_h

#include "ChannelModule.hh"
#include "Pulse.hh"

#include <vector>
#include <set>
/** @class PulseEvaluator
    @brief evaluate pulse/ROI quantities
    @ingroup modules
    */
class ChannelData;

class PulseEvaluator : public ChannelModule
{
public:
  PulseEvaluator();
  ~PulseEvaluator();
  
  int Initialize();
  int Finalize();
  int Process(ChannelData* chdata);
  
  static const std::string GetDefaultName(){ return "PulseEvaluator";}

  //parameters
  bool evaluate_pulses;            /// whether pulse quantities should be evaluated
  bool evaluate_rois;              /// whether roi quantities should be evaluated
  double max_peak_time_us;         /// maximum peak time in us from the beginning of the pulse

private:
  //function to calculate pulse parameters, if max_peak_time>0
  //the peak of the pulse should be within (start time, start time + max_peak_time)
  int EvaluatePulse(Pulse& pulse, ChannelData* chdata, double max_peak_time=-1);

};
#endif
