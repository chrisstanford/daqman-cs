/** @file TOF.hh
    @brief Evaluate the TimeOfFlight
    @author jingke
    @ingroup modules
*/
#ifndef TOF_h
#define TOF_h

#include "ChannelModule.hh"

/** @class TOF
    @brief Searches for valid pulses in specified time windows and combine them for TOF study
    @ingroup modules
*/

class TOF : public ChannelModule{
public:
  TOF();
  ~TOF();
  
  int Initialize();
  int Finalize();
  int Process(ChannelData* chdata);
  
  int ProcessRefChannel(ChannelData* chdata);
  int ProcessDataChannel(ChannelData* chdata);
  
  static const std::string GetDefaultName(){ return "TOF"; }
  
private:
  //parameters
  int ref_ch;
  int signal_ch;
  double ref_ch_offset;
  double signal_start_time;
  double neutron_start_time;
  double search_end_time;
  double pulse_length_us;
};
#endif
