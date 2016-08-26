/** @file EventData.hh
    @brief Defines the EventData storage class
    @author bloer
    @ingroup modules
    @ingroup daqroot
*/

#ifndef EVENTDATA_h
#define EVENTDATA_h

#include <vector>
#include <string>
#include <stdint.h>
#include <iomanip>

#include "Rtypes.h" //has the classdef macro
#include "ChannelData.hh"

//notice: members with comment starting with ! are not saved

/** @class EventData
    @brief Processed data for a trigger which is stored in the ROOT tree
    @ingroup modules
    @ingroup daqroot
*/
class EventData{
  //interface members
public:
  EventData() { Clear(); }
  virtual ~EventData() {} //anything need cleaning up?
  /// Reset all variables to defaults
  void Clear(); //inlined below
  static const char* GetBranchName(){ return "event"; }
  void Print (int verbosity);
public:
  
  enum STATUS_FLAGS { NORMAL=0, ID_MISMATCH=1, /*enter more here*/};
  
  //data members
  int run_id; ///< id of this run
  int event_id; ///< id of this event withing the run
  uint64_t status;             ///< flag denoting event status
  int trigger_count; ///< number of triggers (including unaccepted triggers)
  long timestamp; ///< unix timestamp for the event
  uint64_t dt; ///< time since the last event in ns
  uint64_t event_time; ///< time since run start in ns
  int nchans; ///< physical channels that lit up 
  bool saturated; ///< true if any channel hit the limit of its digitizer range
  bool pulses_aligned; ///< true if the pulse edges on all channels are aligned
  std::set<int> channels_summed; ///< channels that contribute to the sum channel
  std::vector<double> generic; ///< vector of generic results for testing
  std::vector<ChannelData> channels;  ///< results for each channel
  
  /// Return channel info pointer by id, rather than vector index
  ChannelData* GetChannelByID(int id){
    std::vector<ChannelData>::iterator it = channels.begin();
    while(it != channels.end() && it->channel_id != id) it++;
    return (it == channels.end() ? 0 : &(*it) );
  }
  
  ///Get a pointer to the sum s1 pulse
  Pulse* GetPulse(size_t pulsenum, int channel_id = ChannelData::CH_SUM){ 
    ChannelData* ch = GetChannelByID(channel_id);
    if(ch && ch->pulses.size() > pulsenum)
      return &(ch->pulses[pulsenum]);
    return 0;
  }
  Pulse* GetROI(size_t region, int channel_id = ChannelData::CH_SUM){
    ChannelData* ch = GetChannelByID(channel_id);
    if(ch && ch->regions.size() > region)
      return &(ch->regions[region]);
    return 0;
  }
  
  ClassDef(EventData,16)
};

inline void EventData::Clear()
{
  run_id = -1;
  event_id = -1;
  status = NORMAL;
  trigger_count = -1;
  timestamp = 0;
  dt = 0;
  event_time = 0;
  nchans = -1;
  saturated = false;
  pulses_aligned = false;
  generic.clear();
  channels.clear();
}
#endif
