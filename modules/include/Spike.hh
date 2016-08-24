#ifndef SPIKE_h
#define SPIKE_h

#include "Message.hh"
#include <vector>
#include <iostream>

/** @class Spike
    @brief Stores information relevant to a single spike
    @ingroup modules
    @ingroup daqroot
*/
class Spike{
public:
  Spike() { Clear(); }
  virtual ~Spike() {;}
  void Clear();
public:
  double start_time;     ///< time at start of spike in us
  double peak_time;      ///< time at peak of spike in us
  double width;          ///< width of spike from start to end
  double peak_amplitude; ///< amplitude of the spike
  double integral;       ///< integral of the spike
  
  ClassDef(Spike,1);
};

inline void Spike::Clear()
{
  spike_time = 0;
  peak_amplitude = -1;
  integral=0;
};

#endif
