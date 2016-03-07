#ifndef PEAK_h
#define PEAK_h

#include "Message.hh"
#include <vector>
#include <iostream>

/** @class Peak
    @brief Stores information relevant to a single peak
    @ingroup modules
    @ingroup daqroot
*/
class Peak{
public:
  Peak() { Clear(); }
  virtual ~Peak() {;}
  void Clear();
public:
  double peak_time;     ///< time since trigger at start of peak
  double peak_amplitude; ///< amplitude of the peak
  double integral;       ///< integral of the peak
  
  ClassDef(Peak,1);
};

inline void Peak::Clear()
{
  peak_time = 0;
  peak_amplitude = -1;
  integral=0;
};

#endif
