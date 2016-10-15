/** @file Spe.hh
    @brief Defines the Spe storage class
    @ingroup modules
    @ingroup daqroot
    @author rsaldanha
*/

#ifndef SPE_h
#define SPE_h
#include "TObject.h"

/** @class Spe
    @brief Store information about a single photoelectron pulse in the root tree
    @ingroup modules
    @ingroup daqroot
*/
class Spe
{
public:
  Spe()
  {
    Clear();
  }
  virtual ~Spe() {}
  void Clear()
  {
    integral = -1;
    start_time = -1;
    width = -1;
    peak_time=-1;
    amplitude=-1;
    baseline_mean = -1;
    baseline_sigma = -1;
  }
  
  double integral;       ///< Integral in counts*samples of the found charge
  double start_time;     ///< time of the start of the pulse
  int    width;          ///< with of the spe pulse including the edge samples
  double peak_time;      ///< time at which the pulse reached max amplitude
  double amplitude;      ///< Maximum height of the spe pulse
  double baseline_mean;  ///< value of the local baseline mean
  double baseline_sigma; ///< value of the local baseline standard deviation
  ClassDef(Spe,3)
};
  
#endif
