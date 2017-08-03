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
    peak_time = -1;
    end_time = -1;
    width = -1;
    start_index=-1;
    peak_index=-1;
    end_index=-1;
    amplitude=-1;
    amplitude_sigma =-1;
    baseline_mean = -1;
    baseline_sigma = -1;
  }
  
  double integral;       ///< Integral in counts*samples of the found charge
  double start_time;     ///< time of the start of the pulse
  double peak_time;     ///< time of the peak of the pulse
  double end_time;     ///< time of the start of the pulse
  int    width;          ///< with of the spe pulse including the edge samples
  int start_index;      ///< index at which the pulse starts
  int peak_index;      ///< index at which the pulse reached max amplitude
  int end_index;      ///< index at which the pulse ends
  double amplitude;      ///< Maximum height of the spe pulse
  double amplitude_sigma;      ///< Maximum height of the spe pulse in terms oc global baseline sigma
  double baseline_mean;  ///< value of the local baseline mean
  double baseline_sigma; ///< value of the local baseline standard deviation
  ClassDef(Spe,3)
};
  
#endif
