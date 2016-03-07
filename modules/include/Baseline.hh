/** @file Baseline.hh
    @brief Defines the Baseline data class
    @author jingkexu
    @ingroup modules
*/

#ifndef BASELINE_h
#define BASELINE_h

/** @class Baseline
    @brief Variables related to the found baseline for a channel
    @ingroup modules
*/
class Baseline{
public:
  Baseline(){ Clear(); }
  virtual ~Baseline() {}
  /// Reset all to default
  void Clear();
public:
  bool found_baseline;    ///< was the baseline finder successful?
  double mean;            ///< the average baseline found
  double sigma;           ///< standard deviation of the samples in the baseline region
  int length;             ///< num of samples over which baseline was averaged
  bool saturated;         ///< Did the baseline hit the max or min y value?
  
  ClassDef(Baseline,3)

};

inline void Baseline::Clear()
{
  found_baseline = false;
  mean = -1;
  sigma = -1;
  length = -1;
  saturated = false;
}

#endif
