/** @file PSD.hh
    @brief defines the PSD module
    @author jingkexu
    @ingroup modules
*/

#ifndef PSD_h
#define PSD_h

#include "ChannelModule.hh"
#include <vector>
//#include <map>
#include <string>

class TH1F;

/** @class PSD
    @brief calculate the PSD parameters gatti here
    @ingroup modules
*/
class PSD : public ChannelModule{
public:
  PSD();
  ~PSD();

  int Initialize();
  int Finalize();
  int Process(ChannelData* chdata);
  static const std::string GetDefaultName(){ return "PSD";}

  //fparameter calculations
  std::vector<double> fparameter_times;
  //tparameter calculations
  std::vector<double> tparameter_ratios;

  //gatti parameters
  bool calculate_gatti;
  std::set<int> gatti_channels;
  double start_time_us;
  double end_time_us;
  double bin_width_us;
  std::string wfm_file;

private:
  void Cleanup();
  void ProcessPulse(ChannelData *chdata, Pulse & pulse);

  TH1F * weights;
  TH1F * rebined_pulse;

};

#endif
