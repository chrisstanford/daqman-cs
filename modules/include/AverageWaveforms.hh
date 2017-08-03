/** @file AverageWaveforms.hh
    @brief defines the AverageWaveforms module
    @author bloer, jingkexu
    @ingroup modules
*/

#ifndef AVERAGE_WAVEFORMS_h
#define AVERAGE_WAVEFORMS_h

#include "ChannelModule.hh"
#include "TGraph.h"
#include "TTree.h"
#include <map>
#include <string>

class TGraphErrors;

/** @class AverageWaveforms
    @brief Averages the signals for each channel over an entire run (with some basic cuts)
    @ingroup modules
    jingke makes this module a channel module, May 25 2014
*/
class AverageWaveforms : public ChannelModule{
public:
  AverageWaveforms();
  ~AverageWaveforms();

  int Initialize();
  int Finalize();
  int Process(ChannelData* chdata);
  static const std::string GetDefaultName(){ return "AverageWaveforms";}
  bool PassesCuts(const int run_id, ChannelData* channel, const double spe_int_mean, const double spe_height, TString mode, TString option);

  int ref_roi_index;
  int ref_fp_index;
  double aver_time_hr;
  double pre_trigger_time;
  double post_trigger_time;
  double min_pulse_height;
  double max_pulse_height;
  double min_npe;
  double max_npe;
  double min_fprompt;
  double max_fprompt;
  double min_tof;
  double max_tof;
  bool normalize_by_total_waveforms;
  bool one_pulse_only;
  bool bipo_mode;
  int bipo_pulse_id;
  int run_id;
  double spe_height_mean;
  double spe_int_mean;
  double pre_spe_nsamps;
  double post_spe_nsamps;
  double verbose;
  bool afterpulse_cut;
  double afterpulse_height_spe;
  //not set in config file

private:
  void Cleanup();
  std::map<int,TGraphErrors*> _plots;
  std::map<int,double> _sum;
  std::map<int,double> _sum2;

  std::map<int,int> _total_waveforms;
  
  TGraphErrors* waveform; // current waveform
  TGraphErrors* spewf; // current waveform
  TTree* waveforms; // all waveforms to be averaged
  TTree* spewfs; // all spes to be averaged
};

#endif
