/** @file Position.hh
    @brief Defines the Position reconstruction module
    @author jingke
    @ingroup modules
*/
#ifndef POSITION_h
#define POSITION_h

#include "BaseModule.hh"
#include "Pulse.hh"
#include <map>
//#include <vector>
//#include <utility>

/** @class Position
    @brief Calculate position of the pulses
    @ingroup modules
*/
class Position : public BaseModule{
public:
  Position();
  ~Position();
   
  int Initialize();
  int Finalize();
  int Process(EventPtr evt);
  static const std::string GetDefaultName(){ return "Position"; }
  
  typedef std::map<int,double> id_map;
  typedef std::map<int,XYZ> ip_map;

private:

  void GravityCenter(XYZ &pos, id_map hits);
  void LogRatio(XYZ &pos, id_map hits);
  void CalculateLnR(double &lnr, double neg, double pos);

  ///parameters
  int n_pmts;                ///< Total number of PMTs, the positions need to be read from input text file
  ip_map pmts_pos;           ///< positions of the PMTs
  double max_lnr;            ///< maximum valve of log(ratio) when ratio approaches 0 or infinity

};

#endif
