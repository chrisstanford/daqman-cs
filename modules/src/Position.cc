#include "Position.hh"
#include "utilities.hh"
#include "ConvertData.hh"
#include "BaselineFinder.hh"
#include "Integrator.hh"
#include "SumChannels.hh"
#include "PulseUtility.hh"
#include "intarray.hh"
#include "TMath.h"
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include <math.h>
#include <fstream>
#include <cstdio>

Position::Position() : 
  BaseModule(GetDefaultName(), "Search for individual physical scintillation pulses within the hardware trigger")
{
  AddDependency<ConvertData>();
  AddDependency<BaselineFinder>();
  AddDependency<Integrator>();
  AddDependency<SumChannels>();
  
  RegisterParameter("n_pmts", n_pmts = 4,
 		    "Total number of PMTs, the positions need to be read from input text file");
//   RegisterParameter("",  = ,
// 		    "");

}
  
Position::~Position(){;}

int Position::Initialize(){

  if(n_pmts<1) return -1;
  max_lnr = -10;

  const char pos_file[32]="inputs/PMT_Positions.txt";

  std::vector<int> pmt_id;
  std::vector<double> pmt_x;
  std::vector<double> pmt_y;
  std::vector<double> pmt_z;
  pmt_id.resize(n_pmts);
  pmt_x.resize(n_pmts);
  pmt_y.resize(n_pmts);
  pmt_z.resize(n_pmts);
  int num = n_pmts;
  //  num = LoadArray<int>(num, &(pmt_id[0]), pos_file, 1, 1);
  num = LoadArray(num, &(pmt_id[0]), pos_file, 1, 1);
  if(num<n_pmts){
    Message(ERROR)<<pos_file<<" does not contain as many PMT positions! \n";
    return -1;
  }
  num = LoadArray(num, &(pmt_x[0]), pos_file, 1, 2);
  num = LoadArray(num, &(pmt_y[0]), pos_file, 1, 3);
  num = LoadArray(num, &(pmt_z[0]), pos_file, 1, 4);

  for(size_t ii=0; ii<pmt_id.size(); ii++){
    std::cout<<"PMT #"<<pmt_id[ii]<<'\t'<<pmt_x[ii]<<'\t'<<pmt_y[ii]<<'\t'<<pmt_z[ii]<<std::endl;
    XYZ pos;
    pos.x=pmt_x[ii];
    pos.y=pmt_y[ii];
    pos.z=pmt_z[ii];
    pmts_pos[pmt_id[ii]] = pos;
  }//end for 

  return 0;
}

int Position::Finalize(){
  return 0;
}

int Position::Process(EventPtr evt){
  
  EventDataPtr event = evt->GetEventData();
  //only process the sum channel for event positions
  ChannelData *chsum = event->GetChannelByID(ChannelData::CH_SUM);
  if(!chsum) return -1;
  
  std::set<int> &channels_summed = event->channels_summed;
  std::vector<Pulse> & sum_pulses = chsum->pulses;
  id_map pmt_hits_integral;
  id_map pmt_hits_spikes;
  for(size_t jj=0; jj<sum_pulses.size(); jj++){
    if(!(sum_pulses.at(jj).evaluated)) continue;
    for(std::set<int>::iterator it=channels_summed.begin(); it!=channels_summed.end(); it++){
      ChannelData * chdata = event->GetChannelByID(*it);
      if(jj<chdata->pulses.size() && chdata->pulses.at(jj).evaluated){
	pmt_hits_integral[*it] = std::fabs(chdata->pulses.at(jj).integral);
	pmt_hits_spikes[*it] = std::fabs(chdata->pulses.at(jj).nspikes);
      }//end if
      else{
	pmt_hits_integral[*it] = 0;
	pmt_hits_spikes[*it] = 0;
      }//end else
    }//end for it
    //calculate the position parameters
    GravityCenter(sum_pulses.at(jj).pos_com_integral, pmt_hits_integral);
    GravityCenter(sum_pulses.at(jj).pos_com_spikes, pmt_hits_spikes);
    LogRatio(sum_pulses.at(jj).pos_lnr_integral, pmt_hits_integral);
    LogRatio(sum_pulses.at(jj).pos_lnr_spikes, pmt_hits_spikes);

  }//end for jj
  
  return 0;
}

//WARNING: if we have a bottom PMT, this needs to be modified
void Position::GravityCenter(XYZ &pos, id_map hits){

  pos.Clear();
  double all_hits=0;
  for(id_map::iterator it = hits.begin(); it != hits.end(); it++){
    all_hits += it->second;
    ip_map::iterator iter = pmts_pos.find(it->first);
    if(iter == pmts_pos.end()){
      Message(ERROR)<<"Position::GravityCenter PMT position for channel "<<it->first<<" not defined! \n";
      continue;
    }
    pos.x += it->second * iter->second.x;
    pos.y += it->second * iter->second.y;
    pos.z += it->second * iter->second.z;
  }//end for
  
  if(all_hits>0){
    pos.x /= all_hits;
    pos.y /= all_hits;
    pos.z /= all_hits;
  }
  else pos.Clear();
}

//WARNING: this only works for the symmetric 4 pmt layout!!!
void Position::LogRatio(XYZ &pos, id_map hits){
  if(hits.size()!=4) return;
  
  XYZ neg_counts;
  XYZ pos_counts;
  //considering finding the center of the array first
  for(id_map::iterator it = hits.begin(); it != hits.end(); it++){
    ip_map::iterator iter = pmts_pos.find(it->first);
    if(iter == pmts_pos.end()){
      Message(ERROR)<<"Position::GravityCenter PMT position for channel "<<it->first<<" not defined! \n";
      continue;
    }
    if(iter->second.x<0) neg_counts.x += it->second;
    else if(iter->second.x>0) pos_counts.x += it->second;
    if(iter->second.y<0) neg_counts.y += it->second;
    else if(iter->second.y>0) pos_counts.y += it->second;
    if(iter->second.z<0) neg_counts.z += it->second;
    else if(iter->second.z>0) pos_counts.z += it->second;
  }//end for
  
  CalculateLnR(pos.x, neg_counts.x, pos_counts.x);
  CalculateLnR(pos.y, neg_counts.y, pos_counts.y);
  CalculateLnR(pos.z, neg_counts.z, pos_counts.z);
  return;
}

void Position::CalculateLnR(double &lnr, double neg, double pos){
  if(neg>0) neg=log(neg);
  else neg = max_lnr;
  if(pos>0) pos=log(pos);
  else pos = max_lnr;
  lnr = pos-neg;
}
