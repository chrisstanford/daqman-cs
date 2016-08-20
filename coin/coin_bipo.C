#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <exception>
#include <fstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
//#include <numeric>
//#include <map>

#include "EventData.hh"
#include "Pulse.hh"
//#include "/Volumes/Users/jingke/LowBackground/OfficialCodes/CommonFunc.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
//#include "T.h"
//#include "TGraph.h"
//#include "TString.h"
//#include "TROOT.h"
//#include "TSystem.h"

using std::vector;

//for run150
// const int sig_ch = 3;
// const int att_ch = 4;
// const double spe = 28.9;

//for run157
const int sig_ch = 3;
const int att_ch = 4;
const double spe = 29.4;

//const double spe = 42.; //run 192
//const double triplet_t = 1.35;


//const double spe = 34.;


//only allow alphabet, numbers,_, and . to show in filenames
bool ValidString(const char *fname){
  std::string file_name(fname);
  static const char valid_letters[100]="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890/_.";
  return(file_name.find_first_not_of(valid_letters)==std::string::npos);
}

int GetRunNumber(const char *fname){
  std::string file_name(fname);
  size_t pos = file_name.find("Run");
  if(pos==std::string::npos) return -1;
  pos = file_name.find_first_of("1234567890", pos);
  if(pos==std::string::npos) return -1;
  return atoi(file_name.substr(pos,file_name.length()-pos).c_str());

}

Double_t linterp(double x, int n, const double xx[], const double yy[]){
  double y=0;
  if(x<=xx[0]) y = yy[0];
  else if(x>=xx[n-1]) y = yy[n-1];
  else{
    int i;
    for(i=0; i<n; i++) if(x<xx[i]) break;
    double r=(x-xx[i-1])/(xx[i]-xx[i-1]);
    y = r*yy[i]+(1.-r)*yy[i-1];
  }
  return y;
}

bool IsGamma(double npe, double fp){
  const int N = 9;
  //for run 123 on June 12
  // const double npes[N] = {0,   200, 400,  600,  1000, 5000, 10000, 15000, 20000};
  // const double low[N]  = {0,   0,   0.1,  0.1,  0.15, 0.15, 0.15,  0.16,  0.14};
  // const double high[N] = {0.7, 0.6, 0.55, 0.45, 0.45, 0.4,  0.3,   0.26,  0.26};
  //for run 137 on June 19
  // const double npes[N] = {0,   200, 400,  600,  1000, 5000, 10000, 15000, 20000};
  // const double low[N]  = {0,   0.15,0.15, 0.15, 0.15, 0.15, 0.1,   0.1,  0.1};
  // const double high[N] = {0.6, 0.5, 0.45, 0.45, 0.45, 0.35, 0.25,  0.2,  0.2};
  //for run 150 on July 17
  // const double npes[N] = {0,   100,   200,   300,   500,   5000,  10000, 15000, 25000};
  // const double low[N]  = {0,   0,     0,     0,     0,     0.2,   0.2,   0.2,   0.2};   
  // const double high[N] = {0.6, 0.5,   0.48,  0.46,  0.46,  0.45,  0.4,   0.4,   0.4};
  //for run 157 on July 21
  const double npes[N] = {0,   100,   200,   300,   500,   5000,  10000, 15000, 25000};
  const double low[N]  = {0,   0,     0,     0,     0,     0.2,   0.2,   0.2,   0.2};   
  const double high[N] = {0.6, 0.5,   0.48,  0.46,  0.42,  0.32,  0.32,  0.32,  0.32};
//for run 192 on Oct 2
//  const double npes[N] = {0,   100,   200,   300,   500,   5000,  10000, 15000, 25000};
//  const double low[N]  = {0,   0,     0,     0,     0,     0.2,   0.2,   0.2,   0.2};   
//  const double high[N] = {0.6, 0.5,   0.48,  0.46,  0.46,  0.45,  0.4,   0.4,   0.4};
  double glow  = linterp(npe, N, npes, low);
  double ghigh = linterp(npe, N, npes, high);
  return (fp>glow) && (fp<ghigh);
}

bool IsAlpha(double npe, double fp){
  const int N = 9;
  //for run 123 on June 12
  // const double npes[N] = {0,   200, 400, 600, 1000, 5000, 10000, 15000, 20000};
  // const double low[N]  = {0.4, 0.4, 0.4, 0.4, 0.45, 0.4,  0.3,   0.26,  0.26};
  // const double high[N] = {0.9, 0.9, 0.9, 0.9, 0.8,  0.8,  0.8,   0.8,   0.8};
  //for run 137 on June 19
  // const double npes[N] = {0,   200,  400,  600,  1000, 5000, 10000, 15000, 20000};
  // const double low[N]  = {0.2, 0.35, 0.4,  0.4,  0.4,  0.35, 0.25,  0.18,  0.16};
  // const double high[N] = {0.9, 0.85, 0.85, 0.85, 0.85, 0.6,  0.35,  0.35,  0.26};
  //for run 150 on July 17
  // const double npes[N] = {0,   100,   200,   300,   500,   5000,  10000, 15000, 25000};
  // const double low[N]  = {0.3, 0.4,   0.44,  0.46,  0.46,  0.45,  0.4,   0.4,   0.4};
  // const double high[N] = {0.9, 0.9,   0.85,  0.85,  0.85,  0.85,  0.85,  0.85,  0.85};
  //for run 157 on July 21
//   const double npes[N] = {0,   100,   200,   300,   500,   5000,  10000, 15000, 25000};
//   const double low[N]  = {0.3, 0.35,  0.4,   0.42,  0.44,  0.45,  0.4,   0.4,   0.4};
//   const double high[N] = {0.9, 0.9,   0.85,  0.85,  0.85,  0.85,  0.85,  0.85,  0.85};
//for run 192 on Oct 2
  const double npes[N] = {0,   100,   200,   300,   500,   5000,  10000, 15000, 25000};
  const double low[N]  = {0.3, 0.4,  0.45,   0.45,  0.45,  0.45,  0.4,   0.4,   0.4};
  const double high[N] = {0.95,0.95,  0.85,  0.85,  0.85,  0.85,  0.85,  0.85,  0.85};
  double glow  = linterp(npe, N, npes, low);
  double ghigh = linterp(npe, N, npes, high);
  return (fp>glow) && (fp<ghigh);
}

void search(const char *fname){
  
  std::string cmd(fname);
  cmd = "[[ -f " + cmd +" ]]";
  if(system(cmd.c_str())){
    std::cerr<<fname<<" does not exist! Check the filename!"<<std::endl;
    return;
  }

  TFile *fin = new TFile(fname);
  if(!fin || !(fin->IsOpen()) || fin->IsZombie()){
    std::cerr<<fname<<" can not be opened!"<<std::endl;
    return;
  }

  TTree *Events = (TTree *)(fin->Get("Events"));
  if(!Events || Events->IsZombie() || Events->GetEntries()<1){
    std::cout<<"No data available. "<<std::endl;
    return;
  }

  EventData* evt = 0;
  Events->SetBranchAddress("event",&evt);

  TH2F *hgamma = new TH2F("GammaCoin", "Fp vs NPE for Gamma", 1000,0,20000, 220,0,1.1);
  TH2F *halpha = new TH2F("AlphaCoin", "Fp vs NPE for Alpha", 5000,0,25000, 220,0,1.1);
  TH1F *hdelay = new TH1F("Tdelay", "Delay time between coincidence", 1000, 0, 1e6);
  TH1F *ht     = new TH1F("Tcoin", "Time of conincidence", 500, 0, 13000);
  //  TH1F *ha     = new TH1F("AlphaE", "Alpha Energy with corrections", 25000, 0, 25000);
  hgamma->SetDirectory(0);
  halpha->SetDirectory(0);
  hdelay->SetDirectory(0);
  ht->SetDirectory(0);
  //  ha->SetDirectory(0);
  TH2F *hmt = new TH2F("AlphaMT", "MT vs NPE for Alpha", 2500,0,25000, 200,0,0.5);
  hmt->SetDirectory(0);
  TH2F *hmtfp = new TH2F("AlphaMTFp", "Fp vs MT for Alpha", 500,0,10,220,0,1.1);
  hmt->SetDirectory(0);
  
  ChannelData* ch = NULL;
  ChannelData* ch_att = NULL;
  double bad_time = -1e10;
  vector<Pulse> regions;
  vector<double> npes;
  vector<double> fps;
  vector<double> eventtimes;
  vector<double> timesince;
  vector<int> event_ids;
  double max_coin_window = 1e6;//1ms
  int printout = 0;

  for(int entry=0; entry<Events->GetEntries(); entry++){
    //  for(int entry=0; entry<1000; entry++){
    Events->GetEntry(entry);
    //    if(evt->event_time*1e-9<2000) continue;
    //this is for Run157, reject high trigger rate time
    if(evt->event_time*1e-9<200 ||
       (evt->event_time*1e-9>695  && evt->event_time*1e-9<710) ||
       (evt->event_time*1e-9>3995 && evt->event_time*1e-9<4002)) continue;
    
    //process the stored events if there are any
    if(eventtimes.size()>0 && evt->event_time-eventtimes.at(0)>max_coin_window){
      //      if(printout++<50) std::cout<<"Number of saved entries: "<<eventtimes.size()<<'\t'<<evt->event_id<<std::endl;    
      //here comes our golden vent
      //      if(eventtimes.size()==2 && (timesince.at(1)>1e6 || timesince.at(1)==0) && regions.at(1).mean_time<6){	
      if(eventtimes.size()>=2 && (timesince.at(1)>1e6 || timesince.at(1)==0)) {
// 	 && ( npes.at(1)>100 || regions.at(1).mean_time<2)){	
	hgamma->Fill(npes.at(0), fps.at(0));

	//this is not really necessary if we require there are only two events
	bool largest_first = true;
	if(eventtimes.size()>2){
	  for(size_t ii=2; ii<eventtimes.size(); ii++){
	    if(npes.at(ii)>npes.at(1) ) largest_first = false;
	  }
	}
	if(largest_first && eventtimes.at(1)-eventtimes.at(0)>1.5e5 && fps.at(1)<0.98 
	   && regions.at(1).time_since_last>0.95
	   && (npes.at(1)>600 || (regions.at(1).fparameters[5]-regions.at(1).fparameters[2]>0.15 &&
				  regions.at(1).fparameters[5]-regions.at(1).fparameters[2]<0.4 &&
				  fps.at(1)>0.4 && fps.at(1)<0.75))
	   //	   && regions.at(1).tparameters[3]>0.3
	   ){
	  halpha->Fill(npes.at(1), fps.at(1));
	  //	  if(npes.at(1)<50 && printout++<50) 
	  if(npes.at(1)>150 &&npes.at(1)<250 && printout++<50) 
	    std::cout<<event_ids.at(1)<<'\t'<<npes.at(1)<<'\t'<<fps.at(1)<<'\t'
		     <<regions.at(1).mean_time<<'\t'<<regions.at(1).test<<std::endl;

	  //	  hmt->Fill(npes.at(1), regions.at(1).test);
	  //	  hmt->Fill(npes.at(1), regions.at(1).mean_time);
	  //	  hmt->Fill(npes.at(1), regions.at(1).tparameters[2]);
	  hmt->Fill(npes.at(1), regions.at(1).fparameters[5]-regions.at(1).fparameters[2]);
	  //this is to correct the npe based on lifetime
	  // if(npes.at(1)<100){
	  //   if(fps.at(1)>0.74 && fps.at(1)<0.9) ha->Fill(npes.at(1)*fps.at(1) + npes.at(1)*(1.-0.74)*1.55/triplet_t);
	  // }
	  // else{
	  //   ha->Fill(npes.at(1)*fps.at(1) + npes.at(1)*(1.-fps.at(1))*1.55/triplet_t);
	  // }

	  //	  hmt->Fill(npes.at(1), regions.at(1).mean_time);
	  //	  hmt->Fill(npe, region.end_time-region.start_time);
	  if(npes.at(1)<100) hmtfp->Fill(regions.at(1).mean_time, fps.at(1));
	}
	//	if(npes.at(1)<20 && fps.at(1)<0.6) 
	//	if(fps.at(1)>0.98) 
	//	if(regions.at(1).mean_time>2 && regions.at(1).mean_time<3) 
	hdelay->Fill((eventtimes.at(1)-eventtimes.at(0)));
	ht->Fill(eventtimes.at(0)*1e-9);
      }
      regions.clear();
      npes.clear();
      fps.clear();
      eventtimes.clear();
      timesince.clear();
      event_ids.clear();
    }

    ch = evt->GetChannelByID(sig_ch);
    if(!ch){
      std::cout<<"Signal channel does not exist in the channel data!"<<std::endl;
      break;
    }
    if(!(ch->baseline.found_baseline) || ch->tof.size()<1 ){
      //      std::cout<<"failed to find regionss in event "<<evt->event_id<<std::endl;
      continue;
    }
    //apply a cut to remove tails of large pulses
    if(ch->tof[0].tail_peak>400) bad_time = evt->event_time;

    ch_att = evt->GetChannelByID(att_ch);
    if(!ch_att){
      std::cout<<"Attenuated signal channel does not exist in the channel data!"<<std::endl;
      break;
    }
    
    Pulse region = ch->tof.at(0);
    //    if(region.test<80) continue;
    //    if(region.time_since_last<3.8) continue;
    //    if(region.mean_time<0.06 || region.mean_time>3.5) continue;
    //    if(region.mean_time>4) continue;
    
    double npe = region.integral*-1;
    double fp = region.fparameters[1];
    if(region.min<100){
      if(!(ch_att->baseline.found_baseline) || ch_att->tof.size()<1 || ch_att->tof[0].fparameters.size()<2){
	//      std::cout<<"failed to find regionss in event "<<evt->event_id<<std::endl;
	continue;
      }
      Pulse region_att = ch_att->tof.at(0);
      npe = npe*(1.-fp)-region_att.integral*region_att.fparameters[1]*10.271;
      fp = -1.*region_att.integral*region_att.fparameters[1]*10.271/npe;
    }
    npe /= spe;
    
//     halpha->Fill(npe,fp);
//     continue;

//make sure the gamma cuts are strigent and no alphas mistagged
    if(IsGamma(npe, fp) && npe>600){
      if(eventtimes.size()>0){
	regions.clear();
	npes.clear();
	fps.clear();
	eventtimes.clear();
	timesince.clear();
	event_ids.clear();
      }//end if size>0
      if(npe>800){
	regions.push_back(region);
	npes.push_back(npe);
	fps.push_back(fp);
	eventtimes.push_back(evt->event_time);
	timesince.push_back(evt->event_time-bad_time);
	event_ids.push_back(evt->event_id);
      }//end if npe
      continue;
    }//end if gamma

    if(eventtimes.size()==0) continue;
    regions.push_back(region);
    npes.push_back(npe);
    fps.push_back(fp);
    eventtimes.push_back(evt->event_time);
    timesince.push_back(evt->event_time-bad_time);
    event_ids.push_back(evt->event_id);
    
  }
  
  if(evt) delete evt;

  TCanvas * crun=new TCanvas("ccoin", "Coincidence Plot", 1000,700);
  crun->Divide(3,2);
  crun->cd(1);
  //  crun->cd(1)->SetLogz();
  hgamma->Draw("colz");
  //  crun->cd(2)->SetLogz();
  crun->cd(2);
  halpha->Draw("colz");
  crun->cd(3)->SetLogy();
  hdelay->Draw();
  crun->cd(4)->SetLogy();
  ht->Draw();
  crun->cd(5);
  hmt->Draw("colz");
  crun->cd(6);
  hmtfp->Draw("colz");
//   crun->cd(4)->SetLogy();
//   hmt->Draw();//ht

  return ;
}


