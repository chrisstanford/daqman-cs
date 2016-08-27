#include "utilities.hh"
#include "EventData.hh"

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TSeqCollection.h"
#include "TMath.h"
#include "TPad.h"
#include "TBox.h"
#include "TPave.h"
#include "TTimeStamp.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TClass.h"
#include "TList.h"
#include "TClassMenuItem.h"

  

#include "TGFileDialog.h"

#include "Reader.hh"
#include "EventHandler.hh"
#include "ConvertData.hh"

#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;

TFile* OpenFile(const char* filename, const char* option)
{
  //see if it's already open
  TSeqCollection* files = gROOT->GetListOfFiles();
  for(int i=0; i< files->GetSize(); i++){
    if( TString(files->At(i)->GetName()) == TString(filename))
      return (TFile*)(files->At(i));
  }
  
  //otherwise try to open it
  TFile* f = new TFile(filename,option);
  if(!f || !f->IsOpen() || f->IsZombie()){
    cerr<<"Unable to open file "<<filename<<endl;
    return 0;
  }
  return f;
}

TTree* GetEventsTree(const char* filename)
{
  TFile* f = OpenFile(filename);
  if(!f)
    return 0;
  TTree* Events = (TTree*)(f->Get("Events"));
  if(!Events){
    cerr<<"Unable to open Events tree from file "<<filename<<endl;
    return 0;
  }
  return Events;
}
  
TGraph* GetAverageGraph(const char* filename, int channel)
{
  TFile* f = OpenFile(filename);
  if(!f)
    return 0;
  TString name = "average_channel";
  name += channel;
  TGraph* g = (TGraph*)(f->Get(name));
  if(!g){
    cerr<<"Unable to open graph "<<name<<" from file "<<filename<<endl;
    return 0;
  }
  return g;
}

TGraphErrors* GetAverageGraph(int run_no, int channel) {
  std::string run;
  std::stringstream out;
  out << run_no;
  run = out.str();
  std::string filename(6-run.size(), '0');
  filename = "/data/s1waveform/Run" + filename + run + ".root";
  //std::cerr << filename << std::endl;
  return (TGraphErrors*) GetAverageGraph(filename.c_str(), channel);
}

int DividePad(TPad* p, int nplots)
{
  int npadsx = 1, npadsy=1;
  if( nplots < 2)
    {}
  else if( nplots == 2)
    npadsx=2;
  else if (nplots < 5){
    npadsx=2; npadsy=2;
  }
  else if(nplots < 7){
    npadsx=3; npadsy=2;
  }
  else if(nplots < 10){
    npadsx=3; npadsy=3;
  }
  else if(nplots < 13){
    npadsx=4, npadsy=3;
  }
  else if(nplots < 17){
    npadsx=4; npadsy=4;
  }
  else{
    npadsx=(int)TMath::Ceil(sqrt(nplots)); 
    npadsy=(int)TMath::Ceil(sqrt(nplots));
  }
  p->Divide(npadsx,npadsy);
  return npadsx*npadsy;
}

void DrawOperationsBoxes(bool drawbubble, bool drawrecirc)
{
  double dummy, y1, y2;
  gPad->Update();
  gPad->GetRangeAxis(dummy,y1,dummy,y2);
  
  if(drawbubble){
    double x1[] = { TTimeStamp(2010,04,15,11,06,55,0,0).GetSec(),
		    TTimeStamp(2010,05,03,15,51,33,0,0).GetSec(),
		    TTimeStamp(2010,05,05,10,54,19,0,0).GetSec(),
		    TTimeStamp(2010,05,05,15,57,33,0,0).GetSec()
    };
    double x2[] = { TTimeStamp(2010,04,15,16,38,48,0,0).GetSec(),
		    TTimeStamp(2010,05,03,20,40,28,0,0).GetSec(),
		    TTimeStamp(2010,05,05,13,29,44,0,0).GetSec(),
		    TTimeStamp(2010,05,05,18,45,53,0,0).GetSec()
    };
    
    for(size_t i = 0; i < sizeof(x1)/sizeof(double); i++){
      TBox b(x1[i], y1, x2[i], y2);
      b.SetFillStyle(3002);
      b.SetFillColor(kRed);
      b.DrawClone();
    }
  }
  if(drawrecirc){
    double x1[] = { TTimeStamp(2010,05,07,15,40,32,0,0).GetSec(),
		    TTimeStamp(2010,07,14,17,23,00,0,0).GetSec()
    };
    double x2[] = { TTimeStamp(2010,05,07,23,59,59,0,0).GetSec(),
		    TTimeStamp(2010,07,15,16,27,21,0,0).GetSec()
    };
    for(size_t i = 0; i < sizeof(x1)/sizeof(double); i++){
      TBox b(x1[i], y1, x2[i], y2);
      b.SetFillStyle(3002);
      b.SetFillColor(kBlue);
      b.DrawClone();
    }
  }
}

double GetBaseline(TGraph* g, int npts, bool subtract)
{
  double* gy = g->GetY();
  double baseline = accumulate(gy,gy+npts,0.)/(1.*npts);
  if(subtract){
    for(int i=0; i < g->GetN(); i++) gy[i] -= baseline;
  }
  return baseline;
}

TGraph* GetRealEvent(const char* filename, int eventnum, int channel)
{
  Reader reader(filename);
  if(!reader.IsOk()){
    std::cerr<<"Unable to open file "<<filename<<std::endl;
    return 0;
  }
  
  EventHandler* handler = EventHandler::GetInstance();
  handler->AddModule<ConvertData>();
  handler->Initialize();
  RawEventPtr event = reader.GetEventWithID(eventnum);
  if(event == RawEventPtr()){
    std::cerr<<"Unable to load event with id "<<eventnum<<std::endl;
    return 0;
  }
  handler->Process(event);
  EventPtr evt = handler->GetCurrentEvent();
  handler->Finalize();
  
  ChannelData* chdata = evt->GetEventData()->GetChannelByID(channel);
  if(!chdata){
    std::cerr<<"Unable to load data for channel "<<channel<<std::endl;
    return 0;
  }
  return chdata->GetTGraph();
}

ChannelData* GetChannelData(const char* filename, int eventnum, int channel)
{
  Reader reader(filename);
  if(!reader.IsOk()){
    std::cerr<<"Unable to open file "<<filename<<std::endl;
    return 0;
  }

  EventHandler* handler = EventHandler::GetInstance();
  handler->AddModule<ConvertData>();
  handler->Initialize();

  ChannelData* chdata = 0;
  
  while(reader.IsOk() && !reader.eof()){
    RawEventPtr event = reader.GetEventWithID(eventnum++);
    if(event == RawEventPtr()){
      std::cerr<<"Unable to load event with id "<<eventnum<<std::endl;
      return 0;
    }
    if(eventnum%5000 == 0)
      Message(INFO)<<"Processing event "<<eventnum<<std::endl;
    
    handler->Process(event);
  }
  EventPtr evt = handler->GetCurrentEvent();
  handler->Finalize();
  
  chdata = evt->GetEventData()->GetChannelByID(channel);
  if(!chdata){
    std::cerr<<"Unable to load data for channel "<<channel<<std::endl;
    return 0;
  }
  
  return chdata;
}

double CorrelationCoefficient(int npts, double* x, double* y)
{
  double sumx=0, sumy=0, sumx2=0, sumy2=0, sumxy=0;
  for(int i=0; i< npts; i++){
    sumx += x[i];
    sumy += y[i];
    sumx2 += x[i]*x[i];
    sumy2 += y[i]*y[i];
    sumxy += x[i]*y[i];
  }
  //formula from wikipedia Correlation and Dependence
  return (npts * sumxy - sumx*sumy) / 
    ( sqrt(npts*sumx2 - sumx*sumx) * sqrt(npts*sumy2 - sumy*sumy) );
}

void SaveHistoToFile(TObject* c)
{
  if(!c->InheritsFrom("TH1"))
    return;
  TH1* h = (TH1*)c;
  static const char* filetypes[] = {
    "Text files","*.txt",
    0,0
  };
  TGFileInfo fi;
  fi.fFileTypes = filetypes;
  new TGFileDialog(gClient->GetRoot(),0,kFDSave,&fi);
  if(!fi.fFilename)
    return;
  std::string fname = fi.fFilename;
  if(fname.rfind(".txt")==std::string::npos)
    fname.append(".txt");
  std::cout<<"Saving histogram "<<h->GetName()<<" to file "<<fname<<std::endl;
  std::ofstream fout(fname.c_str());
  if(fout.is_open()){
    fout<<"#X\tY\tYerr\n";
    for(int i=1; i<=h->GetNbinsX(); ++i){
      fout<<h->GetBinLowEdge(i)<<"\t"
	  <<h->GetBinContent(i)<<"\t"
	  <<h->GetBinError(i)
	  <<endl;
    }
  }
}

void CustomizeHistogramMenus()
{
  //have to do each child class separately
  //TH1F
  TClass* cl = TH1F::Class();
  TList* l = cl->GetMenuList();
  l->AddFirst(new TClassMenuItem(TClassMenuItem::kPopupUserFunction,cl,
				 "Save as ASCII","SaveHistoToFile",0,
				 "TObject*",0));
  //TH1D
  cl = TH1D::Class();
  l = cl->GetMenuList();
  l->AddFirst(new TClassMenuItem(TClassMenuItem::kPopupUserFunction,cl,
				 "Save as ASCII","SaveHistoToFile",0,
				 "TObject*",0));
  

}


/// Fit the SPE spectrum
void FitChannelSPE(int ch_id, double fit_min, double fit_max, 
		   double max_baseline_sigma, double min_spe_amp, 
		   bool single_gaus){
  
  TTree *Events = (TTree*)(gDirectory->Get("Events"));  
  if(!Events || Events->GetEntries()<1){
    std::cout<<"No data available!"<<std::endl;
    return;
  }
  
  EventData* evt = 0;
  Events->SetBranchAddress("event",&evt);
  
  int run_id = Events->GetMaximum("run_id");
  char chname[30];
  sprintf(chname, "Run%d_Ch%d_SPE",run_id,ch_id);
  TH1F *hspe = (TH1F*)gDirectory->Get(chname);
  if(hspe) delete hspe;
  hspe = new TH1F(chname, chname, 200,-50,550);
  hspe->SetDirectory(0);
  sprintf(chname, "Run%d_Ch%d_SPE_plot",run_id,ch_id);
  //  sprintf(chname, "Channel_%d_SPE_plot",ch_id);
  TCanvas *cspe = (TCanvas *)(gDirectory->Get(chname));
  if(cspe) cspe->cd();
  else cspe=new TCanvas(chname, chname, 800, 600);
  cspe->cd();  
  //  cspe->SetLogy(true);

  ChannelData* ch_signal = NULL;

  for(int entry=0; entry<Events->GetEntries(); entry++){
  //  for(int entry=0; entry<1000; entry++){
    Events->GetEntry(entry);
    //look at the signal channel
    ch_signal = evt->GetChannelByID(ch_id);
    if(!ch_signal){
      std::cout<<ch_id<<" channel does not exist in the channel data!"<<std::endl;
      break;
    }
    if(ch_signal->single_pe.size()<1){
      //      std::cout<<"failed to find correct signal pulses in event "<<evt->event_id<<std::endl;
      continue;
    }
    for(size_t spe_id=0; spe_id<ch_signal->single_pe.size(); spe_id++){
      if(ch_signal->single_pe.at(spe_id).baseline_sigma>max_baseline_sigma ||
	 ch_signal->single_pe.at(spe_id).amplitude<min_spe_amp) continue;
      hspe->Fill(-ch_signal->single_pe.at(spe_id).integral);
    }
  }//end for loop

  hspe->Draw();

  int spe_bin = hspe->GetMaximumBin();
  double spe_peak = hspe->GetBinContent(spe_bin);
  double ped_peak = spe_peak;
  
  TF1 *fspe=(TF1*)(gDirectory->Get("fspe"));
  if(fspe) fspe->SetRange(fit_min, fit_max);
  else fspe=new TF1("fspe", "gaus(0)+gaus(3)", fit_min, fit_max);
  fspe->SetParameters(ped_peak,0,4,spe_peak,200,50);
  fspe->FixParameter(1,0);
  if(single_gaus){
    fspe->FixParameter(0,0);
    fspe->FixParameter(2,1);
  }

  std::cout<<"******** Fitting data: "<<std::endl;
  int fitresult = hspe->Fit("fspe", "REM");
  std::cout<<"******** The Fit result is: "<<fitresult<<std::endl
	   <<"******** The Fit Chi2/NDF is: "<<fspe->GetChisquare()<<'/'<<fspe->GetNDF()<<std::endl
	   <<"******** The Probability of the fit is: "<<TMath::Prob(fspe->GetChisquare(), fspe->GetNDF())
	   <<std::endl<<std::endl;
  ped_peak = fspe->GetParameter(0);
  spe_peak = fspe->GetParameter(3);
  fspe->SetParameter(0,0);
  fspe->SetLineColor(3);
  fspe->DrawCopy("same");
  fspe->SetParameter(0,ped_peak);
  fspe->SetParameter(3,0);
  fspe->SetLineColor(4);
  fspe->DrawCopy("same");
  cspe->Modified();
  cspe->Update();

  while(1){

    std::cout<<"Would you like to change the fit range and redo the fit? y/[n]"<<std::endl;
    std::string result;
    getline(std::cin, result);
    if (result == "" || (result[0]!='y'&&result[0]!='Y')) break;

    std::cout<<"Please enter the new lower bound for SPE fit:"<<std::endl;
    getline(std::cin, result);
    fit_min = atof(result.c_str());
    std::cout<<"Please enter the new upper bound for SPE fit:"<<std::endl;
    getline(std::cin, result);
    fit_max = atof(result.c_str());
    
    if(fit_max<=fit_min || fit_min<hspe->GetBinLowEdge(1) ||
       fit_max>hspe->GetBinLowEdge(hspe->GetNbinsX())+hspe->GetBinWidth(1)){
      std::cout<<"Illegal fit range!"<<std::endl;
    }

    cspe->cd();  
    //    cspe->SetLogy(true);
    hspe->Draw();
    ped_peak = hspe->GetBinContent(hspe->GetMaximumBin());
    spe_peak = hspe->GetBinContent(hspe->FindBin(40));
    fspe->SetParameters(ped_peak,0,4,spe_peak,40,20);
    fspe->SetRange(fit_min, fit_max);
    fspe->SetLineColor(2);
    fitresult = hspe->Fit("fspe", "REM");
    std::cout<<"******** The Fit result is: "<<fitresult<<std::endl
	     <<"******** The Fit Chi2/NDF is: "<<fspe->GetChisquare()<<'/'<<fspe->GetNDF()<<std::endl
	     <<"******** The Probability of the fit is: "<<TMath::Prob(fspe->GetChisquare(), fspe->GetNDF())
	     <<std::endl<<std::endl;
    ped_peak = fspe->GetParameter(0);
    spe_peak = fspe->GetParameter(3);
    fspe->SetParameter(0,0);
    fspe->SetLineColor(3);
    fspe->DrawCopy("same");
    fspe->SetParameter(0,ped_peak);
    fspe->SetParameter(3,0);
    fspe->SetLineColor(4);
    fspe->DrawCopy("same");
    cspe->Modified();
    cspe->Update();

  }

  std::cout<<"Final fit result: ( "<<ch_id<<" , "<<fspe->GetParameter(4)<<" ) "<<std::endl;

  return;

}
