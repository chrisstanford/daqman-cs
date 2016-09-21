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
    double x1[] = { (double)TTimeStamp(2010,04,15,11,06,55,0,0).GetSec(),
		    (double)TTimeStamp(2010,05,03,15,51,33,0,0).GetSec(),
		    (double)TTimeStamp(2010,05,05,10,54,19,0,0).GetSec(),
		    (double)TTimeStamp(2010,05,05,15,57,33,0,0).GetSec()
    };
    double x2[] = { (double)TTimeStamp(2010,04,15,16,38,48,0,0).GetSec(),
		    (double)TTimeStamp(2010,05,03,20,40,28,0,0).GetSec(),
		    (double)TTimeStamp(2010,05,05,13,29,44,0,0).GetSec(),
		    (double)TTimeStamp(2010,05,05,18,45,53,0,0).GetSec()
    };
    
    for(size_t i = 0; i < sizeof(x1)/sizeof(double); i++){
      TBox b(x1[i], y1, x2[i], y2);
      b.SetFillStyle(3002);
      b.SetFillColor(kRed);
      b.DrawClone();
    }
  }
  if(drawrecirc){
    double x1[] = { (double)TTimeStamp(2010,05,07,15,40,32,0,0).GetSec(),
		    (double)TTimeStamp(2010,07,14,17,23,00,0,0).GetSec()
    };
    double x2[] = { (double)TTimeStamp(2010,05,07,23,59,59,0,0).GetSec(),
		    (double)TTimeStamp(2010,07,15,16,27,21,0,0).GetSec()
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
