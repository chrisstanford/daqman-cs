/** @file utilities.hh
    @brief Useful functions for loading trees/files/etc within root
    @author bloer
    @ingroup daqroot
*/

#ifndef UTILITIES_h
#define UTILITIES_h

//have the makefile add these functions to the root lib:
//ClassDef

class TFile;
class TTree;
class TGraph;
class TGraphErrors;
class TPad;
class TChain;
class TObject;

#include "TCanvas.h"
#include "TCut.h"
#include "ChannelData.hh"

#include <string>
#include <vector>
#include <fstream>
#include <sstream>

/// Open a file, and make sure it actually did open properly
TFile* OpenFile(const char* filename, const char* option = "READ");
/// Load the 'Events' tree from a given file
TTree* GetEventsTree(const char* filename);
/// Load the 'average_channel#' graph from a given file for channel n
TGraph* GetAverageGraph(const char* filename, int channel);
/// Load the 'average_channel#' graph from the official file for run run_no and channel n
TGraphErrors* GetAverageGraph(int run_no, int channel); 
/// Get the baseline of a TGraph by averagint the first npts points
double GetBaseline(TGraph* g, int npts=100, bool subtract = false);
/// Decide how to divde a pad based on the number of plots
int DividePad(TPad* p, int nplots);

/// Draw TBoxes on a graph where we did mechanical operations on the detector
void DrawOperationsBoxes(bool drawbubble=true, bool drawrecirc=true);

/// Get the response for the given channel and event in file as a TGraph
TGraph* GetRealEvent(const char* filename, int eventnum, int channel);

/// Get the pointer to ChannelData the the given event and channel
ChannelData* GetChannelData(const char* filename, int eventnum, int channel);

/// Calculate the correlation coefficient between two vectors
double CorrelationCoefficient(int npts, double* x, double* y);

/// Print a TH1 bin contents to a text file
void SaveHistoToFile(TObject* c);
/// Add menus to TH1{F,D}
void CustomizeHistogramMenus();

/// Fit the SPE spectrum
void FitChannelSPE(int ch_id, double fit_min, double fit_max, 
		   double max_baseline_sigma=100, double min_spe_amp=3, 
		   bool single_gaus=true);

//test a file exists or not w/o trying to open it
bool FileExists(const char* filename, bool print=false);

template <typename T>
inline size_t LoadArray(size_t N, T ary[], const char *fname,
			size_t lstart, size_t column){
  if(!FileExists(fname)) return 0;
  size_t count=0, lnum=0, content=0, space=0;
  bool accept=false;
  string line;
  ifstream infile(fname);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile, line);
      if(++lnum<lstart) continue;
      accept=false;
      content=space=0;
      for(size_t i=0; i<column; i++){
        content=line.find_first_not_of(" \t", space);
        if(content==string::npos) break;
        space=line.find_first_of(" \t", content);
        if(space==string::npos) space=line.length();
        if(i==column-1&&space>content) accept=true;
      }
      if(!accept) continue;
      if(count>=N) {
        std::cout<<"More data than wanted, the first "<<N<<" are read!"<<std::endl;
        break;
      }
      //this may need to be improved.
      //may want to use  boost/lexical_cast if we can make sure boost is available 
      T value;
      istringstream istr(line.substr(content,space-content));
      if(!(istr>>value)) value = 0;
      ary[count] = value;
      //      std::cout<<count<<'\t'<<ary[count]<<std::endl;
      count ++;
    }//endwhile()
    infile.close();
  }//endif()
  
  if(count<N){
    //    std::cout<<"Less data acqusited! "<<count<<std::endl;
    for(size_t i=count;i<N;i++) ary[i]=0;
  }
  std::cout<<fname<<" Loaded!"<<std::endl;
  return count;
}

#endif
