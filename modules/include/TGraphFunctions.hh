#ifndef TGRAPHFUNCTIONS_h
#define TGRAPHFUNCTIONS_h

// std includes
#include <algorithm>
#include <iterator>
#include <iomanip>

// ROOT includes
#include "TF1.h"
#include "TString.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TStyle.h"

using namespace std;

void printGraph(TGraphErrors* aw) {
  int an = aw->GetN();
  double* ax = aw->GetX();
  double* ay = aw->GetY();
  cout<<fixed;
  for (int i=0; i<an; i++) {
    cout<<setprecision(9)<<ax[i]<<",";
  }
  cout<<endl;
  for (int i=0; i<an; i++) {
    cout<<setprecision(9)<<ay[i]<<",";
  }
  cout<<endl;
  for (int i=0; i<an; i++) {
    cout<<setprecision(9)<<"{"<<ax[i]<<","<<ay[i]<<"}"<<",";
  }
  cout<<endl;
  return;
}

TGraphErrors* addGraph(TGraphErrors* aw1, TGraphErrors* aw2) {
  // Get values
  int an1 = aw1->GetN();
  double* ax1 = aw1->GetX();
  double* ay1 = aw1->GetY();
  double* aex1 = aw1->GetEX();
  double* aey1 = aw1->GetEY();
  int an2 = aw2->GetN();
  double* ax2 = aw2->GetX();
  double* ay2 = aw2->GetY();
  double* aex2 = aw2->GetEX();
  double* aey2 = aw2->GetEY();  

  if (an1 != an2) {
    cout<<"ERROR (addGraph): Graphs are not the same size."<<endl;
    return 0;
  }
  
  double* ax3 = new double[an1];
  double* ay3 = new double[an1];
  double* aex3 = new double[an1];
  double* aey3 = new double[an1];
  
  for (int i=0; i<an1; i++) {
    if (fabs(ax1[i]-ax2[i])>0.001) {
      cout<<"ERROR (addGraph): Graphs don't have the same axis."<<endl;
      return 0;
    }
    ax3[i] = ax1[i];
    ay3[i] = ay1[i] + ay2[i];
    aex3[i] = sqrt(aex1[i]*aex1[i] + aex2[i]*aex2[i]);
    aey3[i] = sqrt(aey1[i]*aey1[i] + aey2[i]*aey2[i]);
  }
  return new TGraphErrors(an1,ax3,ay3,aex3,aey3);
}

TGraphErrors* normGraph(TGraphErrors* aw) {
  int an = aw->GetN();
  double* ay = aw->GetY();
  //  double* aey = aw->GetEY();
  double sum = 0.;
  for (int i=0;i<an; i++) {
    sum += ay[i];
  }
  double* ny = new double[an];
  double* ney = new double[an];
  for (int i=0;i<an; i++) {
    ney[i] = sqrt(ay[i]) / sum;
    ny[i] = ay[i] / sum;
    //    ney[i] = aey[i] / sum;
  }
  return new TGraphErrors(an,aw->GetX(),ny,aw->GetEX(),ney);
}

TGraphErrors* cutGraph(TGraphErrors* aw, double startTime, double endTime, bool reZero = 0) {
  // Get values
  int an = aw->GetN();
  double* ax = aw->GetX();
  double* ay = aw->GetY();
  double* aex = aw->GetEX();
  double* aey = aw->GetEY();
  int startIndex = 0;
  int endIndex = an;
  for (int i=0; i<an; i++) {
    if (ax[i]>startTime) {
      startIndex = i;
      break;
    }
  }
  for (int i=an-1; i>=0; i--) {
    if (ax[i]<endTime) {
      endIndex = i;
      break;
    }
  }
  const int cn = endIndex-startIndex+1;
  double* cx = new double[cn];
  double* cy = new double[cn];
  double* cex = new double[cn];
  double* cey = new double[cn];
  for (int i=0; i<cn; i++) {
    cx[i]=ax[i+startIndex]-reZero*ax[startIndex];
    cy[i]=ay[i+startIndex];
    cex[i]=aex[i+startIndex];
    cey[i]=aey[i+startIndex];
  }
  TGraphErrors *tge = new TGraphErrors(cn,cx,cy,cex,cey);
  delete[] cx;
  delete[] cy;
  delete[] cex;
  delete[] cey;
  return tge;
}

double halfMaxTime(TGraphErrors* aw) {
  int an = aw->GetN();
  double* ax = aw->GetX();
  double* ay = aw->GetY();
  const double maxy = TMath::MaxElement(an,ay);
  const double halfmax = maxy/2;
  for (int i=0; i<an; i++) {
    if (ay[i]>halfmax) {
      return ax[i];
    }
  }
  return 0;
}

TGraphErrors* preTriggerGraph(TGraphErrors* aw, double halfmaxtime) {
  return cutGraph(aw, TMath::MinElement(aw->GetN(),aw->GetX()), halfMaxTime(aw)-halfmaxtime);
}


TGraphErrors* postTriggerGraph(TGraphErrors* aw, double halfmaxtime, double length=-1.) {
  if (length>0) return cutGraph(aw, halfMaxTime(aw)-halfmaxtime, halfMaxTime(aw)-halfmaxtime+length, true);
  else return cutGraph(aw, halfMaxTime(aw)-halfmaxtime, TMath::MaxElement(aw->GetN(),aw->GetX()), true);
}

double integrateGraph(TGraphErrors* aw, double x1, double x2) {
  // Get values
  int an = aw->GetN();
  double* ax = aw->GetX();
  double* ay = aw->GetY();
  double sum(0.);
  for (int i=0; i<an-1; i++) {
    double dx = ax[i+1]-ax[i];
    if (ax[i]>=x1 && ax[i]<=x2)
      sum+=ay[i]*dx;
  }
  return sum;
}

// Rebin TGraph
TGraphErrors* linearCombineGraph(TGraphErrors* aw, const int nBins=100) {
  // Get values
  int an = aw->GetN();
  double* ax = aw->GetX();
  double* ay = aw->GetY();
  //  double* aex = aw->GetEX();
  double* aey = aw->GetEY();
  // Align peak to 0
  //  double peakXvalue = 0;
  double peakYvalue = 0;
  for (int i=0;i<an;i++) {
    if (ay[i]>peakYvalue) {
      //      peakXvalue = ax[i];
      peakYvalue = ay[i];
    }
  }
  for (int i=0;i<an;i++) {
    //    ax[i] -= peakXvalue; // comment here to disable
  }
  // Set parameters
  const double startT = TMath::MinElement(an,ax);
  const double endT = TMath::MaxElement(an,ax);
  // Build axis
  TH1F* h = new TH1F("h","h",nBins,startT,endT);
  double* binCenters = new double[nBins];
  for (int i=0;i<nBins;i++) {
    binCenters[i] = h->GetBinCenter(i);
  }
  int* timesFilled = new int[nBins];
  // Combined values
  double* aX = binCenters;
  double* aY = new double[nBins];
  double* aeX = new double[nBins];
  double* aeY = new double[nBins];
  for (int i=0;i<nBins;i++) {
    timesFilled[i]=0;
    aY[i]=0.;
    aeX[i]=0.;
    aeY[i]=0.;
  }
  for (int i=0;i<an;i++) {
    double x  = ax[i];
    double y  = ay[i];
    //    double ex = aex[i];
    double ey = aey[i];
    int bin = h->FindBin(x);
    if (bin>=nBins) continue;
    timesFilled[bin]++;
    aY[bin]+=y;
    aeY[bin]+=ey*ey;
  }
  for (int i=0;i<nBins;i++) {
    if (timesFilled[i]<=0) continue;
    aY[i]/=timesFilled[i];
    aeY[i]=sqrt(aeY[i])/timesFilled[i];
  }
  // Make 0-suppressed graph
  int N=0;
  for (int i=0;i<nBins;i++)
    if (timesFilled[i]>0) N++;
  const int nBinsNew = N;
  double* aXnew = new double[nBinsNew];
  double* aYnew = new double[nBinsNew];
  double* aeXnew = new double[nBinsNew];
  double* aeYnew = new double[nBinsNew];
  int I=0;
  for (int i=0;i<nBins;i++) {
    if (timesFilled[i]>0) {
      aXnew[I] = aX[i];
      aYnew[I] = aY[i];
      aeXnew[I] = aeX[i];
      aeYnew[I] = aeY[i];
      //      aeYnew[I] = 0.; // No errors
      I++;
    }
  }
  delete[] binCenters;
  delete[] timesFilled;
  delete h;
  delete[] aY;
  delete[] aeX;
  delete[] aeY;
  return new TGraphErrors(nBinsNew,aXnew,aYnew,aeXnew,aeYnew);
}

TGraphErrors* logCombineGraph(TGraphErrors* aw, const int nBins = 200) {
  // Get values
  int an = aw->GetN();
  double* ax = aw->GetX();
  double* ay = aw->GetY();
  //  double* aex = aw->GetEX();
  double* aey = aw->GetEY();
  // Align peak to 0
  //  double peakXvalue = 0;
  double peakYvalue = 0;
  for (int i=0;i<an;i++) {
    if (ay[i]>peakYvalue) {
      //      peakXvalue = ax[i];
      peakYvalue = ay[i];
    }
  }
  for (int i=0;i<an;i++) {
    //    ax[i] -= peakXvalue; // comment here to disable
  }
  // Set parameters
  const double startT = 1e-4;
  const double endT = TMath::MaxElement(an,ax);
  // Build axis
  TH1F* h = new TH1F("h","h",nBins,log10(startT),log10(endT));
  double* binEdges = new double[nBins];
  for (int i=0;i<nBins;i++) {
    binEdges[i] = TMath::Power(10,h->GetBinCenter(i));
  }
  int* timesFilled = new int[nBins];
  // Combined values
  double* aX = new double[nBins];
  double* aY = new double[nBins];
  double* aeX = new double[nBins];
  double* aeY = new double[nBins];
  for (int i=0;i<nBins;i++) {
    timesFilled[i]=0;
    aX[i]=0.;
    aY[i]=0.;
    aeX[i]=0.;
    aeY[i]=0.;
  }
  for (int i=0;i<an;i++) {
    double x  = ax[i];
    double y  = ay[i];
    //    double ex = aex[i];
    double ey = aey[i];
    if (x<=0) continue;
    //    if (x>0.35&&x<1.5) continue;
    int bin = h->FindBin(log10(x));
    if (bin>=nBins) continue;
    timesFilled[bin]++;
    aX[bin]+=x;
    aY[bin]+=y;
    aeY[bin]+=ey*ey;
  }
  for (int i=0;i<nBins;i++) {
    if (timesFilled[i]<=0) continue;
    aX[i]/=timesFilled[i];
    aY[i]/=timesFilled[i];
    aeY[i]=sqrt(aeY[i])/pow(timesFilled[i],1.);
  }
  // Make 0-suppressed graph
  int N=0;
  for (int i=0;i<nBins;i++)
    if (timesFilled[i]>0) N++;
  const int nBinsNew = N;
  double* aXnew = new double[nBinsNew];
  double* aYnew = new double[nBinsNew];
  double* aeXnew = new double[nBinsNew];
  double* aeYnew = new double[nBinsNew];
  int I=0;
  for (int i=0;i<nBins;i++) {
    if (timesFilled[i]>0) {
      aXnew[I] = aX[i];
      aYnew[I] = aY[i];
      aeXnew[I] = aeX[i];
      aeYnew[I] = aeY[i];
      //      if (aXnew[I]>500) aeYnew[I]/=100;
      //      aeYnew[I] = 0.; // No errors
      I++;
    }
  }
  delete[] binEdges;
  delete[] timesFilled;
  delete h;
  delete[] aX;
  delete[] aY;
  delete[] aeX;
  delete[] aeY;
  TGraphErrors* tge = new TGraphErrors(nBinsNew,aXnew,aYnew,aeXnew,aeYnew);
  delete[] aXnew;
  delete[] aYnew;
  delete[] aeXnew;
  delete[] aeYnew;
  //  printGraph(tge);
  return tge;
}

/* TGraphErrors* logCombineGraph2(TGraphErrors* aw, const int nBins = 200) { */
/*   // Get values */
/*   int an = aw->GetN(); */
/*   double* ax = aw->GetX(); */
/*   double* ay = aw->GetY(); */
/*   //  double* aex = aw->GetEX(); */
/*   double* aey = aw->GetEY(); */

/*   double* aX[an]; */
/*   double* aY[an]; */
/*   double* aeX[an]; */
/*   double* aeY[an]; */
/*   for (int i=0; i<an; ) { */
/*     aX[i]=0; aY[i]=0; aeX[i]=0; aeY[i]=0; */
/*   } */
  
/*   for (int i=0; i<an; ) { */
/*     double x  = ax[i]; */
/*     double y  = ay[i]; */
/*     //    double ex = aex[i]; */
/*     double ey = aey[i]; */
/*     int nConsecutive = x; */
/*     if (nConsecutive<1) nConsecutive=1; */
    

/*   } */

/* } */

TGraph* residualGraph(TGraphErrors* aw, TF1* fn) {
  // Get values
  int an = aw->GetN();
  double* ax = aw->GetX();
  double* ay = aw->GetY();
  double xmin, xmax;
  fn->GetRange(xmin,xmax);
  double res[an];
  for (int i=0; i<an; i++) {
    double x = ax[i];
    if (x<xmin || x>xmax) res[i]=0;
    //    else res[i] = fabs(ay[i]-fn->Eval(x));
    else res[i] = (ay[i]-fn->Eval(x))/ay[i];
  }
  TGraph* tg = new TGraph(an,ax,res);
  return tg;
}

#endif
