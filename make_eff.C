#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TVirtualFitter.h"
#include "TGraph.h"

using std::cout;
using std::endl;

//using namespace std;

float get_nevents(TH1F *hst, float bin_lo, float bin_hi);
//float get_nevents_err(TH1F *hst, float bin_lo, float bin_hi);
void decorate(TGraph *h, TString gtitle, float gmax, float gmin, 
	      int markercolor, int markerstyle, int linecolor, int linewidth);

void make_eff(int opt1=0)
  
{
  //opt1 can be used to add options (add opt2 if necessary)
  TString filename = "output/hst_DY.root";//Add filename with histograms

  //Add the string name according to what you have given
  // We will construct denominator and numerator names by doing
  // den_ + plotname, num_ + plotname
  // So if plotname=effpt, then the histograms in the hst file
  // should be den_effpt and num_effpt
  TString plot0name = "mupt";
  TString plot1name = "ept";
  TString plot2name = "mueta";
  TString plot3name = "eeta";

  for(int i=0; i<4; i++){
    opt1=i;
    if(opt1==0){
      // Determine the bin-edges; efficiency will be calculated in these bins
      float lowed[] =  {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
      float hied[]  = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,200};
      int nbin = 16; 
      
      //Now open the file
      TFile *f1 = new TFile(filename);
      
      //Open the histograms
      TString denname = "den_" + plot0name;
      TString numname = "num_" + plot0name;
      TH1F *h0 = (TH1F*)f1->Get(denname);
      TH1F *h1 = (TH1F*)f1->Get(numname); //h0 is the denominator, h1 is the numerator
      
      //Declare some arrays to store the efficiency
      int n = 0;
      float x[150],y[150],exl[150],exh[150],eyl[150],eyh[150], ex[150];
      for(int i=0; i<150; i++){
	x[i]=y[i]=0; exl[i]=exh[i]=eyl[i]=eyh[i]=ex[i]=0;
      }
      float eff1[3];
      
      // Loop over the nbins and calculate the efficiency in each bin
      // and then store it in the appropriate arrays.
      float et_begin =0;
      float et_end = 0;
      float et_step = 0;
      float et = et_begin;
      double et_lo,et_hi;
      for(int i=0; i<nbin; i++){
	et_lo = lowed[i]; et_hi = hied[i];
	et_step = et_hi - et_lo;
	
	eff1[0]=eff1[1]=eff1[2]=0;
	
	// Next 4 lines get the events in the bin
	// from denominator and numerator histograms
	float den = get_nevents(h0,et_lo,et_hi);
	//	float dene = get_nevents_err(h0,et_lo,et_hi);
	float num = get_nevents(h1,et_lo,et_hi);
	//	float nume = get_nevents_err(h1,et_lo,et_hi);
	
	x[n] = (et_lo+et_hi)/2;
	exh[n] = (et_hi - et_lo)/2.;       exl[n] = (et_hi - et_lo)/2;
	
	// Uncomment next line for debug
	//	cout<<x[n]<<" "<<num<<" ("<<nume<<") / "<<den<<" ("<<dene<<") "<<endl;//debug
	
	if(num>0 && den>0){
	  y[n] = float(num/den);
	  
	  // Different way to calculate uncertainty
	  // float eypc = sqrt( pow(nume/num,2) + pow(dene/den,2) );
	  // float ey = eypc*y[n];
	  // Main way to calculate uncertainty
	  float ey = sqrt( y[n]*(1-y[n])/den );
	  eyl[n] = ey;
	  eyh[n] = ey;
	}
	else
	  y[n]=0.;
	n++;
      }
      
      //Print the results
      cout<<"Bin Efficiency Error"<<endl;
      for(int i=0; i<nbin; i++)
	cout<<lowed[i]<<"-"<<hied[i]<<" "<<y[i]<<" +/- "<<eyl[i]<<endl;
      
      //----------------------------------------------------------------------------------------
      
      //Now we can graph the efficiency as well
      
      // using nbin, x, ex, y, eyl, eyh.
      
      // Open a canvas
      TCanvas *c0 = new TCanvas("c0","c0",800,600);
      
      //Declare the graph
      TGraphErrors *g0 = new TGraphErrors(nbin,x,y,ex,eyl);  
      //decorate the graph
      decorate(g0,";Muon p_{T} [GeV]; Efficiency",0.1,0,kBlack,21,kBlack,2);
      g0->Draw("ap");
      
      // In case needed,
      //    this is how to draw a line
      // TLine line;line.SetLineColor(kBlue); line.SetLineWidth(3); line.SetLineStyle(7);
      // line.DrawLine(0,1,210,1);
      //    this is how to draw latex text
      // TLatex tL;  tL.SetTextColor(kBlack) ; tL.SetTextSize(0.07) ;
      // tL.DrawLatex(140,0.92,"ele: ptcone40/p_{T} < 0.1");
      
      
    }
    if(opt1==1){
      // Determine the bin-edges; efficiency will be calculated in these bins
      float lowed[] =  {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
      float hied[]  = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,200};
      int nbin = 16; 
      
      //Now open the file
      TFile *f1 = new TFile(filename);
      
      //Open the histograms
      TString denname = "den_" + plot1name;
      TString numname = "num_" + plot1name;
      TH1F *h0 = (TH1F*)f1->Get(denname);
      TH1F *h1 = (TH1F*)f1->Get(numname); //h0 is the denominator, h1 is the numerator
      
      //Declare some arrays to store the efficiency
      int n = 0;
      float x[150],y[150],exl[150],exh[150],eyl[150],eyh[150], ex[150];
      for(int i=0; i<150; i++){
	x[i]=y[i]=0; exl[i]=exh[i]=eyl[i]=eyh[i]=ex[i]=0;
      }
      float eff1[3];
      
      // Loop over the nbins and calculate the efficiency in each bin
      // and then store it in the appropriate arrays.
      float et_begin =0;
      float et_end = 0;
      float et_step = 0;
      float et = et_begin;
      double et_lo,et_hi;
      for(int i=0; i<nbin; i++){
	et_lo = lowed[i]; et_hi = hied[i];
	et_step = et_hi - et_lo;
	
	eff1[0]=eff1[1]=eff1[2]=0;
	
	// Next 4 lines get the events in the bin
	// from denominator and numerator histograms
	float den = get_nevents(h0,et_lo,et_hi);
	//	float dene = get_nevents_err(h0,et_lo,et_hi);
	float num = get_nevents(h1,et_lo,et_hi);
	//	float nume = get_nevents_err(h1,et_lo,et_hi);
	
	x[n] = (et_lo+et_hi)/2;
	exh[n] = (et_hi - et_lo)/2.;       exl[n] = (et_hi - et_lo)/2;
	
	// Uncomment next line for debug
	//	cout<<x[n]<<" "<<num<<" ("<<nume<<") / "<<den<<" ("<<dene<<") "<<endl;//debug
	
	if(num>0 && den>0){
	  y[n] = float(num/den);
	  
	  // Different way to calculate uncertainty
	  // float eypc = sqrt( pow(nume/num,2) + pow(dene/den,2) );
	  // float ey = eypc*y[n];
	  // Main way to calculate uncertainty
	  float ey = sqrt( y[n]*(1-y[n])/den );
	  eyl[n] = ey;
	  eyh[n] = ey;
	}
	else
	  y[n]=0.;
	n++;
      }
      
      //Print the results
      cout<<"Bin Efficiency Error"<<endl;
      for(int i=0; i<nbin; i++)
	cout<<lowed[i]<<"-"<<hied[i]<<" "<<y[i]<<" +/- "<<eyl[i]<<endl;
      
      //----------------------------------------------------------------------------------------
      
      //Now we can graph the efficiency as well
      
      // using nbin, x, ex, y, eyl, eyh.
      
      // Open a canvas
      TCanvas *c1 = new TCanvas("c1","c1",800,600);
      
      //Declare the graph
      TGraphErrors *g1 = new TGraphErrors(nbin,x,y,ex,eyl);  
      //decorate the graph
      decorate(g1,";Electron p_{T} [GeV]; Efficiency",0.1,0,kBlack,21,kBlack,2);
      
      g1->Draw("ap");
      
      // In case needed,
      //    this is how to draw a line
      // TLine line;line.SetLineColor(kBlue); line.SetLineWidth(3); line.SetLineStyle(7);
      // line.DrawLine(0,1,210,1);
      //    this is how to draw latex text
      // TLatex tL;  tL.SetTextColor(kBlack) ; tL.SetTextSize(0.07) ;
      // tL.DrawLatex(140,0.92,"ele: ptcone40/p_{T} < 0.1");    
    }
    
    if(opt1==2){
      // Determine the bin-edges; efficiency will be calculated in these bins
      float lowed[] =  {-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5};
      float hied[]  = {-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3};
      int nbin = 12; 
      
      //Now open the file
      TFile *f1 = new TFile(filename);
      
      //Open the histograms
      TString denname = "den_" + plot2name;
      TString numname = "num_" + plot2name;
      TH1F *h0 = (TH1F*)f1->Get(denname);
      TH1F *h1 = (TH1F*)f1->Get(numname); //h0 is the denominator, h1 is the numerator
      
      //Declare some arrays to store the efficiency
      int n = 0;
      float x[150],y[150],exl[150],exh[150],eyl[150],eyh[150], ex[150];
      for(int i=0; i<150; i++){
	x[i]=y[i]=0; exl[i]=exh[i]=eyl[i]=eyh[i]=ex[i]=0;
      }
      float eff1[3];
      
      // Loop over the nbins and calculate the efficiency in each bin
      // and then store it in the appropriate arrays.
      float et_begin =0;
      float et_end = 0;
      float et_step = 0;
      float et = et_begin;
      double et_lo,et_hi;
      for(int i=0; i<nbin; i++){
	et_lo = lowed[i]; et_hi = hied[i];
	et_step = et_hi - et_lo;
	
	eff1[0]=eff1[1]=eff1[2]=0;
	
	// Next 4 lines get the events in the bin
	// from denominator and numerator histograms
	float den = get_nevents(h0,et_lo,et_hi);
	//	float dene = get_nevents_err(h0,et_lo,et_hi);
	float num = get_nevents(h1,et_lo,et_hi);
	//	float nume = get_nevents_err(h1,et_lo,et_hi);
	
	x[n] = (et_lo+et_hi)/2;
	exh[n] = (et_hi - et_lo)/2.;       exl[n] = (et_hi - et_lo)/2;
	
	// Uncomment next line for debug
	//	cout<<x[n]<<" "<<num<<" ("<<nume<<") / "<<den<<" ("<<dene<<") "<<endl;//debug
	
	if(num>0 && den>0){
	  y[n] = float(num/den);
	  
	  // Different way to calculate uncertainty
	  // float eypc = sqrt( pow(nume/num,2) + pow(dene/den,2) );
	  // float ey = eypc*y[n];
	  // Main way to calculate uncertainty
	  float ey = sqrt( y[n]*(1-y[n])/den );
	  eyl[n] = ey;
	  eyh[n] = ey;
	}
	else
	  y[n]=0.;
	n++;
      }
      
      //Print the results
      cout<<"Bin Efficiency Error"<<endl;
      for(int i=0; i<nbin; i++)
	cout<<lowed[i]<<"-"<<hied[i]<<" "<<y[i]<<" +/- "<<eyl[i]<<endl;
      
      //----------------------------------------------------------------------------------------
      
      //Now we can graph the efficiency as well
      
      // using nbin, x, ex, y, eyl, eyh.
      
      // Open a canvas
      TCanvas *c2 = new TCanvas("c2","c2",800,600);
      
      //Declare the graph
      TGraphErrors *g2 = new TGraphErrors(nbin,x,y,ex,eyl);  
      //decorate the graph
      decorate(g2,";Muon eta; Efficiency",0.1,0,kBlack,21,kBlack,2);
      
      g2->Draw("ap");
      
      // In case needed,
      //    this is how to draw a line
      // TLine line;line.SetLineColor(kBlue); line.SetLineWidth(3); line.SetLineStyle(7);
      // line.DrawLine(0,1,210,1);
      //    this is how to draw latex text
      // TLatex tL;  tL.SetTextColor(kBlack) ; tL.SetTextSize(0.07) ;
      // tL.DrawLatex(140,0.92,"ele: ptcone40/p_{T} < 0.1");    
    }
    if(opt1==3){
      // Determine the bin-edges; efficiency will be calculated in these bins
      float lowed[] =  {-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5};
      float hied[]  = {-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3};
      int nbin = 12; 
      
      //Now open the file
      TFile *f1 = new TFile(filename);
      
      //Open the histograms
      TString denname = "den_" + plot3name;
      TString numname = "num_" + plot3name;
      TH1F *h0 = (TH1F*)f1->Get(denname);
      TH1F *h1 = (TH1F*)f1->Get(numname); //h0 is the denominator, h1 is the numerator
      
      //Declare some arrays to store the efficiency
      int n = 0;
      float x[150],y[150],exl[150],exh[150],eyl[150],eyh[150], ex[150];
      for(int i=0; i<150; i++){
	x[i]=y[i]=0; exl[i]=exh[i]=eyl[i]=eyh[i]=ex[i]=0;
      }
      float eff1[3];
      
      // Loop over the nbins and calculate the efficiency in each bin
      // and then store it in the appropriate arrays.
      float et_begin =0;
      float et_end = 0;
      float et_step = 0;
      float et = et_begin;
      double et_lo,et_hi;
      for(int i=0; i<nbin; i++){
	et_lo = lowed[i]; et_hi = hied[i];
	et_step = et_hi - et_lo;
	
	eff1[0]=eff1[1]=eff1[2]=0;
	
	// Next 4 lines get the events in the bin
	// from denominator and numerator histograms
	float den = get_nevents(h0,et_lo,et_hi);
	//	float dene = get_nevents_err(h0,et_lo,et_hi);
	float num = get_nevents(h1,et_lo,et_hi);
	//	float nume = get_nevents_err(h1,et_lo,et_hi);
	
	x[n] = (et_lo+et_hi)/2;
	exh[n] = (et_hi - et_lo)/2.;       exl[n] = (et_hi - et_lo)/2;
	
      // Uncomment next line for debug
	//	cout<<x[n]<<" "<<num<<" ("<<nume<<") / "<<den<<" ("<<dene<<") "<<endl;//debug
	
	if(num>0 && den>0){
	  y[n] = float(num/den);
	  
	  // Different way to calculate uncertainty
	  // float eypc = sqrt( pow(nume/num,2) + pow(dene/den,2) );
	  // float ey = eypc*y[n];
	  // Main way to calculate uncertainty
	  float ey = sqrt( y[n]*(1-y[n])/den );
	  eyl[n] = ey;
	  eyh[n] = ey;
	}
	else
	  y[n]=0.;
	n++;
      }
      
      //Print the results
      cout<<"Bin Efficiency Error"<<endl;
      for(int i=0; i<nbin; i++)
	cout<<lowed[i]<<"-"<<hied[i]<<" "<<y[i]<<" +/- "<<eyl[i]<<endl;
      
      //----------------------------------------------------------------------------------------
      
      //Now we can graph the efficiency as well
      
      // using nbin, x, ex, y, eyl, eyh.
      
      // Open a canvas
      TCanvas *c3 = new TCanvas("c3","c3",800,600);
      
      //Declare the graph
      TGraphErrors *g3 = new TGraphErrors(nbin,x,y,ex,eyl);  
      //decorate the graph
      decorate(g3,";Electron eta; Efficiency",0.1,0,kBlack,21,kBlack,2);
      g3->Draw("ap");
      
      // In case needed,
      //    this is how to draw a line
      // TLine line;line.SetLineColor(kBlue); line.SetLineWidth(3); line.SetLineStyle(7);
      // line.DrawLine(0,1,210,1);
      //    this is how to draw latex text
      // TLatex tL;  tL.SetTextColor(kBlack) ; tL.SetTextSize(0.07) ;
      // tL.DrawLatex(140,0.92,"ele: ptcone40/p_{T} < 0.1");    
    }
  }
}

float get_nevents(TH1F *hst, float bin_lo, float bin_hi)
{
    int bin_width = hst->GetBinWidth(1);
    int ibin_begin = 1;
    float nevents = 0.;
    while ( hst->GetBinCenter(ibin_begin) < bin_lo )
        ibin_begin++;
    int ibin_end = ibin_begin;
    while ( hst->GetBinCenter(ibin_end) < bin_hi )
        ibin_end++;
    for ( int i=ibin_begin; i<ibin_end; i++ )
        nevents += hst->GetBinContent(i);

    return nevents;
}
/*
float get_nevents_err(TH1F *hst, float bin_lo, float bin_hi)
{
    int bin_width = hst->GetBinWidth(1);
    int ibin_begin = 1;
    float nevents = 0.;
    while ( hst->GetBinCenter(ibin_begin) < bin_lo )
        ibin_begin++;
    int ibin_end = ibin_begin;
    while ( hst->GetBinCenter(ibin_end) < bin_hi )
        ibin_end++;
    for ( int i=ibin_begin; i<ibin_end; i++ )
      nevents += pow(hst->GetBinError(i),2);
    nevents = sqrt(nevents);
    return nevents;
}

*/
void decorate(TGraph *h, TString gtitle, float gmax, float gmin, 
	      int markercolor, int markerstyle, int linecolor, int linewidth)
{
  h->SetTitle(gtitle);
  h->SetMarkerColor(markercolor); h->SetMarkerStyle(markerstyle);
  h->SetLineColor(linecolor); h->SetLineWidth(linewidth);
  h->SetMaximum(gmax); h->SetMinimum(gmin);
}

