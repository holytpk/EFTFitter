// -*- C++ -*-
// author: afiq anuar
// short: styling methods for ROOT plotting and histogram manipulation

#ifndef PLOTUTIL_H
#define PLOTUTIL_H

#include "TemplateUtil.h"

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TClass.h"
#include "TKey.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TColor.h"

#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "THStack.h"

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "TAxis.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TLatex.h"

#include "TMarker.h"
#include "TLine.h"

void setH1Style()
{
  // always reset everything first
  gStyle->Reset();
  gROOT->SetStyle("Plain");

  // if color alpha is needed (some install may not have this by default)
  gStyle->SetCanvasPreferGL(true);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);

  //gStyle->SetFrameColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetFillColor(0);

  //gStyle->SetPaperSize(20, 26);
  //gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.085);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadLeftMargin(0.06);
  //gStyle->SetCanvasDefH(800);
  //gStyle->SetCanvasDefW(800);
  //gStyle->SetPadGridX(1);
  //gStyle->SetPadGridY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.03);
  gStyle->SetLabelFont(42, "xyz");
  gStyle->SetTitleFont(42, "xyz");
  gStyle->SetLabelSize(0.025, "xyz");
  gStyle->SetTitleSize(0.027, "xyz");
  gStyle->SetTitleOffset(1.1, "y");
  gStyle->SetTitleOffset(1.1, "x");
    
  gStyle->SetTitleX(0.5); // suit the plot
  gStyle->SetTitleY(0.97);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFontSize(0.035);
  //gStyle->SetTitleStyle(1001);
  //gStyle->SetPadTopMargin(0.05);
  //gStyle->SetPadBottomMargin(0.10);
  //gStyle->SetPadLeftMargin(0.13);
  //gStyle->SetPadRightMargin(0.02);

  // use bold lines and markers
  gStyle->SetMarkerSize(4);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineWidth(2);
  
  //gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);
}



void setH2Style()
{
  // always reset everything first (doesn't seem to really work hm)
  gStyle->Reset();
  setH1Style();

  gStyle->SetPalette( kBird ); // https://root.cern.ch/doc/master/classTColor.html#C06
  gStyle->SetOptStat(0);

  //gStyle->SetPadTopMargin(0.025);
  //gStyle->SetPadBottomMargin(0.19);
  //gStyle->SetPadLeftMargin(0.08);
  //gStyle->SetPadRightMargin(0.015);

  gStyle->SetTextFont(42);
  gStyle->SetTextSizePixels(23);
  gStyle->SetPaintTextFormat(".3g");
}



template <typename Plot> 
void stylePlot(Plot *plot, 
               const int useColor, const double colorAlpha, const int fillStyle, 
               const int markStyle, const double markSize, 
               const int lineStyle, const int lineWidth, 
               const std::string &mainTitle = "") 
{
  if (useColor > -1) {
    plot->SetFillColorAlpha(useColor, colorAlpha);
    plot->SetMarkerColor(useColor);
    plot->SetLineColor(useColor);
  }

  if (fillStyle > -1)
    plot->SetFillStyle(fillStyle);

  if (markStyle > -1)
    plot->SetMarkerStyle(markStyle);

  plot->SetMarkerSize(markSize);

  if (lineStyle > -1)
    plot->SetLineStyle(lineStyle);

  if (lineWidth > -1)
    plot->SetLineWidth(lineWidth);

  plot->SetTitle( mainTitle.c_str() );
}



template <typename Plot>
void axisPlot(Plot *plot,
              const double yMin, const double yMax, const int yDiv,
              const std::string &yTxt, const double ySiz, const double yOff, const double yLab,
              const double xMin, const double xMax, const int xDiv,
              const std::string &xTxt, const double xSiz, const double xOff, const double xLab) 
{
  // if we don't want the range setter, just give an invalid range
  // btw thanks ROOT: Y U NO SAME METHOD FOR AXIS RANGE FOR HISTOGRAM AND GRAPH??? https://root.cern.ch/how/how-set-ranges-axis

  axisRange(*plot, yMin, yMax, xMin, xMax);

  plot->GetYaxis()->SetTitle(yTxt.c_str());
  plot->GetYaxis()->SetTitleSize(ySiz);
  plot->GetYaxis()->SetTitleOffset(yOff);
  plot->GetYaxis()->SetLabelSize(yLab);
  plot->GetYaxis()->SetNdivisions(yDiv);

  plot->GetXaxis()->SetTitle(xTxt.c_str());
  plot->GetXaxis()->SetTitleSize(xSiz);
  plot->GetXaxis()->SetTitleOffset(xOff);
  plot->GetXaxis()->SetLabelSize(xLab);
  plot->GetXaxis()->SetNdivisions(xDiv);
}



void axisRange(TGraph &plot,
               const double yMin, const double yMax, const double xMin, const double xMax)
{
  if (yMin < yMax) {
    // for graphs histograms are just to draw the frame
    // if it has an error then that would draw a line at 0
    // until ROOT 6.14.02 this bug is there
    plot.GetHistogram()->Sumw2(false);
    plot.GetHistogram()->SetMinimum(yMin);
    plot.GetHistogram()->SetMaximum(yMax);
  }
  if (xMin < xMax)
    plot.GetXaxis()->SetLimits(xMin, xMax);
}



void axisRange(TH1 &plot,
               const double yMin, const double yMax, const double xMin, const double xMax)
{
  if (yMin < yMax)
    plot.GetYaxis()->SetRangeUser(yMin, yMax);
  if (xMin < xMax)
    plot.GetXaxis()->SetRangeUser(xMin, xMax);
}



void styleLegend(TLegend *leg, 
                 const int nColumn, const int fillColor, const int borderSize, 
                 const int txtFont, const double txtSize, 
                 const std::string &legHead = "")
{
  leg->SetNColumns(nColumn);
  leg->SetFillStyle(0);
  leg->SetFillColor(fillColor);
  leg->SetBorderSize(borderSize);
  leg->SetTextFont(txtFont);
  leg->SetTextSize(txtSize);
  leg->SetHeader( legHead.c_str() );
}



void putLegend(TLegend *leg, const double x1, const double x2, const double y1, const double y2)
{
  leg->SetX1(x1); leg->SetX2(x2);
  leg->SetY1(y1); leg->SetY2(y2);
}



// to add back the under/overflows to first/last bin
void add_uoflow_bin(TH1 *hist) {
  const int nBin = hist->GetNbinsX();
  hist->SetBinContent(1, hist->GetBinContent(1) + hist->GetBinContent(0));
  hist->SetBinError(1, std::sqrt((hist->GetBinError(1) * hist->GetBinError(1)) + (hist->GetBinError(0) * hist->GetBinError(0))));

  hist->SetBinContent(0, 0.);
  hist->SetBinError(0, 0.);

  hist->SetBinContent(nBin, hist->GetBinContent(nBin) + hist->GetBinContent(nBin + 1));
  hist->SetBinError(nBin, std::sqrt((hist->GetBinError(nBin) * hist->GetBinError(nBin)) + (hist->GetBinError(nBin + 1) * hist->GetBinError(nBin + 1))));

  hist->SetBinContent(nBin + 1, 0.);
  hist->SetBinError(nBin + 1, 0.);
}

#endif
