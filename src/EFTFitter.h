// -*- C++ -*-
// takes templates from MC and data, performs some chi2 fit to set constraints on EFT operators
// includes some methods for plotting output
// updated 31/08/2018

#ifndef EFTFITTER_H
#define EFTFITTER_H

#include "TemplateUtil.h"

//#include <random>
//#include <chrono>

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

#include "TGraph.h"
#include "TGraphAsymmErrors.h"

#include "TAxis.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TLatex.h"

#include "TMatrixD.h"

class EFTFitter {
 public:
  // what kind of fit we're gonna be using
  enum class Fit : char {
    absolute = 'a', // absolute: use the xsec/counts and shape together in the fit
    shape = 's', // shape comparison only
    hybrid = 'h' // hybrid fit where rate and shape are separately fitted
  };

  // what kind of statistic to report
  enum class Stat : char {
    count = 'c', // use event counts
    xsec = 'x' // use xsec (no unit check, so ensure it's consistent throughout)
  };

  /// constructors
  EFTFitter(const std::string &dataName_, const double eftLambda_ = 1., 
            const Fit fitMode_ = Fit::absolute, const Stat statMode_ = Stat::xsec, const double shapeSum_ = 1.);

  // what kind of sample we're working with
  enum class Sample : int {
    all = 0, // sum of linear, quadratic
    linear = 1, // pure interference linear dim6 term
    quadratic = 2 // pure 'resonance' quadratic dim8 term
  };

  /// dump the histogram from file into map after normalization; can also add to existing key instead of overwrite
  /// if fit mode is Fit::hybrid, normalizedSum is exactly that - assigned to count/xsec based on arg sumStat
  /// in this case if normalizedSum error (second entry) is 0 and sumStat is Stat::count, it is set to sqrt of value
  /// note that template integral is completely ignored in this case - it supplies only the shape
  /// otherwise, if fit mode is absolute, the integral of the template is assigned to count/xsec (likewise for its error)
  /// if sumStat is Stat::count and normalizedSum value is not 0, lumi normalization is done (assuming MC template of weighted count)
  /// to not do lumi normalization just set normalizedSum value to 0 (its error is ignored anyway)
  /// Fit::shape acts like absolute, but it doesn't matter either way (only shape and shapeSum are used in fit)
  void addRawInput(const std::string &keyName, const Sample sampleType, 
                   const std::string &fileName, const std::string &histName, const std::string &sumName = "",
                   const int nRebin = 1, const std::array<double, 2> &normalizedSum = {0., 0.}, const Stat sumStat = Stat::count,
                   const bool addIfPresent = false);

  /// draw histograms into canvas
  void drawHistogram(const std::vector< std::tuple<std::string, Sample, std::string> > &vt_keySampleLegend,
                     const std::string &plotName, const std::string &yLabel, const std::string &xLabel,
                     const double histMin, const double histMax, const double ratioMin, const double ratioMax, 
                     const bool drawLogY = false, const std::string &legHeader = "", const bool divideBinWidth = false);

  /// make automatic stat corr matrix of diag(1)
  void autoCovMatStatCorr();

  /// make automatic covariance matrix of data sample base on rate uncertainties (e. g. lumi)
  /// matrix coeff c_ii = sq(rateUnc * content of bin i)
  void autoCovMatRate(const double rateUnc = 0.01);

  /// read covariance matrix from a root file (assumes it's stored as TH2, can also select parts of it)
  void readCovMatRoot(const std::string &keyMat, const std::string &fileName, const std::string &histName,
                      std::array<int, 2> a_binIndex = {-1, -1});

  /// make the final covariance matrix to be used in fit - user needs to provide a function detailing how this is done
  /// example is provided in the sample execution script - said function must not modify other elements in the map!!
  /// technical remark: std::function wrapping lambdas are needed here due to the in-place nature of lambda types
  /// somehow it works for the free templates though; related to namespace blocking? hmm...
  void makeFinalCovMat(const std::function<TMatrixD (const std::map<std::string, TMatrixD> &)> &func);

  /// draw the covariance matrices (if no keys provided, draw all available in the map)
  void drawCovMat(const std::string &dirName = "./", const std::vector<std::string> &v_keyMat = {}) const;

  /// prepare the base content (op = 1) needed for interpolation; follows arXiv:1611.01165
  void prepareInterpolationBase();

  /// compute the test statistic for a given set of keys and fill into map
  void computeFitChi2(const std::map<std::string, std::vector<double>> &m_opGrid);

  /// draw all the 1D dChi2 graphs for each operator
  /// tuple is op text (as to appear in legend), op range to plot, y and x range of plot
  void draw1DChi2(const std::map<std::string, std::tuple<std::string, std::vector<double>, 
                  std::array<double, 2>, std::array<double, 2>>> &mt_opInfo,
                  const std::string &dirName, const std::vector<Sample> &v_sample = {Sample::all, Sample::linear}) const;

  /// draw all the 2D dChi2 contours for each operator pair
  /// tuple is op1-op2 text (as to appear in axes), op1-op2 range to plot
  /// op1 will be on y-axis, op2 on x-axis
  void draw2DChi2(const std::map<std::array<std::string, 2>, 
                  std::array<std::pair<std::string, std::array<double, 2>>, 2>> &mt_opPair,
                  const std::string &dirName, const std::vector<Sample> &v_sample = {Sample::all, Sample::linear}) const;

  /// clear all contents, leaving only those initialized in ctor intact
  void clearContent();

 private:

  /// make automatic covariance matrix of data sample based on stat uncertainties
  /// matrix coeff c_ij = unc of bin i * unc of bin j
  void autoCovMatStatBin();

  /// ensure that the binning of all input hist are compatible
  bool checkInputBin(const std::unique_ptr<TH1D> &hist) const;

  /// quickly assign "in-situ" xsec by scaling from SM for binContents without it
  void assignInSituXsec(std::vector<std::array<double, 2>> &v_binContent);

  /// construct interpolated content for a given key
  std::vector<std::array<double, 2>> interpolateOpValue(const std::string &keyName, const Sample sampleType);

  /// parse the operator name and its coeff value from key
  std::map<std::string, double> parseOpFromKey(const std::string &keyName) const;

  /// and the other way around
  static std::string parseKeyFromOp(const std::map<std::string, double> &m_opVal);

  /// fix the formatting of numbers in key
  std::string fixKeyFormat(const std::string &keyName) const;

  /// ensure that the fit value map cover the same operator set and each operator-value set occurs only once
  bool checkOpSet(const std::string &keyName, const Sample sampleType);

  /// unpack the grid into individual operator-value set and return them as keys
  void unpackOpGrid(std::vector<std::string> &v_opGrid, 
                    const std::map<std::string, std::vector<double>> &m_opGrid, 
                    std::map<std::string, std::vector<double>>::const_iterator i_opGrid, 
                    const std::string &iniStr) const;

  /// convert contents back to histogram for drawing
  std::unique_ptr<TH1D> convertContentToHist(const std::string &keyName, const Sample sampleType, const bool divideBinWidth = false);

  /// members
  /// sample name to be treated as data and fit against
  std::string dataName;

  /// lambda; default is 1 TeV
  double eftLambda;

  /// what fit to do: absolute vs shape
  Fit fitMode; 

  /// what statistic to report: xsec vs count
  Stat statMode;

  /// what should the templates sum up to in case of shape comparison
  double shapeSum;

  bool hasData;

  /// set of operator names
  std::vector<std::string> v_opName;

  /// binning used in the templates
  std::vector<double> v_rawBin;

  /// map of covariance matrices
  /// cheats https://root.cern.ch/root/htmldoc/guides/users-guide/LinearAlgebra.html
  std::map<std::string, TMatrixD> m_covMat;

  /// function to process the inputs in covariance matrix map into the final matrix
  std::function<TMatrixD (const std::map<std::string, TMatrixD> &)> f_finalCov;

  /// map of contents to be fit on
  std::map< std::pair<std::string, Sample>, std::vector<std::array<double, 2>> > m_binContent;

  /// map of contents where each operator had its coeff = 1 and all others are 0
  std::map< std::pair<std::string, Sample>, std::vector<std::array<double, 2>> > m_op1Eq1;

  // a map of bin contents where 2 op had their coeff = 1 and all others are 0
  std::map< std::array<std::string, 2>, std::vector<std::array<double, 2>> > m_op2Eq1;

  /// map of the fit chi2 for all requested points
  std::map< std::pair<std::string, Sample>, double > m_fitChi2;
};



namespace FitUtil {
  /// a namespace for some random stuff that's convenient to have
  /// could be moved to its own file if the fwk gets large enough

  /// style of TH1 - reproducing the TDRStyle macro
  void setH1Style();

  /// style of TH2
  void setH2Style();

  /// extract the binning used in a hist
  std::vector<double> extractBin(const std::unique_ptr<TH1D> &hist);

  /// extract the bin contents and errors of a hist (assumes symmetrical errors)
  std::vector<std::array<double, 2>> extractContentError(const std::unique_ptr<TH1D> &hist);

  /// styling of hist, graph
  template <typename Plot> void stylePlot(const std::unique_ptr<Plot> &plot, 
                                          const int useColor, const double colorAlpha, const int fillStyle, 
                                          const int markStyle, const double markSize, 
                                          const int lineStyle, const int lineWidth, 
                                          const std::string &mainTitle = "");

  /// range and labeling of y, x axes
  template <typename Plot> void axisPlot(const std::unique_ptr<Plot> &plot,
                                         const double yMin, const double yMax,
                                         const std::string &yTxt, const double ySiz, const double yOff, const double yLab,
                                         const double xMin, const double xMax,
                                         const std::string &xTxt, const double xSiz, const double xOff, const double xLab);

  /// best overloads ever; if only we could slug them in axisPlot()...
  void axisRange(TGraph &plot,
                 const double yMin, const double yMax, const double xMin, const double xMax);

  void axisRange(TH1 &plot,
                 const double yMin, const double yMax, const double xMin, const double xMax);

  /// styling, positioning of plot legend
  void styleLegend(const std::unique_ptr<TLegend> &leg, 
                   const int nColumn, const int fillColor, const int borderSize, 
                   const int txtFont, const double txtSize, 
                   const std::string &legHead = "");

  void putLegend(const std::unique_ptr<TLegend> &leg, const double x1, const double x2, const double y1, const double y2);
}



/// this overload is for Sample printouts and string conversion
inline std::ostream& operator<<(std::ostream& out, const EFTFitter::Sample &samp)
{
  if (samp == EFTFitter::Sample::all)
    out << "all";
  if (samp == EFTFitter::Sample::linear)
    out << "linear";
  if (samp == EFTFitter::Sample::quadratic)
    out << "quadratic";

  return out;
}



/// this overload is for printAll to print twin numbers (val, err or x, y what else...) without resorting to printAll
template <typename Content> std::ostream& operator<<(std::ostream& out, const std::array<Content, 2> &array)
{
  out << array.at(0) << ", " << array.at(1);
  return out;
}

#endif
