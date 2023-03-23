// -*- C++ -*-
// author: afiq anuar
// short: takes templates from MC and data, performs some chi2 fit to set constraints on EFT operators
// todo: it may be useful to save the weight vector of each dofs in the chi2

#ifndef EFTFITTER_H
#define EFTFITTER_H

#include "TemplateUtil.h"
#include <random>

#include "TMath.h"
#include "TMatrixD.h"
//#include "TDecompLU.h"

class TH1;
class TH1D;

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
  EFTFitter(const std::string &dataName_, const double &eftLambda_ = 1., 
            const Fit fitMode_ = Fit::absolute, const Stat statMode_ = Stat::xsec, const double &shapeSum_ = 1.);

  /// self-explanatory
  void setShapeSum(const double &shapeSum_ = 1.);

  // what kind of sample we're working with
  enum class Sample : int {
    all = 0, // sum of linear, quadratic
    linear = 1, // pure interference linear dim6 term
    quadratic = 2 // pure 'resonance' quadratic dim8 term
  };

  /// dump the histogram from file into map after normalization; except for data, can also add to existing key instead of overwrite
  /// in Fit::hybrid, this method can't be used to add data, instead a dedicated method is supplied
  /// for MC, it acts like in absolute (see below), except some checks are hybrid-based
  /// in Fit::absolute, the integral of the template is assigned to count/xsec (likewise for its error)
  /// if sumStat is Stat::count and normalizedSum value is not 0, lumi normalization is done (assuming MC template of weighted count)
  /// to not do lumi normalization just set normalizedSum value to 0 (its error is ignored anyway)
  /// in Fit::xsec, normalizedSum is ignored, so templates should already have proper normalization
  /// Fit::shape acts like absolute, but of course only shape and shapeSum are used in fit
  /// keyName assumes the syntax is op1_val1--op2_val2-- ... --opN_valN
  void addRawInput(const std::string &keyName, const Sample sampleType, 
                   const std::string &fileName, const std::string &histName, const std::string &sumName = "",
                   const int nRebin = 1, const std::array<double, 2> &normalizedSum = {0., 0.}, const Stat sumStat = Stat::count,
                   const bool addIfPresent = false);

  /// set the function to transform the binContent templates to the desired form to be fitted
  /// must be from one binContent to another, has no side effects yak yak
  /// basically the same restrictions binding f_finalCov apply here as well
  /// oh yeah, it must also be completely self-contained meaning no accessing class members etc
  void setHybridTransformation(const std::function<std::vector<std::array<double, 2>> (const std::vector<std::array<double, 2>> &)> &func);

  /// dedicated method to add data in Fit::hybrid
  /// assumes that the template is already in a form that is compatible with the result of f_hybridTransform(v_binContent)
  /// or reachable from it through argument func
  /// no check that such is the case, you're on your own
  /// takes the histogram AS IT IS and tack on the normalizedSum as xsec/count
  void addHybridData(const std::string &fileName, const std::string &histName, 
                     const std::array<double, 2> &normalizedSum = {0., 0.}, const Stat sumStat = Stat::count, 
                     const std::function<std::vector<std::array<double, 2>> (const std::vector<std::array<double, 2>> &)> &func = nullptr);

  /// draw histograms into canvas
  /// most options involve controlling the axes (names should be self-explanatory)
  /// ratioMode are either none, simple or covariance - none is none, simple means simple MC/data ratio, while covariance does MC - data / error
  /// covariance mode will only work if final covariance matrix is already available
  void drawHistogram(const std::vector< std::tuple<std::string, Sample, std::string> > &vt_keySampleLegend,
                     const std::string &plotName, const std::string &yLabel, const std::string &xLabel,
                     const double &histMin, const double &histMax, const double &ratioMin, const double &ratioMax, 
                     const bool drawLogY = false, const std::string &legHeader = "", 
                     const bool divideBinWidth = false, const std::string &ratioMode = "simple");

  /// make automatic stat corr matrix of diag(1)
  void autoCovMatStatCorr();

  /// make automatic covariance matrix of data sample base on rate uncertainties (e. g. lumi)
  /// matrix coeff c_ij = sq(rateUnc * rateUnc * content of bin i * content of bin j)
  void autoCovMatRate(const double &rateUnc = 0.01);

  /// read covariance matrix from a root file (assumes it's stored as TH2, can also select parts of it)
  void readCovMatRoot(const std::string &keyMat, const std::string &fileName, const std::string &histName,
                      const std::vector< std::array<int, 2> > &v_endbin = {});

  /// make the final covariance matrix to be used in fit - user needs to provide a function detailing how this is done
  /// example is provided in the sample execution script - said function must not modify other elements in the map!!
  /// technical remark: std::function wrapping lambdas are needed here due to the in-place nature of lambda types
  /// somehow it works for the free templates though; related to namespace blocking? hmm...
  void makeFinalCovMat(const std::function<TMatrixD (const std::map<std::string, TMatrixD> &)> &func);

  /// trivial case of final covmat being the sum of everything else that was read in
  /// why was this not in there right at the beginning???
  void makeFinalCovMat(std::vector<std::string> &&names = {});

  /// draw the covariance matrices (if no keys provided, draw all available in the map)
  /// last flag is when one wants to draw correlation rather than covariance matrix
  void drawCovMat(const std::string &dirName = "./", const std::vector<std::string> &v_keyMat = {}, bool do_correlation = false) const;

  /// get data rate according to the inputs and covMat
  std::array<double, 2> getDataRate() const;

  /// prepare the base content (op = 1) needed for interpolation; follows arXiv:1611.01165
  void prepareInterpolationBase();

  /// assign as data something that is taken from interpolation or what else
  void assignAsData(const std::string &keyName, const Sample sampleType, const bool varyWithinError = false);

  /// make the list of keys to fit on
  /// currently grid building code stops at 2 operators != 0
  /// easily disabled when proper minimization is available
  void listKeyToFit(const std::map<std::string, std::vector<double>> &m_opGrid);

  /// FIXME proper minimization of chi2 such that listKeyToFit() doesnt need to supply billions of points
  /// compute the test statistic for a given set of keys and fill into map
  /// can control if want the fit only for a sample
  /// also which bin to ignore in case of shape fit (for the automatic bin dropper algo)
  void computeFitChi2(const std::vector<Sample> &v_sample = {Sample::all, Sample::linear}, int binToIgnore = 1);

  /// save the computed chi2 into a tree
  void saveFitChi2(const std::string &fileName = "./fit");

  /// FIXME in their current forms, drawNDChi2 methods are extremely inefficient
  /// FIXME handling of output file names is awkward at best
  /// FIXME restriction of drawing each operator once per invocation is also terrible
  /// FIXME methods spend most time doing the minimum and interval/contour search
  /// FIXME but there's no need to do this at every single invocation - do it once per request and cache
  /// FIXME structure needs rethinking/recoding
  /// draw all the 1D dChi2 graphs for each operator
  /// tuple is op text (as to appear in legend), op range to plot, y and x range of plot
  void draw1DChi2(const std::map<std::string, std::tuple<std::string, std::vector<double>, 
                  std::array<double, 2>, std::array<double, 2>>> &mt_opInfo,
                  const std::string &plotName, const std::vector<Sample> &v_sample = {Sample::all, Sample::linear}) const;

  /// draw all the 2D dChi2 contours for each operator pair
  /// array is op1-op2 text (as to appear in axes), op1-op2 range to plot
  /// op1 will be on y-axis, op2 on x-axis
  /// dChi2FracScan is a way to remove some points from the edge scan when drawing the contours
  /// only points within (1 - dChi2FracScan) * dChi2Cut < dChi2 < dChi2Cut will be kept for defining the contour 
  void draw2DChi2(const std::map<std::array<std::string, 2>, 
                  std::array<std::pair<std::string, std::array<double, 2>>, 2>> &mt_opPair,
                  const std::string &plotName, const std::vector<Sample> &v_sample = {Sample::all, Sample::linear},
                  const double &dChi2FracScan = 1.) const;

  /// clear contents
  /// clearLevel 0 clears all leaving only those initialized in ctor intact
  /// clearLevel 1 means only the fitChi2 map and list of keys are cleared up
  void clearContent(const int clearLevel = 0);

 private:

  /// make automatic covariance matrix of data sample based on stat uncertainties
  /// matrix coeff c_ij = unc of bin i * unc of bin j
  void autoCovMatStatBin();

  /// ensure that the binning of all input hist are compatible
  bool checkInputBin(const std::unique_ptr<TH1D> &hist) const;

  /// get the correct integral to use for a given binContent based on current fitter settings
  double getContentSum(const std::vector<std::array<double, 2>> &v_binContent) const;

  /// normalize binContent to the {rate, shape} format
  void normalizeContent(std::vector<std::array<double, 2>> &v_binContent, const int iXs);

  /// quickly assign "in-situ" xsec by scaling with intLumi
  void assignInSituXsec(std::vector<std::array<double, 2>> &v_binContent);

  /// construct interpolated content for a given key
  /// will return the raw binContent if it's available and doesn't need to be Poisson-varied
  std::vector<std::array<double, 2>> interpolateOpValue(const std::string &keyName, const Sample sampleType, const bool varyWithinError = false);

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

  /// get the difference between mean and random sample of Poisson distribution around it
  double getDPoissonVariation(const double &mean) const;

  /// members
  /// sample name to be treated as data and fit against
  std::string dataName;
  bool hasData;

  /// lambda; default is 1 TeV
  double eftLambda;

  /// what fit to do: absolute vs shape
  Fit fitMode; 

  /// what statistic to report: xsec vs count
  Stat statMode;

  /// what should the templates sum up to in case of shape comparison
  double shapeSum;

  /// set of operator names
  std::vector<std::string> v_opName;

  /// binning used in the templates
  std::vector<double> v_rawBin;

  /// hybrid mode shenanigan: function to transform the binContent into some other form to be fitted
  std::function<std::vector<std::array<double, 2>> (const std::vector<std::array<double, 2>> &)> f_hybridTransform;

  /// map of covariance matrices
  /// cheats https://root.cern.ch/root/htmldoc/guides/users-guide/LinearAlgebra.html
  std::map<std::string, TMatrixD> m_covMat;

  /// function to process the inputs in covariance matrix map into the final matrix
  std::function<TMatrixD (const std::map<std::string, TMatrixD> &)> f_finalCov;

  /// map of contents to be fit on
  std::map< std::pair<std::string, Sample>, std::vector<std::array<double, 2>> > m_binContent;

  /// map of contents where each operator had its coeff = 1 and all others are 0
  std::map< std::pair<std::string, Sample>, std::vector<std::array<double, 2>> > m_op1Eq1;

  /// a map of bin contents where 2 op had their coeff = 1 and all others are 0
  std::map< std::array<std::string, 2>, std::vector<std::array<double, 2>> > m_op2Eq1;

  /// set of operator points to be fitted (just keys as fit is done for all and linear always)
  std::vector<std::string> v_keyToFit;

  /// map of the fit chi2 for all requested points
  std::map< std::pair<std::string, Sample>, double > m_fitChi2;
};



namespace FitUtil {
  /// a namespace for some random stuff that's convenient to have
  /// could be moved to its own file if the fwk gets large enough

  /// rng for any random stuff
  static std::mt19937_64 rng(std::random_device{}());

  /// extract the binning used in a hist
  std::vector<double> extractBin(TH1 *hist);

  /// extract the bin contents and errors of a hist (assumes symmetrical errors)
  std::vector<std::array<double, 2>> extractContentError(TH1 *hist);

  /// int lumi for normalization
  static constexpr double intLumi = 35922.;
}



/// this overload is for enum printouts and string conversion
inline std::ostream& operator<<(std::ostream& out, const EFTFitter::Fit &fit)
{
  if (fit == EFTFitter::Fit::absolute)
    out << "absolute";
  if (fit == EFTFitter::Fit::shape)
    out << "shape";
  if (fit == EFTFitter::Fit::hybrid)
    out << "hybrid";

  return out;
}



inline std::ostream& operator<<(std::ostream& out, const EFTFitter::Stat &stat)
{
  if (stat == EFTFitter::Stat::count)
    out << "count";
  if (stat == EFTFitter::Stat::xsec)
    out << "xsec";

  return out;
}



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



/// this overload is to print twin numbers (val, err or x, y what else...) without using printAll
template <typename Content> std::ostream& operator<<(std::ostream& out, const std::array<Content, 2> &array)
{
  out << array.at(0) << ", " << array.at(1);
  return out;
}

#endif
