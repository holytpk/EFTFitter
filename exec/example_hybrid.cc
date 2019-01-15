// -*- C++ -*-
// example code using the EFTFitter plugin
// plugins rely on ROOT6 and C++14 env to work
// utility tool execMacro.sh compiles and executes the macro, to use it do: ./exec/execMacro.sh exec/example_hybrid.cc
// version for doing hybrid fits - two examples given: asymmetry and pure rate

#include "../src/EFTFitter.h"

int main() {
  using EFT = EFTFitter;

  // construction done with the key we want to treat as data, lambda, fit mode and stat mode
  EFT eft("data", 1., EFT::Fit::hybrid, EFT::Stat::xsec);

  // these are just to avoid writing them out repeatedly
  const std::string input_dir = "./example_fit_input/", output_dir = "./example_fit_output/";
  const std::string histName = "TTbarSpinDensityMatrix/some_histogram", sumName = "TTbarSpinDensityMatrix/sumWgt_noCut";
  const int nRebin = 5;
  const double k_nnlo_lo = 1.667296656, br_tt_2l = 0.041062412; // NNLO/LO k-factor

  // tell the framework know how the distributions are to be transformed for MC
  // must be done before any input is added as to perform the consistency checks on function result instead of distributions themselves
  // in this example assume that the distribution has 6 bins
  // asymmetry is then just ratio between difference and sum between bins 1-3 and 4-6
  eft.setHybridTransformation(std::function<std::vector<std::array<double, 2>> (const std::vector<std::array<double, 2>> &)> ([] (const auto &v_binC) {
        const int nBin = v_binC.size();

        // distributions are stored as bin content = {xsec, count, bin fractions}
        // and the hybrid transformation must result in {xsec, count, actual distribution}
        // in case one doesnt care about rates - like in this example - set them to 0
        // non-0 rates are only for when you want to fit with rate too - the code just adds another term for this (assumed uncorrelated to the rest)
        // in this case ensure the uncertainty is correct - rate uncertainty taken from the data content and not from matrix
        std::vector<std::array<double, 2>> v_tmpC = {{0., 0.}, {0., 0.}};

        double sum13 = 0., sum46 = 0.;
        for (const int index : {iB, iB + 1, iB + 2})
          sum13 += v_binC.at(index).at(0);
        for (const int index : {iB + 3, iB + 4, iB + 5})
          sum46 += v_binC.at(index).at(0);

        // evaluate the asymmetry
        const double sum16 = sum13 + sum46, diff61 = sum46 - sum13;
        const double asymm_val = diff61 / sum16;

        // set the uncertainty to 0 - doesnt affect the fit result
        // if one plans to plot these with drawHistogram() then of course it needs to be computed
        v_tmpC.push_back({asymm_val, 0.});

        return v_tmpC;
      }));

  // special case of doing pure rate fits - in this case bin content should be just {xsec, count}
  // {xsec, count, 1.} will not work; in hybrid mode bin content is interpreted as {xsec, count, actual distribution}
  // alternatively {0., 0., xsec/count} will work, but one needs to specify the uncertainty in a matrix
  //eft.setHybridTransformation(std::function<std::vector<std::array<double, 2>> (const std::vector<std::array<double, 2>> &)> ([] (const auto &v_binC) {
  //      // ensure the error is correct - at least for data (here it is 54.37 pb, no count error since it isnt used)
  //      // matrices are ignored in this case, so this error is used directly
  //      std::vector<std::array<double, 2>> v_tmpC = {{v_binC.at(0).at(0), 54.37}, {v_binC.at(1).at(0), 0.}};
  //      return v_tmpC;
  //}));

  // for real data one directly adds using addHybridData() assuming it's already transformed matching the desired format
  // or optionally supply a function to transform it, which is not necessarily the same as in setHybridTransformation()
  // all that matters is that they both give a bin content of the same format
  // useful when data transformation can't be expressed simply - eg asymmetry accounting bin correlations which aren't needed for MC
  // note that this method has no consistency check whatsoever, so user needs to be sure that it is consistent
  //eft.addHybridData(input_dir + "unfolded_data.root", histName, {0., 0.}, Stat::xsec);

  // and MC following the syntax op1_val1--op2_val2-- ... --opN_valN
  // all operators to be considered must be present in key - ie dont write c2_4 when doing c1 c2 c3 fits, write c1_0--c2_4--c3_0
  // only the Sample::all types are considered for interpolation
  // xsec given is some dummy values (with k-factor applied), last arg stands for the kind of histogram: count vs xsec
  // more explanation in header
  // SM
  eft.addRawInput("c1_0", EFT::Sample::all, input_dir + "c1_0.root", 
                  histName, sumName, nRebin, (k_nnlo_lo / br_tt_2l) * 20.4847, EFT::Stat::count);

  // 1D inputs - here c1 = 10 and c1 = 272 is chosen as raw inputs
  eft.addRawInput("c1_10", EFT::Sample::all, input_dir + "c1_10.root", 
                  histName, sumName, nRebin, (k_nnlo_lo / br_tt_2l) * 48.1389, EFT::Stat::count);

  eft.addRawInput("c1_272", EFT::Sample::all, input_dir + "c1_272.root", 
                  histName, sumName, nRebin, (k_nnlo_lo / br_tt_2l) * 20310.85, EFT::Stat::count);

  // prepare the base for interpolation ie compute individual contribution at 1
  eft.prepareInterpolationBase();

  // in case of fitting to MC - the transformation is as specified in setHybridTransformation()
  // assign as data a particular key of choice
  eft.assignAsData("c1_0", EFT::Sample::all);

  /*/ data needs to be in the list to be drawn - in the braces are key, type and legend text
  std::vector<std::tuple<std::string, EFT::Sample, std::string>> vt_keySampleLegend;
  vt_keySampleLegend.push_back({"data", EFT::Sample::all, "Data"});
  vt_keySampleLegend.push_back({"c1_1", EFT::Sample::all, "c1 = 1"});
  vt_keySampleLegend.push_back({"c1_0", EFT::Sample::all, "SM"});

  // args are just filename and customizing the axes and flags for dividing by bin width etc
  // also can control what to put in ratio plot (more detail in header)
  eft.drawHistogram(vt_keySampleLegend, 
                    output_dir + "cQq11_snake_spinCorr", "#sigma [pb]", "Index", 0.0001, 599.9999, 0.5001, 1.9999, false, "", false, "simple");
  */

  // ensure the added matrices are consistent with transformed distributions
  // grab the total stat error matrix - as usual matrix name is such that file->Get("some_matrix_name") works
  // can also partially extract along the diagonal, pass a vector of bin index range eg {{1, 6}, {115, 120}} as last arg
  eft.readCovMatRoot("statErr", inDir + "/unfolded_data.root", "some_matrix_name");

  // can add more matrices if needed
  eft.readCovMatRoot("totalSyst", inDir + "/unfolded_data.root", "other_matrix_name");

  // alternatively if one wants no bbb correlations - just diag(1) - matrix key is "statCorr"
  // the "statBin" taking bbb covariance is created off data automatically when "statCorr" is given
  // useful if the actual matrix is missing and one just wants to "get going"
  //eft.autoCovMatStatCorr();

  // note: this method isnt compulsory; if a final matrix is already available in root format, just give it as finalCov
  // make the final covariance matrix to be used in fit
  // to do this a function returning a TMatrixD object has to be provided with map of covmat as its only arg
  // said function must not modify other elements in the map!!
  // example here is simply summing totalSyst, statErr and rateErr matrices
  //eft.makeFinalCovMat(std::function<TMatrixD (const std::map<std::string, TMatrixD> &)> ([] (const auto &map) {
  //      // this function can't have side effects, so no std::cout :((

  //      TMatrixD outMat(map.at("statErr"), TMatrixD::kPlus, map.at("rateErr"));
  //      outMat += map.at("totalSyst");
  //      return outMat;
  //    }));

  // draw and store all available covmats as TH2 (can also provide a vector of keys to draw only specific ones)
  // also produces the matrix in text format
  //eft.drawCovMat(output_dir);

  // make the range to compute the chi2; in this case [min, max: step]
  const std::vector<double> v_opPoint = makeInterval(-5., 5., 0.00001);
  //fillInterval(v_opPoint, -3., 3., 0.00001); // alternatively for non-uniform intervals - second arg has to be the same as last vector element

  // make the list of keys to fit
  // in this case op1 and op2 over the specified list of points
  // the points are made gridwise considering all operators
  eft.listKeyToFit({ {"c1", v_opPoint} });

  const std::vector<EFT::Sample> v_sample = {EFT::Sample::all, EFT::Sample::linear};
  eft.computeFitChi2(v_sample);

  // now we provide the op for which we wanna draw the dChi2 plot on
  // insert into map: op key, op string in plot, op range (if none, select by dChi2), y axis range, x axis range
  std::map<std::string, std::tuple<std::string, std::vector<double>, std::array<double, 2>, std::array<double, 2>>> m_1D_all;
  m_1D.insert({"c1", { "c1", {/* op range in min, max */}, {0., 9.999}, {-1.499, 1.999} }});

  // in this case we just draw for cQq11 - also include a filename for the resulting plot
  eft.draw1DChi2(m_1D, output_dir + "c1_constraint", v_sample);

  // same thing for 2D - but here only the contour is drawn
  // to brace-init the value, provide an extra pair of braces to obey the brace-elision rule
  // insert into map: op pair key (y-x), op1 string in plot, its axis range, op2 string in plot, its axis range
  //std::map<std::array<std::string, 2>, std::array<std::pair<std::string, std::array<double, 2>>, 2>> m_2D;
  //m_2D.insert({ {"c1", "c2"}, {{{"c1", {-0.799, 0.799}}, {"c2", {-0.799, 1.199}}}} });

  //eft.draw2DChi2(m_2D_all, output_dir + "c1_c2_constraint", v_sample);

  return 0;
}

