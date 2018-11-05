// -*- C++ -*-
// example code using the EFTFitter plugin
// plugins rely on ROOT6 and C++14 env to work
// utility tool execMacro.sh compiles and executes the macro, to use it do: ./exec/execMacro.sh exec/example_fit.cc

#include "../src/EFTFitter.h"

int main() {
  using EFT = EFTFitter;

  // construction done with the key we want to treat as data, lambda, fit mode and event mode
  EFT eft("data", 1., EFT::Fit::absolute, EFT::Stat::xsec);

  // these are just to avoid writing them out repeatedly
  const std::string input_dir = "./example_fit_input/", output_dir = "./example_fit_output/";
  const std::string histName = "TTbarSpinDensityMatrix/snake_spinCorr", sumName = "TTbarSpinDensityMatrix/sumWgt_noCut";
  const int nRebin = 5;
  const double k_nnlo_lo = 22. * 1.667296656, br_tt_2l = 0.041062412; // NNLO/LO k-factor

  // add the input file and hist names (including if necessary the sum of weight hist for normalization)
  // hist name is such that file->Get(name) works
  // please ensure there is exactly 1 input with dataName as in ctor (assign 0 xsec to deactivate normalization)
  eft.addRawInput("data", EFT::Sample::all, inDir + "/unfolded_data.root", 
                  histName, "", 1, 0., EFT::Stat::xsec);

  // and MC following the syntax op1_val1--op2_val2-- ... --opN_valN
  // all operators to be considered must be present in key - eg dont write c2_4 when doing c1 c2 c3 fits, write c1_0--c2_4--c3_0
  // only the Sample::all types are considered for interpolation
  // more explanation in header
  // SM
  eft.addRawInput("cQq11_0", EFT::Sample::all, input_dir + "cQq11_0p0_rwgt_snake.root", 
                  histName, sumName, nRebin, (k_nnlo_lo / br_tt_2l) * 20.4847, EFT::Stat::count);

  // 1D inputs
  eft.addRawInput("cQq11_10", EFT::Sample::all, input_dir + "cQq11_10p0_rwgt_snake.root", 
                  histName, sumName, nRebin, (k_nnlo_lo / br_tt_2l) * 48.1389, EFT::Stat::count);

  eft.addRawInput("cQq11_272", EFT::Sample::all, input_dir + "cQq11_272p0_rwgt_snake.root", 
                  histName, sumName, nRebin, (k_nnlo_lo / br_tt_2l) * 20310.85, EFT::Stat::count);

  // prepare the base for interpolation ie compute individual contribution at 1
  eft.prepareInterpolationBase();

  /*/ data needs to be in the list to be drawn - in the braces are key, type and legend text
  std::vector<std::tuple<std::string, EFT::Sample, std::string>> vt_keySampleLegend;
  vt_keySampleLegend.push_back({"data", EFT::Sample::all, "Data"});
  vt_keySampleLegend.push_back({"cQq13_0--cQq11_1--ctq1_0", EFT::Sample::all, "cQq11 = 1"});
  vt_keySampleLegend.push_back({"cQq13_0--cQq11_0--ctq1_0", EFT::Sample::all, "SM"});

  // args are just filename and customizing the axes and flags for dividing by bin width etc
  // also can control what to put in ratio plot (more detail in header)
  eft.drawHistogram(vt_keySampleLegend, 
                    output_dir + "cQq11_snake_spinCorr", "#sigma [pb]", "Index", 0.0001, 599.9999, 0.5001, 1.9999, false, "", false, "simple");
  */

  // make the rate covmat assuming 2.5% lumi uncertainty - matrix key is "rateErr"
  eft.autoCovMatRate(0.025);

  // grab the total stat error matrix - as usual matrix name is such that file->Get("some_matrix_name") works
  // can also partially extract along the diagonal, pass a vector of bin index range eg {{1, 6}, {115, 120}} as last arg
  eft.readCovMatRoot("statErr", inDir + "/unfolded_data.root", "some_matrix_name");

  // alternatively if one wants no bbb correlations - just diag(1) - matrix key is "statCorr"
  // the "statBin" taking bbb covariance is created off data automatically when "statCorr" is given
  //eft.autoCovMatStatCorr();

  // can add more matrices if needed
  eft.readCovMatRoot("totalSyst", inDir + "/unfolded_data.root", "other_matrix_name");

  // make the final covariance matrix to be used in fit
  // to do this a function returning a TMatrixD object has to be provided with map of covmat as its only arg
  // said function must not modify other elements in the map!!
  // example here is simply summing totalSyst, statErr and rateErr matrices
  eft.makeFinalCovMat(std::function<TMatrixD (const std::map<std::string, TMatrixD> &)> ([] (const auto &map) {
        // this function can't have side effects, so no std::cout :((

        TMatrixD outMat(map.at("statErr"), TMatrixD::kPlus, map.at("rateErr"));
        outMat += map.at("totalSyst");
        return outMat;
      }));

  // draw and store all available covmats as TH2 (can also provide a vector of keys to draw only specific ones)
  // also produces the matrix in text format
  //eft.drawCovMat(output_dir);

  // make the range to compute the chi2; in this case [min, max: step]
  const std::vector<double> v_opPoint = makeInterval(-5., 5., 0.00001);
  //fillInterval(v_opPoint, -3., 3., 0.00001); // alternatively for non-uniform intervals - second arg has to be the same as last vector element

  // make the list of keys to fit
  // in this case op1 and op2 over the specified list of points
  // the points are made gridwise considering all operators
  eft.listKeyToFit({ {"cQq11", v_opPoint}, {"ctq1", v_opPoint} });

  const std::vector<EFT::Sample> v_sample = {EFT::Sample::all, EFT::Sample::linear};
  eft.computeFitChi2(v_sample);

  // now we provide the op for which we wanna draw the dChi2 plot on
  // insert into map: op key, op string in plot, op range (if none, select by dChi2), y axis range, x axis range
  std::map<std::string, std::tuple<std::string, std::vector<double>, std::array<double, 2>, std::array<double, 2>>> m_1D_all;
  m_1D.insert({"cQq11", { "cQq11", {/* op range in min, max */}, {0., 9.999}, {-1.499, 1.999} }});

  // in this case we just draw for cQq11 - also include a filename for the resulting plot
  eft.draw1DChi2(m_1D, output_dir + "cQq11_constraint", v_sample);

  // same thing for 2D - but here only the contour is drawn
  // to brace-init the value, provide an extra pair of braces to obey the brace-elision rule
  // insert into map: op pair key (y-x), op1 string in plot, its axis range, op2 string in plot, its axis range
  std::map<std::array<std::string, 2>, std::array<std::pair<std::string, std::array<double, 2>>, 2>> m_2D;
  m_2D.insert({ {"cQq11", "ctq1"}, {{{"cQq11", {-0.799, 0.799}}, {"ctq1", {-0.799, 1.199}}}} });

  eft.draw2DChi2(m_2D_all, output_dir + "cQq11_ctq1_constraint", v_sample);

  return 0;
}

