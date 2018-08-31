// -*- C++ -*-
// example code using the EFTFitter plugin
// plugins rely on ROOT6 and C++14 env to work
// utility tool execMacro.sh compiles and executes the macro, to use it do: ./exec/execMacro.sh exec/example_fit.cc

#include "../src/EFTFitter.h"

int main() {
  using EFT = EFTFitter;

  // construction done with the key we want to treat as data, lambda, fit mode and event mode
  EFT eft("data", 1., EFT::Fit::absolute, EFT::Stat::xsec);
  const std::string inDir = "./example_fit_root/", outDir = "./example_fit_output/";
  const std::string hName = "TTbarSpinDensityMatrix/snake_spinCorr", sName = "TTbarSpinDensityMatrix/sumWgt_noCut";
  const std::string sufDat = "_unfolded_data", sufSim = "_part_cut_0";
  const int nRebin = 5;
  const double k_nnlo_lo = 22. * 1.667296656, br_tt_2l = 0.041062412; // NNLO/LO k-factor

  // add the input file and hist names (including if necessary the sum of weight hist for normalization)
  // hist name is such that file->Get(name) works
  // please ensure there is exactly 1 input with dataName as in ctor (assign 0 xsec to deactivate normalization)
  eft.addRawInput("data", EFT::Sample::all, "./covMat_180816/unfolded_data.root", 
                  hName + sufDat, "", 1, 0., EFT::Stat::xsec);

  // and MC following the syntax op1_val1--op2_val2-- ... --opN_valN
  // SM
  eft.addRawInput("cQq11_0", EFT::Sample::all, inDir + "cQq11_0p0_rwgt_snake.root", 
                  hName + sufSim, sName, nRebin, (k_nnlo_lo / br_tt_2l) * 20.4847, EFT::Stat::count);

  // 1D inputs
  eft.addRawInput("cQq11_10", EFT::Sample::all, inDir + "cQq11_10p0_rwgt_snake.root", 
                  hName + sufSim, sName, nRebin, (k_nnlo_lo / br_tt_2l) * 48.1389, EFT::Stat::count);

  eft.addRawInput("cQq11_272", EFT::Sample::all, inDir + "cQq11_272p0_rwgt_snake.root", 
                  hName + sufSim, sName, nRebin, (k_nnlo_lo / br_tt_2l) * 20310.85, EFT::Stat::count);

  // prepare the base for interpolation
  //eft.prepareInterpolationBase();

  /*/ data needs to be in the list to be drawn - in the braces are key, type and legend text
  std::vector<std::tuple<std::string, EFT::Sample, std::string>> vt_keySampleLegend1;
  vt_keySampleLegend1.push_back({"data", EFT::Sample::all, "Data"});
  vt_keySampleLegend1.push_back({"cQq13_0--cQq11_1--ctq1_0", EFT::Sample::all, "cQq11 = 1"});
  vt_keySampleLegend1.push_back({"cQq13_0--cQq11_0--ctq1_0", EFT::Sample::all, "SM"});

  // args are just filename and customizing the axes (consult header for more)
  eft.drawHistogram(vt_keySampleLegend1, 
                    outDir + "cQq11_snake_spinCorr_part_2_raw_1", "#sigma [pb]", "Index", 0.0001, 599.9999, 0.5001, 1.9999);
  */

  // make the rate covmat assuming 2.5% lumi uncertainty - matrix key is "rateErr"
  eft.autoCovMatRate(0.025);

  // grab the total stat error matrix 
  // can also partially extract along the diagonal, pass a bin index range eg {115, 120}
  eft.readCovMatRoot("statErr", "./covMat_180816/Systematics_AllVars_nobw.root", "TotalStatCovMatrix_AllVar_rebinnedA");

  // alternatively if one wants no bbb correlations - just diag(1) - matrix key is "statCorr"
  // the "statBin" taking bbb covariance is created off data automatically when "statCorr" is given
  //eft.autoCovMatStatCorr();

  // syst covmats by Jacob
  const std::string systSuffix = "SystCovMatrix_AllVar_rebinnedA";
  eft.readCovMatRoot("totalSyst", "./covMat_180816/Systematics_AllVars_nobw.root", "Total" + systSuffix, covMat_binRange);

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

  // draw and store all available covmats (can also provide a vector of keys for specific ones)
  //eft.drawCovMat(outDir);

  // make the range to compute the chi2; in this case [min, max: step]
  std::vector<double> v_opPoint = {-1.};
  fillInterval(v_opPoint, -1., 1., 0.001);
  eft.computeFitChi2({ {"cQq11", v_opPoint} });

  // now we provide the op for which we wanna draw the dChi2 plot on
  // insert into map: op key, op string in plot, op range (if none, select by dChi2), y axis range, x axis range
  std::map<std::string, std::tuple<std::string, std::vector<double>, std::array<double, 2>, std::array<double, 2>>> m_1D_all;
  m_1D_all.insert({"cQq11", { "cQq11", {/* op range in min, max */}, {0., 9.999}, {-1.499, 1.999} }});

  // default for last arg draws both the all and linear part
  eft.draw1DChi2(m_1D_all, outDir, {EFT::Sample::all});

  // to brace-init the value, provide an extra pair of braces to obey the brace-elision rule
  // insert into map: op pair key (y-x), op1 string in plot, its axis range, op2 string in plot, its axis range
  std::map<std::array<std::string, 2>, std::array<std::pair<std::string, std::array<double, 2>>, 2>> m_2D_all;
  m_2D_all.insert({ {"cQq11", "ctq1"}, {{{"cQq11", {-0.799, 0.799}}, {"ctq1", {-0.799, 1.199}}}} });

  eft.draw2DChi2(m_2D_all, outDir, {EFT::Sample::all});

  return 0;
}

