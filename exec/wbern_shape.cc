// example code using the EFTFitter plugin
// plugins rely on ROOT6 env to work
// utility tool execMacro.sh compiles and executes the macro, to use it do: ./exec/execMacro.sh exec/simple_fit.cc
// exec config expanding the coeffs into shapes and fitting

// correlated Cii, cPrk, cHel: ut cvv c1
// correlated cMnr, cMnk: dt cmm
// correlated Blk, Blr: cva c3
// correlated Blj, Blq: cav c123
// none: cmp (just Bln, but checking with all Bli guys)

#include "../src/EFTFitter.h"
#include "TF1.h"

int main() {
  using EFT = EFTFitter;

  // construction done with the key we want to treat as data, lambda, fit mode and stat mode (optionally sum of shape templates)
  EFT eft("data", 1., EFT::Fit::hybrid, EFT::Stat::xsec, 1.);
  const std::string inDir = "./wbern_0519/root/";

  // dummy rate as usual
  const std::array<double, 2> rate_zero = {0., 0.};

  std::vector<std::unique_ptr<TF1>> v_spinvar;
  v_spinvar.emplace_back( std::make_unique<TF1>("f_bli", "0.5 * (1. + ([0] * x))", -1., 1.) );
  v_spinvar.emplace_back( std::make_unique<TF1>("f_cii", "0.5 * (1. - ([0] * x)) * std::log(1. / std::abs(x))", -1., 1.) );
  v_spinvar.emplace_back( std::make_unique<TF1>("f_cPMij", "0.5 * (1. - (0.5 * [0] * x)) * std::acos(std::abs(x))", -1., 1.) );
  v_spinvar.emplace_back( std::make_unique<TF1>("f_cHel", "0.5 * (1. - ([0] * x))", -1., 1.) );

  const std::vector<std::string> v_operator = {"==op1==", "==op2=="};
  const std::vector<std::string> v_syst = {"==syst_sm==", "==syst_op1==", "==syst_op2=="};

  const std::map<std::vector<std::string>, std::pair<std::vector<double>, std::vector<int>>> m_range_var = {
    // main 2D results (take 1D from here too)
    {{"ut", "cvv"}, {{0.04, 0.08}, {10, 12, 13, 19}}},
    {{"ut", "c1"}, {{0.04, 0.7}, {10, 12, 13, 19}}},
    {{"dt", "cmm"}, {{0.04, 0.08}, {3, 4, 16, 18}}},

    // uncorrelated - for thesis 
    {{"cvv", "cav"}, {{0.06, 0.08}, {8, 9, 12, 19}}},
    {{"dt", "cmp"}, {{0.04, 0.02}, {4, 5, 16, 18}}},

    // cancelling pair for thesis
    {{"cvv", "c1"}, {{0.4, 2.}, {1, 10, 12, 19}}},

    // dopes because cba to split into 2 macros; ut is just a dummy
    // take all other 1D fits here
    {{"ut", "cav"}, {{0.04, 0.08}, {6, 7, 8, 9}}},
    {{"ut", "c123"}, {{0.04, 0.5}, {6, 7, 8, 9}}},
    {{"ut", "cva"}, {{0.04, 0.05}, {1, 3, 10, 15}}},
    {{"ut", "c3"}, {{0.04, 0.5}, {1, 3, 10, 15}}},
    {{"ut", "cmp"}, {{0.04, 0.02}, {4, 5, 8, 17}}}
  };

  /*
  const std::vector<std::string> v_hStr = {"b1k", "b2k", "b1r", "b2r", "b1n", "b2n", "b1j", "b2j", "b1q", "b2q",
                                           "ckk", "crr", "cnn", 
                                           "cP_rk", "cM_rk", "cP_nr", "cM_nr", "cP_nk", "cM_nk",
                                           "cHel",
                                           "cLab", "LL_dPhi"};
  */
  // variable index
  //const std::vector<int> v_variable = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

  for (uint iOp1 = 0; iOp1 < v_operator.size() - 1; ++iOp1) {
    for (uint iOp2 = iOp1 + 1; iOp2 < v_operator.size(); ++iOp2) {
      const std::string op1N = v_operator.at(iOp1), op2N = v_operator.at(iOp2);
      const double op1R = m_range_var.at(v_operator).first.at(0), op2R = m_range_var.at(v_operator).first.at(1);
      const std::vector<int> &v_variable = m_range_var.at(v_operator).second;
      const int binToIgnore = 1;

      const std::string op0S = v_syst.at(0), op1S = v_syst.at(1), op2S = v_syst.at(2);

      const std::string outDir = inDir + "../indep4_bug_itHi2_0726/shape/sm_" + op0S + "_" + op1N + "_" + op1S + "_" + op2N + "_" + op2S + "/";

      // emulate shape fits with arbitrary template selection
      // be very specific about the captures!
      auto f_emulate = [&v_variable, &rate_zero, &binToIgnore, &v_spinvar] (const auto &v_binC) {
        // convert the coeffs into shapes using the known formulas
        std::vector<std::array<double, 2>> v_tmpC = {rate_zero, rate_zero};
        TF1 *f_ref;

        for (int iC = 2; iC < v_binC.size(); ++iC) {
          if (!std::count(std::begin(v_variable), std::end(v_variable), iC - 2)) continue;

          // get the function ref
          if (iC - 2 < 10)
            f_ref = v_spinvar.at(0).get();
          else if (iC - 2 < 13)
            f_ref = v_spinvar.at(1).get();
          else if (iC - 2 < 19)
            f_ref = v_spinvar.at(2).get();
          else if (iC - 2 < 20)
            f_ref = v_spinvar.at(3).get();

          const bool is_cii = (std::string(f_ref->GetName()).find("f_cii") != std::string::npos);

          f_ref->SetParameter(0, v_binC.at(0).at(0) * v_binC.at(iC).at(0));
          const double integral = (is_cii) ? f_ref->Integral(-1., -DBL_MIN) + f_ref->Integral(DBL_MIN, 1.) : f_ref->Integral(-1., 1.);
          for (int iB = 0; iB < 6; ++iB) {
            if (iB == binToIgnore) continue;

            double min = -1. + (iB / 3.), max = -1. + ((iB + 1) / 3.);
            if (is_cii) {
              if (min == 0.)
                min = DBL_MIN;
              if (max == 0.)
                max = -DBL_MIN;
            }

            v_tmpC.push_back({f_ref->Integral(min, max) / integral, 0.});
          }
        }

        return v_tmpC;
      };
      eft.setHybridTransformation(std::function<std::vector<std::array<double, 2>> (const std::vector<std::array<double, 2>> &)> (f_emulate));

      // add the data
      auto f_shape_data = [&v_variable, &rate_zero, &binToIgnore] (const auto &v_binC) {
        // list of variables to include - ensure consistent with main fitter loop below
        const double shapeSum = 22.;

        const int nBin = v_binC.size();
        const int nBinEach = (nBin - 2) / int(shapeSum);

        // actual rate to be used in addHybridData()
        std::vector<std::array<double, 2>> v_tmpC = {{v_binC.at(0).at(0), v_binC.at(0).at(1)}, rate_zero};

        for (int iB = 2; iB < v_binC.size(); ++iB) {
          const int iVar = (iB - 2) / 6;
          if (!std::count(std::begin(v_variable), std::end(v_variable), iVar)) continue;
          if ((iB - 2) % nBinEach == binToIgnore) continue;

          // no shapeSum scaling unlike in MC since addHybridData doesn't normalize the binContent
          v_tmpC.push_back({v_binC.at(iB).at(0), v_binC.at(iB).at(1)});
        }

        return v_tmpC;
      };
      eft.addHybridData("./covariance_matrix/unfolded_data_190114.root", "TTbarSpinDensityMatrix/snake_spinCorr_shape_data_a", 
                        rate_zero, EFT::Stat::xsec, f_shape_data);

      eft.addRawInput(op1N + "_0--" + op2N + "_0", EFT::Sample::all, 
                      inDir + op1N + "_coeff.root", "sm_" + op0S, "", 1, rate_zero, EFT::Stat::xsec);
      eft.addRawInput(op1N + "_1--" + op2N + "_0", EFT::Sample::all, 
                      inDir + op1N + "_coeff.root", "sm_" + op0S + "_plus_" + op1N + "_" + op1S, "", 1, rate_zero, EFT::Stat::xsec);
      eft.addRawInput(op1N + "_-1--" + op2N + "_0", EFT::Sample::all, 
                      inDir + op1N + "_coeff.root", "sm_" + op0S + "_minus_" + op1N + "_" + op1S, "", 1, rate_zero, EFT::Stat::xsec);
      eft.addRawInput(op1N + "_0--" + op2N + "_1", EFT::Sample::all, 
                      inDir + op2N + "_coeff.root", "sm_" + op0S + "_plus_" + op2N + "_" + op2S, "", 1, rate_zero, EFT::Stat::xsec);
      eft.addRawInput(op1N + "_0--" + op2N + "_-1", EFT::Sample::all, 
                      inDir + op2N + "_coeff.root", "sm_" + op0S + "_minus_" + op2N + "_" + op2S, "", 1, rate_zero, EFT::Stat::xsec);
      eft.addRawInput(op1N + "_1--" + op2N + "_1", EFT::Sample::all, 
                      inDir + op1N + "_" + op2N + "_coeff.root", 
                      "sm_" + op0S + "_plus_" + op1N + "_" + op1S + "_plus_" + op2N + "_" + op2S, "", 1, rate_zero, EFT::Stat::xsec);

      // prepare the base for interpolation
      eft.prepareInterpolationBase();

      // assign as data a particular key of choice
      //eft.assignAsData("ut_0--dt_0", EFT::Sample::all, false);

      std::vector<std::tuple<std::string, EFT::Sample, std::string>> vt_keySampleLegend;
      vt_keySampleLegend.push_back({op1N + "_1--" + op2N + "_0", EFT::Sample::linear, op1N + " 1"});
      vt_keySampleLegend.push_back({op1N + "_0--" + op2N + "_1", EFT::Sample::linear, op2N + " 1"});
      vt_keySampleLegend.push_back({op1N + "_1--" + op2N + "_1", EFT::Sample::linear, op1N + " " + op2N + " 1"});
      vt_keySampleLegend.push_back({op1N + "_0--" + op2N + "_0", EFT::Sample::linear, "SM"});
      vt_keySampleLegend.push_back({"data", EFT::Sample::all, "Data"});

      eft.drawHistogram(vt_keySampleLegend, 
                        outDir + op1N + "_" + op2N + "_shape", 
                        "Fraction", "Index", -0.4999, 0.4999, 0.0001, 1.9999,
                        false, "", false, "none");

      // grab the cov matrix
      std::vector<std::array<int, 2>> covMat_binRange;
      for (auto &var : v_variable)
        covMat_binRange.push_back({((var + 1) * 5) - 4, ((var + 1) * 5)});
      eft.readCovMatRoot("finalCov", "./covariance_matrix/covmat_190114.root", "TotalStatSyst_shape_a_drop_bin_1", covMat_binRange);

      // draw and store all available covmats
      eft.drawCovMat(outDir);

      // make the range to interpolate over; in this case [min, max: step]
      const double eps = 0.0001;
      const std::vector<EFT::Sample> v_sample = {EFT::Sample::linear};

      eft.listKeyToFit({ {op1N, {0.}}, {op2N, makeInterval(-op2R, op2R, op2R / 10000.)} });
      eft.computeFitChi2(v_sample);
      eft.draw1DChi2({ {op2N, {op2N, {/* op range in min, max */}, {0., 9.999}, {-op2R + eps, op2R - eps} }} }, outDir, v_sample);
      eft.clearContent(1);

      if (op1N == "ut" and op2N != "cvv" and op2N != "c1")
        eft.clearContent();
      else {
        eft.listKeyToFit({ {op1N, makeInterval(-op1R, op1R, op1R / 10000.)}, {op2N, {0.}} });
        eft.computeFitChi2(v_sample);
        eft.draw1DChi2({ {op1N, {op1N, {/* op range in min, max */}, {0., 9.999}, {-op1R + eps, op1R - eps} }} }, outDir, v_sample);
        eft.clearContent(1);

        eft.listKeyToFit({ {op1N, makeInterval(-op1R, op1R, op1R / 1000.)}, {op2N, makeInterval(-op2R, op2R, op2R / 1000.)} });
        eft.computeFitChi2(v_sample);
        eft.draw2DChi2({{ {op1N, op2N}, {{{op1N, {-op1R + eps, op1R - eps}}, {op2N, {-op2R + eps, op2R - eps}}}} }}, outDir, v_sample, 0.125);
        eft.clearContent();
      }
    }
  }

  return 0;
}

