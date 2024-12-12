// example code using the EFTFitter plugin
// plugins rely on ROOT6 env to work
// utility tool execMacro.sh compiles and executes the macro, to use it do: ./exec/execMacro.sh exec/this_file.cc

// a sane version of evolution plot that does not require too much manual work

#include <filesystem>

#include "../src/EFTFitter.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h" 
#include "TSystem.h"

std::array<std::unique_ptr<TGraphAsymmErrors>, 2> fitResult(const std::string &op = "", const std::string &fileName = "",
                                                            const int &iIter = 0, const int &iVar = 0, 
                                                            const bool useAll = true) {
  using TG = TGraphAsymmErrors;
  std::array<std::unique_ptr<TG>, 2> a_graph = {nullptr, nullptr};
  if (op == "" or fileName == "")
    return a_graph;

  const std::string samp = (useAll) ? "all" : "linear", tag = toStr(iIter) + "_" + toStr(iVar);
  auto file = std::unique_ptr<TFile>(TFile::Open( fileName.c_str() ));

  if (!file->GetListOfKeys()->Contains((op + "_sigma2_" + samp).c_str()))
    return a_graph;

  // read the real best fit and sigmas
  a_graph.at(0) = std::unique_ptr<TG>(dynamic_cast<TG *>(( file->Get( (op + "_sigma1_" + samp).c_str() ))->Clone( (tag + "_s1").c_str() ) ));
  a_graph.at(1) = std::unique_ptr<TG>(dynamic_cast<TG *>(( file->Get( (op + "_sigma2_" + samp).c_str() ))->Clone( (tag + "_s2").c_str() ) ));

  return a_graph;
}

int main(/*int argc, char** argv*/) {
  using EFT = EFTFitter;

  // common flags for an evolution test
  // fit with data or SM, check evolution of all or linear part, 2 sigma or 1 sigma interval width as figure of merit 
  const bool useData = true, useAll = false, useSig2 = true;

  // construction done with the key we want to treat as data, lambda, fit mode and stat mode (optionally sum of shape templates)
  EFT eft("data", 1., EFT::Fit::hybrid, EFT::Stat::xsec, 1.);
  const std::string dName = "TTbarSpinDensityMatrix/";
  const std::string /*hName = "LL_dPhi",*/ sName = dName + "sumWgt_noCut";
  const std::string inDir = "./nlo_ctG_root_density/";
  const std::string sufDat = "_shape_data_a", sufSim = "_absolute_part_cut_0";
  const double k_nnlo_nlo = 1.220251456; // see notes at bottom
  const std::string outDir = "./nlo_ctG_0115/nominal_toppt_default/shape_data/evolution/";

  // make the range to interpolate over; in this case [min, max: step]
  const std::vector<double> v_opPoint = makeInterval(-5., 5., 0.0001);
  std::vector<EFT::Sample> v_sample; 
  if (useAll)
    v_sample.push_back(EFT::Sample::all);
  else
    v_sample.push_back(EFT::Sample::linear);
  const std::string sample = (useAll) ? "all" : "linear";

  const std::vector<std::string> v_hStr = {"b1k", "b2k", "b1r", "b2r", "b1n", "b2n", "b1j", "b2j", "b1q", "b2q",
                                           "ckk", "crr", "cnn", 
                                           "cP_rk", "cM_rk", "cP_nr", "cM_nr", "cP_nk", "cM_nk",
                                           "cHel"/*,
                                           "cLab", "LL_dPhi"*/};

  const std::vector<std::array<int, 2>> v_snake = {{1, 6}, {7, 12}, {13, 18}, {19, 24}, {25, 30}, {31, 36}, {37, 42}, 
                                             {43, 48}, {49, 54}, {55, 60}, {61, 66}, {67, 72}, {73, 78}, {79, 84}, 
                                             {85, 90}, {91, 96}, {97, 102}, {103, 108}, {109, 114}, 
                                             {115, 120}/*, {121, 126}, {127, 132}*/};

  // save the rates in case they're used
  //const std::array<double, 2> rate_17_001 = {803., 33.}, rate_18_006 = {836.925, 51.685};
  const std::array<double, 2> rate_data = {0., 0.}, rate_zero = {0., 0.};

  // variable index
  auto evoFile = std::make_unique<TFile>((outDir + "ctG_evolution.root").c_str(), "recreate");
  std::vector<int> v_all_var = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19/*, 20, 21*/};
  std::vector<int> v_fit_var;

  for (int iVar = 0; iVar < v_hStr.size(); ++iVar) {
    std::cout << "EFTFitter: iteration " << iVar << " starting..." << std::endl;
    std::vector<std::array<std::unique_ptr<TGraphAsymmErrors>, 2>> v_current_result;

    for (int iFit = 0; iFit < v_hStr.size(); ++iFit) {
      if ( std::count(std::begin(v_fit_var), std::end(v_fit_var), iFit) ) {
        v_current_result.emplace_back(fitResult("", ""));
        continue;
      }

      std::vector<int> v_iEB(v_fit_var);
      v_iEB.push_back( iFit );
      std::sort(std::begin(v_iEB), std::end(v_iEB));

      std::vector<std::array<int, 2>> v_endbin;
      for (int iEB = 0; iEB < v_iEB.size(); ++iEB)
        v_endbin.push_back( {((v_iEB.at(iEB) + 1) * 5) - 4, ((v_iEB.at(iEB) + 1) * 5)} ); // hybrid shape

      std::cout << "EFTFitter: starting fit on variable " << v_hStr.at(iFit) << " in iteration " << iVar << std::endl;

      // emulate absolute/shape fits with arbitrary template selection (with snake input and normalization)
      // be very specific about the captures!
      auto f_emulate = [&v_iEB, &rate_data, &rate_zero] (const auto &v_binC) {
        // variables that are included (ensure it's !count in loop)
        const std::vector<int> v_var = v_iEB;
        const double shapeSum = 22.;

        // just use mc rate with data err - don't care about yields in xsec mode
        const std::array<double, 2> rate_mc = {v_binC.at(0).at(0) / shapeSum, rate_data.at(1)};
        std::vector<std::array<double, 2>> v_tmpC = {(rate_data.at(0) > 0.) ? rate_mc : rate_zero, rate_zero};

        // shape fit
        // ensure binToIgnore matches the _drop_bin_N in matrix name - hybrid assumes matrix is good to go
        const int nBin = v_binC.size();
        const int binToIgnore = 1, nBinEach = (nBin - 2) / int(shapeSum);

        for (int iB = 2; iB < v_binC.size(); ++iB) {
          const int iFit = (iB - 2) / 6;
          if (!std::count(std::begin(v_var), std::end(v_var), iFit)) continue;
          if ((iB - 2) % nBinEach == binToIgnore) continue;

          v_tmpC.push_back({shapeSum * v_binC.at(iB).at(0), shapeSum * v_binC.at(iB).at(1)});
        }

        return v_tmpC;
      };

      eft.setShapeSum(v_endbin.size());
      eft.setHybridTransformation(std::function<std::vector<std::array<double, 2>> (const std::vector<std::array<double, 2>> &)> (f_emulate));

      // add the input templates
      // please ensure there is exactly 1 input with dataName as in ctor (real data gets 0 xsec)
      // and MC following the syntax op1_val1--op2_val2-- ... --opN_valN
      // NLO (taken off James' EFT_Fitter)
      // SM
      eft.addRawInput("ctG_0", EFT::Sample::all, inDir + "nlo_ctG_p0_nominal_toppt_default_shape_concatenated.root", 
                      dName + "snake_spinCorr" + sufSim, sName, 1, {22. * k_nnlo_nlo * 681.63, 0.}, EFT::Stat::count);

      // and EFT
      eft.addRawInput("ctG_2", EFT::Sample::all, inDir + "nlo_ctG_p2_nominal_toppt_default_shape_concatenated.root", 
                      dName + "snake_spinCorr" + sufSim, sName, 1, {22. * k_nnlo_nlo * 1260.07, 0.}, EFT::Stat::count);

      eft.addRawInput("ctG_-2", EFT::Sample::all, inDir + "nlo_ctG_m2_nominal_toppt_default_shape_concatenated.root", 
                      dName + "snake_spinCorr" + sufSim, sName, 1, {22. * k_nnlo_nlo * 400.97, 0.}, EFT::Stat::count);

      // prepare the base for interpolation
      eft.prepareInterpolationBase();

      // data
      if (useData) {
        auto f_shape_data = [&v_iEB, &rate_zero] (const auto &v_binC) {
          // list of variables to include - ensure consistent with main fitter loop below
          std::vector<int> v_var = v_iEB;
          const double shapeSum = 22.;

          const int nBin = v_binC.size();
          const int binToIgnore = 1, nBinEach = (nBin - 2) / int(shapeSum);

          // actual rate to be used in addHybridData()
          std::vector<std::array<double, 2>> v_tmpC = {{v_binC.at(0).at(0), v_binC.at(0).at(1)}, rate_zero};

          for (int iB = 2; iB < v_binC.size(); ++iB) {
            const int iFit = (iB - 2) / 6;
            if (!std::count(std::begin(v_var), std::end(v_var), iFit)) continue;
            if ((iB - 2) % nBinEach == binToIgnore) continue;

            // no shapeSum scaling unlike in MC since addHybridData doesn't normalize the binContent
            v_tmpC.push_back({v_binC.at(iB).at(0), v_binC.at(iB).at(1)});
          }

          return v_tmpC;
        };
        eft.addHybridData("./HEPData-concatenated-histogram.root", "concatenated_diff_cross_section;1", 
                          rate_data, EFT::Stat::xsec, f_shape_data);
      }
      else
        eft.assignAsData("ctG_0", EFT::Sample::all, false);

      // grab the stat correlation matrix given by Jacob
      const std::vector<std::array<int, 2>> covMat_binRange = v_endbin;
      eft.readCovMatRoot("finalcov", "./covariance_matrix/statistical_covariance_matrix.root", "covariance_matrix", covMat_binRange);

      eft.listKeyToFit({ {"ctG", v_opPoint} });
      eft.computeFitChi2(v_sample);

      eft.draw1DChi2({ {"ctG", { "ctG", {/* op range in min, max */}, {0., 9.999}, {-4.999, 9.999} }} }, outDir, v_sample);

      eft.clearContent();

      v_current_result.emplace_back(fitResult("ctG", outDir + "ctG_dChi2.root", iVar, iFit, useAll));
      if (v_current_result.back().at(0) == nullptr)
        std::cout << "EFTFitter: fit on variable " << v_hStr.at(iFit) << " in iteration " << iVar << " doesn't result in a constraint!" << std::endl;
      else 
        std::cout << "EFTFitter: fit on variable " << v_hStr.at(iFit) << " in iteration " << iVar << " saved for evolution scan..." << std::endl;

      gSystem->Exec( ("rm " + outDir + "ctG_dChi2.root").c_str() );
    }

    // first find the best variable in a given iter by minimum sig2 or sig1 width
    int iMin = -999;
    double width = 9999.;
    for (int iRes = 0; iRes < v_current_result.size(); ++iRes) {
      if (v_current_result.at(iRes).at(0) == nullptr) continue;
      auto &graph = (useSig2) ? v_current_result.at(iRes).at(1) : v_current_result.at(iRes).at(0);

      if (graph->GetErrorXlow(0) + graph->GetErrorXhigh(0) < width) {
        width = graph->GetErrorXlow(0) + graph->GetErrorXhigh(0);
        iMin = iRes;
      }
    }

    v_fit_var.push_back(iMin);
    std::cout << "EFTFitter: variable " << v_hStr.at(iMin) << " with width " << width << " is the best variable in iteration " << iVar << std::endl;

    evoFile->cd();
    for (int iRes = 0; iRes < v_current_result.size(); ++iRes) {
      if (v_current_result.at(iRes).at(0) == nullptr) continue;

      if (iRes == iMin) {
        v_current_result.at(iRes).at(0)->SetName( ("iter_" + toStr(iVar) + "_best_" + sample + "_" + toStr(iRes) + "_sigma1").c_str() );
        v_current_result.at(iRes).at(0)->Write( ("iter_" + toStr(iVar) + "_best_" + sample + "_" + toStr(iRes) + "_sigma1").c_str() );

        v_current_result.at(iRes).at(1)->SetName( ("iter_" + toStr(iVar) + "_best_" + sample + "_" + toStr(iRes) + "_sigma2").c_str() );
        v_current_result.at(iRes).at(1)->Write( ("iter_" + toStr(iVar) + "_best_" + sample + "_" + toStr(iRes) + "_sigma2").c_str() );
      }
      else {
        v_current_result.at(iRes).at(0)->SetName( ("iter_" + toStr(iVar) + "_fit_" + sample + "_" + toStr(iRes) + "_sigma1").c_str() );
        v_current_result.at(iRes).at(0)->Write( ("iter_" + toStr(iVar) + "_fit_" + sample + "_" + toStr(iRes) + "_sigma1").c_str() );

        v_current_result.at(iRes).at(1)->SetName( ("iter_" + toStr(iVar) + "_fit_" + sample + "_" + toStr(iRes) + "_sigma2").c_str() );
        v_current_result.at(iRes).at(1)->Write( ("iter_" + toStr(iVar) + "_fit_" + sample + "_" + toStr(iRes) + "_sigma2").c_str() );
      }
    }

    std::cout << "EFTFitter: iteration " << iVar << " completed!" << std::endl;
  }

  return 0;
}

/***
    average from log xsecs (pb)
    ctg 2 1260.07 pb
    ctg 0 681.63 pb
    ctg -2 400.97

    tt > 2l BR
    old 0.0493877 (off James' ctg log)

    tt xsec
    NLO 681.63 pb (James setup)
    NNLO 831.76 pb (https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO)

    k factor (= NNLO / NLO SM)
    1.220251456

    theory
    shapes by scales 1/2, 2 variation in uR, uF
    rate by 0.939 and 1.058 multiplied to k-factor
***/
