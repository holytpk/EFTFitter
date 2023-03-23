// -*- C++ -*-
// author: afiq anuar
// short: please refer to header for information

#include "EFTFitter.h"
#include "PlotUtil.h"
#include "TTree.h"

/***
 * EFTFitter source (public)
 ***/



// ctors
EFTFitter::EFTFitter(const std::string &dataName_, const double &eftLambda_, 
                     const Fit fitMode_, const Stat statMode_, const double &shapeSum_) :
  dataName(dataName_),
  hasData(false),
  eftLambda(eftLambda_),
  fitMode(fitMode_),
  statMode(statMode_),
  shapeSum(shapeSum_)
{
  if (fitMode != Fit::absolute and fitMode != Fit::shape and fitMode != Fit::hybrid)
    throw std::invalid_argument( "Invalid fit mode: which fit mode to perform, must be (from namespace EFTFitter) "
                                 "'Fit::absolute' (count/xsec as-is within the template), 'Fit::shape' (shape ie event fraction) " 
                                 "or 'Fit::hybrid' (separately fit on count/xsec and then shape)\n" 
                                 "Note: 'Fit::hybrid' ignores the integral of all templates given to addRawInput(). " 
                                 "The integral is set to the value passed to normalizedSum (can be count/xsec based on sumStat)" );

  if (statMode != Stat::count and statMode != Stat::xsec)
    throw std::invalid_argument( "Invalid event mode: how to count the events, must be (from namespace EFTFitter)"
                                 " 'Stat::count' (use raw event count) or 'Stat::xsec' (use xsec instead)" );

  std::cout << "EFTFitter is initialized with mode " << fitMode << 
    ", event statistics counted as " << statMode << "...\n" << 
    "Expecting as data a template with key " << dataName << "..." << std::endl << std::endl;

  // some ROOT-related stuff
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2(true);
}



void EFTFitter::setShapeSum(const double &shapeSum_)
{
  shapeSum = shapeSum_;
}



void EFTFitter::addRawInput(const std::string &keyName, const Sample sampleType, 
                            const std::string &fileName, const std::string &histName, const std::string &sumName,
                            const int nRebin, const std::array<double, 2> &normalizedSum, const Stat sumStat,
                            const bool addIfPresent) 
{
  if (!checkOpSet(keyName, sampleType))
    throw std::logic_error( "Added key " + keyName + " :: " + toStr(sampleType) + " is not consistent with available input!!" );

  if (fitMode == Fit::hybrid) {
    if (keyName == dataName)
      throw std::invalid_argument( "In Fit::hybrid, data templates should not be added using this method. Please use addHybridData() instead!!" );

    if (!f_hybridTransform)
      throw std::logic_error( "Mode is Fit::hybrid, but f_hybridTransform is undefined. Please provide a transformation function to setHybridTransformation()!!" );
  }

  std::cout << "Adding input with key " << keyName << " of type " << sampleType << "..." << std::endl << std::endl;

  // open the input file, get the hist
  auto file = std::unique_ptr<TFile>(TFile::Open( fileName.c_str() ));
  auto hist = std::unique_ptr<TH1D>(dynamic_cast<TH1D *>(( file->Get( histName.c_str() ))->Clone(fixKeyFormat(keyName).c_str()) ));

  // rebin iff nRebin is compatible
  if (nRebin > 1 and hist->GetNbinsX() % nRebin == 0)
    hist->Rebin(nRebin);

  // define the binning of the templates
  if (v_rawBin.empty()) {
    if (fitMode != Fit::hybrid)
      v_rawBin = FitUtil::extractBin( hist.get() );
    else {
      // Fit::hybrid, so just prepare a dummy binContent to make a binning of just integer interval
      // all that matters is that interval length is what one gets out of f_hybridTransform
      std::vector<std::array<double, 2>> v_pre = FitUtil::extractContentError(hist.get());
      v_pre.push_back({0., 0.});
      v_pre.push_back({0., 0.});

      std::vector<std::array<double, 2>> v_pos = f_hybridTransform(v_pre);
      v_rawBin = makeInterval(0., double(v_pos.size() - 2), 1.);
    }
  }
  else {
    if (!checkInputBin(hist)) {
      std::cout << "Expected binning based on first included template:" << std::endl;
      printAll(v_rawBin);

      std::cout << "Binning of currently provided template:" << std::endl;
      printAll(FitUtil::extractBin( hist.get() ));

      throw std::logic_error( "Given raw input with key " + keyName + " :: " + toStr(sampleType) + 
                              " is binned differently than others!!" );
    }
  }

  // sum of weight histogram for normalization
  std::unique_ptr<TH1D> sumHist(nullptr);

  // assign an xsec for intLumi normalization and if needed interpolation
  const double &xsec = normalizedSum.at(0);
  if (fitMode != Fit::hybrid and xsec != 0. and sumStat == Stat::count) {
    sumHist = std::unique_ptr<TH1D>(dynamic_cast<TH1D *>( (file->Get( sumName.c_str() ))->Clone() ));
    hist->Scale( FitUtil::intLumi * (xsec / sumHist->GetBinContent(1)) );
  }

  // get the integral, error, and normalized bin content
  // "" to get the sum(content), while "width" to get sum(content * width) 
  double hsum, herr;
  hsum = hist->IntegralAndError(-1, -1, herr, "");
  hist->Scale( 1. / hist->Integral() );

  // make the initial content vector - first is xsec, second count
  std::vector<std::array<double, 2>> v_binContent;
  if (sumStat == Stat::count) {
    if (xsec != 0.) {
      // get the stat error of xsec - keep relative stat unc of integral
      const double xerr = xsec * (herr / hsum);

      v_binContent = { {xsec, xerr}, {hsum, herr} };
    }
    else 
      v_binContent = { {hsum / FitUtil::intLumi, herr / FitUtil::intLumi}, {hsum, herr} };
  }
  else if (sumStat == Stat::xsec)
    v_binContent = { {hsum, herr}, {hsum * FitUtil::intLumi, herr * FitUtil::intLumi} };

  const std::vector<std::array<double, 2>> v_binNorm = FitUtil::extractContentError(hist.get());
  v_binContent.insert(std::end(v_binContent), std::begin(v_binNorm), std::end(v_binNorm));

  // if element is already present, erase or add to it as requested
  const std::pair<std::string, Sample> keyPair = {fixKeyFormat(keyName), sampleType};
  if (m_binContent.count(keyPair)) {
    if (!addIfPresent)
      std::cout << "Key " << keyName << " is already present in map and will be overwritten..." << std::endl;
    else {
      std::cout << "Key " << keyName << " is already present in map, adding this input to it..." << std::endl;

      // make copy of the binContent and update the current one
      const std::vector<std::array<double, 2>> v_tmpC = v_binContent;
      for (int iB = 0; iB < v_tmpC.size(); ++iB) {
        if (iB < 2) {
          // xsec, integral are just sum
          v_binContent.at(iB).at(0) = v_tmpC.at(iB).at(0) + m_binContent.at(keyPair).at(iB).at(0);

          // integral error is just quad sum
          if (iB == 1) {
            v_binContent.at(iB).at(1) = (v_tmpC.at(iB).at(1) * v_tmpC.at(iB).at(1)) + 
              (m_binContent.at(keyPair).at(iB).at(1) * m_binContent.at(keyPair).at(iB).at(1));
            v_binContent.at(iB).at(1) = std::sqrt(v_binContent.at(iB).at(1));
          }
        }
        else {
          // in each bin first make the sum of content as above
          v_binContent.at(iB).at(0) = v_tmpC.at(1).at(0) * v_tmpC.at(iB).at(0);
          v_binContent.at(iB).at(0) += m_binContent.at(keyPair).at(1).at(0) * m_binContent.at(keyPair).at(iB).at(0);

          // and also error
          v_binContent.at(iB).at(1) = (v_tmpC.at(1).at(0) * v_tmpC.at(iB).at(1)) * (v_tmpC.at(1).at(0) * v_tmpC.at(iB).at(1));
          v_binContent.at(iB).at(1) += (m_binContent.at(keyPair).at(1).at(0) * m_binContent.at(keyPair).at(iB).at(1)) * 
            (m_binContent.at(keyPair).at(1).at(0) * m_binContent.at(keyPair).at(iB).at(1));
          v_binContent.at(iB).at(1) = std::sqrt(v_binContent.at(iB).at(1));

          // finally normalize them
          v_binContent.at(iB).at(0) /= v_binContent.at(1).at(0);
          v_binContent.at(iB).at(1) /= v_binContent.at(1).at(0);
        }
      }

      // finally xsec error keeping the same proportionality as updated integral
      v_binContent.at(0).at(1) = v_binContent.at(0).at(0) * (v_binContent.at(1).at(1) / v_binContent.at(1).at(0));
    }

    // we're done making the sum, erase the current
    m_binContent.erase(keyPair);
  }

  // and dump the new one in
  m_binContent.insert({keyPair, v_binContent});
}



void EFTFitter::setHybridTransformation(const std::function<std::vector<std::array<double, 2>> (const std::vector<std::array<double, 2>> &)> &func)
{
  if (fitMode != Fit::hybrid)
    throw std::logic_error( "setHybridTransformation() method is only usable in Fit::hybrid!!" );

  if (f_hybridTransform)
    std::cout << "Hybrid transformation function already available, overwriting it..." << std::endl;
  f_hybridTransform = func;
}



void EFTFitter::addHybridData(const std::string &fileName, const std::string &histName, 
                              const std::array<double, 2> &normalizedSum, const Stat sumStat, 
                              const std::function<std::vector<std::array<double, 2>> (const std::vector<std::array<double, 2>> &)> &func)
{
  if (fitMode != Fit::hybrid)
    throw std::logic_error( "addHybridData() method is only usable in Fit::hybrid!!" );

  std::cout << "Adding hybrid data..." << std::endl;

  const std::pair<std::string, Sample> dataKey = {dataName, Sample::all};
  if (m_binContent.count(dataKey)) {
    std::cout << "Data already present in map and will be overwritten..." << std::endl;
    m_binContent.erase(dataKey);
  }

  // open the input file, get the hist
  auto file = std::unique_ptr<TFile>(TFile::Open( fileName.c_str() ));
  auto hist = std::unique_ptr<TH1D>(dynamic_cast<TH1D *>(( file->Get( histName.c_str() ))->Clone(dataName.c_str()) ));

  // assign the xsec of the data
  const double &xsec = normalizedSum.at(0);

  // make the initial content vector - first is xsec, second count
  std::vector<std::array<double, 2>> v_binContent;
  if (sumStat == Stat::count) {
    // for counts just assume Poissonian stats if the user didn't put the error in
    const double xerr = (normalizedSum.at(1) > 0.) ? normalizedSum.at(1) : std::sqrt(xsec);

    v_binContent = { {xsec / FitUtil::intLumi, xerr / FitUtil::intLumi}, {xsec, xerr} };
  }
  else if (sumStat == Stat::xsec) {
    const double xerr = (normalizedSum.at(1) > 0.) ? normalizedSum.at(1) : -1.;

    if (xsec != 0. and xerr < 0.)
      throw std::invalid_argument( "Invalid error of normalizedSum not allowed when using Fit::hybrid "
                                   "with a valid normalizedSum value of Stat::xsec; please assign a valid error!!");

    v_binContent = { {xsec, xerr}, {xsec * FitUtil::intLumi, xerr * FitUtil::intLumi} };
  }

  const std::vector<std::array<double, 2>> v_binNorm = FitUtil::extractContentError(hist.get());
  v_binContent.insert(std::end(v_binContent), std::begin(v_binNorm), std::end(v_binNorm));

  // and in it goes - mark that it has been done
  if (func == nullptr)
    m_binContent.insert({dataKey, v_binContent});
  else
    m_binContent.insert({dataKey, func(v_binContent)});
  hasData = true;

  std::cout << std::endl;
}



void EFTFitter::drawHistogram(const std::vector< std::tuple<std::string, Sample, std::string> > &vt_keySampleLegend,
                              const std::string &plotName, const std::string &yLabel, const std::string &xLabel,
                              const double &histMin, const double &histMax, const double &ratioMin, const double &ratioMax,
                              const bool drawLogY, const std::string &legHeader, 
                              const bool divideBinWidth, const std::string &ratioMode)
{
  if (!hasData)
    throw std::logic_error( "Data input must already be available for this method to be called!!" );

  if (!std::count_if(std::begin(vt_keySampleLegend), std::end(vt_keySampleLegend), [&] (const auto &ksl) 
                     {return (std::get<0>(ksl) == dataName and std::get<1>(ksl) == Sample::all);}))
    throw std::invalid_argument( "This method requires the data histogram be plotted... :((" );

  if (ratioMode != "simple" and ratioMode != "covariance" and ratioMode != "none")
    throw std::logic_error( "Ratio mode not understood: only simple, covariance or none are allowed!" );
  const bool noRatio = (ratioMode == "none"), useCov = (ratioMode == "covariance");
  if (useCov and !m_covMat.count("finalcov"))
    throw std::logic_error( "Final covariance matrix must be available to use covariance ratio mode!!" );

  std::cout << "Drawing the requested key-sample templates " << "..." << std::endl;
  std::for_each(std::begin(vt_keySampleLegend), std::end(vt_keySampleLegend), [] (const auto &ksl) {
      std::cout << std::pair<std::string, Sample>(std::get<0>(ksl), std::get<1>(ksl)) << std::endl;
    });
  std::cout << std::endl;

  // let's face it, you CAN'T write an automatic plotter that always looks right
  // so just be prepared to tune this some...
  // first prepare the color map for our hist
  const int nPlot = std::count_if(std::begin(vt_keySampleLegend), std::end(vt_keySampleLegend), [&] (const auto &ksl) 
                                  {return std::get<0>(ksl) != dataName;});
  gStyle->SetPalette( kRainBow ); // just use the palette, has 255 colors
  std::vector<int> v_kColor;
  for (int iP = 0; iP < nPlot; ++iP)
    v_kColor.push_back(TColor::GetColorPalette( 15. + std::floor(224. * (iP * (1. / nPlot))) ));
  if (v_kColor.empty()) v_kColor.push_back(kBlack);

  // then we make the hist and ratio = hist / data and dress them up; also, legend
  auto iColor = std::begin(v_kColor);
  auto leg = std::make_unique<TLegend>();
  std::map< std::pair<std::string, Sample>, std::unique_ptr<TH1D> > m_hist, m_ratio;

  // first make the data hists and dress em up
  // cant use the init-list version (it implicitly calls copy ctor?), must use make_pair...
  auto hdat = convertContentToHist(dataName, Sample::all, divideBinWidth);
  auto hdra = std::unique_ptr<TH1D>(dynamic_cast<TH1D *>( hdat->Clone(("r_" + dataName).c_str()) ));
  m_hist.insert(std::make_pair( std::make_pair(dataName, Sample::all), std::move(hdat) ));
  m_ratio.insert(std::make_pair( std::make_pair(dataName, Sample::all), std::move(hdra) ));
  m_ratio.at({dataName, Sample::all})->Reset();
  if (useCov)
    m_ratio.at({dataName, Sample::all})->Add(m_hist.at({dataName, Sample::all}).get(), m_hist.at({dataName, Sample::all}).get(), 1., -1.);
  else
    m_ratio.at({dataName, Sample::all})->Divide(m_hist.at({dataName, Sample::all}).get(), m_hist.at({dataName, Sample::all}).get());

  // and then the MCs
  for (const auto &ksl : vt_keySampleLegend) {
    const std::string key = fixKeyFormat(std::get<0>(ksl));
    const Sample &samp = std::get<1>(ksl);

    if (key == dataName) continue;

    auto htmp = convertContentToHist(key, samp, divideBinWidth);
    auto hrat = std::unique_ptr<TH1D>(dynamic_cast<TH1D *>( htmp->Clone(("r_" + key + "_" + toStr(samp)).c_str()) ));
    hrat->Reset();

    if (useCov) {
      hrat->Add(htmp.get(), m_hist.at({dataName, Sample::all}).get(), 1., -1.);
      TMatrixD covmat(m_covMat.at("finalcov"));

      // for shape, drop all non-diagonals (inaccurate, but for visualization it's probably alright)
      if (fitMode == Fit::shape) {
        for (int iR = 0; iR < covmat.GetNrows(); ++iR) {
          for (int iC = 0; iC < covmat.GetNrows(); ++iC)
            if (iR != iC) covmat(iR, iC) = 0.;
        }
      }

      covmat.Invert();

      for (int iB = 1; iB <= hrat->GetNbinsX(); ++iB) {
        hrat->SetBinContent(iB, hrat->GetBinContent(iB) * std::sqrt(covmat(iB - 1, iB - 1)));
        hrat->SetBinError(iB, hrat->GetBinError(iB) * std::sqrt(covmat(iB - 1, iB - 1)));
      }
    }
    else
      hrat->Divide(htmp.get(), m_hist.at({dataName, Sample::all}).get());

    m_hist.insert(std::make_pair( std::make_pair(key, samp), std::move(htmp) ));
    m_ratio.insert(std::make_pair( std::make_pair(key, samp), std::move(hrat) ));

    stylePlot(m_hist.at({key, samp}).get(), *iColor, 1., 0, 21 + std::distance(std::begin(v_kColor), iColor), 1.5, 1, 2);
    stylePlot(m_ratio.at({key, samp}).get(), *iColor, 1., 0, 21 + std::distance(std::begin(v_kColor), iColor), 1.5, 1, 2);
    leg->AddEntry(m_hist.at({key, samp}).get(), std::get<2>(ksl).c_str(), "lp");

    iColor = std::next(iColor);
  }

  stylePlot(m_hist.at({dataName, Sample::all}).get(), kBlack, 1., 0, 20, 1.5, 1, 2);
  stylePlot(m_ratio.at({dataName, Sample::all}).get(), kBlack, 1., 0, 20, 1.5, 1, 2);

  // just to ensure data always comes last in legend
  const auto idat = std::find_if(std::begin(vt_keySampleLegend), std::end(vt_keySampleLegend), 
                                 [&] (const auto &ksl) {return (std::get<0>(ksl) == dataName and std::get<1>(ksl) == Sample::all);});
  const std::string ldat = std::get<2>(*idat);
  leg->AddEntry(m_hist.at({dataName, Sample::all}).get(), ldat.c_str(), "p");

  // the tuning pit; put whatever that "looks right" for axes label size etc
  for (const auto &ksl : vt_keySampleLegend) {
    const std::string key = fixKeyFormat(std::get<0>(ksl));
    const Sample &samp = std::get<1>(ksl);
    const std::string rLabel = (useCov) ? "#frac{MC - " + ldat + "}{Error}" : "MC / " + ldat;
    const double rOffset = (useCov) ? 0.37 : 0.37, rTitle = (useCov) ? 0.108 : 0.121;

    if (noRatio) {
      axisPlot(m_hist.at({key, samp}).get(), 
               histMin, histMax, 505, yLabel, 0.043, 1.21, 0.037, 0., 0., 505, xLabel, 0.043, 1.11, 0.037);
      axisPlot(m_ratio.at({key, samp}).get(), 
               ratioMin, ratioMax, 505, rLabel, 0.121, 0.41, 0.093, 0., 0., 505, xLabel, 0.131, 0.71, 0.107);
    }
    else {
      axisPlot(m_hist.at({key, samp}).get(), 
               histMin, histMax, 505, yLabel, 0.059, 0.73, 0.047, 0., 0., 505, xLabel, 0.037, 1.15, 0.033);
      axisPlot(m_ratio.at({key, samp}).get(), 
               ratioMin, ratioMax, 505, rLabel, rTitle, rOffset, 0.093, 0., 0., 505, xLabel, 0.131, 0.71, 0.107);
    }
  }

  // ok, now let's get to drawing...
  setH1Style();
  auto can = std::make_unique<TCanvas>("can", "", 200, 10, 1000, 1000);

  // has to be declared here to survive the entire scope - always cd to canvas first to ensure correct scaling
  can->cd();
  auto pad1 = std::make_unique<TPad>("pad1", "", 0., 0.305, 1., 1.);
  can->cd();
  auto pad2 = std::make_unique<TPad>("pad2", "", 0., 0., 1., 0.30);

  // for tacking on some text on the plot
  const std::string topMid = "#Lambda = " + toStr(eftLambda) + " TeV";
  TLatex txt;
  txt.SetTextAlign(13);

  const std::string drawOpt = "hist e2 ";

  if (noRatio) {
    can->cd();
    can->SetTopMargin(0.045);
    can->SetBottomMargin(0.11);
    can->SetLeftMargin(0.11);
    can->SetRightMargin(0.025);

    txt.SetTextSize(0.039);

    styleLegend(leg.get(), 2, 0, 0, 42, 0.037, "");
    //putLegend(leg.get(), 0.555, 0.955, 0.615, 0.875); // top right
    //putLegend(leg.get(), 0.135, 0.615, 0.615, 0.875); // top left
    //putLegend(leg.get(), 0.455, 0.955, 0.135, 0.395); // bottom right
    putLegend(leg.get(), 0.135, 0.615, 0.135, 0.395); // bottom left

    if (drawLogY) can->SetLogy();

    m_hist.at({dataName, Sample::all})->Draw("axis");
    leg->Draw();
    txt.DrawLatexNDC(0.461, 0.923, topMid.c_str());

    for (const auto &hist : m_hist) {
      if (hist.first.first == dataName) continue;

      hist.second->Draw((drawOpt + "same").c_str());
    }
    m_hist.at({dataName, Sample::all})->Draw("p e x0 same");

    can->RedrawAxis();
  }
  else {
    can->cd();
    pad1->SetTopMargin(0.05);
    pad1->SetBottomMargin(0.);
    pad1->SetLeftMargin(0.095);
    pad1->SetRightMargin(0.035);
    pad1->Draw();
    pad1->cd();

    txt.SetTextSize(0.051);

    styleLegend(leg.get(), 2, 0, 0, 42, 0.041, legHeader.c_str());
    //putLegend(leg.get(), 0.555, 0.955, 0.615, 0.875); // top right
    putLegend(leg.get(), 0.115, 0.615, 0.615, 0.875); // top left
    //putLegend(leg.get(), 0.455, 0.955, 0.055, 0.315); // bottom right

    if (drawLogY) pad1->SetLogy();

    m_hist.at({dataName, Sample::all})->Draw("axis");
    leg->Draw();
    txt.DrawLatexNDC(0.451, 0.913, topMid.c_str());

    for (const auto &hist : m_hist) {
      if (hist.first.first == dataName) continue;

      hist.second->Draw((drawOpt + "same").c_str());
    }
    m_hist.at({dataName, Sample::all})->Draw("p e x0 same");

    pad1->RedrawAxis();

    can->cd();
    pad2->SetTopMargin(0.);
    pad2->SetBottomMargin(0.21);
    pad2->SetLeftMargin(0.095);
    pad2->SetRightMargin(0.035);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    m_ratio.at({dataName, Sample::all})->Draw("axis");
    for (const auto &ratio : m_ratio) {
      if (ratio.first.first == dataName) continue;

      ratio.second->Draw((drawOpt + "same").c_str());
    }
  }

  // we're done, save the plot
  can->cd();
  can->SaveAs((plotName + ".pdf").c_str());
  //can->SaveAs((plotName + ".C").c_str());

  // save them into a file
  auto file = std::make_unique<TFile>((plotName + ".root").c_str(), "recreate");
  file->cd();
  //can->Write();
  for (auto &hist : m_hist) hist.second->Write();
  for (auto &ratio : m_ratio) ratio.second->Write();

  std::cout << std::endl;
}



void EFTFitter::autoCovMatStatCorr()
{
  if (m_covMat.count("statCorr")) {
    std::cout << "Statistical covariance matrix already available!! Aborting..." << std::endl;
    return;
  }

  std::cout << "Making automatic statistical covariance matrix assuming no correlations..." << std::endl << std::endl;

  const int nBin = m_binContent.at({dataName, Sample::all}).size();

  TMatrixD tmpMat(nBin - 2, nBin - 2);
  for (int iR = 2; iR < nBin; ++iR ) {
    for (int iC = 2; iC < nBin; ++iC )
      tmpMat(iR - 2, iC - 2) = (iC == iR) ? 1. : 0.;
  }

  m_covMat.insert({"statCorr", tmpMat});
  autoCovMatStatBin();
}



void EFTFitter::autoCovMatRate(const double &rateUnc)
{
  if (fitMode == Fit::hybrid) {
    std::cout << "Automatic matrix generation not supported in Fit::hybrid mode. Aborting..." << std::endl;
    return;
  }

  if (!hasData) {
    std::cout << "No data content available, aborting..." << std::endl << std::endl;
    return;
  }

  std::cout << "Making automatic rate uncertainty matrix from data with " << rateUnc << " factor..." << std::endl << std::endl;

  // read off the data content
  const std::vector<std::array<double, 2>> &v_binC = m_binContent.at({dataName, Sample::all});
  const double &integral = (statMode == Stat::xsec) ? v_binC.at(0).at(0) : v_binC.at(1).at(0);
  const int nBin = v_binC.size();

  // make the rate matrix with requested rate uncertainty parameter
  TMatrixD tmpMat(nBin - 2, nBin - 2);
  for (int iR = 2; iR < nBin; ++iR ) {
    for (int iC = 2; iC < nBin; ++iC ) {
      tmpMat(iR - 2, iC - 2) = (fitMode == Fit::absolute) ? 
        (rateUnc * integral) * (rateUnc * integral) * v_binC.at(iR).at(0) * v_binC.at(iC).at(0) : 0.;
    }
  }

  m_covMat.insert({"rateErr", tmpMat});
}



void EFTFitter::readCovMatRoot(const std::string &keyMat, const std::string &fileName, const std::string &histName,
                               const std::vector< std::array<int, 2> > &v_endbin)
{
  std::cout << "Reading matrix " << histName << " from file " << fileName << 
    " and assigning it as " << keyMat << std::endl << std::endl;

  const int nBin = v_rawBin.size() - 1;

  // if bin index range is invalid (like the default), then assume the hist matches raw inputs
  std::vector<int> v_index;
  if (!v_endbin.empty() and std::all_of(std::begin(v_endbin), std::end(v_endbin), 
                                        [] (const auto &endbin) {return endbin.at(1) >= endbin.at(0);})) {
    for (const auto &endbin : v_endbin) {
      v_index.push_back(endbin.at(0));
      fillInterval(v_index, endbin.at(0), endbin.at(1), 1);
    }
  }
  else {
    v_index.push_back(1);
    fillInterval(v_index, 1, nBin, 1);
  }

  if (v_index.size() != nBin) {
    std::cout << "Number of bins in fit: " << nBin << std::endl;
    std::cout << "List of requested index:" << std::endl;
    printAll(v_index);

    throw std::range_error( "Requested bin range to read covMat does not match the number of bins of available input!!!" );
  }

  // open the input file, get the hist
  auto file = std::unique_ptr<TFile>(TFile::Open( fileName.c_str() ));
  auto hist = std::unique_ptr<TH2D>(dynamic_cast<TH2D *>( (file->Get( histName.c_str() ))->Clone() ));

  TMatrixD tmpMat(nBin, nBin);
  for (int iR = 0; iR < v_index.size(); ++iR) {
    for (int iC = 0; iC < v_index.size(); ++iC)
      tmpMat(iR, iC) = hist->GetBinContent(v_index.at(iR), v_index.at(iC));
  }

  // warn if already available and erase
  if (m_covMat.count(keyMat)) {
    std::cout << "Covariance matrix with key " << keyMat << " already available, overwriting it..." << std::endl;
    m_covMat.erase(keyMat);
  }

  m_covMat.insert({keyMat, tmpMat});
  if (keyMat == "statCorr")
    autoCovMatStatBin();
}



void EFTFitter::makeFinalCovMat(const std::function<TMatrixD (const std::map<std::string, TMatrixD> &)> &func)
{
  // warn if already available and erase
  if (m_covMat.count("finalcov")) {
    std::cout << "Final covariance matrix already available, overwriting it..." << std::endl;
    m_covMat.erase("finalcov");
  }

  std::cout << "Making the final covariance matrix with the given function..." << std::endl << std::endl;

  // save the function in case toy fit method is gonna be called
  if (f_finalCov)
    std::cout << "Final covariance function already available, overwriting it..." << std::endl;
  f_finalCov = func;

  // just apply the function passed by the user and tack the result onto the map
  m_covMat.insert({"finalcov", f_finalCov(m_covMat)});
}



void EFTFitter::makeFinalCovMat(std::vector<std::string> &&names)
{
  if (m_covMat.empty()) {
    std::cout << "This method do nothing if no matrices have been read in/created..." << std::endl;
    return;
  }

  if (names.empty())
    std::cout << "Empty list of covariance matrix names given, assuming that final cov is the sum of all existing matrices..." << std::endl;

  // warn if already available and erase
  if (m_covMat.count("finalcov")) {
    std::cout << "Final covariance matrix already available, overwriting it..." << std::endl;
    m_covMat.erase("finalcov");
  }

  TMatrixD finalmat(m_covMat.at(names[0]));
  for (int iname = 1; iname < names.size(); ++iname)
    finalmat += m_covMat.at(names[iname]);
  m_covMat.insert({"finalcov", finalmat});
}



void EFTFitter::drawCovMat(const std::string &dirName, const std::vector<std::string> &v_keyMat, bool do_correlation) const
{
  // quick checks to ensure we don't make unnecessary files
  if (m_covMat.empty()) {
    std::cout << "No point in calling this method when the covMat map is empty, aborting..." << std::endl;
    return;
  }

  std::cout << "Drawing the requested covariance matrices..." << std::endl;
  if (!v_keyMat.empty())
    printAll(v_keyMat);
  else
    printAll(extractKey(m_covMat));

  int matrixToSave = (v_keyMat.empty()) ? 1 : 0;
  for (const auto &keyMat : v_keyMat)
    matrixToSave += m_covMat.count(keyMat);
  if (!matrixToSave) {
    std::cout << "Requested key list results in no matrix to be saved, aborting..." << std::endl;
    return;
  }

  // create a file in which to save the matrices
  auto file = std::make_unique<TFile>((dirName + "cov_all.root").c_str(), "recreate");

  const int nBin = std::begin(m_covMat)->second.GetNrows();

  for (const auto &p_covMat : m_covMat) {
    if (!v_keyMat.empty() and !std::count(std::begin(v_keyMat), std::end(v_keyMat), p_covMat.first)) continue;

    auto h_cMat = std::make_unique<TH2D>(("cov_" + p_covMat.first).c_str(), "", nBin, 0., nBin, nBin, 0., nBin);
    std::ofstream t_cMat(dirName + "cov_" + p_covMat.first + ".txt");

    for (int iR = 0; iR < nBin; ++iR) {
      for (int iC = 0; iC < nBin; ++iC) {
        const auto value = (do_correlation) ? p_covMat.second(iR, iC) / std::sqrt(p_covMat.second(iR, iR) * p_covMat.second(iC, iC)) : p_covMat.second(iR, iC);

        h_cMat->SetBinContent(iR + 1, iC + 1, value);

        t_cMat << std::setprecision(7) << value;
        if (iC != nBin - 1)
          t_cMat << "    ";
      }

      t_cMat << std::endl;
    }

    t_cMat.close();

    stylePlot(h_cMat.get(), kBlack, 1., 0, 20, 1.5, 1, 1); // 3rd last arg is marker size, affects text size too...
    axisPlot(h_cMat.get(), 0., 0., 505, "Column", 0.043, 1.11, 0.037, 0., 0., 505, "Row", 0.043, 1.11, 0.037);

    // apply style - reminder: palette with dark colors in one end likely looks bad for matrices with elements of both sign
    setH2Style();
    gStyle->SetPalette( kAquamarine );
    //TColor::InvertPalette();

    auto can = std::make_unique<TCanvas>(("can_" + p_covMat.first).c_str(), "", 200, 10, 1000, 1000);

    // FIXME incomplete: colz in fact better but needs canvas rearranging
    const std::string drawOpt = (nBin > 10) ? "col" : "coltext";

    can->SetTopMargin(0.045);
    can->SetBottomMargin(0.11);
    can->SetLeftMargin(0.11);
    can->SetRightMargin(0.025);

    // for tacking on some text on the plot
    TLatex txt;
    txt.SetTextSize(0.043);
    txt.SetTextAlign(13);

    can->cd();

    h_cMat->Draw(drawOpt.c_str());

    can->SaveAs((dirName + "cov_" + p_covMat.first + ".pdf").c_str());
    //can->SaveAs((dirName + "cov_" + p_covMat.first + ".C").c_str());

    file->cd();
    //can->Write();
    h_cMat->Write();
  }

  std::cout << std::endl;
}



std::array<double, 2> EFTFitter::getDataRate() const
{
  if (!hasData or fitMode != Fit::absolute or m_covMat.count("finalcov") == 0) {
    std::cout << "Method is useless if not doing absolute fit or "
      "there is no data template or final covariance matrix!! Returning junk!!" << std::endl;
    return {0., 0.};
  }

  const int nBin = v_rawBin.size() - 1, iXs = (statMode == Stat::xsec) ? 0 : 1;

  // prepare the vectors to get the single rate covariance
  TMatrixD rowM(1, nBin);
  for (int iB = 0; iB < nBin; ++iB)
    rowM(0, iB) = 1.;

  // and we do exactly that
  TMatrixD rateCov(rowM, TMatrixD::kMult, TMatrixD(m_covMat.at("finalcov"), TMatrixD::kMultTranspose, rowM));
  //std::cout << rateCov.GetNrows() << " " << rateCov.GetNcols() << std::endl;
  return {m_binContent.at({dataName, Sample::all}).at(iXs).at(0), std::sqrt(rateCov(0, 0))};
}



void EFTFitter::prepareInterpolationBase()
{
  std::cout << "Preparing the base contents for all operators..." << std::endl;

  // copy the bin content map into another - keeping only the ones relevant for interpolation
  std::map< std::vector<double>, std::vector<std::array<double, 2>> > m_interInput;
  uint nOpAll0 = 0, nOp1Non0 = 0, nOp2Non0 = 0;

  for (const auto &p_binContent : m_binContent) {
    const std::string &key = p_binContent.first.first;
    const Sample &samp = p_binContent.first.second;
    const auto &content = p_binContent.second;

    if (key == dataName) continue;

    // interpolation takes only all samples as input
    if (samp != Sample::all) continue;

    // get the operator coeff values
    const auto v_opVal = extractValue( parseOpFromKey(key) );

    // is there a SM ie all operators are 0 case?
    if (std::all_of(std::begin(v_opVal), std::end(v_opVal), [] (const auto &opV) {return opV == 0.;})) {
      m_interInput.insert({v_opVal, content});
      ++nOpAll0;

      // while we're at it we copy the SM case into the content map for linear
      // needed for drawChi2 methods
      m_binContent.insert({{key, Sample::linear}, content});

      // no point going further in the loop at this point
      continue;
    }

    // is there a case where exactly 1 or 2 operator != 0?
    const int nOpNon0 = std::count_if(std::begin(v_opVal), std::end(v_opVal), [] (const auto &opV) {return opV != 0.;});

    if (nOpNon0 == 1 or nOpNon0 == 2) {
      m_interInput.insert({v_opVal, content});

      if (nOpNon0 == 1) 
        ++nOp1Non0;

      if (nOpNon0 == 2) 
        ++nOp2Non0;
    }
  }

  // ensure the minimal input for the interpolation is present
  const int nOp = v_opName.size();
  if (nOpAll0 == 0 or nOp1Non0 < 2 * nOp or nOp2Non0 < nOp - 1) {
    printAll(v_opName);
    throw std::logic_error( "Minimal required input for interpolating the requested operator set is not available!!" );
  }

  // these are gonna pop up often
  // std::pow not as efficient for integer powers, so yeah
  const int nBin = (std::begin(m_binContent)->first.first != dataName) ? std::begin(m_binContent)->second.size() 
    : std::next(std::begin(m_binContent))->second.size();
  const int iXs = (statMode == Stat::xsec) ? 0 : 1;
  const double sqLambda = eftLambda * eftLambda; 
  const double quLambda = sqLambda * sqLambda;

  // ok let's start interpolating - first prepare the bin contents starting from SM
  const auto v_xs_op0 = std::find_if(std::begin(m_interInput), std::end(m_interInput), [] (const auto &p) {
        return std::all_of(std::begin(p.first), std::end(p.first), [] (const auto &op) {return op == 0.;});
      });

  // and we fill this map up - it'll be the base input for our interpolation
  for (int iOp = 0; iOp < nOp; ++iOp) {
    const std::string &opName = v_opName.at(iOp);

    // start by preparing the bin contents of op being considered - first 2 non-0 instance
    const auto f_isOp1Eq1 = [&] (const auto &p) {
      return (std::count_if(std::begin(p.first), std::end(p.first), [] (const auto &c) {return c != 0.;}) == 1 and 
              p.first.at(iOp) != 0.);
    };
    const auto v_xs_op1 = std::find_if(std::begin(m_interInput), std::end(m_interInput), f_isOp1Eq1);
    const auto v_xs_op2 = std::find_if(std::next(v_xs_op1), std::end(m_interInput), f_isOp1Eq1);

    // a quick check ensuring these are correct
    if (v_xs_op1 == std::end(m_interInput) or v_xs_op2 == std::end(m_interInput))
      throw std::logic_error( "There are insufficient input to interpolate operator " + opName + "!!!");

    // get the value of the 1st and 2nd coeff of op
    const double &c1_op = v_xs_op1->first.at(iOp), &c2_op = v_xs_op2->first.at(iOp);

    // remember: binContent vector has the total count at 1 and normalized to 1 shape for the rest
    std::vector<std::array<double, 2>> v_opI(nBin, {0., 0.}), v_opR(nBin, {0., 0.});
    for (int iB = 2; iB < nBin; ++iB) {
      // SM xsec, error in bin iB
      const double xs_op0 = v_xs_op0->second.at(iXs).at(0) * v_xs_op0->second.at(iB).at(0);
      const double xe_op0 = v_xs_op0->second.at(iXs).at(0) * v_xs_op0->second.at(iB).at(1);

      // op xsec, error at c1 and c2
      const double xs_op1 = v_xs_op1->second.at(iXs).at(0) * v_xs_op1->second.at(iB).at(0);
      const double xe_op1 = v_xs_op1->second.at(iXs).at(0) * v_xs_op1->second.at(iB).at(1);

      const double xs_op2 = v_xs_op2->second.at(iXs).at(0) * v_xs_op2->second.at(iB).at(0);
      const double xe_op2 = v_xs_op2->second.at(iXs).at(0) * v_xs_op2->second.at(iB).at(1);

      // compute the terms that go into eq 3.7 (correcting for the presumed typo in the linear term)
      const double c1_x_c2 = c1_op * c2_op;
      const double c1_p_c2 = c1_op + c2_op;
      const double c2_m_c1 = c2_op - c1_op;
      const double cc1 = c1_op * c1_op, cc2 = c2_op * c2_op;

      // the coeffs are exactly known, so compute errors only for the xsec terms
      const double cc1x2 = cc1 * xs_op2, cc2x1 = cc2 * xs_op1;
      const double cc1e2 = cc1 * xe_op2, cc2e1 = cc2 * xe_op1;
      const double cc2x1_m_cc1x2 = cc2x1 - cc1x2;
      const double cc2e1_m_cc1e2 = std::sqrt((cc2e1 * cc2e1) + (cc1e2 * cc1e2));

      const double c1x2 = c1_op * xs_op2, c2x1 = c2_op * xs_op1;
      const double c1e2 = c1_op * xe_op2, c2e1 = c2_op * xe_op1;

      const double c2x1_m_c1x2 = c2x1 - c1x2;
      const double c2e1_m_c1e2 = std::sqrt((c2e1 * c2e1) + (c1e2 * c1e2));

      // then the xsec, error for this op at coeff = 1
      const double sqLambda_o_c1c2 = sqLambda / c1_x_c2;
      const double xs_opI = sqLambda_o_c1c2 * ((cc2x1_m_cc1x2 / c2_m_c1) - (c1_p_c2 * xs_op0));
      const double xe_opI = std::abs(sqLambda_o_c1c2) * std::sqrt(((cc2e1_m_cc1e2 / c2_m_c1) * (cc2e1_m_cc1e2 / c2_m_c1)) + 
                                                                  ((c1_p_c2 * xe_op0) * (c1_p_c2 * xe_op0)));

      const double quLambda_o_c1c2 = quLambda / c1_x_c2;
      const double xs_opR = quLambda_o_c1c2 * (xs_op0 - (c2x1_m_c1x2 / c2_m_c1));
      const double xe_opR = std::abs(quLambda_o_c1c2) * std::sqrt((xe_op0 * xe_op0) + 
                                                                  ((c2e1_m_c1e2 / c2_m_c1) * (c2e1_m_c1e2 / c2_m_c1)));

      // store the xsec of this bin in vector
      v_opI.at(iB).at(0) = xs_opI;
      v_opI.at(iB).at(1) = xe_opI;

      v_opR.at(iB).at(0) = xs_opR;
      v_opR.at(iB).at(1) = xe_opR;
    }

    // here we use the sample tags to refer to the individual op parts
    m_op1Eq1.insert({{opName, Sample::linear}, v_opI});
    m_op1Eq1.insert({{opName, Sample::quadratic}, v_opR});
  }

  // for 2D obtain the op1 - op2 interference term (treated as res since it's quadratic)
  for (int iOp1 = 0; iOp1 < nOp; ++iOp1) {
    for (int iOp2 = iOp1 + 1; iOp2 < nOp; ++iOp2) {
      const std::array<std::string, 2> a_opName = {v_opName.at(iOp1), v_opName.at(iOp2)};

      // get the bin contents where the 2 ops are the only non-0 ones
      const auto v_xs_op12 = std::find_if(std::begin(m_interInput), std::end(m_interInput), [&] (const auto &p) {
          return (std::count_if(std::begin(p.first), std::end(p.first), [] (const auto &c) {return c != 0.;}) == 2 and
                  p.first.at(iOp1) != 0. and p.first.at(iOp2) != 0.);
        });

      // and the coeff values themselves
      const double &c_op1 = v_xs_op12->first.at(iOp1), &c_op2 = v_xs_op12->first.at(iOp2);

      // and it's the usual business again
      std::vector<std::array<double, 2>> v_op12(nBin, {0., 0.});
      for (int iB = 2; iB < nBin; ++iB) {
        // SM xsec, error in bin iB
        const double xs_op0 = v_xs_op0->second.at(iXs).at(0) * v_xs_op0->second.at(iB).at(0);
        const double xe_op0 = v_xs_op0->second.at(iXs).at(0) * v_xs_op0->second.at(iB).at(1);

        // total op12 xsec, error in bin iB
        const double xs_op12 = v_xs_op12->second.at(iXs).at(0) * v_xs_op12->second.at(iB).at(0);
        const double xe_op12 = v_xs_op12->second.at(iXs).at(0) * v_xs_op12->second.at(iB).at(1);

        // op1 and op2 xsec, error for each part
        const double op1ConI = m_op1Eq1.at({a_opName.at(0), Sample::linear}).at(iB).at(0);
        const double op1ErrI = m_op1Eq1.at({a_opName.at(0), Sample::linear}).at(iB).at(1);
        const double op1ConR = m_op1Eq1.at({a_opName.at(0), Sample::quadratic}).at(iB).at(0);
        const double op1ErrR = m_op1Eq1.at({a_opName.at(0), Sample::quadratic}).at(iB).at(1);
        const double cc1 = c_op1 * c_op1;

        const double xs_op1I = (c_op1 / sqLambda) * op1ConI;
        const double xe_op1I = std::abs(c_op1 / sqLambda) * op1ErrI;
        const double xs_op1R = (cc1 / quLambda) * op1ConR;
        const double xe_op1R = std::abs(cc1 / quLambda) * op1ErrR;

        const double op2ConI = m_op1Eq1.at({a_opName.at(1), Sample::linear}).at(iB).at(0);
        const double op2ErrI = m_op1Eq1.at({a_opName.at(1), Sample::linear}).at(iB).at(1);
        const double op2ConR = m_op1Eq1.at({a_opName.at(1), Sample::quadratic}).at(iB).at(0);
        const double op2ErrR = m_op1Eq1.at({a_opName.at(1), Sample::quadratic}).at(iB).at(1);
        const double cc2 = c_op2 * c_op2;

        const double xs_op2I = (c_op2 / sqLambda) * op2ConI;
        const double xe_op2I = std::abs(c_op2 / sqLambda) * op2ErrI;
        const double xs_op2R = (cc2 / quLambda) * op2ConR;
        const double xe_op2R = std::abs(cc2 / quLambda) * op2ErrR;

        // op1 and op2 int in the bin is just the remainder of subtracting the individual parts
        v_op12.at(iB).at(0) = (quLambda / (c_op1 * c_op2)) * (xs_op12 - xs_op0 - xs_op1I - xs_op1R - xs_op2I - xs_op2R);
        v_op12.at(iB).at(1) = std::abs(quLambda / (c_op1 * c_op2)) * std::sqrt((xe_op12 * xe_op12) + (xe_op0 * xe_op0) + 
                                                                               (xe_op1I * xe_op1I) + (xe_op1R * xe_op1R) + 
                                                                               (xe_op2I * xe_op2I) + (xe_op2R * xe_op2R));
      }

      // and into map it goes
      m_op2Eq1.insert({a_opName, v_op12});
    }
  }

  // to find out the contribution of each case
  std::cout << "Contribution of each operator at 1 (in " << statMode << ") is:" << std::endl;
  for (auto &p_op1Eq1 : m_op1Eq1) {
    const double count = std::accumulate(std::begin(p_op1Eq1.second), std::end(p_op1Eq1.second), 0., 
                                         [] (const double &sum, const auto &con) { return sum + con.at(0); });
    const double sqerr = std::accumulate(std::begin(p_op1Eq1.second), std::end(p_op1Eq1.second), 0., 
                                         [] (const double &err, const auto &con) { return err + (con.at(1) * con.at(1)); });

    std::cout << std::make_pair(p_op1Eq1.first, std::array<double, 2>{count, std::sqrt(sqerr)}) << std::endl;
  }
  for (auto &p_op2Eq1 : m_op2Eq1) {
    const double count = std::accumulate(std::begin(p_op2Eq1.second), std::end(p_op2Eq1.second), 0., 
                                         [] (const double &sum, const auto &con) { return sum + con.at(0); });
    const double sqerr = std::accumulate(std::begin(p_op2Eq1.second), std::end(p_op2Eq1.second), 0., 
                                         [] (const double &err, const auto &con) { return err + (con.at(1) * con.at(1)); });

    std::cout << std::make_pair(p_op2Eq1.first, std::array<double, 2>{count, std::sqrt(sqerr)}) << std::endl;
  }

  std::cout << "For reference: SM " << statMode << " is " << v_xs_op0->second.at(iXs) << std::endl << std::endl;
}



void EFTFitter::assignAsData(const std::string &keyName, const Sample sampleType, const bool varyWithinError)
{
  std::cout << "Assigning key " << keyName << " of type " << sampleType << " as " << dataName << "..." << std::endl;
  std::cout << "Note: previous data is overwritten by this method!!" << std::endl << std::endl;
  hasData = true;

  if (varyWithinError)
    std::cout << "Templates will be separately varied within uncertainties for (SM + linear) and quadratic parts..." 
              << std::endl << std::endl;

  // grab the content to be assigned, interpolate if necessary
  const std::pair<std::string, Sample> keyPair = {fixKeyFormat(keyName), sampleType}, dataKey = {dataName, Sample::all};
  const std::vector<std::array<double, 2>> v_binC = interpolateOpValue(keyName, sampleType, varyWithinError);

  if (m_binContent.count(dataKey))
    m_binContent.erase(dataKey);

  m_binContent.insert({dataKey, v_binC});
}



void EFTFitter::listKeyToFit(const std::map<std::string, std::vector<double>> &m_opGrid)
{
  std::cout << "Unpacking the grid of requested operator values into individual keys..." << std::endl << std::endl;

  // check consistency of operator set and extract their names
  if (extractKey(m_opGrid) != v_opName)
    throw std::logic_error( "Requested operator set for interpolation is not consistent with available input!!" );

  // may be large, so reserve the space beforehand
  const ulong nKey = std::accumulate(std::begin(m_opGrid), std::end(m_opGrid), 1, 
                                     [] (ulong nKey, const auto &p) {return nKey * p.second.size();});
  v_keyToFit.reserve( nKey );
  unpackOpGrid(v_keyToFit, m_opGrid, std::begin(m_opGrid), "");
  v_keyToFit.shrink_to_fit();
  std::sort(std::begin(v_keyToFit), std::end(v_keyToFit));
  v_keyToFit.erase(std::unique(std::begin(v_keyToFit), std::end(v_keyToFit)), std::end(v_keyToFit));
  //printAll(v_keyToFit);
}



void EFTFitter::computeFitChi2(const std::vector<Sample> &v_sample, int binToIgnore)
{
  const bool rateFit = (fitMode == Fit::hybrid and v_rawBin.empty());

  if (!hasData or v_keyToFit.empty() or (!m_covMat.count("finalcov") and !rateFit))
    throw std::range_error( "This method shouldn't be called given the insuffiencient input! "
                            "Ensure contents (inserted with addRawInput() or interpolatable with prepareInterpolationBase()), "
                            "list of keys to be fitted on (made with listKeyToFit()) "
                            "and final covmat (made with the auto/readCovMat methods and finally makeFinalCovMat()) are available!" );

  std::cout << "Computing the chi2 for all requested keys..." << std::endl << std::endl;

  const std::vector<std::array<double, 2>> &dataContent = m_binContent.at({dataName, Sample::all});
  const int nBin = dataContent.size(), iXs = (statMode == Stat::xsec) ? 0 : 1;
  const double dataInt = getContentSum(dataContent);

  // ok here we copy and invert the matrix because this is what we actually use
  // a bit of shenanigan needed in case of shape fit since we need to ignore one bin for convertible matrix
  // on the other hand interpolation needs to see the full template, so it must not be done any earlier
  // bin index to ignore means every nth bin in case of stitched templates (assumed to have same nBinEach)
  // number of stitched templates obtained assuming shapes sum up to 1
  const int nBinEach = (nBin - 2) / int(shapeSum);
  const int nHist = (nBin - 2) / nBinEach;
  //std::cout << "ignore " << binToIgnore << " nBin " << nBin << " nBinEach " << nBinEach << " nHist " << nHist << std::endl;

  // reset the bin index in case requested argument doesn't make sense
  // probably fine to do so silently, it's only used to verify the bin dropper works correctly
  if (binToIgnore >= nBinEach or binToIgnore < 0)
    binToIgnore = 1;

  TMatrixD invMat;
  // do nothing for rate fits - no matrix needed as there's no template
  if (rateFit)
    ;
  else if (fitMode != Fit::shape) {
    // must resize before copy-assign: https://root.cern.ch/how/how-create-and-fill-matrix
    invMat.ResizeTo(m_covMat.at("finalcov"));
    invMat = m_covMat.at("finalcov");
  }
  else {
    invMat.ResizeTo(nBin - 2 - nHist, nBin - 2 - nHist);

    // two indices needed as these run off differently :(
    int iShpR = 0, iShpC = 0;

    for (int iAbsR = 0; iAbsR < nBin - 2; ++iAbsR) {
      if (iAbsR % nBinEach == binToIgnore) continue;

      iShpC = 0;
      for (int iAbsC = 0; iAbsC < nBin - 2; ++iAbsC) {
        if (iAbsC % nBinEach == binToIgnore) continue;

        invMat(iShpR, iShpC) = m_covMat.at("finalcov")(iAbsR, iAbsC);
        ++iShpC;
      }

      ++iShpR;
    }
  }

  // std::cout << "" << TDecompLU(invMat).Condition() << std::endl; // condition number of the matrix
  invMat.Invert();

  for (const auto &key : v_keyToFit) {
    for (const auto &samp : v_sample) {
      const std::vector<std::array<double, 2>> opContent = interpolateOpValue(key, samp);
      const double opInt = getContentSum(opContent);

      // here they're used as offset counters; increment every time a bin is ignored
      int iShpR = 0, iShpC = 0;
      double fitChi2 = 0.;
      for (int iAbsR = 2; iAbsR < nBin; ++iAbsR) {
        if (fitMode == Fit::shape and ((iAbsR - 2) % nBinEach == binToIgnore)) {
          ++iShpR;
          continue;
        }

        // the actual matrix index, which is strictly related to the actual row/col indices and offset
        // but doing it like this improves readability slightly
        int iMatR = iAbsR - 2 - iShpR;

        iShpC = 0;
        for (int iAbsC = iAbsR; iAbsC < nBin; ++iAbsC) {
          if (fitMode == Fit::shape and ((iAbsC - 2) % nBinEach == binToIgnore)) {
            ++iShpC;
            continue;
          }
          int iMatC = iAbsC - 2 - iShpC;
          double factor = (iMatR == iMatC) ? 1. : 2.;

          // first compute the bin differences
          const double deltaR = (dataInt * dataContent.at(iAbsR).at(0)) - (opInt * opContent.at(iAbsR).at(0));
          const double deltaC = (dataInt * dataContent.at(iAbsC).at(0)) - (opInt * opContent.at(iAbsC).at(0));

          fitChi2 += deltaR * deltaC * invMat(iMatR, iMatC) * factor;
        }
      }

      if (fitMode == Fit::hybrid and dataContent.at(iXs).at(0) > 0.) {
        const double deltaS = dataContent.at(iXs).at(0) - opContent.at(iXs).at(0);
        fitChi2 += (deltaS * deltaS) / (dataContent.at(iXs).at(1) * dataContent.at(iXs).at(1));
      }

      m_fitChi2.insert({{key, samp}, fitChi2});
    }
  }

  //printAll(m_fitChi2);
}



void EFTFitter::saveFitChi2(const std::string &fileName)
{
  if (m_fitChi2.empty())
    throw std::logic_error( "This method shouldn't be called before the computeFitChi2() method!!" );

  auto file = std::make_unique<TFile>((fileName + "_chi2.root").c_str(), "recreate");
  auto tree = std::make_unique<TTree>("fit_chi2", "", 505);
  tree->SetAutoSave(0);
  tree->SetImplicitMT(false);

  std::vector<double> v_opVal(v_opName.size());
  double chi2 = 0.;
  int sample = -1;
  for (uint iOp = 0; iOp < v_opName.size(); ++iOp)
    tree->Branch(v_opName[iOp].c_str(), &v_opVal[iOp], (v_opName[iOp] + "/D").c_str());
  tree->Branch("chi2", &chi2, "chi2/D");
  tree->Branch("sample", &sample, "sample/I");

  for (const auto &p_fitChi2 : m_fitChi2) {
    const std::vector<double> v_opTmp = extractValue( parseOpFromKey(p_fitChi2.first.first) );
    sample = static_cast<int>(p_fitChi2.first.second);

    chi2 = p_fitChi2.second;
    for (uint iOp = 0; iOp < v_opName.size(); ++iOp)
      v_opVal[iOp] = v_opTmp[iOp];

    tree->Fill();
  }

  file->cd();
  tree->Write();
}



void EFTFitter::draw1DChi2(const std::map<std::string, std::tuple<std::string, std::vector<double>, 
                           std::array<double, 2>, std::array<double, 2>>> &mt_opInfo, 
                           const std::string &plotName, const std::vector<Sample> &v_sample) const
{
  if (m_fitChi2.empty())
    throw std::logic_error( "This method shouldn't be called before the computeFitChi2() method!!" );

  std::cout << "Drawing the 1D dChi2 graph for each operator..." << std::endl << std::endl;

  // bind the operator-chi2 set in a map; consider only those with at most 1 op != 0
  std::map< std::pair<std::vector<double>, Sample>, double> m_opMax1Non0;
  for (const auto &p_fitChi2 : m_fitChi2) {
    const std::vector<double> v_opVal = extractValue( parseOpFromKey(p_fitChi2.first.first) );

    // is there a SM ie all operators are 0 case?
    if (std::all_of(std::begin(v_opVal), std::end(v_opVal), [] (const double &opV) {return opV == 0.;}))
      m_opMax1Non0.insert( {{v_opVal, p_fitChi2.first.second}, p_fitChi2.second} );

    // is there a case where exactly 1 operator != 0?
    if (std::count_if(std::begin(v_opVal), std::end(v_opVal), [] (const double &opV) {return opV != 0.;}) == 1)
      m_opMax1Non0.insert( {{v_opVal, p_fitChi2.first.second}, p_fitChi2.second} );
  }

  // declare here things that are gonna be applicable for all anyway
  const std::array<std::string, 2> a_sampLeg = {" (all)", " (linear)"};
  const std::array<int, 2> a_sampStyL = {1, 2};
  const std::array<double, 2> a_zero = {0., 0.};
  const std::array<double, 2> a_vPnt = {0., 8.};
  const std::array<double, 2> a_vErrD = {0., 4.};
  const std::array<double, 2> a_vErrU = {4., 0.};

  // nBin - (1 for shape bin drop) + (1 if hybrid has rate) - nOp
  const int iXs = (statMode == Stat::xsec) ? 0 : 1;
  const int nDoF = v_rawBin.size() - 1 + v_rawBin.empty()
    - ((fitMode == Fit::shape) * int(shapeSum))
    + (fitMode == Fit::hybrid and m_binContent.at({dataName, Sample::all}).at(iXs).at(0) > 0.)
    - 1;

  // 1 and 2 sigma band finder - always such that it's narrowest around minimum
  const auto f_dChi2Band = [] (std::vector<double> &v_dChi2, const double &dChi2Cut) {
    // first U side
    // define the minimum iterator
    auto fMin = std::min_element(std::begin(v_dChi2), std::end(v_dChi2));
    auto itU = std::find_if(fMin, std::end(v_dChi2), 
                            [&] (const double &dChi2) {return dChi2 > dChi2Cut;});

    // and then the D side - same thing but with reverse iterators
    auto rMin = std::min_element(std::rbegin(v_dChi2), std::rend(v_dChi2));
    auto itD = std::find_if(rMin, std::rend(v_dChi2), 
                            [&] (const double &dChi2) {return dChi2 > dChi2Cut;});

    // if any of them aren't found, return forward end
    std::pair<decltype(fMin), decltype(fMin)> p_band = {std::end(v_dChi2), std::end(v_dChi2)};
    if (itU == std::end(v_dChi2) or itD == std::rend(v_dChi2))
      return p_band;

    // return both as forward iterators, D then U
    p_band = {itD.base() - 1, itU};
    return p_band;
  };

  for (const auto &p_opInfo : mt_opInfo) {
    const auto &opName = p_opInfo.first;
    const int iOp = std::distance(std::begin(v_opName), std::find(std::begin(v_opName), std::end(v_opName), opName));
    const auto &opLeg = std::get<0>(p_opInfo.second);
    const auto &opRange = std::get<1>(p_opInfo.second);
    const auto &yRange = std::get<2>(p_opInfo.second);
    const auto &xRange = std::get<3>(p_opInfo.second);

    // now we make the points to be plotted - separately for all and linear
    std::array<std::vector<double>, 2> av_opVal, av_opChi2, av_dChi2;
    for (const auto &p_op1 : m_opMax1Non0) {
      const auto &v_opVal = p_op1.first.first;
      const auto &opChi2 = p_op1.second;
      const auto &samp = p_op1.first.second;
      const int iSamp = static_cast<int>(samp);

      // skip if it's not the op under consideration
      if (std::count_if(std::begin(v_opVal), std::end(v_opVal), 
                        [] (const double &opV) {return opV != 0.;}) == 1 and v_opVal.at(iOp) == 0.)
        continue;

      // if the range for the corresponding op is empty, select based on chi2 (3 sigma)
      if (opRange.empty() or
          (!opRange.empty() and v_opVal.at(iOp) >= opRange.at(0) and v_opVal.at(iOp) <= opRange.at(1))) {
        av_opVal.at(iSamp).push_back(v_opVal.at(iOp));
        av_opChi2.at(iSamp).push_back(opChi2);
      }
    }

    // and we make the graph
    bool drawSig1 = true, drawSig2 = true;
    std::array<std::unique_ptr<TGraph>, 2> ag_dChi2 = {nullptr, nullptr};
    std::array<std::array<std::unique_ptr<TGraphAsymmErrors>, 3>, 2> ag_sigma = {{{nullptr, nullptr, nullptr},
                                                                                  {nullptr, nullptr, nullptr}}};
    auto leg = std::make_unique<TLegend>();
    for (auto &samp : v_sample) {
      const int iSamp = static_cast<int>(samp);
      if (av_opVal.at(iSamp).empty()) continue;

      // grab the index of the min chi2 and make the dChi2 vector
      // note: due to the realloc below, DO NOT take references for opMin and chi2Min!!!
      const int iMin = std::distance(std::begin(av_opChi2.at(iSamp)), 
                                     std::min_element(std::begin(av_opChi2.at(iSamp)), std::end(av_opChi2.at(iSamp))));
      const double opMin = av_opVal.at(iSamp).at(iMin), chi2Min = av_opChi2.at(iSamp).at(iMin);
      for (double &chi2 : av_opChi2.at(iSamp)) 
        av_dChi2.at(iSamp).push_back(chi2 - chi2Min);

      const double chi2Prob = TMath::Prob(chi2Min, nDoF);
      std::cout << "Minimum for op "<< opName << ", sample " << samp << " found at " << opMin << 
        " with chi2/nDoF " << chi2Min << "/" << nDoF << ", with p-value " << chi2Prob << std::endl;

      // a bit of interlude in case no range is specified and so we do the dChi2-based range
      // note: due to the realloc here, DO NOT take references for opMin and chi2Min!!!
      if (opRange.empty()) {
        // temp copy for storing the filtered vectors
        std::vector<double> tmp_opVal, tmp_opChi2, tmp_dChi2;
        for (int iD = 0; iD < av_dChi2.at(iSamp).size(); ++iD) {
          if (av_dChi2.at(iSamp).at(iD) <= 9.) {
            tmp_opVal.push_back(av_opVal.at(iSamp).at(iD));
            tmp_opChi2.push_back(av_opChi2.at(iSamp).at(iD));
            tmp_dChi2.push_back(av_dChi2.at(iSamp).at(iD));
          }
        }

        // and now we replace the main ones with these temps
        av_opVal.at(iSamp).clear();
        av_opChi2.at(iSamp).clear();
        av_dChi2.at(iSamp).clear();

        av_opVal.at(iSamp) = tmp_opVal;
        av_opChi2.at(iSamp) = tmp_opChi2;
        av_dChi2.at(iSamp) = tmp_dChi2;
      }

      // and make the graph
      ag_dChi2.at(iSamp) = std::make_unique<TGraph>(av_opVal.at(iSamp).size(), av_opVal.at(iSamp).data(), av_dChi2.at(iSamp).data());
      stylePlot(ag_dChi2.at(iSamp).get(), kBlack, 1., 0, 20, 1.5, a_sampStyL.at(iSamp), 2);
      axisPlot(ag_dChi2.at(iSamp).get(), 
               yRange.at(0), yRange.at(1), 505, "#Delta #chi^{2}", 0.043, 1.21, 0.037, 
               xRange.at(0), xRange.at(1), 505, opLeg.c_str(), 0.043, 1.11, 0.037);
      ag_dChi2.at(iSamp)->SetName((opName + "_" + toStr(samp)).c_str());
      ag_dChi2.at(iSamp)->SetTitle("");

      leg->AddEntry(ag_dChi2.at(iSamp).get(), (opLeg + a_sampLeg.at(iSamp)).c_str(), "l");

      // also prepare the stuff for 2 and 1 sigma bands
      // threshold is dChi2 containing 1 and 2 sigma for chi2 with 1 ndf (since we're fitting 1 op)
      // obtained with 1. - TMath::Prob(raw_chi2, ndf)
      const auto itBegin = std::begin(av_dChi2.at(iSamp));

      // 1 sigma D, U band
      const std::pair<decltype(itBegin), decltype(itBegin)> p_it1Sig = f_dChi2Band(av_dChi2.at(iSamp), 1.);

      // check that both iterators are good
      if ( p_it1Sig.first == std::end(av_dChi2.at(iSamp)) or p_it1Sig.second == std::end(av_dChi2.at(iSamp)) ) {
        std::cout << "Requested fit range for operator " << opName << 
          " is insufficient to define the 1 sigma bands! Skipping..." << std::endl;
        drawSig1 = false;
      }

      // 2 sigma D, U band
      const std::pair<decltype(itBegin), decltype(itBegin)> p_it2Sig = f_dChi2Band(av_dChi2.at(iSamp), 4.);

      // check that both iterators are good
      if ( p_it2Sig.first == std::end(av_dChi2.at(iSamp)) or p_it2Sig.second == std::end(av_dChi2.at(iSamp)) ) {
        std::cout << "Requested fit range for operator " << opName << 
          " is insufficient to define the 2 sigma bands! Skipping..." << std::endl;
        drawSig2 = false;
      }

      // ok we've got the iterators, now get the actual errors - 2D, 2U, 1D, 1U
      const int i1SigD = std::distance(itBegin, p_it1Sig.first);
      const int i1SigU = std::distance(itBegin, p_it1Sig.second);
      const int i2SigD = std::distance(itBegin, p_it2Sig.first);
      const int i2SigU = std::distance(itBegin, p_it2Sig.second);

      const double err1SigD = (drawSig1) ? std::abs(opMin - av_opVal.at(iSamp).at(i1SigD)) : 0.;
      const double err1SigU = (drawSig1) ? std::abs(av_opVal.at(iSamp).at(i1SigU) - opMin) : 0.;
      const double err2SigD = (drawSig2) ? std::abs(opMin - av_opVal.at(iSamp).at(i2SigD)) : 0.;
      const double err2SigU = (drawSig2) ? std::abs(av_opVal.at(iSamp).at(i2SigU) - opMin) : 0.;

      const double dChi1SigD = (drawSig1) ? av_dChi2.at(iSamp).at(i1SigD) : 0.;
      const double dChi1SigU = (drawSig1) ? av_dChi2.at(iSamp).at(i1SigU) : 0.;
      const double dChi2SigD = (drawSig2) ? av_dChi2.at(iSamp).at(i2SigD) : 0.;
      const double dChi2SigU = (drawSig2) ? av_dChi2.at(iSamp).at(i2SigU) : 0.;

      // make the needed vectors for the points of the band graphs - zeros, verticals and such
      const std::array<double, 2> a_opMin = {opMin, opMin};
      const std::array<double, 2> a_1SigD = {err1SigD, err1SigD};
      const std::array<double, 2> a_1SigU = {err1SigU, err1SigU};
      const std::array<double, 2> a_2SigD = {err2SigD, err2SigD};
      const std::array<double, 2> a_2SigU = {err2SigU, err2SigU};

      std::cout << "Best fit result " << opName << ", sample " << samp << ": " << opMin << 
        " - " << err1SigD << " + " << err1SigU << " (1 sigma) - " << 
        err2SigD << " + " << err2SigU << " (2 sigma)" << std::endl;

      std::cout << "Best fit dChi2 "<< opName << ", sample " << samp << ": 0 - " << 
        dChi1SigD << " + " << dChi1SigU << " (1 sigma) - " << 
        dChi2SigD << " + " << dChi2SigU << " (2 sigma)" << std::endl;

      if (drawSig1) {
        // guiding line graph
        ag_sigma.at(iSamp).at(0) = std::make_unique<TGraphAsymmErrors>(2, a_opMin.data(), a_vPnt.data(), 
                                                                       a_zero.data(), a_zero.data(), a_zero.data(), a_zero.data());
        stylePlot(ag_sigma.at(iSamp).at(0).get(), kBlack, 1., 0, 20, 1.5, 8, 2);
        axisPlot(ag_sigma.at(iSamp).at(0).get(), 
                 yRange.at(0), yRange.at(1), 505, "#Delta #chi^{2}", 0.043, 1.21, 0.037, 
                 xRange.at(0), xRange.at(1), 505, opLeg.c_str(), 0.043, 1.11, 0.037);
        ag_sigma.at(iSamp).at(0)->SetName((opName + "_sigma0_" + toStr(samp)).c_str());
        ag_sigma.at(iSamp).at(0)->SetTitle("");

        // 1 sigma graph
        ag_sigma.at(iSamp).at(1) = std::make_unique<TGraphAsymmErrors>(2, a_opMin.data(), a_vPnt.data(), 
                                                                       a_1SigD.data(), a_1SigU.data(), a_vErrD.data(), a_vErrU.data());
        stylePlot(ag_sigma.at(iSamp).at(1).get(), kGreen + 1, 1., 1001, 0, 1.5, 1, 2);
        axisPlot(ag_sigma.at(iSamp).at(1).get(), 
                 yRange.at(0), yRange.at(1), 505, "#Delta #chi^{2}", 0.043, 1.21, 0.037, 
                 xRange.at(0), xRange.at(1), 505, opLeg.c_str(), 0.043, 1.11, 0.037);
        ag_sigma.at(iSamp).at(1)->SetName((opName + "_sigma1_" + toStr(samp)).c_str());
        ag_sigma.at(iSamp).at(1)->SetTitle("");

        if (drawSig2) {
          // 2 sigma graph
          ag_sigma.at(iSamp).at(2) = std::make_unique<TGraphAsymmErrors>(2, a_opMin.data(), a_vPnt.data(), 
                                                                         a_2SigD.data(), a_2SigU.data(), a_vErrD.data(), a_vErrU.data());
          axisPlot(ag_sigma.at(iSamp).at(2).get(), 
                   yRange.at(0), yRange.at(1), 505, "#Delta #chi^{2}", 0.043, 1.21, 0.037, 
                   xRange.at(0), xRange.at(1), 505, opLeg.c_str(), 0.043, 1.11, 0.037);
          stylePlot(ag_sigma.at(iSamp).at(2).get(), kOrange, 1., 1001, 0, 1.5, 1, 2);
          ag_sigma.at(iSamp).at(2)->SetName((opName + "_sigma2_" + toStr(samp)).c_str());
          ag_sigma.at(iSamp).at(2)->SetTitle("");

          if (samp == v_sample.front())
            leg->AddEntry(ag_sigma.at(iSamp).at(2).get(), ("#pm 2#sigma" + a_sampLeg.at(iSamp)).c_str(), "f");
        }

        if (samp == v_sample.front())
          leg->AddEntry(ag_sigma.at(iSamp).at(1).get(), ("#pm 1#sigma" + a_sampLeg.at(iSamp)).c_str(), "f");
      }
    }

    setH1Style();
    std::unique_ptr<TCanvas> can = std::make_unique<TCanvas>("can", "can", 200, 10, 1000, 1000);
    can->SetTopMargin(0.045);
    can->SetBottomMargin(0.11);
    can->SetLeftMargin(0.11);
    can->SetRightMargin(0.025);

    // for tacking on some text on the plot
    const std::string topMid = "#Lambda = " + toStr(eftLambda) + " TeV";
    TLatex txt;
    //txt.SetTextSize(0.039); // standard
    txt.SetTextSize(0.043); // cms
    txt.SetTextAlign(13);

    can->cd();

    styleLegend(leg.get(), 1, 0, 0, 42, 0.041, "");
    //putLegend(leg.get(), 0.685, 0.945, 0.655, 0.935); // top right
    //putLegend(leg.get(), 0.125, 0.385, 0.655, 0.935); // top left
    putLegend(leg.get(), 0.685, 0.945, 0.155, 0.435); // bottom right

    ag_dChi2.at(static_cast<int>(v_sample.front()))->Draw("a c");
    leg->Draw();
    txt.DrawLatexNDC(0.461, 0.915, topMid.c_str());

    if (drawSig1) {
      if (drawSig2)
        ag_sigma.at(static_cast<int>(v_sample.front())).at(2)->Draw("2");

      ag_sigma.at(static_cast<int>(v_sample.front())).at(1)->Draw("2");
      ag_sigma.at(static_cast<int>(v_sample.front())).at(0)->Draw("l");
    }

    for (const auto &g_dChi2 : ag_dChi2) {
      if (g_dChi2 == nullptr) continue;
      g_dChi2->Draw("c");
    }

    can->RedrawAxis();
    can->SaveAs((plotName + opName + "_dChi2.pdf").c_str());
    //can->SaveAs((plotName + opName + "_dChi2.C").c_str());

    // save them into a file
    auto file = std::make_unique<TFile>((plotName + opName + "_dChi2.root").c_str(), "recreate");
    file->cd();
    //can->Write();
    for (const auto &g_dChi2 : ag_dChi2) {
      if (g_dChi2 == nullptr) continue; 
      g_dChi2->Write();
    }
    for (const auto &ag_sig : ag_sigma) {
      for (const auto &g_sig : ag_sig) {
        if (g_sig == nullptr) continue; 
        g_sig->Write();
      }
    }

    std::cout << std::endl;
  }
}



void EFTFitter::draw2DChi2(const std::map<std::array<std::string, 2>,
                           std::array<std::pair<std::string, std::array<double, 2>>, 2>> &mt_opPair,
                           const std::string &plotName, const std::vector<Sample> &v_sample,
                           const double &dChi2FracScan) const
{
  if (m_fitChi2.empty())
    throw std::logic_error( "This method shouldn't be called before the computeFitChi2() method!!" );

  if (v_opName.size() < 2) {
    std::cout << "Method does nothing with less than 2 operators, skipping..." << std::endl << std::endl;
    return;
  }

  std::cout << "Drawing the 2D dChi2 graph for each operator pair..." << std::endl << std::endl;

  // bind the operator-chi2 set in a map; consider only those with at most 2 op != 0
  std::map< std::pair<std::vector<double>, Sample>, double> m_opMax2Non0;
  for (const auto &p_fitChi2 : m_fitChi2) {
    const std::vector<double> v_opVal = extractValue( parseOpFromKey(p_fitChi2.first.first) );

    // is there a SM ie all operators are 0 case?
    if (std::all_of(std::begin(v_opVal), std::end(v_opVal), [] (const double &opV) {return opV == 0.;}))
      m_opMax2Non0.insert( {{v_opVal, p_fitChi2.first.second}, p_fitChi2.second} );

    // is there a case where exactly 1 or 2 operator != 0?
    const int nOpNon0 = std::count_if(std::begin(v_opVal), std::end(v_opVal), [] (const auto &opV) {return opV != 0.;});
    if (nOpNon0 == 1 or nOpNon0 == 2)
      m_opMax2Non0.insert( {{v_opVal, p_fitChi2.first.second}, p_fitChi2.second} );
  }

  // declare here things that are gonna be applicable for all anyway
  const std::array<std::string, 2> a_sampLeg = {" (all)", " (lin)"};
  const std::array<int, 2> a_sampCol0 = {kPink - 1, kAzure - 1};
  const std::array<int, 2> a_sampCol1 = {kPink - 3, kAzure - 3}; // brazil is kGreen + 1
  const std::array<int, 2> a_sampCol2 = {kPink - 4, kAzure - 4}; // and kOrange
  const std::array<int, 2> a_sampStyF = {0, 0}; // 0 means no fill, 1001 means full fill
  const std::array<int, 2> a_sampStyL = {1, 1};
  const std::array<int, 2> a_sampStyM = {kFullCrossX, kOpenCrossX};
  const std::array<std::string, 2> a_sampOpt = {"lf", "lf"};
  const std::array<double, 2> a_zero = {0., 0.};

  // nBin - (1 for shape bin drop) + (1 if hybrid has rate) - nOp
  const int iXs = (statMode == Stat::xsec) ? 0 : 1;
  const int nDoF = v_rawBin.size() - 1 + v_rawBin.empty()
    - ((fitMode == Fit::shape) * int(shapeSum))
    + (fitMode == Fit::hybrid and m_binContent.at({dataName, Sample::all}).at(iXs).at(0) > 0.)
    - 2;

  // for the sorting of points in op1 - op2 space - distance calculator
  const auto f_sqDist = [] (const std::array<double, 2> &p1, const std::array<double, 2> &p2) {
    const double d1 = p2.at(0) - p1.at(0), d2 = p2.at(1) - p1.at(1);
    return (d1 * d1) + (d2 * d2);
  };

  // and the actual sorter - such that 2nd element is closest to 1st, N+1'th is closest to N'th...
  const auto f_sortNearestNeighbor = [&f_sqDist] (std::vector<std::array<double, 2>> &v_opPnt) {
    if (v_opPnt.size() < 3) return;

    for (auto iPnt1 = std::begin(v_opPnt); iPnt1 != std::end(v_opPnt) - 2; ++iPnt1) {
      double sqDist = f_sqDist(*iPnt1, *(std::next(iPnt1)));

      for (auto iPnt2 = std::next(iPnt1, 2); iPnt2 != std::end(v_opPnt); ++iPnt2) {
        if (f_sqDist(*iPnt1, *iPnt2) < sqDist) {
          sqDist = f_sqDist(*iPnt1, *iPnt2);
          std::iter_swap(iPnt2, std::next(iPnt1));
        }
      }
    }
  };

  // ok let's get it rolling...
  for (const auto &p_opPair : mt_opPair) {
    const auto &op1Name = p_opPair.first.at(0);
    const int iOp1 = std::distance(std::begin(v_opName), std::find(std::begin(v_opName), std::end(v_opName), op1Name));
    const auto &op1Leg = p_opPair.second.at(0).first;
    const auto &op1Range = p_opPair.second.at(0).second;

    const auto &op2Name = p_opPair.first.at(1);
    const int iOp2 = std::distance(std::begin(v_opName), std::find(std::begin(v_opName), std::end(v_opName), op2Name));
    const auto &op2Leg = p_opPair.second.at(1).first;
    const auto &op2Range = p_opPair.second.at(1).second;

    // now we make the points to be plotted - separately for all and linear
    std::array<std::vector<std::array<double, 2>>, 2> av_opVal;
    std::array<std::vector<double>, 2> av_opChi2, av_dChi2;
    for (const auto &p_op2 : m_opMax2Non0) {
      const auto &v_opVal = p_op2.first.first;
      const auto &opChi2 = p_op2.second;
      const auto &samp = p_op2.first.second;
      const int iSamp = static_cast<int>(samp);

      // skip if it's not the op under consideration
      const int nOpNon0 = std::count_if(std::begin(v_opVal), std::end(v_opVal), [] (const auto &opV) {return opV != 0.;});
      if ((nOpNon0 == 1 and (v_opVal.at(iOp1) == 0 and v_opVal.at(iOp2) == 0)))
        continue;
      if ((nOpNon0 == 2 and (v_opVal.at(iOp1) == 0 or v_opVal.at(iOp2) == 0)))
        continue;

      av_opVal.at(iSamp).push_back({v_opVal.at(iOp1), v_opVal.at(iOp2)});
      av_opChi2.at(iSamp).push_back(opChi2);
    }

    // and we make the graphs
    std::array<std::unique_ptr<TGraph>, 2> ag_sigma0 = {nullptr, nullptr};
    std::array<std::unique_ptr<TGraph>, 2> ag_sigma1 = {nullptr, nullptr};
    std::array<std::unique_ptr<TGraph>, 2> ag_sigma2 = {nullptr, nullptr};
    auto leg = std::make_unique<TLegend>();
    for (auto &samp : v_sample) {
      const int iSamp = static_cast<int>(samp);
      if (av_opVal.at(iSamp).empty()) continue;

      // grab the index of the min chi2 and make the dChi2 vector
      const int iMin = std::distance(std::begin(av_opChi2.at(iSamp)), 
                                     std::min_element(std::begin(av_opChi2.at(iSamp)), std::end(av_opChi2.at(iSamp))));
      const auto opMin = av_opVal.at(iSamp).at(iMin);
      const double chi2Min = av_opChi2.at(iSamp).at(iMin);
      for (double &chi2 : av_opChi2.at(iSamp)) 
        av_dChi2.at(iSamp).push_back(chi2 - chi2Min);

      // make the best fit graph
      const std::vector<double> v_op1Sig0(1, opMin.at(0)), v_op2Sig0(1, opMin.at(1));
      ag_sigma0.at(iSamp) = std::make_unique<TGraph>(1, v_op2Sig0.data(), v_op1Sig0.data());
      stylePlot(ag_sigma0.at(iSamp).get(), a_sampCol0.at(iSamp), 1., 0, a_sampStyM.at(iSamp), 2.5, a_sampStyL.at(iSamp), 3);
      axisPlot(ag_sigma0.at(iSamp).get(), 
               op1Range.at(0), op1Range.at(1), 505, op1Leg.c_str(), 0.043, 1.21, 0.037, 
               op2Range.at(0), op2Range.at(1), 505, op2Leg.c_str(), 0.043, 1.11, 0.037);
      ag_sigma0.at(iSamp)->SetName((op1Name + "_" + op2Name + "_sigma0_" + toStr(samp)).c_str());
      ag_sigma0.at(iSamp)->SetTitle("");
      leg->AddEntry(ag_sigma0.at(iSamp).get(), ("Best fit" + a_sampLeg.at(iSamp)).c_str(), "p");

      const double chi2Prob = TMath::Prob(chi2Min, nDoF);
      std::cout << "Operator pair " << op1Name << " - " << op2Name << ", sample " << samp <<
        ": Best fit result is " << opMin.at(0) << ", " << opMin.at(1) << 
        " with chi2/nDoF " << chi2Min << "/" << nDoF << ", with p-value " << chi2Prob << std::endl;

      // and the points going into 1, 2 sigma contours
      // threshold is dChi2 containing 1 and 2 sigma for chi2 with 2 ndf
      // obtained with 1. - TMath::Prob(raw_chi2, ndf)
      const double dChi2Sig1 = 2.295743, dChi2Sig2 = 6.180063;
      std::vector<std::array<double, 2>> v_opSig1, v_opSig2;
      for (int iD = 0; iD < av_dChi2.at(iSamp).size(); ++iD) {
        if (av_dChi2.at(iSamp).at(iD) <= dChi2Sig2) {
          // the > cut is there because we're interested only in the boundary for plotting
          // and removing points that are certainly not the boundary makes things faster
          const bool scanSig1 = (std::abs(dChi2FracScan) >= 1. or av_dChi2.at(iSamp).at(iD) > (1. - dChi2FracScan) * dChi2Sig1);
          const bool scanSig2 = (std::abs(dChi2FracScan) >= 1. or av_dChi2.at(iSamp).at(iD) > (1. - dChi2FracScan) * dChi2Sig2);

          if (scanSig1 and av_dChi2.at(iSamp).at(iD) <= dChi2Sig1)
            v_opSig1.push_back(av_opVal.at(iSamp).at(iD));
          else if (scanSig2)
            v_opSig2.push_back(av_opVal.at(iSamp).at(iD));
        }
      }

      // ok now we need to define a function to weed out all interior points of the contour
      // idea is to check the op1, op2 neighbors of each point and count those still within the contour
      // an interior point would have 4 (all) neighbors still within, while edge 0 < n < 4
      // we need it to be a unary predicate in order for it to be usable with remove_if
      // and so the threshold is put outside
      // similarly the function needs to be here as it refers to av_opVal.at(iSamp) etc which is scoped here
      // probably can be made faster eg by keeping the points that are neighbors along the edge with the first edge point
      // but seems circular in that neighbors along edge is known only when the edge is known hmm...
      const int nPnt = av_opVal.at(iSamp).size();
      bool prev3Neighbor = false;
      double dChi2Cut = 0.;
      const auto f_isNotEdgePnt = [&] (const auto &opPnt) {
        // exploit the fact that points are sorted lexicographically
        // so loop finding the last low and first high neighbors
        int iHi1 = -1, iHi2 = -1, iLo1 = -1, iLo2 = -1;
        for (int iP = 0; iP < nPnt; ++iP) {
          if (iHi1 != -1 and iHi2 != -1) break;
          if (av_opVal.at(iSamp).at(iP).at(0) != opPnt.at(0) and av_opVal.at(iSamp).at(iP).at(1) != opPnt.at(1)) continue;

          if (iHi1 == -1 and av_opVal.at(iSamp).at(iP).at(1) == opPnt.at(1)) {
            if (av_opVal.at(iSamp).at(iP).at(0) < opPnt.at(0)) {
              iLo1 = iP;
              continue;
            }

            if (av_opVal.at(iSamp).at(iP).at(0) > opPnt.at(0)) {
              iHi1 = iP;
              continue;
            }
          }

          if (iHi2 == -1 and av_opVal.at(iSamp).at(iP).at(0) == opPnt.at(0)) {
            if (av_opVal.at(iSamp).at(iP).at(1) < opPnt.at(1)) {
              iLo2 = iP;
              continue;
            }

            if (av_opVal.at(iSamp).at(iP).at(1) > opPnt.at(1)) {
              iHi2 = iP;
              continue;
            }
          }
        }

        if (iHi1 == -1 or iHi2 == -1 or iLo1 == -1 or iLo2 == -1) {
          std::cout << "Op pair " << op1Name << " - " << op2Name << " type " << samp << " point " << opPnt <<
            ": unable to define all 4 neighbor points! Dropping..." << std::endl;
          return true;
        }

        // now we simply check how many of the neighbors are within the contour
        const int nNeighborIn = (av_dChi2.at(iSamp).at(iHi1) <= dChi2Cut) + (av_dChi2.at(iSamp).at(iLo1) <= dChi2Cut)
        + (av_dChi2.at(iSamp).at(iHi2) <= dChi2Cut) + (av_dChi2.at(iSamp).at(iLo2) <= dChi2Cut);

        // it really shouldn't be 0 since we're running over points within contour, scream if so
        if (!nNeighborIn) {
          std::cout << "Op pair " << op1Name << " - " << op2Name << " type " << samp << " point " << opPnt <<
            " has 0 neighbors within contour! Dropping..." << std::endl;
          return true;
        }

        if (nNeighborIn == 3)
          prev3Neighbor = !prev3Neighbor;

        // only every other point with 3 neighbors inside the contour are considered on-edge
        // since while they also lie on the edge, in this case the edge is purely vertical or horizontal
        // so they don't change the bounded area
        return nNeighborIn == 4 or (nNeighborIn == 3 and prev3Neighbor);
      };

      // ok now actually apply the function to weed out the points
      dChi2Cut = dChi2Sig1;
      std::cout << "dChi2 <= " << dChi2Cut << " contour has " << v_opSig1.size() << " points. Removing the insides..." << std::endl;
      v_opSig1.erase( std::remove_if(std::begin(v_opSig1), std::end(v_opSig1), f_isNotEdgePnt), std::end(v_opSig1) );
      std::cout << "Done. " << v_opSig1.size() << " edge points remain in 1 sigma contour." << std::endl;
      f_sortNearestNeighbor(v_opSig1);
      //printAll(v_opSig1);

      dChi2Cut = dChi2Sig2;
      std::cout << "dChi2 <= " << dChi2Cut << " contour has " << v_opSig2.size() << " points. Removing the insides..." << std::endl;
      v_opSig2.erase( std::remove_if(std::begin(v_opSig2), std::end(v_opSig2), f_isNotEdgePnt), std::end(v_opSig2) );
      std::cout << "Done. " << v_opSig2.size() << " edge points remain in 2 sigma contour." << std::endl;
      f_sortNearestNeighbor(v_opSig2);

      /*/ insert at back the first point so as to have a closed shape
      v_opSig1.push_back(v_opSig1.front());
      v_opSig2.push_back(v_opSig2.front());
      */

      // ok now we split them to be passed to the graph
      std::vector<double> v_op1Sig1, v_op2Sig1, v_op1Sig2, v_op2Sig2;
      for (int iP1 = 0; iP1 < v_opSig1.size(); ++iP1) {
        v_op1Sig1.push_back(v_opSig1.at(iP1).at(0));
        v_op2Sig1.push_back(v_opSig1.at(iP1).at(1));
      }
      for (int iP2 = 0; iP2 < v_opSig2.size(); ++iP2) {
        v_op1Sig2.push_back(v_opSig2.at(iP2).at(0));
        v_op2Sig2.push_back(v_opSig2.at(iP2).at(1));
      }

      ag_sigma1.at(iSamp) = std::make_unique<TGraph>(v_op1Sig1.size(), v_op2Sig1.data(), v_op1Sig1.data());
      stylePlot(ag_sigma1.at(iSamp).get(), a_sampCol1.at(iSamp), 0.43, a_sampStyF.at(iSamp), 0, 2.5, a_sampStyL.at(iSamp), 4);
      axisPlot(ag_sigma1.at(iSamp).get(), 
               op1Range.at(0), op1Range.at(1), 505, op1Leg.c_str(), 0.043, 1.21, 0.037, 
               op2Range.at(0), op2Range.at(1), 505, op2Leg.c_str(), 0.043, 1.11, 0.037);
      ag_sigma1.at(iSamp)->SetName((op1Name + "_" + op2Name + "_sigma1_" + toStr(samp)).c_str());
      ag_sigma1.at(iSamp)->SetTitle("");

      ag_sigma2.at(iSamp) = std::make_unique<TGraph>(v_op1Sig2.size(), v_op2Sig2.data(), v_op1Sig2.data());
      stylePlot(ag_sigma2.at(iSamp).get(), a_sampCol2.at(iSamp), 0.43, a_sampStyF.at(iSamp), 0, 2.5, a_sampStyL.at(iSamp), 4);
      axisPlot(ag_sigma2.at(iSamp).get(), 
               op1Range.at(0), op1Range.at(1), 505, op1Leg.c_str(), 0.043, 1.21, 0.037, 
               op2Range.at(0), op2Range.at(1), 505, op2Leg.c_str(), 0.043, 1.11, 0.037);
      ag_sigma2.at(iSamp)->SetName((op1Name + "_" + op2Name + "_sigma2_" + toStr(samp)).c_str());
      ag_sigma2.at(iSamp)->SetTitle("");
    }

    // and now we make the line to guide the eye along 0 in y and x
    std::array<std::unique_ptr<TGraph>, 2> ag_zero = {nullptr, nullptr};
    ag_zero.at(0) = std::make_unique<TGraph>(2, op2Range.data(), a_zero.data());
    stylePlot(ag_zero.at(0).get(), kGray + 2, 1., 0, 0, 2.5, 1, 1);
    axisPlot(ag_zero.at(0).get(), 
             op1Range.at(0), op1Range.at(1), 505, op1Leg.c_str(), 0.043, 1.21, 0.037, 
             op2Range.at(0), op2Range.at(1), 505, op2Leg.c_str(), 0.043, 1.11, 0.037);
    ag_zero.at(0)->SetName((op1Name + "_" + op2Name + "_xX_y0").c_str());
    ag_zero.at(0)->SetTitle("");

    ag_zero.at(1) = std::make_unique<TGraph>(2, a_zero.data(), op1Range.data());
    stylePlot(ag_zero.at(1).get(), kGray + 2, 1., 0, 0, 2.5, 1, 1);
    axisPlot(ag_zero.at(1).get(), 
             op1Range.at(0), op1Range.at(1), 505, op1Leg.c_str(), 0.043, 1.21, 0.037, 
             op2Range.at(0), op2Range.at(1), 505, op2Leg.c_str(), 0.043, 1.11, 0.037);
    ag_zero.at(1)->SetName((op1Name + "_" + op2Name + "_x0_yY").c_str());
    ag_zero.at(1)->SetTitle("");

    setH1Style();
    std::unique_ptr<TCanvas> can = std::make_unique<TCanvas>("can", "can", 200, 10, 1000, 1000);
    can->SetTopMargin(0.045);
    can->SetBottomMargin(0.11);
    can->SetLeftMargin(0.11);
    can->SetRightMargin(0.025);

    // for tacking on some text on the plot
    const std::string topMid = "#Lambda = " + toStr(eftLambda) + " TeV";
    TLatex txt;
    txt.SetTextSize(0.039);
    txt.SetTextAlign(13);

    can->cd();

    styleLegend(leg.get(), 1, 0, 0, 42, 0.041, "");
    //putLegend(leg.get(), 0.685, 0.945, 0.755, 0.935); // top right
    //putLegend(leg.get(), 0.125, 0.385, 0.755, 0.935); // top left
    putLegend(leg.get(), 0.685, 0.945, 0.155, 0.335); // bottom right

    ag_sigma0.at(static_cast<int>(v_sample.front()))->Draw("a p");
    leg->Draw();
    txt.DrawLatexNDC(0.461, 0.923, topMid.c_str());

    ag_zero.at(0)->Draw("l");
    ag_zero.at(1)->Draw("l");

    for (const auto &samp : v_sample) {
      const int iSamp = static_cast<int>(samp);

      if (ag_sigma2.at(iSamp) != nullptr)
        ag_sigma2.at(iSamp)->Draw(a_sampOpt.at(iSamp).c_str());

      if (ag_sigma1.at(iSamp) != nullptr)
        ag_sigma1.at(iSamp)->Draw(a_sampOpt.at(iSamp).c_str());
    }
    for (const auto &g_sigma0 : ag_sigma0) {
      if (g_sigma0 == nullptr) continue;
      g_sigma0->Draw("p");
    }

    can->RedrawAxis();
    can->SaveAs((plotName + op1Name + "_" + op2Name + "_dChi2.pdf").c_str());
    //can->SaveAs((plotName + op1Name + "_" + op2Name + "_dChi2.C").c_str());

    // save them into a file
    auto file = std::make_unique<TFile>((plotName + op1Name + "_" + op2Name + "_dChi2.root").c_str(), "recreate");
    file->cd();
    //can->Write();
    for (const auto &samp : v_sample) {
      const int iSamp = static_cast<int>(samp);

      if (ag_sigma2.at(iSamp) != nullptr)
        ag_sigma2.at(iSamp)->Write();

      if (ag_sigma1.at(iSamp) != nullptr)
        ag_sigma1.at(iSamp)->Write();
    }
    for (const auto &g_sigma0 : ag_sigma0) {
      if (g_sigma0 == nullptr) continue;
      g_sigma0->Write();
    }

    std::cout << std::endl;
  }
}



void EFTFitter::clearContent(const int clearLevel)
{
  std::cout << "Clearing up contents of the current EFTFitter object..." << std::endl << std::endl;

  m_fitChi2.clear();
  v_keyToFit.clear();

  std::cout << "Chi2 values and corresponding keys cleared." << std::endl << std::endl;
  if (clearLevel > 0) return;

  m_op1Eq1.clear();
  m_op2Eq1.clear();

  m_binContent.clear();
  hasData = false;

  m_covMat.clear();
  f_finalCov = nullptr;

  v_opName.clear();
  v_rawBin.clear();
  f_hybridTransform = nullptr;

  std::cout << "All EFTFitter object content cleared." << std::endl << std::endl;
}




/***
 * EFTFitter source (private)
 ***/



void EFTFitter::autoCovMatStatBin()
{
  if (fitMode == Fit::hybrid) {
    std::cout << "Automatic matrix generation not supported in Fit::hybrid mode. Aborting..." << std::endl;
    return;
  }

  if (!hasData) {
    std::cout << "No data content available, aborting..." << std::endl << std::endl;
    return;
  }

  const std::vector<std::array<double, 2>> &dataC = m_binContent.at({dataName, Sample::all});
  const double &integral = (fitMode == Fit::shape) ? shapeSum : (statMode == Stat::xsec) ? dataC.at(0).at(0) : dataC.at(1).at(0);
  const int nBin = dataC.size();

  TMatrixD tmpMat(nBin - 2, nBin - 2);
  for (int iR = 2; iR < nBin; ++iR ) {
    for (int iC = 2; iC < nBin; ++iC )
      tmpMat(iR - 2, iC - 2) = integral * integral * dataC.at(iR).at(1) * dataC.at(iC).at(1);
  }

  m_covMat.insert({"statBin", tmpMat});

  // this is how to construct (in-place or otherwise) a product of 2 matrices
  //m_covMat.insert({"statErr", TMatrixD(m_covMat.at("statBin"), TMatrixD::kMult, m_covMat.at("statCorr"))});

  // but in this case we want just element-wise multiplication
  TMatrixD eleMat = m_covMat.at("statBin");
  ElementMult(eleMat, m_covMat.at("statCorr"));
  m_covMat.insert({"statErr", eleMat});
}



bool EFTFitter::checkInputBin(const std::unique_ptr<TH1D> &hist) const
{
  // compare the binning against the known template binning (obtained at first feeding into the map)
  if (fitMode != Fit::hybrid)
    return (hist->GetNbinsX() == v_rawBin.size() - 1 and FitUtil::extractBin( hist.get() ) == v_rawBin);
  else {
    // yes Fit::hybrid is AWFUL
    std::vector<std::array<double, 2>> v_pre = FitUtil::extractContentError(hist.get());
    v_pre.push_back({0., 0.});
    v_pre.push_back({0., 0.});

    std::vector<std::array<double, 2>> v_pos = f_hybridTransform(v_pre);
    const auto v_bin = makeInterval(0., double(v_pos.size() - 2), 1.);

    return v_bin == v_rawBin;
  }
}



double EFTFitter::getContentSum(const std::vector<std::array<double, 2>> &v_binContent) const
{
  if (fitMode == Fit::shape)
    return shapeSum;

  if (fitMode == Fit::absolute) {
    if (statMode == Stat::xsec)
      return v_binContent.at(0).at(0);

    return v_binContent.at(1).at(0);
  }

  // Fit::hybrid, assume nothing was/needs to be done
  return 1.;
}



void EFTFitter::normalizeContent(std::vector<std::array<double, 2>> &v_binContent, const int iXs)
{
  // first check if it's already normalized
  const double normalizedSum = std::accumulate(std::next(std::begin(v_binContent), 2), std::end(v_binContent), 0., 
                                               [] (const double &sum, const auto &con) { return sum + con.at(0); });
  if (normalizedSum == 1.) return;

  const int nBin = v_binContent.size();
  for (int iB = 2; iB < nBin; ++iB) {
    v_binContent.at(iXs).at(0) += v_binContent.at(iB).at(0);
    v_binContent.at(iXs).at(1) += v_binContent.at(iB).at(1) * v_binContent.at(iB).at(1);
  }
  v_binContent.at(iXs).at(1) = std::sqrt(v_binContent.at(iXs).at(1));

  for (int iB = 2; iB < nBin; ++iB) {
    v_binContent.at(iB).at(0) = v_binContent.at(iB).at(0) / v_binContent.at(iXs).at(0);
    v_binContent.at(iB).at(1) = v_binContent.at(iB).at(1) / v_binContent.at(iXs).at(0);
  }
}



void EFTFitter::assignInSituXsec(std::vector<std::array<double, 2>> &v_binContent)
{
  // nothing to scale if it has both numbers already
  if (v_binContent.at(0).at(0) != 0. and v_binContent.at(1).at(0) != 0.) return;

  if (v_binContent.at(0).at(0) == 0.) {
    v_binContent.at(0).at(0) = v_binContent.at(1).at(0) / FitUtil::intLumi;
    v_binContent.at(0).at(1) = v_binContent.at(1).at(1) / FitUtil::intLumi;
  }

  if (v_binContent.at(1).at(0) == 0.) {
    v_binContent.at(1).at(0) = v_binContent.at(0).at(0) * FitUtil::intLumi;
    v_binContent.at(1).at(1) = v_binContent.at(0).at(1) * FitUtil::intLumi;
  }
}



std::vector<std::array<double, 2>> EFTFitter::interpolateOpValue(const std::string &keyName, const Sample sampleType, 
                                                                 const bool varyWithinError)
{
  // return the raw input if it is available and don't need to be varied
  if (m_binContent.count({keyName, sampleType}) and !varyWithinError)
    return (fitMode != Fit::hybrid or keyName == dataName) ? 
      m_binContent.at({keyName, sampleType}) : f_hybridTransform(m_binContent.at({keyName, sampleType}));

  // otherwise check if the ingredients are there and hack away
  if (m_op1Eq1.empty() or (v_opName.size() > 1 and m_op2Eq1.empty()))
    throw std::logic_error( "This method shouldn't be called given the insuffiencient input! "
                            "Ensure prepareInterpolationBase() has been called!!!");

  // these are gonna pop up often
  const int nBin = (std::begin(m_binContent)->first.first != dataName) ? 
    std::begin(m_binContent)->second.size() : std::next(std::begin(m_binContent))->second.size();
  const int iXs = (statMode == Stat::xsec) ? 0 : 1;
  const double sqLambda = eftLambda * eftLambda; 
  const double quLambda = sqLambda * sqLambda;

  // ok let's start interpolating - first prepare the bin contents starting from SM
  const auto v_xs_op0 = std::find_if(std::begin(m_binContent), std::end(m_binContent), [this] (const auto &p) {
      const auto smOp = extractValue( this->parseOpFromKey(p.first.first) );
      return !smOp.empty() and std::all_of(std::begin(smOp), std::end(smOp), [] (const auto &op) {return op == 0.;});
    });

  // get the op-coeff set represented by this key
  const std::map<std::string, double> m_opVal = parseOpFromKey(keyName);

  std::vector<std::string> v_opNon0Name;
  std::vector<double> v_opNon0Val;
  for (const auto& p_opV : m_opVal) {
    if (p_opV.second == 0.) continue;

    v_opNon0Name.push_back(p_opV.first);
    v_opNon0Val.push_back(p_opV.second);
  }

  // and we build the bin content for this key
  std::vector<std::array<double, 2>> v_intResult(nBin, {0., 0.});
  for (int iB = 2; iB < nBin; ++iB) {
    // SM xsec, error in bin iB
    const double xs_op0 = v_xs_op0->second.at(iXs).at(0) * v_xs_op0->second.at(iB).at(0);
    const double xe_op0 = v_xs_op0->second.at(iXs).at(0) * v_xs_op0->second.at(iB).at(1);

    // get the int, res xsec, error for the given op-coeff set
    double xs_opI = 0., xe_opI = 0., xs_opR = 0., xe_opR = 0.;
    for (int iOp1 = 0; iOp1 < v_opNon0Val.size(); ++iOp1) {
      const std::string &op1Name = v_opNon0Name.at(iOp1);
      const double &op1Val = v_opNon0Val.at(iOp1);
      const double &op1ConI = m_op1Eq1.at({op1Name, Sample::linear}).at(iB).at(0);
      const double &op1ErrI = m_op1Eq1.at({op1Name, Sample::linear}).at(iB).at(1);
      const double &op1ConR = m_op1Eq1.at({op1Name, Sample::quadratic}).at(iB).at(0);
      const double &op1ErrR = m_op1Eq1.at({op1Name, Sample::quadratic}).at(iB).at(1);

      xs_opI += op1Val * op1ConI;
      xe_opI += (op1Val * op1ErrI) * (op1Val * op1ErrI);

      xs_opR += op1Val * op1Val * op1ConR;
      xe_opR += (op1Val * op1ErrR) * (op1Val * op1ErrR);

      for (int iOp2 = iOp1 + 1; iOp2 < v_opNon0Val.size(); ++iOp2) {
        const std::string &op2Name = v_opNon0Name.at(iOp2);
        const double &op2Val = v_opNon0Val.at(iOp2);
        const double &op12Sum = m_op2Eq1.at({op1Name, op2Name}).at(iXs).at(0);
        const double &op12Con = op12Sum * m_op2Eq1.at({op1Name, op2Name}).at(iB).at(0);
        const double &op12Err = op12Sum * m_op2Eq1.at({op1Name, op2Name}).at(iB).at(1);

        xs_opR += op1Val * op2Val * op12Con;
        xe_opR += (op1Val * op2Val * op12Err) * (op1Val * op2Val * op12Err);
      }
    }
    xe_opI = std::sqrt(xe_opI);
    xe_opR = std::sqrt(xe_opR);

    // and this becomes the content of this bin (int, res are both sum of itself and SM, ignoring the other)
    if (sampleType == Sample::all) {
      v_intResult.at(iB).at(0) = xs_op0 + (xs_opI / sqLambda) + (xs_opR / quLambda);
      v_intResult.at(iB).at(1) = std::sqrt((xe_op0 * xe_op0) + 
                                           ((xe_opI / sqLambda) * (xe_opI / sqLambda)) + 
                                           ((xe_opR / quLambda) * (xe_opR / quLambda)));

      if (varyWithinError) {
        // get variation of SM + linear and quadratic separately
        // done this way since linear might be negative, so can't take it independently
        const double scale = (statMode == Stat::xsec) ? FitUtil::intLumi : 1.;
        const double count0I = scale * (xs_op0 + (xs_opI / sqLambda));
        const double countR = scale * (xs_opR / quLambda);

        v_intResult.at(iB).at(0) -= (getDPoissonVariation(count0I) + getDPoissonVariation(countR)) / scale;
      }
    }
    else if (sampleType == Sample::linear) {
      v_intResult.at(iB).at(0) = xs_op0 + (xs_opI / sqLambda);
      v_intResult.at(iB).at(1) = std::sqrt((xe_op0 * xe_op0) + ((xe_opI / sqLambda) * (xe_opI / sqLambda)));

      if (varyWithinError) {
        // comment at all part
        const double scale = (statMode == Stat::xsec) ? FitUtil::intLumi : 1.;
        const double count0I = scale * (xs_op0 + (xs_opI / sqLambda));

        v_intResult.at(iB).at(0) -= getDPoissonVariation(count0I) / scale;
      }
    }
    else if (sampleType == Sample::quadratic) {
      v_intResult.at(iB).at(0) = xs_op0 + (xs_opR / quLambda);
      v_intResult.at(iB).at(1) = std::sqrt((xe_op0 * xe_op0) + ((xe_opR / quLambda) * (xe_opR / quLambda)));

      if (varyWithinError) {
        // comment at all part
        const double scale = (statMode == Stat::xsec) ? FitUtil::intLumi : 1.;
        const double count0 = scale * xs_op0;
        const double countR = scale * (xs_opR / quLambda);

        v_intResult.at(iB).at(0) -= (getDPoissonVariation(count0) + getDPoissonVariation(countR)) / scale;
      }
    }
  }

  this->normalizeContent(v_intResult, iXs);
  this->assignInSituXsec(v_intResult);
  return (fitMode != Fit::hybrid) ? v_intResult : f_hybridTransform(v_intResult);
}



std::map<std::string, double> EFTFitter::parseOpFromKey(const std::string &keyName) const
{
  // assumes the syntax is op1_val1--op2_val2-- ... --opN_valN
  std::map<std::string, double> m_opVal;
  if (keyName == dataName) return m_opVal;

  int iBegin = 0, iDelim = keyName.find("--");

  while (true) {
    const std::string opStr = keyName.substr(iBegin, iDelim - iBegin);
    const int iUsc = opStr.find('_');
    const std::string opName = opStr.substr(0, iUsc);
    const double opVal = std::stod( opStr.substr(iUsc + 1) );

    m_opVal.insert( {opName, opVal} );
    if (iDelim > keyName.length()) break;
    else {
      iBegin = iDelim + 2;
      iDelim = keyName.find("--", iBegin);
    }
  }

  return m_opVal;
}



std::string EFTFitter::parseKeyFromOp(const std::map<std::string, double> &m_opVal)
{
  std::string keyName = "";
  for (auto &p_opVal : m_opVal) {
    keyName += (keyName == "") ? "" : "--";

    keyName += p_opVal.first + "_" + toStr(p_opVal.second);
  }

  return keyName;
}



std::string EFTFitter::fixKeyFormat(const std::string &keyName) const
{
  return (keyName == dataName) ? dataName : parseKeyFromOp( parseOpFromKey(keyName) );
}




bool EFTFitter::checkOpSet(const std::string &keyName, const Sample sampleType)
{
  // data checks
  if (keyName == dataName) {
    if (hasData) return false;
    if (sampleType != Sample::all) return false;
    hasData = true;
    return true;
  }

  // remove those trivially not following the syntax
  if (keyName.find("_") == std::string::npos)
    return false;

  const std::map<std::string, double> m_opVal = parseOpFromKey(keyName);

  // first time a not-data key is added, make list of operator name from it
  if (v_opName.empty()) {
    v_opName = extractKey( m_opVal );
    return true;
  }

  // check for name set match
  if (extractKey( m_opVal ) == v_opName)
    return true;

  // shouldnt reach this point; veto everything that does
  return false;
}



void EFTFitter::unpackOpGrid(std::vector<std::string> &v_opGrid, 
                             const std::map<std::string, std::vector<double>> &m_opGrid,
                             std::map<std::string, std::vector<double>>::const_iterator i_opGrid, 
                             const std::string &iniStr) const
{
  const int nOp = i_opGrid->second.size();

  for (int iOp = 0; iOp < nOp; ++iOp) {
    // keys with e- is likely a rounding issue of 0; fix that
    std::string valStr = toStr(i_opGrid->second.at(iOp));
    if (valStr.find("e-") != std::string::npos)
      valStr = "0";

    std::string varStr = i_opGrid->first + "_" + valStr;
    std::string binStr = (iniStr == "") ? varStr : iniStr + "--" + varStr;

    const auto v_opVal = extractValue( parseOpFromKey(binStr) );
    const int nOpNon0 = std::count_if(std::begin(v_opVal), std::end(v_opVal), [] (const auto &op) {return op != 0.;});

    // terminating condition - we're at the end
    if (std::next(i_opGrid) == std::end(m_opGrid))
      v_opGrid.push_back(binStr);
    // if we ever need to allow only N op != 0; early termination
    else if (nOpNon0 > 1) {
      for (auto iOpG = std::next(i_opGrid); iOpG != std::end(m_opGrid); ++iOpG)
        binStr = binStr + "--" + iOpG->first + "_0";
      v_opGrid.push_back(binStr);
    }
    else
      unpackOpGrid(v_opGrid, m_opGrid, std::next(i_opGrid), binStr);
  }
}



std::unique_ptr<TH1D> EFTFitter::convertContentToHist(const std::string &keyName, const Sample sampleType, const bool divideBinWidth)
{
  // grab the content to be converted
  const std::vector<std::array<double, 2>> v_binC = interpolateOpValue(keyName, sampleType);

  // and make the histograms
  const int nBin = v_rawBin.size() - 1;
  auto hist = std::make_unique<TH1D>((keyName + "_" + toStr(sampleType)).c_str(), "", nBin, v_rawBin.data());
  //std::cout << keyName << " " << nBin << " " << v_binC.size() << std::endl;

  const double integral = getContentSum(v_binC);

  for (int iB = 1; iB <= nBin; ++iB) {
    // a scaled histogram just have every bin content, error scaled down by a factor
    // so we simply have to unscale it back
    const double binCon = integral * v_binC.at(iB + 1).at(0);
    const double binErr = integral * v_binC.at(iB + 1).at(1);

    const double binWid = divideBinWidth ? hist->GetXaxis()->GetBinUpEdge(iB) - hist->GetXaxis()->GetBinLowEdge(iB) : 1.;

    hist->SetBinContent( iB, binCon / binWid);
    hist->SetBinError( iB, binErr / binWid );
  }

  return hist;
}



double EFTFitter::getDPoissonVariation(const double &mean) const
{
  if (!(mean > 0.)) return 0.;
  std::poisson_distribution<int> pois(mean);
  return mean - pois(FitUtil::rng);
}



/***
 * FitUtil source
 ***/



std::vector<double> FitUtil::extractBin(TH1 *hist)
{
  const int nBin = hist->GetNbinsX();
  std::vector<double> v_bin = { hist->GetXaxis()->GetBinLowEdge(1) };
  for (int iB = 1; iB <= nBin; ++iB)
    v_bin.push_back( hist->GetXaxis()->GetBinUpEdge( iB ) );

  return v_bin;
}



std::vector<std::array<double, 2>> FitUtil::extractContentError(TH1 *hist)
{
  const int nBin = hist->GetNbinsX();
  std::vector<std::array<double, 2>> v_conErr;
  for (int iB = 1; iB <= nBin; ++iB)
    v_conErr.push_back( {hist->GetBinContent( iB ), hist->GetBinError( iB )} );

  return v_conErr;
}
