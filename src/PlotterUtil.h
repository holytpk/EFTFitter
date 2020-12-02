// -*- C++ -*-
// an extension of PlotUtil, but now it directly makes the canvas and such

#ifndef PLOTTERUTIL_H
#define PLOTTERUTIL_H

#include "PlotUtil.h"
#include "TPaletteAxis.h"

// a simple wrapper - the plot itself, the legend text, legend type and draw options
template <typename P>
struct Plot {
  std::unique_ptr<P> plot;
  std::string legend_txt, legend_opt, draw_opt;
};

// do not put extension in filename
void standard_plot(std::vector<Plot<TH1>> &v_hist, const std::vector<Plot<TGraph>> &v_graph,
                   const std::string &filename, const bool add_uoflow, const bool normalize_shape, const bool logX, const bool logY, 
                   const int leg_ncolumn, const int leg_fill_color, const int leg_border_size, 
                   const int leg_font, const double leg_txt_size, const std::string &leg_header, 
                   const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                   const double axy_min, const double axy_max, const int axy_ndiv, 
                   const std::string &axy_txt, const double axy_txt_size, const double axy_offset, const double axy_label,
                   const double axx_min, const double axx_max, const int axx_ndiv, 
                   const std::string &axx_txt, const double axx_txt_size, const double axx_offset, const double axx_label,
                   const double can_top_margin, const double can_bottom_margin, 
                   const double can_left_margin, const double can_right_margin,
                   const std::string &plotformat = ".pdf") {
  // add underflow to bin 0, and overflow to bin N, if requested
  if (add_uoflow) {
    for (auto &hist: v_hist)
      add_uoflow_bin(hist.plot.get());
  }

  // check validity of histograms
  // if not, then drop it from the list
  // validity is defined simply as having at least one non-uoflow bin with either nonzero content or error
  for (auto &hist: v_hist) {
    const int nbin = hist.plot->GetNbinsX();
    bool onenonzero = false, allnonnan = true;
    for (int ib = 1; ib <= nbin; ++ib) {
      onenonzero = onenonzero or hist.plot->GetBinContent(ib) != 0. or hist.plot->GetBinError(ib) != 0.;
      allnonnan = allnonnan and !std::isnan(hist.plot->GetBinContent(ib)) and !std::isnan(hist.plot->GetBinError(ib));
    }

    if (!onenonzero or !allnonnan) {
      auto ptr = hist.plot.release();
      hist.legend_txt = ""; 
      hist.legend_opt = "";
      delete ptr;
    }
  }

  // normalize histograms if requested
  if (normalize_shape) {
    for (auto &hist: v_hist) {
      if (hist.plot != nullptr and hist.plot->Integral() != 0.)
        hist.plot->Scale(1. / hist.plot->Integral());
    }
  }

  // and axis setting
  for (auto &hist: v_hist) {
    if (hist.plot != nullptr)
      axisPlot(hist.plot.get(), 
               axy_min, axy_max, axy_ndiv, axy_txt, axy_txt_size, axy_offset, axy_label, 
               axx_min, axx_max, axx_ndiv, axx_txt, axx_txt_size, axx_offset, axx_label);
  }
  for (auto &graph: v_graph) {
    axisPlot(graph.plot.get(), 
             axy_min, axy_max, axy_ndiv, axy_txt, axy_txt_size, axy_offset, axy_label, 
             axx_min, axx_max, axx_ndiv, axx_txt, axx_txt_size, axx_offset, axx_label);
  }

  auto leg = std::make_unique<TLegend>();
  if (leg_x2 > leg_x1 and leg_y2 > leg_y1) {
    for (auto &hist: v_hist) {
      if (hist.plot != nullptr and hist.legend_txt != "" and hist.legend_opt != "")
        leg->AddEntry(hist.plot.get(), hist.legend_txt.c_str(), hist.legend_opt.c_str());
    }
    for (auto &graph: v_graph) {
      if (graph.legend_txt != "" and graph.legend_opt != "")
        leg->AddEntry(graph.plot.get(), graph.legend_txt.c_str(), graph.legend_opt.c_str());
    }
  }

  // making the canvas
  setH1Style();

  auto can = std::make_unique<TCanvas>("canvas", "canvas", 200, 10, 1000, 1000);
  can->SetFrameFillColor(0); // transparent plot
  can->SetFrameFillStyle(0); // transparent plot
  can->SetFrameBorderMode(0); // transparent plot
  can->SetFillColor(0); // transparent plot
  can->SetFillStyle(0); // transparent plot
  can->SetTopMargin(can_top_margin);
  can->SetBottomMargin(can_bottom_margin);
  can->SetLeftMargin(can_left_margin);
  can->SetRightMargin(can_right_margin);

  can->cd();
  if (logX)
    can->SetLogx();
  if (logY)
    can->SetLogy();

  if (leg_x2 > leg_x1 and leg_y2 > leg_y1) {
    styleLegend(leg.get(), leg_ncolumn, leg_fill_color, leg_border_size, 
                leg_font, leg_txt_size, leg_header);
    putLegend(leg.get(), leg_x1, leg_x2, leg_y1, leg_y2);
  }

  if (!v_hist.empty() and v_hist.front().plot != nullptr)
    v_hist.front().plot->Draw( v_hist.front().draw_opt.c_str() );
  else if (!v_graph.empty())
    v_graph.front().plot->Draw( ("a " + v_graph.front().draw_opt).c_str() );

  if (leg_x2 > leg_x1 and leg_y2 > leg_y1)
    leg->Draw();

  for (auto &hist: v_hist)
    if (hist.plot != nullptr)
      hist.plot->Draw( (hist.draw_opt + " same").c_str() );
  for (auto &graph: v_graph) 
    graph.plot->Draw( graph.draw_opt.c_str() );

  can->RedrawAxis();
  can->SaveAs( (filename + plotformat).c_str() );

  auto file = std::make_unique<TFile>((filename + ".root").c_str(), "recreate");
  file->cd();
  for (auto &hist : v_hist)
    if (hist.plot != nullptr)
      hist.plot->Write();
  for (auto &graph : v_graph)
    graph.plot->Write();
}



// do not put extension in filename
void standard_ratio(std::vector<Plot<TH1>> &v_hist,
                    const std::string &filename, const bool add_uoflow, const bool normalize_shape, const bool binomial_error,
                    const bool logX, const bool logY, 
                    const int leg_ncolumn, const int leg_fill_color, const int leg_border_size, 
                    const int leg_font, const double leg_txt_size, const std::string &leg_header, 
                    const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                    const double axy_min, const double axy_max, const int axy_ndiv,
                    const std::string &axy_txt, const double axy_txt_size, const double axy_offset, const double axy_label,
                    const double axr_min, const double axr_max, const int axr_ndiv,
                    const std::string &axr_txt, const double axr_txt_size, const double axr_offset, const double axr_label,
                    const double axx_min, const double axx_max, const int axx_ndiv,
                    const std::string &axx_txt, const double axx_txt_size, const double axx_offset, const double axx_label,
                    const double pad1_top_end, const double pad1_bottom_end, 
                    const double pad1_top_margin, const double pad1_bottom_margin,
                    const double pad2_top_end, const double pad2_bottom_end, 
                    const double pad2_top_margin, const double pad2_bottom_margin, 
                    const double can_left_margin, const double can_right_margin,
                    const std::string &plotformat = ".pdf") {
  // add underflow to bin 0, and overflow to bin N, if requested
  if (add_uoflow) {
    for (auto &hist: v_hist)
      add_uoflow_bin(hist.plot.get());
  }

  // check validity of histograms
  // if not, then drop it from the list
  // validity is defined simply as having at least one non-uoflow bin with either nonzero content or error
  for (auto &hist: v_hist) {
    const int nbin = hist.plot->GetNbinsX();
    bool onenonzero = false, allnonnan = true;
    for (int ib = 1; ib <= nbin; ++ib) {
      onenonzero = onenonzero or hist.plot->GetBinContent(ib) != 0. or hist.plot->GetBinError(ib) != 0.;
      allnonnan = allnonnan and !std::isnan(hist.plot->GetBinContent(ib)) and !std::isnan(hist.plot->GetBinError(ib));
    }

    if (!onenonzero or !allnonnan) {
      auto ptr = hist.plot.release();
      hist.legend_txt = ""; 
      hist.legend_opt = "";
      delete ptr;
    }
  }

  // normalize histograms if requested
  if (normalize_shape) {
    for (auto &hist: v_hist)
      if (hist.plot != nullptr and hist.plot->Integral() != 0.)
        hist.plot->Scale(1. / hist.plot->Integral());
  }

  // and axis setting
  for (auto &hist: v_hist) {
    if (hist.plot != nullptr)
      axisPlot(hist.plot.get(), 
               axy_min, axy_max, axy_ndiv, axy_txt, axy_txt_size, axy_offset, axy_label, 
               axx_min, axx_max, axx_ndiv, axx_txt, axx_txt_size, axx_offset, axx_label);
  }

  // make the ratio plots
  std::vector<Plot<TH1>> v_rhist(v_hist.size());
  for (uint iR = 0; iR < v_hist.size(); ++iR) {
    if (v_hist[iR].plot == nullptr)
      continue;

    const std::string name = v_hist[iR].plot->GetName();
    v_rhist[iR].draw_opt = v_hist[iR].draw_opt;
    v_rhist[iR].plot = std::unique_ptr<TH1D>(dynamic_cast<TH1D *>( v_hist[iR].plot->Clone( (name + "_ratio").c_str() )));
    v_rhist[iR].plot->Reset();

    if (binomial_error)
      v_rhist[iR].plot->Divide(v_hist[iR].plot.get(), v_hist[0].plot.get(), 1., 1., "B");
    else
      v_rhist[iR].plot->Divide(v_hist[iR].plot.get(), v_hist[0].plot.get(), 1., 1.);
  }

  for (auto &rhist: v_rhist) {
    if (rhist.plot != nullptr)
      axisPlot(rhist.plot.get(), 
               axr_min, axr_max, axr_ndiv, axr_txt, axr_txt_size, axr_offset, axr_label, 
               axx_min, axx_max, axx_ndiv, axx_txt, axx_txt_size, axx_offset, axx_label);
  }

  auto leg = std::make_unique<TLegend>();
  if (leg_x2 > leg_x1 and leg_y2 > leg_y1) {
    for (auto &hist: v_hist) {
      if (hist.plot != nullptr and hist.legend_txt != "" and hist.legend_opt != "")
        leg->AddEntry(hist.plot.get(), hist.legend_txt.c_str(), hist.legend_opt.c_str());
    }
  }

  // making the canvas
  setH1Style();

  auto can = std::make_unique<TCanvas>("canvas", "canvas", 200, 10, 1000, 1000);
  can->SetFrameFillColor(0); // transparent plot
  can->SetFrameFillStyle(0); // transparent plot
  can->SetFrameBorderMode(0); // transparent plot
  can->SetFillColor(0); // transparent plot
  can->SetFillStyle(0); // transparent plot
  can->cd();

  auto pad1 = std::make_unique<TPad>("pad1", "", 0., pad1_bottom_end, 1., pad1_top_end);
  pad1->SetFrameFillColor(0); // transparent plot
  pad1->SetFrameFillStyle(0); // transparent plot
  pad1->SetFrameBorderMode(0); // transparent plot
  pad1->SetFillColor(0); // transparent plot
  pad1->SetFillStyle(0); // transparent plot
  pad1->SetTopMargin(pad1_top_margin);
  pad1->SetBottomMargin(pad1_bottom_margin);
  pad1->SetLeftMargin(can_left_margin);
  pad1->SetRightMargin(can_right_margin);
  pad1->Draw();

  pad1->cd();
  if (logX)
    pad1->SetLogx();
  if (logY)
    pad1->SetLogy();

  if (leg_x2 > leg_x1 and leg_y2 > leg_y1) {
    styleLegend(leg.get(), leg_ncolumn, leg_fill_color, leg_border_size, 
                leg_font, leg_txt_size, leg_header);
    putLegend(leg.get(), leg_x1, leg_x2, leg_y1, leg_y2);
  }

  if (!v_hist.empty() and v_hist.front().plot != nullptr)
    v_hist.front().plot->Draw( v_hist.front().draw_opt.c_str() );

  if (leg_x2 > leg_x1 and leg_y2 > leg_y1)
    leg->Draw();

  for (auto &hist: v_hist) 
    if (hist.plot != nullptr)
      hist.plot->Draw( (hist.draw_opt + " same").c_str() );

  pad1->RedrawAxis();
  can->cd();

  auto pad2 = std::make_unique<TPad>("pad2", "", 0., pad2_bottom_end, 1., pad2_top_end);
  pad2->SetFrameFillColor(0); // transparent plot
  pad2->SetFrameFillStyle(0); // transparent plot
  pad2->SetFrameBorderMode(0); // transparent plot
  pad2->SetFillColor(0); // transparent plot
  pad2->SetFillStyle(0); // transparent plot
  pad2->SetTopMargin(pad2_top_margin);
  pad2->SetBottomMargin(pad2_bottom_margin);
  pad2->SetLeftMargin(can_left_margin);
  pad2->SetRightMargin(can_right_margin);
  pad2->Draw();

  pad2->cd();

  if (!v_rhist.empty() and v_rhist.front().plot != nullptr)
    v_rhist.front().plot->Draw( v_rhist.front().draw_opt.c_str() );

  for (auto &rhist: v_rhist)
    if (rhist.plot != nullptr) 
      rhist.plot->Draw( (rhist.draw_opt + " same").c_str() );

  pad2->RedrawAxis();
  can->cd();

  can->SaveAs( (filename + plotformat).c_str() );

  auto file = std::make_unique<TFile>((filename + ".root").c_str(), "recreate");
  file->cd();
  for (auto &hist : v_hist)
    hist.plot->Write();
  for (auto &rhist : v_rhist)
    rhist.plot->Write();
}



// do not put extension in filename
// FIXME still to be finalized; many args are dummy
void standard_colormap(const Plot<TH2> &hist,
                       const std::string &filename, const bool add_uoflow, const bool normalize_shape, 
                       const bool logX, const bool logY, const bool logZ, const int color_palette,
                       const double axz_min, const double axz_max, const int axz_ndiv,
                       const std::string &axz_txt, const double axz_txt_size, const double axz_offset, const double axz_label,
                       const double axz_begin, const double axz_width,
                       const double axy_min, const double axy_max, const int axy_ndiv,
                       const std::string &axy_txt, const double axy_txt_size, const double axy_offset, const double axy_label,
                       const double axx_min, const double axx_max, const int axx_ndiv,
                       const std::string &axx_txt, const double axx_txt_size, const double axx_offset, const double axx_label,
                       const double can_top_margin, const double can_bottom_margin, 
                       const double can_left_margin, const double can_right_margin,
                       const std::string &plotformat = ".pdf") {
  // will be implemented in the fwk plotter
  if (add_uoflow) {
    //add_uoflow_bin(hist.plot.get());
  }

  // normalize histograms if requested
  if (normalize_shape) {
    if (hist.plot->Integral() != 0.)
      hist.plot->Scale(1. / hist.plot->Integral());
  }

  // and axis setting
  axisPlot(hist.plot.get(), 
           axy_min, axy_max, axy_ndiv, axy_txt, axy_txt_size, axy_offset, axy_label, 
           axx_min, axx_max, axx_ndiv, axx_txt, axx_txt_size, axx_offset, axx_label);
  hist.plot->SetMinimum(axz_min);
  hist.plot->SetMaximum(axz_max);
  hist.plot->GetZaxis()->SetTitle(axz_txt.c_str());
  hist.plot->GetZaxis()->SetNdivisions(axz_ndiv);
  hist.plot->GetZaxis()->SetTitleSize(axz_txt_size);
  hist.plot->GetZaxis()->SetTitleOffset(axz_offset);
  hist.plot->SetContour(100);

  // making the canvas
  setH2Style();
  gStyle->SetPalette( color_palette );

  auto can = std::make_unique<TCanvas>("canvas", "canvas", 200, 10, 2000, 2000);
  can->SetFrameFillColor(0); // transparent plot
  can->SetFrameFillStyle(0); // transparent plot
  can->SetFrameBorderMode(0); // transparent plot
  can->SetFillColor(0); // transparent plot
  can->SetFillStyle(0); // transparent plot
  can->SetTopMargin(can_top_margin);
  can->SetBottomMargin(can_bottom_margin);
  can->SetLeftMargin(can_left_margin);
  can->SetRightMargin(can_right_margin);

  can->cd();
  if (logX)
    can->SetLogx();
  if (logY)
    can->SetLogy();
  if (logZ)
    can->SetLogz();

  hist.plot->Draw("colz");
  gPad->Update();
  can->SaveAs( (filename + plotformat).c_str() );

  // for tuning z axis size
  //auto palette = std::unique_ptr<TPaletteAxis>(dynamic_cast<TPaletteAxis *>( hist.plot->GetListOfFunctions()->FindObject("palette") ));
  TPaletteAxis *palette = (TPaletteAxis*) hist.plot->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(axz_begin);
  palette->SetX2NDC(axz_begin + axz_width);
  gPad->Modified();
  gPad->Update();

  can->RedrawAxis();
  can->SaveAs( (filename + plotformat).c_str() );

  auto file = std::make_unique<TFile>((filename + ".root").c_str(), "recreate");
  file->cd();
  hist.plot->Write();
}



void legend_canvas(const std::vector<Plot<TH1>> &v_hist, const std::vector<Plot<TGraph>> &v_graph, 
                   const std::string &filename,
                   const int leg_ncolumn, const int leg_fill_color, const int leg_border_size, 
                   const int leg_font, const double leg_txt_size, const std::string &leg_header, 
                   const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                   const std::string &plotformat = ".pdf") {
  if (!(leg_x2 > leg_x1 and leg_y2 > leg_y1))
    return;

  auto leg = std::make_unique<TLegend>();
  for (auto &hist: v_hist) {
    if (hist.legend_txt != "" and hist.legend_opt != "")
      leg->AddEntry(hist.plot.get(), hist.legend_txt.c_str(), hist.legend_opt.c_str());
  }
  for (auto &graph: v_graph) {
    if (graph.legend_txt != "" and graph.legend_opt != "")
      leg->AddEntry(graph.plot.get(), graph.legend_txt.c_str(), graph.legend_opt.c_str());
  }

  // making the canvas
  setH1Style();
  auto can = std::make_unique<TCanvas>("canvas", "canvas", 200, 10, 1000, 1000);
  can->SetFrameFillColor(0); // transparent plot
  can->SetFrameFillStyle(0); // transparent plot
  can->SetFrameBorderMode(0); // transparent plot
  can->SetFillColor(0); // transparent plot
  can->SetFillStyle(0); // transparent plot
  can->SetTopMargin(0.);
  can->SetBottomMargin(0.);
  can->SetLeftMargin(0.);
  can->SetRightMargin(0.);

  can->cd();
  styleLegend(leg.get(), leg_ncolumn, leg_fill_color, leg_border_size, 
              leg_font, leg_txt_size, leg_header);
  putLegend(leg.get(), leg_x1, leg_x2, leg_y1, leg_y2);

  // draw and save
  leg->Draw();
  can->SaveAs( (filename + plotformat).c_str() );
}



// written without safety to be faster
std::array<double, 2> cumulative_sum(TH1* hist, int iBin) {
  std::array<double, 2> a_sum = {0., 0.};
  for (int iB = iBin; iB > 0; iB--) {
    a_sum[0] += hist->GetBinContent(iB);
    a_sum[1] += hist->GetBinError(iB) * hist->GetBinError(iB);
  }

  a_sum[1] = std::sqrt(a_sum[1]);
  return a_sum;
}



// last arg: true for efficiency, false for rejection (1 - eff)
Plot<TGraph> efficiency_profile(const Plot<TH1> &hist, const bool efficiency = true, const std::string &draw_opt = "lx") {
  if (hist.plot == nullptr) {
    std::cout << "ERROR: efficiency_profile(): histogram is null. Returning a null plot!!" << std::endl;
    return Plot<TGraph>();
  }

  // replace returns a bool indicating whether the replacement succeeds
  // i.e. the string has a from substring to be replaced
  std::string name = hist.plot->GetName(), tag = (efficiency) ? "_eff" : "_rej";
  if (!replace(name, "_hist", tag))
    name = name + tag;

  const int nBin = hist.plot->GetNbinsX();
  const double integral = hist.plot->Integral();
  std::vector<double> v_axx(nBin + 2), v_axe(nBin + 2), v_eff(nBin + 2), v_err(nBin + 2);
  v_axx[0] = hist.plot->GetBinLowEdge(1);
  v_axe[0] = hist.plot->GetBinWidth(nBin) / 2.;
  v_eff[0] = (efficiency) ? 0. : 1.;
  v_err[0] = 0.;

  for (int iB = 1; iB <= nBin; ++iB) {
    v_axx[iB] = hist.plot->GetBinLowEdge(iB) + (hist.plot->GetBinWidth(iB) / 2.);
    v_axe[iB] = hist.plot->GetBinWidth(iB) / 2.;

    auto sum = cumulative_sum(hist.plot.get(), iB);
    v_eff[iB] = (efficiency) ? sum[0] / integral : 1. - (sum[0] / integral);
    v_err[iB] = sum[1] / integral;
  }
  v_axx[nBin + 1] = hist.plot->GetBinLowEdge(nBin) + hist.plot->GetBinWidth(nBin);
  v_axe[nBin + 1] = hist.plot->GetBinWidth(nBin) / 2.;
  v_eff[nBin + 1] = v_eff[nBin];
  v_err[nBin + 1] = v_err[nBin];

  const int nPnt = (count_substring(draw_opt, "3") or count_substring(draw_opt, "4") or count_substring(draw_opt, "x")) 
    ? v_axx.size() : v_axx.size() - 1;
  Plot<TGraph> pl;
  pl.plot = std::make_unique<TGraphErrors>(nPnt, v_axx.data(), v_eff.data(), v_axe.data(), v_err.data());
  pl.plot->SetName(name.c_str());
  pl.draw_opt = draw_opt;
  pl.legend_opt = "l";
  pl.legend_txt = hist.legend_txt;
  stylePlot(pl.plot.get(), hist.plot->GetLineColor(), 1., 1001, 1, 1.5, 1, 3);

  return std::move(pl);
}



// for the moment explicitly ignore the uncertainties in sig and bkg 
Plot<TGraph> discrimination_profile(const std::vector<Plot<TH1>> &v_hist, uint isig = 0u, uint ibkg = 1u,
                                    double start_sig = 0., double start_bkg = 1., bool interpolate = false,
                                    int color = kBlack, bool print_integral = false) {
  if (v_hist.size() < 2u or isig > v_hist.size() - 1u or ibkg > v_hist.size() - 1u) {
    std::cout << "ERROR: discrimination_profile(): unsuitable argument set. Returning a null plot!!" << std::endl;
    return Plot<TGraph>();
  }

  if (start_sig < 0. or start_sig > 1.)
    start_sig = 0.;
  if (start_bkg < 0. or start_bkg > 1.)
    start_bkg = 1.;

  auto eff_sig = efficiency_profile(v_hist[isig], true, "lx");
  auto rej_bkg = efficiency_profile(v_hist[ibkg], false, "lx");

  if (eff_sig.plot->GetN() != rej_bkg.plot->GetN()) {
    std::cout << "ERROR: discrimination_profile(): signal and background graphs do not match. Returning a null plot!!" << std::endl;
    return Plot<TGraph>();
  }

  std::string name = eff_sig.plot->GetName();
  replace(name, "_eff", "_roc");

  Plot<TGraph> pl;
  int npnt = eff_sig.plot->GetN(), iinterp = 1;
  std::vector<double> sigpnt, bkgpnt;
  double area = 0., sig1 = 0., bkg1 = 0., sx1 = 0., sy1 = 0., bx1 = 0., by1 = 0.;

  if (start_sig == 0. and start_bkg == 1.) {
    pl.plot = std::make_unique<TGraph>(npnt, eff_sig.plot->GetY(), rej_bkg.plot->GetY());
    area = 0.5;
  }
  else {
    for (int ipnt = 0; ipnt < npnt; ++ipnt) {
      if (eff_sig.plot->GetPoint(ipnt, sx1, sy1) == -1 or rej_bkg.plot->GetPoint(ipnt, bx1, by1) == -1)
        continue;

      if (sy1 >= start_sig and by1 <= start_bkg) {
        if (interpolate and iinterp and sy1 != start_sig and by1 != start_bkg) {
          double sx0 = 0., sy0 = 0., bx0 = 0., by0 = 0.;

          if (eff_sig.plot->GetPoint(ipnt - 1, sx0, sy0) == -1 or rej_bkg.plot->GetPoint(ipnt - 1, bx0, by0) == -1)
            continue;

          double slope = (by1 - by0) / (sy1 - sy0);
          double offs = (start_sig > sy0 and start_sig < sy1) ? (start_sig - sy0) / (sy1 - sy0) : -9999.;
          double offb = (start_bkg < by0 and start_bkg > by1) ? (start_bkg - by0) / (by1 - by0) : -9999.;
          double syi = 0., byi = 0.;

          if (offs > offb) {
            syi = start_sig;
            byi = (slope * (syi - sy0)) + by0;
          }
          else {
            byi = start_bkg;
            syi = ((byi - by0) / slope) + sy0;
          }

          sigpnt.emplace_back(syi);
          bkgpnt.emplace_back(byi);

          iinterp = 0;
        }

        sigpnt.emplace_back(sy1);
        bkgpnt.emplace_back(by1);

        if (sig1 == 0.) {
          sig1 = sigpnt.front();
          bkg1 = bkgpnt.front();
        }
      }
    }

    pl.plot = std::make_unique<TGraph>(sigpnt.size(), sigpnt.data(), bkgpnt.data());
    area = (1. - sig1) * bkg1 / 2.;
  }

  pl.plot->SetName(name.c_str());
  pl.draw_opt = "l";
  pl.legend_opt = "l";
  pl.legend_txt = "";
  stylePlot(pl.plot.get(), color, 1., 0, 24, 1.5, 1, 3);

  // the + area is because of how ROOT::TGraph::Integral() works
  // what it really calculates is the area formed by a polygon formed by the set of points
  // where the polygon is closed by connecting the last and first points with a straight line
  // which for our ROCs is the line from (0, 1) to (1, 1) i.e. the worst possible ROC 
  if (print_integral)
    std::cout << "Integral of ROC curve " << name << ": " << pl.plot->Integral() + area << std::endl;

  return std::move(pl);
}

#endif
