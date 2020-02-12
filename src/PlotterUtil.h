// -*- C++ -*-
// an extension of PlotUtil, but now it directly makes the canvas and such

#ifndef PLOTTERUTIL_H
#define PLOTTERUTIL_H

#include "PlotUtil.h"

// a simple wrapper - the plot itself, the legend text, legend type and draw options
template <typename P>
struct Plot {
  std::unique_ptr<P> plot;
  std::string legend_txt, legend_opt, draw_opt;
};

// do not put extension in filename
void standard_plot(const std::vector<Plot<TH1>> &v_hist, const std::vector<Plot<TGraph>> &v_graph,
                   const std::string &filename, const bool add_uoflow, const bool normalize_shape, const bool logX, const bool logY, 
                   const int leg_ncolumn, const int leg_fill_color, const int leg_border_size, 
                   const int leg_font, const double leg_txt_size, const std::string &leg_header, 
                   const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                   const double axy_min, const double axy_max,
                   const std::string &axy_txt, const double axy_txt_size, const double axy_offset, const double axy_label,
                   const double axx_min, const double axx_max,
                   const std::string &axx_txt, const double axx_txt_size, const double axx_offset, const double axx_label,
                   const double can_top_margin, const double can_bottom_margin, 
                   const double can_left_margin, const double can_right_margin,
                   const std::string &plotformat = ".pdf") {
  if (add_uoflow) {
    for (auto &hist: v_hist)
      add_uoflow_bin(hist.plot.get());
  }

  // normalize histograms if requested
  if (normalize_shape) {
    for (auto &hist: v_hist)
      hist.plot->Scale(1. / hist.plot->Integral());
  }

  // and axis setting
  for (auto &hist: v_hist) {
    axisPlot(hist.plot.get(), 
             axy_min, axy_max, axy_txt, axy_txt_size, axy_offset, axy_label, 
             axx_min, axx_max, axx_txt, axx_txt_size, axx_offset, axx_label);
  }
  for (auto &graph: v_graph) {
    axisPlot(graph.plot.get(), 
             axy_min, axy_max, axy_txt, axy_txt_size, axy_offset, axy_label, 
             axx_min, axx_max, axx_txt, axx_txt_size, axx_offset, axx_label);
  }

  auto leg = std::make_unique<TLegend>();
  if (leg_x2 > leg_x1 and leg_y2 > leg_y1) {
    for (auto &hist: v_hist) {
      if (hist.legend_txt != "" and hist.legend_opt != "")
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

  if (!v_hist.empty())
    v_hist.front().plot->Draw( v_hist.front().draw_opt.c_str() );
  else if (!v_graph.empty())
    v_graph.front().plot->Draw( ("a " + v_graph.front().draw_opt).c_str() );

  if (leg_x2 > leg_x1 and leg_y2 > leg_y1)
    leg->Draw();

  for (auto &hist: v_hist) 
    hist.plot->Draw( (hist.draw_opt + " same").c_str() );
  for (auto &graph: v_graph) 
    graph.plot->Draw( graph.draw_opt.c_str() );

  can->RedrawAxis();
  can->SaveAs( (filename + plotformat).c_str() );
  can->SaveAs( (filename + ".root").c_str() );
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
Plot<TGraph> discrimination_profile(const std::vector<Plot<TH1>> &v_hist, const int sig = 0, const int bkg = 1,
                                    const int color = kBlack, const std::string &legend_txt = "", const bool print_integral = false) {
  if (v_hist.size() < 2 or sig < 0 or bkg < 0 or sig > v_hist.size() - 1 or bkg > v_hist.size() - 1) {
    std::cout << "ERROR: discrimination_profile(): unsuitable argument set. Returning a null plot!!" << std::endl;
    return Plot<TGraph>();
  }

  auto eff_sig = efficiency_profile(v_hist[sig], true, "lx");
  auto rej_bkg = efficiency_profile(v_hist[bkg], false, "lx");

  if (eff_sig.plot->GetN() != rej_bkg.plot->GetN()) {
    std::cout << "ERROR: discrimination_profile(): signal and background graphs do not match. Returning a null plot!!" << std::endl;
    return Plot<TGraph>();
  }

  std::string name = eff_sig.plot->GetName();
  replace(name, "_eff", "_roc");

  Plot<TGraph> pl;
  pl.plot = std::make_unique<TGraph>(eff_sig.plot->GetN(), eff_sig.plot->GetY(), rej_bkg.plot->GetY());
  pl.plot->SetName(name.c_str());
  pl.draw_opt = "l";
  pl.legend_opt = "l";
  pl.legend_txt = legend_txt;
  stylePlot(pl.plot.get(), color, 1., 0, 1, 1.5, 1, 3);

  // the + 0.5 is because of how ROOT::TGraph::Integral() works
  // what it really calculates is the area formed by a polygon formed by the set of points
  // where the polygon is closed by connecting the last and first points with a straight line
  // which for our ROCs is the line from (0, 1) to (1, 1) i.e. the worst possible ROC 
  if (print_integral) {
    std::cout << "Integral of ROC curve " << name << " from signal/background index of " << sig << "/" << bkg 
              << ": " << pl.plot->Integral() + 0.5 << std::endl;
  }

  return std::move(pl);
}

#endif
