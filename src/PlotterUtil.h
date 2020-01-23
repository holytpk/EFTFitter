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
                   const std::string &filename, const bool normalize_shape, 
                   const int leg_ncolumn, const int leg_fill_color, const int leg_border_size, 
                   const int leg_font, const double leg_txt_size, const std::string &leg_header, 
                   const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2,
                   const double axy_min, const double axy_max,
                   const std::string &axy_txt, const double axy_txt_size, const double axy_offset, const double axy_label,
                   const double axx_min, const double axx_max,
                   const std::string &axx_txt, const double axx_txt_size, const double axx_offset, const double axx_label,
                   const double can_top_margin, const double can_bottom_margin, 
                   const double can_left_margin, const double can_right_margin) {
  // normalize histograms if requested
  if (normalize_shape) {
    for (auto &hist: v_hist)
      hist.plot->Scale(1. / hist.plot->Integral());
  }

  // and axis setting
  if (!v_hist.empty()) {
    for (auto &hist: v_hist) {
      axisPlot(hist.plot.get(), 
               axy_min, axy_max, axy_txt, axy_txt_size, axy_offset, axy_label, 
               axx_min, axx_max, axx_txt, axx_txt_size, axx_offset, axx_label);
    }
  }
  if (!v_graph.empty()) {
    for (auto &graph: v_hist) {
      axisPlot(graph.plot.get(), 
               axy_min, axy_max, axy_txt, axy_txt_size, axy_offset, axy_label, 
               axx_min, axx_max, axx_txt, axx_txt_size, axx_offset, axx_label);
    }
  }

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
  can->SetTopMargin(can_top_margin);
  can->SetBottomMargin(can_bottom_margin);
  can->SetLeftMargin(can_left_margin);
  can->SetRightMargin(can_right_margin);

  can->cd();
  styleLegend(leg.get(), leg_ncolumn, leg_fill_color, leg_border_size, 
              leg_font, leg_txt_size, leg_header);
  putLegend(leg.get(), leg_x1, leg_x2, leg_y1, leg_y2); // bottom center

  if (!v_hist.empty())
    v_hist.front().plot->Draw( v_hist.front().draw_opt.c_str() );
  else if (!v_graph.empty())
    v_graph.front().plot->Draw( v_graph.front().draw_opt.c_str() );
  leg->Draw();

  for (auto &hist: v_hist) 
    hist.plot->Draw( (hist.draw_opt + " same").c_str() );
  for (auto &graph: v_graph) 
    graph.plot->Draw( graph.draw_opt.c_str() );

  can->RedrawAxis();
  can->SaveAs( (filename + ".pdf").c_str() );
  can->SaveAs( (filename + ".root").c_str() );
}

#endif
