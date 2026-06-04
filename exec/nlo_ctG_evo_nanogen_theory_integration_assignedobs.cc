// ./execMacro.sh nlo_ctG_evo_nanogen_theory_integration_assignedobs.cc

// -----------------------------------------------------------------------------
// Example: how the per-bin theory histogram is constructed for gen_c_kk
// at ctG = 2.0 using the inclusive SMEFT theory coefficients.
//
// Embedded theory table inputs (LO example):
//
//   sigma_SM             = 29.12
//   sigma_ctG_lin        =  9.32
//   sigma_ctG_quad       =  1.63
//
//   sigma_Ckk_SM         =  9.42
//   sigma_Ckk_ctG_lin    =  4.60
//   sigma_Ckk_ctG_quad   =  0.97
//
// Step 1: compute inclusive total cross section
//
//   den = sigma_SM
//       + ctG * sigma_ctG_lin
//       + ctG^2 * sigma_ctG_quad
//
// For ctG = 2:
//
//   den = 29.12 + 2*(9.32) + 4*(1.63)
//       = 54.28
//
// Step 2: compute inclusive spin-correlation numerator
//
//   num = sigma_Ckk_SM
//       + ctG * sigma_Ckk_ctG_lin
//       + ctG^2 * sigma_Ckk_ctG_quad
//
//   num = 9.42 + 2*(4.60) + 4*(0.97)
//       = 22.50
//
// Step 3: compute inclusive coefficient
//
//   Ckk(ctG=2) = num / den
//               = 22.50 / 54.28
//               = 0.4145
//
// Step 4: convert coefficient into forward-backward asymmetry
//
// For spin-correlation observables:
//
//   AFB = coeff / asymmetry_factor(obs)
//
// and for gen_c_kk:
//
//   asymmetry_factor("gen_c_kk") = -4
//
// therefore:
//
//   AFB = 0.4145 / (-4)
//       = -0.1036
//
// Step 5: redistribute the reference histogram halves
//
// The reference histogram already contains the differential shape:
//
//   bins 0,1,2 = backward half
//   bins 3,4,5 = forward half
//
// Define:
//
//   total = total integral of reference histogram
//   Bsum  = sum of bins 0,1,2
//   Fsum  = sum of bins 3,4,5
//
// Target backward/forward normalizations become:
//
//   target_B = 0.5 * total * (1 - AFB)
//            = 0.5518 * total
//
//   target_F = 0.5 * total * (1 + AFB)
//            = 0.4482 * total
//
// Each bin is then rescaled proportionally:
//
//   for bins 0,1,2:
//       out[i] = ref[i] * target_B / Bsum
//
//   for bins 3,4,5:
//       out[i] = ref[i] * target_F / Fsum
//
// Example: 2nd bin of gen_c_kk (C++ index i=1)
//
//   out[1]
//     = ref[1] * target_B / Bsum
//
// so the theory table only fixes the inclusive coefficient,
// while the detailed differential shape comes from the reference histogram.
// ----------------------------------------------------------------------------- 

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TSystem.h>
#include <TMarker.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TBox.h>
#include <TGraphAsymmErrors.h>
#include <TH2D.h>
#include <TColor.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TPaveText.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <cstdio>

static const int BINS_PER_OBS = 6;
static const int N_OBS_TOTAL = 38;
static const std::string THEORY_ORDER = "LO";  // compare LO MC with LO theory

std::vector<int> parse_obs(const std::string& s) {
  std::vector<int> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (!item.empty()) out.push_back(std::stoi(item));
  }
  return out;
}

std::vector<int> build_keep_indices(const std::vector<int>& obs, int drop_bin_idx) {
  std::vector<int> keep;
  for (int iobs : obs) {
    for (int ibin = 0; ibin < BINS_PER_OBS; ++ibin) {
      if (ibin != drop_bin_idx) keep.push_back(iobs * BINS_PER_OBS + ibin);
    }
  }
  return keep;
}

TH1* get_hist1(TFile* f, const std::string& key = "diff_cross_section") {
  TH1* h = dynamic_cast<TH1*>(f->Get(key.c_str()));
  if (h) return h;

  TIter next(f->GetListOfKeys());
  TObject* k = nullptr;
  while ((k = next())) {
    TObject* obj = f->Get(k->GetName());
    h = dynamic_cast<TH1*>(obj);
    if (h) return h;
  }
  return nullptr;
}

std::vector<double> load_values(const std::string& path) {
  TFile f(path.c_str(), "READ");
  if (f.IsZombie()) throw std::runtime_error("Cannot open " + path);

  TH1* h = get_hist1(&f);
  if (!h) throw std::runtime_error("No TH1 found in " + path);

  std::vector<double> v(h->GetNbinsX());
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    v[i - 1] = h->GetBinContent(i);
  }
  return v;
}

TMatrixD load_cov_tree_one(const std::string& path,
                           const std::string& tree_name,
                           const std::string& branch_hint = "CovMatrix_Norm_bin_width") {
  TFile f(path.c_str(), "READ");
  if (f.IsZombie()) throw std::runtime_error("Cannot open covariance file: " + path);

  TTree* t = dynamic_cast<TTree*>(f.Get(tree_name.c_str()));
  if (!t) {
    f.ls();
    throw std::runtime_error("Cannot find TTree " + tree_name + " in " + path);
  }

  TBranch* br = t->GetBranch(branch_hint.c_str());
  if (!br) {
    t->Print();
    throw std::runtime_error("Cannot find branch " + branch_hint + " in tree " + tree_name);
  }

  const Long64_t nentries = t->GetEntries();
  const int n = (int)nentries;

  double row[228];
  for (int j = 0; j < 228; ++j) row[j] = 0.0;

  t->SetBranchAddress(branch_hint.c_str(), row);

  TMatrixD cov(n, n);

  for (Long64_t i = 0; i < nentries; ++i) {
    t->GetEntry(i);
    for (int j = 0; j < n; ++j) {
      cov((int)i, j) = row[j];
    }
  }

  std::cout << "[OK] loaded covariance array branch "
            << tree_name << "/" << branch_hint
            << " as " << n << "x" << n << std::endl;

  return cov;
}

TMatrixD load_covariance(const std::string& stat_path, const std::string& syst_path) {
  TMatrixD cov_stat = load_cov_tree_one(stat_path, "stat_1D_variables");
  TMatrixD cov_syst = load_cov_tree_one(syst_path, "syst_1D_variables");

  if (cov_stat.GetNrows() != cov_syst.GetNrows()) {
    throw std::runtime_error("Stat/syst covariance dimensions do not match");
  }

  TMatrixD cov = cov_stat;
  cov += cov_syst;
  return cov;
}

std::vector<double> select_vec(const std::vector<double>& v, const std::vector<int>& keep) {
  std::vector<double> out;
  out.reserve(keep.size());
  for (int idx : keep) out.push_back(v.at(idx));
  return out;
}

TMatrixD select_cov(const TMatrixD& cov_full,
                    const std::vector<int>& keep228,
                    int drop_bin_idx) {
  std::vector<int> allobs;
  for (int i = 0; i < N_OBS_TOTAL; ++i) allobs.push_back(i);

  std::vector<int> full_keep = build_keep_indices(allobs, drop_bin_idx);

  std::map<int, int> remap;
  for (int i = 0; i < (int)full_keep.size(); ++i) remap[full_keep[i]] = i;

  std::vector<int> keep_cov;

  if (cov_full.GetNrows() == (int)full_keep.size()) {
    for (int idx : keep228) keep_cov.push_back(remap.at(idx));
  } else if (cov_full.GetNrows() == N_OBS_TOTAL * BINS_PER_OBS) {
    keep_cov = keep228;
  } else {
    std::stringstream ss;
    ss << "Unexpected covariance size: " << cov_full.GetNrows()
       << ", expected " << full_keep.size()
       << " or " << N_OBS_TOTAL * BINS_PER_OBS;
    throw std::runtime_error(ss.str());
  }

  int n = keep_cov.size();
  TMatrixD out(n, n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      out(i, j) = cov_full(keep_cov[i], keep_cov[j]);
    }
  }
  return out;
}


void save_matrix_inspection_plot(const TMatrixD& m,
                                 const std::string& outdir,
                                 const std::string& tag,
                                 const std::string& ztitle = "matrix value") {
  gStyle->SetOptStat(0);
  std::string od = outdir + "/covariance_inspection";
  gSystem->mkdir(od.c_str(), true);

  const int nr = m.GetNrows();
  const int nc = m.GetNcols();
  TH2D h(("h_" + tag).c_str(), "", nc, 0, nc, nr, 0, nr);
  double zmax_abs = 0.0;
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      const double v = m(i, j);
      h.SetBinContent(j + 1, nr - i, v);
      zmax_abs = std::max(zmax_abs, std::fabs(v));
    }
  }
  h.SetTitle("");
  h.GetXaxis()->SetTitle("matrix column");
  h.GetYaxis()->SetTitle("matrix row");
  h.GetZaxis()->SetTitle(ztitle.c_str());
  h.GetXaxis()->SetTitleSize(0.045);
  h.GetYaxis()->SetTitleSize(0.045);
  h.GetZaxis()->SetTitleSize(0.040);
  h.GetXaxis()->SetLabelSize(0.030);
  h.GetYaxis()->SetLabelSize(0.030);
  h.GetZaxis()->SetLabelSize(0.030);
  if (zmax_abs > 0.0) h.GetZaxis()->SetRangeUser(-zmax_abs, zmax_abs);

  TCanvas c(("c_" + tag).c_str(), ("c_" + tag).c_str(), 900, 760);
  
  gStyle->SetOptStat(0);c.SetLeftMargin(0.12);
  c.SetRightMargin(0.18);
  c.SetBottomMargin(0.11);
  c.SetTopMargin(0.08);
  h.SetStats(0);
  h.Draw("COLZ");

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextSize(0.032);
  lat.DrawLatex(0.12, 0.945, tag.c_str());

  c.SaveAs((od + "/" + tag + ".png").c_str());
  c.SaveAs((od + "/" + tag + ".pdf").c_str());

  std::ofstream csv(od + "/" + tag + ".csv");
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      if (j) csv << ",";
      csv << std::setprecision(12) << m(i, j);
    }
    csv << "\n";
  }
}

void dump_covariance_inspection_plots(const TMatrixD& cov_full,
                                      const std::string& outdir,
                                      int drop_bin_idx) {
  std::cout << "[INFO cov] writing covariance inspection plots to "
            << outdir << "/covariance_inspection" << std::endl;

  save_matrix_inspection_plot(cov_full, outdir,
                              "cov_00_raw_loaded_before_reduction",
                              "covariance");

  std::vector<int> allobs;
  for (int i = 0; i < N_OBS_TOTAL; ++i) allobs.push_back(i);
  std::vector<int> keep = build_keep_indices(allobs, drop_bin_idx);

  TMatrixD cov_reduced = select_cov(cov_full, keep, drop_bin_idx);
  save_matrix_inspection_plot(cov_reduced, outdir,
                              "cov_01_reduced_before_inversion",
                              "covariance");

  TDecompSVD svd(cov_reduced);
  TMatrixD cov_inv = svd.Invert();
  save_matrix_inspection_plot(cov_inv, outdir,
                              "cov_02_reduced_inverse_after_inversion",
                              "inverse covariance");
}

std::string replace_all(std::string s, const std::string& a, const std::string& b) {
  size_t pos = 0;
  while ((pos = s.find(a, pos)) != std::string::npos) {
    s.replace(pos, a.size(), b);
    pos += b.size();
  }
  return s;
}

std::string valstr(int v) {
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%+d", v);
  return std::string(buf);
}

std::string make_template_path(std::string pattern,
                               const std::string& wc,
                               int val) {
  pattern = replace_all(pattern, "{wc}", wc);
  pattern = replace_all(pattern, "{val}", valstr(val));
  return pattern;
}

std::vector<double> block_renorm(std::vector<double> pred,
                                 const std::vector<double>& ref,
                                 int bins_eff) {
  int nblock = pred.size() / bins_eff;

  for (int ib = 0; ib < nblock; ++ib) {
    double sp = 0.0;
    double sr = 0.0;

    for (int k = 0; k < bins_eff; ++k) {
      sp += pred[ib * bins_eff + k];
      sr += ref [ib * bins_eff + k];
    }

    if (std::abs(sp) > 0.0 && std::abs(sr) > 0.0) {
      double scale = sr / sp;
      for (int k = 0; k < bins_eff; ++k) {
        pred[ib * bins_eff + k] *= scale;
      }
    }
  }

  return pred;
}

struct FitResult1D {
  std::string name;
  std::string label;
  double best = 0.0;
  double chi2min = 0.0;
  double lo68 = 0.0, hi68 = 0.0;
  double lo95 = 0.0, hi95 = 0.0;
  int nbins = 0;
};

std::string pretty_label(const std::string& wc) {
  // ROOT/TLatex labels for the WC summary y-axis.
  // Keep these close to CMS/TOP publication notation.
  if (wc == "ctGRe") return "c_{tG}^{Re}";
  if (wc == "ctGIm") return "c_{tG}^{Im}";
  if (wc == "cQj11") return "c_{Qq}^{(1,1)}";
  if (wc == "cQj31") return "c_{Qq}^{(3,1)}";
  if (wc == "cQj18") return "c_{Qq}^{(1,8)}";
  if (wc == "cQj38") return "c_{Qq}^{(3,8)}";
  if (wc == "cQu1")  return "c_{Qu}^{(1)}";
  if (wc == "cQu8")  return "c_{Qu}^{(8)}";
  if (wc == "cQd1")  return "c_{Qd}^{(1)}";
  if (wc == "cQd8")  return "c_{Qd}^{(8)}";
  if (wc == "ctu1")  return "c_{tu}^{(1)}";
  if (wc == "ctu8")  return "c_{tu}^{(8)}";
  if (wc == "ctd1")  return "c_{td}^{(1)}";
  if (wc == "ctd8")  return "c_{td}^{(8)}";
  if (wc == "ctj1")  return "c_{tq}^{(1)}";
  if (wc == "ctj8")  return "c_{tq}^{(8)}";
  if (wc == "mut") return "#hat{#mu}_{t}";
  if (wc == "#hat{d}_{t}") return "#hat{d}_{t}";
  return wc;
}

std::vector<int> range_obs(int a, int b_inclusive) {
  std::vector<int> v;
  for (int i = a; i <= b_inclusive; ++i) v.push_back(i);
  return v;
}

std::pair<double,double> interval_crossing(const std::vector<double>& xs,
                                           const std::vector<double>& ys,
                                           double level,
                                           double best) {
  double lo = best, hi = best;

  for (int i = 1; i < (int)xs.size(); ++i) {
    if ((xs[i-1] <= best && best <= xs[i]) || (xs[i] <= best && best <= xs[i-1])) {
      break;
    }
  }

  for (int i = 1; i < (int)xs.size(); ++i) {
    if (xs[i] > best) break;
    if ((ys[i-1] - level) * (ys[i] - level) <= 0.0) {
      double t = (level - ys[i-1]) / (ys[i] - ys[i-1] + 1e-30);
      lo = xs[i-1] + t * (xs[i] - xs[i-1]);
    }
  }

  for (int i = 1; i < (int)xs.size(); ++i) {
    if (xs[i-1] < best) continue;
    if ((ys[i-1] - level) * (ys[i] - level) <= 0.0) {
      double t = (level - ys[i-1]) / (ys[i] - ys[i-1] + 1e-30);
      hi = xs[i-1] + t * (xs[i] - xs[i-1]);
      break;
    }
  }

  return {lo, hi};
}

void save_chi2_plot(const std::string& outpdf,
                    const std::string& outpng,
                    const std::vector<double>& xs,
                    const std::vector<double>& dchi,
                    const FitResult1D& r,
                    double xscale = 1.0,
                    const std::string& xlabel = "") {
  gStyle->SetOptStat(0);

  std::vector<double> xplot(xs.size());
  for (size_t i = 0; i < xs.size(); ++i) xplot[i] = xs[i] * xscale;

  double xmin = r.lo95 * xscale;
  double xmax = r.hi95 * xscale;
  double pad = 0.25 * std::max(1e-9, xmax - xmin);
  xmin -= pad;
  xmax += pad;

  TCanvas c("c", "c", 760, 660);
  
  gStyle->SetOptStat(0);c.SetLeftMargin(0.13);
  c.SetRightMargin(0.05);
  c.SetTopMargin(0.08);
  c.SetBottomMargin(0.13);

  TGraph gr(xs.size(), xplot.data(), dchi.data());
  gr.SetLineWidth(3);
  gr.SetLineColor(kBlack);
  gr.SetTitle("");
  gr.GetXaxis()->SetTitle(xlabel.empty() ? r.label.c_str() : xlabel.c_str());
  gr.GetYaxis()->SetTitle("#Delta#chi^{2}");
  gr.GetXaxis()->SetLimits(xmin, xmax);
  gr.GetYaxis()->SetRangeUser(0.0, 10.0);
  gr.GetXaxis()->SetTitleSize(0.055);
  gr.GetYaxis()->SetTitleSize(0.060);
  gr.GetXaxis()->SetLabelSize(0.045);
  gr.GetYaxis()->SetLabelSize(0.045);
  gr.Draw("AL");

  const int green = TColor::GetColor("#2ca25f");
  const int yellow = TColor::GetColor("#ffd92f");

  TBox b95(r.lo95*xscale, 0, r.hi95*xscale, 10);
  b95.SetFillColorAlpha(yellow, 0.85);
  b95.SetLineColor(0);
  b95.Draw("same");

  TBox b68(r.lo68*xscale, 0, r.hi68*xscale, 10);
  b68.SetFillColorAlpha(green, 0.90);
  b68.SetLineColor(0);
  b68.Draw("same");

  gr.Draw("L same");

  TLine l1(xmin, 1.0, xmax, 1.0); l1.SetLineStyle(2); l1.SetLineColor(kGray+2); l1.Draw();
  TLine l4(xmin, 4.0, xmax, 4.0); l4.SetLineStyle(3); l4.SetLineColor(kGray+2); l4.Draw();

  TLine l0(0.0, 0.0, 0.0, 10.0); l0.SetLineStyle(1); l0.SetLineColor(kBlack); l0.SetLineWidth(2); l0.Draw();
  TLine lb(r.best*xscale, 0.0, r.best*xscale, 10.0); lb.SetLineStyle(2); lb.SetLineColor(kBlack); lb.SetLineWidth(2); lb.Draw();

  TLegend leg(0.68, 0.64, 0.92, 0.82);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.030);
  leg.SetMargin(0.22);
  leg.SetEntrySeparation(0.15);
  TBox dummy95; dummy95.SetFillColor(yellow);
  TBox dummy68; dummy68.SetFillColor(green);
  TLine dummyNom; dummyNom.SetLineColor(kBlack); dummyNom.SetLineWidth(2);
  TLine dummyBest; dummyBest.SetLineColor(kBlack); dummyBest.SetLineStyle(2); dummyBest.SetLineWidth(2);
  leg.AddEntry(&dummy95, "95% CL", "f");
  leg.AddEntry(&dummy68, "68% CL", "f");
  leg.AddEntry(&dummyNom, "Nominal", "l");
  leg.AddEntry(&dummyBest, "Best fit", "l");
  leg.Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.038);
  latex.DrawLatex(0.13, 0.945, "#bf{CMS} #it{Simulation Work in Progress}");
  latex.DrawLatex(0.70, 0.945, "138 fb^{-1} (13 TeV)");

  c.SaveAs(outpdf.c_str());
  c.SaveAs(outpng.c_str());
}

FitResult1D fit_one_wc(const std::string& wc,
                       const std::vector<int>& obs,
                       const std::string& tag,
                       const std::string& outdir,
                       const std::string& data_root,
                       const std::string& eft_template_pattern,
                       const TMatrixD& cov_full,
                       int drop_bin_idx,
                       double scan_min,
                       double scan_max,
                       int scan_n,
                       double xscale = 1.0,
                       const std::string& xlabel_override = "") {
  std::vector<int> keep = build_keep_indices(obs, drop_bin_idx);

  auto data_full = load_values(data_root);
  auto data = select_vec(data_full, keep);

  TMatrixD cov = select_cov(cov_full, keep, drop_bin_idx);
  TDecompSVD svd(cov);
  TMatrixD cov_inv = svd.Invert();

  std::vector<int> wc_vals = {-8, -4, -2, 0, 2, 4, 8};
  std::map<int, std::vector<double>> tmpl;

  for (int v : wc_vals) {
    std::string path = make_template_path(eft_template_pattern, wc, v);
    tmpl[v] = select_vec(load_values(path), keep);
  }

  const int nb = data.size();

  TMatrixD X(wc_vals.size(), 3);
  for (int i = 0; i < (int)wc_vals.size(); ++i) {
    double x = wc_vals[i];
    X(i, 0) = 1.0;
    X(i, 1) = x;
    X(i, 2) = x * x;
  }

  TDecompSVD xsvd(X);
  TMatrixD Xpinv = xsvd.Invert();

  std::vector<double> A(nb), B(nb), C(nb);
  for (int ib = 0; ib < nb; ++ib) {
    TVectorD y(wc_vals.size());
    for (int i = 0; i < (int)wc_vals.size(); ++i) y(i) = tmpl[wc_vals[i]][ib];
    TVectorD beta = Xpinv * y;
    A[ib] = beta(0);
    B[ib] = beta(1);
    C[ib] = beta(2);
  }

  auto pred_at = [&](double x) {
    std::vector<double> p(nb);
    for (int i = 0; i < nb; ++i) p[i] = A[i] + B[i] * x + C[i] * x * x;
    return block_renorm(p, tmpl[0], BINS_PER_OBS - 1);
  };

  auto chi2_at = [&](double x) {
    auto p = pred_at(x);
    TVectorD d(nb);
    for (int i = 0; i < nb; ++i) d(i) = data[i] - p[i];
    TVectorD tmp = cov_inv * d;
    return d * tmp;
  };

  std::vector<double> xs(scan_n), chi(scan_n), dchi(scan_n);
  double best = 0.0, chi_min = 1e300;

  for (int i = 0; i < scan_n; ++i) {
    double x = scan_min + (scan_max - scan_min) * double(i) / double(scan_n - 1);
    xs[i] = x;
    chi[i] = chi2_at(x);
    if (chi[i] < chi_min) {
      chi_min = chi[i];
      best = x;
    }
  }

  for (int i = 0; i < scan_n; ++i) dchi[i] = chi[i] - chi_min;

  auto i68 = interval_crossing(xs, dchi, 1.0, best);
  auto i95 = interval_crossing(xs, dchi, 4.0, best);

  FitResult1D r;
  r.name = wc;
  r.label = pretty_label(wc);
  r.best = best;
  r.chi2min = chi_min;
  r.lo68 = i68.first;
  r.hi68 = i68.second;
  r.lo95 = i95.first;
  r.hi95 = i95.second;
  r.nbins = nb;

  const std::string od = outdir + "/" + tag;
  gSystem->mkdir(od.c_str(), true);

  // ------------------------------------------------------------
  // Parabola-basis sanity check for the concatenated histogram fit.
  // The official/result fit above remains the global all-template fit
  // using {-8,-4,-2,0,2,4,8}.  This diagnostic only compares how the
  // chi2 parabola moves if the bin-wise EFT interpolation is rebuilt
  // from {-4,-2,0,2,4} or {-2,0,2}.
  // ------------------------------------------------------------
  auto solve_quadratic_basis = [&](const std::vector<int>& vals,
                                   std::vector<double>& AA,
                                   std::vector<double>& BB,
                                   std::vector<double>& CC) {
    TMatrixD XX(vals.size(), 3);
    for (int i = 0; i < (int)vals.size(); ++i) {
      const double xx = vals[i];
      XX(i, 0) = 1.0;
      XX(i, 1) = xx;
      XX(i, 2) = xx * xx;
    }
    TMatrixD XXpinv = TDecompSVD(XX).Invert();
    AA.assign(nb, 0.0);
    BB.assign(nb, 0.0);
    CC.assign(nb, 0.0);
    for (int ib = 0; ib < nb; ++ib) {
      TVectorD yy(vals.size());
      for (int i = 0; i < (int)vals.size(); ++i) yy(i) = tmpl[vals[i]][ib];
      TVectorD beta = XXpinv * yy;
      AA[ib] = beta(0);
      BB[ib] = beta(1);
      CC[ib] = beta(2);
    }
  };

  auto pred_from_basis = [&](const std::vector<double>& AA,
                             const std::vector<double>& BB,
                             const std::vector<double>& CC,
                             double xx) {
    std::vector<double> p(nb);
    for (int i = 0; i < nb; ++i) p[i] = AA[i] + BB[i] * xx + CC[i] * xx * xx;
    return block_renorm(p, tmpl[0], BINS_PER_OBS - 1);
  };

  auto chi2_from_basis = [&](const std::vector<double>& AA,
                             const std::vector<double>& BB,
                             const std::vector<double>& CC,
                             double xx) {
    auto p = pred_from_basis(AA, BB, CC, xx);
    TVectorD d(nb);
    for (int i = 0; i < nb; ++i) d(i) = data[i] - p[i];
    TVectorD tmp = cov_inv * d;
    return d * tmp;
  };

  std::vector<int> vals_all7  = {-8, -4, -2, 0, 2, 4, 8};
  std::vector<int> vals_mid5  = {-4, -2, 0, 2, 4};
  std::vector<int> vals_pm2   = {-2, 0, 2};
  std::vector<double> A_all7=A, B_all7=B, C_all7=C;
  std::vector<double> A_mid5, B_mid5, C_mid5, A_pm2, B_pm2, C_pm2;
  solve_quadratic_basis(vals_mid5, A_mid5, B_mid5, C_mid5);
  solve_quadratic_basis(vals_pm2,  A_pm2,  B_pm2,  C_pm2);

  auto build_scan_for_basis = [&](const std::vector<double>& AA,
                                  const std::vector<double>& BB,
                                  const std::vector<double>& CC,
                                  std::vector<double>& cc,
                                  std::vector<double>& dc,
                                  double& bbest,
                                  double& cmin) {
    cc.assign(scan_n, 0.0);
    dc.assign(scan_n, 0.0);
    bbest = 0.0;
    cmin = 1e300;
    for (int i = 0; i < scan_n; ++i) {
      const double xnow = xs[i];
      cc[i] = chi2_from_basis(AA, BB, CC, xnow);
      if (cc[i] < cmin) {
        cmin = cc[i];
        bbest = xnow;
      }
    }
    for (int i = 0; i < scan_n; ++i) dc[i] = cc[i] - cmin;
  };

  std::vector<double> chi_all7, dchi_all7, chi_mid5, dchi_mid5, chi_pm2, dchi_pm2;
  double best_all7=0.0, min_all7=0.0, best_mid5=0.0, min_mid5=0.0, best_pm2=0.0, min_pm2=0.0;
  build_scan_for_basis(A_all7, B_all7, C_all7, chi_all7, dchi_all7, best_all7, min_all7);
  build_scan_for_basis(A_mid5, B_mid5, C_mid5, chi_mid5, dchi_mid5, best_mid5, min_mid5);
  build_scan_for_basis(A_pm2,  B_pm2,  C_pm2,  chi_pm2,  dchi_pm2,  best_pm2,  min_pm2);

  {
    std::ofstream pcsv(od + "/parabola_basis_compare_" + wc + ".csv");
    pcsv << wc << ",chi2_all7,dchi2_all7,chi2_mid5,dchi2_mid5,chi2_pm2,dchi2_pm2\n";
    for (int i = 0; i < scan_n; ++i) {
      pcsv << xs[i] << "," << chi_all7[i] << "," << dchi_all7[i]
           << "," << chi_mid5[i] << "," << dchi_mid5[i]
           << "," << chi_pm2[i] << "," << dchi_pm2[i] << "\n";
    }
  }
  {
    std::ofstream psum(od + "/parabola_basis_summary_" + wc + ".csv");
    psum << "basis,template_points,best,chi2min\n";
    psum << "all7,{-8,-4,-2,0,2,4,8}," << best_all7 << "," << min_all7 << "\n";
    psum << "mid5,{-4,-2,0,2,4}," << best_mid5 << "," << min_mid5 << "\n";
    psum << "pm2,{-2,0,2}," << best_pm2 << "," << min_pm2 << "\n";
  }

  {
    std::vector<double> xplot(scan_n);
    for (int i = 0; i < scan_n; ++i) xplot[i] = xs[i] * xscale;

    const double yplot_min = 0.0;
    const double yplot_max = 8.0;

    auto find_interval_68 = [&](const std::vector<double>& yy,
                                double best,
                                double& lo,
                                double& hi) {
      lo = best;
      hi = best;

      int ibest = 0;
      double best_dist = 1e300;
      for (int i = 0; i < scan_n; ++i) {
        const double dist = std::fabs(xs[i] - best);
        if (dist < best_dist) {
          best_dist = dist;
          ibest = i;
        }
      }

      // Left crossing of Delta chi2 = 1
      bool found_lo = false;
      for (int i = ibest; i > 0; --i) {
        const double y1 = yy[i];
        const double y0 = yy[i - 1];
        if ((y1 - 1.0) * (y0 - 1.0) <= 0.0 && std::fabs(y1 - y0) > 1e-14) {
          const double t = (1.0 - y0) / (y1 - y0);
          lo = xs[i - 1] + t * (xs[i] - xs[i - 1]);
          found_lo = true;
          break;
        }
      }
      if (!found_lo) lo = xs.front();

      // Right crossing of Delta chi2 = 1
      bool found_hi = false;
      for (int i = ibest; i + 1 < scan_n; ++i) {
        const double y0 = yy[i];
        const double y1 = yy[i + 1];
        if ((y0 - 1.0) * (y1 - 1.0) <= 0.0 && std::fabs(y1 - y0) > 1e-14) {
          const double t = (1.0 - y0) / (y1 - y0);
          hi = xs[i] + t * (xs[i + 1] - xs[i]);
          found_hi = true;
          break;
        }
      }
      if (!found_hi) hi = xs.back();
    };

    double lo_all7 = 0.0, hi_all7 = 0.0;
    double lo_mid5 = 0.0, hi_mid5 = 0.0;
    double lo_pm2  = 0.0, hi_pm2  = 0.0;
    find_interval_68(dchi_all7, best_all7, lo_all7, hi_all7);
    find_interval_68(dchi_mid5, best_mid5, lo_mid5, hi_mid5);
    find_interval_68(dchi_pm2,  best_pm2,  lo_pm2,  hi_pm2);

    {
      std::ofstream psum2(od + "/parabola_basis_summary_" + wc + ".csv", std::ios::app);
      psum2 << "#68CL_summary,basis,best,lo68,hi68,minus,plus\n";
      psum2 << "68CL,all7," << best_all7 << "," << lo_all7 << "," << hi_all7 << "," << (best_all7 - lo_all7) << "," << (hi_all7 - best_all7) << "\n";
      psum2 << "68CL,mid5," << best_mid5 << "," << lo_mid5 << "," << hi_mid5 << "," << (best_mid5 - lo_mid5) << "," << (hi_mid5 - best_mid5) << "\n";
      psum2 << "68CL,pm2,"  << best_pm2  << "," << lo_pm2  << "," << hi_pm2  << "," << (best_pm2  - lo_pm2)  << "," << (hi_pm2  - best_pm2)  << "\n";
    }

    TCanvas cpar("cpar", "cpar", 1120, 760);
    cpar.SetLeftMargin(0.12);
    cpar.SetBottomMargin(0.12);
    cpar.SetRightMargin(0.38);
    cpar.SetTopMargin(0.09);
    cpar.SetTicks(1, 1);

    TGraph gAll(scan_n, xplot.data(), dchi_all7.data());
    TGraph gMid(scan_n, xplot.data(), dchi_mid5.data());
    TGraph gPm2(scan_n, xplot.data(), dchi_pm2.data());
    gAll.SetName(("g_parabola_all7_" + wc).c_str());
    gMid.SetName(("g_parabola_mid5_" + wc).c_str());
    gPm2.SetName(("g_parabola_pm2_" + wc).c_str());

    const int lineAll = kBlack;
    const int lineMid = kBlue + 1;
    const int linePm2 = kRed + 1;
    const int markAll = kGray + 2;
    const int markMid = kAzure + 7;
    const int markPm2 = kOrange + 7;

    gAll.SetLineColor(lineAll); gAll.SetLineWidth(2); gAll.SetLineStyle(1);
    gMid.SetLineColor(lineMid); gMid.SetLineWidth(2); gMid.SetLineStyle(1);
    gPm2.SetLineColor(linePm2); gPm2.SetLineWidth(2); gPm2.SetLineStyle(1);
    gAll.SetTitle("");
    gAll.GetXaxis()->SetTitle(xlabel_override.empty() ? pretty_label(wc).c_str() : xlabel_override.c_str());
    gAll.GetYaxis()->SetTitle("#Delta#chi^{2}");
    gAll.GetYaxis()->SetRangeUser(yplot_min, yplot_max);
    gAll.GetYaxis()->SetTitleOffset(0.95);
    gAll.GetXaxis()->SetTitleOffset(1.05);
    gAll.Draw("AL");
    gMid.Draw("L SAME");
    gPm2.Draw("L SAME");

    TLine l68(scan_min * xscale, 1.0, scan_max * xscale, 1.0);
    l68.SetLineStyle(3); l68.SetLineColor(kGray+2); l68.SetLineWidth(1); l68.Draw("SAME");
    TLine l95(scan_min * xscale, 4.0, scan_max * xscale, 4.0);
    l95.SetLineStyle(3); l95.SetLineColor(kGray+2); l95.SetLineWidth(1); l95.Draw("SAME");

    // Do NOT draw vertical best-fit / 68% CL interval lines here.
    // They were visually mistaken for extra parabolas because there are
    // three basis choices and each would add three vertical lines.  The
    // best fit and 68% CL interval are kept in the legend and CSV summary.
    // Draw only a small marker at the minimum of each parabola.
    auto draw_best_marker = [&](double best, int color, int style) {
      TMarker* m = new TMarker(best * xscale, 0.0, style);
      m->SetMarkerColor(color);
      m->SetMarkerSize(1.65);
      m->Draw("SAME");
    };
    draw_best_marker(best_all7, lineAll, 29);
    draw_best_marker(best_mid5, lineMid, 29);
    draw_best_marker(best_pm2,  linePm2, 29);

    auto draw_markers = [&](const std::vector<int>& vals,
                            const std::vector<double>& AA,
                            const std::vector<double>& BB,
                            const std::vector<double>& CC,
                            double cmin,
                            int color,
                            int style) {
      for (int vv : vals) {
        const double yy = chi2_from_basis(AA, BB, CC, (double)vv) - cmin;
        if (yy < yplot_min || yy > yplot_max) continue;
        TMarker* m = new TMarker(vv * xscale, yy, style);
        m->SetMarkerColor(color);
        m->SetMarkerSize(1.45);
        m->Draw("SAME");
      }
    };
    // Draw markers last so template points stay above the curves and vertical interval lines.
    draw_markers(vals_all7, A_all7, B_all7, C_all7, min_all7, markAll, 20);
    draw_markers(vals_mid5, A_mid5, B_mid5, C_mid5, min_mid5, markMid, 21);
    draw_markers(vals_pm2,  A_pm2,  B_pm2,  C_pm2,  min_pm2,  markPm2, 22);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.034);
    latex.DrawLatex(0.12, 0.955, ("#bf{CMS} #it{Simulation Work in Progress}  " + tag).c_str());
    latex.SetTextAlign(31);
    latex.DrawLatex(0.96, 0.955, "138 fb^{-1} (13 TeV)");
    latex.SetTextAlign(11);
    latex.SetTextSize(0.022);
    // latex.DrawLatex(0.125, 0.83, "95% CL");
    // latex.DrawLatex(0.125, 0.46, "68% CL");

    auto fit_label = [](const char* set_label, double best, double lo, double hi) {
      return Form("%s: %.3g^{+%.2g}_{-%.2g}",
                  set_label, best, hi - best, best - lo);
    };

    TGraph legAll(1), legMid(1), legPm2(1);
    legAll.SetLineColor(lineAll); legAll.SetLineWidth(2); legAll.SetLineStyle(1);
    legAll.SetMarkerColor(markAll); legAll.SetMarkerStyle(20); legAll.SetMarkerSize(1.35);
    legMid.SetLineColor(lineMid); legMid.SetLineWidth(2); legMid.SetLineStyle(1);
    legMid.SetMarkerColor(markMid); legMid.SetMarkerStyle(21); legMid.SetMarkerSize(1.35);
    legPm2.SetLineColor(linePm2); legPm2.SetLineWidth(2); legPm2.SetLineStyle(1);
    legPm2.SetMarkerColor(markPm2); legPm2.SetMarkerStyle(22); legPm2.SetMarkerSize(1.35);

    TLegend leg(0.705, 0.55, 0.995, 0.89);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.019);
    leg.AddEntry(&legAll, fit_label("{-8,-4,-2,0,2,4,8}", best_all7, lo_all7, hi_all7), "lp");
    leg.AddEntry(&legMid, fit_label("{-4,-2,0,2,4}", best_mid5, lo_mid5, hi_mid5), "lp");
    leg.AddEntry(&legPm2, fit_label("{-2,0,2}", best_pm2, lo_pm2, hi_pm2), "lp");
    leg.Draw();

    cpar.RedrawAxis();
    cpar.SaveAs((od + "/parabola_basis_compare_" + wc + ".png").c_str());
    cpar.SaveAs((od + "/parabola_basis_compare_" + wc + ".pdf").c_str());
  }

  std::cout << "[SANITY parabola] " << tag << " " << wc
            << " global/all7 best=" << best_all7
            << " mid5 best=" << best_mid5
            << " pm2 best=" << best_pm2
            << " (reported FIT result remains global/all7)"
            << std::endl;

  
std::ofstream csv(od + "/scan_" + wc + ".csv");
csv << wc << ",chi2,delta_chi2\n";
for (int i = 0; i < scan_n; ++i) {
  csv << xs[i] << "," << chi[i] << "," << dchi[i] << "\n";
}
csv.close();


  
std::ofstream js(od + "/fit_result_" + wc + ".json");
js << "{\n";
js << "  \"wc\": \"" << wc << "\",\n";
js << "  \"tag\": \"" << tag << "\",\n";
js << "  \"best\": " << r.best << ",\n";
js << "  \"chi2_min\": " << r.chi2min << ",\n";
js << "  \"lo68\": " << r.lo68 << ",\n";
js << "  \"hi68\": " << r.hi68 << ",\n";
js << "  \"lo95\": " << r.lo95 << ",\n";
js << "  \"hi95\": " << r.hi95 << ",\n";
js << "  \"n_bins_fit\": " << r.nbins << "\n";
js << "}\n";
js.close();


  save_chi2_plot(
    od + "/deltaChi2_" + wc + ".pdf",
    od + "/deltaChi2_" + wc + ".png",
    xs, dchi, r, xscale, xlabel_override
  );

  std::cout << "[FIT] " << tag << " " << wc
            << " best=" << best
            << " 68=[" << r.lo68 << "," << r.hi68 << "]"
            << " 95=[" << r.lo95 << "," << r.hi95 << "]"
            << " nbins=" << nb << std::endl;

  return r;
}

void fit_2d_pair_grid(const std::string& wc1,
                      const std::string& wc2,
                      const std::vector<int>& obs,
                      const std::string& tag,
                      const std::string& outdir,
                      const std::string& data_root,
                      const std::string& eft_template_pattern,
                      const TMatrixD& cov_full,
                      int drop_bin_idx,
                      double scan_min,
                      double scan_max,
                      int ngrid) {
  std::vector<int> keep = build_keep_indices(obs, drop_bin_idx);
  auto data = select_vec(load_values(data_root), keep);

  TMatrixD cov = select_cov(cov_full, keep, drop_bin_idx);
  TDecompSVD svd(cov);
  TMatrixD cov_inv = svd.Invert();

  auto load_quad = [&](const std::string& wc, std::vector<double>& A, std::vector<double>& B, std::vector<double>& C, std::vector<double>& ref) {
    std::vector<int> vals = {-8,-4,-2,0,2,4,8};
    std::map<int, std::vector<double>> tmpl;
    for (int v : vals) tmpl[v] = select_vec(load_values(make_template_path(eft_template_pattern, wc, v)), keep);
    ref = tmpl[0];

    const int nb = ref.size();
    A.assign(nb, 0.0); B.assign(nb, 0.0); C.assign(nb, 0.0);

    TMatrixD X(vals.size(), 3);
    for (int i = 0; i < (int)vals.size(); ++i) {
      double x = vals[i];
      X(i,0)=1; X(i,1)=x; X(i,2)=x*x;
    }
    TMatrixD Xpinv = TDecompSVD(X).Invert();

    for (int ib = 0; ib < nb; ++ib) {
      TVectorD y(vals.size());
      for (int i = 0; i < (int)vals.size(); ++i) y(i) = tmpl[vals[i]][ib];
      TVectorD beta = Xpinv * y;
      A[ib]=beta(0); B[ib]=beta(1); C[ib]=beta(2);
    }
  };

  std::vector<double> A1,B1,C1,R1,A2,B2,C2,R2;
  load_quad(wc1,A1,B1,C1,R1);
  load_quad(wc2,A2,B2,C2,R2);

  const int nb = data.size();

  auto pred = [&](double x, double y) {
    std::vector<double> p(nb);
    for (int i = 0; i < nb; ++i) {
      double p1 = A1[i] + B1[i]*x + C1[i]*x*x;
      double p2 = A2[i] + B2[i]*y + C2[i]*y*y;
      p[i] = R1[i] + (p1 - R1[i]) + (p2 - R2[i]);
    }
    return block_renorm(p, R1, BINS_PER_OBS - 1);
  };

  auto chi2 = [&](double x, double y) {
    auto p = pred(x,y);
    TVectorD d(nb);
    for (int i = 0; i < nb; ++i) d(i)=data[i]-p[i];
    TVectorD tmp = cov_inv*d;
    return d*tmp;
  };

  TH2D h("h","",ngrid,scan_min,scan_max,ngrid,scan_min,scan_max);
  double cmin = 1e300, bx = 0, by = 0;

  for (int ix = 1; ix <= ngrid; ++ix) {
    double x = h.GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= ngrid; ++iy) {
      double y = h.GetYaxis()->GetBinCenter(iy);
      double c = chi2(x,y);
      h.SetBinContent(ix,iy,c);
      if (c < cmin) { cmin = c; bx=x; by=y; }
    }
  }

  for (int ix = 1; ix <= ngrid; ++ix)
    for (int iy = 1; iy <= ngrid; ++iy)
      h.SetBinContent(ix,iy,h.GetBinContent(ix,iy)-cmin);

  std::string od = outdir + "/" + tag + "_2D";
  gSystem->mkdir(od.c_str(), true);

  TCanvas c("c2d","c2d",760,700);
  
  gStyle->SetOptStat(0);c.SetLeftMargin(0.13);
  c.SetRightMargin(0.15);
  c.SetTopMargin(0.11);
  c.SetBottomMargin(0.13);
  
  h.SetTitle("");
  h.GetXaxis()->SetTitle(pretty_label(wc1).c_str());
  h.GetYaxis()->SetTitle(pretty_label(wc2).c_str());
  h.GetZaxis()->SetTitle("#Delta#chi^{2}");
  h.SetMinimum(0);
  h.SetMaximum(10);
  h.Draw("AXIS");

  double levels[2] = {2.30, 5.99};
  h.SetContour(2, levels);
  h.SetContour(2);
  h.SetContourLevel(0, 2.30);
  h.SetContourLevel(1, 5.99);
  h.SetLineColor(TColor::GetColor("#0057ff"));
  h.SetLineWidth(3);
  h.Draw("CONT3 same");

  TMarker m(bx,by,29);
  m.SetMarkerColor(kBlack);
  m.SetMarkerSize(1.8);
  m.Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.032);
  latex.DrawLatex(0.13,0.955,"#bf{CMS} #it{Simulation Work in Progress}");
  latex.DrawLatex(0.70,0.955,"138 fb^{-1} (13 TeV)");

  c.SaveAs((od + "/contour_" + wc1 + "_vs_" + wc2 + ".pdf").c_str());
  c.SaveAs((od + "/contour_" + wc1 + "_vs_" + wc2 + ".png").c_str());

  std::ofstream js(od + "/fit2d_" + wc1 + "_vs_" + wc2 + ".json");
  js << "{\n";
  js << "  \"wc1\": \"" << wc1 << "\",\n";
  js << "  \"wc2\": \"" << wc2 << "\",\n";
  js << "  \"best1\": " << bx << ",\n";
  js << "  \"best2\": " << by << ",\n";
  js << "  \"chi2_min\": " << cmin << ",\n";
  js << "  \"n_bins_fit\": " << nb << "\n";
  js << "}\n";
  js.close();

  std::cout << "[2D] " << tag << " " << wc1 << " vs " << wc2
            << " best=(" << bx << "," << by << ")"
            << " chi2min=" << cmin << std::endl;
}


struct SummaryEntry {
  std::string key;
  std::string label;
  double best = 0.0, lo68 = 0.0, hi68 = 0.0;
};

struct PubEntry {
  std::string key;
  double best = 0.0;
  double err68 = 0.0;
};

PubEntry pub_from_interval(const std::string& key, double lo, double hi) {
  PubEntry p;
  p.key = key;
  p.best = 0.5 * (lo + hi);
  p.err68 = 0.5 * (hi - lo);
  return p;
}

std::vector<PubEntry> top22_wc_pub_entries() {
  // TOP-22-006, Table 5, profiled 1 sigma intervals.
  // Only direct same-basis/approximately same-name WCs from the current 16-WC list are included.
  // The TOP-22-006 paper uses SMEFTsim-style names: ctG, c_tq^{1,8}, c_Qq^{11,18,31,38}.
  return {
    pub_from_interval("ctGRe", -0.15, 0.12),
    pub_from_interval("cQj18", -0.47, -0.00),
    pub_from_interval("cQj38", -0.09, 0.08),
    pub_from_interval("cQj11", -0.10, 0.10),
    pub_from_interval("cQj31", -0.04, 0.03),
    pub_from_interval("ctj8",  -0.45, 0.03),
    pub_from_interval("ctj1",  -0.11, 0.11)
  };
}

SummaryEntry make_summary_entry(const std::string& key, const std::string& label,
                                const FitResult1D& r, double scale = 1.0) {
  SummaryEntry e;
  e.key = key; e.label = label;
  e.best = r.best * scale; e.lo68 = r.lo68 * scale; e.hi68 = r.hi68 * scale;
  return e;
}

bool template_exists(const std::string& eft_template_pattern, const std::string& wc, int val = 0) {
  const std::string path = make_template_path(eft_template_pattern, wc, val);
  return !gSystem->AccessPathName(path.c_str());
}

void save_summary_plot(const std::vector<SummaryEntry>& entries,
                       const std::string& outpdf, const std::string& outpng,
                       const std::string& xtitle, double xmin, double xmax,
                       const std::string& cms_extra = "Simulation Work in Progress",
                       const std::vector<PubEntry>& pub_entries = std::vector<PubEntry>()) {
  if (entries.empty()) return;
  gStyle->SetOptStat(0);
  const int n = entries.size();
  TH2D frame("summary_frame", "", 1, xmin, xmax, n, 0, n);
  frame.SetTitle("");
  frame.GetXaxis()->SetTitle(xtitle.c_str());
  frame.GetYaxis()->SetTitle("");
  frame.GetXaxis()->SetTitleSize(0.055);
  frame.GetXaxis()->SetLabelSize(0.045);
  frame.GetYaxis()->SetLabelSize(0.050);
  frame.GetYaxis()->SetTickLength(0.0);
  for (int i = 0; i < n; ++i) frame.GetYaxis()->SetBinLabel(n - i, entries[i].label.c_str());
  TCanvas c("c_summary", "c_summary", 820, std::max(620, 52*n + 170));
  
  gStyle->SetOptStat(0);c.SetLeftMargin(0.30); c.SetRightMargin(0.04); c.SetTopMargin(0.16); c.SetBottomMargin(0.13);
  frame.Draw("AXIS");
  const int green = TColor::GetColor("#2ca25f");
  TLine sm(0.0, 0.0, 0.0, n); sm.SetLineColor(green); sm.SetLineWidth(3); sm.Draw("same");
  std::vector<double> x(n), y(n), exlo(n), exhi(n), eylo(n,0.0), eyhi(n,0.0);
  for (int i = 0; i < n; ++i) {
    const auto& e = entries[i];
    x[i] = e.best; y[i] = n - i - 0.5;
    exlo[i] = std::max(0.0, e.best - e.lo68);
    exhi[i] = std::max(0.0, e.hi68 - e.best);
  }
  TGraphAsymmErrors gr(n, x.data(), y.data(), exlo.data(), exhi.data(), eylo.data(), eyhi.data());
  gr.SetMarkerStyle(20); gr.SetMarkerSize(1.05); gr.SetMarkerColor(kBlack); gr.SetLineColor(kBlack); gr.SetLineWidth(2);

  std::map<std::string, PubEntry> pub_map;
  for (const auto& p : pub_entries) pub_map[p.key] = p;
  std::vector<double> px, py, pexlo, pexhi, peylo, peyhi;
  for (int i = 0; i < n; ++i) {
    auto it = pub_map.find(entries[i].key);
    if (it == pub_map.end()) continue;
    px.push_back(it->second.best);
    py.push_back(n - i - 0.5);
    pexlo.push_back(it->second.err68);
    pexhi.push_back(it->second.err68);
    peylo.push_back(0.0);
    peyhi.push_back(0.0);
  }

  TGraphAsymmErrors gr_pub;
  if (!px.empty()) {
    gr_pub = TGraphAsymmErrors(px.size(), px.data(), py.data(), pexlo.data(), pexhi.data(), peylo.data(), peyhi.data());
    gr_pub.SetMarkerStyle(24);
    gr_pub.SetMarkerSize(1.05);
    gr_pub.SetMarkerColor(TColor::GetColor("#5b7fd9"));
    gr_pub.SetLineColor(TColor::GetColor("#5b7fd9"));
    gr_pub.SetLineWidth(2);
    gr_pub.Draw("P same");
  }

  gr.Draw("P same");

  // Place the anomalous-coupling legend slightly lower than the WC-summary legend.
  const bool isAnomSummary = (xtitle.find("Anomalous") != std::string::npos);
  const double legY1 = isAnomSummary ? 0.48 : 0.58;
  const double legY2 = isAnomSummary ? 0.62 : 0.72;
  TLegend leg(0.72, legY1, 0.965, legY2);
  leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.020); leg.SetMargin(0.16); leg.SetEntrySeparation(0.03);

  TLine dummySM; dummySM.SetLineColor(green); dummySM.SetLineWidth(3);
  TGraph dummyPub; dummyPub.SetMarkerStyle(24); dummyPub.SetMarkerColor(TColor::GetColor("#5b7fd9")); dummyPub.SetLineColor(TColor::GetColor("#5b7fd9")); dummyPub.SetLineWidth(2);
  TGraph dummyFit; dummyFit.SetMarkerStyle(20); dummyFit.SetMarkerColor(kBlack); dummyFit.SetLineColor(kBlack); dummyFit.SetLineWidth(2);
  leg.AddEntry(&dummySM, "Standard model", "l");

  const bool isTheoryLO  = (cms_extra.find("Theory LO")  != std::string::npos);
  const bool isTheoryNLO = (cms_extra.find("Theory NLO") != std::string::npos);
  const bool isWilsonSummary = (xtitle.find("Wilson coefficient") != std::string::npos);
  const char* pubLabel = isWilsonSummary ? "TOP-22-006 (1#sigma)" : "TOP-18-006 (68% CL)";
  if (!pub_entries.empty()) leg.AddEntry(&dummyPub, pubLabel, "lep");

  std::string fitLabel = "Dim6Top MC";
  if (isTheoryLO)  fitLabel = "Theory LO";
  if (isTheoryNLO) fitLabel = "Theory NLO";
  fitLabel += " (68% CL)";
  leg.AddEntry(&dummyFit, fitLabel.c_str(), "lep");
  leg.Draw();

  TLatex latex; latex.SetNDC(); latex.SetTextFont(42); latex.SetTextSize(0.034);
  latex.DrawLatex(0.30, 0.965, "#bf{CMS}");
  latex.SetTextSize(0.028);
  latex.DrawLatex(0.43, 0.965, "#it{Simulation Work in Progress}");
  latex.SetTextAlign(31);
  latex.DrawLatex(0.94, 0.965, "138 fb^{-1} (13 TeV)");
  latex.SetTextAlign(11);
  c.SaveAs(outpdf.c_str()); c.SaveAs(outpng.c_str());
}

void write_summary_csv(const std::vector<SummaryEntry>& entries, const std::string& path) {
  std::ofstream f(path);
  f << "key,label,best,lo68,hi68\n";
  for (const auto& e : entries) f << e.key << ",\"" << e.label << "\"," << e.best << "," << e.lo68 << "," << e.hi68 << "\n";
}

struct PairFitSpec {
  std::string wc1, wc2, tag;
  std::vector<int> obs;
  double xmin, xmax;
};

void fit_2d_pair_grid_if_available(const PairFitSpec& spec,
                                   const std::string& outdir,
                                   const std::string& data_root,
                                   const std::string& eft_template_pattern,
                                   const TMatrixD& cov_full,
                                   int drop_bin_idx,
                                   int scan2d_n) {
  if (!template_exists(eft_template_pattern, spec.wc1, 0) || !template_exists(eft_template_pattern, spec.wc2, 0)) {
    std::cout << "[SKIP 2D] missing template(s) for " << spec.wc1 << " vs " << spec.wc2 << std::endl;
    return;
  }
  fit_2d_pair_grid(spec.wc1, spec.wc2, spec.obs, spec.tag, outdir, data_root, eft_template_pattern,
                   cov_full, drop_bin_idx, spec.xmin, spec.xmax, scan2d_n);
}

// ============================================================
// Theory-table integration helpers
//   - The parse_latex_tables.py output is embedded below as CSV columns:
//       order,op,observable,kind,central,up,down
//     where kind is one of sm,lin,quad.
//   - We turn the inclusive functional form
//       (N0 + c N1 + c^2 N2)/(D0 + c D1 + c^2 D2)
//     into a 6-bin toy/asymmetry template using the same A_FB convention
//     as anom_chi2_fit_nanogen_Alex_test_integrated.py.
// ============================================================

struct HistChunk1D {
  std::vector<double> y;
  std::vector<double> e;
};

struct CoeffSummary {
  bool ok = false;
  double afb = 0.0;
  double afb_err = 0.0;
  double coeff = 0.0;
  double coeff_err = 0.0;
  double factor = 0.0;
};

struct TheoryTriplet {
  double sm = 0.0, lin = 0.0, quad = 0.0;
  double sm_err = 0.0, lin_err = 0.0, quad_err = 0.0;
  bool has_sm=false, has_lin=false, has_quad=false;
};
using TheoryTable = std::map<std::string, std::map<std::string, std::map<std::string, TheoryTriplet>>>;

// Embedded output of parse_latex_tables.py.
// This removes the runtime dependency on a separate Python pre-step.
static const char* EMBEDDED_THEORY_CSV = R"THEORYCSV(
order,op,observable,kind,central,up,down
LO,SM,sigma,sm,29.12,9.1,6.49
LO,SM,sigma_B_n+barB_n,sm,-0.09,0.02,0.03
LO,SM,sigma_B_r+barB_r,sm,0.04,0.01,0.01
LO,SM,sigma_B_k+barB_k,sm,0.11,0.03,0.02
LO,SM,sigma_C_nn,sm,9.63,2.96,2.13
LO,SM,sigma_C_rr,sm,0.19,0.19,0.13
LO,SM,sigma_C_kk,sm,9.42,3.12,2.22
LO,SM,sigma_C_rk+C_kr,sm,-6.81,1.45,2.01
LO,SM,sigma_C_nr+C_rn,sm,-0.06,0.01,0.02
LO,SM,sigma_C_nk+C_kn,sm,0.11,0.04,0.03
LO,c_tG,sigma,lin,9.32,2.9,2.08
LO,c_tG,sigma,quad,1.63,0.52,0.37
LO,c_tG,sigma_B_n+barB_n,lin,-0.01,0.0,0.0
LO,c_tG,sigma_B_n+barB_n,quad,-0.0,0.0,0.0
LO,c_tG,sigma_B_r+barB_r,lin,-0.01,0.0,0.0
LO,c_tG,sigma_B_r+barB_r,quad,-0.0,0.0,0.0
LO,c_tG,sigma_B_k+barB_k,lin,-0.02,0.01,0.01
LO,c_tG,sigma_B_k+barB_k,quad,-0.01,0.0,0.0
LO,c_tG,sigma_C_nn,lin,6.68,2.13,1.52
LO,c_tG,sigma_C_nn,quad,0.4,0.13,0.09
LO,c_tG,sigma_C_rr,lin,4.37,1.44,1.02
LO,c_tG,sigma_C_rr,quad,0.65,0.21,0.15
LO,c_tG,sigma_C_kk,lin,4.6,1.47,1.05
LO,c_tG,sigma_C_kk,quad,0.97,0.31,0.22
LO,c_tG,sigma_C_rk+C_kr,lin,-0.65,0.11,0.15
LO,c_tG,sigma_C_rk+C_kr,quad,-0.03,0.01,0.01
LO,c_tG,sigma_C_nr+C_rn,lin,-0.02,0.01,0.01
LO,c_tG,sigma_C_nr+C_rn,quad,-0.0,0.0,0.0
LO,c_tG,sigma_C_nk+C_kn,lin,-0.04,0.01,0.02
LO,c_tG,sigma_C_nk+C_kn,quad,-0.01,0.0,0.0
LO,cQq83,sigma,lin,0.073,0.008,0.007
LO,cQq83,sigma,quad,0.075,0.004,0.004
LO,cQq83,sigma_B_n+barB_n,lin,0.002,0.0,0.0
LO,cQq83,sigma_B_n+barB_n,quad,-0.003,0.001,0.001
LO,cQq83,sigma_B_r+barB_r,lin,-0.008,0.001,0.001
LO,cQq83,sigma_B_r+barB_r,quad,-0.011,0.001,0.0
LO,cQq83,sigma_B_k+barB_k,lin,-0.079,0.008,0.011
LO,cQq83,sigma_B_k+barB_k,quad,-0.133,0.011,0.013
LO,cQq83,sigma_C_nn,lin,0.019,0.002,0.001
LO,cQq83,sigma_C_nn,quad,0.002,0.0,0.0
LO,cQq83,sigma_C_rr,lin,-0.031,0.003,0.003
LO,cQq83,sigma_C_rr,quad,-0.003,0.002,0.001
LO,cQq83,sigma_C_kk,lin,-0.036,0.004,0.005
LO,cQq83,sigma_C_kk,quad,-0.055,0.003,0.003
LO,cQq83,sigma_C_rk+C_kr,lin,0.001,0.007,0.004
LO,cQq83,sigma_C_rk+C_kr,quad,-0.006,0.002,0.002
LO,cQq83,sigma_C_nr+C_rn,lin,0.006,0.002,0.001
LO,cQq83,sigma_C_nr+C_rn,quad,0.005,0.002,0.001
LO,cQq83,sigma_C_nk+C_kn,lin,-0.008,0.002,0.002
LO,cQq83,sigma_C_nk+C_kn,quad,-0.002,0.0,0.0
LO,cQq13,sigma,lin,0.006,0.002,0.001
LO,cQq13,sigma,quad,0.334,0.019,0.017
LO,cQq13,sigma_B_n+barB_n,lin,-0.002,0.002,0.001
LO,cQq13,sigma_B_n+barB_n,quad,-0.006,0.001,0.001
LO,cQq13,sigma_B_r+barB_r,lin,-0.002,0.003,0.005
LO,cQq13,sigma_B_r+barB_r,quad,-0.05,0.003,0.003
LO,cQq13,sigma_B_k+barB_k,lin,0.036,0.007,0.008
LO,cQq13,sigma_B_k+barB_k,quad,-0.584,0.038,0.043
LO,cQq13,sigma_C_nn,lin,0.011,0.001,0.001
LO,cQq13,sigma_C_nn,quad,-0.006,0.002,0.002
LO,cQq13,sigma_C_rr,lin,-0.001,0.001,0.001
LO,cQq13,sigma_C_rr,quad,-0.029,0.002,0.001
LO,cQq13,sigma_C_kk,lin,0.085,0.005,0.004
LO,cQq13,sigma_C_kk,quad,-0.246,0.013,0.015
LO,cQq13,sigma_C_rk+C_kr,lin,0.052,0.011,0.008
LO,cQq13,sigma_C_rk+C_kr,quad,-0.059,0.003,0.002
LO,cQq13,sigma_C_nr+C_rn,lin,-0.009,0.001,0.002
LO,cQq13,sigma_C_nr+C_rn,quad,0.003,0.001,0.001
LO,cQq13,sigma_C_nk+C_kn,lin,-0.011,0.003,0.004
LO,cQq13,sigma_C_nk+C_kn,quad,-0.003,0.0,0.0
LO,cQq81,sigma,lin,0.385,0.046,0.039
LO,cQq81,sigma,quad,0.075,0.004,0.004
LO,cQq81,sigma_B_n+barB_n,lin,-0.013,0.002,0.003
LO,cQq81,sigma_B_n+barB_n,quad,0.001,0.001,0.0
LO,cQq81,sigma_B_r+barB_r,lin,-0.062,0.004,0.005
LO,cQq81,sigma_B_r+barB_r,quad,-0.012,0.001,0.001
LO,cQq81,sigma_B_k+barB_k,lin,-0.441,0.054,0.068
LO,cQq81,sigma_B_k+barB_k,quad,-0.123,0.009,0.011
LO,cQq81,sigma_C_nn,lin,0.089,0.01,0.009
LO,cQq81,sigma_C_nn,quad,0.002,0.0,0.001
LO,cQq81,sigma_C_rr,lin,-0.201,0.019,0.022
LO,cQq81,sigma_C_rr,quad,-0.003,0.002,0.001
LO,cQq81,sigma_C_kk,lin,-0.202,0.02,0.022
LO,cQq81,sigma_C_kk,quad,-0.047,0.003,0.003
LO,cQq81,sigma_C_rk+C_kr,lin,-0.194,0.022,0.028
LO,cQq81,sigma_C_rk+C_kr,quad,-0.011,0.001,0.001
LO,cQq81,sigma_C_nr+C_rn,lin,0.024,0.007,0.005
LO,cQq81,sigma_C_nr+C_rn,quad,-0.002,0.001,0.001
LO,cQq81,sigma_C_nk+C_kn,lin,0.002,0.002,0.001
LO,cQq81,sigma_C_nk+C_kn,quad,0.001,0.0,0.0
LO,cQq11,sigma,lin,-0.002,0.0,0.0
LO,cQq11,sigma,quad,0.335,0.019,0.017
LO,cQq11,sigma_B_n+barB_n,lin,-0.006,0.001,0.001
LO,cQq11,sigma_B_n+barB_n,quad,0.001,0.0,0.0
LO,cQq11,sigma_B_r+barB_r,lin,0.009,0.004,0.002
LO,cQq11,sigma_B_r+barB_r,quad,-0.054,0.002,0.002
LO,cQq11,sigma_B_k+barB_k,lin,-0.004,0.0,0.001
LO,cQq11,sigma_B_k+barB_k,quad,-0.557,0.035,0.04
LO,cQq11,sigma_C_nn,lin,0.004,0.003,0.002
LO,cQq11,sigma_C_nn,quad,-0.002,0.001,0.002
LO,cQq11,sigma_C_rr,lin,0.019,0.006,0.004
LO,cQq11,sigma_C_rr,quad,-0.026,0.003,0.002
LO,cQq11,sigma_C_kk,lin,0.004,0.001,0.001
LO,cQq11,sigma_C_kk,quad,-0.219,0.012,0.013
LO,cQq11,sigma_C_rk+C_kr,lin,0.007,0.003,0.002
LO,cQq11,sigma_C_rk+C_kr,quad,-0.046,0.005,0.005
LO,cQq11,sigma_C_nr+C_rn,lin,-0.025,0.004,0.006
LO,cQq11,sigma_C_nr+C_rn,quad,0.003,0.0,0.0
LO,cQq11,sigma_C_nk+C_kn,lin,-0.01,0.002,0.003
LO,cQq11,sigma_C_nk+C_kn,quad,0.009,0.002,0.002
LO,cQu8,sigma,lin,0.23,0.027,0.023
LO,cQu8,sigma,quad,0.046,0.003,0.002
LO,cQu8,sigma_B_n+barB_n,lin,-0.003,0.001,0.002
LO,cQu8,sigma_B_n+barB_n,quad,0.003,0.001,0.001
LO,cQu8,sigma_B_r+barB_r,lin,-0.033,0.002,0.003
LO,cQu8,sigma_B_r+barB_r,quad,-0.008,0.001,0.001
LO,cQu8,sigma_B_k+barB_k,lin,-0.272,0.034,0.043
LO,cQu8,sigma_B_k+barB_k,quad,-0.08,0.006,0.008
LO,cQu8,sigma_C_nn,lin,0.056,0.007,0.006
LO,cQu8,sigma_C_nn,quad,0.0,0.0,0.001
LO,cQu8,sigma_C_rr,lin,-0.123,0.011,0.013
LO,cQu8,sigma_C_rr,quad,-0.006,0.0,0.0
LO,cQu8,sigma_C_kk,lin,-0.14,0.016,0.019
LO,cQu8,sigma_C_kk,quad,-0.033,0.002,0.002
LO,cQu8,sigma_C_rk+C_kr,lin,-0.116,0.012,0.015
LO,cQu8,sigma_C_rk+C_kr,quad,-0.007,0.001,0.001
LO,cQu8,sigma_C_nr+C_rn,lin,0.014,0.007,0.004
LO,cQu8,sigma_C_nr+C_rn,quad,0.001,0.0,0.0
LO,cQu8,sigma_C_nk+C_kn,lin,0.006,0.003,0.002
LO,cQu8,sigma_C_nk+C_kn,quad,-0.002,0.001,0.001
LO,cQu1,sigma,lin,0.001,0.0,0.0
LO,cQu1,sigma,quad,0.207,0.012,0.011
LO,cQu1,sigma_B_n+barB_n,lin,0.006,0.001,0.001
LO,cQu1,sigma_B_n+barB_n,quad,-0.001,0.0,0.0
LO,cQu1,sigma_B_r+barB_r,lin,-0.0,0.001,0.001
LO,cQu1,sigma_B_r+barB_r,quad,-0.03,0.001,0.001
LO,cQu1,sigma_B_k+barB_k,lin,0.019,0.002,0.002
LO,cQu1,sigma_B_k+barB_k,quad,-0.357,0.023,0.026
LO,cQu1,sigma_C_nn,lin,-0.006,0.001,0.0
LO,cQu1,sigma_C_nn,quad,-0.002,0.001,0.001
LO,cQu1,sigma_C_rr,lin,-0.004,0.002,0.002
LO,cQu1,sigma_C_rr,quad,-0.016,0.001,0.001
LO,cQu1,sigma_C_kk,lin,0.048,0.006,0.004
LO,cQu1,sigma_C_kk,quad,-0.147,0.008,0.009
LO,cQu1,sigma_C_rk+C_kr,lin,0.002,0.002,0.002
LO,cQu1,sigma_C_rk+C_kr,quad,-0.034,0.002,0.001
LO,cQu1,sigma_C_nr+C_rn,lin,0.029,0.008,0.006
LO,cQu1,sigma_C_nr+C_rn,quad,0.001,0.001,0.001
LO,cQu1,sigma_C_nk+C_kn,lin,0.017,0.005,0.003
LO,cQu1,sigma_C_nk+C_kn,quad,0.004,0.002,0.001
LO,ctq8,sigma,lin,0.38,0.045,0.038
LO,ctq8,sigma,quad,0.074,0.004,0.004
LO,ctq8,sigma_B_n+barB_n,lin,-0.002,0.0,0.001
LO,ctq8,sigma_B_n+barB_n,quad,-0.007,0.001,0.002
LO,ctq8,sigma_B_r+barB_r,lin,0.049,0.003,0.002
LO,ctq8,sigma_B_r+barB_r,quad,0.008,0.0,0.0
LO,ctq8,sigma_B_k+barB_k,lin,0.446,0.073,0.057
LO,ctq8,sigma_B_k+barB_k,quad,0.123,0.012,0.01
LO,ctq8,sigma_C_nn,lin,0.099,0.012,0.01
LO,ctq8,sigma_C_nn,quad,0.001,0.001,0.001
LO,ctq8,sigma_C_rr,lin,-0.191,0.016,0.019
LO,ctq8,sigma_C_rr,quad,-0.0,0.002,0.002
LO,ctq8,sigma_C_kk,lin,-0.235,0.027,0.034
LO,ctq8,sigma_C_kk,quad,-0.045,0.003,0.003
LO,ctq8,sigma_C_rk+C_kr,lin,-0.167,0.015,0.018
LO,ctq8,sigma_C_rk+C_kr,quad,-0.013,0.001,0.001
LO,ctq8,sigma_C_nr+C_rn,lin,0.016,0.004,0.003
LO,ctq8,sigma_C_nr+C_rn,quad,0.01,0.003,0.002
LO,ctq8,sigma_C_nk+C_kn,lin,-0.016,0.003,0.004
LO,ctq8,sigma_C_nk+C_kn,quad,0.0,0.0,0.0
LO,ctq1,sigma,lin,-0.001,0.0,0.0
LO,ctq1,sigma,quad,0.335,0.019,0.017
LO,ctq1,sigma_B_n+barB_n,lin,-0.013,0.003,0.005
LO,ctq1,sigma_B_n+barB_n,quad,0.0,0.001,0.001
LO,ctq1,sigma_B_r+barB_r,lin,-0.006,0.001,0.001
LO,ctq1,sigma_B_r+barB_r,quad,0.05,0.002,0.002
LO,ctq1,sigma_B_k+barB_k,lin,0.012,0.002,0.002
LO,ctq1,sigma_B_k+barB_k,quad,0.562,0.042,0.036
LO,ctq1,sigma_C_nn,lin,0.007,0.005,0.003
LO,ctq1,sigma_C_nn,quad,-0.006,0.001,0.002
LO,ctq1,sigma_C_rr,lin,0.003,0.002,0.002
LO,ctq1,sigma_C_rr,quad,-0.022,0.003,0.002
LO,ctq1,sigma_C_kk,lin,-0.03,0.004,0.004
LO,ctq1,sigma_C_kk,quad,-0.217,0.012,0.013
LO,ctq1,sigma_C_rk+C_kr,lin,-0.011,0.001,0.002
LO,ctq1,sigma_C_rk+C_kr,quad,-0.051,0.003,0.003
LO,ctq1,sigma_C_nr+C_rn,lin,0.011,0.003,0.002
LO,ctq1,sigma_C_nr+C_rn,quad,-0.006,0.001,0.001
LO,ctq1,sigma_C_nk+C_kn,lin,0.028,0.011,0.007
LO,ctq1,sigma_C_nk+C_kn,quad,-0.004,0.001,0.001
LO,cQd8,sigma,lin,0.158,0.018,0.016
LO,cQd8,sigma,quad,0.028,0.002,0.002
LO,cQd8,sigma_B_n+barB_n,lin,-0.007,0.002,0.002
LO,cQd8,sigma_B_n+barB_n,quad,-0.003,0.001,0.001
LO,cQd8,sigma_B_r+barB_r,lin,-0.02,0.001,0.002
LO,cQd8,sigma_B_r+barB_r,quad,-0.005,0.001,0.001
LO,cQd8,sigma_B_k+barB_k,lin,-0.195,0.026,0.034
LO,cQd8,sigma_B_k+barB_k,quad,-0.049,0.004,0.005
LO,cQd8,sigma_C_nn,lin,0.05,0.009,0.007
LO,cQd8,sigma_C_nn,quad,0.008,0.002,0.002
LO,cQd8,sigma_C_rr,lin,-0.083,0.007,0.008
LO,cQd8,sigma_C_rr,quad,-0.006,0.001,0.001
LO,cQd8,sigma_C_kk,lin,-0.106,0.013,0.017
LO,cQd8,sigma_C_kk,quad,-0.022,0.002,0.002
LO,cQd8,sigma_C_rk+C_kr,lin,-0.06,0.005,0.005
LO,cQd8,sigma_C_rk+C_kr,quad,-0.004,0.0,0.0
LO,cQd8,sigma_C_nr+C_rn,lin,0.005,0.002,0.002
LO,cQd8,sigma_C_nr+C_rn,quad,0.002,0.001,0.0
LO,cQd8,sigma_C_nk+C_kn,lin,0.012,0.005,0.004
LO,cQd8,sigma_C_nk+C_kn,quad,0.003,0.001,0.001
LO,cQd1,sigma,lin,0.001,0.0,0.0
LO,cQd1,sigma,quad,0.128,0.007,0.007
LO,cQd1,sigma_B_n+barB_n,lin,-0.0,0.001,0.001
LO,cQd1,sigma_B_n+barB_n,quad,-0.001,0.0,0.0
LO,cQd1,sigma_B_r+barB_r,lin,0.011,0.003,0.002
LO,cQd1,sigma_B_r+barB_r,quad,-0.02,0.001,0.001
LO,cQd1,sigma_B_k+barB_k,lin,0.003,0.002,0.002
LO,cQd1,sigma_B_k+barB_k,quad,-0.217,0.015,0.017
LO,cQd1,sigma_C_nn,lin,0.011,0.004,0.003
LO,cQd1,sigma_C_nn,quad,-0.0,0.0,0.0
LO,cQd1,sigma_C_rr,lin,-0.002,0.001,0.002
LO,cQd1,sigma_C_rr,quad,-0.011,0.001,0.001
LO,cQd1,sigma_C_kk,lin,0.007,0.002,0.001
LO,cQd1,sigma_C_kk,quad,-0.078,0.005,0.005
LO,cQd1,sigma_C_rk+C_kr,lin,-0.016,0.004,0.005
LO,cQd1,sigma_C_rk+C_kr,quad,-0.015,0.002,0.001
LO,cQd1,sigma_C_nr+C_rn,lin,0.006,0.002,0.002
LO,cQd1,sigma_C_nr+C_rn,quad,-0.008,0.001,0.001
LO,cQd1,sigma_C_nk+C_kn,lin,-0.002,0.002,0.003
LO,cQd1,sigma_C_nk+C_kn,quad,0.007,0.002,0.001
LO,ctu8,sigma,lin,0.22,0.025,0.021
LO,ctu8,sigma,quad,0.045,0.003,0.002
LO,ctu8,sigma_B_n+barB_n,lin,0.007,0.002,0.001
LO,ctu8,sigma_B_n+barB_n,quad,-0.002,0.0,0.001
LO,ctu8,sigma_B_r+barB_r,lin,0.031,0.001,0.001
LO,ctu8,sigma_B_r+barB_r,quad,0.009,0.001,0.001
LO,ctu8,sigma_B_k+barB_k,lin,0.262,0.04,0.032
LO,ctu8,sigma_B_k+barB_k,quad,0.081,0.008,0.007
LO,ctu8,sigma_C_nn,lin,0.036,0.001,0.001
LO,ctu8,sigma_C_nn,quad,0.0,0.0,0.0
LO,ctu8,sigma_C_rr,lin,-0.132,0.013,0.016
LO,ctu8,sigma_C_rr,quad,-0.001,0.001,0.001
LO,ctu8,sigma_C_kk,lin,-0.127,0.012,0.013
LO,ctu8,sigma_C_kk,quad,-0.032,0.002,0.002
LO,ctu8,sigma_C_rk+C_kr,lin,-0.105,0.01,0.012
LO,ctu8,sigma_C_rk+C_kr,quad,-0.008,0.0,0.0
LO,ctu8,sigma_C_nr+C_rn,lin,-0.017,0.004,0.006
LO,ctu8,sigma_C_nr+C_rn,quad,-0.006,0.002,0.002
LO,ctu8,sigma_C_nk+C_kn,lin,-0.006,0.001,0.002
LO,ctu8,sigma_C_nk+C_kn,quad,0.004,0.001,0.001
LO,ctu1,sigma,lin,0.002,0.0,0.0
LO,ctu1,sigma,quad,0.206,0.012,0.011
LO,ctu1,sigma_B_n+barB_n,lin,-0.005,0.0,0.001
LO,ctu1,sigma_B_n+barB_n,quad,-0.002,0.0,0.0
LO,ctu1,sigma_B_r+barB_r,lin,-0.006,0.002,0.003
LO,ctu1,sigma_B_r+barB_r,quad,0.032,0.001,0.001
LO,ctu1,sigma_B_k+barB_k,lin,0.026,0.004,0.003
LO,ctu1,sigma_B_k+barB_k,quad,0.36,0.026,0.023
LO,ctu1,sigma_C_nn,lin,0.012,0.002,0.002
LO,ctu1,sigma_C_nn,quad,-0.004,0.001,0.002
LO,ctu1,sigma_C_rr,lin,-0.007,0.003,0.004
LO,ctu1,sigma_C_rr,quad,-0.019,0.001,0.001
LO,ctu1,sigma_C_kk,lin,-0.032,0.004,0.005
LO,ctu1,sigma_C_kk,quad,-0.151,0.008,0.01
LO,ctu1,sigma_C_rk+C_kr,lin,0.018,0.005,0.004
LO,ctu1,sigma_C_rk+C_kr,quad,-0.034,0.002,0.002
LO,ctu1,sigma_C_nr+C_rn,lin,-0.008,0.002,0.003
LO,ctu1,sigma_C_nr+C_rn,quad,0.001,0.0,0.0
LO,ctu1,sigma_C_nk+C_kn,lin,-0.02,0.005,0.007
LO,ctu1,sigma_C_nk+C_kn,quad,0.012,0.003,0.002
LO,ctd8,sigma,lin,0.161,0.02,0.017
LO,ctd8,sigma,quad,0.03,0.002,0.002
LO,ctd8,sigma_B_n+barB_n,lin,-0.005,0.001,0.001
LO,ctd8,sigma_B_n+barB_n,quad,-0.002,0.0,0.001
LO,ctd8,sigma_B_r+barB_r,lin,0.026,0.002,0.002
LO,ctd8,sigma_B_r+barB_r,quad,0.005,0.0,0.0
LO,ctd8,sigma_B_k+barB_k,lin,0.188,0.029,0.023
LO,ctd8,sigma_B_k+barB_k,quad,0.05,0.005,0.004
LO,ctd8,sigma_C_nn,lin,0.023,0.0,0.001
LO,ctd8,sigma_C_nn,quad,0.001,0.0,0.0
LO,ctd8,sigma_C_rr,lin,-0.073,0.006,0.007
LO,ctd8,sigma_C_rr,quad,-0.011,0.002,0.003
LO,ctd8,sigma_C_kk,lin,-0.095,0.009,0.01
LO,ctd8,sigma_C_kk,quad,-0.022,0.001,0.001
LO,ctd8,sigma_C_rk+C_kr,lin,-0.099,0.012,0.014
LO,ctd8,sigma_C_rk+C_kr,quad,-0.002,0.001,0.001
LO,ctd8,sigma_C_nr+C_rn,lin,0.006,0.002,0.001
LO,ctd8,sigma_C_nr+C_rn,quad,0.006,0.002,0.001
LO,ctd8,sigma_C_nk+C_kn,lin,0.011,0.004,0.003
LO,ctd8,sigma_C_nk+C_kn,quad,0.0,0.0,0.0
LO,ctd1,sigma,lin,0.001,0.0,0.0
LO,ctd1,sigma,quad,0.129,0.007,0.007
LO,ctd1,sigma_B_n+barB_n,lin,0.005,0.002,0.002
LO,ctd1,sigma_B_n+barB_n,quad,-0.002,0.001,0.001
LO,ctd1,sigma_B_r+barB_r,lin,0.006,0.002,0.002
LO,ctd1,sigma_B_r+barB_r,quad,0.021,0.001,0.001
LO,ctd1,sigma_B_k+barB_k,lin,-0.005,0.001,0.002
LO,ctd1,sigma_B_k+barB_k,quad,0.215,0.016,0.014
LO,ctd1,sigma_C_nn,lin,-0.009,0.002,0.002
LO,ctd1,sigma_C_nn,quad,0.001,0.0,0.0
LO,ctd1,sigma_C_rr,lin,0.008,0.005,0.003
LO,ctd1,sigma_C_rr,quad,-0.007,0.002,0.001
LO,ctd1,sigma_C_kk,lin,-0.028,0.005,0.007
LO,ctd1,sigma_C_kk,quad,-0.08,0.005,0.005
LO,ctd1,sigma_C_rk+C_kr,lin,0.014,0.005,0.003
LO,ctd1,sigma_C_rk+C_kr,quad,-0.025,0.001,0.001
LO,ctd1,sigma_C_nr+C_rn,lin,0.002,0.0,0.001
LO,ctd1,sigma_C_nr+C_rn,quad,-0.006,0.002,0.002
LO,ctd1,sigma_C_nk+C_kn,lin,0.02,0.005,0.004
LO,ctd1,sigma_C_nk+C_kn,quad,-0.009,0.002,0.003
NLO,SM,sigma,sm,36.23,1.7,3.14
NLO,SM,sigma_B_n+barB_n,sm,-0.04,0.0,0.01
NLO,SM,sigma_B_r+barB_r,sm,0.03,0.03,0.02
NLO,SM,sigma_B_k+barB_k,sm,0.04,0.03,0.02
NLO,SM,sigma_C_nn,sm,11.27,0.43,0.8
NLO,SM,sigma_C_rr,sm,1.66,0.65,0.45
NLO,SM,sigma_C_kk,sm,11.83,0.68,1.03
NLO,SM,sigma_C_rk+C_kr,sm,-7.36,0.44,0.13
NLO,SM,sigma_C_nr+C_rn,sm,-0.03,0.04,0.03
NLO,SM,sigma_C_nk+C_kn,sm,0.16,0.05,0.04
NLO,c_tG,sigma,lin,11.5,0.51,0.97
NLO,c_tG,sigma,quad,2.01,0.09,0.17
NLO,c_tG,sigma_B_n+barB_n,lin,0.0,0.01,0.0
NLO,c_tG,sigma_B_n+barB_n,quad,0.01,0.0,0.0
NLO,c_tG,sigma_B_r+barB_r,lin,-0.01,0.01,0.01
NLO,c_tG,sigma_B_r+barB_r,quad,-0.0,0.0,0.0
NLO,c_tG,sigma_B_k+barB_k,lin,0.03,0.01,0.01
NLO,c_tG,sigma_B_k+barB_k,quad,0.0,0.0,0.0
NLO,c_tG,sigma_C_nn,lin,8.3,0.42,0.72
NLO,c_tG,sigma_C_nn,quad,0.56,0.05,0.06
NLO,c_tG,sigma_C_rr,lin,5.75,0.37,0.57
NLO,c_tG,sigma_C_rr,quad,0.81,0.04,0.07
NLO,c_tG,sigma_C_kk,lin,5.66,0.29,0.47
NLO,c_tG,sigma_C_kk,quad,1.17,0.05,0.1
NLO,c_tG,sigma_C_rk+C_kr,lin,-1.07,0.14,0.18
NLO,c_tG,sigma_C_rk+C_kr,quad,-0.04,0.0,0.0
NLO,c_tG,sigma_C_nr+C_rn,lin,-0.03,0.01,0.02
NLO,c_tG,sigma_C_nr+C_rn,quad,-0.01,0.0,0.0
NLO,c_tG,sigma_C_nk+C_kn,lin,0.01,0.03,0.01
NLO,c_tG,sigma_C_nk+C_kn,quad,0.0,0.0,0.0
NLO,cQq83,sigma,lin,0.059,0.002,0.005
NLO,cQq83,sigma,quad,0.072,0.002,0.001
NLO,cQq83,sigma_B_n+barB_n,lin,-0.009,0.004,0.007
NLO,cQq83,sigma_B_n+barB_n,quad,0.0,0.001,0.0
NLO,cQq83,sigma_B_r+barB_r,lin,-0.011,0.003,0.003
NLO,cQq83,sigma_B_r+barB_r,quad,-0.005,0.001,0.001
NLO,cQq83,sigma_B_k+barB_k,lin,-0.083,0.004,0.001
NLO,cQq83,sigma_B_k+barB_k,quad,-0.118,0.002,0.001
NLO,cQq83,sigma_C_nn,lin,-0.002,0.008,0.015
NLO,cQq83,sigma_C_nn,quad,-0.004,0.001,0.002
NLO,cQq83,sigma_C_rr,lin,-0.027,0.008,0.005
NLO,cQq83,sigma_C_rr,quad,-0.004,0.0,0.001
NLO,cQq83,sigma_C_kk,lin,-0.08,0.012,0.017
NLO,cQq83,sigma_C_kk,quad,-0.041,0.002,0.001
NLO,cQq83,sigma_C_rk+C_kr,lin,-0.005,0.02,0.01
NLO,cQq83,sigma_C_rk+C_kr,quad,0.001,0.002,0.002
NLO,cQq83,sigma_C_nr+C_rn,lin,-0.013,0.004,0.007
NLO,cQq83,sigma_C_nr+C_rn,quad,0.001,0.0,0.0
NLO,cQq83,sigma_C_nk+C_kn,lin,-0.009,0.007,0.013
NLO,cQq83,sigma_C_nk+C_kn,quad,-0.002,0.001,0.001
NLO,cQq13,sigma,lin,-0.002,0.001,0.0
NLO,cQq13,sigma,quad,0.476,0.025,0.022
NLO,cQq13,sigma_B_n+barB_n,lin,0.022,0.01,0.006
NLO,cQq13,sigma_B_n+barB_n,quad,0.003,0.001,0.001
NLO,cQq13,sigma_B_r+barB_r,lin,-0.005,0.006,0.004
NLO,cQq13,sigma_B_r+barB_r,quad,-0.059,0.001,0.002
NLO,cQq13,sigma_B_k+barB_k,lin,-0.004,0.007,0.004
NLO,cQq13,sigma_B_k+barB_k,quad,-0.761,0.031,0.034
NLO,cQq13,sigma_C_nn,lin,0.006,0.002,0.003
NLO,cQq13,sigma_C_nn,quad,0.006,0.001,0.001
NLO,cQq13,sigma_C_rr,lin,-0.01,0.003,0.006
NLO,cQq13,sigma_C_rr,quad,-0.04,0.001,0.002
NLO,cQq13,sigma_C_kk,lin,-0.002,0.015,0.007
NLO,cQq13,sigma_C_kk,quad,-0.294,0.008,0.008
NLO,cQq13,sigma_C_rk+C_kr,lin,-0.015,0.007,0.013
NLO,cQq13,sigma_C_rk+C_kr,quad,-0.076,0.005,0.008
NLO,cQq13,sigma_C_nr+C_rn,lin,0.013,0.009,0.004
NLO,cQq13,sigma_C_nr+C_rn,quad,0.004,0.003,0.002
NLO,cQq13,sigma_C_nk+C_kn,lin,0.026,0.014,0.009
NLO,cQq13,sigma_C_nk+C_kn,quad,-0.004,0.001,0.002
NLO,cQq81,sigma,lin,0.345,0.003,0.013
NLO,cQq81,sigma,quad,0.078,0.002,0.001
NLO,cQq81,sigma_B_n+barB_n,lin,0.005,0.003,0.001
NLO,cQq81,sigma_B_n+barB_n,quad,0.001,0.001,0.001
NLO,cQq81,sigma_B_r+barB_r,lin,-0.074,0.004,0.004
NLO,cQq81,sigma_B_r+barB_r,quad,-0.006,0.002,0.001
NLO,cQq81,sigma_B_k+barB_k,lin,-0.469,0.014,0.003
NLO,cQq81,sigma_B_k+barB_k,quad,-0.129,0.003,0.003
NLO,cQq81,sigma_C_nn,lin,0.073,0.007,0.013
NLO,cQq81,sigma_C_nn,quad,-0.003,0.002,0.003
NLO,cQq81,sigma_C_rr,lin,-0.206,0.003,0.001
NLO,cQq81,sigma_C_rr,quad,-0.004,0.0,0.001
NLO,cQq81,sigma_C_kk,lin,-0.286,0.022,0.03
NLO,cQq81,sigma_C_kk,quad,-0.047,0.001,0.002
NLO,cQq81,sigma_C_rk+C_kr,lin,-0.167,0.009,0.005
NLO,cQq81,sigma_C_rk+C_kr,quad,-0.002,0.003,0.002
NLO,cQq81,sigma_C_nr+C_rn,lin,-0.012,0.002,0.003
NLO,cQq81,sigma_C_nr+C_rn,quad,-0.001,0.001,0.001
NLO,cQq81,sigma_C_nk+C_kn,lin,0.005,0.01,0.005
NLO,cQq81,sigma_C_nk+C_kn,quad,-0.001,0.001,0.002
NLO,cQq11,sigma,lin,-0.019,0.004,0.006
NLO,cQq11,sigma,quad,0.478,0.026,0.022
NLO,cQq11,sigma_B_n+barB_n,lin,-0.005,0.003,0.005
NLO,cQq11,sigma_B_n+barB_n,quad,-0.002,0.001,0.002
NLO,cQq11,sigma_B_r+barB_r,lin,0.031,0.008,0.006
NLO,cQq11,sigma_B_r+barB_r,quad,-0.055,0.001,0.001
NLO,cQq11,sigma_B_k+barB_k,lin,-0.044,0.011,0.015
NLO,cQq11,sigma_B_k+barB_k,quad,-0.767,0.032,0.035
NLO,cQq11,sigma_C_nn,lin,-0.03,0.007,0.01
NLO,cQq11,sigma_C_nn,quad,0.004,0.001,0.002
NLO,cQq11,sigma_C_rr,lin,0.017,0.01,0.006
NLO,cQq11,sigma_C_rr,quad,-0.032,0.001,0.002
NLO,cQq11,sigma_C_kk,lin,-0.026,0.01,0.015
NLO,cQq11,sigma_C_kk,quad,-0.305,0.01,0.012
NLO,cQq11,sigma_C_rk+C_kr,lin,0.04,0.02,0.012
NLO,cQq11,sigma_C_rk+C_kr,quad,-0.065,0.002,0.003
NLO,cQq11,sigma_C_nr+C_rn,lin,0.012,0.009,0.006
NLO,cQq11,sigma_C_nr+C_rn,quad,-0.005,0.001,0.001
NLO,cQq11,sigma_C_nk+C_kn,lin,-0.001,0.008,0.013
NLO,cQq11,sigma_C_nk+C_kn,quad,0.003,0.001,0.0
NLO,cQu8,sigma,lin,0.193,0.005,0.012
NLO,cQu8,sigma,quad,0.034,0.001,0.002
NLO,cQu8,sigma_B_n+barB_n,lin,0.003,0.004,0.005
NLO,cQu8,sigma_B_n+barB_n,quad,-0.002,0.001,0.002
NLO,cQu8,sigma_B_r+barB_r,lin,0.002,0.015,0.009
NLO,cQu8,sigma_B_r+barB_r,quad,-0.003,0.001,0.001
NLO,cQu8,sigma_B_k+barB_k,lin,-0.203,0.021,0.009
NLO,cQu8,sigma_B_k+barB_k,quad,-0.055,0.004,0.003
NLO,cQu8,sigma_C_nn,lin,0.055,0.002,0.003
NLO,cQu8,sigma_C_nn,quad,0.001,0.0,0.001
NLO,cQu8,sigma_C_rr,lin,-0.111,0.011,0.004
NLO,cQu8,sigma_C_rr,quad,-0.0,0.0,0.0
NLO,cQu8,sigma_C_kk,lin,-0.138,0.004,0.004
NLO,cQu8,sigma_C_kk,quad,-0.015,0.002,0.002
NLO,cQu8,sigma_C_rk+C_kr,lin,-0.09,0.009,0.005
NLO,cQu8,sigma_C_rk+C_kr,quad,0.002,0.004,0.003
NLO,cQu8,sigma_C_nr+C_rn,lin,0.012,0.004,0.003
NLO,cQu8,sigma_C_nr+C_rn,quad,0.002,0.0,0.001
NLO,cQu8,sigma_C_nk+C_kn,lin,-0.01,0.002,0.002
NLO,cQu8,sigma_C_nk+C_kn,quad,0.001,0.002,0.001
NLO,cQu1,sigma,lin,-0.02,0.004,0.006
NLO,cQu1,sigma,quad,0.296,0.016,0.013
NLO,cQu1,sigma_B_n+barB_n,lin,-0.004,0.003,0.004
NLO,cQu1,sigma_B_n+barB_n,quad,-0.003,0.001,0.001
NLO,cQu1,sigma_B_r+barB_r,lin,0.011,0.009,0.006
NLO,cQu1,sigma_B_r+barB_r,quad,-0.027,0.001,0.001
NLO,cQu1,sigma_B_k+barB_k,lin,0.062,0.017,0.012
NLO,cQu1,sigma_B_k+barB_k,quad,-0.462,0.018,0.02
NLO,cQu1,sigma_C_nn,lin,0.014,0.003,0.002
NLO,cQu1,sigma_C_nn,quad,-0.001,0.001,0.001
NLO,cQu1,sigma_C_rr,lin,0.001,0.001,0.001
NLO,cQu1,sigma_C_rr,quad,-0.021,0.002,0.004
NLO,cQu1,sigma_C_kk,lin,0.044,0.011,0.009
NLO,cQu1,sigma_C_kk,quad,-0.179,0.005,0.005
NLO,cQu1,sigma_C_rk+C_kr,lin,0.044,0.015,0.009
NLO,cQu1,sigma_C_rk+C_kr,quad,-0.037,0.002,0.003
NLO,cQu1,sigma_C_nr+C_rn,lin,-0.001,0.003,0.004
NLO,cQu1,sigma_C_nr+C_rn,quad,-0.003,0.001,0.001
NLO,cQu1,sigma_C_nk+C_kn,lin,-0.003,0.009,0.015
NLO,cQu1,sigma_C_nk+C_kn,quad,-0.006,0.001,0.002
NLO,ctq8,sigma,lin,0.34,0.005,0.015
NLO,ctq8,sigma,quad,0.056,0.002,0.003
NLO,ctq8,sigma_B_n+barB_n,lin,0.011,0.007,0.004
NLO,ctq8,sigma_B_n+barB_n,quad,0.001,0.001,0.001
NLO,ctq8,sigma_B_r+barB_r,lin,0.014,0.01,0.015
NLO,ctq8,sigma_B_r+barB_r,quad,-0.001,0.002,0.003
NLO,ctq8,sigma_B_k+barB_k,lin,0.332,0.016,0.043
NLO,ctq8,sigma_B_k+barB_k,quad,0.089,0.004,0.006
NLO,ctq8,sigma_C_nn,lin,0.065,0.009,0.02
NLO,ctq8,sigma_C_nn,quad,-0.005,0.002,0.003
NLO,ctq8,sigma_C_rr,lin,-0.206,0.003,0.002
NLO,ctq8,sigma_C_rr,quad,-0.005,0.001,0.002
NLO,ctq8,sigma_C_kk,lin,-0.232,0.012,0.017
NLO,ctq8,sigma_C_kk,quad,-0.027,0.004,0.003
NLO,ctq8,sigma_C_rk+C_kr,lin,-0.166,0.007,0.01
NLO,ctq8,sigma_C_rk+C_kr,quad,-0.009,0.001,0.003
NLO,ctq8,sigma_C_nr+C_rn,lin,0.002,0.004,0.009
NLO,ctq8,sigma_C_nr+C_rn,quad,0.0,0.0,0.001
NLO,ctq8,sigma_C_nk+C_kn,lin,0.003,0.005,0.009
NLO,ctq8,sigma_C_nk+C_kn,quad,0.001,0.001,0.001
NLO,ctq1,sigma,lin,-0.03,0.005,0.008
NLO,ctq1,sigma,quad,0.475,0.025,0.022
NLO,ctq1,sigma_B_n+barB_n,lin,0.009,0.002,0.001
NLO,ctq1,sigma_B_n+barB_n,quad,0.002,0.001,0.001
NLO,ctq1,sigma_B_r+barB_r,lin,0.032,0.005,0.003
NLO,ctq1,sigma_B_r+barB_r,quad,0.043,0.003,0.004
NLO,ctq1,sigma_B_k+barB_k,lin,-0.062,0.014,0.019
NLO,ctq1,sigma_B_k+barB_k,quad,0.747,0.031,0.029
NLO,ctq1,sigma_C_nn,lin,0.009,0.004,0.003
NLO,ctq1,sigma_C_nn,quad,0.004,0.001,0.001
NLO,ctq1,sigma_C_rr,lin,-0.008,0.006,0.007
NLO,ctq1,sigma_C_rr,quad,-0.042,0.002,0.005
NLO,ctq1,sigma_C_kk,lin,0.026,0.016,0.014
NLO,ctq1,sigma_C_kk,quad,-0.295,0.01,0.011
NLO,ctq1,sigma_C_rk+C_kr,lin,-0.013,0.006,0.005
NLO,ctq1,sigma_C_rk+C_kr,quad,-0.07,0.004,0.006
NLO,ctq1,sigma_C_nr+C_rn,lin,-0.011,0.006,0.007
NLO,ctq1,sigma_C_nr+C_rn,quad,0.001,0.002,0.002
NLO,ctq1,sigma_C_nk+C_kn,lin,-0.002,0.013,0.006
NLO,ctq1,sigma_C_nk+C_kn,quad,-0.003,0.004,0.005
NLO,cQd8,sigma,lin,0.148,0.001,0.004
NLO,cQd8,sigma,quad,0.023,0.001,0.0
NLO,cQd8,sigma_B_n+barB_n,lin,0.005,0.006,0.003
NLO,cQd8,sigma_B_n+barB_n,quad,-0.001,0.001,0.001
NLO,cQd8,sigma_B_r+barB_r,lin,0.003,0.006,0.004
NLO,cQd8,sigma_B_r+barB_r,quad,-0.001,0.001,0.001
NLO,cQd8,sigma_B_k+barB_k,lin,-0.15,0.017,0.008
NLO,cQd8,sigma_B_k+barB_k,quad,-0.033,0.003,0.002
NLO,cQd8,sigma_C_nn,lin,0.022,0.005,0.008
NLO,cQd8,sigma_C_nn,quad,0.0,0.0,0.001
NLO,cQd8,sigma_C_rr,lin,-0.103,0.008,0.013
NLO,cQd8,sigma_C_rr,quad,-0.0,0.001,0.001
NLO,cQd8,sigma_C_kk,lin,-0.095,0.002,0.001
NLO,cQd8,sigma_C_kk,quad,-0.008,0.002,0.002
NLO,cQd8,sigma_C_rk+C_kr,lin,-0.078,0.005,0.007
NLO,cQd8,sigma_C_rk+C_kr,quad,-0.004,0.002,0.002
NLO,cQd8,sigma_C_nr+C_rn,lin,0.001,0.002,0.001
NLO,cQd8,sigma_C_nr+C_rn,quad,-0.0,0.001,0.001
NLO,cQd8,sigma_C_nk+C_kn,lin,0.002,0.003,0.002
NLO,cQd8,sigma_C_nk+C_kn,quad,0.005,0.002,0.001
NLO,cQd1,sigma,lin,0.002,0.001,0.001
NLO,cQd1,sigma,quad,0.182,0.009,0.008
NLO,cQd1,sigma_B_n+barB_n,lin,-0.003,0.004,0.007
NLO,cQd1,sigma_B_n+barB_n,quad,0.003,0.001,0.001
NLO,cQd1,sigma_B_r+barB_r,lin,-0.008,0.004,0.006
NLO,cQd1,sigma_B_r+barB_r,quad,-0.022,0.001,0.001
NLO,cQd1,sigma_B_k+barB_k,lin,0.036,0.008,0.007
NLO,cQd1,sigma_B_k+barB_k,quad,-0.282,0.011,0.012
NLO,cQd1,sigma_C_nn,lin,0.006,0.006,0.007
NLO,cQd1,sigma_C_nn,quad,0.005,0.001,0.001
NLO,cQd1,sigma_C_rr,lin,-0.006,0.002,0.002
NLO,cQd1,sigma_C_rr,quad,-0.016,0.001,0.002
NLO,cQd1,sigma_C_kk,lin,0.029,0.007,0.006
NLO,cQd1,sigma_C_kk,quad,-0.106,0.004,0.005
NLO,cQd1,sigma_C_rk+C_kr,lin,0.009,0.008,0.006
NLO,cQd1,sigma_C_rk+C_kr,quad,-0.022,0.001,0.001
NLO,cQd1,sigma_C_nr+C_rn,lin,-0.01,0.004,0.005
NLO,cQd1,sigma_C_nr+C_rn,quad,-0.002,0.001,0.002
NLO,cQd1,sigma_C_nk+C_kn,lin,-0.002,0.003,0.003
NLO,cQd1,sigma_C_nk+C_kn,quad,0.005,0.002,0.001
NLO,ctu8,sigma,lin,0.205,0.002,0.008
NLO,ctu8,sigma,quad,0.047,0.001,0.001
NLO,ctu8,sigma_B_n+barB_n,lin,0.0,0.003,0.002
NLO,ctu8,sigma_B_n+barB_n,quad,-0.002,0.001,0.001
NLO,ctu8,sigma_B_r+barB_r,lin,0.029,0.004,0.005
NLO,ctu8,sigma_B_r+barB_r,quad,0.006,0.002,0.001
NLO,ctu8,sigma_B_k+barB_k,lin,0.268,0.002,0.008
NLO,ctu8,sigma_B_k+barB_k,quad,0.078,0.002,0.001
NLO,ctu8,sigma_C_nn,lin,0.042,0.005,0.01
NLO,ctu8,sigma_C_nn,quad,-0.002,0.001,0.001
NLO,ctu8,sigma_C_rr,lin,-0.125,0.005,0.006
NLO,ctu8,sigma_C_rr,quad,-0.001,0.002,0.001
NLO,ctu8,sigma_C_kk,lin,-0.155,0.01,0.014
NLO,ctu8,sigma_C_kk,quad,-0.027,0.001,0.0
NLO,ctu8,sigma_C_rk+C_kr,lin,-0.112,0.008,0.012
NLO,ctu8,sigma_C_rk+C_kr,quad,-0.009,0.001,0.002
NLO,ctu8,sigma_C_nr+C_rn,lin,0.011,0.015,0.009
NLO,ctu8,sigma_C_nr+C_rn,quad,-0.002,0.001,0.002
NLO,ctu8,sigma_C_nk+C_kn,lin,-0.001,0.012,0.007
NLO,ctu8,sigma_C_nk+C_kn,quad,0.008,0.006,0.003
NLO,ctu1,sigma,lin,-0.0,0.001,0.001
NLO,ctu1,sigma,quad,0.296,0.015,0.013
NLO,ctu1,sigma_B_n+barB_n,lin,0.012,0.005,0.003
NLO,ctu1,sigma_B_n+barB_n,quad,-0.004,0.001,0.002
NLO,ctu1,sigma_B_r+barB_r,lin,-0.004,0.001,0.004
NLO,ctu1,sigma_B_r+barB_r,quad,0.032,0.0,0.001
NLO,ctu1,sigma_B_k+barB_k,lin,0.024,0.011,0.007
NLO,ctu1,sigma_B_k+barB_k,quad,0.471,0.023,0.02
NLO,ctu1,sigma_C_nn,lin,-0.02,0.005,0.007
NLO,ctu1,sigma_C_nn,quad,0.001,0.001,0.002
NLO,ctu1,sigma_C_rr,lin,0.017,0.006,0.003
NLO,ctu1,sigma_C_rr,quad,-0.021,0.001,0.002
NLO,ctu1,sigma_C_kk,lin,0.001,0.004,0.006
NLO,ctu1,sigma_C_kk,quad,-0.177,0.006,0.006
NLO,ctu1,sigma_C_rk+C_kr,lin,-0.011,0.013,0.023
NLO,ctu1,sigma_C_rk+C_kr,quad,-0.043,0.003,0.004
NLO,ctu1,sigma_C_nr+C_rn,lin,0.009,0.007,0.003
NLO,ctu1,sigma_C_nr+C_rn,quad,0.005,0.003,0.002
NLO,ctu1,sigma_C_nk+C_kn,lin,0.013,0.014,0.008
NLO,ctu1,sigma_C_nk+C_kn,quad,0.004,0.0,0.0
NLO,ctd8,sigma,lin,0.139,0.003,0.009
NLO,ctd8,sigma,quad,0.03,0.002,0.001
NLO,ctd8,sigma_B_n+barB_n,lin,0.006,0.001,0.001
NLO,ctd8,sigma_B_n+barB_n,quad,-0.001,0.0,0.0
NLO,ctd8,sigma_B_r+barB_r,lin,0.011,0.002,0.005
NLO,ctd8,sigma_B_r+barB_r,quad,0.002,0.001,0.001
NLO,ctd8,sigma_B_k+barB_k,lin,0.191,0.005,0.007
NLO,ctd8,sigma_B_k+barB_k,quad,0.044,0.001,0.0
NLO,ctd8,sigma_C_nn,lin,0.029,0.005,0.007
NLO,ctd8,sigma_C_nn,quad,0.001,0.001,0.001
NLO,ctd8,sigma_C_rr,lin,-0.082,0.003,0.001
NLO,ctd8,sigma_C_rr,quad,-0.001,0.0,0.0
NLO,ctd8,sigma_C_kk,lin,-0.122,0.013,0.019
NLO,ctd8,sigma_C_kk,quad,-0.017,0.0,0.001
NLO,ctd8,sigma_C_rk+C_kr,lin,-0.063,0.011,0.006
NLO,ctd8,sigma_C_rk+C_kr,quad,0.013,0.004,0.004
NLO,ctd8,sigma_C_nr+C_rn,lin,-0.019,0.009,0.015
NLO,ctd8,sigma_C_nr+C_rn,quad,-0.004,0.0,0.0
NLO,ctd8,sigma_C_nk+C_kn,lin,0.015,0.016,0.007
NLO,ctd8,sigma_C_nk+C_kn,quad,-0.002,0.003,0.006
NLO,ctd1,sigma,lin,-0.018,0.004,0.005
NLO,ctd1,sigma,quad,0.183,0.01,0.008
NLO,ctd1,sigma_B_n+barB_n,lin,-0.011,0.005,0.007
NLO,ctd1,sigma_B_n+barB_n,quad,-0.002,0.0,0.0
NLO,ctd1,sigma_B_r+barB_r,lin,0.02,0.007,0.005
NLO,ctd1,sigma_B_r+barB_r,quad,0.02,0.0,0.0
NLO,ctd1,sigma_B_k+barB_k,lin,0.02,0.005,0.004
NLO,ctd1,sigma_B_k+barB_k,quad,0.288,0.015,0.013
NLO,ctd1,sigma_C_nn,lin,-0.016,0.004,0.004
NLO,ctd1,sigma_C_nn,quad,-0.002,0.001,0.002
NLO,ctd1,sigma_C_rr,lin,-0.0,0.003,0.004
NLO,ctd1,sigma_C_rr,quad,-0.013,0.002,0.001
NLO,ctd1,sigma_C_kk,lin,-0.04,0.01,0.014
NLO,ctd1,sigma_C_kk,quad,-0.106,0.003,0.003
NLO,ctd1,sigma_C_rk+C_kr,lin,0.003,0.0,0.002
NLO,ctd1,sigma_C_rk+C_kr,quad,-0.026,0.001,0.002
NLO,ctd1,sigma_C_nr+C_rn,lin,-0.027,0.006,0.007
NLO,ctd1,sigma_C_nr+C_rn,quad,-0.001,0.001,0.002
NLO,ctd1,sigma_C_nk+C_kn,lin,-0.017,0.008,0.012
NLO,ctd1,sigma_C_nk+C_kn,quad,0.002,0.003,0.002
)THEORYCSV";

std::vector<std::string> observable_names() {
  return {
    "gen_b1k", "gen_b2k", "gen_b1r", "gen_b2r", "gen_b1n", "gen_b2n",
    "gen_b1j", "gen_b2j", "gen_b1q", "gen_b2q",
    "gen_c_kk", "gen_c_rr", "gen_c_nn",
    "gen_c_kj", "gen_c_rq",
    "gen_c_Prk", "gen_c_Mrk",
    "gen_c_Pnr", "gen_c_Mnr",
    "gen_c_Pnk", "gen_c_Mnk",
    "gen_c_Prj", "gen_c_Mrj",
    "gen_c_han", "gen_c_tra", "gen_c_sca",
    "gen_c_kjL", "gen_c_rqL",
    "gen_c_rkP", "gen_c_rkM",
    "gen_c_nrP", "gen_c_nrM",
    "gen_c_nkP", "gen_c_nkM",
    "gen_ll_cHel", "gen_ll_cLab",
    "gen_llbar_delta_phi", "gen_llbar_delta_eta"
  };
}

std::string obs_root_label(const std::string& obs) {
  if (obs == "gen_b1k") return "B_{1k}";
  if (obs == "gen_b2k") return "B_{2k}";
  if (obs == "gen_b1r") return "B_{1r}";
  if (obs == "gen_b2r") return "B_{2r}";
  if (obs == "gen_b1n") return "B_{1n}";
  if (obs == "gen_b2n") return "B_{2n}";
  if (obs == "gen_b1j") return "B_{1k^{*}}";
  if (obs == "gen_b2j") return "B_{2k^{*}}";
  if (obs == "gen_b1q") return "B_{1r^{*}}";
  if (obs == "gen_b2q") return "B_{2r^{*}}";
  if (obs == "gen_c_kk") return "C_{kk}";
  if (obs == "gen_c_rr") return "C_{rr}";
  if (obs == "gen_c_nn") return "C_{nn}";
  if (obs == "gen_c_kj") return "C_{kj}";
  if (obs == "gen_c_rq") return "C_{rq}";
  if (obs == "gen_c_Prk") return "C_{rk}^{+}";
  if (obs == "gen_c_Mrk") return "C_{rk}^{-}";
  if (obs == "gen_c_Pnr") return "C_{nr}^{+}";
  if (obs == "gen_c_Mnr") return "C_{nr}^{-}";
  if (obs == "gen_c_Pnk") return "C_{nk}^{+}";
  if (obs == "gen_c_Mnk") return "C_{nk}^{-}";
  if (obs == "gen_c_Prj") return "C_{rj}^{+}";
  if (obs == "gen_c_Mrj") return "C_{rj}^{-}";
  if (obs == "gen_c_han") return "C_{han}";
  if (obs == "gen_c_tra") return "C_{tra}";
  if (obs == "gen_c_sca") return "C_{sca}";
  if (obs == "gen_c_kjL") return "C_{kj}^{L}";
  if (obs == "gen_c_rqL") return "C_{rq}^{L}";
  if (obs == "gen_c_rkP") return "C_{rk}^{P}";
  if (obs == "gen_c_rkM") return "C_{rk}^{M}";
  if (obs == "gen_c_nrP") return "C_{nr}^{P}";
  if (obs == "gen_c_nrM") return "C_{nr}^{M}";
  if (obs == "gen_c_nkP") return "C_{nk}^{P}";
  if (obs == "gen_c_nkM") return "C_{nk}^{M}";
  if (obs == "gen_ll_cHel") return "cos#varphi_{ll}^{hel}";
  if (obs == "gen_ll_cLab") return "cos#varphi_{ll}^{lab}";
  if (obs == "gen_llbar_delta_phi") return "#Delta#phi_{ll}";
  if (obs == "gen_llbar_delta_eta") return "#Delta#eta_{ll}";
  return obs;
}
std::string theory_obs_key(const std::string& obs) {
  if (obs == "gen_b1n" || obs == "gen_b2n") return "sigma_B_n+barB_n";
  if (obs == "gen_b1r" || obs == "gen_b2r") return "sigma_B_r+barB_r";
  if (obs == "gen_b1k" || obs == "gen_b2k") return "sigma_B_k+barB_k";
  if (obs == "gen_c_nn") return "sigma_C_nn";
  if (obs == "gen_c_rr") return "sigma_C_rr";
  if (obs == "gen_c_kk") return "sigma_C_kk";
  if (obs == "gen_c_Prk" || obs == "gen_c_Mrk") return "sigma_C_rk+C_kr";
  if (obs == "gen_c_Pnr" || obs == "gen_c_Mnr") return "sigma_C_nr+C_rn";
  if (obs == "gen_c_Pnk" || obs == "gen_c_Mnk") return "sigma_C_nk+C_kn";
  return "";
}

std::string theory_op_from_wc(const std::string& wc) {
  if (wc == "ctGRe" || wc == "ctGIm") return "c_tG";
  if (wc == "cQj38") return "cQq83";
  if (wc == "cQj31") return "cQq13";
  if (wc == "cQj18") return "cQq81";
  if (wc == "cQj11") return "cQq11";
  if (wc == "cQu8") return "cQu8";
  if (wc == "cQu1") return "cQu1";
  if (wc == "cQd8") return "cQd8";
  if (wc == "cQd1") return "cQd1";
  if (wc == "ctu8") return "ctu8";
  if (wc == "ctu1") return "ctu1";
  if (wc == "ctd8") return "ctd8";
  if (wc == "ctd1") return "ctd1";
  if (wc == "ctj8") return "ctq8";
  if (wc == "ctj1") return "ctq1";
  return wc;
}

double asymmetry_factor(const std::string& obs) {
  if (obs.find("gen_b") == 0) return 2.0;
  if (obs == "gen_c_Prk" || obs == "gen_c_Mrk" || obs == "gen_c_Pnr" || obs == "gen_c_Mnr" ||
      obs == "gen_c_Pnk" || obs == "gen_c_Mnk" || obs == "gen_c_Prj" || obs == "gen_c_Mrj") return -16.0 / M_PI;
  if (obs == "gen_c_han" || obs == "gen_c_sca" || obs == "gen_c_tra" ||
      obs == "gen_c_kjL" || obs == "gen_c_rqL" || obs == "gen_c_rkP" || obs == "gen_c_rkM" ||
      obs == "gen_c_nrP" || obs == "gen_c_nrM" || obs == "gen_c_nkP" || obs == "gen_c_nkM") return -2.0;
  if (obs == "gen_ll_cHel" || obs == "gen_ll_cLab") return -4.0;
  if (obs.find("gen_c_") == 0) return -4.0;
  return 0.0;
}

std::map<std::string, HistChunk1D> split_values_into_chunks(const std::vector<double>& vals) {
  // Fallback helper for code paths that only have bin contents.
  // Do NOT assign sqrt(content) errors for normalized differential spectra:
  // those contents are not raw Poisson counts and this produced O(1-3)
  // fake uncertainties in the per-observable plots.
  std::map<std::string, HistChunk1D> out;
  const auto names = observable_names();
  const int nobs = std::min((int)names.size(), (int)vals.size() / BINS_PER_OBS);
  for (int iobs = 0; iobs < nobs; ++iobs) {
    HistChunk1D c;
    for (int ib = 0; ib < BINS_PER_OBS; ++ib) {
      const double y = vals[iobs * BINS_PER_OBS + ib];
      c.y.push_back(y);
      c.e.push_back(0.0);
    }
    out[names[iobs]] = c;
  }
  return out;
}

std::vector<double> obs_axis_edges(const std::string& obs);

std::map<std::string, HistChunk1D> load_hist_chunks_from_root(const std::string& path,
                                                       bool print_mc_error_debug = false,
                                                       const std::string& err_tag = "") {
  TFile f(path.c_str(), "READ");
  if (f.IsZombie()) throw std::runtime_error("Cannot open " + path);

  TH1* h = get_hist1(&f);
  if (!h) throw std::runtime_error("No TH1 found in " + path);

  std::map<std::string, HistChunk1D> out;
  const auto names = observable_names();
  const int nobs = std::min((int)names.size(), h->GetNbinsX() / BINS_PER_OBS);

  // IMPORTANT:
  // The concatenation script writes a hist.storage.Weight() histogram, i.e. each
  // flattened bin stores both "value" and "variance".  When ROOT reads this back
  // as a TH1, h->GetBinContent(gbin) is the normalized density value and
  // h->GetBinError(gbin) should be sqrt(saved variance).  Do NOT reconstruct the
  // MC errors from GetEntries(): after uproot/hist export, Entries is not a
  // reliable per-observable event count for this flattened density histogram.
  // If the printed errors below are not consistent with the generated-event count, then the
  // problem is upstream in concatenate_histograms.py or in the ROOT writing step.

  const double entries_total = h->GetEntries();
  const int sumw2_n = h->GetSumw2N();

  for (int iobs = 0; iobs < nobs; ++iobs) {
    const std::string& obs = names[iobs];
    HistChunk1D c;
    std::vector<double> edges = obs_axis_edges(obs);

    double integral_density = 0.0;
    double sum_var_density = 0.0;
    double min_err = 1e300;
    double max_err = 0.0;
    double neff_sum = 0.0;
    int neff_n = 0;

    for (int ib = 0; ib < BINS_PER_OBS; ++ib) {
      const int gbin = iobs * BINS_PER_OBS + ib + 1;
      const double y = h->GetBinContent(gbin);
      const double e = std::max(0.0, h->GetBinError(gbin));
      const double dx = (edges.size() == BINS_PER_OBS + 1) ? std::abs(edges[ib+1] - edges[ib]) : 1.0;

      c.y.push_back(y);
      c.e.push_back(e);

      integral_density += y * dx;
      sum_var_density += e * e;
      min_err = std::min(min_err, e);
      max_err = std::max(max_err, e);

      // Back-calculate an approximate multinomial N_eff from the saved density
      // error: p_i = y_i * dx, sigma(p_i) = e_i * dx,
      // sigma^2(p_i) ~= p_i(1-p_i)/N_eff.
      const double pbin = y * dx;
      const double ep = e * dx;
      if (pbin > 0.0 && pbin < 1.0 && ep > 0.0) {
        neff_sum += pbin * (1.0 - pbin) / (ep * ep);
        ++neff_n;
      }
    }

    if (print_mc_error_debug &&
        (obs == "gen_b1k" || obs == "gen_c_kk" || obs == "gen_ll_cHel")) {
      std::cout << std::scientific << std::setprecision(6)
                << "[EFTMC-ERR-READBACK] tag=" << err_tag
                << " obs=" << obs
                << " path=" << path
                << " entries_total=" << entries_total
                << " sumw2N=" << sumw2_n
                << " density_integral=" << integral_density
                << " min_bin_err=" << min_err
                << " max_bin_err=" << max_err
                << " rms_bin_err=" << std::sqrt(sum_var_density / std::max(1, BINS_PER_OBS));
      if (neff_n > 0) {
        std::cout << " approx_Neff_from_errors=" << (neff_sum / double(neff_n));
      }
      std::cout << std::endl;

      for (int ib = 0; ib < BINS_PER_OBS; ++ib) {
        const double dx = (edges.size() == BINS_PER_OBS + 1) ? std::abs(edges[ib+1] - edges[ib]) : 1.0;
        const double pbin = c.y[ib] * dx;
        const double ep = c.e[ib] * dx;
        double neff_bin = -1.0;
        if (pbin > 0.0 && pbin < 1.0 && ep > 0.0) {
          neff_bin = pbin * (1.0 - pbin) / (ep * ep);
        }
        std::cout << std::scientific << std::setprecision(6)
                  << "[EFTMC-ERR-BIN] tag=" << err_tag
                  << " obs=" << obs
                  << " ib=" << ib
                  << " y_density=" << c.y[ib]
                  << " err_density=" << c.e[ib]
                  << " p_width=" << pbin
                  << " err_p_width=" << ep
                  << " approx_Neff=" << neff_bin
                  << std::endl;
      }
    }

    out[obs] = c;
  }
  return out;
}

CoeffSummary coefficient_from_chunk(const std::string& obs, const HistChunk1D& c) {
  CoeffSummary s;
  const double f = asymmetry_factor(obs);
  if (f == 0.0 || c.y.size() != BINS_PER_OBS) return s;

  // Histograms are normalized differential densities.  For an asymmetry we
  // must integrate density over bin width first:
  //   p_i = density_i * width_i
  //   sigma(p_i) = sigma(density_i) * width_i
  // The old implementation summed raw density values directly.  For equal-width
  // bins the central AFB accidentally cancels the common width, but the stored
  // errors and any non-uniform binning should be propagated in probability space.
  std::vector<double> edges = obs_axis_edges(obs);
  double F = 0.0, B = 0.0, vF = 0.0, vB = 0.0;
  for (int i = 0; i < BINS_PER_OBS; ++i) {
    const double dx = (edges.size() == BINS_PER_OBS + 1) ? std::abs(edges[i+1] - edges[i]) : 1.0;
    const double pbin = c.y[i] * dx;
    const double epbin = c.e[i] * dx;
    if (i < BINS_PER_OBS / 2) {
      B += pbin;
      vB += epbin * epbin;
    } else {
      F += pbin;
      vF += epbin * epbin;
    }
  }

  if (std::abs(F + B) < 1e-15) return s;
  s.afb = (F - B) / (F + B);
  s.afb_err = (2.0 / ((F + B) * (F + B))) * std::sqrt(std::max(0.0, B * B * vF + F * F * vB));
  s.coeff = f * s.afb;
  s.coeff_err = std::abs(f) * s.afb_err;
  s.factor = f;
  s.ok = true;
  return s;
}

std::vector<std::string> csv_split_simple(const std::string& line) {
  std::vector<std::string> out; std::stringstream ss(line); std::string item;
  while (std::getline(ss, item, ',')) out.push_back(item);
  return out;
}

TheoryTable load_theory_csv_stream(std::istream& f) {
  TheoryTable tab;
  std::string line;
  std::getline(f, line); // header
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    auto c = csv_split_simple(line);
    if (c.size() < 6) continue;
    const std::string order = c[0], op = c[1], obs = c[2], kind = c[3];
    const double val = std::atof(c[4].c_str());
    double err = 0.0;
    if (c.size() >= 7) {
      const double up = std::abs(std::atof(c[5].c_str()));
      const double dn = std::abs(std::atof(c[6].c_str()));
      err = std::max(up, dn);
    }
    TheoryTriplet& t = tab[order][op][obs];
    if (kind == "sm")   { t.sm = val;   t.sm_err = err;   t.has_sm = true; }
    if (kind == "lin")  { t.lin = val;  t.lin_err = err;  t.has_lin = true; }
    if (kind == "quad") { t.quad = val; t.quad_err = err; t.has_quad = true; }
  }
  return tab;
}

TheoryTable load_embedded_theory_csv() {
  std::stringstream ss(EMBEDDED_THEORY_CSV);
  return load_theory_csv_stream(ss);
}

TheoryTable load_theory_csv(const std::string& csv_path) {
  TheoryTable tab;
  std::ifstream f(csv_path);
  if (!f.good()) return tab;
  std::string line;
  std::getline(f, line);
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    auto c = csv_split_simple(line);
    if (c.size() < 6) continue;
    const std::string order = c[0], op = c[1], obs = c[2], kind = c[3];
    const double val = std::atof(c[4].c_str());
    double err = 0.0;
    if (c.size() >= 7) {
      const double up = std::abs(std::atof(c[5].c_str()));
      const double dn = std::abs(std::atof(c[6].c_str()));
      err = std::max(up, dn);
    }
    TheoryTriplet& t = tab[order][op][obs];
    if (kind == "sm")   { t.sm = val;   t.sm_err = err;   t.has_sm = true; }
    if (kind == "lin")  { t.lin = val;  t.lin_err = err;  t.has_lin = true; }
    if (kind == "quad") { t.quad = val; t.quad_err = err; t.has_quad = true; }
  }
  return tab;
}

bool ensure_theory_csv(const std::string& parser, const std::string& csv_path) {
  std::ifstream test(csv_path);
  if (test.good()) return true;
  std::string cmd = "python3 " + parser + " --dump-csv " + csv_path;
  std::cout << "[INFO theory] trying: " << cmd << std::endl;
  int ret = gSystem->Exec(cmd.c_str());
  std::ifstream test2(csv_path);
  if (ret != 0 || !test2.good()) {
    std::cout << "[WARN theory] could not auto-create " << csv_path
              << "; create it with parse_latex_tables.py or pass an existing CSV." << std::endl;
    return false;
  }
  return true;
}

double theory_value(const TheoryTable& tab, const std::string& order, const std::string& op,
                    const std::string& obs_key, double c) {
  if (!tab.count(order)) throw std::runtime_error("missing theory order " + order);
  const auto& byop = tab.at(order);
  if (!byop.count("SM") || !byop.count(op)) throw std::runtime_error("missing theory op " + op);
  const auto& sm = byop.at("SM");
  const auto& oo = byop.at(op);
  if (!sm.count("sigma") || !oo.count("sigma") || !sm.count(obs_key) || !oo.count(obs_key)) {
    throw std::runtime_error("missing theory observable " + obs_key + " for " + op);
  }
  const auto& D0 = sm.at("sigma");
  const auto& D = oo.at("sigma");
  const auto& N0 = sm.at(obs_key);
  const auto& N = oo.at(obs_key);
  const double den = D0.sm + c * D.lin + c * c * D.quad;
  const double num = N0.sm + c * N.lin + c * c * N.quad;
  if (std::abs(den) < 1e-15) return 0.0;
  return num / den;
}

double theory_value_error(const TheoryTable& tab, const std::string& order, const std::string& op,
                          const std::string& obs_key, double c) {
  if (!tab.count(order)) throw std::runtime_error("missing theory order " + order);
  const auto& byop = tab.at(order);
  if (!byop.count("SM") || !byop.count(op)) throw std::runtime_error("missing theory op " + op);
  const auto& sm = byop.at("SM");
  const auto& oo = byop.at(op);
  if (!sm.count("sigma") || !oo.count("sigma") || !sm.count(obs_key) || !oo.count(obs_key)) {
    throw std::runtime_error("missing theory observable " + obs_key + " for " + op);
  }
  const auto& D0 = sm.at("sigma");
  const auto& D = oo.at("sigma");
  const auto& N0 = sm.at(obs_key);
  const auto& N = oo.at(obs_key);

  const double cc = c * c;
  const double den = D0.sm + c * D.lin + cc * D.quad;
  const double num = N0.sm + c * N.lin + cc * N.quad;
  if (std::abs(den) < 1e-15) return 0.0;

  // Propagate the up/down theory uncertainties from the table into the
  // inclusive coefficient shown in the coefficient box. Correlations among
  // numerator and denominator entries are unknown, so this is a conservative
  // uncorrelated propagation using max(up,down) for each table entry.
  const double sig_num = std::sqrt(
      N0.sm_err * N0.sm_err +
      c * c * N.lin_err * N.lin_err +
      cc * cc * N.quad_err * N.quad_err);
  const double sig_den = std::sqrt(
      D0.sm_err * D0.sm_err +
      c * c * D.lin_err * D.lin_err +
      cc * cc * D.quad_err * D.quad_err);

  const double dfdnum = 1.0 / den;
  const double dfdden = -num / (den * den);
  return std::sqrt(std::max(0.0,
      dfdnum * dfdnum * sig_num * sig_num +
      dfdden * dfdden * sig_den * sig_den));
}

HistChunk1D theory_asymmetry_hist_from_reference(const HistChunk1D& ref, const std::string& obs, double coeff) {
  // The theory table gives inclusive spin-correlation numerators. It has no
  // 6-bin shape information. To sync with NanoGEN/data, keep the reference
  // 6-bin intra-half shape and only rescale the backward/forward halves so
  // that the resulting histogram has the requested A_FB.
  HistChunk1D out;
  out.y.assign(BINS_PER_OBS, 0.0);
  out.e.assign(BINS_PER_OBS, 0.0);

  double bsum = 0.0, fsum = 0.0;
  for (int i = 0; i < BINS_PER_OBS / 2; ++i) bsum += ref.y[i];
  for (int i = BINS_PER_OBS / 2; i < BINS_PER_OBS; ++i) fsum += ref.y[i];
  double total = bsum + fsum;
  if (std::abs(total) < 1e-15) total = 1.0;

  const double fac = asymmetry_factor(obs);
  double afb = (fac == 0.0) ? 0.0 : coeff / fac;
  afb = std::max(-0.95, std::min(0.95, afb));

  const double target_b = 0.5 * total * (1.0 - afb);
  const double target_f = 0.5 * total * (1.0 + afb);

  if (std::abs(bsum) > 1e-15) {
    for (int i = 0; i < BINS_PER_OBS / 2; ++i) out.y[i] = ref.y[i] * target_b / bsum;
  } else {
    for (int i = 0; i < BINS_PER_OBS / 2; ++i) out.y[i] = target_b / (BINS_PER_OBS / 2);
  }

  if (std::abs(fsum) > 1e-15) {
    for (int i = BINS_PER_OBS / 2; i < BINS_PER_OBS; ++i) out.y[i] = ref.y[i] * target_f / fsum;
  } else {
    for (int i = BINS_PER_OBS / 2; i < BINS_PER_OBS; ++i) out.y[i] = target_f / (BINS_PER_OBS / 2);
  }

  return out;
}

std::vector<double> linspace_edges(double lo, double hi, int nbins) {
  std::vector<double> edges(nbins + 1);
  for (int i = 0; i <= nbins; ++i) {
    edges[i] = lo + (hi - lo) * double(i) / double(nbins);
  }
  return edges;
}

std::vector<double> obs_axis_edges(const std::string& obs) {
  if (obs == "gen_llbar_delta_phi") return linspace_edges(0.0, TMath::Pi(), 6);
  if (obs == "gen_llbar_delta_eta") return linspace_edges(0.0, 5.0, 6);
  return linspace_edges(-1.0, 1.0, 6);
}

TGraphAsymmErrors make_theory_band_graph(const TH1D& hcen, const TH1D& herr, const std::string& name) {
  const int n = hcen.GetNbinsX();
  std::vector<double> x(n), y(n), exlo(n), exhi(n), eylo(n), eyhi(n);
  for (int i = 1; i <= n; ++i) {
    x[i-1] = hcen.GetBinCenter(i);
    y[i-1] = hcen.GetBinContent(i);
    exlo[i-1] = 0.5 * hcen.GetBinWidth(i);
    exhi[i-1] = 0.5 * hcen.GetBinWidth(i);
    const double e = std::max(0.0, herr.GetBinError(i));
    eylo[i-1] = e;
    eyhi[i-1] = e;
  }
  TGraphAsymmErrors gr(n, x.data(), y.data(), exlo.data(), exhi.data(), eylo.data(), eyhi.data());
  gr.SetName(name.c_str());
  return gr;
}


void save_individual_distribution_plot(const std::string& outpng,
                                       const std::string& obs,
                                       const HistChunk1D& data,
                                       const HistChunk1D& eft,
                                       const HistChunk1D* theory,
                                       const std::string& title,
                                       double theory_coeff_err_override = -1.0) {
  TCanvas c("c_ind", "c_ind", 850, 750);
  
  gStyle->SetOptStat(0);TPad p1("p1", "p1", 0.0, 0.30, 1.0, 1.0);
  TPad p2("p2", "p2", 0.0, 0.0, 1.0, 0.30);
  p1.SetBottomMargin(0.03); p1.SetLeftMargin(0.13); p1.SetRightMargin(0.04);
  p2.SetTopMargin(0.04); p2.SetBottomMargin(0.33); p2.SetLeftMargin(0.13); p2.SetRightMargin(0.04);
  p1.Draw(); p2.Draw();

  std::vector<double> xedges = obs_axis_edges(obs);
  TH1D hdata("hdata", "", BINS_PER_OBS, xedges.data());
  TH1D heft("heft", "", BINS_PER_OBS, xedges.data());
  TH1D hthy("hthy", "", BINS_PER_OBS, xedges.data());
  TH1D hthy_band("hthy_band", "", BINS_PER_OBS, xedges.data());
  for (int i = 0; i < BINS_PER_OBS; ++i) {
    hdata.SetBinContent(i+1, data.y[i]); hdata.SetBinError(i+1, data.e[i]);
    heft.SetBinContent(i+1, eft.y[i]); heft.SetBinError(i+1, eft.e[i]);
    if (theory) { hthy.SetBinContent(i+1, theory->y[i]); hthy.SetBinError(i+1, theory->e[i]); }
  }

  p1.cd();
  hdata.SetMarkerStyle(20); hdata.SetMarkerColor(kBlack); hdata.SetLineColor(kBlack);
  heft.SetLineColor(TColor::GetColor("#0072B2")); heft.SetLineWidth(3);
  const int theoryRed = TColor::GetColor("#D62728");
  hthy.SetLineColor(theoryRed); hthy.SetLineWidth(3); hthy.SetLineStyle(2);
  hthy_band.SetFillColorAlpha(theoryRed, 0.30);
  hthy_band.SetFillStyle(1001);
  hthy_band.SetLineColor(theoryRed);
  hthy_band.SetLineWidth(0);
  hthy_band.SetMarkerSize(0);

  if (theory) {
    const CoeffSummary ct_for_band = coefficient_from_chunk(obs, *theory);
    const double coeff_band_err = (theory_coeff_err_override >= 0.0)
                                  ? theory_coeff_err_override
                                  : (ct_for_band.ok ? ct_for_band.coeff_err : 0.0);
    HistChunk1D thy_up = theory_asymmetry_hist_from_reference(data, obs, ct_for_band.coeff + coeff_band_err);
    HistChunk1D thy_dn = theory_asymmetry_hist_from_reference(data, obs, ct_for_band.coeff - coeff_band_err);
    for (int i = 0; i < BINS_PER_OBS; ++i) {
      const double y0 = theory->y[i];
      const double eup = std::abs(thy_up.y[i] - y0);
      const double edn = std::abs(y0 - thy_dn.y[i]);
      hthy_band.SetBinContent(i+1, y0);
      hthy_band.SetBinError(i+1, std::max(eup, edn));
    }
  }
  hdata.GetYaxis()->SetTitle("normalized diff. xsec");
  hdata.GetYaxis()->SetTitleSize(0.055); hdata.GetYaxis()->SetLabelSize(0.045);
  hdata.GetXaxis()->SetLabelSize(0.0);
  double ymax = 0.0;
  for (int i=1;i<=BINS_PER_OBS;++i) {
    ymax = std::max(ymax, hdata.GetBinContent(i) + hdata.GetBinError(i));
    ymax = std::max(ymax, heft.GetBinContent(i) + heft.GetBinError(i));
    if (theory) ymax = std::max(ymax, hthy.GetBinContent(i) + hthy_band.GetBinError(i));
  }
  hdata.SetMinimum(0.0);
  hdata.SetMaximum(1.8 * std::max(1e-9, ymax));
  TGraphAsymmErrors gr_theory_band;
  if (theory) {
    gr_theory_band = make_theory_band_graph(hthy, hthy_band, "gr_theory_band");
    gr_theory_band.SetFillColorAlpha(theoryRed, 0.30);
    gr_theory_band.SetFillStyle(1001);
    gr_theory_band.SetLineColor(theoryRed);
    gr_theory_band.SetLineWidth(0);
    gr_theory_band.SetMarkerSize(0);
  }
  hdata.Draw("E1");
  if (theory) gr_theory_band.Draw("2 same");
  heft.Draw("hist same");
  if (theory) hthy.Draw("hist same");
  hdata.Draw("E1 same");
  // Main-panel legend: move slightly to the right to avoid the coefficient box.
  TLegend leg(0.66, 0.70, 0.965, 0.88);
  leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.034);
  leg.AddEntry(&hdata, "Data", "lep");
  leg.AddEntry(&heft, "Dim6Top MC", "l");
  std::string theoryLegend = "Theory ref";
  if (title.find("Theory LO") != std::string::npos) theoryLegend = "Theory LO";
  if (title.find("Theory NLO") != std::string::npos) theoryLegend = "Theory NLO";
  if (theory) {
    leg.AddEntry(&hthy, theoryLegend.c_str(), "l");
    leg.AddEntry(&hthy_band, (theoryLegend + " unc.").c_str(), "f");
  }
  leg.Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextFont(42);
  lat.SetTextSize(0.040);
  lat.SetTextAlign(11);
  lat.DrawLatex(0.13, 0.93, "#bf{CMS} #it{Simulation Work in Progress}");
  lat.SetTextAlign(31);
  lat.SetTextSize(0.034);
  lat.DrawLatex(0.96, 0.93, "138 fb^{-1} (13 TeV)");
  lat.SetTextAlign(11);

  CoeffSummary cd = coefficient_from_chunk(obs, data);
  CoeffSummary ce = coefficient_from_chunk(obs, eft);
  CoeffSummary ct;
  if (theory) ct = coefficient_from_chunk(obs, *theory);
  if (cd.ok || ce.ok || ct.ok) {
    TPaveText *coeffBox = new TPaveText(0.13, 0.72, 0.48, 0.88, "NDC");
    coeffBox->SetFillColor(0);
    coeffBox->SetFillStyle(0);
    coeffBox->SetBorderSize(1);
    coeffBox->SetTextFont(42);
    coeffBox->SetTextAlign(12);
    coeffBox->SetTextSize(0.027);
    const std::string lab = obs_root_label(obs);
    if (cd.ok) coeffBox->AddText(Form("%s [Data] = %.4f #pm %.4f", lab.c_str(), cd.coeff, cd.coeff_err));
    if (ce.ok) coeffBox->AddText(Form("%s [EFT MC] = %.4f #pm %.4f", lab.c_str(), ce.coeff, ce.coeff_err));
    if (ct.ok) {
      const double theory_box_err = (theory_coeff_err_override >= 0.0) ? theory_coeff_err_override : ct.coeff_err;
      coeffBox->AddText(Form("%s [%s] = %.4f #pm %.4f", lab.c_str(), theoryLegend.c_str(), ct.coeff, theory_box_err));
    }
    coeffBox->Draw("same");
  }

  p2.cd();
  TH1D hratio_data("hratio_data", "", BINS_PER_OBS, xedges.data());
  TH1D hratio_mc("hratio_mc", "", BINS_PER_OBS, xedges.data());
  TH1D hratio_theory_band("hratio_theory_band", "", BINS_PER_OBS, xedges.data());

  const bool useTheoryDenom = (theory != nullptr);
  const std::string denomLabel = useTheoryDenom ? theoryLegend : std::string("Data");

  for (int i = 1; i <= BINS_PER_OBS; ++i) {
    const double d  = hdata.GetBinContent(i);
    const double ed = hdata.GetBinError(i);
    const double m  = heft.GetBinContent(i);
    const double em = heft.GetBinError(i);

    const double den  = useTheoryDenom ? hthy.GetBinContent(i) : d;
    const double eden = useTheoryDenom ? hthy.GetBinError(i)   : ed;

    auto set_ratio = [&](TH1D& h, double num, double enum_) {
      const double r = std::abs(den) > 0.0 ? num / den : 0.0;
      double er = 0.0;
      if (std::abs(den) > 0.0) {
        er = std::sqrt((enum_/den)*(enum_/den) + (num*eden/(den*den))*(num*eden/(den*den)));
      }
      h.SetBinContent(i, r);
      h.SetBinError(i, er);
    };

    set_ratio(hratio_data, d, ed);
    set_ratio(hratio_mc,   m, em);

    if (useTheoryDenom) {
      hratio_theory_band.SetBinContent(i, 1.0);
      hratio_theory_band.SetBinError(i, std::abs(den) > 0.0 ? hthy_band.GetBinError(i) / std::abs(den) : 0.0);
    }
  }

  hratio_data.SetMarkerStyle(20);
  hratio_data.SetMarkerColor(kBlack);
  hratio_data.SetLineColor(kBlack);
  hratio_data.GetYaxis()->SetTitle(useTheoryDenom ? "ratio to theory" : "ratio");
  hratio_data.GetYaxis()->SetRangeUser(0.8, 1.2);
  hratio_data.GetYaxis()->SetTitleSize(0.10);
  hratio_data.GetYaxis()->SetLabelSize(0.09);
  hratio_data.GetYaxis()->SetTitleOffset(0.55);
  hratio_data.GetXaxis()->SetTitle(obs_root_label(obs).c_str());
  hratio_data.GetXaxis()->SetNdivisions(506);
  hratio_data.GetXaxis()->SetTitleSize(0.12);
  hratio_data.GetXaxis()->SetLabelSize(0.10);

  hratio_mc.SetMarkerStyle(20);
  hratio_mc.SetMarkerColor(TColor::GetColor("#0072B2"));
  hratio_mc.SetLineColor(TColor::GetColor("#0072B2"));

  TGraphAsymmErrors gr_ratio_theory_band;
  if (useTheoryDenom) {
    gr_ratio_theory_band = make_theory_band_graph(hratio_theory_band, hratio_theory_band, "gr_ratio_theory_band");
    gr_ratio_theory_band.SetFillColorAlpha(theoryRed, 0.30);
    gr_ratio_theory_band.SetFillStyle(1001);
    gr_ratio_theory_band.SetLineColor(theoryRed);
    gr_ratio_theory_band.SetLineWidth(0);
    gr_ratio_theory_band.SetMarkerSize(0);
  }

  hratio_data.Draw("E1");
  if (useTheoryDenom) gr_ratio_theory_band.Draw("2 same");
  hratio_data.Draw("E1 same");
  hratio_mc.Draw("E1 same");
  TLine l(xedges.front(), 1.0, xedges.back(), 1.0); l.SetLineStyle(2); l.SetLineColor(kGray+2); l.Draw();

  TLegend rleg(0.68, 0.70, 0.965, 0.92);
  rleg.SetBorderSize(0);
  rleg.SetFillStyle(0);
  rleg.SetTextFont(42);
  rleg.SetTextSize(0.074);
  rleg.AddEntry(&hratio_data, (std::string("Data/") + denomLabel).c_str(), "lep");
  rleg.AddEntry(&hratio_mc,   (std::string("MC/")   + denomLabel).c_str(), "lep");
  rleg.Draw();

  c.SaveAs(outpng.c_str());
}


std::vector<double> theory_full_values_for_wc_order(const std::map<std::string, HistChunk1D>& data_chunks,
                                                    const TheoryTable& theory,
                                                    const std::string& order,
                                                    const std::string& wc,
                                                    double cval);

std::vector<double> flatten_chunks_in_observable_order(const std::map<std::string, HistChunk1D>& chunks, bool errors=false) {
  std::vector<double> out;
  out.reserve(N_OBS_TOTAL * BINS_PER_OBS);
  const auto names = observable_names();
  for (const auto& obs : names) {
    auto it = chunks.find(obs);
    for (int ib = 0; ib < BINS_PER_OBS; ++ib) {
      if (it != chunks.end() && ib < (int)it->second.y.size()) {
        out.push_back(errors ? it->second.e[ib] : it->second.y[ib]);
      } else {
        out.push_back(0.0);
      }
    }
  }
  return out;
}

void make_concatenated_inspection_plot_root(const std::string& data_root,
                                            const std::string& eft_template_pattern,
                                            const std::string& wc,
                                            const std::string& outdir,
                                            const TheoryTable* theory = nullptr,
                                            const std::string& theory_order = "LO") {
  const std::string od = outdir + "/inspect_" + wc;
  gSystem->mkdir(od.c_str(), true);

  const auto obs_names = observable_names();
  const int nbins = (int)obs_names.size() * BINS_PER_OBS;
  auto data_chunks = load_hist_chunks_from_root(data_root);
  const std::vector<double> data_y = flatten_chunks_in_observable_order(data_chunks, false);
  const std::vector<double> data_e = flatten_chunks_in_observable_order(data_chunks, true);

  TH1D hdata("hdata_concat", "", nbins, 0.0, (double)nbins);
  for (int i = 0; i < nbins; ++i) {
    hdata.SetBinContent(i + 1, data_y[i]);
    hdata.SetBinError(i + 1, data_e[i]);
  }

  const std::vector<int> vals = {-8, -4, -2, 0, 2, 4, 8};
  const int colors[] = {
    kOrange+7, kAzure+7, kGreen+2, kBlack, kAzure-3, kOrange+1, kMagenta-4
  };

  std::map<int, TH1D*> hmc;
  std::map<int, TH1D*> hratio;
  std::map<int, TH1D*> hthy;
  double ymax = 0.0;
  for (int i = 0; i < nbins; ++i) ymax = std::max(ymax, data_y[i] + data_e[i]);

  for (size_t iv = 0; iv < vals.size(); ++iv) {
    const int v = vals[iv];
    const std::string path = make_template_path(eft_template_pattern, wc, v);
    try {
      auto mc_chunks = load_hist_chunks_from_root(path, true, "fit-mc");
      const std::vector<double> mc_y = flatten_chunks_in_observable_order(mc_chunks, false);
      const std::vector<double> mc_e = flatten_chunks_in_observable_order(mc_chunks, true);
      TH1D* hm = new TH1D(Form("hmc_concat_%s_%d", wc.c_str(), v), "", nbins, 0.0, (double)nbins);
      TH1D* hr = new TH1D(Form("hratio_concat_%s_%d", wc.c_str(), v), "", nbins, 0.0, (double)nbins);
      hm->SetDirectory(nullptr);
      hr->SetDirectory(nullptr);
      for (int i = 0; i < nbins; ++i) {
        hm->SetBinContent(i + 1, mc_y[i]);
        hm->SetBinError(i + 1, mc_e[i]);
        ymax = std::max(ymax, mc_y[i] + mc_e[i]);
        const double d = data_y[i];
        const double ed = data_e[i];
        const double n = mc_y[i];
        const double en = mc_e[i];
        const double r = std::abs(d) > 0 ? n / d : 0.0;
        double er = 0.0;
        if (std::abs(d) > 0) er = std::sqrt((en/d)*(en/d) + (n*ed/(d*d))*(n*ed/(d*d)));
        hr->SetBinContent(i + 1, r);
        hr->SetBinError(i + 1, er);
      }
      hm->SetLineColor(colors[iv]);
      hm->SetLineWidth(v == 0 ? 4 : 2);
      hm->SetLineStyle(1);
      hr->SetLineColor(colors[iv]);
      hr->SetMarkerColor(colors[iv]);
      hr->SetLineWidth(2);
      hr->SetMarkerStyle(20);
      hr->SetMarkerSize(0.35);
      hmc[v] = hm;
      hratio[v] = hr;

      if (theory) {
        const std::vector<double> ty = theory_full_values_for_wc_order(data_chunks, *theory, theory_order, wc, (double)v);
        TH1D* ht = new TH1D(Form("hthy_concat_%s_%s_%d", theory_order.c_str(), wc.c_str(), v), "", nbins, 0.0, (double)nbins);
        ht->SetDirectory(nullptr);
        for (int i = 0; i < nbins && i < (int)ty.size(); ++i) ht->SetBinContent(i + 1, ty[i]);
        ht->SetLineColor(colors[iv]);
        ht->SetLineWidth(1);
        ht->SetLineStyle(2);
        hthy[v] = ht;
      }
    } catch (const std::exception& e) {
      std::cout << "[WARN concat inspect] skip " << path << ": " << e.what() << std::endl;
    }
  }

  TCanvas c("c_concat_inspect", "c_concat_inspect", 1500, 980);
  
  gStyle->SetOptStat(0);c.SetTopMargin(0.03);
  c.SetBottomMargin(0.03);
  TPad p1("p1_concat", "p1_concat", 0.0, 0.30, 1.0, 1.0);
  TPad p2("p2_concat", "p2_concat", 0.0, 0.0, 1.0, 0.30);
  p1.SetBottomMargin(0.02);
  p1.SetLeftMargin(0.08);
  p1.SetRightMargin(0.03);
  p1.SetTopMargin(0.08);
  p2.SetTopMargin(0.02);
  p2.SetBottomMargin(0.36);
  p2.SetLeftMargin(0.08);
  p2.SetRightMargin(0.03);
  p1.Draw();
  p2.Draw();

  p1.cd();
  hdata.SetMarkerStyle(20);
  hdata.SetMarkerSize(0.6);
  hdata.SetLineColor(kBlack);
  hdata.SetMarkerColor(kBlack);
  hdata.GetYaxis()->SetTitle("Normalized\nDifferential\nCross-section");
  hdata.GetYaxis()->SetTitleSize(0.055);
  hdata.GetYaxis()->SetTitleOffset(0.63);
  hdata.GetYaxis()->SetLabelSize(0.040);
  hdata.GetXaxis()->SetLabelSize(0.0);
  hdata.SetMaximum(1.8 * std::max(1e-9, ymax));
  hdata.SetMinimum(0.0);
  hdata.Draw("E1");
  for (const auto& kv : hmc) kv.second->Draw("hist same");
  for (const auto& kv : hthy) kv.second->Draw("hist same");
  hdata.Draw("E1 same");

  for (int iobs = 1; iobs < (int)obs_names.size(); ++iobs) {
    const double x = iobs * BINS_PER_OBS;
    TLine* l = new TLine(x, 0.0, x, hdata.GetMaximum());
    l->SetLineColor(kGray+1);
    l->SetLineStyle(3);
    l->SetLineWidth(1);
    l->Draw("same");
  }

  TLegend leg(0.66, 0.63, 0.94, 0.91);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.030);
  leg.SetNColumns(2);
  for (int v : vals) {
    if (hmc.count(v)) leg.AddEntry(hmc[v], Form("%s=%+d", wc.c_str(), v), "l");
  }
  leg.AddEntry(&hdata, "Data", "lep");
  if (theory && !hthy.empty()) leg.AddEntry(hthy.begin()->second, ("Theory " + theory_order).c_str(), "l");
  leg.Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextFont(42);
  lat.SetTextSize(0.034);
  lat.SetTextAlign(11);
  lat.DrawLatex(0.08, 0.955, "#bf{CMS} #it{Simulation Work in Progress}");
  lat.SetTextAlign(31);
  lat.DrawLatex(0.97, 0.955, "138 fb^{-1} (13 TeV)");
  lat.SetTextAlign(11);

  p2.cd();
  TH1D frame_ratio("frame_ratio_concat", "", nbins, 0.0, (double)nbins);
  frame_ratio.GetYaxis()->SetTitle(theory ? ("MC/Theory " + theory_order).c_str() : "MC/Data");
  frame_ratio.GetYaxis()->SetRangeUser(0.5, 1.5);
  frame_ratio.GetYaxis()->SetTitleSize(0.105);
  frame_ratio.GetYaxis()->SetTitleOffset(0.34);
  frame_ratio.GetYaxis()->SetLabelSize(0.075);
  frame_ratio.GetXaxis()->SetTitle("Observable");
  frame_ratio.GetXaxis()->SetTitleSize(0.105);
  frame_ratio.GetXaxis()->SetLabelSize(0.060);
  frame_ratio.GetXaxis()->SetTickLength(0.08);
  for (int iobs = 0; iobs < (int)obs_names.size(); ++iobs) {
    frame_ratio.GetXaxis()->SetBinLabel(iobs * BINS_PER_OBS + BINS_PER_OBS / 2 + 1,
                                        obs_root_label(obs_names[iobs]).c_str());
  }
  frame_ratio.LabelsOption("v", "X");
  frame_ratio.Draw("AXIS");
  TLine unity(0.0, 1.0, (double)nbins, 1.0);
  unity.SetLineStyle(2);
  unity.SetLineColor(kGray+2);
  unity.Draw("same");
  for (const auto& kv : hratio) kv.second->Draw("hist same");
  for (int iobs = 1; iobs < (int)obs_names.size(); ++iobs) {
    const double x = iobs * BINS_PER_OBS;
    TLine* l = new TLine(x, 0.5, x, 1.5);
    l->SetLineColor(kGray+1);
    l->SetLineStyle(3);
    l->SetLineWidth(1);
    l->Draw("same");
  }

  const std::string suffix = theory ? ("_" + theory_order) : "";
  const std::string pdf = od + "/inspect_" + wc + "_concatenated" + suffix + ".pdf";
  const std::string png = od + "/inspect_" + wc + "_concatenated" + suffix + ".png";
  c.SaveAs(pdf.c_str());
  c.SaveAs(png.c_str());
  std::cout << "[SAVED] " << pdf << std::endl;
  std::cout << "[SAVED] " << png << std::endl;
}

void make_individual_distribution_plots_and_coefficients_root(const std::string& data_root,
                                                              const std::string& eft_template_pattern,
                                                              const std::string& wc,
                                                              const std::string& outdir,
                                                              const TheoryTable* theory = nullptr,
                                                              const std::string& theory_order = "NLO") {
  std::string od = outdir + "/inspect_" + wc;
  std::string pd = od + "/individual";
  gSystem->mkdir(od.c_str(), true); gSystem->mkdir(pd.c_str(), true);
  auto data_chunks = load_hist_chunks_from_root(data_root);
  std::ofstream csv(od + "/coefficients_" + wc + ".csv");
  csv << "wc,wc_value,observable,observable_label,data_coefficient,data_coefficient_err,eft_coefficient,eft_coefficient_err,theory_coefficient,theory_coefficient_err\n";

  const std::vector<int> vals = {-8, -4, -2, 0, 2, 4, 8};
  for (int v : vals) {
    std::string path = make_template_path(eft_template_pattern, wc, v);
    try {
      auto eft_chunks = load_hist_chunks_from_root(path, true, "plot-eft");
      for (const auto& kv : data_chunks) {
        const std::string& obs = kv.first;
        if (!eft_chunks.count(obs)) continue;
        HistChunk1D* theory_ptr = nullptr;
        HistChunk1D theory_hist;
        double theory_coeff = 0.0;
        double theory_coeff_err = -1.0;
        if (theory) {
          const std::string op = theory_op_from_wc(wc);
          const std::string tkey = theory_obs_key(obs);
          if (!tkey.empty()) {
            try {
              theory_coeff = theory_value(*theory, theory_order, op, tkey, (double)v);
              theory_coeff_err = theory_value_error(*theory, theory_order, op, tkey, (double)v);
              theory_hist = theory_asymmetry_hist_from_reference(kv.second, obs, theory_coeff);
              theory_ptr = &theory_hist;
            } catch (...) {}
          }
        }
        std::stringstream name;
        name << pd << "/" << obs << "_" << wc << "_" << v << ".png";
        std::stringstream title;
        title << wc << "=" << v;
        if (theory_ptr) title << " Theory " << theory_order;
        save_individual_distribution_plot(name.str(), obs, kv.second, eft_chunks[obs], theory_ptr, title.str(), theory_coeff_err);
        auto cd = coefficient_from_chunk(obs, kv.second);
        auto ce = coefficient_from_chunk(obs, eft_chunks[obs]);
        csv << wc << "," << v << "," << obs << ",\"" << obs_root_label(obs) << "\",";
        csv << (cd.ok ? cd.coeff : 0.0) << "," << (cd.ok ? cd.coeff_err : 0.0) << ",";
        csv << (ce.ok ? ce.coeff : 0.0) << "," << (ce.ok ? ce.coeff_err : 0.0) << "," << theory_coeff << "," << (theory_coeff_err >= 0.0 ? theory_coeff_err : 0.0) << "\n";
      }
    } catch (const std::exception& e) {
      std::cout << "[WARN inspect] skip " << path << ": " << e.what() << std::endl;
    }
  }
  csv.close();
  std::cout << "[SAVED] " << od << "/coefficients_" << wc << ".csv" << std::endl;
}


std::vector<double> theory_full_values_for_wc_order(const std::map<std::string, HistChunk1D>& data_chunks,
                                                    const TheoryTable& theory,
                                                    const std::string& order,
                                                    const std::string& wc,
                                                    double cval) {
  std::vector<double> out;
  out.reserve(N_OBS_TOTAL * BINS_PER_OBS);
  const auto names = observable_names();
  const std::string op = theory_op_from_wc(wc);
  for (const auto& obs : names) {
    HistChunk1D chunk;
    auto itref = data_chunks.find(obs);
    if (itref != data_chunks.end()) chunk = itref->second;
    else { chunk.y.assign(BINS_PER_OBS, 0.0); chunk.e.assign(BINS_PER_OBS, 0.0); }

    const std::string tkey = theory_obs_key(obs);
    if (!tkey.empty()) {
      try {
        const double coeff = theory_value(theory, order, op, tkey, cval);
        chunk = theory_asymmetry_hist_from_reference(chunk, obs, coeff);
      } catch (const std::exception& e) {
        // No inclusive theory table entry for this WC/observable/order.
        // Keep the reference chunk so this bin does not contribute fake theory dependence.
      }
    }
    for (int ib = 0; ib < BINS_PER_OBS; ++ib) out.push_back(chunk.y[ib]);
  }
  return out;
}

FitResult1D fit_one_wc_theory(const std::string& wc,
                              const std::vector<int>& obs,
                              const std::string& tag,
                              const std::string& outdir,
                              const std::string& data_root,
                              const TheoryTable& theory,
                              const std::string& order,
                              const TMatrixD& cov_full,
                              int drop_bin_idx,
                              double scan_min,
                              double scan_max,
                              int scan_n,
                              double xscale = 1.0,
                              const std::string& xlabel_override = "") {
  std::vector<int> keep = build_keep_indices(obs, drop_bin_idx);
  auto data_chunks = load_hist_chunks_from_root(data_root);
  auto data_full = load_values(data_root);
  auto data = select_vec(data_full, keep);

  TMatrixD cov = select_cov(cov_full, keep, drop_bin_idx);
  TDecompSVD svd(cov);
  TMatrixD cov_inv = svd.Invert();
  const int nb = data.size();

  auto pred_at = [&](double x) {
    return select_vec(theory_full_values_for_wc_order(data_chunks, theory, order, wc, x), keep);
  };

  auto chi2_at = [&](double x) {
    auto p = pred_at(x);
    TVectorD d(nb);
    for (int i = 0; i < nb; ++i) d(i) = data[i] - p[i];
    TVectorD tmp = cov_inv * d;
    return d * tmp;
  };

  std::vector<double> xs(scan_n), chi(scan_n), dchi(scan_n);
  double best = 0.0, chi_min = 1e300;
  for (int i = 0; i < scan_n; ++i) {
    double x = scan_min + (scan_max - scan_min) * double(i) / double(scan_n - 1);
    xs[i] = x;
    chi[i] = chi2_at(x);
    if (chi[i] < chi_min) { chi_min = chi[i]; best = x; }
  }
  for (int i = 0; i < scan_n; ++i) dchi[i] = chi[i] - chi_min;

  auto i68 = interval_crossing(xs, dchi, 1.0, best);
  auto i95 = interval_crossing(xs, dchi, 4.0, best);

  FitResult1D r;
  r.name = wc;
  r.label = pretty_label(wc);
  r.best = best;
  r.chi2min = chi_min;
  r.lo68 = i68.first;
  r.hi68 = i68.second;
  r.lo95 = i95.first;
  r.hi95 = i95.second;
  r.nbins = nb;

  const std::string od = outdir + "/" + tag;
  gSystem->mkdir(od.c_str(), true);
  std::ofstream csv(od + "/scan_" + wc + ".csv");
  csv << wc << ",chi2,delta_chi2\n";
  for (int i = 0; i < scan_n; ++i) csv << xs[i] << "," << chi[i] << "," << dchi[i] << "\n";
  csv.close();

  std::ofstream js(od + "/fit_result_" + wc + ".json");
  js << "{\n";
  js << "  \"wc\": \"" << wc << "\",\n";
  js << "  \"tag\": \"" << tag << "\",\n";
  js << "  \"theory_order\": \"" << order << "\",\n";
  js << "  \"best\": " << r.best << ",\n";
  js << "  \"chi2_min\": " << r.chi2min << ",\n";
  js << "  \"lo68\": " << r.lo68 << ",\n";
  js << "  \"hi68\": " << r.hi68 << ",\n";
  js << "  \"lo95\": " << r.lo95 << ",\n";
  js << "  \"hi95\": " << r.hi95 << ",\n";
  js << "  \"n_bins_fit\": " << r.nbins << "\n";
  js << "}\n";
  js.close();

  save_chi2_plot(od + "/deltaChi2_" + wc + ".pdf", od + "/deltaChi2_" + wc + ".png",
                 xs, dchi, r, xscale, xlabel_override);

  std::cout << "[THEORY " << order << "] " << tag << " " << wc
            << " best=" << best
            << " 68=[" << r.lo68 << "," << r.hi68 << "]"
            << " 95=[" << r.lo95 << "," << r.hi95 << "]"
            << " nbins=" << nb << std::endl;
  return r;
}

void fit_2d_pair_grid_theory(const std::string& wc1,
                             const std::string& wc2,
                             const std::vector<int>& obs,
                             const std::string& tag,
                             const std::string& outdir,
                             const std::string& data_root,
                             const TheoryTable& theory,
                             const std::string& order,
                             const TMatrixD& cov_full,
                             int drop_bin_idx,
                             double scan_min,
                             double scan_max,
                             int ngrid) {
  std::vector<int> keep = build_keep_indices(obs, drop_bin_idx);
  auto data_chunks = load_hist_chunks_from_root(data_root);
  auto data_full = load_values(data_root);
  auto data = select_vec(data_full, keep);

  TMatrixD cov = select_cov(cov_full, keep, drop_bin_idx);
  TDecompSVD svd(cov);
  TMatrixD cov_inv = svd.Invert();
  const int nb = data.size();

  auto sm1_full = theory_full_values_for_wc_order(data_chunks, theory, order, wc1, 0.0);
  auto sm2_full = theory_full_values_for_wc_order(data_chunks, theory, order, wc2, 0.0);
  auto sm_ref = select_vec(sm1_full, keep);

  auto pred = [&](double x, double y) {
    auto p1 = select_vec(theory_full_values_for_wc_order(data_chunks, theory, order, wc1, x), keep);
    auto p2 = select_vec(theory_full_values_for_wc_order(data_chunks, theory, order, wc2, y), keep);
    auto r2 = select_vec(sm2_full, keep);
    std::vector<double> p(nb, 0.0);
    for (int i = 0; i < nb; ++i) p[i] = sm_ref[i] + (p1[i] - sm_ref[i]) + (p2[i] - r2[i]);
    return p;
  };

  auto chi2 = [&](double x, double y) {
    auto p = pred(x, y);
    TVectorD d(nb);
    for (int i = 0; i < nb; ++i) d(i) = data[i] - p[i];
    TVectorD tmp = cov_inv * d;
    return d * tmp;
  };

  TH2D h("h_theory2d", "", ngrid, scan_min, scan_max, ngrid, scan_min, scan_max);
  double cmin = 1e300, bx = 0.0, by = 0.0;
  for (int ix = 1; ix <= ngrid; ++ix) {
    double x = h.GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= ngrid; ++iy) {
      double y = h.GetYaxis()->GetBinCenter(iy);
      double c = chi2(x, y);
      h.SetBinContent(ix, iy, c);
      if (c < cmin) { cmin = c; bx = x; by = y; }
    }
  }
  for (int ix = 1; ix <= ngrid; ++ix)
    for (int iy = 1; iy <= ngrid; ++iy)
      h.SetBinContent(ix, iy, h.GetBinContent(ix, iy) - cmin);

  std::string od = outdir + "/" + tag + "_2D";
  gSystem->mkdir(od.c_str(), true);

  TCanvas c("c2d_theory", "c2d_theory", 760, 700);
  
  gStyle->SetOptStat(0);c.SetLeftMargin(0.13); c.SetRightMargin(0.15); c.SetTopMargin(0.11); c.SetBottomMargin(0.13);
  h.SetTitle("");
  h.GetXaxis()->SetTitle(pretty_label(wc1).c_str());
  h.GetYaxis()->SetTitle(pretty_label(wc2).c_str());
  h.GetZaxis()->SetTitle("#Delta#chi^{2}");
  h.SetMinimum(0); h.SetMaximum(10);
  h.Draw("AXIS");
  h.SetContour(2);
  h.SetContourLevel(0, 2.30);
  h.SetContourLevel(1, 5.99);
  h.SetLineColor(TColor::GetColor("#D55E00"));
  h.SetLineWidth(3);
  h.Draw("CONT3 same");
  TMarker m(bx, by, 29); m.SetMarkerColor(kBlack); m.SetMarkerSize(1.8); m.Draw();
  TLatex latex; latex.SetNDC(); latex.SetTextFont(42); latex.SetTextSize(0.032);
  latex.DrawLatex(0.13, 0.955, "#bf{CMS} #it{Work in Progress}");
  latex.SetTextAlign(31); latex.DrawLatex(0.94, 0.955, "138 fb^{-1} (13 TeV)"); latex.SetTextAlign(11);
  c.SaveAs((od + "/contour_" + wc1 + "_vs_" + wc2 + ".pdf").c_str());
  c.SaveAs((od + "/contour_" + wc1 + "_vs_" + wc2 + ".png").c_str());

  std::ofstream js(od + "/fit2d_" + wc1 + "_vs_" + wc2 + ".json");
  js << "{\n";
  js << "  \"wc1\": \"" << wc1 << "\",\n";
  js << "  \"wc2\": \"" << wc2 << "\",\n";
  js << "  \"theory_order\": \"" << order << "\",\n";
  js << "  \"best1\": " << bx << ",\n";
  js << "  \"best2\": " << by << ",\n";
  js << "  \"chi2_min\": " << cmin << ",\n";
  js << "  \"n_bins_fit\": " << nb << "\n";
  js << "}\n";
  js.close();

  std::cout << "[THEORY 2D " << order << "] " << tag << " " << wc1 << " vs " << wc2
            << " best=(" << bx << "," << by << ") chi2min=" << cmin << std::endl;
}

void add_linear_anom_entry(std::vector<SummaryEntry>& anom_summary,
                           const std::map<std::string, FitResult1D>& result_by_key,
                           const std::string& key,
                           const std::string& label,
                           const std::vector<std::pair<std::string,double>>& terms,
                           double scale) {
  double best = 0.0, varlo = 0.0, varhi = 0.0;
  bool ok = true;
  for (const auto& kv : terms) {
    auto it = result_by_key.find(kv.first);
    if (it == result_by_key.end()) { ok = false; break; }
    const auto& r = it->second;
    best += kv.second * r.best;
    const double elo = std::max(0.0, r.best - r.lo68);
    const double ehi = std::max(0.0, r.hi68 - r.best);
    varlo += (kv.second * elo) * (kv.second * elo);
    varhi += (kv.second * ehi) * (kv.second * ehi);
  }
  if (!ok) { std::cout << "[SKIP theory anom] missing WC input for " << key << std::endl; return; }
  SummaryEntry e;
  e.key = key; e.label = label;
  e.best = best * scale;
  e.lo68 = (best - std::sqrt(varlo)) * scale;
  e.hi68 = (best + std::sqrt(varhi)) * scale;
  anom_summary.push_back(e);
}

void run_theory_fit_suite(const std::string& order,
                          const std::string& suffix,
                          const std::string& outdir,
                          const std::string& data_root,
                          const TheoryTable& theory,
                          const TMatrixD& cov_full,
                          int drop_bin_idx,
                          double scan_min,
                          double scan_max,
                          int scan_n,
                          int scan2d_n,
                          const std::vector<std::string>& wc_list,
                          const std::map<std::string, std::vector<int>>& obs_sets,
                          const std::map<std::string, std::vector<int>>& fig16_obs,
                          const std::vector<PairFitSpec>& fig18_pairs) {
  std::cout << "\n[INFO theory] starting " << order << " data-vs-theory 1D/2D fits" << std::endl;
  std::vector<FitResult1D> results;
  const double mt = 0.1725;
  const double mu_scale = 2.0 * mt * mt;

  // Single source of truth: every 1D WC fit uses the assigned observable set.
  // For WCs explicitly listed in fig16_obs, use that AN-22-028 subset.
  // For the remaining WCs, the assigned fallback is all_0_35.
  std::map<std::string, FitResult1D> assigned_result_by_key;

  auto assigned_obs_for_wc = [&](const std::string& wc) -> std::vector<int> {
    auto it = fig16_obs.find(wc);
    if (it != fig16_obs.end()) return it->second;
    return obs_sets.at("all_0_35");
  };

  for (const auto& wc : wc_list) {
    FitResult1D r = fit_one_wc_theory(wc, assigned_obs_for_wc(wc),
                                      "assigned_perWC" + suffix, outdir,
                                      data_root, theory, order, cov_full, drop_bin_idx,
                                      scan_min, scan_max, scan_n);
    results.push_back(r);
    assigned_result_by_key[wc] = r;
  }

  // Keep the explicitly scaled mut/d_t scans as named diagnostic plots,
  // but do not use them as a separate fit source. The summaries below use
  // assigned_result_by_key, so the WC summary and anomalous translation are consistent.
  FitResult1D r_mu = fit_one_wc_theory("ctGRe", obs_sets.at("AN22_028_Fig16_mu_t_1D"),
                                       "AN22_028_Fig16_mu_t_1D" + suffix, outdir,
                                       data_root, theory, order, cov_full, drop_bin_idx,
                                       scan_min, scan_max, scan_n, mu_scale, "#hat{#mu}_{t}");
  results.push_back(r_mu);

  FitResult1D r_dt = fit_one_wc_theory("ctGIm", obs_sets.at("AN22_028_Fig18_mu_t_vs_d_t"),
                                       "AN22_028_d_t_CPodd_1D" + suffix, outdir,
                                       data_root, theory, order, cov_full, drop_bin_idx,
                                       scan_min, scan_max, scan_n, mu_scale, "#hat{d}_{t}");
  results.push_back(r_dt);

  std::vector<SummaryEntry> wc_summary;
  for (const auto& wc : wc_list) {
    if (assigned_result_by_key.count(wc)) wc_summary.push_back(make_summary_entry(wc, pretty_label(wc), assigned_result_by_key[wc], 1.0));
  }
  const std::vector<PubEntry> top22_wc_pub = top22_wc_pub_entries();
  save_summary_plot(wc_summary,
                    outdir + "/summary_SMEFTsim_WCs" + suffix + ".pdf",
                    outdir + "/summary_SMEFTsim_WCs" + suffix + ".png",
                    "Wilson coefficient / #Lambda^{2} [TeV^{-2}]", -10.0, 10.0,
                    "Theory " + order + " Work in Progress", top22_wc_pub);
  write_summary_csv(wc_summary, outdir + "/summary_SMEFTsim_WCs" + suffix + ".csv");

  const double MT = 0.1725;
  const double GS = 1.1666;
  const double norm = MT * MT / (GS * GS);
  std::vector<SummaryEntry> anom_summary;
  if (assigned_result_by_key.count("ctGRe")) anom_summary.push_back(make_summary_entry("mut", "#hat{#mu}_{t}", assigned_result_by_key["ctGRe"], 2.0*MT*MT));
  if (assigned_result_by_key.count("ctGIm")) anom_summary.push_back(make_summary_entry("#hat{d}_{t}", "#hat{d}_{t}", assigned_result_by_key["ctGIm"], 2.0*MT*MT));
  add_linear_anom_entry(anom_summary, assigned_result_by_key, "#hat{c}_{VV}", "#hat{c}_{VV}", {{"ctj8",0.5},{"cQj18",0.5},{"ctu8",0.25},{"ctd8",0.25},{"cQu8",0.25},{"cQd8",0.25}}, norm);
  add_linear_anom_entry(anom_summary, assigned_result_by_key, "#hat{c}_{VA}", "#hat{c}_{VA}", {{"ctj8",0.5},{"cQj18",-0.5},{"ctu8",0.25},{"ctd8",0.25},{"cQu8",-0.25},{"cQd8",-0.25}}, norm);
  add_linear_anom_entry(anom_summary, assigned_result_by_key, "#hat{c}_{AV}", "#hat{c}_{AV}", {{"ctj8",-0.5},{"cQj18",-0.5},{"ctu8",0.25},{"ctd8",0.25},{"cQu8",0.25},{"cQd8",0.25}}, norm);
  add_linear_anom_entry(anom_summary, assigned_result_by_key, "#hat{c}_{AA}", "#hat{c}_{AA}", {{"ctj8",-0.5},{"cQj18",0.5},{"ctu8",0.25},{"ctd8",0.25},{"cQu8",-0.25},{"cQd8",-0.25}}, norm);
  add_linear_anom_entry(anom_summary, assigned_result_by_key, "#hat{c}_{1}", "#hat{c}_{1}", {{"ctu8",0.5},{"ctd8",-0.5},{"cQu8",0.5},{"cQd8",-0.5},{"cQj38",1.0}}, norm);
  add_linear_anom_entry(anom_summary, assigned_result_by_key, "#hat{c}_{3}", "#hat{c}_{3}", {{"ctu8",0.5},{"ctd8",-0.5},{"cQu8",-0.5},{"cQd8",0.5},{"cQj38",-1.0}}, norm);
  add_linear_anom_entry(anom_summary, assigned_result_by_key, "c1_minus_c2_plus_c3", "#hat{c}_{1}-#hat{c}_{2}+#hat{c}_{3}", {{"ctu8",0.5},{"ctd8",-0.5},{"cQu8",0.5},{"cQd8",-0.5},{"cQj38",-1.0}}, norm);

  std::vector<PubEntry> cms_pub = {
    {"mut", -0.005, 0.005}, {"#hat{d}_{t}", -0.004, 0.008}, {"#hat{c}_{VV}", 0.016, 0.013},
    {"#hat{c}_{VA}", -0.009, 0.018}, {"#hat{c}_{AV}", -0.001, 0.017}, {"#hat{c}_{AA}", 0.000, 0.020},
    {"#hat{c}_{1}", 0.13, 0.11}, {"#hat{c}_{3}", -0.07, 0.14}, {"c1_minus_c2_plus_c3", -0.01, 0.08}
  };
  save_summary_plot(anom_summary,
                    outdir + "/summary_anomalous_couplings" + suffix + ".pdf",
                    outdir + "/summary_anomalous_couplings" + suffix + ".png",
                    "Anomalous coupling", -0.4, 0.4,
                    "Theory " + order + " Work in Progress", cms_pub);
  write_summary_csv(anom_summary, outdir + "/summary_anomalous_couplings" + suffix + ".csv");

  for (const auto& spec : fig18_pairs) {
    fit_2d_pair_grid_theory(spec.wc1, spec.wc2, spec.obs, spec.tag + suffix,
                            outdir, data_root, theory, order, cov_full,
                            drop_bin_idx, spec.xmin, spec.xmax, scan2d_n);
  }

  std::ofstream summary(outdir + "/summary_1D" + suffix + ".csv");
  summary << "tag_or_group,wc,theory_order,best,lo68,hi68,lo95,hi95,chi2min,nbins\n";
  for (const auto& r : results) {
    summary << "theory," << r.name << "," << order << "," << r.best << "," << r.lo68 << "," << r.hi68 << ","
            << r.lo95 << "," << r.hi95 << "," << r.chi2min << "," << r.nbins << "\n";
  }
  summary.close();
}


int main() {
  // IMPORTANT:
  // This must be the unfolded DATA reference, not the MC/Nominal template.
  // The previous version accidentally pointed data_root to concatenated_histograms_mc.root,
  // so the black "Data" points and coefficient errors were actually MC-stat errors.
  // You can override this at runtime with:
  //   DATA_ROOT=/path/to/concatenated_histograms_data.root ./execMacro.sh ...
  const char* env_data_root = std::getenv("DATA_ROOT");
  const std::string data_root = (env_data_root && std::string(env_data_root).size())
    ? std::string(env_data_root)
    : std::string("/depot/cms/top/he614/notebooks/EFT_FullRun2/histogram_output_nanogen_ttbbllnunu_run0_test/concatenated_histograms_data.root");

  if (data_root.find("concatenated_histograms_mc.root") != std::string::npos ||
      data_root.find("concatenated_histograms_Nominal.root") != std::string::npos) {
    std::cerr << "[FATAL] DATA_ROOT points to an MC/Nominal template, not unfolded data: "
              << data_root << std::endl;
    std::cerr << "        Build/use concatenated_histograms_data.root from the gigantic-matrix central values "
              << "or set DATA_ROOT explicitly." << std::endl;
    return 2;
  }

  std::cout << "[INPUT] DATA_ROOT = " << data_root << std::endl;

  const std::string eft_template_pattern =
    "/depot/cms/top/he614/notebooks/EFT_FullRun2/histogram_output_nanogen_ttbbllnunu_run0_test/concatenated_histograms_{wc}_{val}.root";

  // const std::string cov_stat =
  //   "/depot/cms/top/dawoodo/fullRun2_UL_September2024_unfolding/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/stat_gigantic_matrix_fullRun2.root";

  // const std::string cov_syst =
  //   "/depot/cms/top/dawoodo/fullRun2_UL_September2024_unfolding/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/gigantic_matrices/syst_gigantic_matrix_fullRun2.root";
    
  const std::string cov_stat =
    "/depot/cms/top/dawoodo/fullRun2_UL_September2024_unfolding/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/gigantic_matrices/stat_gigantic_matrix_fullRun2.root";

  const std::string cov_syst =
    "/depot/cms/top/dawoodo/fullRun2_UL_September2024_unfolding/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/gigantic_matrices/syst_gigantic_matrix_fullRun2.root";

  const int drop_bin_idx = 1;
  const double scan_min = -20.0;
  const double scan_max = 20.0;
  const int scan_n = 10000;
  const int scan2d_n = 121;

  const std::string outdir = "nanogen_fits_root_ttbbllnunu_run0_test";
  gSystem->mkdir(outdir.c_str(), true);

  TheoryTable theory_table = load_embedded_theory_csv();
  bool have_theory = !theory_table.empty();
  std::cout << "[INFO theory] embedded theory comparison "
            << (have_theory ? "enabled" : "disabled") << std::endl;

  // Optional audit dump: keep the exact embedded table next to the plots.
  {
    std::ofstream theory_dump(outdir + "/theory_coefficients_embedded.csv");
    theory_dump << EMBEDDED_THEORY_CSV;
  }

  TMatrixD cov_full = load_covariance(cov_stat, cov_syst);
  dump_covariance_inspection_plots(cov_full, outdir, drop_bin_idx);

  const std::vector<std::string> wc_list = {
    "ctGRe", "ctGIm",
    "cQj11", "cQj31", "cQj18", "cQj38",
    "cQu1", "cQu8", "cQd1", "cQd8",
    "ctu1", "ctu8", "ctd1", "ctd8",
    "ctj1", "ctj8"
  };

  // Observable index convention follows your 0..37 screenshot and NanoGEN concatenation.
  // Known anchors from your screenshot: 34=cHel, 25=csca, 0=b1k, 2=b1r, 4=b1n, 10=ckk, 12=cnn.
  // AN-22-028 Fig.18 CP-odd cross terms use the standard AN order here: 18=cnr-crn, 20=cnk-ckn.
  // If your local concatenation order differs, only edit these constants.
  const int OBS_b1k = 0, OBS_b1r = 2, OBS_b1n = 4;
  const int OBS_ckk = 10, OBS_cnn = 12;
  const int OBS_cnrM = 18, OBS_cnkM = 20;
  const int OBS_csca = 25, OBS_cHel = 34;

  std::map<std::string, std::vector<int>> obs_sets;
  obs_sets["AN22_028_Fig16_mu_t_1D"] = {OBS_cHel, OBS_csca, OBS_b1k, OBS_ckk};
  obs_sets["AN22_028_Fig18_mu_t_vs_d_t"] = {OBS_cHel, OBS_csca, OBS_cnrM, OBS_cnkM};
  obs_sets["AN22_028_Fig18_mu_t_vs_cVV"] = {OBS_cHel, OBS_csca, OBS_ckk, OBS_b1r};
  obs_sets["AN22_028_Fig18_d_t_vs_cMinusMinus"] = {OBS_cnrM, OBS_cnkM, OBS_b1n, OBS_b1k};
  obs_sets["AN22_028_Fig18_cVV_vs_c1"] = {OBS_cHel, OBS_cnn, OBS_ckk, OBS_b1r};
  obs_sets["all_0_35"] = range_obs(0, 35);

  std::vector<FitResult1D> results;
  // AN-style anomalous ctG directions: mut = 2 m_t^2 ctGRe and d_t = 2 m_t^2 ctGIm, mt=0.1725 TeV.
  const double mt = 0.1725;
  const double mu_scale = 2.0 * mt * mt;

  // --- Assigned observable choices used by *all* 1D fits and summaries. ---
  // WCs listed here use the AN-22-028-motivated subset.
  // Any WC not listed here uses the assigned fallback: all_0_35.
  std::map<std::string, std::vector<int>> FIG16_OBS = {
    {"ctGRe", {OBS_cHel, OBS_csca, OBS_b1k, OBS_ckk}},
    {"ctGIm", {OBS_cHel, OBS_csca, OBS_cnrM, OBS_cnkM}},
    {"cQj18", {OBS_cHel, OBS_csca, OBS_ckk, OBS_b1r}},   // cVV proxy
    {"ctj8",  {OBS_cnrM, OBS_cnkM, OBS_b1n, OBS_b1k}},   // c-- proxy
    {"cQj38", {OBS_cHel, OBS_cnn, OBS_ckk, OBS_b1r}}     // c1 proxy
  };

  auto assigned_obs_for_wc = [&](const std::string& wc) -> std::vector<int> {
    auto it = FIG16_OBS.find(wc);
    if (it != FIG16_OBS.end()) return it->second;
    return obs_sets["all_0_35"];
  };

  // Single source of truth: the 16-WC summary and the anomalous-coupling
  // translations both consume these same assigned-observable fits.
  std::map<std::string, FitResult1D> assigned_result_by_key;

  for (const auto& wc : wc_list) {
    FitResult1D r = fit_one_wc(wc, assigned_obs_for_wc(wc),
                               "assigned_perWC", outdir,
                               data_root, eft_template_pattern, cov_full,
                               drop_bin_idx, scan_min, scan_max, scan_n);

    results.push_back(r);
    assigned_result_by_key[wc] = r;
  }

  // Keep the explicitly scaled mut/d_t scans as named diagnostic plots only.
  // The summary values are produced from assigned_result_by_key below, so there
  // is no second/competing fit definition.
  FitResult1D r_mu = fit_one_wc("ctGRe", obs_sets["AN22_028_Fig16_mu_t_1D"],
                                "AN22_028_Fig16_mu_t_1D", outdir,
                                data_root, eft_template_pattern, cov_full,
                                drop_bin_idx, scan_min, scan_max, scan_n,
                                mu_scale, "#hat{#mu}_{t}");
  results.push_back(r_mu);

  FitResult1D r_dt = fit_one_wc("ctGIm", obs_sets["AN22_028_Fig18_mu_t_vs_d_t"],
                                "AN22_028_d_t_CPodd_1D", outdir,
                                data_root, eft_template_pattern, cov_full,
                                drop_bin_idx, scan_min, scan_max, scan_n,
                                mu_scale, "#hat{d}_{t}");
  results.push_back(r_dt);

  std::vector<SummaryEntry> wc_summary;
  for (const auto& wc : wc_list) {
    if (assigned_result_by_key.count(wc)) wc_summary.push_back(make_summary_entry(wc, pretty_label(wc), assigned_result_by_key[wc], 1.0));
  }
  const std::vector<PubEntry> top22_wc_pub = top22_wc_pub_entries();
  save_summary_plot(wc_summary, outdir + "/summary_SMEFTsim_WCs.pdf", outdir + "/summary_SMEFTsim_WCs.png",
                    "Wilson coefficient / #Lambda^{2} [TeV^{-2}]", -10.0, 10.0,
                    "Simulation Work in Progress", top22_wc_pub);
  write_summary_csv(wc_summary, outdir + "/summary_SMEFTsim_WCs.csv");
  std::vector<SummaryEntry> anom_summary;


  auto add_linear_anom = [&anom_summary, &assigned_result_by_key](const std::string& key, const std::string& label,
                            const std::vector<std::pair<std::string,double>>& terms,
                            double scale) {
    double best = 0.0, varlo = 0.0, varhi = 0.0;
    bool ok = true;
    for (const auto& kv : terms) {
      if (!assigned_result_by_key.count(kv.first)) { ok = false; break; }
      const auto& r = assigned_result_by_key[kv.first];
      best += kv.second * r.best;
      const double elo = std::max(0.0, r.best - r.lo68);
      const double ehi = std::max(0.0, r.hi68 - r.best);
      varlo += (kv.second * elo) * (kv.second * elo);
      varhi += (kv.second * ehi) * (kv.second * ehi);
    }
    if (!ok) { std::cout << "[SKIP anom] missing WC input for " << key << std::endl; return; }
    SummaryEntry e;
    e.key = key; e.label = label;
    e.best = best * scale;
    e.lo68 = (best - std::sqrt(varlo)) * scale;
    e.hi68 = (best + std::sqrt(varhi)) * scale;
    anom_summary.push_back(e);
  };

  const double MT = 0.1725;
  const double GS = 1.1666;
  const double norm = MT*MT/(GS*GS);
  anom_summary.push_back(make_summary_entry("mut", "#hat{#mu}_{t}", assigned_result_by_key["ctGRe"], 2.0*MT*MT));
  anom_summary.push_back(make_summary_entry("#hat{d}_{t}", "#hat{d}_{t}", assigned_result_by_key["ctGIm"], 2.0*MT*MT));
  add_linear_anom("#hat{c}_{VV}", "#hat{c}_{VV}", {{"ctj8",0.5},{"cQj18",0.5},{"ctu8",0.25},{"ctd8",0.25},{"cQu8",0.25},{"cQd8",0.25}}, norm);
  add_linear_anom("#hat{c}_{VA}", "#hat{c}_{VA}", {{"ctj8",0.5},{"cQj18",-0.5},{"ctu8",0.25},{"ctd8",0.25},{"cQu8",-0.25},{"cQd8",-0.25}}, norm);
  add_linear_anom("#hat{c}_{AV}", "#hat{c}_{AV}", {{"ctj8",-0.5},{"cQj18",-0.5},{"ctu8",0.25},{"ctd8",0.25},{"cQu8",0.25},{"cQd8",0.25}}, norm);
  add_linear_anom("#hat{c}_{AA}", "#hat{c}_{AA}", {{"ctj8",-0.5},{"cQj18",0.5},{"ctu8",0.25},{"ctd8",0.25},{"cQu8",-0.25},{"cQd8",-0.25}}, norm);
  add_linear_anom("#hat{c}_{1}", "#hat{c}_{1}", {{"ctu8",0.5},{"ctd8",-0.5},{"cQu8",0.5},{"cQd8",-0.5},{"cQj38",1.0}}, norm);
  add_linear_anom("#hat{c}_{3}", "#hat{c}_{3}", {{"ctu8",0.5},{"ctd8",-0.5},{"cQu8",-0.5},{"cQd8",0.5},{"cQj38",-1.0}}, norm);
  add_linear_anom("c1_minus_c2_plus_c3", "#hat{c}_{1}-#hat{c}_{2}+#hat{c}_{3}", {{"ctu8",0.5},{"ctd8",-0.5},{"cQu8",0.5},{"cQd8",-0.5},{"cQj38",-1.0}}, norm);
  std::vector<PubEntry> cms_pub = {
    {"mut", -0.005, 0.005},
    {"#hat{d}_{t}", -0.004, 0.008},
    {"#hat{c}_{VV}", 0.016, 0.013},
    {"#hat{c}_{VA}", -0.009, 0.018},
    {"#hat{c}_{AV}", -0.001, 0.017},
    {"#hat{c}_{AA}", 0.000, 0.020},
    {"#hat{c}_{1}", 0.13, 0.11},
    {"#hat{c}_{3}", -0.07, 0.14},
    {"c1_minus_c2_plus_c3", -0.01, 0.08}
  };
  save_summary_plot(anom_summary, outdir + "/summary_anomalous_couplings.pdf", outdir + "/summary_anomalous_couplings.png",
                    "Anomalous coupling", -0.4, 0.4,
                    "Work in Progress", cms_pub);
  write_summary_csv(anom_summary, outdir + "/summary_anomalous_couplings.csv");

  // AN-22-028 Fig.18 2D observable choices:
  //   mut,d_t: chel,csca,cnr-crn,cnk-ckn
  //   mut,cVV: chel,csca,ckk,b1r
  //   d_t,c--:  cnr-crn,cnk-ckn,b1n,b1k
  //   cVV,c1:   chel,cnn,ckk,b1r
  // cVV/c1/c-- below are NanoGEN WC-template proxies. Replace proxy WCs with direct experiment-basis templates if available.
  std::vector<PairFitSpec> fig18_pairs = {
    {"ctGRe", "ctGIm", "AN22_028_Fig18_mu_t_vs_d_t", obs_sets["AN22_028_Fig18_mu_t_vs_d_t"], -8.0, 8.0},
    {"ctGRe", "cQj18", "AN22_028_Fig18_ctGRe_vs_cQj18_cVV_proxy", obs_sets["AN22_028_Fig18_mu_t_vs_cVV"], -8.0, 8.0},
    {"ctGIm", "ctj8", "AN22_028_Fig18_ctGIm_vs_ctj8_cMinusMinus_proxy", obs_sets["AN22_028_Fig18_d_t_vs_cMinusMinus"], -8.0, 8.0},
    {"cQj18", "cQj38", "AN22_028_Fig18_cQj18_cVV_proxy_vs_cQj38_c1_proxy", obs_sets["AN22_028_Fig18_cVV_vs_c1"], -8.0, 8.0}
  };
  for (const auto& spec : fig18_pairs) {
    fit_2d_pair_grid_if_available(spec, outdir, data_root, eft_template_pattern, cov_full, drop_bin_idx, scan2d_n);
  }

  if (have_theory) {
    run_theory_fit_suite("LO", "_theory_lo", outdir, data_root, theory_table, cov_full,
                         drop_bin_idx, scan_min, scan_max, scan_n, scan2d_n,
                         wc_list, obs_sets, FIG16_OBS, fig18_pairs);
    run_theory_fit_suite("NLO", "_theory_nlo", outdir, data_root, theory_table, cov_full,
                         drop_bin_idx, scan_min, scan_max, scan_n, scan2d_n,
                         wc_list, obs_sets, FIG16_OBS, fig18_pairs);
  }

  for (const auto& wc : wc_list) {
    make_concatenated_inspection_plot_root(
      data_root, eft_template_pattern, wc, outdir,
      have_theory ? &theory_table : nullptr, "LO"
    );
    if (have_theory) {
      make_concatenated_inspection_plot_root(
        data_root, eft_template_pattern, wc, outdir,
        &theory_table, "NLO"
      );
    }
    make_individual_distribution_plots_and_coefficients_root(
      data_root, eft_template_pattern, wc, outdir,
      have_theory ? &theory_table : nullptr, "LO"
    );
  }

  std::ofstream summary(outdir + "/summary_1D.csv");
  summary << "tag_or_group,wc,best,lo68,hi68,lo95,hi95,chi2min,nbins\n";
  for (const auto& r : results) {
    summary << "mixed," << r.name << "," << r.best << "," << r.lo68 << "," << r.hi68 << ","
            << r.lo95 << "," << r.hi95 << "," << r.chi2min << "," << r.nbins << "\n";
  }
  summary.close();

  std::cout << "\nDone. Outputs in " << outdir << std::endl;
  return 0;
}
