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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static const int BINS_PER_OBS = 6;
static const int N_OBS_TOTAL = 38;

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
  if (wc == "ctGRe") return "c_{tG}";
  if (wc == "ctGIm") return "c^{I}_{tG}";
  if (wc == "mu_t") return "#hat{#mu}_{t}";
  if (wc == "d_t") return "#hat{d}_{t}";
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
  c.SetLeftMargin(0.13);
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

  TLegend leg(0.66, 0.68, 0.93, 0.91);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
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
  c.SetLeftMargin(0.13);
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
  c.SetLeftMargin(0.30); c.SetRightMargin(0.04); c.SetTopMargin(0.16); c.SetBottomMargin(0.13);
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
  TLegend leg(0.58, 0.70, 0.94, 0.84); leg.SetBorderSize(0); leg.SetFillStyle(0);
  TLine dummySM; dummySM.SetLineColor(green); dummySM.SetLineWidth(3);
  TGraph dummyPub; dummyPub.SetMarkerStyle(24); dummyPub.SetMarkerColor(TColor::GetColor("#5b7fd9")); dummyPub.SetLineColor(TColor::GetColor("#5b7fd9")); dummyPub.SetLineWidth(2);
  TGraph dummyFit; dummyFit.SetMarkerStyle(20); dummyFit.SetMarkerColor(kBlack); dummyFit.SetLineColor(kBlack); dummyFit.SetLineWidth(2);
  leg.AddEntry(&dummySM, "Standard model", "l");
  if (!pub_entries.empty()) leg.AddEntry(&dummyPub, "TOP-18-006 (68% CL)", "lep");
  leg.AddEntry(&dummyFit, "This fit (68% CL)", "lep");
  leg.Draw();
  TLatex latex; latex.SetNDC(); latex.SetTextFont(42); latex.SetTextSize(0.034);
  latex.DrawLatex(0.30, 0.965, "#bf{CMS}");
  latex.SetTextSize(0.028);
  latex.DrawLatex(0.43, 0.965, ("#it{" + cms_extra + "}").c_str());
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

int main() {
  const std::string data_root =
    "/depot/cms/top/he614/notebooks/EFT_FullRun2/histogram_output_nanogen/concatenated_histograms_data.root";

  const std::string eft_template_pattern =
    "/depot/cms/top/he614/notebooks/EFT_FullRun2/histogram_output_nanogen/concatenated_histograms_{wc}_{val}.root";

  const std::string cov_stat =
    "/depot/cms/top/dawoodo/fullRun2_UL_September2024_unfolding/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/gigantic_matrices/stat_gigantic_matrix_fullRun2.root";

  const std::string cov_syst =
    "/depot/cms/top/dawoodo/fullRun2_UL_September2024_unfolding/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/gigantic_matrices/syst_gigantic_matrix_fullRun2.root";

  const int drop_bin_idx = 1;
  const double scan_min = -20.0;
  const double scan_max = 20.0;
  const int scan_n = 10000;
  const int scan2d_n = 121;

  const std::string outdir = "nanogen_fits_root";
  gSystem->mkdir(outdir.c_str(), true);

  TMatrixD cov_full = load_covariance(cov_stat, cov_syst);

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
  std::map<std::string, FitResult1D> result_by_key;

  // AN-style anomalous ctG directions: mu_t = 2 m_t^2 ctGRe and d_t = 2 m_t^2 ctGIm, mt=0.1725 TeV.
  const double mt = 0.1725;
  const double mu_scale = 2.0 * mt * mt;

  FitResult1D r_mu = fit_one_wc("ctGRe", obs_sets["AN22_028_Fig16_mu_t_1D"],
                                "AN22_028_Fig16_mu_t_1D", outdir,
                                data_root, eft_template_pattern, cov_full,
                                drop_bin_idx, scan_min, scan_max, scan_n,
                                mu_scale, "#hat{#mu}_{t}");
  results.push_back(r_mu); result_by_key["mu_t"] = r_mu;

  FitResult1D r_dt = fit_one_wc("ctGIm", obs_sets["AN22_028_Fig18_mu_t_vs_d_t"],
                                "AN22_028_d_t_CPodd_1D", outdir,
                                data_root, eft_template_pattern, cov_full,
                                drop_bin_idx, scan_min, scan_max, scan_n,
                                mu_scale, "#hat{d}_{t}");
  results.push_back(r_dt); result_by_key["d_t"] = r_dt;

  // --- Fig.16 per-WC observable fits ---
  std::map<std::string, std::vector<int>> FIG16_OBS = {
    {"ctGRe", {OBS_cHel, OBS_csca, OBS_b1k, OBS_ckk}},
    {"ctGIm", {OBS_cnrM, OBS_cnkM, OBS_b1n, OBS_b1k}},
    {"cQj18", {OBS_cHel, OBS_cnn, OBS_ckk, OBS_b1r}},   // cVV proxy
    {"ctj8",  {OBS_cnrM, OBS_cnkM, OBS_b1n, OBS_b1k}},  // c-- proxy
    {"cQj38", {OBS_cHel, OBS_cnn, OBS_ckk, OBS_b1r}}    // c1 proxy
  };

  for (const auto& wc : wc_list) {
    std::vector<int> obs = obs_sets["all_0_35"];
    if (FIG16_OBS.count(wc)) obs = FIG16_OBS[wc];

    FitResult1D r = fit_one_wc(wc, obs,
                               "FIG16_perWC", outdir,
                               data_root, eft_template_pattern, cov_full,
                               drop_bin_idx, scan_min, scan_max, scan_n);

    results.push_back(r);
    result_by_key[wc] = r;
  }

  std::vector<SummaryEntry> wc_summary;
  for (const auto& wc : wc_list) {
    if (result_by_key.count(wc)) wc_summary.push_back(make_summary_entry(wc, pretty_label(wc), result_by_key[wc], 1.0));
  }
  save_summary_plot(wc_summary, outdir + "/summary_SMEFTsim_WCs.pdf", outdir + "/summary_SMEFTsim_WCs.png",
                    "Wilson coefficient / #Lambda^{2} [TeV^{-2}]", -10.0, 10.0);
  write_summary_csv(wc_summary, outdir + "/summary_SMEFTsim_WCs.csv");
  std::vector<SummaryEntry> anom_summary;


  auto add_linear_anom = [&anom_summary, &result_by_key](const std::string& key, const std::string& label,
                            const std::vector<std::pair<std::string,double>>& terms,
                            double scale) {
    double best = 0.0, varlo = 0.0, varhi = 0.0;
    bool ok = true;
    for (const auto& kv : terms) {
      if (!result_by_key.count(kv.first)) { ok = false; break; }
      const auto& r = result_by_key[kv.first];
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
  anom_summary.push_back(make_summary_entry("mu_t", "#hat{#mu}_{t}", result_by_key["ctGRe"], 2.0*MT*MT));
  anom_summary.push_back(make_summary_entry("d_t", "#hat{d}_{t}", result_by_key["ctGIm"], 2.0*MT*MT));
  add_linear_anom("cVV", "#hat{c}_{VV}", {{"ctj8",0.5},{"cQj18",0.5},{"ctu8",0.25},{"ctd8",0.25},{"cQu8",0.25},{"cQd8",0.25}}, norm);
  add_linear_anom("cVA", "#hat{c}_{VA}", {{"ctj8",0.5},{"cQj18",-0.5},{"ctu8",0.25},{"ctd8",0.25},{"cQu8",-0.25},{"cQd8",-0.25}}, norm);
  add_linear_anom("cAV", "#hat{c}_{AV}", {{"ctj8",-0.5},{"cQj18",-0.5},{"ctu8",0.25},{"ctd8",0.25},{"cQu8",0.25},{"cQd8",0.25}}, norm);
  add_linear_anom("cAA", "#hat{c}_{AA}", {{"ctj8",-0.5},{"cQj18",0.5},{"ctu8",0.25},{"ctd8",0.25},{"cQu8",-0.25},{"cQd8",-0.25}}, norm);
  add_linear_anom("c1", "#hat{c}_{1}", {{"ctu8",0.5},{"ctd8",-0.5},{"cQu8",0.5},{"cQd8",-0.5},{"cQj38",1.0}}, norm);
  add_linear_anom("c3", "#hat{c}_{3}", {{"ctu8",0.5},{"ctd8",-0.5},{"cQu8",-0.5},{"cQd8",0.5},{"cQj38",-1.0}}, norm);
  add_linear_anom("c1_minus_c2_plus_c3", "#hat{c}_{1}-#hat{c}_{2}+#hat{c}_{3}", {{"ctu8",0.5},{"ctd8",-0.5},{"cQu8",0.5},{"cQd8",-0.5},{"cQj38",-1.0}}, norm);
  std::vector<PubEntry> cms_pub = {
    {"mu_t", -0.005, 0.005},
    {"d_t", -0.004, 0.008},
    {"cVV", 0.016, 0.013},
    {"cVA", -0.009, 0.018},
    {"cAV", -0.001, 0.017},
    {"cAA", 0.000, 0.020},
    {"c1", 0.13, 0.11},
    {"c3", -0.07, 0.14},
    {"c1_minus_c2_plus_c3", -0.01, 0.08}
  };
  save_summary_plot(anom_summary, outdir + "/summary_anomalous_couplings.pdf", outdir + "/summary_anomalous_couplings.png",
                    "Anomalous coupling", -0.4, 0.4,
                    "Work in Progress", cms_pub);
  write_summary_csv(anom_summary, outdir + "/summary_anomalous_couplings.csv");

  // AN-22-028 Fig.18 2D observable choices:
  //   mu_t,d_t: chel,csca,cnr-crn,cnk-ckn
  //   mu_t,cVV: chel,csca,ckk,b1r
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
