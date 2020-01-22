// -*- C++ -*-
#include "../src/PlotUtil.h"

// root -l -b wbern_input.cxx++

void addHistogram(std::vector<std::pair<std::string, std::vector<std::unique_ptr<TH1D>>>>::const_iterator iOp, 
                  const std::vector<std::pair<std::string, std::vector<std::unique_ptr<TH1D>>>> &vp_hist,
                  const std::vector<std::string> v_syst, std::unique_ptr<TFile> &file,
                  std::vector<std::pair<std::string, uint>> vp_index, const double &iAdd = 0.)
{
  const uint nHist = iOp->second.size();

  // terminating condition - write to file, reset to null
  if (std::next(iOp) == std::end(vp_hist)) {
    for (uint iH = 0; iH < iOp->second.size(); ++iH) {
      auto vp_iH = vp_index;
      vp_iH.push_back({iOp->first, iH});
      //printAll(vp_iH);

      // actually make the histogram
      const std::string sAdd = (iAdd > 0.) ? "_plus_" : "_minus_";
      std::string name = vp_iH.at(0).first + "_" + v_syst.at(vp_iH.at(0).second);

      // sm is of course only added
      auto hist = std::make_unique<TH1D>("xxx", "", 22, 0., 22.);
      auto iHist0 = std::find_if(std::begin(vp_hist), std::end(vp_hist), [&vp_iH] (auto &p) { return p.first == vp_iH.at(0).first; });
      hist->Add(iHist0->second.at(vp_iH.at(0).second).get(), 1.);

      for (uint iS = 1; iS < vp_iH.size(); ++iS) {
        name = name + sAdd + vp_iH.at(iS).first + "_" + v_syst.at(vp_iH.at(iS).second);

        auto iHist = std::find_if(std::begin(vp_hist), std::end(vp_hist), [&vp_iH, &iS] (auto &p) { return p.first == vp_iH.at(iS).first; });
        hist->Add(hist.get(), iHist->second.at(vp_iH.at(iS).second).get(), 1., iAdd);
      }

      hist->SetName(name.c_str());
      file->cd();
      hist->Write();
      //std::cout << name << std::endl;
    }
  }
  else {
    for (uint iH = 0; iH < iOp->second.size(); ++iH) {
      auto vp_iH = vp_index;
      vp_iH.push_back({iOp->first, iH});
      addHistogram(std::next(iOp), vp_hist, v_syst, file, vp_iH, iAdd);
    }
  }
}

void produceFile(const std::map<std::pair<std::string, std::string>, std::vector<double>> &m_coeff,
                 const std::vector<std::string> &v_opStr = {}) {
  // it should be 22
  for (const auto &pair : m_coeff) {
    if (pair.second.size() != 22) {
      std::cout << "Ok that's weird, somehow there's an element in coeff map that doesn't have 22 elements..." << std::endl;
      std::cout << "First faulty element is " << pair.first << std::endl;
      return;
    }
  }

  if (v_opStr.empty()) return;

  std::string opStr = v_opStr.at(0);
  for (uint iOp = 1; iOp < v_opStr.size(); ++iOp) {
    opStr = opStr + "_" + v_opStr.at(iOp);
  }

  const std::string outName = "/nfs/dust/cms/user/afiqaize/cms/rand/eftRivet_290118/EFTFitter/wbern_0314/root/" + opStr + "_coeff.root";
  const std::vector<std::string> v_syst = {"nominal", "up", "down"};

  // to dump all the histograms - vec<pair> rather than map to keep push order
  std::vector<std::pair<std::string, std::vector<std::unique_ptr<TH1D>>>> vp_hist;
  vp_hist.emplace_back( std::make_pair("sm", std::vector<std::unique_ptr<TH1D>>()) );
  for (const auto &op : v_opStr)
    vp_hist.emplace_back( std::make_pair(op, std::vector<std::unique_ptr<TH1D>>()) );

  // first make the 'pure' histogram - sm and op - separately nominal, up, down
  // sm is sm, for op just make the +1 case (without sm contribution at this stage)
  for (auto &p_hist : vp_hist) {
    auto &v_hist = p_hist.second;

    for (const std::string &syst : v_syst) {
      v_hist.emplace_back(std::make_unique<TH1D>((p_hist.first + "_" + syst).c_str(), "", 22, 0., 22.));

      for (uint iC = 0; iC < std::begin(m_coeff)->second.size(); ++iC) {
        double coeff = m_coeff.at({p_hist.first, "nominal"}).at(iC);
        if (syst != "nominal")
          coeff += m_coeff.at({p_hist.first, syst}).at(iC);

        v_hist.back()->SetBinContent(iC + 1, coeff);
      }
    }
  }

  // create the file and save histograms
  auto file = std::make_unique<TFile>(outName.c_str(), "recreate");
  file->cd();

  // for pure cases it's straightforward, just save as-is for all syst
  for (auto &p_hist : vp_hist) {
    p_hist.second.at(0)->Write();
    p_hist.second.at(1)->Write();
    p_hist.second.at(2)->Write();
  }

  // for the op, well recursion is our (awkward and difficult) friend...
  std::unique_ptr<TH1D> hist = nullptr;
  addHistogram(std::begin(vp_hist), vp_hist, v_syst, file, {}, 1.);
  addHistogram(std::begin(vp_hist), vp_hist, v_syst, file, {}, -1.);
}

void wbern_input() {
  gROOT->Reset();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2(true);

  // put the information in (do it only up to CPMij, cHel and so on coded in)
  // NLO pred: table 3 TOP-18-006
  // nominal scale is mt: up is 2mt, down is mt/2
  // signs are as in WB's mail 12/03/2019 for SM (blj and blq assumed same as blk and blr)
  // EFT numbers as in WB's mail 15/03/2019
  std::map<std::pair<std::string, std::string>, std::vector<double>> m_coeff;
  m_coeff.insert({{"sm", "nominal"}, 
        {0.004, 0.004, 0.0016, 0.0016, 0.0057, 0.0057, 0., 0., 0., 0., 0.331, 0.071, 0.326, -0.206, 0., 0.00106, 0., 0.00215, 0.}});
  m_coeff.insert({{"sm", "up"}, 
        {0.0017, 0.0017, 0.0012, 0.0012, -0.0004, -0.0004, 0.0005, 0.0005, 0.0005, 0.0005, 0.002, -0.006, 0.002, -0.002, 0., -0.00001, 0., -0.00007, 0.}});
  m_coeff.insert({{"sm", "down"}, 
        {-0.0012, -0.0012, -0.0009, -0.0009, 0.0005, 0.0005, -0.0005, -0.0005, -0.0005, -0.0005, -0.002, 0.008, -0.002, 0.002, 0., 0.00001, 0., 0.00004, 0.}});

  m_coeff.insert({{"ut", "nominal"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.917, 2.475, 2.025, 0.74, 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"ut", "up"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.006, -0.019, -0.024, 0.001, 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"ut", "down"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.006, 0.02, 0.025, -0.002, 0., 0., 0., 0., 0.}});

  m_coeff.insert({{"cvv", "nominal"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1.218, -0.697, -0.0799, -0.306, 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"cvv", "up"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.04, -0.027, -0.0085, -0.014, 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"cvv", "down"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.039, 0.027, 0.008, 0.014, 0., 0., 0., 0., 0.}});

  m_coeff.insert({{"c1", "nominal"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.151, -0.0846, -0.00821, -0.0358, 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"c1", "up"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.007, -0.0043, -0.00108, -0.002, 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"c1", "down"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.007, 0.0041, 0.00099, 0.002, 0., 0., 0., 0., 0.}});

  m_coeff.insert({{"dt", "nominal"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -4.143, 0., -0.8}});
  m_coeff.insert({{"dt", "up"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.053, 0., 0.006}});
  m_coeff.insert({{"dt", "down"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.056, 0., -0.003}});

  m_coeff.insert({{"cmm", "nominal"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1.226, 0., -2.157}});
  m_coeff.insert({{"cmm", "up"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.004, 0., 0.108}});
  m_coeff.insert({{"cmm", "down"}, 
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.007, 0., -0.087}});

  // for the sums and differences of CP-even, split as B1 = B2
  m_coeff.insert({{"cva", "nominal"}, 
        {0.8035, 0.8035, 0.105, 0.105, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"cva", "up"}, 
        {0.0255, 0.0255, 0.0045, 0.0045, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"cva", "down"}, 
        {-0.026, -0.026, -0.0045, -0.0045, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}});

  m_coeff.insert({{"c3", "nominal"}, 
        {0.1005, 0.1005, 0.01275, 0.01275, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"c3", "up"}, 
        {0.0045, 0.0045, 0.0007, 0.0007, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"c3", "down"}, 
        {-0.0045, -0.0045, -0.0007, -0.0007, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}});

  m_coeff.insert({{"cav", "nominal"}, 
        {0., 0., 0., 0., 0., 0., 0.4415, 0.4415, 0.396, 0.396, 0., 0., 0., 0., 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"cav", "up"}, 
        {0., 0., 0., 0., 0., 0., 0.017, 0.017, 0.02, 0.02, 0., 0., 0., 0., 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"cav", "down"}, 
        {0., 0., 0., 0., 0., 0., -0.017, -0.017, -0.0195, -0.0195, 0., 0., 0., 0., 0., 0., 0., 0., 0.}});

  m_coeff.insert({{"c123", "nominal"}, 
        {0., 0., 0., 0., 0., 0., 0.0915, 0.0915, 0.082, 0.082, 0., 0., 0., 0., 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"c123", "up"}, 
        {0., 0., 0., 0., 0., 0., 0.004, 0.004, 0.0045, 0.0045, 0., 0., 0., 0., 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"c123", "down"}, 
        {0., 0., 0., 0., 0., 0., -0.004, -0.004, -0.0045, -0.0045, 0., 0., 0., 0., 0., 0., 0., 0., 0.}});

  // this last guy that is P-even but CP-odd, assume maximal CPV => B1 = -B2
  m_coeff.insert({{"cmp", "nominal"}, 
        {0., 0., 0., 0., 2.438, -2.438, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"cmp", "up"}, 
        {0., 0., 0., 0., -0.062, 0.062, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}});
  m_coeff.insert({{"cmp", "down"}, 
        {0., 0., 0., 0., 0.054, -0.054, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}});

  // fill up cHel and labs - labs needed only so snake is there, but it's never used
  // so wont bother with correlation in scales
  for (auto &p_coeff : m_coeff) {
    const double cHel = -1. * (p_coeff.second.at(10) + p_coeff.second.at(11) + p_coeff.second.at(12)) / 3.;
    p_coeff.second.push_back(cHel);

    if (p_coeff.first.first != "sm") {
      p_coeff.second.push_back(0.);
      p_coeff.second.push_back(0.);
    }
    else {
      if (p_coeff.first.second == "nominal") {
        p_coeff.second.push_back(0.181);
        p_coeff.second.push_back(0.108);
      }
      if (p_coeff.first.second == "up") {
        p_coeff.second.push_back(0.004);
        p_coeff.second.push_back(0.009);
      }
      if (p_coeff.first.second == "down") {
        p_coeff.second.push_back(-0.003);
        p_coeff.second.push_back(-0.012);
      }
    }
  }

  const std::vector<std::string> v_op = {"ut", "cvv", "c1", "dt", "cmm", "cva", "c3", "cav", "c123", "cmp"};
  //const std::vector<std::string> v_op = {"ut", "cvv"};

  // only go 2D since first: EFTFitter won't like the huge dimensions and second: these are the highest correlated ops (ignoring degeneracies)
  for (uint iOp0 = 0; iOp0 < v_op.size(); ++iOp0) {
    produceFile(m_coeff, {v_op.at(iOp0)});

    for (uint iOp1 = iOp0 + 1; iOp1 < v_op.size(); ++iOp1)
      produceFile(m_coeff, {v_op.at(iOp0), v_op.at(iOp1)});
  }

  gROOT->ProcessLine(".q");
}
