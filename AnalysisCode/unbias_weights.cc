// unbias_weights.cc
// Apply information-theoretic reweighting (Assi et al., arXiv:2501.17219)
// to unbias a single-observable spectrum (here: jet pT).
//
// Builds: c++ -O2 unbias_weights.cc $(root-config --cflags --libs) -o unbias_weights
//
// Example:
//   ./unbias_weights \
//     --input UnbiasingTest_PYTHIApp_pthatmin50_121825.root \
//     --tree tgenBefore --n-branch nJets --pt-branch pt --weight-branch weight \
//     --target-hist target.root:hpT \
//     --moments 1,2,3,4 --pt0 1.0 --pt-min 50 --pt-max 200 \
//     --out unbias_weights.root

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

static void die(const std::string &msg) {
  std::cerr << "error: " << msg << std::endl;
  std::exit(1);
}

static bool starts_with(const std::string &s, const std::string &pfx) {
  return s.rfind(pfx, 0) == 0;
}

static std::string get_arg(int argc, char* argv[], const std::string &flag, const std::string &def = "") {
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == flag && i + 1 < argc) return argv[i + 1];
  }
  return def;
}

static bool has_flag(int argc, char* argv[], const std::string &flag) {
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == flag) return true;
  }
  return false;
}

static std::vector<int> parse_moments(const std::string &s) {
  std::vector<int> out;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    if (tok.empty()) continue;
    out.push_back(std::stoi(tok));
  }
  if (out.empty()) die("empty --moments list");
  return out;
}

static void split_target(const std::string &s, std::string &file, std::string &hist) {
  auto pos = s.find(':');
  if (pos == std::string::npos) die("--target-hist must be file.root:histname");
  file = s.substr(0, pos);
  hist = s.substr(pos + 1);
  if (file.empty() || hist.empty()) die("--target-hist must be file.root:histname");
}

static double safe_log(double x, double x0) {
  if (x <= 0.0 || x0 <= 0.0) return std::numeric_limits<double>::quiet_NaN();
  return std::log(x / x0);
}

static bool solve_linear(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double> &x) {
  // Gaussian elimination with partial pivoting.
  const int n = static_cast<int>(b.size());
  x.assign(n, 0.0);
  for (int i = 0; i < n; ++i) {
    int piv = i;
    double maxv = std::fabs(A[i][i]);
    for (int r = i + 1; r < n; ++r) {
      double v = std::fabs(A[r][i]);
      if (v > maxv) { maxv = v; piv = r; }
    }
    if (maxv == 0.0) return false;
    if (piv != i) {
      std::swap(A[i], A[piv]);
      std::swap(b[i], b[piv]);
    }
    double diag = A[i][i];
    for (int c = i; c < n; ++c) A[i][c] /= diag;
    b[i] /= diag;
    for (int r = 0; r < n; ++r) {
      if (r == i) continue;
      double f = A[r][i];
      for (int c = i; c < n; ++c) A[r][c] -= f * A[i][c];
      b[r] -= f * b[i];
    }
  }
  for (int i = 0; i < n; ++i) x[i] = b[i];
  return true;
}

int main(int argc, char* argv[]) {
  const std::string input = get_arg(argc, argv, "--input");
  const std::string tree_name = get_arg(argc, argv, "--tree", "tgenBefore");
  const std::string n_branch = get_arg(argc, argv, "--n-branch", "nJets");
  const std::string pt_branch = get_arg(argc, argv, "--pt-branch", "pt");
  const std::string weight_branch = get_arg(argc, argv, "--weight-branch", "weight");
  const std::string target_spec = get_arg(argc, argv, "--target-hist");
  const std::string moments_spec = get_arg(argc, argv, "--moments", "1,2,3,4");
  const std::string out_name = get_arg(argc, argv, "--out", "unbias_weights.root");

  const double pt0 = std::stod(get_arg(argc, argv, "--pt0", "1.0"));
  const double pt_min = std::stod(get_arg(argc, argv, "--pt-min", "0.0"));
  const double pt_max = std::stod(get_arg(argc, argv, "--pt-max", "1e9"));

  const int max_iter = std::stoi(get_arg(argc, argv, "--max-iter", "200"));
  const double tol = std::stod(get_arg(argc, argv, "--tol", "1e-6"));
  const double cov_reg = std::stod(get_arg(argc, argv, "--cov-reg", "1e-12"));

  if (input.empty()) die("--input is required");
  if (target_spec.empty()) die("--target-hist is required");

  std::vector<int> moments = parse_moments(moments_spec);
  const int nm = static_cast<int>(moments.size());

  std::string target_file, target_hist;
  split_target(target_spec, target_file, target_hist);

  // Load target histogram and compute target moments c_j.
  TFile tf(target_file.c_str(), "READ");
  if (tf.IsZombie()) die("failed to open target file: " + target_file);
  TH1 *h = dynamic_cast<TH1*>(tf.Get(target_hist.c_str()));
  if (!h) die("failed to find target hist: " + target_hist);
  double norm = 0.0;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    norm += h->GetBinContent(b) * h->GetBinWidth(b);
  }
  if (norm <= 0.0) die("target hist has non-positive normalization");
  std::vector<double> c(nm, 0.0);
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    double x = h->GetBinCenter(b);
    if (x < pt_min || x > pt_max) continue;
    double w = h->GetBinContent(b) * h->GetBinWidth(b);
    double l = safe_log(x, pt0);
    if (!std::isfinite(l)) continue;
    for (int j = 0; j < nm; ++j) {
      c[j] += w * std::pow(l, moments[j]);
    }
  }
  for (int j = 0; j < nm; ++j) c[j] /= norm;
  tf.Close();

  // Load input tree and flatten jet pT list.
  TFile in(input.c_str(), "READ");
  if (in.IsZombie()) die("failed to open input file: " + input);
  TTree *tree = dynamic_cast<TTree*>(in.Get(tree_name.c_str()));
  if (!tree) die("failed to find tree: " + tree_name);

  int nJets = 0;
  float pt[100] = {0};
  float weight = 1.0f;

  tree->SetBranchAddress(n_branch.c_str(), &nJets);
  tree->SetBranchAddress(pt_branch.c_str(), pt);
  tree->SetBranchAddress(weight_branch.c_str(), &weight);

  std::vector<double> pts;
  std::vector<double> base_w;
  std::vector<std::vector<double>> gvals;

  const Long64_t nentries = tree->GetEntries();
  pts.reserve(nentries * 2);
  base_w.reserve(nentries * 2);

  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);
    for (int j = 0; j < nJets; ++j) {
      double x = pt[j];
      if (x < pt_min || x > pt_max) continue;
      double l = safe_log(x, pt0);
      if (!std::isfinite(l)) continue;
      pts.push_back(x);
      base_w.push_back(static_cast<double>(weight));
    }
  }

  const size_t n = pts.size();
  if (n == 0) die("no jets passed the pT selection");

  gvals.assign(n, std::vector<double>(nm, 0.0));
  for (size_t i = 0; i < n; ++i) {
    double l = safe_log(pts[i], pt0);
    for (int j = 0; j < nm; ++j) {
      gvals[i][j] = std::pow(l, moments[j]);
    }
  }

  // Solve for lambda_j using Newton iterations.
  std::vector<double> lambda(nm, 0.0);
  auto compute_loss = [&](const std::vector<double> &mean_g) {
    double loss = 0.0;
    for (int j = 0; j < nm; ++j) {
      double num = c[j] - mean_g[j];
      double den = c[j] + mean_g[j];
      if (std::fabs(den) < 1e-12) continue;
      double r = num / den;
      loss += r * r;
    }
    return loss;
  };

  for (int iter = 0; iter < max_iter; ++iter) {
    double sumW = 0.0;
    std::vector<double> sumG(nm, 0.0);
    std::vector<std::vector<double>> sumGG(nm, std::vector<double>(nm, 0.0));

    for (size_t i = 0; i < n; ++i) {
      double dot = 0.0;
      for (int j = 0; j < nm; ++j) dot += lambda[j] * gvals[i][j];
      double w = std::exp(-dot) * base_w[i];
      sumW += w;
      for (int j = 0; j < nm; ++j) {
        sumG[j] += w * gvals[i][j];
      }
      for (int j = 0; j < nm; ++j) {
        for (int k = 0; k < nm; ++k) {
          sumGG[j][k] += w * gvals[i][j] * gvals[i][k];
        }
      }
    }
    if (sumW == 0.0) die("sum of weights is zero during optimization");

    std::vector<double> mean_g(nm, 0.0);
    for (int j = 0; j < nm; ++j) mean_g[j] = sumG[j] / sumW;

    std::vector<double> f(nm, 0.0);
    double max_rel = 0.0;
    for (int j = 0; j < nm; ++j) {
      f[j] = mean_g[j] - c[j];
      double rel = std::fabs(f[j]) / (std::fabs(c[j]) + 1e-12);
      if (rel > max_rel) max_rel = rel;
    }

    double loss = compute_loss(mean_g);
    std::cout << "iter " << iter
              << " loss=" << loss
              << " max_rel=" << max_rel
              << std::endl;
    if (max_rel < tol) break;

    std::vector<std::vector<double>> cov(nm, std::vector<double>(nm, 0.0));
    for (int j = 0; j < nm; ++j) {
      for (int k = 0; k < nm; ++k) {
        cov[j][k] = sumGG[j][k] / sumW - mean_g[j] * mean_g[k];
      }
      cov[j][j] += cov_reg;
    }

    std::vector<double> delta;
    if (!solve_linear(cov, f, delta)) {
      die("failed to solve linear system for lambda update");
    }

    double step = 1.0;
    while (step > 1e-3) {
      std::vector<double> lambda_try = lambda;
      for (int j = 0; j < nm; ++j) lambda_try[j] += step * delta[j];

      double sumW2 = 0.0;
      std::vector<double> sumG2(nm, 0.0);
      for (size_t i = 0; i < n; ++i) {
        double dot = 0.0;
        for (int j = 0; j < nm; ++j) dot += lambda_try[j] * gvals[i][j];
        double w = std::exp(-dot) * base_w[i];
        sumW2 += w;
        for (int j = 0; j < nm; ++j) sumG2[j] += w * gvals[i][j];
      }
      std::vector<double> mean_g2(nm, 0.0);
      for (int j = 0; j < nm; ++j) mean_g2[j] = sumG2[j] / sumW2;
      double loss2 = compute_loss(mean_g2);
      if (loss2 <= loss) {
        lambda = lambda_try;
        break;
      }
      step *= 0.5;
    }
  }

  // Write per-jet weights.
  TFile out(out_name.c_str(), "RECREATE");
  TTree tw("tweights", "per-jet unbiasing weights");
  float out_pt = 0.0f;
  double w_unbias = 1.0;
  double w_total = 1.0;
  double w_base = 1.0;

  tw.Branch("pt", &out_pt, "pt/F");
  tw.Branch("w_unbias", &w_unbias, "w_unbias/D");
  tw.Branch("w_base", &w_base, "w_base/D");
  tw.Branch("w_total", &w_total, "w_total/D");

  for (size_t i = 0; i < n; ++i) {
    double dot = 0.0;
    for (int j = 0; j < nm; ++j) dot += lambda[j] * gvals[i][j];
    w_unbias = std::exp(-dot);
    w_base = base_w[i];
    w_total = w_base * w_unbias;
    out_pt = static_cast<float>(pts[i]);
    tw.Fill();
  }

  // Save lambda values as a small tree for bookkeeping.
  TTree tmeta("meta", "fit metadata");
  std::vector<double> lambda_out = lambda;
  std::vector<int> moments_out = moments;
  tmeta.Branch("lambda", &lambda_out);
  tmeta.Branch("moments", &moments_out);
  tmeta.Branch("pt0", const_cast<double*>(&pt0), "pt0/D");
  tmeta.Branch("pt_min", const_cast<double*>(&pt_min), "pt_min/D");
  tmeta.Branch("pt_max", const_cast<double*>(&pt_max), "pt_max/D");
  tmeta.Fill();

  tw.Write();
  tmeta.Write();
  out.Close();
  in.Close();

  std::cout << "wrote " << out_name << std::endl;
  return 0;
}
