#include "JunoDetectorResponse.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "TAxis.h"
#include "TH2D.h"
#include "TString.h"

namespace juno {

namespace {
constexpr double kSqrt2 = 1.4142135623730950488;
constexpr double kPi = 3.14159265358979323846;

std::string Trim(const std::string& s) {
  const std::string ws = " \t\r\n";
  const std::size_t b = s.find_first_not_of(ws);
  if (b == std::string::npos) return "";
  const std::size_t e = s.find_last_not_of(ws);
  return s.substr(b, e - b + 1);
}

std::vector<std::string> Split(const std::string& s, char delim) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) out.push_back(Trim(item));
  return out;
}
} // namespace

JunoDetectorResponse::JunoDetectorResponse() {
  fEnergyEdges = {
    1.08508435, 1.41421357, 1.84317469, 2.40224887, 3.13090217,
    4.08057156, 5.31829591, 6.93144844, 9.03390453, 11.77408038, 20.0
  };
  fCosthEdges = {
    -1.0, -0.845, -0.645, -0.45, -0.225,
     0.0,  0.225,  0.45,  0.645,  0.845, 1.0
  };

  for (auto& mat : fPIDMatrices) {
    for (int t = 0; t < 5; ++t) {
      for (int r = 0; r < 5; ++r) mat[t][r] = (t == r ? 1.0 : 0.0);
    }
  }
}

JunoDetectorResponse::~JunoDetectorResponse() = default;

void JunoDetectorResponse::SetConfig(const Config& cfg) { fConfig = cfg; }

const JunoDetectorResponse::Config& JunoDetectorResponse::GetConfig() const {
  return fConfig;
}

void JunoDetectorResponse::SetAnalysisBinning(
    const std::vector<double>& e_edges,
    const std::vector<double>& costh_edges) {
  if (e_edges.size() < 2 || costh_edges.size() < 2) {
    throw std::runtime_error(
        "JunoDetectorResponse::SetAnalysisBinning: need at least 2 edges per axis");
  }
  fEnergyEdges = e_edges;
  fCosthEdges = costh_edges;
}

bool JunoDetectorResponse::LoadFCCoefficientsText(const std::string& filename) {
  std::ifstream in(filename.c_str());
  if (!in.is_open()) {
    std::cerr << "JunoDetectorResponse: cannot open FC coefficient file: "
              << filename << std::endl;
    return false;
  }

  std::string line;
  while (std::getline(in, line)) {
    line = Trim(line);
    if (line.empty() || line[0] == '#') continue;

    const std::size_t eq = line.find('=');
    if (eq == std::string::npos) continue;

    const std::string key = Trim(line.substr(0, eq));
    const std::string vals = line.substr(eq + 1);
    std::vector<double> coeffs;
    for (const std::string& tok : Split(vals, ',')) {
      if (!tok.empty()) coeffs.push_back(std::atof(tok.c_str()));
    }

    try {
      SetFCCoefficients(FlavorFromString(key), coeffs);
    } catch (...) {
      std::cerr << "JunoDetectorResponse: unknown FC flavor key: "
                << key << std::endl;
      return false;
    }
  }

  return true;
}

bool JunoDetectorResponse::LoadFCCoefficientsText(
    const std::string& filename,
    Flavor flavor)
{
  std::ifstream in(filename.c_str());

  if (!in.is_open()) {
    std::cerr
      << "JunoDetectorResponse: cannot open FC coefficient file: "
      << filename
      << std::endl;
    return false;
  }

  std::vector<double> coeffs;
  std::string line;

  while (std::getline(in, line)) {
    line = Trim(line);

    if (line.empty()) {
      continue;
    }

    if (line[0] == '#') {
      continue;
    }

    // Accept comma-separated or whitespace-separated coefficients.
    std::replace(line.begin(), line.end(), ',', ' ');

    std::stringstream ss(line);
    double x = 0.0;

    while (ss >> x) {
      coeffs.push_back(x);
    }
  }

  if (coeffs.empty()) {
    std::cerr
      << "JunoDetectorResponse: no FC coefficients read from: "
      << filename
      << std::endl;
    return false;
  }

  SetFCCoefficients(flavor, coeffs);

  std::cerr
    << "[JUNO FC] loaded "
    << coeffs.size()
    << " coefficients for "
    << FlavorName(flavor)
    << " from "
    << filename
    << std::endl;

  return true;
}

void JunoDetectorResponse::SetFCCoefficients(
    Flavor flavor,
    const std::vector<double>& coeffs) {
  fFCCoeffs[flavor] = coeffs;
}

bool JunoDetectorResponse::LoadPIDMatrixCSV(const std::string& filename) {
  std::ifstream in(filename.c_str());
  if (!in.is_open()) {
    std::cerr << "JunoDetectorResponse: cannot open PID matrix CSV: "
              << filename << std::endl;
    return false;
  }

  for (auto& mat : fPIDMatrices) {
    for (int t = 0; t < 5; ++t) {
      for (int r = 0; r < 5; ++r) mat[t][r] = 0.0;
    }
  }

  std::string line;
  int n_loaded = 0;
  while (std::getline(in, line)) {
    line = Trim(line);
    if (line.empty() || line[0] == '#') continue;

    const std::vector<std::string> cols = Split(line, ',');
    if (cols.size() < 7) continue;
    if (cols[0] == "group") continue;

    const int group = std::atoi(cols[0].c_str());
    const int true_label = std::atoi(cols[1].c_str());
    if (group < 1 || group > 8 || true_label < 0 || true_label > 4) continue;

    for (int reco_label = 0; reco_label < 5; ++reco_label) {
      fPIDMatrices[group - 1][true_label][reco_label] =
          std::atof(cols[2 + reco_label].c_str());
    }
    ++n_loaded;
  }

  if (n_loaded <= 0) {
    std::cerr << "JunoDetectorResponse: PID matrix CSV contained no usable rows"
              << std::endl;
    return false;
  }

  fHavePID = true;
  fUseSinglePID = false;
  return true;
}

void JunoDetectorResponse::SetSinglePIDMatrix(
    const std::array<std::array<double, 5>, 5>& matrix) {
  for (auto& mat : fPIDMatrices) mat = matrix;
  fHavePID = true;
  fUseSinglePID = true;
}

JunoDetectorResponse::HistMap
JunoDetectorResponse::Apply(const std::vector<Component>& components) const {
  HistMap output;
  ApplyInto(components, output);
  return output;
}

bool JunoDetectorResponse::ApplyInto(
    const std::vector<Component>& components,
    HistMap& output) const {
  if (fEnergyEdges.size() < 2 || fCosthEdges.size() < 2) {
    std::cerr << "JunoDetectorResponse: analysis binning is invalid" << std::endl;
    return false;
  }

  if (fConfig.apply_pid && !fHavePID) {
    std::cerr << "JunoDetectorResponse: apply_pid requested but no PID matrix loaded"
              << std::endl;
    return false;
  }

  for (const Component& comp : components) {
    if (!comp.true_spectrum) continue;

    const TH2D* htrue = comp.true_spectrum;
    const int nETrue = htrue->GetNbinsX();
    const int nCTrue = htrue->GetNbinsY();

    for (int ie = 1; ie <= nETrue; ++ie) {
      const double etrue = htrue->GetXaxis()->GetBinCenter(ie);
      const std::vector<std::pair<int, double> > e_probs =
          EnergyRecoProbs(etrue, comp.final_flavor);
      if (e_probs.empty()) continue;

      for (int ic = 1; ic <= nCTrue; ++ic) {
        const double costh_true = htrue->GetYaxis()->GetBinCenter(ic);
        double rate = htrue->GetBinContent(ie, ic);
        if (rate <= 0.0) continue;

        if (fConfig.apply_fc) rate *= EvalFC(comp.final_flavor, etrue);
        if (rate <= 0.0) continue;

        const std::vector<std::pair<int, double> > c_probs =
            CosthRecoProbs(costh_true, etrue);
        if (c_probs.empty()) continue;

        for (const auto& ep : e_probs) {
          const int ie_reco = ep.first;
          const double pE = ep.second;
          if (pE <= 0.0) continue;
          const double ereco = 0.5 * (fEnergyEdges[ie_reco - 1] + fEnergyEdges[ie_reco]);

          for (const auto& cp : c_probs) {
            const int ic_reco = cp.first;
            const double pC = cp.second;
            if (pC <= 0.0) continue;
            const double creco = 0.5 * (fCosthEdges[ic_reco - 1] + fCosthEdges[ic_reco]);
            const double smeared = rate * pE * pC;
            if (smeared <= 0.0) continue;

            if (!fConfig.apply_pid) {
              TH2D* hout = EnsureHist(output, comp.input_sample, comp.reac_mode);
              hout->Fill(ereco, creco, smeared);
              continue;
            }

            const int true_label = TrueLabel(comp.final_flavor);
            const int group = fUseSinglePID ? 0 : PIDGroupFromRecoEnergy(ereco);

            for (int reco_label = 0; reco_label < 5; ++reco_label) {
              const double pid_prob = fPIDMatrices[group][true_label][reco_label];
              if (pid_prob <= 0.0) continue;
              const int reco_sample = SampleFromRecoLabel(reco_label);
              if (reco_sample < 0) continue;

              TH2D* hout = EnsureHist(output, reco_sample, comp.reac_mode);
              hout->Fill(ereco, creco, smeared * pid_prob);
            }
          }
        }
      }
    }
  }

  return true;
}

JunoDetectorResponse::Flavor
JunoDetectorResponse::FlavorFromString(const std::string& s) {
  const std::string x = Trim(s);
  if (x == "nue" || x == "NuE") return Flavor::kNuE;
  if (x == "numu" || x == "NuMu") return Flavor::kNuMu;
  if (x == "nuebar" || x == "antinue" || x == "NuEBar") return Flavor::kNuEBar;
  if (x == "numubar" || x == "antinumu" || x == "NuMuBar") return Flavor::kNuMuBar;
  if (x == "nc" || x == "NC") return Flavor::kNC;
  throw std::runtime_error("unknown flavor string: " + s);
}

int JunoDetectorResponse::TrueLabel(Flavor flavor) {
  switch (flavor) {
    case Flavor::kNuE:     return 0;
    case Flavor::kNuMu:    return 1;
    case Flavor::kNuEBar:  return 2;
    case Flavor::kNuMuBar: return 3;
    case Flavor::kNC:      return 4;
  }
  return 4;
}

int JunoDetectorResponse::SampleFromRecoLabel(int reco_label) {
  if (reco_label == 0) return 52; // nueCC_fc
  if (reco_label == 1) return 50; // numuCC_fc
  if (reco_label == 2) return 53; // nuebarCC_fc
  if (reco_label == 3) return 51; // numubarCC_fc
  return -1; // NC ignored in current 4CC setup
}

const char* JunoDetectorResponse::FlavorName(Flavor flavor) {
  switch (flavor) {
    case Flavor::kNuE:     return "nue";
    case Flavor::kNuMu:    return "numu";
    case Flavor::kNuEBar:  return "nuebar";
    case Flavor::kNuMuBar: return "numubar";
    case Flavor::kNC:      return "nc";
  }
  return "unknown";
}

double JunoDetectorResponse::EvalFC(Flavor flavor, double energy) const {
  const auto it = fFCCoeffs.find(flavor);
  if (it == fFCCoeffs.end() || it->second.empty()) return 1.0;

  double v = EvalPoly(it->second, energy);
  if (v < 0.0) v = 0.0;
  if (v > 1.0) v = 1.0;
  return v;
}

double JunoDetectorResponse::EvalPoly(
    const std::vector<double>& coeffs,
    double x) const {
  if (coeffs.empty()) return 1.0;
  double y = 0.0;
  for (double a : coeffs) y = y * x + a;
  return y;
}

bool JunoDetectorResponse::IsEFlavor(Flavor flavor) const {
  return flavor == Flavor::kNuE || flavor == Flavor::kNuEBar;
}

double JunoDetectorResponse::SigmaEnergy(double energy, Flavor flavor) const {
  const double e = std::max(1.0, energy);
  if (IsEFlavor(flavor)) {
    return 0.0444 + 0.003280 * std::sqrt(e) + 0.07683 / e - 0.00218 * e;
  }
  return 0.0535 + 0.005934 * std::sqrt(e) + 0.06803 / e - 0.00267 * e;
}

double JunoDetectorResponse::MuEnergyBias(double energy, Flavor flavor) const {
  const double e = std::max(1.0, energy);
  if (IsEFlavor(flavor)) {
    return -0.4557 + 0.235222 * std::sqrt(e) + 0.29617 / e - 0.03182 * e;
  }
  return -0.2657 + 0.160940 * std::sqrt(e) + 0.08111 / e - 0.02438 * e;
}

double JunoDetectorResponse::SigmaThetaDeg(double energy) const {
  const double e = std::max(1.0, energy);
  return 0.593 + 1.040 * std::sqrt(e) + 26.175 / e - 0.140 * e;
}

double JunoDetectorResponse::GaussianBinProb(
    double low,
    double high,
    double mean,
    double sigma) const {
  if (sigma <= 0.0) return 0.0;
  const double z1 = (low - mean) / (kSqrt2 * sigma);
  const double z2 = (high - mean) / (kSqrt2 * sigma);
  return 0.5 * (std::erf(z2) - std::erf(z1));
}

int JunoDetectorResponse::FindBin(const std::vector<double>& edges, double x) const {
  const int nbins = static_cast<int>(edges.size()) - 1;
  for (int i = 1; i <= nbins; ++i) {
    if (x >= edges[i - 1] && x < edges[i]) return i;
  }
  if (nbins >= 1 && x == edges[nbins]) return nbins;
  return -1;
}

int JunoDetectorResponse::PIDGroupFromRecoEnergy(double ereco) const {
  const int nE = static_cast<int>(fEnergyEdges.size()) - 1;
  int ebin = nE - 1;
  for (int i = 0; i < nE; ++i) {
    if (ereco < fEnergyEdges[i + 1]) {
      ebin = i;
      break;
    }
  }

  if (ebin <= 2) return 0;
  if (ebin == 3) return 1;
  if (ebin == 4) return 2;
  if (ebin == 5) return 3;
  if (ebin == 6) return 4;
  if (ebin == 7) return 5;
  if (ebin == 8) return 6;
  return 7;
}

std::vector<std::pair<int, double> >
JunoDetectorResponse::EnergyRecoProbs(double etrue, Flavor flavor) const {
  std::vector<std::pair<int, double> > probs;

  if (!fConfig.apply_energy_smearing) {
    const int ibin = FindBin(fEnergyEdges, etrue);
    if (ibin > 0) probs.emplace_back(ibin, 1.0);
    return probs;
  }

  const double sigma = std::max(1.0e-6, SigmaEnergy(etrue, flavor));
  const double mean = etrue + MuEnergyBias(etrue, flavor);
  double norm = 0.0;

  const int nbins = static_cast<int>(fEnergyEdges.size()) - 1;
  for (int ibin = 1; ibin <= nbins; ++ibin) {
    const double low = fEnergyEdges[ibin - 1];
    const double high = fEnergyEdges[ibin];
    if (high < mean - fConfig.smear_nsigma * sigma) continue;
    if (low  > mean + fConfig.smear_nsigma * sigma) continue;

    const double p = GaussianBinProb(low, high, mean, sigma);
    if (p > 0.0) {
      probs.emplace_back(ibin, p);
      norm += p;
    }
  }

  if (fConfig.preserve_normalization && norm > 0.0) {
    for (auto& x : probs) x.second /= norm;
  }

  return probs;
}

// std::vector<std::pair<int, double> >
// JunoDetectorResponse::CosthRecoProbs(double costh_true, double etrue) const {
//   std::vector<std::pair<int, double> > probs;

//   if (!fConfig.apply_costh_smearing) {
//     const int ibin = FindBin(fCosthEdges, costh_true);
//     if (ibin > 0) probs.emplace_back(ibin, 1.0);
//     return probs;
//   }

//   const double c = std::max(-1.0, std::min(1.0, costh_true));
//   const double theta = std::acos(c);
//   const double sigma_theta_rad =
//       std::max(1.0e-6, SigmaThetaDeg(etrue)) * kPi / 180.0;
//   const double sigma_c =
//       std::max(1.0e-6, std::abs(std::sin(theta)) * sigma_theta_rad);
//   const double mean = c;
//   double norm = 0.0;

//   const int nbins = static_cast<int>(fCosthEdges.size()) - 1;
//   for (int ibin = 1; ibin <= nbins; ++ibin) {
//     const double low = fCosthEdges[ibin - 1];
//     const double high = fCosthEdges[ibin];
//     if (high < mean - fConfig.smear_nsigma * sigma_c) continue;
//     if (low  > mean + fConfig.smear_nsigma * sigma_c) continue;

//     const double p = GaussianBinProb(low, high, mean, sigma_c);
//     if (p > 0.0) {
//       probs.emplace_back(ibin, p);
//       norm += p;
//     }
//   }

//   if (fConfig.preserve_normalization && norm > 0.0) {
//     for (auto& x : probs) x.second /= norm;
//   }

//   return probs;
// }

std::vector<std::pair<int, double> >
JunoDetectorResponse::CosthRecoProbs(double costh_true, double etrue) const {
  std::vector<std::pair<int, double> > probs;

  if (!fConfig.apply_costh_smearing) {
    const int ibin = FindBin(fCosthEdges, costh_true);
    if (ibin > 0) probs.emplace_back(ibin, 1.0);
    return probs;
  }

  const double c =
    std::max(-1.0, std::min(1.0, costh_true));

  const double theta =
    std::acos(c);

  const double sigma_theta_rad =
    std::max(1.0e-6, SigmaThetaDeg(etrue)) * kPi / 180.0;

  const double sigma_c =
    std::max(1.0e-6, std::abs(std::sin(theta)) * sigma_theta_rad);

  const double mean = c;

  const int nbins =
    static_cast<int>(fCosthEdges.size()) - 1;

  if (nbins <= 0) {
    return probs;
  }

  const double support_low =
    mean - fConfig.smear_nsigma * sigma_c;

  const double support_high =
    mean + fConfig.smear_nsigma * sigma_c;

  auto add_prob = [&](int ibin, double p) {
    if (p <= 0.0) return;

    for (auto& x : probs) {
      if (x.first == ibin) {
        x.second += p;
        return;
      }
    }

    probs.emplace_back(ibin, p);
  };

  auto overlap_prob = [&](double a, double b) {
    const double lo = std::max(a, support_low);
    const double hi = std::min(b, support_high);

    if (hi <= lo) {
      return 0.0;
    }

    return GaussianBinProb(lo, hi, mean, sigma_c);
  };

  /*
    Reflection folding on [-1, 1]:

    For each integer k:

      identity branch:
        x in [-1 + 4k, 1 + 4k]
        c_folded = x - 4k

        reco bin [low, high] receives x in:
        [low + 4k, high + 4k]

      reflected branch:
        x in [1 + 4k, 3 + 4k]
        c_folded = 2 + 4k - x

        reco bin [low, high] receives x in:
        [2 + 4k - high, 2 + 4k - low]

    This automatically handles values beyond -2, +2, etc.
  */

  const int k_min =
    static_cast<int>(std::floor((support_low - 3.0) / 4.0)) - 1;

  const int k_max =
    static_cast<int>(std::ceil((support_high + 3.0) / 4.0)) + 1;

  for (int ibin = 1; ibin <= nbins; ++ibin) {
    const double low =
      fCosthEdges[ibin - 1];

    const double high =
      fCosthEdges[ibin];

    double p_total = 0.0;

    for (int k = k_min; k <= k_max; ++k) {
      const double shift = 4.0 * k;

      // Identity branch: x -> x - 4k
      const double id_low =
        low + shift;

      const double id_high =
        high + shift;

      p_total += overlap_prob(id_low, id_high);

      // Reflected branch: x -> 2 + 4k - x
      const double ref_low =
        2.0 + shift - high;

      const double ref_high =
        2.0 + shift - low;

      p_total += overlap_prob(ref_low, ref_high);
    }

    add_prob(ibin, p_total);
  }

  double norm = 0.0;
  for (const auto& x : probs) {
    norm += x.second;
  }

  /*
    With mirror reflection, total probability inside [-1, 1]
    should already be close to the Gaussian probability inside
    mean ± nsigma*sigma. If preserve_normalization is true,
    renormalize the truncated support to exactly 1.
  */
  if (fConfig.preserve_normalization && norm > 0.0) {
    for (auto& x : probs) {
      x.second /= norm;
    }
  }

  return probs;
}

TH2D* JunoDetectorResponse::EnsureHist(
    HistMap& output,
    int sample,
    int reac_mode) const {
  const HistKey key(sample, reac_mode);
  auto it = output.find(key);
  if (it != output.end()) return it->second;

  TString name = TString::Format("h2_sample%d_mode%d", sample, reac_mode);
  TString title = TString::Format(
      "sample %d mode %d;E_{reco};cos#theta_{reco}", sample, reac_mode);

  TH2D* h = new TH2D(name, title,
                     static_cast<int>(fEnergyEdges.size()) - 1, &fEnergyEdges[0],
                     static_cast<int>(fCosthEdges.size()) - 1, &fCosthEdges[0]);
  h->SetDirectory(nullptr);
  output[key] = h;
  return h;
}

bool JunoDetectorResponse::ApplyToVALORArrays(
    std::vector<valor::VALORArray*> & mct_fin_ob,
    const std::vector<bool> & mct_null_ob,
    const std::vector<valor::conv::ReacMode_t> & modes,
    const valor::SingleSampleKineBinning & binning,
    valor::conv::Sample_t sample,
    Option_t * opt) const
{
  (void)opt;

  if (fConfig.apply_pid) {
    std::cerr
      << "JunoDetectorResponse::ApplyToVALORArrays: PID is requested, "
      << "but PID requires cross-sample migration and cannot be applied "
      << "inside one JunoSingleSampleFitModel. Disable apply_pid here "
      << "or handle PID in CompoundFitModel."
      << std::endl;
    return false;
  }

  Flavor flavor = Flavor::kNC;

  if      (sample == valor::conv::ksJUNO_FC_CC_Nue)    flavor = Flavor::kNuE;
  else if (sample == valor::conv::ksJUNO_FC_CC_Numu)   flavor = Flavor::kNuMu;
  else if (sample == valor::conv::ksJUNO_FC_CC_Nuebar) flavor = Flavor::kNuEBar;
  else if (sample == valor::conv::ksJUNO_FC_CC_Numubar)flavor = Flavor::kNuMuBar;
  else {
    std::cerr
      << "JunoDetectorResponse::ApplyToVALORArrays: unsupported sample = "
      << sample << std::endl;
    return false;
  }

  const valor::VALORBinning * reco_ob = binning.RecoOptBinning();

  if (!reco_ob) {
    std::cerr
      << "JunoDetectorResponse::ApplyToVALORArrays: null reco binning"
      << std::endl;
    return false;
  }

  const int n_reco_bins = reco_ob->NBins();

  for (size_t imode = 0; imode < modes.size(); ++imode) {

    if (modes[imode] == valor::conv::kmUndefined) {
      continue;
    }

    if (imode < mct_null_ob.size() && mct_null_ob[imode]) {
      continue;
    }

    if (imode >= mct_fin_ob.size()) {
      continue;
    }

    valor::VALORArray * arr = mct_fin_ob[imode];

    if (!arr) {
      continue;
    }

    const valor::VALORBinning * true_sb =
      binning.TrueSmpBinning(modes[imode]);

    if (!true_sb) {
      continue;
    }

    const int n_true_bins = true_sb->NBins();

    valor::VALORArray out(n_reco_bins, n_true_bins);
    out.SetAll(0.0);

    for (int itrue = 0; itrue < n_true_bins; ++itrue) {

      const double etrue =
        true_sb->BinCentre(itrue, valor::conv::kvEtrue);

      const double costh_true =
        true_sb->BinCentre(itrue, valor::conv::kvQ2true);

      double true_rate = 0.0;

      for (int ireco = 0; ireco < n_reco_bins; ++ireco) {
        true_rate += arr->GetEntry(ireco, itrue);
      }

      if (true_rate <= 0.0) {
        continue;
      }

      if (fConfig.apply_fc) {
        true_rate *= EvalFC(flavor, etrue);
      }

      if (true_rate <= 0.0) {
        continue;
      }

      const std::vector<std::pair<int, double> > e_probs =
        EnergyRecoProbs(etrue, flavor);

      const std::vector<std::pair<int, double> > c_probs =
        CosthRecoProbs(costh_true, etrue);

      if (e_probs.empty() || c_probs.empty()) {
        continue;
      }

      for (const auto & ep : e_probs) {
        const int ie_reco = ep.first;
        const double pE = ep.second;

        if (ie_reco <= 0 || pE <= 0.0) {
          continue;
        }

        for (const auto & cp : c_probs) {
          const int ic_reco = cp.first;
          const double pC = cp.second;

          if (ic_reco <= 0 || pC <= 0.0) {
            continue;
          }

          const double ereco =
            0.5 * (fEnergyEdges[ie_reco - 1] + fEnergyEdges[ie_reco]);

          const double creco =
            0.5 * (fCosthEdges[ic_reco - 1] + fCosthEdges[ic_reco]);

          std::vector<double> reco_vars;
          reco_vars.push_back(ereco);
          reco_vars.push_back(creco);

          const int ireco_ob = reco_ob->FindBin(reco_vars);

          if (ireco_ob < 0) {
            continue;
          }

          const double old = out.GetEntry(ireco_ob, itrue);

          out.SetEntry(
              old + true_rate * pE * pC,
              ireco_ob,
              itrue);
        }
      }
    }

    for (int itrue = 0; itrue < n_true_bins; ++itrue) {
      for (int ireco = 0; ireco < n_reco_bins; ++ireco) {
        arr->SetEntry(out.GetEntry(ireco, itrue), ireco, itrue);
      }
    }
  }

  return true;
}

} // namespace juno
