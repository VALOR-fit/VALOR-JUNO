#ifndef JUNO_DETECTOR_RESPONSE_H
#define JUNO_DETECTOR_RESPONSE_H

#include <array>
#include <map>
#include <string>
#include <utility>
#include <vector>

// VALOR FWMK
#include "valor/Conventions/Sample.h"
#include "valor/Conventions/ReacMode.h"
#include "valor/VArray/VALORArray.h"
#include "valor/Utilities/SingleSampleKineBinning.h"

class TH2D;

namespace juno {

class JunoDetectorResponse {
public:
  enum class Flavor {
    kNuE,
    kNuMu,
    kNuEBar,
    kNuMuBar,
    kNC
  };

  struct Config {
    bool apply_fc = false;
    bool apply_energy_smearing = false;
    bool apply_costh_smearing = false;
    bool apply_pid = false;
    bool preserve_normalization = false;
    double smear_nsigma = 5.0;
  };

  struct Component {
    int input_sample = 0;
    int reac_mode = 0;
    Flavor final_flavor = Flavor::kNC;
    const TH2D* true_spectrum = nullptr;
  };

  using HistKey = std::pair<int, int>; // (sample, reac_mode)
  using HistMap = std::map<HistKey, TH2D*>;

  JunoDetectorResponse();
  ~JunoDetectorResponse();

  void SetConfig(const Config& cfg);
  const Config& GetConfig() const;

  void SetAnalysisBinning(const std::vector<double>& e_edges,
                          const std::vector<double>& costh_edges);

  bool LoadFCCoefficientsText(const std::string& filename);
  void SetFCCoefficients(Flavor flavor, const std::vector<double>& coeffs);

  bool LoadPIDMatrixCSV(const std::string& filename);
  void SetSinglePIDMatrix(const std::array<std::array<double, 5>, 5>& matrix);

  HistMap Apply(const std::vector<Component>& components) const;
  bool ApplyInto(const std::vector<Component>& components, HistMap& output) const;

  static Flavor FlavorFromString(const std::string& s);
  static int TrueLabel(Flavor flavor);
  static int SampleFromRecoLabel(int reco_label);
  static const char* FlavorName(Flavor flavor);

  bool ApplyToVALORArrays(
    std::vector<valor::VALORArray*>& mct_fin_ob,
    const std::vector<bool>& mct_null_ob,
    const std::vector<valor::conv::ReacMode_t>& reaction_modes,
    const valor::SingleSampleKineBinning& binning,
    valor::conv::Sample_t sample,
    Option_t* opt = 0) const;

private:
  Config fConfig;
  std::vector<double> fEnergyEdges;
  std::vector<double> fCosthEdges;

  std::map<Flavor, std::vector<double> > fFCCoeffs;

  // 8 energy groups, matrix[true][reco]
  std::array<std::array<std::array<double, 5>, 5>, 8> fPIDMatrices;
  bool fHavePID = false;
  bool fUseSinglePID = true;

  double EvalFC(Flavor flavor, double energy) const;
  double EvalPoly(const std::vector<double>& coeffs, double x) const;

  double SigmaEnergy(double energy, Flavor flavor) const;
  double MuEnergyBias(double energy, Flavor flavor) const;
  double SigmaThetaDeg(double energy) const;
  double GaussianBinProb(double low, double high, double mean, double sigma) const;

  bool IsEFlavor(Flavor flavor) const;
  int FindBin(const std::vector<double>& edges, double x) const;
  int PIDGroupFromRecoEnergy(double ereco) const;

  std::vector<std::pair<int, double> > EnergyRecoProbs(double etrue,
                                                       Flavor flavor) const;
  std::vector<std::pair<int, double> > CosthRecoProbs(double costh_true,
                                                       double etrue) const;

  TH2D* EnsureHist(HistMap& output, int sample, int reac_mode) const;
};

} // namespace juno

#endif // JUNO_DETECTOR_RESPONSE_H
