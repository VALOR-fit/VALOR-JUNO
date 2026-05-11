#!/usr/bin/env bash
set -euo pipefail

OUTPUT_FILE="data/valor_juno_grid_dst_IO_compact_smear.root"
FLUX_NUMU=""
FLUX_NUE=""
FLUX_NUMUBAR=""
FLUX_NUEBAR=""
XSEC_C12="data/input/tn_13581/xcheck_xsec_tot_cc_C12.root"
XSEC_H1="data/input/tn_13581/xcheck_xsec_tot_cc_H1.root"

APPLY_OSC=0
OSC_NU_FILE=""
OSC_NB_FILE=""
INCLUDE_CROSS_COMPONENTS=0
APPLY_FC_RATIO=0
FC_COEFF_DIR="data/input/fc_ratio/Total_FC_V2"

APPLY_SMEARING=0
APPLY_ENERGY_SMEARING=0
APPLY_COSTH_SMEARING=0
SMEAR_NSIGMA=5.0

APPLY_PID=0
PID_MATRIX_FILE="data/input/detector_response/pid_matrices_8bin_from_user.csv"
PRESERVE_NORMALIZATION=0
EVENT_STYLE_DST=0

RUNTIME_YEARS=6.0
MASS_KG=20000000.0
SOLID_ANGLE=0.0314159265359
CARBON_RATIO=0.88
HYDROGEN_RATIO=0.12

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/build_smear_from_maps.sh [options]

Required flux inputs:
  --flux-numu <file.root>
  --flux-nue <file.root>
  --flux-numubar <file.root>
  --flux-nuebar <file.root>

Options:
  --output <file.root>
  --xsec-c12 <file.root>
  --xsec-h1 <file.root>

  --apply-osc
  --osc-nu-file <file.root>
  --osc-nb-file <file.root>
  --include-cross-components

  --apply-fc-ratio
  --fc-coeff-dir <dir>

  --apply-smearing              Apply both energy and cos(theta) smearing
  --apply-energy-smearing       Apply only energy smearing
  --apply-costh-smearing        Apply only cos(theta) smearing
  --smear-nsigma <float>        Truncate Gaussian support to mean ± N sigma

  --apply-pid
  --pid-matrix-file <file.csv>

  --preserve-normalization
  --compact-template-dst        Write compact/template-style DST (enu_true=enu_reco)
  --event-style-dst             Write event-style DST (truth preserved; much larger)

  --runtime-years <float>
  --mass-kg <float>
  --solid-angle <float>
  --carbon-ratio <float>
  --hydrogen-ratio <float>
  --help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --output) OUTPUT_FILE="$2"; shift 2 ;;
    --flux-numu) FLUX_NUMU="$2"; shift 2 ;;
    --flux-nue) FLUX_NUE="$2"; shift 2 ;;
    --flux-numubar) FLUX_NUMUBAR="$2"; shift 2 ;;
    --flux-nuebar) FLUX_NUEBAR="$2"; shift 2 ;;
    --xsec-c12) XSEC_C12="$2"; shift 2 ;;
    --xsec-h1) XSEC_H1="$2"; shift 2 ;;
    --apply-osc) APPLY_OSC=1; shift ;;
    --osc-nu-file) OSC_NU_FILE="$2"; shift 2 ;;
    --osc-nb-file) OSC_NB_FILE="$2"; shift 2 ;;
    --include-cross-components) INCLUDE_CROSS_COMPONENTS=1; shift ;;
    --apply-fc-ratio) APPLY_FC_RATIO=1; shift ;;
    --fc-coeff-dir) FC_COEFF_DIR="$2"; shift 2 ;;
    --apply-smearing) APPLY_SMEARING=1; APPLY_ENERGY_SMEARING=1; APPLY_COSTH_SMEARING=1; shift ;;
    --apply-energy-smearing) APPLY_ENERGY_SMEARING=1; shift ;;
    --apply-costh-smearing) APPLY_COSTH_SMEARING=1; shift ;;
    --smear-nsigma) SMEAR_NSIGMA="$2"; shift 2 ;;
    --apply-pid) APPLY_PID=1; shift ;;
    --pid-matrix-file) PID_MATRIX_FILE="$2"; shift 2 ;;
    --preserve-normalization) PRESERVE_NORMALIZATION=1; shift ;;
    --compact-template-dst) EVENT_STYLE_DST=0; shift ;;
    --event-style-dst) EVENT_STYLE_DST=1; shift ;;
    --runtime-years) RUNTIME_YEARS="$2"; shift 2 ;;
    --mass-kg) MASS_KG="$2"; shift 2 ;;
    --solid-angle) SOLID_ANGLE="$2"; shift 2 ;;
    --carbon-ratio) CARBON_RATIO="$2"; shift 2 ;;
    --hydrogen-ratio) HYDROGEN_RATIO="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *)
      echo "[ERROR] Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

if [[ "$APPLY_ENERGY_SMEARING" == "1" || "$APPLY_COSTH_SMEARING" == "1" ]]; then
  APPLY_SMEARING=1
fi

if [[ "$APPLY_FC_RATIO" == "1" ]]; then
  for c in coeffs_numu.npy coeffs_nue.npy coeffs_antinumu.npy coeffs_antinue.npy; do
    if [[ ! -f "$FC_COEFF_DIR/$c" ]]; then
      echo "[ERROR] Missing FC coefficient file: $FC_COEFF_DIR/$c"
      exit 1
    fi
  done
fi

for req in "$FLUX_NUMU" "$FLUX_NUE" "$FLUX_NUMUBAR" "$FLUX_NUEBAR" "$XSEC_C12" "$XSEC_H1"; do
  if [[ -z "$req" || ! -f "$req" ]]; then
    echo "[ERROR] Missing required file: $req"
    exit 1
  fi
done

if [[ "$APPLY_OSC" == "1" ]]; then
  if [[ -z "$OSC_NU_FILE" || -z "$OSC_NB_FILE" ]]; then
    echo "[ERROR] --apply-osc requires --osc-nu-file and --osc-nb-file"
    exit 1
  fi
  if [[ ! -f "$OSC_NU_FILE" || ! -f "$OSC_NB_FILE" ]]; then
    echo "[ERROR] Missing oscillation probability file(s)."
    exit 1
  fi
fi

if [[ "$APPLY_PID" == "1" && ! -f "$PID_MATRIX_FILE" ]]; then
  echo "[ERROR] --apply-pid requested but PID matrix file not found: $PID_MATRIX_FILE"
  exit 1
fi

FC_COEFF_NUMU=""
FC_COEFF_NUE=""
FC_COEFF_NUMUBAR=""
FC_COEFF_NUEBAR=""

if [[ "$APPLY_FC_RATIO" == "1" ]]; then
  if ! command -v python3 >/dev/null 2>&1; then
    echo "[ERROR] --apply-fc-ratio requires python3 with numpy available."
    exit 1
  fi

  FC_TMP="/tmp/fc_coeff_dump_$$.txt"
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 - <<PY > "$FC_TMP"
import numpy as np
base = r"$FC_COEFF_DIR"
mapping = {
    "numu": "coeffs_numu.npy",
    "nue": "coeffs_nue.npy",
    "numubar": "coeffs_antinumu.npy",
    "nuebar": "coeffs_antinue.npy",
}
for key, fn in mapping.items():
    coeff = np.load(f"{base}/{fn}", allow_pickle=True).astype(float).ravel()
    print(f"{key}=" + ",".join(f"{x:.18e}" for x in coeff))
PY

  FC_COEFF_NUMU=$(grep '^numu=' "$FC_TMP" | cut -d= -f2-)
  FC_COEFF_NUE=$(grep '^nue=' "$FC_TMP" | cut -d= -f2-)
  FC_COEFF_NUMUBAR=$(grep '^numubar=' "$FC_TMP" | cut -d= -f2-)
  FC_COEFF_NUEBAR=$(grep '^nuebar=' "$FC_TMP" | cut -d= -f2-)
  rm -f "$FC_TMP"
fi

TMP_MACRO="/tmp/build_grid_dst_compact_from_maps_$$.C"
cat > "$TMP_MACRO" <<ROOTMACRO
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <map>
#include <utility>
#include <cstdlib>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TString.h"

void build_grid_dst_compact_from_maps() {
  TH1::AddDirectory(kFALSE);

  const std::string outFile = "${OUTPUT_FILE}";
  const std::string fluxNumuFile = "${FLUX_NUMU}";
  const std::string fluxNueFile = "${FLUX_NUE}";
  const std::string fluxNumubarFile = "${FLUX_NUMUBAR}";
  const std::string fluxNuebarFile = "${FLUX_NUEBAR}";
  const std::string xsecC12File = "${XSEC_C12}";
  const std::string xsecH1File = "${XSEC_H1}";

  const bool applyOsc = (${APPLY_OSC} != 0);
  const bool includeCrossComponents = (${INCLUDE_CROSS_COMPONENTS} != 0);
  const bool applyFCRatio = (${APPLY_FC_RATIO} != 0);

  const bool applySmearing = (${APPLY_SMEARING} != 0);
  const bool applyEnergySmearing = (${APPLY_ENERGY_SMEARING} != 0);
  const bool applyCosthSmearing  = (${APPLY_COSTH_SMEARING} != 0);
  const double smearNSigma = ${SMEAR_NSIGMA};

  const bool applyPID = (${APPLY_PID} != 0);
  const bool preserveNormalization = (${PRESERVE_NORMALIZATION} != 0);
  const bool eventStyleDST = (${EVENT_STYLE_DST} != 0);

  const std::string pidMatrixFile = "${PID_MATRIX_FILE}";
  const std::string oscNuFile = "${OSC_NU_FILE}";
  const std::string oscNbFile = "${OSC_NB_FILE}";

  const std::vector<double> fcCoeffNumu = {${FC_COEFF_NUMU}};
  const std::vector<double> fcCoeffNue = {${FC_COEFF_NUE}};
  const std::vector<double> fcCoeffNumubar = {${FC_COEFF_NUMUBAR}};
  const std::vector<double> fcCoeffNuebar = {${FC_COEFF_NUEBAR}};

  auto evalPoly = [](const std::vector<double>& c, double x) {
    if (c.empty()) return 1.0;
    double y = 0.0;
    for (double a : c) y = y * x + a;
    return y;
  };

  auto getFCRatio = [&](const std::string& finalFlavor, double energy) {
    if (!applyFCRatio) return 1.0;
    double v = 1.0;
    if (finalFlavor == "numu") v = evalPoly(fcCoeffNumu, energy);
    else if (finalFlavor == "nue") v = evalPoly(fcCoeffNue, energy);
    else if (finalFlavor == "numubar") v = evalPoly(fcCoeffNumubar, energy);
    else if (finalFlavor == "nuebar") v = evalPoly(fcCoeffNuebar, energy);
    if (v < 0.0) v = 0.0;
    return v;
  };

  auto isEFlavor = [](const std::string& finalFlavor) {
    return (finalFlavor == "nue" || finalFlavor == "nuebar");
  };

  auto sigmaEnergy = [&](double E, const std::string& finalFlavor) {
    const double e = std::max(1.0, E);
    if (isEFlavor(finalFlavor)) {
      return 0.0444 + 0.003280 * std::sqrt(e) + 0.07683 / e - 0.00218 * e;
    }
    return 0.0535 + 0.005934 * std::sqrt(e) + 0.06803 / e - 0.00267 * e;
  };

  auto muEnergyBias = [&](double E, const std::string& finalFlavor) {
    const double e = std::max(1.0, E);
    if (isEFlavor(finalFlavor)) {
      return -0.4557 + 0.235222 * std::sqrt(e) + 0.29617 / e - 0.03182 * e;
    }
    return -0.2657 + 0.160940 * std::sqrt(e) + 0.08111 / e - 0.02438 * e;
  };

  auto sigmaThetaDeg = [&](double E) {
    const double e = std::max(1.0, E);
    return 0.593 + 1.040 * std::sqrt(e) + 26.175 / e - 0.140 * e;
  };

  const double runtimeYears = ${RUNTIME_YEARS};
  const double massKg = ${MASS_KG};
  const double solidAngle = ${SOLID_ANGLE};
  const double carbonRatio = ${CARBON_RATIO};
  const double hydrogenRatio = ${HYDROGEN_RATIO};

  const double runtime = runtimeYears * 365.25 * 24.0 * 3600.0;
  const double unitConversion = 1e-42;
  const double avogadro = 6.02214076e23;
  const double massCarbonG = massKg * carbonRatio * 1e3;
  const double massHydrogenG = massKg * hydrogenRatio * 1e3;
  const double atomsC12 = (massCarbonG / 12.0) * avogadro;
  const double atomsH1  = (massHydrogenG / 1.0) * avogadro;

  const double neutrinoAtmoHeight = 15.0;
  const double earthRadius = 6341.0;
  const double sqrt2 = std::sqrt(2.0);

  const int nE_ana = 10;
  const int nC_ana = 10;
  const double e_edges_ana[nE_ana + 1] = {
    1.08508435, 1.41421357, 1.84317469, 2.40224887, 3.13090217,
    4.08057156, 5.31829591, 6.93144844, 9.03390453, 11.77408038, 20.0
  };
  const double c_edges_ana[nC_ana + 1] = {
    -1.0, -0.845, -0.645, -0.45, -0.225,
     0.0,  0.225,  0.45,  0.645,  0.845, 1.0
  };

  auto gaussianBinProb = [&](double low, double high, double mean, double sigma) {
    if (sigma <= 0.0) return 0.0;
    const double z1 = (low - mean) / (sqrt2 * sigma);
    const double z2 = (high - mean) / (sqrt2 * sigma);
    return 0.5 * (std::erf(z2) - std::erf(z1));
  };

  TFile fNumu(fluxNumuFile.c_str(), "READ");
  TFile fNue(fluxNueFile.c_str(), "READ");
  TFile fNumubar(fluxNumubarFile.c_str(), "READ");
  TFile fNuebar(fluxNuebarFile.c_str(), "READ");
  TFile fXC(xsecC12File.c_str(), "READ");
  TFile fXH(xsecH1File.c_str(), "READ");

  TH2D* fluxNumu = dynamic_cast<TH2D*>(fNumu.Get("custom_flux_hist"));
  TH2D* fluxNue = dynamic_cast<TH2D*>(fNue.Get("custom_flux_hist"));
  TH2D* fluxNumubar = dynamic_cast<TH2D*>(fNumubar.Get("custom_flux_hist"));
  TH2D* fluxNuebar = dynamic_cast<TH2D*>(fNuebar.Get("custom_flux_hist"));

  if (!fluxNumu || !fluxNue || !fluxNumubar || !fluxNuebar) {
    std::cerr << "[ERROR] Missing custom_flux_hist in one or more flux files." << std::endl;
    return;
  }

  TH1D* xNuMuC = dynamic_cast<TH1D*>(fXC.Get("nu_mu_tot_cc"));
  TH1D* xNuEC = dynamic_cast<TH1D*>(fXC.Get("nu_e_tot_cc"));
  TH1D* xNuMuBarC = dynamic_cast<TH1D*>(fXC.Get("nu_mu_bar_tot_cc"));
  TH1D* xNuEBarC = dynamic_cast<TH1D*>(fXC.Get("nu_e_bar_tot_cc"));

  TH1D* xNuMuH = dynamic_cast<TH1D*>(fXH.Get("nu_mu_tot_cc"));
  TH1D* xNuEH = dynamic_cast<TH1D*>(fXH.Get("nu_e_tot_cc"));
  TH1D* xNuMuBarH = dynamic_cast<TH1D*>(fXH.Get("nu_mu_bar_tot_cc"));
  TH1D* xNuEBarH = dynamic_cast<TH1D*>(fXH.Get("nu_e_bar_tot_cc"));

  if (!xNuMuC || !xNuEC || !xNuMuBarC || !xNuEBarC ||
      !xNuMuH || !xNuEH || !xNuMuBarH || !xNuEBarH) {
    std::cerr << "[ERROR] Missing one or more CC xsec histograms." << std::endl;
    return;
  }

  TFile* pOscNu = nullptr;
  TFile* pOscNb = nullptr;
  TH2D *PmmNu = nullptr, *PemNu = nullptr, *PeeNu = nullptr, *PmeNu = nullptr;
  TH2D *PmmNb = nullptr, *PemNb = nullptr, *PeeNb = nullptr, *PmeNb = nullptr;

  if (applyOsc) {
    pOscNu = new TFile(oscNuFile.c_str(), "READ");
    pOscNb = new TFile(oscNbFile.c_str(), "READ");
    if (!pOscNu || pOscNu->IsZombie() || !pOscNb || pOscNb->IsZombie()) {
      std::cerr << "[ERROR] Failed to open oscillation probability file(s)." << std::endl;
      return;
    }

    PmmNu = dynamic_cast<TH2D*>(pOscNu->Get("NuMuToNuMu3f"));
    PemNu = dynamic_cast<TH2D*>(pOscNu->Get("NuEToNuMu3f"));
    PeeNu = dynamic_cast<TH2D*>(pOscNu->Get("NuEToNuE3f"));
    PmeNu = dynamic_cast<TH2D*>(pOscNu->Get("NuMuToNuE3f"));

    PmmNb = dynamic_cast<TH2D*>(pOscNb->Get("NuMuToNuMu3f"));
    PemNb = dynamic_cast<TH2D*>(pOscNb->Get("NuEToNuMu3f"));
    PeeNb = dynamic_cast<TH2D*>(pOscNb->Get("NuEToNuE3f"));
    PmeNb = dynamic_cast<TH2D*>(pOscNb->Get("NuMuToNuE3f"));

    if (!PmmNu || !PemNu || !PeeNu || !PmeNu || !PmmNb || !PemNb || !PeeNb || !PmeNb) {
      std::cerr << "[ERROR] Missing required oscillation TH2D in osc files." << std::endl;
      return;
    }
  }

  std::vector<std::array<double, 25>> pidMats(8);
  for (auto& m : pidMats) m.fill(0.0);

  if (applyPID) {
    std::ifstream pin(pidMatrixFile);
    if (!pin.is_open()) {
      std::cerr << "[ERROR] Cannot open PID matrix file: " << pidMatrixFile << std::endl;
      return;
    }
    std::string line;
    int nLoaded = 0;
    while (std::getline(pin, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::stringstream ss(line);
      std::string tok;
      std::vector<std::string> cols;
      while (std::getline(ss, tok, ',')) cols.push_back(tok);
      if (cols.size() < 7) continue;
      if (cols[0] == "group") continue;

      int g = std::atoi(cols[0].c_str());
      int tlabel = std::atoi(cols[1].c_str());
      if (g < 1 || g > 8 || tlabel < 0 || tlabel > 4) continue;
      auto& mat = pidMats[g - 1];
      for (int r = 0; r < 5; ++r) {
        mat[tlabel * 5 + r] = std::atof(cols[2 + r].c_str());
      }
      ++nLoaded;
    }
    if (nLoaded <= 0) {
      std::cerr << "[ERROR] PID matrix file contained no usable rows." << std::endl;
      return;
    }
  }

  auto pidGroupFromRecoEnergy = [&](double Ereco) {
    int ebin = nE_ana - 1;
    for (int i = 0; i < nE_ana; ++i) {
      if (Ereco < e_edges_ana[i+1]) {
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
  };

  auto trueFlavorLabelFromFinalFlavor = [&](const std::string& finalFlavor) {
    if (finalFlavor == "nue") return 0;
    if (finalFlavor == "numu") return 1;
    if (finalFlavor == "nuebar") return 2;
    if (finalFlavor == "numubar") return 3;
    return 4;
  };

  auto sampleFromRecoFlavorLabel = [&](int recoLabel) {
    if (recoLabel == 0) return 52;
    if (recoLabel == 1) return 50;
    if (recoLabel == 2) return 53;
    if (recoLabel == 3) return 51;
    return -1;
  };

  auto makeHist = [&](const std::string& name, const std::string& title) {
    TH2D* h = new TH2D(name.c_str(), title.c_str(), nE_ana, e_edges_ana, nC_ana, c_edges_ana);
    h->SetDirectory(nullptr);
    return h;
  };

  std::map<std::pair<int,int>, TH2D*> outHists;

  auto sampleTitle = [&](int sample) -> std::string {
    if (sample == 50) return "numuCC_fc";
    if (sample == 51) return "numubarCC_fc";
    if (sample == 52) return "nueCC_fc";
    if (sample == 53) return "nuebarCC_fc";
    return "unknown";
  };

  auto ensureHist = [&](int sample, int reac_mode) -> TH2D* {
    std::pair<int,int> key(sample, reac_mode);
    if (outHists.find(key) == outHists.end()) {
      TString name = TString::Format("h2_sample%d_mode%d", sample, reac_mode);
      outHists[key] = makeHist(name.Data(), sampleTitle(sample));
    }
    return outHists[key];
  };

  TFile fout(outFile.c_str(), "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "[ERROR] Failed to create output file: " << outFile << std::endl;
    return;
  }

  TTree* t = new TTree("valor_dst", "Compact grid-built VALOR DST");
  t->SetDirectory(&fout);

  int sample = 0;
  int beam = 8;
  int reac_mode = 0;
  double enu_true = 0.0;
  double enu_reco = 0.0;
  double costh_true = 0.0;
  double costh_reco = 0.0;
  double Theta_reco = 0.0;
  double Q2_true = 0.0;
  double baseline = 0.0;
  double wght = 0.0;
  double wght_pot = 1.0;

  t->Branch("sample", &sample, "sample/I");
  t->Branch("beam", &beam, "beam/I");
  t->Branch("reac_mode", &reac_mode, "reac_mode/I");
  t->Branch("enu_true", &enu_true, "enu_true/D");
  t->Branch("enu_reco", &enu_reco, "enu_reco/D");
  t->Branch("costh_true", &costh_true, "costh_true/D");
  t->Branch("costh_reco", &costh_reco, "costh_reco/D");
  t->Branch("Theta_reco", &Theta_reco, "Theta_reco/D");
  t->Branch("Q2_true", &Q2_true, "Q2_true/D");
  t->Branch("baseline", &baseline, "baseline/D");
  t->Branch("wght", &wght, "wght/D");
  t->Branch("wght_pot", &wght_pot, "wght_pot/D");

  auto fillTruthComponentIntoAnalysis = [&](int sampleBase,
                                            int reacMode,
                                            TH2D* fluxHist,
                                            TH2D* probHist,
                                            TH1D* XC,
                                            TH1D* XH,
                                            const std::string& finalFlavor) {
    const int nE_fine = fluxHist->GetNbinsX();
    const int nC_fine = fluxHist->GetNbinsY();

    for (int ie_fine = 1; ie_fine <= nE_fine; ++ie_fine) {
      const double E_true = fluxHist->GetXaxis()->GetBinCenter(ie_fine);
      const int ie_xsec = std::min(ie_fine, XC->GetNbinsX());
      const double xsecC = XC->GetBinContent(ie_xsec);
      const double xsecH = XH->GetBinContent(ie_xsec);

      for (int ic_fine = 1; ic_fine <= nC_fine; ++ic_fine) {
        const double C_true = fluxHist->GetYaxis()->GetBinCenter(ic_fine);

        double fluxVal = fluxHist->GetBinContent(ie_fine, ic_fine);
        if (probHist) fluxVal *= probHist->GetBinContent(ie_fine, ic_fine);

        const double fc = getFCRatio(finalFlavor, E_true);
        const double rate = fluxVal * (xsecC * atomsC12 + xsecH * atomsH1)
                          * runtime * solidAngle * unitConversion * fc;
        if (rate <= 0.0) continue;

        std::vector<std::pair<int,double>> eProbs;
        std::vector<std::pair<int,double>> cProbs;
        eProbs.reserve(nE_ana);
        cProbs.reserve(nC_ana);

        if (!applySmearing) {
          int ieReco = -1;
          int icReco = -1;

          for (int i = 1; i <= nE_ana; ++i) {
            if (E_true >= e_edges_ana[i-1] && E_true < e_edges_ana[i]) { ieReco = i; break; }
          }
          if (E_true == e_edges_ana[nE_ana]) ieReco = nE_ana;

          for (int i = 1; i <= nC_ana; ++i) {
            if (C_true >= c_edges_ana[i-1] && C_true < c_edges_ana[i]) { icReco = i; break; }
          }
          if (C_true == c_edges_ana[nC_ana]) icReco = nC_ana;

          if (ieReco > 0 && icReco > 0) {
            eProbs.emplace_back(ieReco, 1.0);
            cProbs.emplace_back(icReco, 1.0);
          }
        } else {
          if (!applyEnergySmearing) {
            int ieReco = -1;
            for (int i = 1; i <= nE_ana; ++i) {
              if (E_true >= e_edges_ana[i-1] && E_true < e_edges_ana[i]) { ieReco = i; break; }
            }
            if (E_true == e_edges_ana[nE_ana]) ieReco = nE_ana;
            if (ieReco > 0) eProbs.emplace_back(ieReco, 1.0);
          } else {
            const double eSigma = std::max(1.0e-6, sigmaEnergy(E_true, finalFlavor));
            const double eMean  = E_true + muEnergyBias(E_true, finalFlavor);

            double eNorm = 0.0;
            for (int ieReco = 1; ieReco <= nE_ana; ++ieReco) {
              const double low = e_edges_ana[ieReco-1];
              const double high = e_edges_ana[ieReco];
              if (high < eMean - smearNSigma * eSigma) continue;
              if (low  > eMean + smearNSigma * eSigma) continue;

              const double p = gaussianBinProb(low, high, eMean, eSigma);
              if (p > 0.0) {
                eProbs.emplace_back(ieReco, p);
                eNorm += p;
              }
            }
            if (preserveNormalization && eNorm > 0.0) {
              for (auto& x : eProbs) x.second /= eNorm;
            }
          }

          if (!applyCosthSmearing) {
            int icReco = -1;
            for (int i = 1; i <= nC_ana; ++i) {
              if (C_true >= c_edges_ana[i-1] && C_true < c_edges_ana[i]) { icReco = i; break; }
            }
            if (C_true == c_edges_ana[nC_ana]) icReco = nC_ana;
            if (icReco > 0) cProbs.emplace_back(icReco, 1.0);
          } else {
            const double cClamped = std::max(-1.0, std::min(1.0, C_true));
            const double theta = std::acos(cClamped);
            const double sigmaThetaRad = std::max(1.0e-6, sigmaThetaDeg(E_true)) * M_PI / 180.0;
            const double cSigma = std::max(1.0e-6, std::abs(std::sin(theta)) * sigmaThetaRad);
            const double cMean  = cClamped;

            double cNorm = 0.0;
            for (int icReco = 1; icReco <= nC_ana; ++icReco) {
              const double low = c_edges_ana[icReco-1];
              const double high = c_edges_ana[icReco];
              if (high < cMean - smearNSigma * cSigma) continue;
              if (low  > cMean + smearNSigma * cSigma) continue;

              const double p = gaussianBinProb(low, high, cMean, cSigma);
              if (p > 0.0) {
                cProbs.emplace_back(icReco, p);
                cNorm += p;
              }
            }
            if (preserveNormalization && cNorm > 0.0) {
              for (auto& x : cProbs) x.second /= cNorm;
            }
          }
        }

        if (eProbs.empty() || cProbs.empty()) continue;

        for (const auto& eitem : eProbs) {
          const int ieReco = eitem.first;
          const double pE = eitem.second;
          if (pE <= 0.0) continue;

          const double eRecoCenter = 0.5 * (e_edges_ana[ieReco - 1] + e_edges_ana[ieReco]);

          for (const auto& citem : cProbs) {
            const int icReco = citem.first;
            const double pC = citem.second;
            if (pC <= 0.0) continue;

            const double cRecoCenter = 0.5 * (c_edges_ana[icReco - 1] + c_edges_ana[icReco]);
            const double smearedWeight = rate * pE * pC;
            if (smearedWeight <= 0.0) continue;

            if (!applyPID) {
              TH2D* h = ensureHist(sampleBase, reacMode);
              h->Fill(eRecoCenter, cRecoCenter, smearedWeight);
              if (eventStyleDST) {
                sample = sampleBase;
                reac_mode = reacMode;
                enu_true = E_true;
                costh_true = C_true;
                enu_reco = eRecoCenter;
                costh_reco = cRecoCenter;
                Theta_reco = costh_reco;
                Q2_true = costh_true;
                baseline = std::sqrt(std::pow(earthRadius * costh_true, 2.0)
                                   + std::pow(neutrinoAtmoHeight, 2.0)
                                   + 2.0 * earthRadius * neutrinoAtmoHeight)
                         - earthRadius * costh_true;
                wght = smearedWeight;
                t->Fill();
              }
            } else {
              const int tlabel = trueFlavorLabelFromFinalFlavor(finalFlavor);
              const int g = pidGroupFromRecoEnergy(eRecoCenter);
              const auto& mat = pidMats[g];

              for (int rlabel = 0; rlabel < 5; ++rlabel) {
                const double pidProb = mat[tlabel * 5 + rlabel];
                if (pidProb <= 0.0) continue;
                const int recoSample = sampleFromRecoFlavorLabel(rlabel);
                if (recoSample < 0) continue;

                TH2D* h = ensureHist(recoSample, reacMode);
                h->Fill(eRecoCenter, cRecoCenter, smearedWeight * pidProb);
                if (eventStyleDST) {
                  sample = recoSample;
                  reac_mode = reacMode;
                  enu_true = E_true;
                  costh_true = C_true;
                  enu_reco = eRecoCenter;
                  costh_reco = cRecoCenter;
                  Theta_reco = costh_reco;
                  Q2_true = costh_true;
                  baseline = std::sqrt(std::pow(earthRadius * costh_true, 2.0)
                                     + std::pow(neutrinoAtmoHeight, 2.0)
                                     + 2.0 * earthRadius * neutrinoAtmoHeight)
                           - earthRadius * costh_true;
                  wght = smearedWeight * pidProb;
                  t->Fill();
                }
              }
            }
          }
        }
      }
    }
  };

  if (!includeCrossComponents) {
    if (applyOsc) {
      TH2D* fluxNumuEff = dynamic_cast<TH2D*>(fluxNumu->Clone("fluxNumuEff"));
      TH2D* fluxNueEff = dynamic_cast<TH2D*>(fluxNue->Clone("fluxNueEff"));
      TH2D* fluxNumubarEff = dynamic_cast<TH2D*>(fluxNumubar->Clone("fluxNumubarEff"));
      TH2D* fluxNuebarEff = dynamic_cast<TH2D*>(fluxNuebar->Clone("fluxNuebarEff"));

      fluxNumuEff->SetDirectory(nullptr);
      fluxNueEff->SetDirectory(nullptr);
      fluxNumubarEff->SetDirectory(nullptr);
      fluxNuebarEff->SetDirectory(nullptr);

      for (int ie = 1; ie <= fluxNumu->GetNbinsX(); ++ie) {
        for (int ic = 1; ic <= fluxNumu->GetNbinsY(); ++ic) {
          fluxNumuEff->SetBinContent(ie, ic,
            fluxNumu->GetBinContent(ie, ic) * PmmNu->GetBinContent(ie, ic) +
            fluxNue->GetBinContent(ie, ic) * PemNu->GetBinContent(ie, ic));
          fluxNueEff->SetBinContent(ie, ic,
            fluxNue->GetBinContent(ie, ic) * PeeNu->GetBinContent(ie, ic) +
            fluxNumu->GetBinContent(ie, ic) * PmeNu->GetBinContent(ie, ic));
          fluxNumubarEff->SetBinContent(ie, ic,
            fluxNumubar->GetBinContent(ie, ic) * PmmNb->GetBinContent(ie, ic) +
            fluxNuebar->GetBinContent(ie, ic) * PemNb->GetBinContent(ie, ic));
          fluxNuebarEff->SetBinContent(ie, ic,
            fluxNuebar->GetBinContent(ie, ic) * PeeNb->GetBinContent(ie, ic) +
            fluxNumubar->GetBinContent(ie, ic) * PmeNb->GetBinContent(ie, ic));
        }
      }

      fillTruthComponentIntoAnalysis(50, 13, fluxNumuEff, nullptr, xNuMuC, xNuMuH, "numu");
      fillTruthComponentIntoAnalysis(52, 60, fluxNueEff, nullptr, xNuEC, xNuEH, "nue");
      fillTruthComponentIntoAnalysis(51, 37, fluxNumubarEff, nullptr, xNuMuBarC, xNuMuBarH, "numubar");
      fillTruthComponentIntoAnalysis(53, 84, fluxNuebarEff, nullptr, xNuEBarC, xNuEBarH, "nuebar");
    } else {
      fillTruthComponentIntoAnalysis(50, 13, fluxNumu, nullptr, xNuMuC, xNuMuH, "numu");
      fillTruthComponentIntoAnalysis(52, 60, fluxNue, nullptr, xNuEC, xNuEH, "nue");
      fillTruthComponentIntoAnalysis(51, 37, fluxNumubar, nullptr, xNuMuBarC, xNuMuBarH, "numubar");
      fillTruthComponentIntoAnalysis(53, 84, fluxNuebar, nullptr, xNuEBarC, xNuEBarH, "nuebar");
    }
  } else {
    fillTruthComponentIntoAnalysis(50, 13,  fluxNumu,    applyOsc ? PmmNu : nullptr, xNuMuC,    xNuMuH,    "numu");
    fillTruthComponentIntoAnalysis(50, 135, fluxNue,     applyOsc ? PemNu : nullptr, xNuMuC,    xNuMuH,    "numu");
    fillTruthComponentIntoAnalysis(52, 60,  fluxNue,     applyOsc ? PeeNu : nullptr, xNuEC,     xNuEH,     "nue");
    fillTruthComponentIntoAnalysis(52, 108, fluxNumu,    applyOsc ? PmeNu : nullptr, xNuEC,     xNuEH,     "nue");
    fillTruthComponentIntoAnalysis(51, 37,  fluxNumubar, applyOsc ? PmmNb : nullptr, xNuMuBarC, xNuMuBarH, "numubar");
    fillTruthComponentIntoAnalysis(51, 188, fluxNuebar,  applyOsc ? PemNb : nullptr, xNuMuBarC, xNuMuBarH, "numubar");
    fillTruthComponentIntoAnalysis(53, 84,  fluxNuebar,  applyOsc ? PeeNb : nullptr, xNuEBarC,  xNuEBarH,  "nuebar");
    fillTruthComponentIntoAnalysis(53, 162, fluxNumubar, applyOsc ? PmeNb : nullptr, xNuEBarC,  xNuEBarH,  "nuebar");
  }

  if (!eventStyleDST) {
    for (const auto& kv : outHists) {
      sample = kv.first.first;
      reac_mode = kv.first.second;
      TH2D* h = kv.second;

      for (int ie = 1; ie <= nE_ana; ++ie) {
        for (int ic = 1; ic <= nC_ana; ++ic) {
          const double v = h->GetBinContent(ie, ic);
          if (v <= 0.0) continue;

          enu_reco = h->GetXaxis()->GetBinCenter(ie);
          costh_reco = h->GetYaxis()->GetBinCenter(ic);
          enu_true = enu_reco;
          costh_true = costh_reco;
          Theta_reco = costh_reco;
          Q2_true = costh_true;
          baseline = std::sqrt(std::pow(earthRadius * costh_true, 2.0)
                             + std::pow(neutrinoAtmoHeight, 2.0)
                             + 2.0 * earthRadius * neutrinoAtmoHeight)
                   - earthRadius * costh_true;
          wght = v;
          t->Fill();
        }
      }
    }
  }

  fout.cd();
  t->Write("", TObject::kOverwrite);
  fout.Write();
  fout.Close();

  if (pOscNu) { pOscNu->Close(); delete pOscNu; }
  if (pOscNb) { pOscNb->Close(); delete pOscNb; }

  std::cout << "[INFO] Wrote compact grid DST: " << outFile << std::endl;
  std::cout << "[INFO] DST mode: "
            << (eventStyleDST ? "event-style (truth preserved)" : "compact template (truth=reco)")
            << std::endl;
  std::cout << "[INFO] Smearing: "
            << (applySmearing ? "ON" : "OFF")
            << " [energy=" << (applyEnergySmearing ? "on" : "off")
            << ", costh=" << (applyCosthSmearing ? "on" : "off")
            << ", nsigma=" << smearNSigma << "]"
            << std::endl;
  std::cout << "[INFO] Normalization preservation: "
            << (preserveNormalization ? "ON" : "OFF (leakage allowed)") << std::endl;
}
ROOTMACRO

root -l -b <<ROOTRUN
.L ${TMP_MACRO}
build_grid_dst_compact_from_maps();
.q
ROOTRUN

rm -f "$TMP_MACRO"
