#!/usr/bin/env bash
set -euo pipefail

# Build a VALOR DST directly from 400x400 flux/xsec maps.
#
# Key ROOT ownership fixes in this version:
#   - TH1::AddDirectory(kFALSE) to stop auto-registration
#   - all input histograms are cloned into detached local copies
#   - input files are closed immediately after cloning
#   - only detached local histograms are used and deleted manually
#
# Two modes:
# 1) Flux-already-oscillated mode (default): use provided flux maps directly.
# 2) Apply-osc mode: combine unoscillated flux maps with oscillation probability maps on 400x400.

OUTPUT_FILE="data/valor_juno_grid_dst.root"
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
FINE_GRID_MODE=0
APPLY_FC_RATIO=0
FC_COEFF_DIR="data/input/fc_ratio/Total_FC_V2"
APPLY_SMEARING=0
SMEAR_NSIGMA=5.0
APPLY_PID=0
PID_MATRIX_FILE="data/input/detector_response/pid_matrices_8bin_from_user.csv"

RUNTIME_YEARS=6.0
MASS_KG=20000000.0
SOLID_ANGLE=0.0314159265359
CARBON_RATIO=0.88
HYDROGEN_RATIO=0.12

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/build_grid_dst_from_maps.sh [options]

Required flux inputs:
  --flux-numu <file.root>
  --flux-nue <file.root>
  --flux-numubar <file.root>
  --flux-nuebar <file.root>

Options:
  --output <file.root>          Output DST file (default: data/valor_juno_grid_dst.root)
  --xsec-c12 <file.root>        CC xsec C12 ROOT file
  --xsec-h1 <file.root>         CC xsec H1 ROOT file

  --apply-osc                   Apply osc probabilities on 400x400 before rate weighting
  --osc-nu-file <file.root>     Osc prob file for neutrinos (contains NuMuToNuMu3f, NuEToNuMu3f, NuEToNuE3f, NuMuToNuE3f)
  --osc-nb-file <file.root>     Osc prob file for antineutrinos (same histogram names)
  --include-cross-components    Write crossed mode components as separate reac_mode entries
  --fine-grid-mode              Preserve input fine-grid bins in DST output (no 10x10 rebin)
  --apply-fc-ratio              Apply FC-ratio polynomial scaling vs energy
  --fc-coeff-dir <dir>          Directory containing coeffs_*.npy (default: data/input/fc_ratio/Total_FC_V2)
  --apply-smearing              Apply Gaussian smearing in energy and cos(zenith)
  --smear-nsigma <float>        Truncation window in sigma units (default: 5.0)
  --apply-pid                   Apply energy-dependent 5x5 PID matrices (last step)
  --pid-matrix-file <file.csv>  CSV with 8 matrices (groups 1..8), format:
                                group,true_label,p00,p01,p02,p03,p04

  --runtime-years <float>       Default: 6.0
  --mass-kg <float>             Default: 2.0e7
  --solid-angle <float>         Default: 0.0314159265359
  --carbon-ratio <float>        Default: 0.88
  --hydrogen-ratio <float>      Default: 0.12
  --help                        Show this message

Notes:
  - Output tree name is valor_dst with standard VALOR branches.
  - sample is set to JUNO 4CC sample IDs (numu/nue and bars).
  - reac_mode is set to *CCoth coarse modes, consistent with summed modelists.
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
    --fine-grid-mode) FINE_GRID_MODE=1; shift ;;
    --apply-fc-ratio) APPLY_FC_RATIO=1; shift ;;
    --fc-coeff-dir) FC_COEFF_DIR="$2"; shift 2 ;;
    --apply-smearing) APPLY_SMEARING=1; shift ;;
    --smear-nsigma) SMEAR_NSIGMA="$2"; shift 2 ;;
    --apply-pid) APPLY_PID=1; shift ;;
    --pid-matrix-file) PID_MATRIX_FILE="$2"; shift 2 ;;
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

if [[ "$APPLY_SMEARING" == "1" && "$FINE_GRID_MODE" != "1" ]]; then
  echo "[ERROR] --apply-smearing is currently supported only with --fine-grid-mode"
  exit 1
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

TMP_MACRO="/tmp/build_grid_dst_from_maps_$$.C"
cat > "$TMP_MACRO" <<ROOTMACRO
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <utility>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TString.h"

void build_grid_dst_from_maps() {
  TH1::AddDirectory(kFALSE);
  gROOT->cd();

  const std::string outFile = "${OUTPUT_FILE}";
  const std::string fluxNumuFile = "${FLUX_NUMU}";
  const std::string fluxNueFile = "${FLUX_NUE}";
  const std::string fluxNumubarFile = "${FLUX_NUMUBAR}";
  const std::string fluxNuebarFile = "${FLUX_NUEBAR}";
  const std::string xsecC12File = "${XSEC_C12}";
  const std::string xsecH1File = "${XSEC_H1}";

  const bool applyOsc = (${APPLY_OSC} != 0);
  const bool includeCrossComponents = (${INCLUDE_CROSS_COMPONENTS} != 0);
  const bool fineGridMode = (${FINE_GRID_MODE} != 0);
  const bool applyFCRatio = (${APPLY_FC_RATIO} != 0);
  const bool applySmearing = (${APPLY_SMEARING} != 0);
  const double smearNSigma = ${SMEAR_NSIGMA};
  const bool applyPID = (${APPLY_PID} != 0);
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
    if (finalFlavor == "numu") {
      v = evalPoly(fcCoeffNumu, energy);
    } else if (finalFlavor == "nue") {
      v = evalPoly(fcCoeffNue, energy);
    } else if (finalFlavor == "numubar") {
      v = evalPoly(fcCoeffNumubar, energy);
    } else if (finalFlavor == "nuebar") {
      v = evalPoly(fcCoeffNuebar, energy);
    }
    if (v < 0.0) v = 0.0;
    return v;
  };

  auto isEFlavor = [](const std::string& finalFlavor) {
    return (finalFlavor == "nue" || finalFlavor == "nuebar");
  };

  auto sigmaEnergy = [&](double E, const std::string& finalFlavor) {
    const double e = std::max(1.0e-6, E);
    if (isEFlavor(finalFlavor)) {
      return 0.0444 + 0.003280 * std::sqrt(e) + 0.07683 / e - 0.00218 * e;
    }
    return 0.0535 + 0.005934 * std::sqrt(e) + 0.06803 / e - 0.00267 * e;
  };

  auto muEnergyBias = [&](double E, const std::string& finalFlavor) {
    const double e = std::max(1.0e-6, E);
    if (isEFlavor(finalFlavor)) {
      return -0.4557 + 0.235222 * std::sqrt(e) + 0.29617 / e - 0.03182 * e;
    }
    return -0.2657 + 0.160940 * std::sqrt(e) + 0.08111 / e - 0.02438 * e;
  };

  auto sigmaThetaDeg = [&](double E) {
    const double e = std::max(1.0e-6, E);
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
  const double atomsH1 = (massHydrogenG / 1.0) * avogadro;

  auto cloneTH2 = [&](TFile& f, const char* key, const char* newname) -> TH2D* {
    TH2D* h = dynamic_cast<TH2D*>(f.Get(key));
    if (!h) return nullptr;
    TH2D* c = dynamic_cast<TH2D*>(h->Clone(newname));
    if (!c) return nullptr;
    c->SetDirectory(nullptr);
    return c;
  };

  auto cloneTH1 = [&](TFile& f, const char* key, const char* newname) -> TH1D* {
    TH1D* h = dynamic_cast<TH1D*>(f.Get(key));
    if (!h) return nullptr;
    TH1D* c = dynamic_cast<TH1D*>(h->Clone(newname));
    if (!c) return nullptr;
    c->SetDirectory(nullptr);
    return c;
  };

  TFile fNumu(fluxNumuFile.c_str(), "READ");
  TFile fNue(fluxNueFile.c_str(), "READ");
  TFile fNumubar(fluxNumubarFile.c_str(), "READ");
  TFile fNuebar(fluxNuebarFile.c_str(), "READ");
  TFile fXC(xsecC12File.c_str(), "READ");
  TFile fXH(xsecH1File.c_str(), "READ");

  if (fNumu.IsZombie() || fNue.IsZombie() || fNumubar.IsZombie() || fNuebar.IsZombie() || fXC.IsZombie() || fXH.IsZombie()) {
    std::cerr << "[ERROR] Failed to open one or more input files." << std::endl;
    return;
  }

  TH2D* fluxNumu = cloneTH2(fNumu, "custom_flux_hist", "fluxNumu_local");
  TH2D* fluxNue = cloneTH2(fNue, "custom_flux_hist", "fluxNue_local");
  TH2D* fluxNumubar = cloneTH2(fNumubar, "custom_flux_hist", "fluxNumubar_local");
  TH2D* fluxNuebar = cloneTH2(fNuebar, "custom_flux_hist", "fluxNuebar_local");

  if (!fluxNumu || !fluxNue || !fluxNumubar || !fluxNuebar) {
    std::cerr << "[ERROR] Missing custom_flux_hist in one or more flux files." << std::endl;
    delete fluxNumu; delete fluxNue; delete fluxNumubar; delete fluxNuebar;
    return;
  }

  TH1D* xNuMuC = cloneTH1(fXC, "nu_mu_tot_cc", "xNuMuC_local");
  TH1D* xNuEC = cloneTH1(fXC, "nu_e_tot_cc", "xNuEC_local");
  TH1D* xNuMuBarC = cloneTH1(fXC, "nu_mu_bar_tot_cc", "xNuMuBarC_local");
  TH1D* xNuEBarC = cloneTH1(fXC, "nu_e_bar_tot_cc", "xNuEBarC_local");

  TH1D* xNuMuH = cloneTH1(fXH, "nu_mu_tot_cc", "xNuMuH_local");
  TH1D* xNuEH = cloneTH1(fXH, "nu_e_tot_cc", "xNuEH_local");
  TH1D* xNuMuBarH = cloneTH1(fXH, "nu_mu_bar_tot_cc", "xNuMuBarH_local");
  TH1D* xNuEBarH = cloneTH1(fXH, "nu_e_bar_tot_cc", "xNuEBarH_local");

  if (!xNuMuC || !xNuEC || !xNuMuBarC || !xNuEBarC ||
      !xNuMuH || !xNuEH || !xNuMuBarH || !xNuEBarH) {
    std::cerr << "[ERROR] Missing one or more CC xsec histograms." << std::endl;
    delete fluxNumu; delete fluxNue; delete fluxNumubar; delete fluxNuebar;
    delete xNuMuC; delete xNuEC; delete xNuMuBarC; delete xNuEBarC;
    delete xNuMuH; delete xNuEH; delete xNuMuBarH; delete xNuEBarH;
    return;
  }

  // Close inputs early; from here on only detached clones are used.
  fNumu.Close();
  fNue.Close();
  fNumubar.Close();
  fNuebar.Close();
  fXC.Close();
  fXH.Close();

  TH2D *PmmNu = nullptr, *PemNu = nullptr, *PeeNu = nullptr, *PmeNu = nullptr;
  TH2D *PmmNb = nullptr, *PemNb = nullptr, *PeeNb = nullptr, *PmeNb = nullptr;

  if (applyOsc) {
    TFile pOscNu(oscNuFile.c_str(), "READ");
    TFile pOscNb(oscNbFile.c_str(), "READ");
    if (pOscNu.IsZombie() || pOscNb.IsZombie()) {
      std::cerr << "[ERROR] Failed to open oscillation probability file(s)." << std::endl;
      delete fluxNumu; delete fluxNue; delete fluxNumubar; delete fluxNuebar;
      delete xNuMuC; delete xNuEC; delete xNuMuBarC; delete xNuEBarC;
      delete xNuMuH; delete xNuEH; delete xNuMuBarH; delete xNuEBarH;
      return;
    }

    PmmNu = cloneTH2(pOscNu, "NuMuToNuMu3f", "PmmNu_local");
    PemNu = cloneTH2(pOscNu, "NuEToNuMu3f", "PemNu_local");
    PeeNu = cloneTH2(pOscNu, "NuEToNuE3f", "PeeNu_local");
    PmeNu = cloneTH2(pOscNu, "NuMuToNuE3f", "PmeNu_local");

    PmmNb = cloneTH2(pOscNb, "NuMuToNuMu3f", "PmmNb_local");
    PemNb = cloneTH2(pOscNb, "NuEToNuMu3f", "PemNb_local");
    PeeNb = cloneTH2(pOscNb, "NuEToNuE3f", "PeeNb_local");
    PmeNb = cloneTH2(pOscNb, "NuMuToNuE3f", "PmeNb_local");

    if (!PmmNu || !PemNu || !PeeNu || !PmeNu || !PmmNb || !PemNb || !PeeNb || !PmeNb) {
      std::cerr << "[ERROR] Missing required oscillation TH2D in osc files." << std::endl;
      delete fluxNumu; delete fluxNue; delete fluxNumubar; delete fluxNuebar;
      delete xNuMuC; delete xNuEC; delete xNuMuBarC; delete xNuEBarC;
      delete xNuMuH; delete xNuEH; delete xNuMuBarH; delete xNuEBarH;
      delete PmmNu; delete PemNu; delete PeeNu; delete PmeNu;
      delete PmmNb; delete PemNb; delete PeeNb; delete PmeNb;
      return;
    }

    pOscNu.Close();
    pOscNb.Close();
  }

  TFile* fout = TFile::Open(outFile.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "[ERROR] Failed to create output file: " << outFile << std::endl;
    delete fluxNumu; delete fluxNue; delete fluxNumubar; delete fluxNuebar;
    delete xNuMuC; delete xNuEC; delete xNuMuBarC; delete xNuEBarC;
    delete xNuMuH; delete xNuEH; delete xNuMuBarH; delete xNuEBarH;
    delete PmmNu; delete PemNu; delete PeeNu; delete PmeNu;
    delete PmmNb; delete PemNb; delete PeeNb; delete PmeNb;
    return;
  }

  fout->cd();
  TTree* t = new TTree("valor_dst", "Grid-built VALOR DST");

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

  const double neutrinoAtmoHeight = 15.0;
  const double earthRadius = 6341.0;

  const int nE_ana = 10;
  const int nC_ana = 10;
  const double e_edges_ana[nE_ana + 1] = {
    1.08508435, 1.41421357, 1.84317469, 2.40224887, 3.13090217,
    4.08057156, 5.31829591, 6.93144844, 9.03390453, 11.77408038, 20.0
  };
  const double c_edges_ana[nC_ana + 1] = {
    -1.0, -0.845, -0.645, -0.45, -0.225, 0.0, 0.225, 0.45, 0.645, 0.845, 1.0
  };

  std::vector<std::array<double, 25>> pidMats(8);
  for (auto& m : pidMats) m.fill(0.0);
  if (applyPID) {
    std::ifstream pin(pidMatrixFile);
    if (!pin.is_open()) {
      std::cerr << "[ERROR] Cannot open PID matrix file: " << pidMatrixFile << std::endl;
      delete fout; // file not fully owned yet
      delete fluxNumu; delete fluxNue; delete fluxNumubar; delete fluxNuebar;
      delete xNuMuC; delete xNuEC; delete xNuMuBarC; delete xNuEBarC;
      delete xNuMuH; delete xNuEH; delete xNuMuBarH; delete xNuEBarH;
      delete PmmNu; delete PemNu; delete PeeNu; delete PmeNu;
      delete PmmNb; delete PemNb; delete PeeNb; delete PmeNb;
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
      std::cerr << "[ERROR] PID matrix file contained no usable rows: " << pidMatrixFile << std::endl;
      delete fout;
      delete fluxNumu; delete fluxNue; delete fluxNumubar; delete fluxNuebar;
      delete xNuMuC; delete xNuEC; delete xNuMuBarC; delete xNuEBarC;
      delete xNuMuH; delete xNuEH; delete xNuMuBarH; delete xNuEBarH;
      delete PmmNu; delete PemNu; delete PeeNu; delete PmeNu;
      delete PmmNb; delete PemNb; delete PeeNb; delete PmeNb;
      return;
    }
  }

  auto pidGroupFromRecoEnergy = [&](double Ereco) {
    int ebin = 10;
    for (int i = 1; i <= nE_ana; ++i) {
      if (Ereco < e_edges_ana[i]) {
        ebin = i;
        break;
      }
    }
    if (ebin <= 3) return 0;
    if (ebin == 4) return 1;
    if (ebin == 5) return 2;
    if (ebin == 6) return 3;
    if (ebin == 7) return 4;
    if (ebin == 8) return 5;
    if (ebin == 9) return 6;
    return 7;
  };

  auto trueFlavorLabelFromFinalFlavor = [&](const std::string& finalFlavor) {
    if (finalFlavor == "nue") return 0;
    if (finalFlavor == "nuebar") return 1;
    if (finalFlavor == "numu") return 2;
    if (finalFlavor == "numubar") return 3;
    return 4;
  };

  auto sampleFromRecoFlavorLabel = [&](int recoLabel) {
    if (recoLabel == 0) return 52;
    if (recoLabel == 1) return 53;
    if (recoLabel == 2) return 50;
    if (recoLabel == 3) return 51;
    return -1;
  };

  int anaCounter = 0;
  auto buildAnalysisHistogram = [&](TH2D* fluxHist, TH2D* probHist, TH1D* XC, TH1D* XH, const std::string& finalFlavor) {
    TH2D* ana = new TH2D(Form("ana_%d", anaCounter++), "ana", nE_ana, e_edges_ana, nC_ana, c_edges_ana);
    ana->SetDirectory(nullptr);

    const int nE_fine = fluxHist->GetNbinsX();
    const int nC_fine = fluxHist->GetNbinsY();

    for (int ie_fine = 1; ie_fine <= nE_fine; ++ie_fine) {
      double E_center = fluxHist->GetXaxis()->GetBinCenter(ie_fine);
      for (int ic_fine = 1; ic_fine <= nC_fine; ++ic_fine) {
        double C_center = fluxHist->GetYaxis()->GetBinCenter(ic_fine);

        double fluxVal = fluxHist->GetBinContent(ie_fine, ic_fine);
        if (probHist) {
          fluxVal *= probHist->GetBinContent(ie_fine, ic_fine);
        }

        int ie_ana = ana->GetXaxis()->FindBin(E_center);
        int ic_ana = ana->GetYaxis()->FindBin(C_center);

        if (ie_ana >= 1 && ie_ana <= nE_ana && ic_ana >= 1 && ic_ana <= nC_ana) {
          const int ie_xsec = std::min(ie_fine, XC->GetNbinsX());
          const double xsecC = XC->GetBinContent(ie_xsec);
          const double xsecH = XH->GetBinContent(ie_xsec);
          const double fc = getFCRatio(finalFlavor, E_center);
          const double rate = fluxVal * (xsecC * atomsC12 + xsecH * atomsH1) * runtime * solidAngle * unitConversion * fc;
          int globalBin = ana->GetBin(ie_ana, ic_ana);
          ana->AddBinContent(globalBin, rate);
        }
      }
    }
    return ana;
  };

  auto fillFromAnalysis = [&](int s, int rm, TH2D* ana) {
    for (int ie = 1; ie <= nE_ana; ++ie) {
      for (int ic = 1; ic <= nC_ana; ++ic) {
        double rate = ana->GetBinContent(ie, ic);
        if (rate <= 0.0) continue;

        sample = s;
        reac_mode = rm;
        enu_true = ana->GetXaxis()->GetBinCenter(ie);
        costh_true = ana->GetYaxis()->GetBinCenter(ic);
        enu_reco = enu_true;
        costh_reco = costh_true;
        Theta_reco = costh_reco;
        Q2_true = costh_true;
        baseline = std::sqrt(std::pow(earthRadius * costh_true, 2.0)
                           + std::pow(neutrinoAtmoHeight, 2.0)
                           + 2.0 * earthRadius * neutrinoAtmoHeight)
                 - earthRadius * costh_true;
        wght = rate;
        t->Fill();
      }
    }
  };

  auto fillFromFineGrid = [&](int s, int rm, TH2D* fluxHist, TH2D* probHist, TH1D* XC, TH1D* XH, const std::string& finalFlavor) {
    const int nE_fine = fluxHist->GetNbinsX();
    const int nC_fine = fluxHist->GetNbinsY();
    const double sqrt2 = std::sqrt(2.0);

    auto gaussianBinProb = [&](double low, double high, double mean, double sigma) {
      if (sigma <= 0.0) return 0.0;
      const double z1 = (low - mean) / (sqrt2 * sigma);
      const double z2 = (high - mean) / (sqrt2 * sigma);
      return 0.5 * (std::erf(z2) - std::erf(z1));
    };

    if (!applySmearing) {
      for (int ie_fine = 1; ie_fine <= nE_fine; ++ie_fine) {
        const double E_center = fluxHist->GetXaxis()->GetBinCenter(ie_fine);
        const int ie_xsec = std::min(ie_fine, XC->GetNbinsX());
        const double xsecC = XC->GetBinContent(ie_xsec);
        const double xsecH = XH->GetBinContent(ie_xsec);

        for (int ic_fine = 1; ic_fine <= nC_fine; ++ic_fine) {
          const double C_center = fluxHist->GetYaxis()->GetBinCenter(ic_fine);
          double fluxVal = fluxHist->GetBinContent(ie_fine, ic_fine);
          if (probHist) {
            fluxVal *= probHist->GetBinContent(ie_fine, ic_fine);
          }

          const double fc = getFCRatio(finalFlavor, E_center);
          const double rate = fluxVal * (xsecC * atomsC12 + xsecH * atomsH1) * runtime * solidAngle * unitConversion * fc;
          if (rate <= 0.0) continue;

          sample = s;
          reac_mode = rm;
          enu_true = E_center;
          costh_true = C_center;
          enu_reco = enu_true;
          costh_reco = costh_true;
          Theta_reco = costh_reco;
          Q2_true = costh_true;
          baseline = std::sqrt(std::pow(earthRadius * costh_true, 2.0)
                             + std::pow(neutrinoAtmoHeight, 2.0)
                             + 2.0 * earthRadius * neutrinoAtmoHeight)
                   - earthRadius * costh_true;
          wght = rate;
          t->Fill();
        }
      }
      return;
    }

    const int blockE = 8;
    const int blockC = 8;
    const int nE_coarse = (nE_fine + blockE - 1) / blockE;
    const int nC_coarse = (nC_fine + blockC - 1) / blockC;

    auto coarseBinEdgesE = [&](int iCoarse) {
      const int start = (iCoarse - 1) * blockE + 1;
      const int end = std::min(nE_fine, iCoarse * blockE);
      return std::pair<double, double>(
        fluxHist->GetXaxis()->GetBinLowEdge(start),
        fluxHist->GetXaxis()->GetBinUpEdge(end)
      );
    };

    auto coarseBinEdgesC = [&](int iCoarse) {
      const int start = (iCoarse - 1) * blockC + 1;
      const int end = std::min(nC_fine, iCoarse * blockC);
      return std::pair<double, double>(
        fluxHist->GetYaxis()->GetBinLowEdge(start),
        fluxHist->GetYaxis()->GetBinUpEdge(end)
      );
    };

    for (int ie_coarse = 1; ie_coarse <= nE_coarse; ++ie_coarse) {
      const int ie_start = (ie_coarse - 1) * blockE + 1;
      const int ie_end = std::min(nE_fine, ie_coarse * blockE);

      for (int ic_coarse = 1; ic_coarse <= nC_coarse; ++ic_coarse) {
        const int ic_start = (ic_coarse - 1) * blockC + 1;
        const int ic_end = std::min(nC_fine, ic_coarse * blockC);

        double rate = 0.0;
        double wE = 0.0;
        double wC = 0.0;

        for (int ie_fine = ie_start; ie_fine <= ie_end; ++ie_fine) {
          const double E_center_f = fluxHist->GetXaxis()->GetBinCenter(ie_fine);
          const int ie_xsec = std::min(ie_fine, XC->GetNbinsX());
          const double xsecC = XC->GetBinContent(ie_xsec);
          const double xsecH = XH->GetBinContent(ie_xsec);
          const double fc = getFCRatio(finalFlavor, E_center_f);

          for (int ic_fine = ic_start; ic_fine <= ic_end; ++ic_fine) {
            const double C_center_f = fluxHist->GetYaxis()->GetBinCenter(ic_fine);
            double fluxVal = fluxHist->GetBinContent(ie_fine, ic_fine);
            if (probHist) {
              fluxVal *= probHist->GetBinContent(ie_fine, ic_fine);
            }
            const double cellRate = fluxVal * (xsecC * atomsC12 + xsecH * atomsH1) * runtime * solidAngle * unitConversion * fc;
            if (cellRate <= 0.0) continue;
            rate += cellRate;
            wE += cellRate * E_center_f;
            wC += cellRate * C_center_f;
          }
        }

        if (rate <= 0.0) continue;

        const double E_center = wE / rate;
        const double C_center = wC / rate;

        const double eSigma = std::max(1.0e-6, sigmaEnergy(E_center, finalFlavor));
        const double eMean = E_center + muEnergyBias(E_center, finalFlavor);

        const double cClamped = std::max(-1.0, std::min(1.0, C_center));
        const double theta = std::acos(cClamped);
        const double sigmaThetaRad = std::max(1.0e-6, sigmaThetaDeg(E_center)) * M_PI / 180.0;
        const double cSigma = std::max(1.0e-6, std::abs(std::sin(theta)) * sigmaThetaRad);
        const double cMean = cClamped;

        std::vector<std::pair<int, double>> eProbs;
        std::vector<std::pair<int, double>> cProbs;
        eProbs.reserve(16);
        cProbs.reserve(16);

        int ieMinFine = fluxHist->GetXaxis()->FindBin(eMean - smearNSigma * eSigma);
        int ieMaxFine = fluxHist->GetXaxis()->FindBin(eMean + smearNSigma * eSigma);
        ieMinFine = std::max(1, std::min(nE_fine, ieMinFine));
        ieMaxFine = std::max(1, std::min(nE_fine, ieMaxFine));
        if (ieMinFine > ieMaxFine) std::swap(ieMinFine, ieMaxFine);
        int ieMinCoarse = std::max(1, std::min(nE_coarse, (ieMinFine - 1) / blockE + 1));
        int ieMaxCoarse = std::max(1, std::min(nE_coarse, (ieMaxFine - 1) / blockE + 1));

        int icMinFine = fluxHist->GetYaxis()->FindBin(std::max(-1.0, cMean - smearNSigma * cSigma));
        int icMaxFine = fluxHist->GetYaxis()->FindBin(std::min(1.0, cMean + smearNSigma * cSigma));
        icMinFine = std::max(1, std::min(nC_fine, icMinFine));
        icMaxFine = std::max(1, std::min(nC_fine, icMaxFine));
        if (icMinFine > icMaxFine) std::swap(icMinFine, icMaxFine);
        int icMinCoarse = std::max(1, std::min(nC_coarse, (icMinFine - 1) / blockC + 1));
        int icMaxCoarse = std::max(1, std::min(nC_coarse, (icMaxFine - 1) / blockC + 1));

        double eNorm = 0.0;
        for (int ieRecoCoarse = ieMinCoarse; ieRecoCoarse <= ieMaxCoarse; ++ieRecoCoarse) {
          auto edges = coarseBinEdgesE(ieRecoCoarse);
          const double p = gaussianBinProb(edges.first, edges.second, eMean, eSigma);
          if (p > 0.0) {
            eProbs.emplace_back(ieRecoCoarse, p);
            eNorm += p;
          }
        }

        double cNorm = 0.0;
        for (int icRecoCoarse = icMinCoarse; icRecoCoarse <= icMaxCoarse; ++icRecoCoarse) {
          auto edges = coarseBinEdgesC(icRecoCoarse);
          const double p = gaussianBinProb(edges.first, edges.second, cMean, cSigma);
          if (p > 0.0) {
            cProbs.emplace_back(icRecoCoarse, p);
            cNorm += p;
          }
        }

        if (eNorm <= 0.0 || cNorm <= 0.0) continue;

        for (const auto& eitem : eProbs) {
          const bool preserveNormalization = false;

          const int ieRecoCoarse = eitem.first;
          const double pE = preserveNormalization ? (eitem.second / eNorm) : eitem.second;
          if (pE <= 0.0) continue;

          auto eEdges = coarseBinEdgesE(ieRecoCoarse);
          const double eRecoCenter = 0.5 * (eEdges.first + eEdges.second);

          for (const auto& citem : cProbs) {
            const int icRecoCoarse = citem.first;
            const double pC = preserveNormalization ? (citem.second / cNorm) : citem.second;
            if (pC <= 0.0) continue;

            auto cEdges = coarseBinEdgesC(icRecoCoarse);
            const double cRecoCenter = 0.5 * (cEdges.first + cEdges.second);

            const double smearedWeight = rate * pE * pC;
            if (smearedWeight <= 0.0) continue;

            reac_mode = rm;
            enu_true = E_center;
            costh_true = C_center;
            enu_reco = eRecoCenter;
            costh_reco = cRecoCenter;
            Theta_reco = costh_reco;
            Q2_true = costh_true;
            baseline = std::sqrt(std::pow(earthRadius * costh_true, 2.0)
                               + std::pow(neutrinoAtmoHeight, 2.0)
                               + 2.0 * earthRadius * neutrinoAtmoHeight)
                     - earthRadius * costh_true;

            if (!applyPID) {
              sample = s;
              wght = smearedWeight;
              t->Fill();
            } else {
              const int tlabel = trueFlavorLabelFromFinalFlavor(finalFlavor);
              const int g = pidGroupFromRecoEnergy(enu_reco);
              const auto& mat = pidMats[g];

              for (int rlabel = 0; rlabel < 5; ++rlabel) {
                const double pidProb = mat[tlabel * 5 + rlabel];
                if (pidProb <= 0.0) continue;

                const int recoSample = sampleFromRecoFlavorLabel(rlabel);
                if (recoSample < 0) continue;

                sample = recoSample;
                wght = smearedWeight * pidProb;
                t->Fill();
              }
            }
          }
        }
      }
    }
  };

  // sample enums: nue=52, nuebar=53, numu=50, numubar=51
  // intrinsic mode enums: numuCCoth=13, numubarCCoth=37, nueCCoth=60, nuebarCCoth=84
  // crossed mode enums:   nuetonumuCCoth=135, numutonueCCoth=108,
  //                       nuebartonumubarCCoth=188, numubartonuebarCCoth=162

  if (!includeCrossComponents) {
    TH2D* ana50 = nullptr;
    TH2D* ana52 = nullptr;
    TH2D* ana51 = nullptr;
    TH2D* ana53 = nullptr;

    if (applyOsc) {
      TH2D* fluxNumuEff = dynamic_cast<TH2D*>(fluxNumu->Clone("fluxNumuEff"));
      TH2D* fluxNueEff = dynamic_cast<TH2D*>(fluxNue->Clone("fluxNueEff"));
      TH2D* fluxNumubarEff = dynamic_cast<TH2D*>(fluxNumubar->Clone("fluxNumubarEff"));
      TH2D* fluxNuebarEff = dynamic_cast<TH2D*>(fluxNuebar->Clone("fluxNuebarEff"));

      if (!fluxNumuEff || !fluxNueEff || !fluxNumubarEff || !fluxNuebarEff) {
        std::cerr << "[ERROR] Failed to clone effective flux histograms." << std::endl;
        delete fluxNumuEff; delete fluxNueEff; delete fluxNumubarEff; delete fluxNuebarEff;
        delete fluxNumu; delete fluxNue; delete fluxNumubar; delete fluxNuebar;
        delete xNuMuC; delete xNuEC; delete xNuMuBarC; delete xNuEBarC;
        delete xNuMuH; delete xNuEH; delete xNuMuBarH; delete xNuEBarH;
        delete PmmNu; delete PemNu; delete PeeNu; delete PmeNu;
        delete PmmNb; delete PemNb; delete PeeNb; delete PmeNb;
        fout->Write(); fout->Close(); delete fout;
        return;
      }

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

      if (fineGridMode) {
        fillFromFineGrid(50, 13, fluxNumuEff, nullptr, xNuMuC, xNuMuH, "numu");
        fillFromFineGrid(52, 60, fluxNueEff, nullptr, xNuEC, xNuEH, "nue");
        fillFromFineGrid(51, 37, fluxNumubarEff, nullptr, xNuMuBarC, xNuMuBarH, "numubar");
        fillFromFineGrid(53, 84, fluxNuebarEff, nullptr, xNuEBarC, xNuEBarH, "nuebar");
      } else {
        ana50 = buildAnalysisHistogram(fluxNumuEff, nullptr, xNuMuC, xNuMuH, "numu");
        ana52 = buildAnalysisHistogram(fluxNueEff, nullptr, xNuEC, xNuEH, "nue");
        ana51 = buildAnalysisHistogram(fluxNumubarEff, nullptr, xNuMuBarC, xNuMuBarH, "numubar");
        ana53 = buildAnalysisHistogram(fluxNuebarEff, nullptr, xNuEBarC, xNuEBarH, "nuebar");
      }

      delete fluxNumuEff;
      delete fluxNueEff;
      delete fluxNumubarEff;
      delete fluxNuebarEff;
    } else {
      if (fineGridMode) {
        fillFromFineGrid(50, 13, fluxNumu, nullptr, xNuMuC, xNuMuH, "numu");
        fillFromFineGrid(52, 60, fluxNue, nullptr, xNuEC, xNuEH, "nue");
        fillFromFineGrid(51, 37, fluxNumubar, nullptr, xNuMuBarC, xNuMuBarH, "numubar");
        fillFromFineGrid(53, 84, fluxNuebar, nullptr, xNuEBarC, xNuEBarH, "nuebar");
      } else {
        ana50 = buildAnalysisHistogram(fluxNumu, nullptr, xNuMuC, xNuMuH, "numu");
        ana52 = buildAnalysisHistogram(fluxNue, nullptr, xNuEC, xNuEH, "nue");
        ana51 = buildAnalysisHistogram(fluxNumubar, nullptr, xNuMuBarC, xNuMuBarH, "numubar");
        ana53 = buildAnalysisHistogram(fluxNuebar, nullptr, xNuEBarC, xNuEBarH, "nuebar");
      }
    }

    if (!fineGridMode) {
      fillFromAnalysis(50, 13, ana50);
      fillFromAnalysis(52, 60, ana52);
      fillFromAnalysis(51, 37, ana51);
      fillFromAnalysis(53, 84, ana53);
    }

    delete ana50;
    delete ana52;
    delete ana51;
    delete ana53;
  } else {
    if (fineGridMode) {
      fillFromFineGrid(50, 13, fluxNumu, applyOsc ? PmmNu : nullptr, xNuMuC, xNuMuH, "numu");
      fillFromFineGrid(50, 135, fluxNue, applyOsc ? PemNu : nullptr, xNuMuC, xNuMuH, "numu");

      fillFromFineGrid(52, 60, fluxNue, applyOsc ? PeeNu : nullptr, xNuEC, xNuEH, "nue");
      fillFromFineGrid(52, 108, fluxNumu, applyOsc ? PmeNu : nullptr, xNuEC, xNuEH, "nue");

      fillFromFineGrid(51, 37, fluxNumubar, applyOsc ? PmmNb : nullptr, xNuMuBarC, xNuMuBarH, "numubar");
      fillFromFineGrid(51, 188, fluxNuebar, applyOsc ? PemNb : nullptr, xNuMuBarC, xNuMuBarH, "numubar");

      fillFromFineGrid(53, 84, fluxNuebar, applyOsc ? PeeNb : nullptr, xNuEBarC, xNuEBarH, "nuebar");
      fillFromFineGrid(53, 162, fluxNumubar, applyOsc ? PmeNb : nullptr, xNuEBarC, xNuEBarH, "nuebar");
    } else {
      TH2D* ana50_self = buildAnalysisHistogram(fluxNumu, applyOsc ? PmmNu : nullptr, xNuMuC, xNuMuH, "numu");
      TH2D* ana50_cross = buildAnalysisHistogram(fluxNue, applyOsc ? PemNu : nullptr, xNuMuC, xNuMuH, "numu");

      TH2D* ana52_self = buildAnalysisHistogram(fluxNue, applyOsc ? PeeNu : nullptr, xNuEC, xNuEH, "nue");
      TH2D* ana52_cross = buildAnalysisHistogram(fluxNumu, applyOsc ? PmeNu : nullptr, xNuEC, xNuEH, "nue");

      TH2D* ana51_self = buildAnalysisHistogram(fluxNumubar, applyOsc ? PmmNb : nullptr, xNuMuBarC, xNuMuBarH, "numubar");
      TH2D* ana51_cross = buildAnalysisHistogram(fluxNuebar, applyOsc ? PemNb : nullptr, xNuMuBarC, xNuMuBarH, "numubar");

      TH2D* ana53_self = buildAnalysisHistogram(fluxNuebar, applyOsc ? PeeNb : nullptr, xNuEBarC, xNuEBarH, "nuebar");
      TH2D* ana53_cross = buildAnalysisHistogram(fluxNumubar, applyOsc ? PmeNb : nullptr, xNuEBarC, xNuEBarH, "nuebar");

      fillFromAnalysis(50, 13, ana50_self);
      fillFromAnalysis(50, 135, ana50_cross);

      fillFromAnalysis(52, 60, ana52_self);
      fillFromAnalysis(52, 108, ana52_cross);

      fillFromAnalysis(51, 37, ana51_self);
      fillFromAnalysis(51, 188, ana51_cross);

      fillFromAnalysis(53, 84, ana53_self);
      fillFromAnalysis(53, 162, ana53_cross);

      delete ana50_self; delete ana50_cross;
      delete ana52_self; delete ana52_cross;
      delete ana51_self; delete ana51_cross;
      delete ana53_self; delete ana53_cross;
    }
  }

  std::cout << "[INFO] Entries in valor_dst = " << t->GetEntries() << std::endl;

  fout->cd();
  t->Write("valor_dst", TObject::kOverwrite);
  fout->Close();

  delete fluxNumu; delete fluxNue; delete fluxNumubar; delete fluxNuebar;
  delete xNuMuC; delete xNuEC; delete xNuMuBarC; delete xNuEBarC;
  delete xNuMuH; delete xNuEH; delete xNuMuBarH; delete xNuEBarH;
  delete PmmNu; delete PemNu; delete PeeNu; delete PmeNu;
  delete PmmNb; delete PemNb; delete PeeNb; delete PmeNb;
  delete fout;

  std::cout << "[INFO] Wrote grid DST: " << outFile << std::endl;
}
ROOTMACRO

root -l -b <<ROOTRUN
.L ${TMP_MACRO}+
build_grid_dst_from_maps();
.q
ROOTRUN

rm -f "$TMP_MACRO"
