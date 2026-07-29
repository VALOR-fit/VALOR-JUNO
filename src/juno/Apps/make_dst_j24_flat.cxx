//____________________________________________________________________________
/*!

\app       valor-juno_make_dst_j24_flat

\brief     Generate VALOR/JUNO DSTs from the J24 flat truth/reco sample

\options   --input-flat        : J24 flat ROOT file containing truth/reco trees [mandatory]
           --input-mode-friend : mode_friend.root containing modeFriend tree [mandatory]
           --output            : Output DST file [default: valor_juno_j24_flat_dst.root]
           --beam              : Beam configuration, e.g. JUNO_Atmo [mandatory]
           --samples           : Samples to keep, or all [mandatory]
           --require-valid-mode: true/false, skip entries without valid modeFriend match [default: true]
           --mesg-thresholds   : Message thresholds [default: verbosity.xml]

\created   June 2026

*/
//____________________________________________________________________________

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TVector3.h>

// VALOR
#include "valor/Utilities/Env.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/NameRegistry.h"
#include "valor/Conventions/Constants.h"
#include "valor/Conventions/ReacMode.h"
#include "valor/Conventions/Sample.h"
#include "valor/Conventions/BeamConf.h"
#include "valor/Messenger/Messenger.h"
#include "valor/CmdLine/CmdLnArgList.h"
#include "valor/Utilities/UtlNumStat.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Conventions/ModeList.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

using namespace valor;

namespace {

const double neutrino_atmo_height = 15.0; // km
const double earth_radius         = 6341.0; // km
const int    gPrintEveryNEvents   = 2000;

//--------------------------------------------------------------------------------------------------
bool HasBranch(TTree* t, const char* name)
{
  return t && t->GetBranch(name);
}

//--------------------------------------------------------------------------------------------------
template <typename T>
bool BindFirst(TTree* t, T& var, const std::vector<const char*>& names, const char* logical_name)
{
  for (const char* n : names) {
    if (HasBranch(t, n)) {
      t->SetBranchAddress(n, &var);
      VLOG("valor_juno", pNOTICE) << "Using branch " << t->GetName() << ":" << n
                                  << " for " << logical_name;
      return true;
    }
  }

  VLOG("valor_juno", pWARN) << "Could not bind " << logical_name << " from tree "
                            << (t ? t->GetName() : "<null>");
  return false;
}

//--------------------------------------------------------------------------------------------------
double Clamp(double x, double lo, double hi)
{
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

//--------------------------------------------------------------------------------------------------
double GetBaseline(double cosTheta, double R, double h)
{
  // Same convention as the previous JUNO make_dst.
  return std::sqrt(std::pow(R * cosTheta, 2.0) + 2.0 * h * R + h * h) - R * cosTheta;
}

//--------------------------------------------------------------------------------------------------
bool ParseBool(const std::string& s_in, bool default_value = true)
{
  std::string s = s_in;
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);

  if (s == "true"  || s == "1" || s == "yes" || s == "y") return true;
  if (s == "false" || s == "0" || s == "no"  || s == "n") return false;

  return default_value;
}

//--------------------------------------------------------------------------------------------------
valor::conv::Sample_t InferSampleType(bool is_cc, bool is_nc, bool is_nue, bool is_antineutrino)
{
  if (is_cc) {
    if ( is_nue && !is_antineutrino) return valor::conv::ksJUNO_FC_CC_Nue;
    if ( is_nue &&  is_antineutrino) return valor::conv::ksJUNO_FC_CC_Nuebar;
    if (!is_nue && !is_antineutrino) return valor::conv::ksJUNO_FC_CC_Numu;
    if (!is_nue &&  is_antineutrino) return valor::conv::ksJUNO_FC_CC_Numubar;
  }

  if (is_nc) {
    if ( is_nue && !is_antineutrino) return valor::conv::ksJUNO_FC_NC_Nue;
    if ( is_nue &&  is_antineutrino) return valor::conv::ksJUNO_FC_NC_Nuebar;
    if (!is_nue && !is_antineutrino) return valor::conv::ksJUNO_FC_NC_Numu;
    if (!is_nue &&  is_antineutrino) return valor::conv::ksJUNO_FC_NC_Numubar;
  }

  return valor::conv::ksUndefined;
}

//--------------------------------------------------------------------------------------------------
int J24CompactModeCode(int interType, int ctag)
{
  // Matches build_mode_friend.C convention:
  //   CC: 1 QEL, 2 RES, 3 COH, 4 DIS, 5 MEC/2p2h
  //   NC: 101 QEL, 102 RES, 103 COH, 104 DIS, 105 MEC/2p2h
  if (interType < 1 || interType > 5) return 0;
  if (ctag == 2) return interType;
  if (ctag == 1) return 100 + interType;
  return 0;
}

} // namespace

//--------------------------------------------------------------------------------------------------
int processEvents(
    const std::string& inputFlatFile,
    const std::string& inputModeFriendFile,
    const std::vector<conv::Sample_t>& inputSamples,
    const conv::BeamConf_t& inputBeam,
    const std::string& outputFileName,
    bool requireValidMode)
{
  TFile* fin = TFile::Open(inputFlatFile.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    VLOG("valor_juno", pERROR) << "Cannot open input flat file: " << inputFlatFile;
    return 1;
  }

  TTree* truth = dynamic_cast<TTree*>(fin->Get("truth"));
  TTree* reco  = dynamic_cast<TTree*>(fin->Get("reco"));

  if (!truth) truth = dynamic_cast<TTree*>(fin->Get("Truth"));
  if (!reco)  reco  = dynamic_cast<TTree*>(fin->Get("Reco"));

  if (!truth) {
    VLOG("valor_juno", pERROR) << "Cannot find truth/Truth tree in " << inputFlatFile;
    fin->Close();
    return 1;
  }

  if (!reco) {
    VLOG("valor_juno", pWARN) << "Cannot find reco/Reco tree. Will use truth tree for reco fallback.";
    reco = truth;
  }

  TFile* fmode = TFile::Open(inputModeFriendFile.c_str(), "READ");
  if (!fmode || fmode->IsZombie()) {
    VLOG("valor_juno", pERROR) << "Cannot open input mode friend file: " << inputModeFriendFile;
    fin->Close();
    return 1;
  }

  TTree* modeFriend = dynamic_cast<TTree*>(fmode->Get("modeFriend"));
  if (!modeFriend) {
    VLOG("valor_juno", pERROR) << "Cannot find modeFriend tree in " << inputModeFriendFile;
    fmode->Close();
    fin->Close();
    return 1;
  }

  const Long64_t n_truth = truth->GetEntries();
  const Long64_t n_reco  = reco->GetEntries();
  const Long64_t n_mode  = modeFriend->GetEntries();
  const Long64_t n       = std::min(n_truth, std::min(n_reco, n_mode));

  if (n_truth != n_reco || n_truth != n_mode) {
    VLOG("valor_juno", pWARN) << "Entry count mismatch: truth=" << n_truth
                              << ", reco=" << n_reco
                              << ", modeFriend=" << n_mode
                              << ". Will process " << n << " aligned entries.";
  }

  // Input branches: truth.
  int pid = 0;
  int isNC = 0;
  int evtid = 0;
  int file_num = -1;
  int evt = -1;
  double e_true_in = 0.0;
  double theta_true = 0.0;
  double costh_true_in = 999.0;
  double inPx = 0.0, inPy = 0.0, inPz = 0.0;

  bool ok_pid   = BindFirst(truth, pid,       {"pid", "nu_pdg", "pPdg"}, "nu_pdg");
  bool ok_etrue = BindFirst(truth, e_true_in, {"energy", "enu_true", "E_true", "pEn"}, "enu_true");
  bool ok_isNC  = BindFirst(truth, isNC,      {"isNC", "nNC"}, "isNC");

  BindFirst(truth, evtid,     {"evtid", "event_id", "evtID"}, "event_id");
  BindFirst(truth, file_num,  {"file_num", "fileno"}, "file_num");
  BindFirst(truth, evt,       {"evt", "evt_no"}, "evt");
  BindFirst(truth, inPx,      {"inPx", "nuPx"}, "nuPx");
  BindFirst(truth, inPy,      {"inPy", "nuPy"}, "nuPy");
  BindFirst(truth, inPz,      {"inPz", "nuPz"}, "nuPz");

  const bool has_costh_true = BindFirst(truth, costh_true_in,
                                        {"costh_true", "costheta", "cosTheta", "costh"},
                                        "costh_true");
  const bool has_theta_true = has_costh_true ? false :
                              BindFirst(truth, theta_true, {"theta", "theta_true"}, "theta_true");

  if (!ok_pid || !ok_etrue || !ok_isNC) {
    VLOG("valor_juno", pERROR) << "Missing required truth branches. Need pid/energy/isNC equivalents.";
    fmode->Close();
    fin->Close();
    return 1;
  }

  // Input branches: reco.
  // J24 flat reco energy is enu_pn, as used in your resolution scripts.
  // J24 flat reco direction is reco_xyz_pn, a 3-vector. We derive
  // theta_reco = reco_dir.Theta() and costh_reco = cos(theta_reco).
  double e_reco_in = 0.0;
  double theta_reco_in = 0.0;
  double costh_reco_in = 999.0;
  std::vector<double>* reco_xyz_pn = nullptr;

  const bool has_e_reco = BindFirst(reco, e_reco_in,
                                    {"enu_pn", "energy", "enu_reco", "E_reco", "Evis", "evis", "recE"},
                                    "enu_reco");

  const bool has_reco_xyz_pn = HasBranch(reco, "reco_xyz_pn");
  if (has_reco_xyz_pn) {
    reco->SetBranchAddress("reco_xyz_pn", &reco_xyz_pn);
    VLOG("valor_juno", pNOTICE) << "Using branch " << reco->GetName()
                                << ":reco_xyz_pn for reco direction";
  } else {
    VLOG("valor_juno", pWARN) << "Could not bind reco_xyz_pn from tree " << reco->GetName()
                              << ". Will try scalar reco angle/costh branches.";
  }

  const bool has_costh_reco = has_reco_xyz_pn ? false :
                              BindFirst(reco, costh_reco_in,
                                        {"cos_theta_reco", "costh_reco", "costheta", "cosTheta", "costh", "recCosth"},
                                        "costh_reco");
  const bool has_theta_reco = (has_reco_xyz_pn || has_costh_reco) ? false :
                              BindFirst(reco, theta_reco_in,
                                        {"theta_reco", "theta", "recTheta"},
                                        "theta_reco");

  // Input branches: mode friend.
  int valid_mode = 0;
  int interType = 0;
  int ctag = 0;
  int mode_code = 0;

  BindFirst(modeFriend, valid_mode, {"valid"}, "valid_mode");
  BindFirst(modeFriend, interType,  {"interType"}, "interType");
  BindFirst(modeFriend, ctag,       {"ctag"}, "ctag");
  BindFirst(modeFriend, mode_code,  {"mode"}, "mode_code");

  // Output tree: keep the same external structure as the old make_dst.
  TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    VLOG("valor_juno", pERROR) << "Cannot create output file: " << outputFileName;
    fmode->Close();
    fin->Close();
    return 1;
  }

  TTree* outputTree = new TTree("juno_true_tree", "Processed Data");

  double Q2 = 0., W = 0., E_reco = 0., E_true = 0., costh_true = 0., costh_reco = 0., wght = 0.;
  double theta_true_out = 0.0;
  double theta_reco_out = 0.0;
  double L_true = 0., L_reco = 0.;
  bool isnc = 0, iscc = 0;
  int reac_mod_wo_PID = -1, reac_mod = -1, nu_pdg = -1, event_id = 0, sample = -1, beam = -1;

  outputTree->Branch("Q2", &Q2);
  outputTree->Branch("W", &W);
  outputTree->Branch("enu_true", &E_true);
  outputTree->Branch("enu_reco", &E_reco);
  outputTree->Branch("costh_true", &costh_true);
  outputTree->Branch("costh_reco", &costh_reco);
  outputTree->Branch("theta_true", &theta_true_out);
  outputTree->Branch("theta_reco", &theta_reco_out);
  outputTree->Branch("L_true", &L_true);
  outputTree->Branch("L_reco", &L_reco);
  outputTree->Branch("nCC", &iscc);
  outputTree->Branch("nNC", &isnc);
  outputTree->Branch("reac_mode_wo_PID", &reac_mod_wo_PID);
  outputTree->Branch("reac_mode", &reac_mod);
  outputTree->Branch("nu_pdg", &nu_pdg);
  outputTree->Branch("event_id", &event_id);
  outputTree->Branch("wght", &wght);
  outputTree->Branch("sample", &sample);
  outputTree->Branch("beam", &beam);

  // Debug/provenance branches. Existing VALOR code can ignore them.
  int flat_file_num = -1;
  int flat_evt = -1;
  int j24_interType = 0;
  int j24_ctag = 0;
  int j24_mode = 0;
  int mode_valid = 0;

  outputTree->Branch("flat_file_num", &flat_file_num);
  outputTree->Branch("flat_evt", &flat_evt);
  outputTree->Branch("j24_interType", &j24_interType);
  outputTree->Branch("j24_ctag", &j24_ctag);
  outputTree->Branch("j24_mode", &j24_mode);
  outputTree->Branch("mode_valid", &mode_valid);

  Long64_t savedEvents = 0;
  Long64_t n_bad_mode = 0;
  Long64_t n_bad_flavour = 0;
  Long64_t n_bad_reco = 0;
  Long64_t n_bad_sample = 0;

  for (Long64_t i = 0; i < n; ++i) {
    truth->GetEntry(i);
    reco->GetEntry(i);
    modeFriend->GetEntry(i);

    if (gPrintEveryNEvents && !(i % gPrintEveryNEvents)) {
      VLOG("valor_juno", pNOTICE) << "** Reading event: " << i
                                  << " - so far accepted " << savedEvents;
    }

    if (!(pid == 12 || pid == -12 || pid == 14 || pid == -14)) {
      ++n_bad_flavour;
      continue;
    }

    if (requireValidMode && valid_mode == 0) {
      ++n_bad_mode;
      continue;
    }

    E_true = e_true_in;
    if (has_costh_true) {
      costh_true = costh_true_in;
    } else if (has_theta_true) {
      costh_true = std::cos(theta_true);
    } else if (E_true > 0.0) {
      // Same convention as the old GHEP make_dst: costh_true = -Pz/E.
      costh_true = -inPz / E_true;
    } else {
      costh_true = 0.0;
    }
    costh_true = Clamp(costh_true, -1.0, 1.0);

    if (has_theta_true) {
      theta_true_out = theta_true;
    } else {
      theta_true_out = std::acos(Clamp(costh_true, -1.0, 1.0));
    }

    if (has_e_reco) E_reco = e_reco_in;
    else            E_reco = E_true;

    if (has_reco_xyz_pn && reco_xyz_pn && reco_xyz_pn->size() >= 3) {
      const double x_reco = reco_xyz_pn->at(0);
      const double y_reco = reco_xyz_pn->at(1);
      const double z_reco = reco_xyz_pn->at(2);

      if (std::isfinite(x_reco) && std::isfinite(y_reco) && std::isfinite(z_reco)) {
        TVector3 reco_dir(x_reco, y_reco, z_reco);
        if (reco_dir.Mag() > 0.0) {
          theta_reco_out = reco_dir.Theta();
          costh_reco = std::cos(theta_reco_out);
        } else {
          theta_reco_out = std::acos(Clamp(costh_true, -1.0, 1.0));
          costh_reco = costh_true;
        }
      } else {
        theta_reco_out = std::acos(Clamp(costh_true, -1.0, 1.0));
        costh_reco = costh_true;
      }
    } else if (has_costh_reco) {
      costh_reco = costh_reco_in;
      theta_reco_out = std::acos(Clamp(costh_reco, -1.0, 1.0));
    } else if (has_theta_reco) {
      theta_reco_out = theta_reco_in;
      costh_reco = std::cos(theta_reco_out);
    } else {
      costh_reco = costh_true;
      theta_reco_out = std::acos(Clamp(costh_reco, -1.0, 1.0));
    }
    costh_reco = Clamp(costh_reco, -1.0, 1.0);

    if (!std::isfinite(E_true) || !std::isfinite(E_reco) || !std::isfinite(costh_reco) || E_true <= 0.0 || E_reco < 0.0) {
      ++n_bad_reco;
      continue;
    }

    // Leave Q2 and W at zero for now. They are not stored directly in the J24 flat/detsim files.
    Q2 = 0.0;
    W  = 0.0;

    nu_pdg = pid;
    isnc = (isNC != 0);
    iscc = !isnc;

    j24_interType = interType;
    j24_ctag = ctag;
    j24_mode = (mode_code != 0) ? mode_code : J24CompactModeCode(interType, ctag);
    mode_valid = valid_mode;

    // For now, store the J24 compact mode code in the existing reac_mode branches.
    // This preserves channel information even before mapping to a VALOR/GENIE ReacMode_t enum.
    reac_mod = j24_mode;
    reac_mod_wo_PID = reac_mod;

    bool is_nue = (std::abs(nu_pdg) == 12);
    bool is_antineutrino = (nu_pdg < 0);
    sample = static_cast<int>(InferSampleType(iscc, isnc, is_nue, is_antineutrino));

    int cnt = count(inputSamples.begin(), inputSamples.end(), static_cast<conv::Sample_t>(sample));
    if (sample == static_cast<int>(conv::ksUndefined)) {
      ++n_bad_sample;
      continue;
    } else if (cnt == 0) {
      ++n_bad_sample;
      continue;
    }

    event_id = evtid;
    flat_file_num = file_num;
    flat_evt = evt;

    L_true = GetBaseline(costh_true, earth_radius, neutrino_atmo_height);
    L_reco = GetBaseline(costh_reco, earth_radius, neutrino_atmo_height);

    beam = static_cast<int>(inputBeam);

    // The flat sample is already produced with detector response. Do not downsample here.
    // Put flux/xsec/exposure normalisation in the VALOR configuration/weights instead.
    wght = 1.0;

    outputTree->Fill();
    ++savedEvents;
  }

  outputFile->cd();
  outputTree->Write();
  outputFile->Close();

  fmode->Close();
  fin->Close();
  delete outputFile;

  VLOG("valor_juno", pNOTICE) << "Processed aligned entries : " << n;
  VLOG("valor_juno", pNOTICE) << "Saved events              : " << savedEvents;
  VLOG("valor_juno", pNOTICE) << "Skipped invalid mode      : " << n_bad_mode;
  VLOG("valor_juno", pNOTICE) << "Skipped bad flavour       : " << n_bad_flavour;
  VLOG("valor_juno", pNOTICE) << "Skipped bad reco/energy   : " << n_bad_reco;
  VLOG("valor_juno", pNOTICE) << "Skipped bad sample        : " << n_bad_sample;
  VLOG("valor_juno", pNOTICE) << "Wrote                     : " << outputFileName;

  return 0;
}

//--------------------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  time_t scriptstartrawtime;
  struct tm* scriptstarttimeinfo;
  time(&scriptstartrawtime);
  scriptstarttimeinfo = localtime(&scriptstartrawtime);

  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Script start local time and date:  " << asctime(scriptstarttimeinfo);
  std::cout << "-----------------------------------------------------------" << std::endl;

  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-h") {
      VLOG("valor_juno", pFATAL)
        << "Usage: valor-juno_make_dst_j24_flat "
        << " --input-flat path_to_j24_flat.root "
        << " --input-mode-friend path_to_mode_friend.root "
        << " --output path_to_output.root "
        << " --beam JUNO_Atmo "
        << " --samples all "
        << " --require-valid-mode true "
        << " --mesg-thresholds verbosity.xml";
      return 0;
    }
  }

  CmdLnArgList arg;
  arg.Read(argc, argv, "$VALOR_EXPT/config/apps/make_dst_j24_flat.arg");
  VLOG("valor_juno", pNOTICE) << arg;

  string arg_input_flat              = arg.GetString("input-flat");
  string arg_input_mode_friend       = arg.GetString("input-mode-friend");
  string arg_output                  = arg.GetString("output");
  string arg_input_beam              = arg.GetString("beam");
  vector<string> arg_input_samples   = arg.GetStringVec("samples");
  string arg_require_valid_mode      = arg.GetString("require-valid-mode");

  std::string thresholds = util::CompleteConfigPath(arg.GetString("mesg-thresholds"));
  valor::Messenger::Instance()->Configure(thresholds);

  bool requireValidMode = ParseBool(arg_require_valid_mode, true);

  vector<string> samples;
  vector<conv::Sample_t> inputSamples;

  if (std::find(arg_input_samples.begin(), arg_input_samples.end(), "all") != arg_input_samples.end()) {
    samples = {
      "numuCC_fc", "numubarCC_fc", "nueCC_fc", "nuebarCC_fc",
      "numuNC_fc", "numubarNC_fc", "nueNC_fc", "nuebarNC_fc"
    };

    for (string sample_name : samples) {
      inputSamples.push_back(util::SampleFromString(sample_name));
    }
  } else {
    for (string sample_name : arg_input_samples) {
      conv::Sample_t curr_sample = util::SampleFromString(sample_name);
      assert(curr_sample != conv::ksUndefined);
      inputSamples.push_back(curr_sample);
    }
  }

  conv::BeamConf_t inputBeam = util::BeamConfFromString(arg_input_beam);
  assert(inputBeam != conv::kbUndefined);

  int ret = processEvents(
      arg_input_flat,
      arg_input_mode_friend,
      inputSamples,
      inputBeam,
      arg_output,
      requireValidMode);

  VLOG("valor_juno", pNOTICE) << "Done!";

  time_t endrawtime;
  struct tm* endtimeinfo;
  time(&endrawtime);
  endtimeinfo = localtime(&endrawtime);
  time_t total_time = endrawtime - scriptstartrawtime;

  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End local time and date  :  " << asctime(endtimeinfo);
  std::cout << " Total time taken         :  "
            << total_time / static_cast<double>(60) << " minutes" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return ret;
}