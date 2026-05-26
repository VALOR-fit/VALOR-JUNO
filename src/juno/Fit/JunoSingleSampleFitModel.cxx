#include "JunoSingleSampleFitModel.h"

#include "Rtypes.h"

#include "valor/Messenger/Messenger.h"
#include "valor/Utilities/UtlXML.h"
#include "valor/Utilities/UtlSystem.h"

using namespace valor;

JunoSingleSampleFitModel::JunoSingleSampleFitModel()
  : GenericSingleSampleFitModel(),
    fUseJunoDetectorResponse(false)
{
}

JunoSingleSampleFitModel::~JunoSingleSampleFitModel()
{
}

bool JunoSingleSampleFitModel::ReadExperimentSpecificConfig(
    TXMLEngine* xml,
    XMLNodePointer_t root_node)
{
  fUseJunoDetectorResponse = false;

  const std::string name = "detector_response/enabled";

  if (util::NodeExists(xml, root_node, name)) {
    fUseJunoDetectorResponse =
      util::XML2Bool(xml, root_node, name);
  }

  if (!fUseJunoDetectorResponse) {
    return true;
  }

  juno::JunoDetectorResponse::Config cfg;

  cfg.apply_fc =
    util::NodeExists(xml, root_node, "detector_response/apply_fc") ?
    util::XML2Bool(xml, root_node, "detector_response/apply_fc") : false;

  cfg.apply_energy_smearing =
  util::NodeExists(xml, root_node, "detector_response/apply_energy_smearing") ?
  util::XML2Bool(xml, root_node, "detector_response/apply_energy_smearing") : false;

  cfg.apply_costh_smearing =
    util::NodeExists(xml, root_node, "detector_response/apply_costh_smearing") ?
    util::XML2Bool(xml, root_node, "detector_response/apply_costh_smearing") : false;

  cfg.apply_pid =
    util::NodeExists(xml, root_node, "detector_response/apply_pid") ?
    util::XML2Bool(xml, root_node, "detector_response/apply_pid") : false;

  cfg.preserve_normalization =
    util::NodeExists(xml, root_node, "detector_response/preserve_normalization") ?
    util::XML2Bool(xml, root_node, "detector_response/preserve_normalization") : false;

  cfg.smear_nsigma =
    util::NodeExists(xml, root_node, "detector_response/smear_nsigma") ?
    util::XML2Dbl(xml, root_node, "detector_response/smear_nsigma") : 5.0;

  fJunoDetectorResponse.SetConfig(cfg);

  if (cfg.apply_fc &&
    util::NodeExists(xml, root_node, "detector_response/fc_coeff_dir")) {
    std::string fcfile =
    util::CompleteConfigPath(
      util::XML2String(xml, root_node, "detector_response/fc_coeff_dir"));

    if (!fJunoDetectorResponse.LoadFCCoefficientsText(fcfile)) {
        VLOG("JunoSingleSampleFitModel", pERROR)
        << "[JUNO DETRESP] Failed to load FC coefficients from: " << fcfile;
        return false;
    }
    }

    if (cfg.apply_pid &&
        util::NodeExists(xml, root_node, "detector_response/pid_matrix_file")) {
    std::string pidfile =
        util::CompleteConfigPath(
        util::XML2String(xml, root_node, "detector_response/pid_matrix_file"));

    if (!fJunoDetectorResponse.LoadPIDMatrixCSV(pidfile)) {
        VLOG("JunoSingleSampleFitModel", pERROR)
        << "[JUNO DETRESP] Failed to load PID file: " << pidfile;
        return false;
    }
  }

  VLOG("JunoSingleSampleFitModel", pNOTICE)
    << "[JUNO DETRESP] Enabled detector response";

  return true;
}

bool JunoSingleSampleFitModel::ApplyExperimentSpecificResponse(
    Option_t* opt)
{
  std::cerr << "[DEBUG JUNO] ApplyExperimentSpecificResponse enabled="
          << fUseJunoDetectorResponse << std::endl;
          
  if (!fUseJunoDetectorResponse) {
    return true;
  }

  return fJunoDetectorResponse.ApplyToVALORArrays(
    this->MCTFinObInternal(),
    this->MCTNullObInternal(),
    this->ReactionModesInternal(),
    this->BinningInternal(),
    this->SampleInternal(),
    opt);
}