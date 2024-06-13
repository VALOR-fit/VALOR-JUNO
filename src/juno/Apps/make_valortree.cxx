#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TSystem.h>
#include <TRandom.h>
#include "valor/Utilities/Env.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/NameRegistry.h"
#include "valor/Conventions/Constants.h"
#include "valor/Conventions/ReacMode.h"
#include "valor/Messenger/Messenger.h"
#include "valor/CmdLine/CmdLnArgList.h"
#include "valor/Utilities/UtlNumStat.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Conventions/ModeList.h"
#include "valor/Conventions/ParameterPool.h"
#include "valor/GENIEUtilities/GENIEIncludes.h"
#include "valor/GENIEUtilities/UtlEventRec.h"
#include "valor/Utilities/SingleSampleKineBinning.h"
#include "valor/Utilities/RFL.h"
#include "valor/Utilities/KinematicBin.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/RunOpt.h"
#include "valor/Oscillation/NuOsc3Flavour.h"
#include "valor/Conventions/Sample.h"
#include "valor/Conventions/BeamConf.h"

using std::string;
using namespace genie;
using namespace valor; 

int processEvents(const std::vector<std::string>& inputFiles, const std::string& outputFileName, const string& sampleName, bool otherSample, bool calculateAppearance) {
  TChain* ghep_tree = new TChain("gtree");
  static const int gPrintEveryNEvents = 500000;
  ghep_tree->SetMaxVirtualSize(-1);
  for (const auto& file : inputFiles) {
    ghep_tree->Add(file.c_str());
  }
  
  genie::NtpMCEventRecord *mcrec = nullptr;
  ghep_tree->SetBranchAddress("gmcrec", &mcrec);
  
  TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
  TTree* outputTree = new TTree("juno_true_tree", "Processed Data");
  double Q2, W, E_true, costh_true, costh_reco, wght;
  int isnc = 0, iscc = 0, nkaon = 0, npip = 0, npim = 0, npi0 = 0, reac_mod = -1, nu_pdg = -1, event_id = 0, sample_branch = -1, beam = -1;
  outputTree->Branch("Q2", &Q2);
  outputTree->Branch("W", &W);
  outputTree->Branch("enu_true", &E_true);
  outputTree->Branch("costh_true", &costh_true);
  outputTree->Branch("costh_reco", &costh_reco);
  outputTree->Branch("nCC", &isnc);
  outputTree->Branch("nNC", &iscc);
  outputTree->Branch("nKaon", &nkaon);
  outputTree->Branch("nPionp", &npip);
  outputTree->Branch("nPionm", &npim);
  outputTree->Branch("nPion0", &npi0);
  outputTree->Branch("reac_mode", &reac_mod);
  outputTree->Branch("nu_PDG", &nu_pdg);
  outputTree->Branch("event_id", &event_id);
  outputTree->Branch("wght", &wght); 
  outputTree->Branch("sample", &sample_branch);  
  outputTree->Branch("beam", &beam);  

  NuOsc3Flavour osc; // Create NuOsc3Flavour instance
  osc.SetBaseline(12742); // Set baseline distance in kilometers
  osc.SetDensity(2.6); // Set density in g/cm^3
  osc.SetMassHierarchy(valor::conv::NuMassHierarchy_t::kNmhNormal); // Set mass hierarchy
  osc.SetOscParameters(0.846, 0.093, 0.528, 7.53e-5, 2.44e-3, 1.36); // Set oscillation parameters

  TRandom randomGen; // Create random number generator instance

  Long64_t nentries = ghep_tree->GetEntries();

  conv::Sample_t sample = util::SampleFromString(sampleName);
 
  sample_branch = sample;
  beam = valor::conv::kbAtmo_JUNO; 

  int savedEvents = 0; 
  for (Long64_t i = 0; i < nentries; ++i) {
    ghep_tree->GetEntry(i);
    if (gPrintEveryNEvents && !(i % gPrintEveryNEvents)) {
      VLOG("valor_juno", pNOTICE) << "** Reading event: " << i;
    }

    if (mcrec && mcrec->event) {
      const genie::EventRecord& event = *(mcrec->event);
      const TLorentzVector& nuP4 = *(event.Probe()->P4());
      const genie::Interaction* interaction = event.Summary();
      const genie::GHepParticle* neutrino = event.Probe();
      iscc = 0; isnc = 0; npip = 0; npim = 0; npi0 = 0; nkaon = 0; reac_mod = 0; event_id = 0;

      if (neutrino != nullptr) {
        nu_pdg = neutrino->Pdg();
      }
      if (interaction->ProcInfo().IsWeakCC()) {
        iscc = 1;
        switch (nu_pdg) {
          case 12:
            event_id = 121;
            break;
          case -12:
            event_id = -121;
            break;
          case 14:
            event_id = 141;
            break;
          case -14:
            event_id = -141;
            break;
        }
      } else if (interaction->ProcInfo().IsWeakNC()) {
        isnc = 1;
        switch (nu_pdg) {
          case 12:
            event_id = 120;
            break;
          case -12:
            event_id = -120;
            break;
          case 14:
            event_id = 140;
            break;
          case -14:
            event_id = -140;
            break;
        }
      }

      TIter iter(event.MakeIterator());
      genie::GHepParticle* p = nullptr;
      while ((p = dynamic_cast<genie::GHepParticle*>(iter.Next()))) {
        if (p->Status() == genie::kIStStableFinalState) {
          int pdg = p->Pdg();
          if (pdg == genie::kPdgPiP) npip++;
          else if (pdg == genie::kPdgPiM) npim++;
          else if (pdg == genie::kPdgPi0) npi0++;
          if (pdg == genie::kPdgKP || pdg == genie::kPdgKM || pdg == genie::kPdgK0 || pdg == genie::kPdgAntiK0) {
            nkaon++;
          }
        }
      }

      valor::conv::ReacMode_t reac_mode = valor::util::DetermineReactionMode(&event, false);
      reac_mod = reac_mode;

      E_true = nuP4.Energy();
      VLOG("valor_juno", pDEBUG) << "Read event: "
                                 << "Etrue = " << E_true << " GeV; ";

      costh_true = -nuP4.Pz() / E_true;
      Q2 = mcrec->event->Summary()->Kine().Q2(true);
      W = mcrec->event->Summary()->Kine().W(true);

      // Calculate survival and transition probabilities
      double osc_probs[3][3];
      if (nu_pdg > 0) { // Neutrino
        osc.CalcMatterOscProbAll(false, E_true, osc_probs);
      } else { // Antineutrino
        osc.CalcMatterOscProbAll(true, E_true, osc_probs);
      }

      double prob_threshold;
      if (calculateAppearance) { // Calculate appearance probability
        if (nu_pdg == 14) {
          prob_threshold = osc_probs[0][1]; // P(ν_μ -> ν_e)
        } else if (nu_pdg == -14) {
          prob_threshold = osc_probs[0][1]; // P(ν̄_μ -> ν̄_e)
        } else {
          prob_threshold = 0.0;
        }
      } else { // Calculate survival probability
        switch (nu_pdg) {
          case 12:
            prob_threshold = osc_probs[0][0]; // P(ν_e -> ν_e)
            break;
          case -12:
            prob_threshold = osc_probs[0][0]; // P(ν̄_e -> ν̄_e)
            break;
          case 14:
            prob_threshold = osc_probs[1][1]; // P(ν_μ -> ν_μ)
            break;
          case -14:
            prob_threshold = osc_probs[1][1]; // P(ν̄_μ -> ν̄_μ)
            break;
          default:
            prob_threshold = 0.0;
            break;
        }
      }

      double random_value = randomGen.Uniform(0.0, 1.0); // Generate a random number between 0 and 1
      if (random_value <= prob_threshold) {
        wght = 1.0;  
        outputTree->Fill(); 
        savedEvents++; 
      }
    }
  }

  outputTree->Write();
  outputFile->Close();

  delete ghep_tree;
  delete outputFile; // Also deletes outputTree
  VLOG("valor_juno", pNOTICE) << "Done!";
  std::cout << "Total saved events: " << savedEvents << std::endl; // 输出保存的事件数量

  return 0;
}

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cout << "Usage: " << argv[0] << " -s|-a <sample_name> <input files> <output file name>" << std::endl;
    return 1;
  }

  bool calculateAppearance = false;
  bool otherSample = false;
  int inputStartIndex = 2;

  if (string(argv[1]) == "-a") {
    calculateAppearance = true;
    inputStartIndex = 3;
  } else if (string(argv[1]) == "-s") {
    calculateAppearance = false;
    inputStartIndex = 3;
  } else {
    std::cout << "Usage: " << argv[0] << " -s|-a <sample_name> <input files> <output file name>" << std::endl;
    return 1;
  }

  std::string sampleName = argv[2];
  std::vector<std::string> inputFiles;
  for (int i = inputStartIndex; i < argc - 1; ++i) {
    inputFiles.emplace_back(argv[i]);
  }
  std::string outputFileName = argv[argc - 1];
  return processEvents(inputFiles, outputFileName, sampleName, otherSample, calculateAppearance);
}
