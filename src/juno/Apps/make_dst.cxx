#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include <sstream>
#include <map>
#include <utility>
#include <fstream>
#include <cstdlib>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <TColor.h>

//#include <ROOT/RDF/RInterface.hxx>
//#include <ROOT/RDFHelpers.hxx>
//#include <ROOT/RDataFrame.hxx>
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

struct NormEntry {
    double energyLow;
    double costhLow;
    double numu;
    double numubar;
    double nue;
    double nuebar;
    double nc;
};

std::vector<NormEntry> gNormTable;

// Function to read normalization table from txt file
bool readNormTable(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error opening normalization file: " << filename << std::endl;
        return false;
    }
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        // Skip header or comment lines
        if (line[0] == '#' || line.find("----") != std::string::npos) continue;
        std::istringstream iss(line);
        NormEntry entry;
        if (!(iss >> entry.energyLow >> entry.costhLow >> entry.numu >> entry.numubar >> entry.nue >> entry.nuebar >> entry.nc)) {
            continue;
        }
        gNormTable.push_back(entry);
    }
    return true;
}

std::string SampleTypeToString(valor::conv::Sample_t sample) {
    switch (sample) {
        case valor::conv::ksNuMuCC: return "numuCC";
        case valor::conv::ksNuMuCC1piC: return "numuCC1piC";
        case valor::conv::ksNuMuCC1pi0: return "numuCC1pi0";
        case valor::conv::ksNuMuCC1piC1pi0: return "numuCC1piC1pi0";
        case valor::conv::ksNuMuCCoth: return "numuCCoth";
        case valor::conv::ksNuECC0pi: return "nueCC0pi";
        case valor::conv::ksNuECC1piC: return "nueCC1piC";
        case valor::conv::ksNuECC1pi0: return "nueCC1pi0";
        case valor::conv::ksNuECCoth: return "nueCCoth";
        case valor::conv::ksNuMuBarCChigh: return "numubarCChigh";
        case valor::conv::ksNuMuBarCClow: return "numubarCClow";
        case valor::conv::ksNuEBarCChigh: return "nuebarCChigh";
        case valor::conv::ksNuEBarCClow: return "nuebarCClow";
        default: return "undefined";
    }
}

valor::conv::Sample_t InferSampleType(int npi, int npi0, bool is_nue, bool is_antineutrino) {
    if (is_nue) {
        if (is_antineutrino) {
            return valor::conv::ksNuEBarCChigh; 
        } else {
            if (npi == 0 && npi0 == 0) return valor::conv::ksNuECC;
            if (npi == 1 && npi0 == 0) return valor::conv::ksNuECC1piC;
            if (npi == 0 && npi0 == 1) return valor::conv::ksNuECC1pi0;
            return valor::conv::ksNuECCoth;
        }
    } else {
        if (is_antineutrino) {
            return valor::conv::ksNuMuBarCChigh; 
        } else {
            if (npi == 0 && npi0 == 0) return valor::conv::ksNuMuCC;
            if (npi == 1 && npi0 == 0) return valor::conv::ksNuMuCC1piC;
            if (npi == 0 && npi0 == 1) return valor::conv::ksNuMuCC1pi0;
            if (npi == 1 && npi0 == 1) return valor::conv::ksNuMuCC1piC1pi0;
            return valor::conv::ksNuMuCCoth;
        }
    }
}

#include <unordered_map>
// ... (other includes)

int processEvents(const std::vector<std::string>& inputFiles, const std::string& outputFileName) {
    TChain* ghep_tree = new TChain("gtree");
    static const int gPrintEveryNEvents = 200000;
    ghep_tree->SetMaxVirtualSize(-1);
    for (const auto& file : inputFiles) {
        ghep_tree->Add(file.c_str());
    }
    genie::NtpMCEventRecord* mcrec = nullptr;
    ghep_tree->SetBranchAddress("gmcrec", &mcrec);
    Long64_t nentries = ghep_tree->GetEntries();

    std::vector<int> binCounts(gNormTable.size(), 0);
    for (Long64_t i = 0; i < nentries; ++i) {
        ghep_tree->GetEntry(i);
        if (mcrec && mcrec->event) {
            const TLorentzVector& nuP4 = *(mcrec->event->Probe()->P4());
            double E_true_tmp = nuP4.Energy();
            double costh_true_tmp = -nuP4.Pz() / E_true_tmp;
            int binIndex = -1;
            double bestEnergy = -1;
            double bestCosth = -100;
            for (size_t j = 0; j < gNormTable.size(); ++j) {
                const NormEntry& entry = gNormTable[j];
                if (E_true_tmp >= entry.energyLow && costh_true_tmp >= entry.costhLow) {
                    if (entry.energyLow > bestEnergy || (entry.energyLow == bestEnergy && entry.costhLow > bestCosth)) {
                        bestEnergy = entry.energyLow;
                        bestCosth = entry.costhLow;
                        binIndex = j;
                    }
                }
            }
            if (binIndex != -1) {
                binCounts[binIndex]++;
            }
        }
    }

    TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    TTree* outputTree = new TTree("juno_true_tree", "Processed Data");
    double Q2, W, E_true, costh_true, costh_reco, wght, LoE_true;
    int isnc = 0, iscc = 0, nkaon = 0, npip = 0, npim = 0, npi0 = 0, reac_mod_wo_PID = -1;
    int reac_mod = -1, nu_pdg = -1, event_id = 0, sample_branch = -1, beam = -1;
    outputTree->Branch("Q2", &Q2);
    outputTree->Branch("W", &W);
    outputTree->Branch("enu_true", &E_true);
    outputTree->Branch("LoverE_true", &LoE_true);
    outputTree->Branch("costh_true", &costh_true);
    outputTree->Branch("costh_reco", &costh_reco);
    outputTree->Branch("nCC", &iscc);
    outputTree->Branch("nNC", &isnc);
    outputTree->Branch("nKaon", &nkaon);
    outputTree->Branch("nPionp", &npip);
    outputTree->Branch("nPionm", &npim);
    outputTree->Branch("nPion0", &npi0);
    outputTree->Branch("reac_mode_wo_PID", &reac_mod_wo_PID);
    outputTree->Branch("reac_mode", &reac_mod);
    outputTree->Branch("nu_pdg", &nu_pdg);
    outputTree->Branch("event_id", &event_id);
    outputTree->Branch("wght", &wght);
    outputTree->Branch("sample", &sample_branch);
    outputTree->Branch("beam", &beam);

    // Map to store counts for each reaction mode
    std::unordered_map<int, int> modeCounts;

    int savedEvents = 0;
    TRandom3 randGen;
    for (Long64_t i = 0; i < nentries; ++i) {
        ghep_tree->GetEntry(i);
        if (gPrintEveryNEvents && !(i % gPrintEveryNEvents)) {
            VLOG("valor_juno", pNOTICE) << "** Reading event: " << i;
        }
        if (mcrec && mcrec->event) {
            const TLorentzVector& nuP4 = *(mcrec->event->Probe()->P4());
            E_true = nuP4.Energy();
            costh_true = -nuP4.Pz() / E_true;
            int binIndex = -1;
            double bestEnergy = -1;
            double bestCosth = -100;
            for (size_t j = 0; j < gNormTable.size(); ++j) {
                const NormEntry& entry = gNormTable[j];
                if (E_true >= entry.energyLow && costh_true >= entry.costhLow) {
                    if (entry.energyLow > bestEnergy || (entry.energyLow == bestEnergy && entry.costhLow > bestCosth)) {
                        bestEnergy = entry.energyLow;
                        bestCosth = entry.costhLow;
                        binIndex = j;
                    }
                }
            }
            double desiredCount = 0.0;
            if (binIndex != -1) {
                const NormEntry& bestEntry = gNormTable[binIndex];
                const genie::GHepParticle* neutrino = mcrec->event->Probe();
                if (neutrino != nullptr) {
                    nu_pdg = neutrino->Pdg();
                }
                const genie::Interaction* interaction = mcrec->event->Summary();
                if (interaction->ProcInfo().IsWeakCC()) {
                    iscc = 1; isnc = 0;
                } else if (interaction->ProcInfo().IsWeakNC()) {
                    iscc = 0; isnc = 1;
                }
                if (isnc == 1) {
                    desiredCount = bestEntry.nc;
                } else if (std::abs(nu_pdg) == 12) {
                    desiredCount = (nu_pdg > 0) ? bestEntry.nue : bestEntry.nuebar;
                } else if (std::abs(nu_pdg) == 14) {
                    desiredCount = (nu_pdg > 0) ? bestEntry.numu : bestEntry.numubar;
                }
            }
            double p = 1.0;
            if (binIndex != -1 && binCounts[binIndex] > 0) {
                p = desiredCount / binCounts[binIndex];
                if (p > 1.0) p = 1.0;
            }
            if (randGen.Uniform() < p) {
                double E_reco = E_true;
                costh_reco = costh_true;
                LoE_true = 12742 / E_true;
                Q2 = mcrec->event->Summary()->Kine().Q2(true);
                W = mcrec->event->Summary()->Kine().W(true);
                
                // Apply the reaction mode determination:
                const auto& event = *(mcrec->event);
                valor::conv::ReacMode_t reac_mode = valor::util::DetermineReactionMode(&event, false);
                reac_mod = reac_mode;
                
                // Update reaction mode count
                modeCounts[reac_mod]++;
                
                bool is_nue = (std::abs(nu_pdg) == 12);
                bool is_antineutrino = (nu_pdg < 0);
                if (isnc == 1) {
                    sample_branch = valor::conv::ksJFCAtmoNC;       // JUNO Full contained Atmospheric NC (53)
                } else if (std::abs(nu_pdg) == 12) {
                    if (nu_pdg > 0)
                        sample_branch = valor::conv::ksJFCAtmoNuECC;    // JUNO Full contained Atmo nue CC (51)
                    else
                        sample_branch = valor::conv::ksJFCAtmoNuEBarCC; // JUNO Full contained Atmo nuebar CC (52)
                } else if (std::abs(nu_pdg) == 14) {
                    if (nu_pdg > 0)
                        sample_branch = valor::conv::ksJFCAtmoNuMuCC;   // JUNO Full contained Atmo numu CC (49)
                    else
                        sample_branch = valor::conv::ksJFCAtmoNuMuBarCC; // JUNO Full contained Atmo numubar CC (50)
                } else {
                    sample_branch = -1; 
                }
                beam = valor::conv::kbAtmo_Jiangmen;
                wght = p;
                outputTree->Fill();
                savedEvents++;
            }
        }
    }

    // Write out the tree and close the file
    outputTree->Write();
    outputFile->Close();
    //delete ghep_tree;
    delete outputFile;
    VLOG("valor_juno", pNOTICE) << "Done!";
    std::cout << "Total saved events: " << savedEvents << std::endl;

    // Output the event counts for each reaction mode
    std::cout << "Reaction Mode Event Counts:" << std::endl;
    for (const auto& kv : modeCounts) {
        std::cout << "Mode " << kv.first << ": " << kv.second << " events." << std::endl;
    }
    return 0;
}



int main(int argc, char** argv) {
  /*
    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " <normalization txt file> <input files> <output file name>" << std::endl;
        return 1;
    }
    std::string normFile = argv[1];
    if (!readNormTable(normFile)) {
        return 1;
    }
    int inputStartIndex = 2;
    std::vector<std::string> inputFiles;
    for (int i = inputStartIndex; i < argc - 1; ++i) {
        inputFiles.emplace_back(argv[i]);
    }
    std::string outputFileName = argv[argc - 1];
    return processEvents(inputFiles, outputFileName);
  */

  CmdLnArgList arg;
  arg.Read(argc,argv, "$VALOR_EXPT/config/apps/make_dst.arg");
  VLOG("valor_juno", pNOTICE) << arg;
  
  string arg_input_normalisation = arg.GetString("input-normalisation");
  string arg_input_directory     = arg.GetString("input-directory");
  string arg_output              = arg.GetString("output");
  string arg_sample              = arg.GetString("sample");

  std::string thresholds = util::CompleteConfigPath(arg.GetString("mesg-thresholds"));
  valor::Messenger::Instance()->Configure(thresholds);

  if (!readNormTable(arg_input_normalisation)) {
    return 1;
  }

  std::vector<string> inputFiles;
  for( const auto& entry : std::filesystem::directory_iterator(arg_input_directory) ) {
    if( entry.path().extension() == ".root" ) inputFiles.emplace_back( entry.path() );
  }

  return processEvents(inputFiles, arg_output);
}

