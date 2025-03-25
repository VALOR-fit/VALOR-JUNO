#include <vector>
#include <string>
#include <iostream>
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

bool readNormTable(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error opening normalization file: " << filename << std::endl;
        return false;
    }
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
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

struct EventBranchData {
    int nu_pdg;
    bool is_CC;
    bool is_NC;
    bool is_charm;
    bool is_strange;
    bool is_QE;
    bool is_MEC;
    bool is_NuEEl;
    bool is_IMD;
    bool is_IMDanh;
    bool is_Coherent;
    std::vector<int> particlePDG;
};

EventBranchData ExtractEventBranchData(const genie::EventRecord* event, bool use_final_state) {
    EventBranchData data;
    const genie::GHepParticle* neutrino = event->Probe();
    data.nu_pdg = neutrino ? neutrino->Pdg() : -1;
    
    const genie::Interaction* interaction = event->Summary();
    if (interaction->ProcInfo().IsWeakCC()) {
        data.is_CC = true;
        data.is_NC = false;
    } else if (interaction->ProcInfo().IsWeakNC()) {
        data.is_CC = false;
        data.is_NC = true;
    } else {
        data.is_CC = data.is_NC = false;
    }
    
    data.is_charm   = util::EventContainsQuarkFlavour(event, genie::kPdgCQuark);
    data.is_strange = util::EventContainsQuarkFlavour(event, genie::kPdgSQuark);
    data.is_QE      = interaction->ProcInfo().IsQuasiElastic() && !data.is_strange && !data.is_charm;
    data.is_MEC     = interaction->ProcInfo().IsMEC();
    data.is_NuEEl   = interaction->ProcInfo().IsNuElectronElastic();
    data.is_IMD     = interaction->ProcInfo().IsInverseMuDecay();
    data.is_IMDanh  = interaction->ProcInfo().IsIMDAnnihilation();
#if __GENIE_RELEASE_CODE__ >= GRELCODE(3, 4, 0)
    data.is_Coherent = interaction->ProcInfo().IsCoherentProduction();
#else
    data.is_Coherent = interaction->ProcInfo().IsCoherent();
#endif

    std::vector<int> primary_hadr_system, final_hadr_system;
    valor::util::GetHadronicSystems(event, primary_hadr_system, final_hadr_system);
    const std::vector<int>& hadr_system = (use_final_state ? final_hadr_system : primary_hadr_system);
    for (size_t i = 0; i < hadr_system.size(); ++i) {
        const genie::GHepParticle* p = event->Particle(hadr_system[i]);
        data.particlePDG.push_back(p->Pdg());
    }
    
    return data;
}

// -------------------- processEvents --------------------
int processEvents(const std::vector<std::string>& inputFiles, const std::string& outputFileName) {
    TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    outputFile->cd();
    TTree* outputTree = new TTree("juno_true_tree", "Processed Data");
    TChain* ghep_tree = new TChain("gtree");
    static const int gPrintEveryNEvents = 200000;
    ghep_tree->SetMaxVirtualSize(-1);
    for (const auto& file : inputFiles) {
        ghep_tree->Add(file.c_str());
    }
    genie::NtpMCEventRecord* mcrec = nullptr;
    ghep_tree->SetBranchAddress("gmcrec", &mcrec);
    Long64_t nentries = ghep_tree->GetEntries();

    // First pass: count events per bin and per category
    // Categories: "nc", "nue", "nuebar", "numu", "numubar"
    std::vector<std::unordered_map<std::string, int>> binCategoryCounts(gNormTable.size());
    for (Long64_t i = 0; i < nentries; ++i) {
        ghep_tree->GetEntry(i);
        if (mcrec && mcrec->event) {
            const TLorentzVector& nuP4 = *(mcrec->event->Probe()->P4());
            double E_true_tmp = nuP4.Energy();
            double costh_true_tmp = -nuP4.Pz() / E_true_tmp;            
            if (E_true_tmp > 20.0) continue;
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
                const genie::GHepParticle* neutrino = mcrec->event->Probe();
                int nu_pdg = (neutrino ? neutrino->Pdg() : -1);
                const genie::Interaction* interaction = mcrec->event->Summary();
                std::string category;
                if (interaction->ProcInfo().IsWeakNC()) {
                    category = "nc";
                } else if (std::abs(nu_pdg) == 12) {
                    category = (nu_pdg > 0) ? "nue" : "nuebar";
                } else if (std::abs(nu_pdg) == 14) {
                    category = (nu_pdg > 0) ? "numu" : "numubar";
                } else {
                    category = "unknown";
                }
                if (category != "unknown") {
                    binCategoryCounts[binIndex][category]++;
                }
            }
        }
    }

    // Branch variable declarations
    double Q2, W, E_true, E_reco, costh_true, costh_reco, wght, LoE_true;
    int isnc = 0, iscc = 0, nkaon = 0, npip = 0, npim = 0, npi0 = 0;
    int reac_mod_wo_PID = -1;
    int reac_mod = -1, nu_pdg = -1, event_id = 0, sample_branch = -1, beam = -1;
    outputTree->Branch("Q2", &Q2);
    outputTree->Branch("W", &W);
    outputTree->Branch("enu_true", &E_true);
    outputTree->Branch("enu_reco", &E_reco);
    outputTree->Branch("LoverE_true", &LoE_true);
    outputTree->Branch("costh_true", &costh_true);
    outputTree->Branch("costh_reco", &costh_reco);
    outputTree->Branch("nCC", &iscc);
    outputTree->Branch("nNC", &isnc);
    outputTree->Branch("nKaon", &nkaon);
    outputTree->Branch("nPionp", &npip);
    outputTree->Branch("nPionm", &npim);
    outputTree->Branch("nPi0", &npi0);
    outputTree->Branch("reac_mode_wo_PID", &reac_mod_wo_PID);
    outputTree->Branch("reac_mode", &reac_mod);
    outputTree->Branch("nu_pdg", &nu_pdg);
    outputTree->Branch("event_id", &event_id);
    outputTree->Branch("wght", &wght);
    outputTree->Branch("sample", &sample_branch);
    outputTree->Branch("beam", &beam);

    int nP = 0, nN = 0, nPip = 0, nPim = 0;
    int nKp = 0, nKm = 0, nK0 = 0, nEM = 0, nOtherHad = 0, nNucleus = 0;
    outputTree->Branch("nP", &nP);
    outputTree->Branch("nN", &nN);
    outputTree->Branch("nPip", &nPip);
    outputTree->Branch("nPim", &nPim);
    outputTree->Branch("nPi0_count", &npi0);
    outputTree->Branch("nKp", &nKp);
    outputTree->Branch("nKm", &nKm);
    outputTree->Branch("nK0", &nK0);
    outputTree->Branch("nEM", &nEM);
    outputTree->Branch("nOtherHad", &nOtherHad);
    outputTree->Branch("nNucleus", &nNucleus);
    bool isCC_flag = false, isNC_flag = false, isCharm_flag = false, isStrange_flag = false;
    bool isQE_flag = false, isMEC_flag = false, isNuEEl_flag = false;
    bool isIMD_flag = false, isIMDAnnihilation_flag = false, isCoherent_flag = false;
    outputTree->Branch("isCC", &isCC_flag);
    outputTree->Branch("isNC", &isNC_flag);
    outputTree->Branch("isCharm", &isCharm_flag);
    outputTree->Branch("isStrange", &isStrange_flag);
    outputTree->Branch("isQE", &isQE_flag);
    outputTree->Branch("isMEC", &isMEC_flag);
    outputTree->Branch("isNuEEl", &isNuEEl_flag);
    outputTree->Branch("isIMD", &isIMD_flag);
    outputTree->Branch("isIMDAnnihilation", &isIMDAnnihilation_flag);
    outputTree->Branch("isCoherent", &isCoherent_flag);

    std::unordered_map<int, int> modeCounts;

    int savedEvents = 0;
    TRandom3 randGen;

    // Prepare maps to track processed and accepted counts per bin and category
    std::vector<std::unordered_map<std::string, int>> binCategoryProcessed(gNormTable.size());
    std::vector<std::unordered_map<std::string, int>> binCategoryAccepted(gNormTable.size());

    for (Long64_t i = 0; i < nentries; ++i) {
        ghep_tree->GetEntry(i);
        if (gPrintEveryNEvents && !(i % gPrintEveryNEvents)) {
            // Logging can be added here if needed
        }
        if (mcrec && mcrec->event) {
            const TLorentzVector& nuP4 = *(mcrec->event->Probe()->P4());
            E_true = nuP4.Energy();
            if (E_true > 20.0) continue;
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
            if (binIndex != -1) {
                const genie::GHepParticle* neutrino = mcrec->event->Probe();
                nu_pdg = (neutrino ? neutrino->Pdg() : -1);
                const genie::Interaction* interaction = mcrec->event->Summary();
                bool currentIsNC = false;
                if (interaction->ProcInfo().IsWeakCC()) {
                    iscc = 1; isnc = 0;
                    currentIsNC = false;
                } else if (interaction->ProcInfo().IsWeakNC()) {
                    iscc = 0; isnc = 1;
                    currentIsNC = true;
                }
                std::string category;
                if (currentIsNC) {
                    category = "nc";
                } else if (std::abs(nu_pdg) == 12) {
                    category = (nu_pdg > 0) ? "nue" : "nuebar";
                } else if (std::abs(nu_pdg) == 14) {
                    category = (nu_pdg > 0) ? "numu" : "numubar";
                } else {
                    category = "unknown";
                }
                if (category == "unknown") continue;

                const NormEntry& bestEntry = gNormTable[binIndex];
                double desiredCountDouble = 0.0;
                if (category == "nc") {
                    desiredCountDouble = bestEntry.nc;
                } else if (category == "nue") {
                    desiredCountDouble = bestEntry.nue;
                } else if (category == "nuebar") {
                    desiredCountDouble = bestEntry.nuebar;
                } else if (category == "numu") {
                    desiredCountDouble = bestEntry.numu;
                } else if (category == "numubar") {
                    desiredCountDouble = bestEntry.numubar;
                }
                desiredCountDouble *= 100.0;
                int targetCount = static_cast<int>(desiredCountDouble + 0.5);


                int totalCount = binCategoryCounts[binIndex][category];
                int processed = binCategoryProcessed[binIndex][category];
                int accepted = binCategoryAccepted[binIndex][category];
                int remaining = totalCount - processed;
                int needed = targetCount - accepted;
                binCategoryProcessed[binIndex][category] = processed + 1;
                bool select = false;
                if (remaining <= needed) {
                    select = true;
                } else {
                    double prob = static_cast<double>(needed) / remaining;
                    if (randGen.Uniform() < prob) {
                        select = true;
                    }
                }
                if (select) {
                    binCategoryAccepted[binIndex][category] = accepted + 1;
                    E_reco = E_true;
                    costh_reco = costh_true;
                    LoE_true = 12742 / E_true;
                    Q2 = mcrec->event->Summary()->Kine().Q2(true);
                    W  = mcrec->event->Summary()->Kine().W(true);
                    const auto& event = *(mcrec->event);
                    valor::conv::ReacMode_t reac_mode = valor::util::DetermineReactionMode(&event, false);
                    reac_mod = reac_mode;
                    modeCounts[reac_mod]++;
                    EventBranchData branchData = ExtractEventBranchData(mcrec->event, false);
                    nP = 0; nN = 0; nPip = 0; nPim = 0; npi0 = 0;
                    nKp = 0; nKm = 0; nK0 = 0; nEM = 0; nOtherHad = 0; nNucleus = 0;
                    for (size_t j = 0; j < branchData.particlePDG.size(); ++j) {
                        int pdg = branchData.particlePDG[j];
                        if (pdg == genie::kPdgProton || pdg == genie::kPdgAntiProton)
                            nP++;
                        else if (pdg == genie::kPdgNeutron || pdg == genie::kPdgAntiNeutron)
                            nN++;
                        else if (pdg == genie::kPdgPiP)
                            nPip++;
                        else if (pdg == genie::kPdgPiM)
                            nPim++;
                        else if (pdg == genie::kPdgPi0) {
                            npi0++;
                        }
                        else if (pdg == genie::kPdgKP)
                            nKp++;
                        else if (pdg == genie::kPdgKM)
                            nKm++;
                        else if (pdg == genie::kPdgK0 || pdg == genie::kPdgAntiK0)
                            nK0++;
                        else if (pdg == genie::kPdgGamma || pdg == genie::kPdgElectron || pdg == genie::kPdgPositron)
                            nEM++;
                        else if (genie::pdg::IsHadron(pdg))
                            nOtherHad++;
                        else if (genie::pdg::IsIon(pdg))
                            nNucleus++;
                    }
                    isCC_flag         = branchData.is_CC;
                    isNC_flag         = branchData.is_NC;
                    isCharm_flag      = branchData.is_charm;
                    isStrange_flag    = branchData.is_strange;
                    isQE_flag         = branchData.is_QE;
                    isMEC_flag        = branchData.is_MEC;
                    isNuEEl_flag      = branchData.is_NuEEl;
                    isIMD_flag        = branchData.is_IMD;
                    isIMDAnnihilation_flag = branchData.is_IMDanh;
                    isCoherent_flag   = branchData.is_Coherent;
                    bool is_nue = (std::abs(nu_pdg) == 12);
                    bool is_antineutrino = (nu_pdg < 0);
                    if (isNC_flag) {
                        sample_branch = valor::conv::ksJFCAtmoNC;
                    } else if (std::abs(nu_pdg) == 12) {
                        if (nu_pdg > 0)
                            sample_branch = valor::conv::ksJFCAtmoNuECC;
                        else
                            sample_branch = valor::conv::ksJFCAtmoNuEBarCC;
                    } else if (std::abs(nu_pdg) == 14) {
                        if (nu_pdg > 0)
                            sample_branch = valor::conv::ksJFCAtmoNuMuCC;
                        else
                            sample_branch = valor::conv::ksJFCAtmoNuMuBarCC;
                    } else {
                        sample_branch = -1;
                    }
                    beam = valor::conv::kbAtmo_Jiangmen;
                    wght = 1.0;
                    outputTree->Fill();
                    savedEvents++;
                }
            }
        }
    }

    outputTree->Write();
    outputFile->Close();
    delete ghep_tree;
    std::cout << "Total saved events: " << savedEvents << std::endl;

    std::cout << "Reaction Mode Event Counts:" << std::endl;
    for (const auto& kv : modeCounts) {
        std::cout << "Mode " << kv.first << ": " << kv.second << " events." << std::endl;
    }
    return 0;
}

int main(int argc, char** argv) {
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
}
