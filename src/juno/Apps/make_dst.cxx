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
            return valor::conv::ksNuEBarCChigh; // Assuming high resolution, adjust as necessary
        } else {
            if (npi == 0 && npi0 == 0) return valor::conv::ksNuECC;
            if (npi == 1 && npi0 == 0) return valor::conv::ksNuECC1piC;
            if (npi == 0 && npi0 == 1) return valor::conv::ksNuECC1pi0;
            return valor::conv::ksNuECCoth;
        }
    } else {
        if (is_antineutrino) {
            return valor::conv::ksNuMuBarCChigh; // Assuming high resolution, adjust as necessary
        } else {
            if (npi == 0 && npi0 == 0) return valor::conv::ksNuMuCC;
            if (npi == 1 && npi0 == 0) return valor::conv::ksNuMuCC1piC;
            if (npi == 0 && npi0 == 1) return valor::conv::ksNuMuCC1pi0;
            if (npi == 1 && npi0 == 1) return valor::conv::ksNuMuCC1piC1pi0;
            return valor::conv::ksNuMuCCoth;
        }
    }
}

int processEvents(const std::vector<std::string>& inputFiles, const std::string& outputFileName) {
    TChain* ghep_tree = new TChain("gtree");
    static const int gPrintEveryNEvents = 50000;
    ghep_tree->SetMaxVirtualSize(-1);
    for (const auto& file : inputFiles) {
        ghep_tree->Add(file.c_str());
    }

    genie::NtpMCEventRecord *mcrec = nullptr;
    ghep_tree->SetBranchAddress("gmcrec", &mcrec);

    TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    TTree* outputTree = new TTree("juno_true_tree", "Processed Data");

    double Q2, W, E_true, costh_true, costh_reco, wght, LoE_true;
    int isnc = 0, iscc = 0, nkaon = 0, npip = 0, npim = 0, npi0 = 0,reac_mod_wo_PID = -1, reac_mod = -1, nu_pdg = -1, event_id = 0, sample_branch = -1, beam = -1;
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
    outputTree->Branch("nu_PDG", &nu_pdg);
    outputTree->Branch("event_id", &event_id);
    outputTree->Branch("wght", &wght);
    outputTree->Branch("sample", &sample_branch);
    outputTree->Branch("beam", &beam);
    
    // Confusion matrix for efficiency
    double efficiencyMatrix[5][5] = {
        {0.53, 0.29, 0.04, 0.13, 0.01},
        {0.08, 0.71, 0.15, 0.07, 0.00},
        {0.08, 0.38, 0.54, 0.01, 0.00},
        {0.06, 0.07, 0.00, 0.74, 0.13},
        {0.06, 0.03, 0.00, 0.34, 0.57}
    };

    int confusionMatrix[5][5] = {0};

    TRandom3 randGen;


    auto flavorToIndex = [](int nu_pdg) {
        switch (nu_pdg) {
            case 12:  return 1; // nue
            case -12: return 2; // nuebar
            case 14:  return 3; // numu
            case -14: return 4; // numubar
            default:  return 0;  // nc
        }
    };


    Long64_t nentries = ghep_tree->GetEntries();

    

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
            iscc = 0; isnc = 0; npip = 0; npim = 0; npi0 = 0; nkaon = 0; reac_mod = 0; event_id = 0; reac_mod_wo_PID = 0;

            if (neutrino != nullptr) {
                nu_pdg = neutrino->Pdg();
            }
            if (interaction->ProcInfo().IsWeakCC()) {
                iscc = 1;
                isnc = 0; 
            } else if (interaction->ProcInfo().IsWeakNC()) {
                isnc = 1;
                iscc = 0; 
            }
            switch(nu_pdg){
                case 12: beam = valor::conv::kbAtmo_JUNONuE;
                case -12: beam = valor::conv::kbAtmo_JUNONuEBar;
                case 14: beam = valor::conv::kbAtmo_JUNONuMu;
                case -14: beam = valor::conv::kbAtmo_JUNONuMuBar;
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
            reac_mod_wo_PID = reac_mode;

            // Determine final flavor using efficiency matrix
            int initialIndex = flavorToIndex(nu_pdg);
            double randomValue = randGen.Uniform();
            double cumulativeProbability = 0.0;
            int finalIndex = 0; // Default to nc

            for (int j = 0; j < 5; ++j) {
                cumulativeProbability += efficiencyMatrix[initialIndex][j];
                if (randomValue < cumulativeProbability) {
                    finalIndex = j;
                    break;
                }
            }

            // Map matrix index back to PDG code
            switch (finalIndex) {
                case 1: nu_pdg = 12; break;
                case 2: nu_pdg = -12; break;
                case 3: nu_pdg = 14; break;
                case 4: nu_pdg = -14; break;
                default: nu_pdg = 0; break;
            }

            // Update confusion matrix
            confusionMatrix[initialIndex][finalIndex]++;

            // If a CC event is misidentified as an NC event
            if (iscc == 1 && finalIndex == 0) {
                iscc = 0;
                isnc = 1;
            }
            valor::conv::ReacMode_t reac_mode_PID = valor::util::DetermineReactionModeRead(nu_pdg, iscc, isnc,&event, false);
            reac_mod = reac_mode_PID;

            E_true = nuP4.Energy();
            LoE_true = 12742 / E_true; // L in km
            VLOG("valor_juno", pDEBUG) << "Read event: "
                                       << "Etrue = " << E_true << " GeV; ";

            costh_true = -nuP4.Pz() / E_true;
            Q2 = mcrec->event->Summary()->Kine().Q2(true);
            W = mcrec->event->Summary()->Kine().W(true);

            bool is_nue = (nu_pdg == 12 || nu_pdg == -12);
            bool is_antineutrino = (nu_pdg < 0);
            valor::conv::Sample_t sample = InferSampleType(npip + npim, npi0, is_nue, is_antineutrino);
            sample_branch = sample;

            wght = 1.0;
            outputTree->Fill();
            savedEvents++;
        }
    }

    outputTree->Write();
    outputFile->Close();

    delete ghep_tree;
    delete outputFile;

    std::cout << "Confusion Matrix:\n";
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            std::cout << confusionMatrix[i][j] << " ";
        }
        std::cout << "\n";
    }
    VLOG("valor_juno", pNOTICE) << "Done!";
    std::cout << "Total saved events: " << savedEvents << std::endl; 

    return 0;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << "<input files> <output file name>" << std::endl;
        return 1;
    }

    int inputStartIndex = 1;

    std::vector<std::string> inputFiles;
    for (int i = inputStartIndex; i < argc - 1; ++i) {
        inputFiles.emplace_back(argv[i]);
    }
    std::string outputFileName = argv[argc - 1];
    return processEvents(inputFiles, outputFileName);
}



