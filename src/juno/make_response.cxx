#include <cmath>
#include <vector>
#include <string>
#include <TRandom3.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include "valor/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"

int processEvents(float efficiency, const std::vector<std::string>& inputFiles, const std::string& outputFileName) {
    TRandom3 randomGen(0);  // Random number generator

    TChain* ghep_tree = new TChain("gtree");  // GHEP event chain
    static const int gPrintEveryNEvents = 500000;
    ghep_tree->SetMaxVirtualSize(-1);
    for (const auto& file : inputFiles) {
        ghep_tree->Add(file.c_str());
    }

    int discarded = 0;
    genie::NtpMCEventRecord *mcrec = nullptr;
    ghep_tree->SetBranchAddress("gmcrec", &mcrec);

    TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    TTree* outputTree = new TTree("valor_tree", "Processed Data");
    double E_true, E_reco, cos_true, cos_reco, theta_true, theta_reco;
    double lepton_px, lepton_py, lepton_pz;  // Lepton momentum components
    outputTree->Branch("E_true", &E_true);
    outputTree->Branch("E_reco", &E_reco);
    outputTree->Branch("costh_true", &cos_true);
    outputTree->Branch("costh_reco", &cos_reco);
    outputTree->Branch("lepton_px", &lepton_px);
    outputTree->Branch("lepton_py", &lepton_py);
    outputTree->Branch("lepton_pz", &lepton_pz);

    Long64_t nentries = ghep_tree->GetEntries();
    Long64_t processed = 0;

    for (Long64_t i = 0; i < nentries; ++i) {
        ghep_tree->GetEntry(i);
        if (gPrintEveryNEvents && !(i % gPrintEveryNEvents)) {
            VLOG("valor_juno", pNOTICE) << "** Reading event: " << i;
        }

        float randValue = randomGen.Uniform();
        if (randValue > efficiency) {
            discarded++;
            continue; 
        }

        processed++;

        if (mcrec && mcrec->event) {
    const genie::EventRecord& event = *(mcrec->event);
    const TLorentzVector& nuP4 = *(event.Probe()->P4());
    E_true = nuP4.Energy();
    cos_true = -nuP4.Pz() / E_true;

    genie::GHepParticle *particle = nullptr;
    TLorentzVector *leptonP4 = nullptr;

    // Use MakeIterator() if that's what suggested in your GENIE version
    TIter iter = event.MakeIterator();
    while ((particle = dynamic_cast<genie::GHepParticle *>(iter.Next()))) {
        if (particle->Status() == genie::kIStStableFinalState && 
            (abs(particle->Pdg()) == 11 || abs(particle->Pdg()) == 13)) {  // Electron or muon
            leptonP4 = particle->P4();
            break;
        }
    }


    if (leptonP4) {
        lepton_px = leptonP4->Px();
        lepton_py = leptonP4->Py();
        lepton_pz = leptonP4->Pz();

        // Calculate cos_reco based on the angle between neutrino and lepton momentum vectors
        double lepton_mag = sqrt(lepton_px*lepton_px + lepton_py*lepton_py + lepton_pz*lepton_pz);
        cos_reco = -lepton_pz / lepton_mag;
    }

    outputTree->Fill();
        }
    }

    outputTree->Write();
    outputFile->Close();

    delete ghep_tree;
    delete outputFile;

    VLOG("valor_juno", pNOTICE) << "Done!";
    VLOG("valor_juno", pNOTICE) << "Total events: " << nentries;
    VLOG("valor_juno", pNOTICE) << "Processed events: " << processed;
    VLOG("valor_juno", pNOTICE) << "Discarded events: " << discarded;

    return 0;
}

int main(int argc, char **argv) {
    if (argc < 4) { 
        std::cout << "valor-juno_make_response Usage: " << argv[0] << " <efficiency> <input files> <output file name>" << std::endl;
        return 1;
    }
    
    float efficiency = std::stof(argv[1]);
    if (efficiency < 0.0 || efficiency > 1.0) {
        std::cerr << "Efficiency must be between 0 and 1." << std::endl;
        return 1;
    }

    std::vector<std::string> inputFiles;
    for (int i = 2; i < argc - 1; ++i) { 
        inputFiles.emplace_back(argv[i]);
    }

    std::string outputFileName = argv[argc - 1];

    return processEvents(efficiency, inputFiles, outputFileName);
}
