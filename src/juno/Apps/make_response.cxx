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

int processEvents(float efficiency, const std::vector<std::string>& inputFiles, const std::string& outputFileName) {
  TRandom3 randomGen(0); 

  TChain* ghep_tree = new TChain("gtree");
  static const int gPrintEveryNEvents = 500000;
  ghep_tree->SetMaxVirtualSize(-1);
  for (const auto& file : inputFiles) {
    ghep_tree->Add(file.c_str());
  }

  int discarded = 0; // Fixed missing semicolon
  genie::NtpMCEventRecord *mcrec = nullptr;
  ghep_tree->SetBranchAddress("gmcrec", &mcrec);

  TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
  TTree* outputTree = new TTree("valor_tree", "Processed Data");
  double E_true, E_reco, cos_true, cos_reco;
  outputTree->Branch("E_true", &E_true);
  outputTree->Branch("E_reco", &E_reco);
  outputTree->Branch("cos_theta_true", &cos_true);
  outputTree->Branch("cos_theta_reco", &cos_reco);

  Long64_t nentries = ghep_tree->GetEntries();
  Long64_t processed = 0; // To track processed events

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

    processed++; // Increment for each processed event

    if (mcrec && mcrec->event) {
      const genie::EventRecord& event = *(mcrec->event);
      const TLorentzVector& nuP4 = *(event.Probe()->P4());
      E_true = nuP4.Energy();
      cos_true = -nuP4.Pz() / E_true;

      double E_sigma = 3 / sqrt(E_true);
      E_reco = randomGen.Gaus(E_true, E_sigma);
      
      double cos_sigma = 1 / 2.35482; // Assuming sigma for Gaussian distribution of cos(theta)
      cos_reco = randomGen.Gaus(cos_true, cos_sigma); // Fixed missing semicolon
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
