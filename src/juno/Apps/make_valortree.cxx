#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>

// System
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
// ROOT
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TSystem.h>
// VALOR/SDK
#include "valor/Conventions/Constants.h"
#include "valor/Conventions/ReacMode.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Utilities/Env.h"
#include "valor/CmdLine/CmdLnArgList.h"
#include "valor/Utilities/UtlNumStat.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/Rnd.h"

// GENIE
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/RunOpt.h"

using std::string;
using namespace genie;

int processEvents(const std::vector<std::string>& inputFiles, const std::string& outputFileName) {
  TChain* ghep_tree = new TChain("gtree");
  ghep_tree->SetMaxVirtualSize(-1);
  for (const auto& file : inputFiles) {
    ghep_tree->Add(file.c_str());
  }
  
  genie::NtpMCEventRecord *mcrec = nullptr;
  ghep_tree->SetBranchAddress("gmcrec", &mcrec);
  
  TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
  TTree* outputTree = new TTree("valor_tree", "Processed Data");
  double Q2, W, E_nu, cos_theta;
  outputTree->Branch("Q2", &Q2);
  outputTree->Branch("W", &W);
  outputTree->Branch("E_nu", &E_nu);
  outputTree->Branch("cos_theta", &cos_theta);
  
 Long64_t nentries = ghep_tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    ghep_tree->GetEntry(i);
    
    if (mcrec && mcrec->event) {
      const genie::EventRecord& event = *(mcrec->event);
      const TLorentzVector& nuP4 = *(event.Probe()->P4()); 
      E_nu = nuP4.Energy();
      if (E_nu < 3.0 || E_nu > 3.4) {
        continue;}
      cos_theta = -nuP4.Pz()/E_nu;       
      Q2 = mcrec->event->Summary()->Kine().Q2(true);
      W = mcrec->event->Summary()->Kine().W(true);

      outputTree->Fill();
    }
  }
  
  outputTree->Write();
  outputFile->Close();
  
  delete ghep_tree;
  delete outputFile; // Also deletes outputTree
  
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 3) { // Change here to require at least one input file and an output file name
    std::cout << "Usage: " << argv[0] << " <input files> <output file name>" << std::endl;
    return 1;
  }
  std::vector<std::string> inputFiles;
  for (int i = 1; i < argc - 1; ++i) { // Adjust here to reserve the last argument for the output file name
    inputFiles.emplace_back(argv[i]);
  }
  std::string outputFileName = argv[argc - 1]; // Get the last argument as the output file name
  return processEvents(inputFiles, outputFileName);
}
