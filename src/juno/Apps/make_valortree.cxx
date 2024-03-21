// System
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
// ROOT
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TSystem.h>
// VALOR-SDK
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

using std::string;
using namespace genie;

typedef enum NuMuCCExclusiveTopology{

  kUnknownTopology              = -1,
  kCC0Pi                        =  0,
  kCC1Pi                        =  1,
  kCC2Pi                        =  2,
  kCCPi0                        =  3,
  kCCOther                      =  4

} ExclusiveTopology_t;

int processEvents(const std::vector<std::string>& inputFiles, const std::string& outputFileName) {
  TChain* ghep_tree = new TChain("gtree");
  static const int gPrintEveryNEvents   = 500000;
  ghep_tree->SetMaxVirtualSize(-1);
  for (const auto& file : inputFiles) {
    ghep_tree->Add(file.c_str());
  }
  
  genie::NtpMCEventRecord *mcrec = nullptr;
  ghep_tree->SetBranchAddress("gmcrec", &mcrec);
  
  TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
  TTree* outputTree = new TTree("valor_tree", "Processed Data");
  double Q2, W, E_true, E_reco, cos_theta;
  int isnc=0;int iscc=0;int nkaon=0;int npip=0;int npim=0;int npi0=0;int reac_mode=-1;int nu_pdg=-1;
  outputTree->Branch("Q2", &Q2);
  outputTree->Branch("W", &W);
  outputTree->Branch("E_true", &E_true);
  outputTree->Branch("E_reco", &E_reco);
  outputTree->Branch("cos_theta", &cos_theta);
  outputTree->Branch("nCC", &isnc);
  outputTree->Branch("nNC", &iscc);
  outputTree->Branch("nKaon", &nkaon);
  outputTree->Branch("nPionp", &npip);
  outputTree->Branch("nPionm", &npim);
  outputTree->Branch("nPion0", &npi0);
  outputTree->Branch("reac_mode", &reac_mode);
  outputTree->Branch("nu_PDG", &nu_pdg);

  Long64_t nentries = ghep_tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    ghep_tree->GetEntry(i);
    if (gPrintEveryNEvents && !(i%gPrintEveryNEvents)) {
      VLOG("valor_juno", pNOTICE) << "** Reading event: " << i;
    }
    else {
      VLOG("valor_juno", pDEBUG) << "** Reading event: " << i;
    }
    
    if (mcrec && mcrec->event) {
      const genie::EventRecord& event = *(mcrec->event);
      const genie::Interaction* interaction = event.Summary();
      const TLorentzVector& nuP4 = *(event.Probe()->P4()); 
      const genie::GHepParticle *neutrino = event.Probe();

      if (neutrino != nullptr) {
        nu_pdg = neutrino->Pdg();}

      iscc = 0;  isnc = 0;  npip = 0;  npim = 0;  npi0 = 0;  nkaon = 0;
      if (interaction->ProcInfo().IsWeakCC()) {
        iscc = 1;
      } else if (interaction->ProcInfo().IsWeakNC()) {
        isnc = 1;
      }
      TIter iter(event.MakeIterator()); 
      genie::GHepParticle* p = nullptr;
      while ((p = dynamic_cast<genie::GHepParticle*>(iter.Next()))) {
        if (p->Status() == genie::kIStStableFinalState) {
          int pdg = p->Pdg();
          E_reco += p->Energy();
          if (pdg == genie::kPdgPiP) npip++;
          else if (pdg == genie::kPdgPiM) npim++;
          else if (pdg == genie::kPdgPi0) npi0++;
          if (pdg == genie::kPdgKP || pdg == genie::kPdgKM || pdg == genie::kPdgK0 || pdg == genie::kPdgAntiK0) {
            nkaon++;
          }
        }
      }
      valor::conv::ReacType_t reac = valor::conv::kmtUndefined; 
    if (iscc) { 
        if      (nkaon )                      {reac = valor::conv::kmtCCoth;}
        else if (npip+npim == 1 && npi0 == 0)  {reac = valor::conv::kmtCC1piC;}
        else if (npip+npim == 0 && npi0 == 1)  {reac = valor::conv::kmtCC1pi0;}
        else if (npip+npim == 2 && npi0 == 0)  {reac = valor::conv::kmtCC2piC;}
        else if (npip+npim == 0 && npi0 == 2)  {reac = valor::conv::kmtCC2pi0;}
        else if (npip+npim == 1 && npi0 == 1)  {reac = valor::conv::kmtCCpiCpi0;}
        else                                   {reac = valor::conv::kmtCCoth;}
    }
    else if (isnc) { // start counting pions
        if      (nkaon)                       {reac = valor::conv::kmtNCoth;}
        else if (npip+npim == 1 && npi0 == 0)  {reac = valor::conv::kmtNC1piC;}
        else if (npip+npim == 0 && npi0 == 1)  {reac = valor::conv::kmtNC1pi0;}
        else if (npip+npim == 2 && npi0 == 0)  {reac = valor::conv::kmtNC2piC;}
        else if (npip+npim == 0 && npi0 == 2)  {reac = valor::conv::kmtNC2pi0;}
        else if (npip+npim == 1 && npi0 == 1)  {reac = valor::conv::kmtNCpiCpi0;}
        else                                   {reac = valor::conv::kmtNCoth;}
    }
    else {
        reac = valor::conv::kmtOther;
    }

    valor::conv::NuFlavour_t final_flav = valor::util::FlavourFromPDG(nu_pdg);
    valor::conv::NuFlavour_t init_flav = final_flav;

    // Get the real valor mode
    valor::conv::ReacMode_t  reac_mode  = valor::util::ReacTypeToMode(reac,init_flav,final_flav);

    E_true = nuP4.Energy();
    VLOG("valor_juno",pDEBUG) << "Read event: "
      << "Etrue         = "<< E_true  << " GeV; ";
  
      /*if (E_true < 3.0 || E_true > 3.4) {
        continue;}*/
      
    cos_theta = -nuP4.Pz()/E_true;       
    Q2 = mcrec->event->Summary()->Kine().Q2(true);
    W = mcrec->event->Summary()->Kine().W(true);

    outputTree->Fill();
    }
  }
  
  outputTree->Write();
  outputFile->Close();
  
  delete ghep_tree;
  delete outputFile; // Also deletes outputTree
  VLOG("valor_juno", pNOTICE)  << "Done!";
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
