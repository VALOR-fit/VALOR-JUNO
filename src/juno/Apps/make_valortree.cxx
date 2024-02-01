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
#include "valor/GENIEUtilities/GENIEIncludes.h"
#include "valor/GENIEUtilities/UtlEventRec.h"

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"




int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <input files>" << std::endl;
    return 1;
  }
  std::vector<std::string> inputFiles;
  for (int i = 1; i < argc; ++i) {
    inputFiles.emplace_back(argv[i]);
  }

  ROOT::EnableImplicitMT();

  ROOT::RDataFrame df("gtree", inputFiles);
  auto f = df.Define("Q2",
                     [](const genie::NtpMCEventRecord &mcrec) {
                       return mcrec.event->Summary()->Kine().Q2(true);
                     },
                     {"gmcrec"})
               .Define("W",
                       [](const genie::NtpMCEventRecord &mcrec) {
                         return mcrec.event->Summary()->Kine().W(true);
                       },
                       {"gmcrec"})
               .Snapshot("mytree", "test.root", {"Q2", "W"});
  return 0;
}
