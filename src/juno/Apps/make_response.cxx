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

int processEvents(const std::vector<std::string>& inputFiles, const std::string& outputFileName) {
    time_t scriptstartrawtime;
    struct tm * scriptstarttimeinfo;
    time (&scriptstartrawtime);
    scriptstarttimeinfo = localtime (&scriptstartrawtime);
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << " Script start local time and date:  " << asctime(scriptstarttimeinfo);
    std::cout << "-----------------------------------------------------------" << std::endl;
    // Create a chain to manage input files and branches
    TChain ghep_tree("juno_true_tree");
    for (const auto& file : inputFiles) {
        ghep_tree.Add(file.c_str());
    }

    // Variables for branch data
    double E_true, costh_true, wght;
    int nu_pdg, isnc, iscc, nkaon, npip, npim, npi0, reac_mod, event_id, sample_branch, beam;
    double Q2, W;

    // Set branch addresses
    ghep_tree.SetBranchAddress("enu_true", &E_true);
    ghep_tree.SetBranchAddress("nu_pdg", &nu_pdg);
    ghep_tree.SetBranchAddress("costh_true", &costh_true);
    ghep_tree.SetBranchAddress("Q2", &Q2);
    ghep_tree.SetBranchAddress("W", &W);
    ghep_tree.SetBranchAddress("nCC", &iscc);
    ghep_tree.SetBranchAddress("nNC", &isnc);
    ghep_tree.SetBranchAddress("nKaon", &nkaon);
    ghep_tree.SetBranchAddress("nPionp", &npip);
    ghep_tree.SetBranchAddress("nPionm", &npim);
    ghep_tree.SetBranchAddress("nPion0", &npi0);
    ghep_tree.SetBranchAddress("reac_mode", &reac_mod);
    ghep_tree.SetBranchAddress("event_id", &event_id);
    ghep_tree.SetBranchAddress("wght", &wght);
    ghep_tree.SetBranchAddress("sample", &sample_branch);
    ghep_tree.SetBranchAddress("beam", &beam);

    // Prepare output file and tree
    TFile outputFile(outputFileName.c_str(), "RECREATE");
    TTree outputTree("juno_reco_dst", "Processed Data");

    // Create new branches and copy old branches to the output tree
    double E_reco, costh_reco;
    outputTree.Branch("enu_reco", &E_reco);
    outputTree.Branch("costh_reco", &costh_reco);
    outputTree.Branch("Q2", &Q2);
    outputTree.Branch("W", &W);
    outputTree.Branch("enu_true", &E_true);
    outputTree.Branch("costh_true", &costh_true);
    outputTree.Branch("nCC", &iscc);
    outputTree.Branch("nNC", &isnc);
    outputTree.Branch("nKaon", &nkaon);
    outputTree.Branch("nPionp", &npip);
    outputTree.Branch("nPionm", &npim);
    outputTree.Branch("nPion0", &npi0);
    outputTree.Branch("reac_mode", &reac_mod);
    outputTree.Branch("nu_pdg", &nu_pdg);
    outputTree.Branch("event_id", &event_id);
    outputTree.Branch("wght", &wght);
    outputTree.Branch("sample", &sample_branch);
    outputTree.Branch("beam", &beam);

    // Random number generator
    TRandom3 randGen;

    Long64_t nentries = ghep_tree.GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        ghep_tree.GetEntry(i);

        double ratio = 1.0, efficiency = 0.0;
        if (nu_pdg == 14 || nu_pdg == -14) {
            ratio = (E_true <= 3.5) ? (nu_pdg == 14 ? exp(-0.02956 * pow(E_true, 2) - 0.005532 * E_true - 0.004352) : exp(-0.06302 * pow(E_true, 2) + 0.06550 * E_true - 0.03779))
                                    : (nu_pdg == 14 ? 0.747 * exp(-0.0351 * E_true) : 0.674 * exp(-0.0633 * E_true));
        } else if (nu_pdg == 12 || nu_pdg == -12) {
            ratio = -0.076 * log(E_true) + 0.9561;
        }
        efficiency = ((nu_pdg == 14 ? 0.70 : (nu_pdg == -14 ? 0.60 : (nu_pdg == 12 ? 0.80 : (nu_pdg == -12 ? 0.55 : 0.0)))) * ratio);

        if (randGen.Uniform() < efficiency) {
            double E_sigma = 0.3 / sqrt(E_true);
            E_reco = randGen.Gaus(E_true, E_sigma);

            double cos_sigma = 1 / 2.35482; // Assuming sigma for Gaussian distribution of cos(theta)
            costh_reco = randGen.Gaus(costh_true, cos_sigma);

            // Ensure E_reco is non-negative
            if (E_reco < 0) {
                E_reco = 0;
            }

            // Ensure costh_reco is within [-1, 1]
            if (costh_reco < -1) {
                costh_reco = -1;
            } else if (costh_reco > 1) {
                costh_reco = 1;
            }

            outputTree.Fill();  // Fill the tree with all branches only for efficient events
        }
    }

    outputTree.Write();
    outputFile.Close();

    std::cout << "Total events processed: " << nentries << std::endl;
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input files> <output file name>" << std::endl;
        return 1;
    }

    std::vector<std::string> inputFiles;
    for (int i = 1; i < argc - 1; ++i) {
        inputFiles.push_back(argv[i]);
    }

    std::string outputFileName = argv[argc - 1];

    return processEvents(inputFiles, outputFileName);
}
