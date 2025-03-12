#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
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

int processEvents(const std::vector<std::string>& inputFiles, const std::string& outputFileName) {
    time_t scriptstartrawtime;
    struct tm* scriptstarttimeinfo;
    time(&scriptstartrawtime);
    scriptstarttimeinfo = localtime(&scriptstartrawtime);
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << " Script start local time and date:  " << asctime(scriptstarttimeinfo);
    std::cout << "-----------------------------------------------------------" << std::endl;

    TChain ghep_tree("juno_true_tree");
    for (const auto& file : inputFiles) {
        ghep_tree.Add(file.c_str());
    }
    
    double E_true, costh_true, wght;
    int nu_pdg, isnc, iscc, nkaon, npip, npim, npi0, reac_mod, event_id, sample_branch, beam;
    double Q2, W;
    bool isCharm, isStrange, isQE, isMEC, isNuEEl, isIMD, isIMDAnnihilation, isCoherent;
    int nP, nN, nPip, nPim, nKp, nKm, nK0, nEM, nOtherHad, nNucleus;
    
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
    ghep_tree.SetBranchAddress("nPi0", &npi0);
    ghep_tree.SetBranchAddress("reac_mode", &reac_mod);
    ghep_tree.SetBranchAddress("event_id", &event_id);
    ghep_tree.SetBranchAddress("wght", &wght);
    ghep_tree.SetBranchAddress("sample", &sample_branch);
    ghep_tree.SetBranchAddress("beam", &beam);
    ghep_tree.SetBranchAddress("isCharm", &isCharm);
    ghep_tree.SetBranchAddress("isStrange", &isStrange);
    ghep_tree.SetBranchAddress("isQE", &isQE);
    ghep_tree.SetBranchAddress("isMEC", &isMEC);
    ghep_tree.SetBranchAddress("isNuEEl", &isNuEEl);
    ghep_tree.SetBranchAddress("isIMD", &isIMD);
    ghep_tree.SetBranchAddress("isIMDAnnihilation", &isIMDAnnihilation);
    ghep_tree.SetBranchAddress("isCoherent", &isCoherent);
    ghep_tree.SetBranchAddress("nP", &nP);
    ghep_tree.SetBranchAddress("nN", &nN);
    ghep_tree.SetBranchAddress("nPip", &nPip);
    ghep_tree.SetBranchAddress("nPim", &nPim);
    ghep_tree.SetBranchAddress("nKp", &nKp);
    ghep_tree.SetBranchAddress("nKm", &nKm);
    ghep_tree.SetBranchAddress("nK0", &nK0);
    ghep_tree.SetBranchAddress("nEM", &nEM);
    ghep_tree.SetBranchAddress("nOtherHad", &nOtherHad);
    ghep_tree.SetBranchAddress("nNucleus", &nNucleus);

    TFile outputFile(outputFileName.c_str(), "RECREATE");
    TTree outputTree("juno_reco_dst", "Processed Data");
    double E_reco, costh_reco;
    int new_sample;
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
    outputTree.Branch("nPi0", &npi0);
    outputTree.Branch("reac_mode", &reac_mod);
    outputTree.Branch("nu_pdg", &nu_pdg);
    outputTree.Branch("event_id", &event_id);
    outputTree.Branch("wght", &wght);
    outputTree.Branch("sample", &new_sample);
    outputTree.Branch("beam", &beam);
    
    int count_NuECC = 0, count_NuEBarCC = 0, count_NuMuCC = 0, count_NuMuBarCC = 0;
    
    TRandom3 randGen;
    Long64_t nentries = ghep_tree.GetEntries();
    int savedEvents = 0;
    for (Long64_t i = 0; i < nentries; ++i) {
        ghep_tree.GetEntry(i);
        if (i % 200000 == 0)
            std::cout << "** Reading event: " << i << std::endl;
        
        double ratio = 1.0;
        double efficiency = 0.0;
        if (std::abs(nu_pdg) == 14) {
            ratio = (E_true <= 3.5) ? (nu_pdg == 14 ? exp(-0.02956 * pow(E_true, 2) - 0.005532 * E_true - 0.004352)
                                                : exp(-0.06302 * pow(E_true, 2) + 0.06550 * E_true - 0.03779))
                                    : (nu_pdg == 14 ? 0.747 * exp(-0.0351 * E_true) : 0.674 * exp(-0.0633 * E_true));
        } else if (std::abs(nu_pdg) == 12) {
            ratio = -0.076 * log(E_true) + 0.9561;
        }
        if (std::abs(nu_pdg) == 14)
            efficiency = (nu_pdg == 14 ? 0.70 : 0.60) * ratio;
        else if (std::abs(nu_pdg) == 12)
            efficiency = (nu_pdg == 12 ? 0.80 : 0.55) * ratio;
        else
            efficiency = 0.0;
        

        double E_sigma = 0.3 / sqrt(E_true);
        E_reco = randGen.Gaus(E_true, E_sigma);
        double cos_sigma = 0.1;
        costh_reco = randGen.Gaus(costh_true, cos_sigma);
        if (E_reco < 0)
            E_reco = 0;
        if (costh_reco < -1)
            costh_reco = -1;
        else if (costh_reco > 1)
            costh_reco = 1;
        
        if (randGen.Uniform() < efficiency) {

            int original_nu_pdg = nu_pdg;
            double conf_r = randGen.Uniform(0, 100);
            if (conf_r < 60.0) {
                nu_pdg = original_nu_pdg;
            } else if (conf_r < 80.0) {
                nu_pdg = -original_nu_pdg;
            } else if (conf_r < 95.0) {
                if (std::abs(original_nu_pdg) == 14)
                    nu_pdg = (original_nu_pdg > 0) ? 12 : -12;
                else if (std::abs(original_nu_pdg) == 12)
                    nu_pdg = (original_nu_pdg > 0) ? 14 : -14;
                else
                    nu_pdg = original_nu_pdg;
            } else {
                if (std::abs(original_nu_pdg) == 14)
                    nu_pdg = (original_nu_pdg > 0) ? -12 : 12;
                else if (std::abs(original_nu_pdg) == 12)
                    nu_pdg = (original_nu_pdg > 0) ? -14 : 14;
                else
                    nu_pdg = original_nu_pdg;
            }
            

            if (isnc == 1) {
                sample_branch = valor::conv::ksJFCAtmoNC;
            } else if (std::abs(nu_pdg) == 12) {
                sample_branch = (nu_pdg > 0) ? valor::conv::ksJFCAtmoNuECC : valor::conv::ksJFCAtmoNuEBarCC;
            } else if (std::abs(nu_pdg) == 14) {
                sample_branch = (nu_pdg > 0) ? valor::conv::ksJFCAtmoNuMuCC : valor::conv::ksJFCAtmoNuMuBarCC;
            } else {
                sample_branch = -1;
            }
            new_sample = sample_branch;
            beam = valor::conv::kbAtmo_Jiangmen;
            wght = ratio;
            
            reac_mod = valor::util::DetermineReactionModeReadJuno(
                nu_pdg, iscc, isnc,
                nP, nN, nPip, nPim, npi0,
                nKp, nKm, nK0, nEM, nOtherHad, nNucleus,
                isCharm, isStrange, isQE, isMEC,
                isNuEEl, isIMD, isIMDAnnihilation, isCoherent
            );
            
            outputTree.Fill();
            savedEvents++;
            
            if (new_sample == valor::conv::ksJFCAtmoNuECC)
                count_NuECC++;
            else if (new_sample == valor::conv::ksJFCAtmoNuEBarCC)
                count_NuEBarCC++;
            else if (new_sample == valor::conv::ksJFCAtmoNuMuCC)
                count_NuMuCC++;
            else if (new_sample == valor::conv::ksJFCAtmoNuMuBarCC)
                count_NuMuBarCC++;
        }
    }
    outputTree.Write();
    outputFile.Close();
    std::cout << "Total events saved: " << savedEvents << std::endl;
    std::cout << "Sample counts:" << std::endl;
    std::cout << "  ksJFCAtmoNuECC:    " << count_NuECC << std::endl;
    std::cout << "  ksJFCAtmoNuEBarCC: " << count_NuEBarCC << std::endl;
    std::cout << "  ksJFCAtmoNuMuCC:   " << count_NuMuCC << std::endl;
    std::cout << "  ksJFCAtmoNuMuBarCC:" << count_NuMuBarCC << std::endl;
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
