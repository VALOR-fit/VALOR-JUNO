//____________________________________________________________________________
/*!

\app       valor-juno_make_dst

\brief     Generate DSTs for JUNO software

\options     --input-normalisation : Location of normalisation file (for event rates) [mandatory]
             --input-directory     : Path to directory containing ghep input files [mandatory]
	     --output              : Path to output [default: valor_juno_dst.root]
	     --mesg-thresholds     : Message thresholds [default: verbosity.xml]

\author    John Plows <kplows \at liverpool.ac.uk>
           Beth Slater <b.slater2 \at liverpool.ac.uk>
           University of Liverpool

           Yaoqi Cao <yaoqi.cao \at wariwck.ac.uk>
	   Ziou He <zh316 \at cantab.ac.uk>
	   University of Warwick and Chinese Academy of Sciences

\created   March 12, 2025

\copyright Copyright (c) 2010-2025, The VALOR group. All rights reserved

*/
//____________________________________________________________________________

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <fstream>
#include <cstdlib>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TRandom.h>

// VALOR
#include "valor/Utilities/Env.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/NameRegistry.h"
#include "valor/Conventions/Constants.h"
#include "valor/Conventions/ReacMode.h"
#include "valor/Conventions/Sample.h"
#include "valor/Conventions/BeamConf.h"
#include "valor/Messenger/Messenger.h"
#include "valor/CmdLine/CmdLnArgList.h"
#include "valor/Utilities/UtlNumStat.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Conventions/ModeList.h"
#include "valor/GENIEUtilities/GENIEIncludes.h"
#include "valor/GENIEUtilities/UtlEventRec.h"

// GENIE
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/ProcessInfo.h"

// memory tracking
#include "TObjectTable.h"
#include <filesystem>

using namespace genie;
using namespace valor;
namespace fs = std::filesystem;

// A struct to let one read the normalisation table
struct NormEntry {
  double energyLow;
  double costhLow;
  double numu;
  double numubar;
  double nue;
  double nuebar;
  double nc;
};

// Global variables
const double neutrino_atmo_height = 15.0; // km
const double earth_radius = 6341.0; // km

std::vector<NormEntry> gNormTable;
const int gPrintEveryNEvents = 2000;

//--------------------------------------------------------------------------------------------------
// Function to read normalization table from txt file
bool readNormTable(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) {
      VLOG("valor_juno-make_dst", pERROR) << "Error opening normalization file: " << filename;
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
//--------------------------------------------------------------------------------------------------
// Read out some information...
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

//--------------------------------------------------------------------------------------------------
// SampleTypeToString is Sample.cxx::SampleTypeToString (with the correct JUNO samples)
//--------------------------------------------------------------------------------------------------
valor::conv::Sample_t InferSampleType(bool is_cc, bool is_nc, bool is_nue, bool is_antineutrino) {
  if( is_cc ) {
    if ( is_nue && !is_antineutrino) return valor::conv::ksJUNO_FC_CC_Nue;
    if ( is_nue &&  is_antineutrino) return valor::conv::ksJUNO_FC_CC_Nuebar;
    if (!is_nue && !is_antineutrino) return valor::conv::ksJUNO_FC_CC_Numu;
    if (!is_nue &&  is_antineutrino) return valor::conv::ksJUNO_FC_CC_Numubar;
  } else if( is_nc ) {
    if ( is_nue && !is_antineutrino) return valor::conv::ksJUNO_FC_NC_Nue;
    if ( is_nue &&  is_antineutrino) return valor::conv::ksJUNO_FC_NC_Nuebar;
    if (!is_nue && !is_antineutrino) return valor::conv::ksJUNO_FC_NC_Numu;
    if (!is_nue &&  is_antineutrino) return valor::conv::ksJUNO_FC_NC_Numubar;
  }

  return valor::conv::ksUndefined;
}
//--------------------------------------------------------------------------------------------------
double GetBaseline(double cosTheta, double R, double h) {
  // production height = h, earth radius = R, neutrino incidence angle = theta
  return std::sqrt( std::pow( R * cosTheta, 2.0 ) + 2.0 * h * R + h * h ) - R * cosTheta;
}
//--------------------------------------------------------------------------------------------------
int processEvents( const std::vector<std::string> & inputFiles,const std::vector<conv::Sample_t> &  inputSamples, const conv::BeamConf_t & inputBeam, const std::string & outputFileName ) {
    TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    outputFile->cd();
    TTree* outputTree = new TTree("juno_true_tree", "Processed Data");
    TChain* ghep_tree = new TChain("gtree");
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
            TLorentzVector& nuP4 = *mcrec->event->Probe()->GetP4();
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
    double Q2, W, E_true, E_reco, costh_true, costh_reco, wght, L_true,L_reco;
    int isnc = 0, iscc = 0, nkaon = 0, npip = 0, npim = 0, npi0 = 0;
    int reac_mod_wo_PID = -1;
    int reac_mod = -1, nu_pdg = -1, event_id = 0, sample_branch = -1, beam = -1;
    outputTree->Branch("Q2", &Q2);
    outputTree->Branch("W", &W);
    outputTree->Branch("enu_true", &E_true);
    outputTree->Branch("enu_reco", &E_reco);
    outputTree->Branch("L_true", &L_true);
    outputTree->Branch("L_reco", &L_reco);
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
        genie::EventRecord  * event       = mcrec->event;
        genie::GHepParticle * neutrino    = event->Probe();
        TLorentzVector      &nuP4         = *(neutrino->GetP4());
        genie::Interaction  * interaction = event->Summary();
        genie::Kinematics   * kinematics  = interaction->KinePtr();
        genie::ProcessInfo  * procInfo    = interaction->ProcInfoPtr();
        if (gPrintEveryNEvents && !(i % gPrintEveryNEvents)) {
            // Logging can be added here if needed
        }
        if (mcrec && mcrec->event) {
            //const TLorentzVector& nuP4 = *(mcrec->event->Probe()->P4());
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
                nu_pdg = (neutrino ? neutrino->Pdg() : -1);
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
                    L_reco = GetBaseline( costh_reco, earth_radius, neutrino_atmo_height );
                    L_true = GetBaseline( costh_true, earth_radius, neutrino_atmo_height );
                    Q2 = mcrec->event->Summary()->Kine().Q2(true);
                    W  = mcrec->event->Summary()->Kine().W(true);
                    valor::conv::ReacMode_t reac_mode = valor::util::DetermineReactionMode(event, false);
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
                    sample_branch = static_cast<int>( InferSampleType( iscc, isnc, is_nue, is_antineutrino ) );
                    int cnt = count(inputSamples.begin(), inputSamples.end(), sample_branch);
                    if(sample_branch == -1){
                        VLOG("valor_juno",pINFO) << "Current event has an undefined sample so skipping";
	                continue; 
                    }else if(cnt == 0){
                        VLOG("valor_juno",pINFO) << "Current event has an undesired sample so skipping";
	                continue;
                    }
                    beam = static_cast<int>(inputBeam);
                    wght = 1.0;
                    outputTree->Fill();
                    savedEvents++;
                }
            }
        }
    }
 
    double f = 100. * (double)savedEvents / (double)nentries;
    VLOG("valor_juno", pNOTICE) <<
      "Used " << f << "% (" << savedEvents << "/" << nentries
      << ") of input events.";

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

    time_t scriptstartrawtime;
    struct tm * scriptstarttimeinfo;
    time (&scriptstartrawtime);
    scriptstarttimeinfo = localtime (&scriptstartrawtime);
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << " Script start local time and date:  " << asctime(scriptstarttimeinfo);
    std::cout << "-----------------------------------------------------------" << std::endl;
  
    // Return usage if user asks for -h
    for( int i = 1 ; i < argc ; i++ ) {
      if( std::string(argv[i]) == "-h" ) {
        VLOG("valor_juno", pFATAL) << "Usage:" << " valor-juno_make_dst "
                   << " --input-normalisation path_to_norm_table.txt "
                   << " --input-directory path_to_ghep_files"
                   << " [--output path_to_output] [--mesg-thresholds path_to_verbosity]";
        return 0;
      }
    }
    
    CmdLnArgList arg;
    arg.Read(argc,argv, "$VALOR_EXPT/config/apps/make_dst.arg");
    VLOG("valor_juno", pNOTICE) << arg;
    
    string arg_input_normalisation         = arg.GetString("input-normalisation");
    string arg_input_directory             = arg.GetString("input-directory");
    string arg_output                      = arg.GetString("output");
    string arg_input_beam                  = arg.GetString("beam");
    std::vector<string>  arg_input_samples = arg.GetStringVec("samples");

  
    std::string thresholds = util::CompleteConfigPath(arg.GetString("mesg-thresholds"));
    valor::Messenger::Instance()->Configure(thresholds);
  
    if (!readNormTable(arg_input_normalisation)) {
      return 1;
    }
  
    // Get the sample and beam configuration
    std::vector<std::string> samples;
    std::vector<conv::Sample_t> inputSamples;
    if(std::find(arg_input_samples.begin(), arg_input_samples.end(), "all") != arg_input_samples.end()) {
      //samples = {"CCincl_fc", "numuCC_fc", "numubarCC_fc", "nueCC_fc", "nuebarCC_fc", "NCincl_fc", "numuNC_fc", "numubarNC_fc", "nueNC_fc", "nuebarNC_fc"};
      samples = {"numuCC_fc", "numubarCC_fc", "nueCC_fc", "nuebarCC_fc", "numuNC_fc", "numubarNC_fc", "nueNC_fc", "nuebarNC_fc"};
      for(string sample : samples){
        inputSamples.push_back(util::SampleFromString(sample));
      }
    }else{
      for(string sample : arg_input_samples){
        conv::Sample_t curr_sample = util::SampleFromString(sample);
        assert(curr_sample != conv::ksUndefined);
        inputSamples.push_back(curr_sample);
      }
    }
    conv::BeamConf_t inputBeam = util::BeamConfFromString(arg_input_beam);
    assert(inputBeam != conv::kbUndefined);
  
    std::vector<std::string> inputFiles;
    for (const auto& entry : fs::recursive_directory_iterator(arg_input_directory)) {
      if (entry.is_regular_file() && entry.path().extension() == ".root")
          inputFiles.emplace_back(entry.path().string());
    }
  
    processEvents(inputFiles, inputSamples, inputBeam, arg_output);
  
    VLOG("valor_juno", pNOTICE)  << "Done!";
  
    time_t endrawtime;
    struct tm * endtimeinfo;
    time (&endrawtime);
    endtimeinfo = localtime (&endrawtime);
    time_t total_time = endrawtime-scriptstartrawtime;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << " End local time and date  :  " << asctime(endtimeinfo);
    std::cout << " Total time taken         :  " << total_time/static_cast<double>(60) << " minutes" << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
  
    return 0;
  }