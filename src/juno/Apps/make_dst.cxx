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

using namespace genie;
using namespace valor;

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

  // First setup output
  TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
  TTree* outputTree = new TTree("juno_true_tree", "Processed Data");

  double Q2 = 0., W = 0., E_reco = 0., E_true = 0., costh_true = 0., costh_reco = 0., wght = 0.;
  double L_true = 0., L_reco = 0.;
  bool isnc = 0, iscc = 0;
  int reac_mod_wo_PID = -1, reac_mod = -1, nu_pdg = -1, event_id = 0, sample = -1, beam = -1;
  outputTree->Branch("Q2", &Q2);
  outputTree->Branch("W", &W);
  outputTree->Branch("enu_true", &E_true);
  outputTree->Branch("enu_reco", &E_reco);
  outputTree->Branch("costh_true", &costh_true);
  outputTree->Branch("costh_reco", &costh_reco);
  outputTree->Branch("L_true", &L_true);
  outputTree->Branch("L_reco", &L_reco);
  outputTree->Branch("nCC", &iscc);
  outputTree->Branch("nNC", &isnc);
  outputTree->Branch("reac_mode_wo_PID", &reac_mod_wo_PID);
  outputTree->Branch("reac_mode", &reac_mod);
  outputTree->Branch("nu_pdg", &nu_pdg);
  outputTree->Branch("event_id", &event_id);
  outputTree->Branch("wght", &wght);
  outputTree->Branch("sample", &sample);
  outputTree->Branch("beam", &beam);

  TChain* ghep_tree = new TChain("gtree");
  for (const auto& file : inputFiles) {
    ghep_tree->Add(file.c_str());
  }

  // Setup the DST to read GHEP trees.
  genie::NtpMCEventRecord* mcrec = nullptr;
  ghep_tree->SetBranchAddress("gmcrec", &mcrec);
  Long64_t nentries = ghep_tree->GetEntries();

  // Read the norm table and figure out expected numbers of entries
  std::vector<int> binCounts(gNormTable.size(), 0);
  for (Long64_t i = 0; i < nentries; ++i) {
    ghep_tree->GetEntry(i);
    if( !mcrec ) continue;

    TLorentzVector& nuP4 = *(mcrec->event->Probe()->GetP4()); // GetP4() NOT P4(), it's critical!
    double E_true_tmp = nuP4.Energy();
    double costh_true_tmp = -nuP4.Pz() / E_true_tmp;
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

      if (binIndex != -1) {
	binCounts[binIndex]++;
      }
    }

    // Crucially, delete the mcrec INSIDE HERE
    delete mcrec; mcrec = 0;
    
  } // loop over entries

  //
  // event loop
  //

  // Map to store counts for each reaction mode
  std::unordered_map<int, int> modeCounts;
  std::vector<Long64_t> nev_accept_by_mode;
  for (int im = 0 ; im < valor::conv::kNReacModes; im++) {
    nev_accept_by_mode.push_back(0);
  }

  Long64_t savedEvents = 0;
  TRandom3 randGen = valor::Rnd::Instance()->RndOther();
  for (Long64_t i = 0; i < nentries; ++i) {
    ghep_tree->GetEntry(i);
    if (gPrintEveryNEvents && !(i % gPrintEveryNEvents)) {
      VLOG("valor_juno", pNOTICE) << "** Reading event: " << i;
    } else {
      VLOG("valor_juno", pDEBUG) << "** Reading event: " << i
        <<" - so far accepted " << savedEvents;
    }

    if( !mcrec ) continue;

    genie::EventRecord  * event       = mcrec->event;
    genie::GHepParticle * neutrino    = event->Probe();
    TLorentzVector      * nuP4        = neutrino->GetP4();
    genie::Interaction  * interaction = event->Summary();
    genie::Kinematics   * kinematics  = interaction->KinePtr();
    genie::ProcessInfo  * procInfo    = interaction->ProcInfoPtr();
    
    E_true = nuP4->Energy();
    costh_true = -nuP4->Pz() / E_true;
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
    double desiredCount = 0.0;
    if (binIndex != -1) {
      const NormEntry& bestEntry = gNormTable[binIndex];
      if (neutrino != nullptr) {
	nu_pdg = neutrino->Pdg();
      }
      
      if (procInfo->IsWeakCC()) {
	iscc = 1; isnc = 0;
      } else if (procInfo->IsWeakNC()) {
	iscc = 0; isnc = 1;
      }
      if (isnc == 1) {
	desiredCount = bestEntry.nc;
      } else if (std::abs(nu_pdg) == 12) {
	desiredCount = (nu_pdg > 0) ? bestEntry.nue : bestEntry.nuebar;
      } else if (std::abs(nu_pdg) == 14) {
	desiredCount = (nu_pdg > 0) ? bestEntry.numu : bestEntry.numubar;
      }
    }
    double p = 1.0;
    if (binIndex != -1 && binCounts[binIndex] > 0) {
      p = desiredCount / binCounts[binIndex];
      if (p > 1.0) p = 1.0;
    }
    if (randGen.Uniform() < p) {
      E_reco = E_true;
      costh_reco = costh_true;

      /*
       * I will make some assumptions. 
       * Assume the Earth is a 2d circle, (i.e. the plane of neutrino propagation contains a great circle - rotational invariance implies that is always true, as phi is defined wrt the ray connecting JUNO to the centre of the Earth.
       * Assuming the radius of the earth to be R, and the neutrino to be emitted at height h in the atmosphere, (i.e. sitting on a sphere of radius R+h), Mathematica gives us the length as
       * L(h, R, theta) = sqrt( (R*cos(theta))^2 + h^2 + 2Rh ) - Rcos(theta)
       * TODO - extend to different baselines? Keep all height = 15 km now
       */

      L_reco = GetBaseline( costh_reco, earth_radius, neutrino_atmo_height );
      L_true = GetBaseline( costh_true, earth_radius, neutrino_atmo_height );

      Q2 = kinematics->Q2(true);
      W = kinematics->W(true);
                
      // Apply the reaction mode determination:
      valor::conv::ReacMode_t reac_mode = valor::util::DetermineReactionMode(mcrec->event, false);
      reac_mod = reac_mode;
                
      // Update reaction mode count
      modeCounts[reac_mod]++;
                
      bool is_nue = (std::abs(nu_pdg) == 12);
      bool is_antineutrino = (nu_pdg < 0);
      sample = static_cast<int>( InferSampleType( iscc, isnc, is_nue, is_antineutrino ) );
      int cnt = count(inputSamples.begin(), inputSamples.end(), sample);
      if(sample == -1){
        VLOG("valor_juno",pINFO) << "Current event has an undefined sample so skipping";
	continue; 
      }else if(cnt == 0){
        VLOG("valor_juno",pINFO) << "Current event has an undesired sample so skipping";
	continue;
      }
      beam = static_cast<int>(inputBeam);
      //beam = static_cast<int>( valor::conv::kbJUNO_Atmo );
      wght = p;
      outputTree->Fill();
      savedEvents++;
      nev_accept_by_mode[reac_mod] += 1;
    } // save this event

    // Clean up each event. It'll handle destruction and freeing of nested classes
    delete event;
    // But TLorentzVectors are not nested, so they need to be deleted
    delete nuP4;
  } // loop over events

  for (unsigned int im = 0 ; im < nev_accept_by_mode.size(); im++) {
    valor::conv::ReacMode_t m = (valor::conv::ReacMode_t)im;
    if(nev_accept_by_mode[im] > 0) VLOG("valor_sbn", pNOTICE) << "Selected "<<nev_accept_by_mode[im]<<" events of mode "<<valor::util::AsString(m);
  }

  double f = 100. * (double)savedEvents / (double)nentries;
  VLOG("valor_juno", pNOTICE) <<
    "Used " << f << "% (" << savedEvents << "/" << nentries
    << ") of input events.";

  // Write out the tree and close the file
  outputTree->Write();
  outputFile->Close();
  //delete ghep_tree;
  delete outputFile;
  
  /*// Output the event counts for each reaction mode
  std::cout << "Reaction Mode Event Counts:" << std::endl;
  for (const auto& kv : modeCounts) {
    std::cout << "Mode " << kv.first << ": " << kv.second << " events." << std::endl;
  }*/
  
  return 0;
}
//--------------------------------------------------------------------------------------------------

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
  for( const auto& entry : std::filesystem::recursive_directory_iterator(arg_input_directory) ) {
    if( entry.is_regular_file() && entry.path().extension() == ".root" )
      inputFiles.emplace_back( entry.path() );
  }

  processEvents(inputFiles, inputSamples, inputBeam, arg_output);

  VLOG("valor_sbn", pNOTICE)  << "Done!";

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
