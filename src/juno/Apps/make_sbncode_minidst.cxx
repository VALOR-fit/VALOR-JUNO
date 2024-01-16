//____________________________________________________________________________
/*!

  \app       valor_sbn-make_sbncode_minidst

  \brief     Converts input SBNCode analysis DSTs to the VALOR mini-DST format.

  \syntax    valor_sbn-make_sbncode_minidst
                --input  input_analysis_file.root
                --output output_file_containing_valor_minidst.root
               [--write-universes] optionally write out the systematic variations
               [--disable-checks]
               [--off-axis-min off_axis_lower_edge] optionally give offaxis angle
                                                    range in deg; default 0.0
               [--off-axis-max off_axis_upper_edge] default 90.0

  \author    Rhiannon Jones <rhiannon.jones \at liverpool.ac.uk>
             University of Liverpool

  \created   March 2019

  \copyright Copyright (c) 2010-2020, The VALOR group. All rights reserved.

*/
//____________________________________________________________________________

// System
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <utility>
#include <algorithm>
// ROOT
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TSystem.h>
// VALOR/SDK
#include "valor/Conventions/Constants.h"
#include "valor/Conventions/BeamConf.h"
#include "valor/Conventions/ReacMode.h"
#include "valor/Conventions/Sample.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Utilities/Env.h"
#include "valor/CmdLine/CmdLnArgList.h"
#include "valor/Utilities/UtlNumStat.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/Rnd.h"
#include "valor/GENIEUtilities/GENIEIncludes.h"
#include "valor/GENIEUtilities/GENIEIncludes.h"
#include "valor/GENIEUtilities/UtlEventRec.h"
// SBNcode - See for example Gray or Rhiannon's branch (not the default)
#include "core/Event.hh"

using std::string;
using std::ostringstream;

void GetCommandLineArgs (int /*argc*/, char ** /*argv*/);
static const string gTreeName       = "sbnreco";
static const string gSubRunTreeName = "sbnsubrun";

static const int gPrintEveryNEvents   = 50000;

static const int max_n_universes = 10000;

// Should move this to a header file (eg DSTTypes or something).
//~ typedef Double_t variation_output_type;
//~ const char variation_output_type_root = 'D';
typedef Float_t variation_output_type;
const char variation_output_type_root = 'F';

// The enum simb::int_type_ from LarSoft (copied from MCNeutrino.h, in v06_05_00)
typedef enum ELarMode{
  kLArUnknownInteraction = -1,
  kLArQE                 =  0,
  kLArRes                =  1,
  kLArDIS                =  2,
  kLArCoh                =  3,
  kLArCohElastic         =  4,
  kLArElectronScattering =  5,
  kLArIMDAnnihilation    =  6,
  kLArInverseBetaDecay   =  7,
  kLArGlashowResonance   =  8,
  kLArAMNuGamma          =  9,
  kLArMEC                = 10,
  kLArDiffractive        = 11,
  kLArEM                 = 12,
  kLArWeakMix            = 13
} LArMode_t;

std::vector<TString> nonResVariations = 
{
"genie_NonResRvpCC1pi_Genie",
"genie_NonResRvnCC1pi_Genie",
"genie_NonResRvbarnCC1pi_Genie",
"genie_NonResRvbarpCC1pi_Genie",
"genie_NonResRvpCC2pi_Genie",
"genie_NonResRvnCC2pi_Genie",
"genie_NonResRvbarnCC2pi_Genie",
"genie_NonResRvbarpCC2pi_Genie",
"genie_NonResRvpNC1pi_Genie",
"genie_NonResRvnNC1pi_Genie",
"genie_NonResRvbarnNC1pi_Genie",
"genie_NonResRvbarpNC1pi_Genie",
"genie_NonResRvpNC2pi_Genie",
"genie_NonResRvnNC2pi_Genie",
"genie_NonResRvbarnNC2pi_Genie",
"genie_NonResRvbarpNC2pi_Genie",
};

bool CheckSample(const event::RecoInteraction &reco, const bool &numu, const valor::conv::Sample_t &sample, const bool &other);
valor::conv::Sample_t GetOtherSample(const valor::conv::Sample_t &s, const bool &numu);
void FillVariationBranches(const TString &vname,
                           const size_t &n_universes,
                           TObjArray &variation_names,
                           std::map<TString,int> &variation_nentries,
                           std::map<TString,variation_output_type *> &variation_entries,
                           TTree *valor_minidst);
void GetCorrectVariationName(const int nuc_pdg, const bool iscc, TString &vname);
void GetPair(const TString &vname, std::pair<TString,TString> &new_vnames);

//___________________________________________________________________
int main(int argc, char ** argv)
{
  time_t scriptstartrawtime;
  struct tm * scriptstarttimeinfo;
  time (&scriptstartrawtime);
  scriptstarttimeinfo = localtime (&scriptstartrawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Script start local time and date:  " << asctime(scriptstarttimeinfo);
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Set GENIE print level
  genie::GHepRecord::SetPrintLevel(genie::RunOpt::Instance()->EventRecordPrintLevel());

  using namespace valor;

  CmdLnArgList arg;
  arg.Read(argc,argv, "$VALOR_EXPT/config/apps/make_sbncode_minidst.arg");
  LOG("valor_sbn", pNOTICE) << arg;

  std::string thresholds = util::CompleteConfigPath(arg.GetString("mesg-thresholds"));
  Messenger::Instance()->Configure(thresholds);

  string opt_beam_name      = arg.GetString("beam");
  string opt_sample_name    = arg.GetString("sample");
  string opt_output_name    = arg.GetString("output");
  string opt_input_name     = arg.GetString("input");
  string opt_oscillation    = arg.GetString("oscillation");
  bool   opt_disable_checks = arg.GetBool("disable-checks");
  bool   write_universes    = arg.GetBool("write-universes");
  bool   apply_weights      = !(arg.GetBool("dont-apply-weights"));
  bool   no_truth           = arg.GetBool("no-truth");
  bool   other_sample       = arg.GetBool("other-sample");
  bool   fix_nonres         = arg.GetBool("fix-nonres");
  bool   status_code_set    = arg.GetBool("status-code-set");
  int    max_n_events       = arg.GetInt("max-n-events");
  double off_axis_bin_min   = arg.GetDouble("off-axis-min");
  double off_axis_bin_max   = arg.GetDouble("off-axis-max");

  // Get the sample and beam configuration
  conv::Sample_t sample = util::SampleFromString(opt_sample_name);
  assert(sample != conv::ksUndefined);

  conv::BeamConf_t beam = util::BeamConfFromString(opt_beam_name);
  assert(beam != conv::kbUndefined);

  if(no_truth && write_universes)
    VLOG("valor_sbn", pFATAL) << " Error: Cannot write universes if there is no truth information supplied ";

  if(apply_weights)
    VLOG("valor_sbn", pNOTICE) << " Chosen to apply valor_minidst_wght ";
  else
    VLOG("valor_sbn", pNOTICE) << " Chosen not to apply valor_minidst_wght";

  // If we want to define 'other' as anything but the sample defined in the arguments, set this to true
  if(other_sample)
    VLOG("valor_sbn", pNOTICE) << " Chosen to look for events which do NOT fall under the definition of the sample, "
                               << util::AsString(sample);

  // Check if we are looking at nue_disappearance from the opt_oscillation argument
  // so that we reject the nue appearance oscillated events later
  string nue_dis = "nue-disappearance";
  size_t found_nue_dis = opt_oscillation.find(nue_dis);
  bool nue_disappearance = false;
  if(found_nue_dis != string::npos){
    nue_disappearance = true;
    VLOG("valor_sbn", pINFO)
        << "Looking at a nue disappearance event from the oscillation option: " << opt_oscillation;
  }
  else
    VLOG("valor_sbn", pINFO)
        << "NOT looking at a nue disappearance event from the oscillation option: " << opt_oscillation;

  // Determine sample-types
  bool is_nueint        = false;
  bool is_nueosc        = false;
  bool is_nuenumu       = false;
  bool is_numu          = false;
  bool is_nue           = false;
  bool is_dirt          = false;
  bool is_cosmic        = false;
  string str_numu       = "Numu";
  string str_nue        = "Nue";
  string str_int        = "_Int";
  string str_osc        = "_Osc";
  string str_dirt       = "_Dirt";
  string str_cosmic     = "_Nu_Cosmic_Dirt";

  //Find substrings and check the selected nuE sample-type if necessary
  size_t found_nue  = opt_input_name.find(str_nue);
  size_t found_numu = opt_input_name.find(str_numu);

  // Get selected nuMu from sample name
  if(found_nue == string::npos && found_numu != string::npos)
    is_numu = true;

  if(found_nue != string::npos){
    is_nue = true;
    // We are looking at a nue sample, check if osc, intrinsic or numu
    size_t found_nueosc  = opt_input_name.find(str_osc);
    size_t found_nueint  = opt_input_name.find(str_int);
    size_t found_nuenumu = opt_input_name.find(str_numu);
    size_t found_nuedirt = opt_input_name.find(str_dirt);                                                                                                          
    size_t found_cosmic  = opt_input_name.find(str_cosmic); 
    
    if(found_cosmic != string::npos){                                                                                                                                    
      is_cosmic = true;  
    } 
    if(found_nueosc != string::npos && found_nueint == string::npos && found_nuenumu == string::npos)
      is_nueosc = true;
    else if(found_nueosc != string::npos){
      VLOG("valor_sbn", pFATAL)
        << "Error: Found selected nuE oscillated file with hint of intrinsic and/or numu in the name - misleading: " << opt_input_name;
      exit(1);
    }
    if(found_nueint != string::npos && found_nueosc == string::npos && found_nuenumu == string::npos)
      is_nueint = true;
    else if(found_nueint != string::npos){
      VLOG("valor_sbn", pFATAL)
        << "Error: Found selected nuE intrinsic file with hint of oscillated and/or numu in the name - misleading: " << opt_input_name;
      exit(1);
    }
    if(found_nuenumu != string::npos && found_nueosc == string::npos && found_nueint == string::npos)
      is_nuenumu = true;
    else if(found_nuenumu != string::npos){
      VLOG("valor_sbn", pFATAL)
        << "Error: Found selected nuE numu file with hint of oscillated and/or intrinsic in the name - misleading: " << opt_input_name;
      exit(1);
    }
    if(found_nuedirt != string::npos && found_nueosc == string::npos && found_nueint == string::npos && found_nuenumu == string::npos)                                
      is_dirt = true;                                                                                                                                                 
    else if(found_nuedirt != string::npos){                                                                                                                           
      VLOG("make_sbncode_minidst", pFATAL)                                                                                                                               
        << "Error: Found appearance dirt file with hint of oscillated and/or intrinsic and/or numu in the name - misleading: " << opt_input_name;                     
      exit(1);                                                                                                                                                           
    }
  }
  

  bool found_nues[] = {is_nueosc, is_nueint, is_nuenumu, is_dirt};
  if(found_nue && std::count(found_nues, found_nues+sizeof(found_nues), 1) > 1){
    VLOG("valor_sbn", pFATAL)
      << "Error: Found multiple selected nuE types in filename - misleading: " << opt_input_name;
    exit(1);
  }

  // Make sure we are finding at least 1 of numu or nue in the name
  bool found[] = {is_numu, is_nueosc, is_nueint, is_nuenumu, is_dirt, is_cosmic};
  if(std::count(found, found+sizeof(found), 1) == 0){
    VLOG("valor_sbn", pFATAL)
      << "Error: Found hint of neither selected nuE (" << str_nue << ", " << str_int << ", " << str_osc << ", " << str_dirt << "," << str_cosmic << ") nor selected nuMu ( " << str_numu << ") in filename - misleading: " << opt_input_name;
    exit(1);
  }

  // CHECKING OUTPUT OF NEW CODE
  std::cout << " Input filename : " << opt_input_name << " - ";
  if(is_nueosc)
    std::cout << " Selected NuE, Oscillated " << std::endl;
  if(is_nueint)
    std::cout << " Selected NuE, Intrinsic " << std::endl;
  if(is_nuenumu)
    std::cout << " Selected NuE, NuMu " << std::endl;
  if(is_numu)
    std::cout << " Selected NuMu, NuMu " << std::endl;
  if(is_dirt && !is_cosmic)
    std::cout << " Dirt " << std::endl;
  if(is_dirt && is_cosmic)
    std::cout << " Cosmics with Dirt" << std::endl;

  // Get the POT of the sample from the subrun tree
  double mcsample_pot = -1.;
  TChain * subrun_tree = new TChain(gSubRunTreeName.c_str());
  subrun_tree->Add(opt_input_name.c_str());
  if(!subrun_tree || !subrun_tree->GetEntries()) {
    VLOG("valor_sbn", pFATAL)
      << "Error: Unable to read " << gSubRunTreeName
      <<" from "<< opt_input_name;
    exit(1);
  }

  Long64_t nev_subrun = subrun_tree->GetEntries();
  VLOG("valor_sbn", pNOTICE)
    << "Input subrun tree has " << nev_subrun << ((nev_subrun==1) ? " entry." : " entries.");

  TBranch *b = subrun_tree->GetBranch("totgoodpot");
  unsigned int npot = subrun_tree->GetEntries();
  for(unsigned int i = 0; i < npot; ++i){
    b->GetEntry(i);
    mcsample_pot += b->GetLeaf("totgoodpot")->GetValue();
  }

  if(mcsample_pot == -1.){
    VLOG("valor_sbn", pFATAL)
      << "Error: Not able to access MCSample POT from " << gSubRunTreeName;
    exit(1);
  }
  VLOG("valor_sbn", pNOTICE)
    << "POT in the input file is " << mcsample_pot;


  int    valor_minidst_sample       = sample;
  int    valor_minidst_beam         = -1;  // see Conventions/BeamConf.h
  int    valor_minidst_reac_mode    = -1;  // see Conventions/ReacMode.h
  int    valor_minidst_reac_type    = -1;  // see Conventions/ReacMode.h
  int    valor_minidst_target       = -1;  // nuclear target -  GENIE PDG convention
  double valor_minidst_baseline     = -1.; // true baseline of the event
  int    valor_minidst_init_nupdg   = -1;  // Initial PDG
  int    valor_minidst_final_nupdg  = -1;  // Final PDG
  double valor_minidst_enu_true     = -1.; // true neutrino energy, in GeV
  double valor_minidst_Q2_true      = -1.; // true neutrino energy, in GeV
  double valor_minidst_enu_reco     = -1.; // reco neutrino energy, in GeV
  double valor_minidst_q0_reco      = -1.; // reco lepton energy, in GeV
  double valor_minidst_Theta_reco   = -1.; // reco lepton opening angle
  double valor_minidst_y_reco       = -1.; // inelasticity
  double valor_minidst_wght_pot     = -1.; // POT weight
  double valor_minidst_wght         = -1.; // misc weight
  double valor_minidst_posx_true    = -999.; // true vertex x positon
  double valor_minidst_posy_true    = -999.; // true vertex y positon
  double valor_minidst_posz_true    = -999.; // true vertex z positon
  double valor_minidst_offaxis_true = -1.0; // true offaxis angle to beam centre

  TFile fout(opt_output_name.c_str(),"recreate");
  fout.cd();

  TTree * valor_minidst =
    new TTree("valor_sbn_minidst","Mini-DST input for VALOR/SBN fit");

  valor_minidst->Branch("beam"        , &valor_minidst_beam        ,"beam/I"        );
  valor_minidst->Branch("sample"      , &valor_minidst_sample      ,"sample/I"      );
  valor_minidst->Branch("reac_mode"   , &valor_minidst_reac_mode   ,"reac_mode/I"   );
  valor_minidst->Branch("reac_type"   , &valor_minidst_reac_type   ,"reac_type/I"   );
  valor_minidst->Branch("init_pdg"    , &valor_minidst_init_nupdg  ,"init_pdg/I"    );
  valor_minidst->Branch("final_pdg"   , &valor_minidst_final_nupdg ,"final_pdg/I"   );
  valor_minidst->Branch("baseline"    , &valor_minidst_baseline    ,"baseline/D"    );
  valor_minidst->Branch("target"      , &valor_minidst_target      ,"target/I"      );
  valor_minidst->Branch("enu_true"    , &valor_minidst_enu_true    ,"enu_true/D"    );
  valor_minidst->Branch("Q2_true"     , &valor_minidst_Q2_true     ,"q2_true/D"     );
  valor_minidst->Branch("enu_reco"    , &valor_minidst_enu_reco    ,"enu_reco/D"    );
  valor_minidst->Branch("q0_reco"     , &valor_minidst_q0_reco     ,"q0_reco/D"     );
  valor_minidst->Branch("Theta_reco"  , &valor_minidst_Theta_reco  ,"Theta_reco/D"  );
  valor_minidst->Branch("y_reco"      , &valor_minidst_y_reco      ,"y_reco/D"      );
  valor_minidst->Branch("wght_pot"    , &valor_minidst_wght_pot    ,"wght_pot/D"    );
  valor_minidst->Branch("wght"        , &valor_minidst_wght        ,"wght/D"        );
  valor_minidst->Branch("posx_true"   , &valor_minidst_posx_true   ,"posx_true/D"   );
  valor_minidst->Branch("posy_true"   , &valor_minidst_posy_true   ,"posy_true/D"   );
  valor_minidst->Branch("posz_true"   , &valor_minidst_posz_true   ,"posz_true/D"   );
  valor_minidst->Branch("offaxis_true", &valor_minidst_offaxis_true,"offaxis_true/D");

  // open the ROOT file and get the TTree & its header

  TChain * analysis_tree = new TChain(gTreeName.c_str());
  analysis_tree->Add(opt_input_name.c_str());
  if(!analysis_tree || !analysis_tree->GetEntries()) {
    VLOG("valor_sbn", pFATAL)
      << "Error: Unable to read " << gTreeName
      <<" from "<< opt_input_name;
    exit(1);
  }

  Long64_t nev = analysis_tree->GetEntries();
  VLOG("valor_sbn", pNOTICE)
    << "Input analysis tree has " << nev << ((nev==1) ? " entry." : " entries.");
  if (max_n_events >= 0) {
    nev = std::min((Long64_t)max_n_events,nev);
    VLOG("valor_sbn", pNOTICE) <<
      "Reading "<< nev << ((nev==1) ? " entry." : " entries.");
  }

  event::RecoEvent *pevent = nullptr;
  analysis_tree->SetBranchAddress("reco_events", &pevent);

  std::map<TString,variation_output_type *> variation_entries;
  std::map<TString,int>     variation_nentries;
  //~ std::map<std::string,TBranch *> variation_branches;
  // First, let's read one true event map to get the universe details
  if (!no_truth && write_universes) {
    analysis_tree->GetEntry(0);
    event::Interaction * ptruth = pevent->GetTruth();
    assert(ptruth);
    const auto & weights = ptruth->weights; // map<string,whatever>
    TObjArray variation_names;

    // Boolean for if we have already added the replacement non-res branches
    bool replacedNonRes = false;
    for (const auto & variation : weights) {
      // make a list of the variation names
      const TString & vname(variation.first);
      
      // Get the number of universes
      const size_t n_universes = variation.second.size();
     
      // If the variation is a non-resonant background parameter, add 4 branches to the tree
      if(vname.Contains("NonRes") && fix_nonres){
        if(!replacedNonRes){
          for(const TString &new_vname : nonResVariations){
            // Now translate to the new list of parameters, conversion table is here: https://github.com/VALOR-fit/VALOR-SBN/wiki/XSec-systematics-2020a
            FillVariationBranches(new_vname,n_universes,variation_names,variation_nentries,variation_entries,valor_minidst);
            replacedNonRes = true;
            LOG("valor_sbn",pINFO) << "Adding new parameter: " << new_vname << " branches in replacement of the current non-resonant branch definitions";
          }
        }
        else continue;
      }
      // Otherwise add the branches as-named
      else{
        FillVariationBranches(vname,n_universes,variation_names,variation_nentries,variation_entries,valor_minidst);
      }
    }
    variation_names.Write("variation_names",TObject::kSingleKey);
  }

  //
  // event loop
  //

  Long64_t nev_accept = 0;

  std::vector<int> nev_accept_by_mode;
  for (int im = 0 ; im < valor::conv::kNReacModes; im++) {
    nev_accept_by_mode.push_back(0);
  }

  // Commented out now we only have one tree per sample. But remember if that changes!
  //~ int last_tree_number = -1; // store which number file we're in, so we can count each tree's weight once

  for(Long64_t iev = 0; iev < nev; iev++) {

    analysis_tree->GetEntry(iev);

    if (gPrintEveryNEvents && !(iev%gPrintEveryNEvents)) {
      VLOG("valor_sbn", pNOTICE) << "** Reading event: " << iev
        <<" - so far accepted " << nev_accept;
    }
    else {
      VLOG("valor_sbn", pDEBUG) << "** Reading event: " << iev
        <<" - so far accepted " << nev_accept;
    }

    // reset the systematic universes
    if (!no_truth && write_universes) {
      for (const auto & entry : variation_entries) {
        variation_output_type * arr = entry.second;
        for (int ii = 0 ; ii < max_n_universes ; ii++) {
          arr[ii] = static_cast<variation_output_type>(conv::kBadResult);
        }
      }
      for (auto & entry : variation_nentries) {
        int & nentries = entry.second;
        nentries = -1;
      }
    } // Variations

    double enu_reco    =  -999.;
    double q0_reco     =  -999.;
    double Theta_reco  =  -999.;
    double reco_weight = 1.;

    analysis_tree->GetEntry(iev);

    event::RecoEvent & event = *pevent;
    event::RecoInteraction &reco = event.reco;
    event::Interaction * ptruth = event.GetTruth();
    event::Interaction & truth = *ptruth;

    // Before we do anything, make sure the current event fits the description of the sample
    if(!CheckSample(reco,is_numu,sample,other_sample)){
      LOG("valor_sbn",pINFO)   << "Current event has \n" 
                               << reco.npion << " pions, \n " 
                               << reco.npi0 << " pi0s and \n " 
                               << reco.nproton << " protons and does not match the desired sample type: " 
                               << util::AsString(sample) << "\n Skipping.";
      continue;
    }

    // If there is truth in the file, access the variables
    if(!no_truth){
      int iscc           = false;
      int isnc           = false;
      int nu_pdg         = -1;
      int init_pdg       = -1;
      int nuc_pdg        = -1; // PDG of the struck nucleon
      double baseline    = -999.;
      double enu_true    = -999.;
      double Q2s         = -999.;

      int nipip    = 0;
      int nipim    = 0;
      int nipi0    = 0;
      int nikp     = 0;
      int nikm     = 0;
      int nik0     = 0;
      int niak0    = 0;
      int nin      = 0;
      int nian     = 0;
      int nip      = 0;
      int niap     = 0;
      int niem     = 0;
      int niep     = 0;
      int ninu     = 0;
      int nimu     = 0;
      int nimup    = 0;
      int nig      = 0;
      int niother  = 0;
      int ninuc    = 0;

      int larsoft_mode = -1;

      double posx_true = -999.0;
      double posy_true = -999.0;
      double posz_true = -999.0;

      double offaxis_true = -1.0;

      /*
         analysis_tree->GetEntry(iev);

         event::RecoEvent & event = *pevent;
         event::RecoInteraction &reco = event.reco;
         */
      if(ptruth==nullptr) {
        VLOG("valor_sbn", pNOTICE)
          << "The number of truth entries is 0, skipping event "<<iev;
        continue;
      }

      std::vector<event::FinalStateParticle> &fsps = truth.finalstate;
      iscc         = truth.neutrino.iscc;
      isnc         = truth.neutrino.isnc;
      baseline     = truth.neutrino.baseline;
      nu_pdg       = truth.neutrino.pdg;
      enu_true     = truth.neutrino.energy;
      Q2s          = truth.neutrino.Q2;
      larsoft_mode = truth.neutrino.genie_intcode;
      posx_true    = truth.neutrino.position.x();
      posy_true    = truth.neutrino.position.y();
      posz_true    = truth.neutrino.position.z();

      //Calculate the off angle position
      double beamx = 45.7; //read this from somewhere sensible!
      double beamy = 0.0;
      double exp_baseline = 11000.0;
      double offaxis_true_cm = sqrt( (posx_true-beamx)*(posx_true-beamx) + (posy_true-beamy)*(posy_true-beamy));
      offaxis_true = atan(offaxis_true_cm/exp_baseline) * (180/M_PI);


      // Particle counting
      bool nucleonPDGSet = false;
      for(unsigned int p = 0; p < fsps.size(); ++p){
        if(status_code_set){
          if(fsps[p].status_code == 11){
            if(nucleonPDGSet){
              LOG("valor_sbn",pNOTICE) << "\nIn event " << iev
                << "\nThe struck nucleon pdg has already been set as: " << nuc_pdg << ", will take the first answer."
                << "\nThis nucleon has pdg: " << fsps[p].pdg;
            }
            else{
              LOG("valor_sbn",pINFO) << "\nIn event " << iev
                << "\nThe struck nucleon has pdg: " << fsps[p].pdg;
              nuc_pdg = fsps[p].pdg;
              nucleonPDGSet = true;
            }
          }
        }
        // Make sure we are counting initial states for the interaction definition 
        if (status_code_set && fsps[p].status_code != 14 && fsps[p].status_code != 12) {continue;}
        if     (fsps[p].pdg == genie::kPdgProton      ) nip++;
        else if(fsps[p].pdg == genie::kPdgAntiProton  ) niap++;
        else if(fsps[p].pdg == genie::kPdgNeutron     ) nin++;
        else if(fsps[p].pdg == genie::kPdgAntiNeutron ) nian++;
        else if(fsps[p].pdg == genie::kPdgPiP         ) nipip++;
        else if(fsps[p].pdg == genie::kPdgPiM         ) nipim++;
        else if(fsps[p].pdg == genie::kPdgPi0         ) nipi0++;
        else if(fsps[p].pdg == genie::kPdgKP          ) nikp++;
        else if(fsps[p].pdg == genie::kPdgKM          ) nikm++;
        else if(fsps[p].pdg == genie::kPdgK0          ) nik0++;
        else if(fsps[p].pdg == genie::kPdgAntiK0      ) niak0++;
        else if(fsps[p].pdg == genie::kPdgGamma       ) nig++;
        else if(fsps[p].pdg == genie::kPdgElectron    ) niem++;
        else if(fsps[p].pdg == genie::kPdgPositron    ) niep++;
        else if(fsps[p].pdg == genie::kPdgMuon        ) nimu++;
        else if(fsps[p].pdg == genie::kPdgAntiMuon    ) nimup++;
        else if(fsps[p].pdg == genie::kPdgAntiMuon    ) nimup++;
        else if(fsps[p].pdg == 14 || fsps[p].pdg == -14 || fsps[p].pdg == 12 || fsps[p].pdg == -12) ninu++;
        else if(fsps[p].pdg > 1e4) ninuc++;
        else niother++;
      }

      if(nuc_pdg > 5000 && (nipip+nipim+nipi0 == 1 || nipip+nipim+nipi0 == 2)){
        LOG("valor_sbn",pERROR) << "The number of pions in this event is: " << nipip+nipim+nipi0 << " and the struck nucleon is neither a proton nor a neutron: " << nuc_pdg;
        util::Die(-1);
      }

      // Construct the reaction mode
      conv::ReacType_t reac = conv::kmtUndefined;
      if (is_dirt){
        reac = conv::kmtDirt; 

        if(is_cosmic){                                                                                                                                                   
          reac = conv::kmtCosmic;                                                                                                                                        
        }                                                                                                                                                                
      }   
      else if (iscc) {
        if      (larsoft_mode==kLArQE                ) { reac = conv::kmtCCQE;    }
        else if (larsoft_mode==kLArMEC               ) { reac = conv::kmtCCMEC;   }
        else if (larsoft_mode==kLArCoh               ) { reac = conv::kmtCCcoh;   }
        else if (larsoft_mode==kLArIMDAnnihilation   ) { reac = conv::kmtCCELScat;}
        else if (larsoft_mode==kLArElectronScattering) { reac = conv::kmtCCELScat;}

        else { // start counting pions
          if      (nikp || nik0 || nikm           )  { reac = conv::kmtCCoth;}
          else if (nipip+nipim == 1 && nipi0 == 0 )  { reac = conv::kmtCC1piC;}
          else if (nipip+nipim == 0 && nipi0 == 1 )  { reac = conv::kmtCC1pi0;}
          else if (nipip+nipim == 2 && nipi0 == 0 )  { reac = conv::kmtCC2piC;}
          else if (nipip+nipim == 0 && nipi0 == 2 )  { reac = conv::kmtCC2pi0;}
          else if (nipip+nipim == 1 && nipi0 == 1 )  { reac = conv::kmtCCpiCpi0;}
          else                                       { reac = conv::kmtCCoth;}
        }
      }
      else if (isnc) {
        if      (larsoft_mode==kLArQE                ) { reac = conv::kmtNCEL;    }
        else if (larsoft_mode==kLArMEC               ) { reac = conv::kmtNCMEC;   }
        else if (larsoft_mode==kLArCoh               ) { reac = conv::kmtNCcoh;   }
        else if (larsoft_mode==kLArElectronScattering) { reac = conv::kmtNCELScat;}

        else { // start counting pions
          if      (nikp || nik0 || nikm           )  { reac = conv::kmtNCoth;}
          else if (nipip+nipim == 1 && nipi0 == 0 )  { reac = conv::kmtNC1piC;}
          else if (nipip+nipim == 0 && nipi0 == 1 )  { reac = conv::kmtNC1pi0;}
          else if (nipip+nipim == 2 && nipi0 == 0 )  { reac = conv::kmtNC2piC;}
          else if (nipip+nipim == 0 && nipi0 == 2 )  { reac = conv::kmtNC2pi0;}
          else if (nipip+nipim == 1 && nipi0 == 1 )  { reac = conv::kmtNCpiCpi0;}
          else                                       { reac = conv::kmtNCoth;}
        }
      }
      else if (larsoft_mode == kLArElectronScattering) {
        reac = conv::kmtELScat;
      }
      
      else if(is_cosmic){ reac = conv::kmtCosmic; 
      }
      else reac = conv::kmtOther;

      conv::NuFlavour_t final_flav = valor::util::FlavourFromPDG(nu_pdg);
      conv::NuFlavour_t init_flav;

      // If we are looking at the selected nuE samples, get the initial neutrino pdg from the file
      init_pdg  = truth.neutrino.initpdg;
      init_flav = valor::util::FlavourFromPDG(init_pdg);

      // Get the real valor mode - this will only be defined for modes turned on right now
      conv::ReacMode_t  reac_mode  = valor::util::ReacTypeToMode(reac,init_flav,final_flav);
      // If we have any events that are in disabled modes, put them into other
      if (reac_mode==conv::kmUndefined) {
        if (isnc) {reac = conv::kmtNCoth;}
        if (iscc) {reac = conv::kmtCCoth;}
        reac_mode  = valor::util::ReacTypeToMode(reac,init_flav,final_flav);
      }

      if(!opt_disable_checks) {
        //
        // ignore events that can not be used in the VALOR analysis
        //
        if(reac_mode == conv::kmUndefined) {
          VLOG("valor_sbn", pWARN) << "- Rejecting event [bad mode] "<<iev;
          VLOG("valor_sbn", pNOTICE)
            <<"\n larsoft_mode = "<<larsoft_mode
            <<"\n baseline = "<< baseline << " m"
            <<"\n iscc="<<iscc<<" and isnc="<<isnc
            <<"\n initpdg="<<init_pdg<<" and finalpdg="<<nu_pdg
            <<"\n valor reaction="<<reac<<" which is "<<util::AsString(reac);
          continue;
        }
        if(isnc && init_flav != final_flav){
          VLOG("valor_sbn", pWARN) << "- Rejecting event [oscillated nc] "<<iev;
          VLOG("valor_sbn", pINFO)
            <<"\n larsoft_mode = "<<larsoft_mode
            <<"\n baseline = "<< baseline << " m"
            <<"\n iscc="<<iscc<<" and isnc="<<isnc
            <<"\n initpdg="<<init_pdg<<" and finalpdg="<<nu_pdg
            <<"\n valor mode="<<reac_mode<<" which is "<<util::AsString(reac_mode);
          continue;
        }
        if(is_nue && (init_pdg == 12 || init_pdg == -12) && (!is_nueint)){
          VLOG("valor_sbn", pWARN) << "- Rejecting event [nue not from intrinsic] "<<iev;
          VLOG("valor_sbn", pINFO)
            <<"\n is_nueint = "<<is_nueint
            <<"\n larsoft_mode = "<<larsoft_mode
            <<"\n baseline = "<< baseline << " m"
            <<"\n iscc="<<iscc<<" and isnc="<<isnc
            <<"\n initpdg="<<init_pdg<<" and finalpdg="<<nu_pdg
            <<"\n valor mode="<<reac_mode<<" which is "<<util::AsString(reac_mode);
          continue;
        }
        if(is_nueosc && init_flav == final_flav){
          VLOG("valor_sbn", pWARN) << "- Rejecting event [non-oscillated event in oscillated file] "<<iev;
          VLOG("valor_sbn", pINFO)
            <<"\n larsoft_mode = "<<larsoft_mode
            <<"\n baseline = "<< baseline << " m"
            <<"\n iscc="<<iscc<<" and isnc="<<isnc
            <<"\n initpdg="<<init_pdg<<" and finalpdg="<<nu_pdg
            <<"\n valor mode="<<reac_mode<<" which is "<<util::AsString(reac_mode);
          continue;
        }
        if(is_numu && init_flav != final_flav){
          VLOG("valor_sbn", pWARN) << "- Rejecting event [oscillated event in non-oscillated file] "<<iev;
          VLOG("valor_sbn", pINFO)
            <<"\n larsoft_mode = "<<larsoft_mode
            <<"\n baseline = "<< baseline << " m"
            <<"\n iscc="<<iscc<<" and isnc="<<isnc
            <<"\n initpdg="<<init_pdg<<" and finalpdg="<<nu_pdg
            <<"\n valor mode="<<reac_mode<<" which is "<<util::AsString(reac_mode);
          continue;
        }
        if(nue_disappearance && is_nueosc){
          VLOG("valor_sbn", pWARN) << "- Rejecting event [appearance oscillated event in nue disappearance template construction] "<<iev;
          VLOG("valor_sbn", pINFO)
            <<"\n larsoft_mode = "<<larsoft_mode
            <<"\n baseline = "<< baseline << " m"
            <<"\n iscc="<<iscc<<" and isnc="<<isnc
            <<"\n initpdg="<<init_pdg<<" and finalpdg="<<nu_pdg
            <<"\n valor mode="<<reac_mode<<" which is "<<util::AsString(reac_mode);
          continue;

        }
      }

      VLOG("valor_sbn",pDEBUG) << "Read event: "
        << "Etrue         = "<< enu_true  << " GeV; "
        << "LArSoft mode  = "<< larsoft_mode  << "; ";
      valor_minidst_reac_mode  = reac_mode;
      valor_minidst_enu_true   = enu_true;
      valor_minidst_baseline   = baseline/1000.;
      valor_minidst_Q2_true    = Q2s;
      valor_minidst_reac_type    = reac;
      valor_minidst_init_nupdg   = init_pdg;
      valor_minidst_final_nupdg  = nu_pdg;
      valor_minidst_target       = 1000180400; // Argon40
      valor_minidst_posx_true  = posx_true;
      valor_minidst_posy_true  = posy_true;
      valor_minidst_posz_true  = posz_true;
      valor_minidst_offaxis_true = offaxis_true;

      // Copy the systematic universes across
      if (write_universes) {
        // Sort out the NonRes parameters
        for (const auto & entry : truth.weights) {
          TString vname = entry.first;
          TString input_vname(vname);
          if(!vname.Contains("NonRes") || !fix_nonres) continue;
          LOG("valor_sbn",pINFO) << " Current input parameter: " << vname;

          const auto & in_vec = entry.second;
          const int nentries = in_vec.size();
          assert(nentries<=max_n_universes);

          // If we are looking at a CC event and an Alt paramater, set all the NC NonRes param entries to be 1
          // Vice-versa for NC !Alt
          if((iscc && input_vname.Contains("Alt")) || (!iscc && !input_vname.Contains("Alt"))){
            for(const TString &new_vname : nonResVariations){
              if((iscc && new_vname.Contains("CC")) || (!iscc && new_vname.Contains("NC"))) continue;
              LOG("valor_sbn",pINFO) << " Setting all NC parameter entries to be 1.0 in a CC event and vice versa";
              for (size_t ii = 0 ; ii < in_vec.size() ; ii++) {
                variation_entries[new_vname][ii] = 1.f;
              }
              variation_nentries[new_vname] = nentries;
            }
          }
          else{ 
            // Get the corresponding variation name if possible
            // Get also the pair of variations this parameter should be split into
            std::pair<TString,TString> new_vnames;
            GetPair(vname,new_vnames);
            std::cout << vname << std::endl;
            std::cout << nuc_pdg << std::endl;
            GetCorrectVariationName(nuc_pdg, iscc, vname);
            if(vname != ""){
              if(vname != new_vnames.first && vname != new_vnames.second){
                LOG("valor_sbn",pERROR) << "Associated parameter name does not match one of the expected names:" 
                  << "\n Parameter : " << vname 
                  << "\n Expected names: " << new_vnames.first << ", " << new_vnames.second;
                util::Die(-1);
              }
              // Now get the other name from the pair
              TString other_name = (vname == new_vnames.first) ? new_vnames.second : new_vnames.first;
              for (size_t ii = 0 ; ii < in_vec.size() ; ii++) {
                LOG("valor_sbn",pINFO) << " New parameter: " << vname << ", other parameter option: " << other_name << " (setting entries to 1.f)";
                // Set the second element of the pair to be 1.f
                variation_entries[other_name][ii] = 1.f;

                // Set the matched name properly
                variation_output_type val = static_cast<variation_output_type>(in_vec[ii]);
                variation_entries[vname][ii] = val;
              }
              variation_nentries[vname] = nentries;
              variation_nentries[other_name] = nentries;
            }
            // Otherwise, get the two parameters it should be split into and set the entries to be 1.0 for this event
            else{
              for (size_t ii = 0 ; ii < in_vec.size() ; ii++) {
                if(abs(in_vec.at(ii)-1)>std::numeric_limits<float>::epsilon()){
                  LOG("valor_sbn",pERROR) << " Setting alternate parameters for " << input_vname << " to be 1 when this parameter isn't always 1, this value: " << in_vec.at(ii);
                  util::Die(-1);
                }
                variation_entries[new_vnames.first][ii]  = 1.f;
                variation_entries[new_vnames.second][ii] = 1.f;
              }
              variation_nentries[new_vnames.first]  = nentries;
              variation_nentries[new_vnames.second] = nentries;
            }
          }
        }

        // For the rest of the parameters, copy the input values to the ouput TTree
        for (const auto & entry : truth.weights) {
          TString vname = entry.first;
          if(vname.Contains("NonRes") && fix_nonres) continue;

          const auto & in_vec = entry.second;
          const int nentries = in_vec.size();
          assert(nentries<=max_n_universes);

          // Now fill the branch
          for (size_t ii = 0 ; ii < in_vec.size() ; ii++) {
            variation_output_type val = static_cast<variation_output_type>(in_vec[ii]);
            variation_output_type & destination = variation_entries[vname][ii];
            if (val < 0.0) {
              LOG("valor_sbn",pINFO) << "Setting negative event weight to zero: "
                <<vname<<"\tEvent="<<iev<<" ; Universe="<<ii
                <<" ; wght = "<<std::scientific << std::setprecision(6)<<val;
              val = 0.0;
            }
            destination = val;
          }
          variation_nentries[vname] = nentries;
        }
      }
    }// !no-truth

    // If we are finding everything 'other' than the sample defined, label the DST sample 'other'
    if(other_sample)
      valor_minidst_sample = GetOtherSample(sample,is_numu);
    else
      valor_minidst_sample     = sample;
    valor_minidst_beam       = beam;
    valor_minidst_wght_pot   = 1./mcsample_pot;
    //valor_minidst_y_reco    = 0.5; // FORCE y=0.5 for now

    VLOG("valor_sbn", pINFO)
      << "\n Information copied to VALOR/SBN Mini-DST: "
      << "\n Beam          = " << valor_minidst_beam
      << "\n Sample        = " << valor_minidst_sample
      << "\n Reaction mode = " << valor_minidst_reac_mode
      << "\n Target        = " << valor_minidst_target
      << "\n Baseline      = " << valor_minidst_baseline << " m"
      << "\n Q2 (true)     = " << valor_minidst_Q2_true << " GeV^2"
      << "\n Ev (true)     = " << valor_minidst_enu_true << " GeV"
      << "\n y (reco)      = " << valor_minidst_y_reco
      << "\n POT weight    = " << valor_minidst_wght_pot
      << "\n Vtx x (true)  = " << valor_minidst_posx_true << " cm"
      << "\n Vtx y (true)  = " << valor_minidst_posy_true << " cm"
      << "\n Vtx z (true)  = " << valor_minidst_posz_true << " cm"
      << "\n Off-axis distance (true)  = " << valor_minidst_offaxis_true << " deg";


    enu_reco    = reco.reco_energy;
    q0_reco     = reco.lepton_energy;
    Theta_reco  = reco.lepton_costh;
    reco_weight = reco.weight;

    VLOG("valor_sbn",pDEBUG) << "Read event: "
      << "Ereco         = "<< enu_reco  << " GeV; "
      << "Off-axis         = "<< valor_minidst_offaxis_true << " deg; "
      << "\n Reco lepton E     = "<< q0_reco  << " GeV; "
      << "\n Reco lepton angle = "<< Theta_reco;
    if(!opt_disable_checks) {

      //
      // ignore events that can not be used in the VALOR analysis
      //
      if(enu_reco <= 0) {
        VLOG("valor_sbn", pINFO) << "- Rejecting event [Ev (reco) <= 0]: "<<enu_reco;
        continue;
      }
    }


    //Check if the event falls into the correct off axis angle range
    if (valor_minidst_offaxis_true<off_axis_bin_min || valor_minidst_offaxis_true>off_axis_bin_max){
      VLOG("valor_sbn", pINFO) << "- Rejecting event [Off-axis angle (true) out of bounds]: "<<valor_minidst_offaxis_true << " deg";
        continue;
      }

    VLOG("valor_sbn", pINFO) << "+ Accepting event";

    valor_minidst_enu_reco    = enu_reco;
    valor_minidst_q0_reco     = q0_reco;
    valor_minidst_Theta_reco  = Theta_reco;
    if(apply_weights)
      valor_minidst_wght = reco_weight;
    else
      valor_minidst_wght = 1.;

    valor_minidst->Fill();

    nev_accept += 1;
    if(!no_truth)
      nev_accept_by_mode[valor_minidst_reac_mode] += 1;
  } // event loop

  if(!no_truth){
    for (unsigned int im = 0 ; im < nev_accept_by_mode.size(); im++) {
      valor::conv::ReacMode_t m = (valor::conv::ReacMode_t)im;
      VLOG("valor_sbn", pNOTICE) << "Selected "<<nev_accept_by_mode[im]<<" events of mode "<<valor::util::AsString(m);
    }
  }

  double f = 100. * (double)nev_accept / (double)nev;
  VLOG("valor_sbn", pNOTICE) <<
    "Used " << f << "% (" << nev_accept << "/" << nev
    << ") of input events.";

  //~ VLOG("valor_sbn", pNOTICE) <<"Total POT read: "<<total_pot_read;

  // Saving VALOR mini-FST to output file
  VLOG("valor_sbn", pNOTICE)
    << "Saving VALOR/SBN mini-DST to: " << opt_output_name;
  valor_minidst->Write();
  fout.Close();

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

bool CheckSample(const event::RecoInteraction &reco, const bool &numu, const valor::conv::Sample_t &sample, const bool &other){
  if(numu){
    if(sample == valor::conv::ksNuMuCC || sample == valor::conv::ksNuMuLikeCChigh|| 
       sample == valor::conv::ksNuMuCCoa00 || sample == valor::conv::ksNuMuCCoa01 || 
       sample == valor::conv::ksNuMuCCoa02 || sample == valor::conv::ksNuMuCCoa03 || 
       sample == valor::conv::ksNuMuCCoa04 || sample == valor::conv::ksNuMuCCoa05 || 
       sample == valor::conv::ksNuMuCCoa06 || sample == valor::conv::ksNuMuCCoa07) return true;
    else if(sample == valor::conv::ksNuMuCC0piNtrk){
      if((reco.npion+reco.npi0) == 0)
        return ((other) ? false : true);
      else return ((other) ? true : false);
    }
    else if(sample == valor::conv::ksNuMuCC0pi1trk){
      if((reco.npion+reco.npi0) == 0 && reco.nproton == 0)
        return ((other) ? false : true);
      else return ((other) ? true : false);
    }
    else if(sample == valor::conv::ksNuMuCC0pi2trk){
      if((reco.npion+reco.npi0) == 0 && reco.nproton == 1)
        return ((other) ? false : true);
      else return ((other) ? true : false);
    }
    else if(sample == valor::conv::ksNuMuCC1piC){
      if(reco.npion == 1 && reco.npi0 == 0)
        return ((other) ? false : true);
      else return ((other) ? true : false);
    }
    else if(sample == valor::conv::ksNuMuCC1pi0){
      if(reco.npion == 0 && reco.npi0 == 1)
        return ((other) ? false : true);
      else return ((other) ? true : false);
    }
    else if(sample == valor::conv::ksNuMuCC1piC1pi0){
      if(reco.npion == 1 && reco.npi0 == 1)
        return ((other) ? false : true);
      else return ((other) ? true : false);
    }
    else if(sample == valor::conv::ksNuMuCC2piC){
      if(reco.npion == 2 && reco.npi0 == 0)
        return ((other) ? false : true);
      else return ((other) ? true : false);
    }
    else if(sample == valor::conv::ksNuMuCCoth){
      if((reco.npion + reco.npi0) >= 2)
        return ((other) ? false : true);
      else return ((other) ? true : false);
    }
    else{
      LOG("valor_sbn",pFATAL) << "Error: Desired sample " << valor::util::AsString(sample) << " is not distinguishable from the information in this set of events ";
      return false;
    }
  }
  if(!numu){
    if(sample == valor::conv::ksNuECC || sample == valor::conv::ksNuELikeCChigh|| 
       sample == valor::conv::ksNuECCoa00 || sample == valor::conv::ksNuECCoa01 || 
       sample == valor::conv::ksNuECCoa02 || sample == valor::conv::ksNuECCoa03 || 
       sample == valor::conv::ksNuECCoa04 || sample == valor::conv::ksNuECCoa05 || 
       sample == valor::conv::ksNuECCoa06 || sample == valor::conv::ksNuECCoa07) return true;
    if(sample == valor::conv::ksNuECC || sample == valor::conv::ksNuELikeCChigh) return true;
    else if(sample == valor::conv::ksNuECC0pi){
      if((reco.npion+reco.npi0) == 0)
        return ((other) ? false : true);
      else return ((other) ? true : false);
    }
    else if(sample == valor::conv::ksNuECC1piC){
      if(reco.npion == 1 && reco.npi0 == 0)
        return ((other) ? false : true);
      else return ((other) ? true : false);
    }
    else if(sample == valor::conv::ksNuECCoth){
      if((reco.npion+reco.npi0) >= 2)
        return ((other) ? false : true);
      else return ((other) ? true : false);
    }
    else{
      LOG("valor_sbn",pFATAL) << "Error: Desired sample " << valor::util::AsString(sample) << " is not distinguishable from the information in this set of events ";
      return false;
    }
  }
  return false; 
}
valor::conv::Sample_t GetOtherSample(const valor::conv::Sample_t &s, const bool &numu){

  // Numu CC
  if(valor::util::IsCC(s) && numu){
    valor::conv::Sample_t new_s = valor::util::SampleFromString("numuCCoth");
    assert(new_s != valor::conv::ksUndefined);
    return new_s;
  }

  // Numu NC
  else if(!valor::util::IsCC(s) && numu){
    valor::conv::Sample_t new_s = valor::util::SampleFromString("numuNCoth");
    assert(new_s != valor::conv::ksUndefined);
    return new_s;
  }

  // Nue CC
  else if(valor::util::IsCC(s) && !numu){
    valor::conv::Sample_t new_s = valor::util::SampleFromString("nueCCoth");
    assert(new_s != valor::conv::ksUndefined);
    return new_s;
  }
  
  // Nue NC
  else if(!valor::util::IsCC(s) && !numu){
    valor::conv::Sample_t new_s = valor::util::SampleFromString("nueNCoth");
    assert(new_s != valor::conv::ksUndefined);
    return new_s;
  }
  else{
    LOG("valor_sbn",pNOTICE) << "The 'other' version of this sample cannot be determined, returning the input, " << valor::util::AsString(s);
    return s;
  }
}
void FillVariationBranches(const TString &vname,
                           const size_t &n_universes,
                           TObjArray &variation_names,
                           std::map<TString,int> &variation_nentries,
                           std::map<TString,variation_output_type *> &variation_entries,
                           TTree *valor_minidst){

  // Fill the variation names object for writing
  variation_names.AddLast(new TObjString(vname.Data()));

  // Make a branch to store the number of universes
  std::string branch_name_n (Form("variation_n_%s"  ,vname.Data()));
  std::string branch_name_nI(Form("variation_n_%s/I",vname.Data()));
  variation_nentries[vname] = -999;
  valor_minidst->Branch(branch_name_n.c_str(),&variation_nentries[vname],branch_name_nI.c_str());
  LOG("valor_sbn",pNOTICE) << "Branch constructed with name: " << branch_name_n;

  assert(n_universes <= max_n_universes);

  // and a branch to store them all in
  variation_entries[vname] = new variation_output_type[max_n_universes];
  for (size_t ii = 0 ; ii < n_universes ; ii++) {
    variation_entries[vname][ii] = static_cast<variation_output_type>(valor::conv::kBadResult);
  }
  std::string branch_name (Form("variation_%s"  ,vname.Data()));
  std::string branch_nameF(Form("variation_%s[%s]/%c",vname.Data(),branch_name_n.c_str(),variation_output_type_root));
  LOG("valor_sbn",pNOTICE) << "Branch constructed with name: " << branch_name;
  valor_minidst->Branch(branch_name.c_str(),variation_entries[vname],branch_nameF.c_str());
}
void GetCorrectVariationName(const int nuc_pdg, const bool iscc, TString &vname){

  TString name("");
  if(iscc){
    if(!vname.Contains("Alt")){
      if(!vname.Contains("bar")){
        if(vname.Contains("p1pi")){
          if(nuc_pdg == 2212)
            name = "genie_NonResRvpCC1pi_Genie";
          else if(nuc_pdg == 2112)
            name = "genie_NonResRvbarnCC1pi_Genie";
        }
        else if(vname.Contains("p2pi")){
          if(nuc_pdg == 2212)
            name = "genie_NonResRvpCC2pi_Genie";
          else if(nuc_pdg == 2112)
            name = "genie_NonResRvbarnCC2pi_Genie";
        }
      }
      else if(vname.Contains("bar")){
        if(vname.Contains("p1pi")){
          if(nuc_pdg == 2212)
            name = "genie_NonResRvbarpCC1pi_Genie";
          else if(nuc_pdg == 2112)
            name = "genie_NonResRvnCC1pi_Genie";
        }
        else if(vname.Contains("p2pi")){
          if(nuc_pdg == 2212)
            name = "genie_NonResRvbarpCC2pi_Genie";
          else if(nuc_pdg == 2112)
            name = "genie_NonResRvnCC2pi_Genie";
        }
      }
    }
  } 
  else{
    if(vname.Contains("Alt")){
      if(!vname.Contains("bar")){
        if(vname.Contains("p1pi")){
          if(nuc_pdg == 2212)
            name = "genie_NonResRvpNC1pi_Genie";
          else if(nuc_pdg == 2112)
            name = "genie_NonResRvbarnNC1pi_Genie";
        }
        else if(vname.Contains("p2pi")){
          if(nuc_pdg == 2212)
            name = "genie_NonResRvpNC2pi_Genie";
          else if(nuc_pdg == 2112)
            name = "genie_NonResRvbarnNC2pi_Genie";
        }
      }
      else if(vname.Contains("bar")){
        if(vname.Contains("p1pi")){
          if(nuc_pdg == 2212)
            name = "genie_NonResRvbarpNC1pi_Genie";
          else if(nuc_pdg == 2112)
            name = "genie_NonResRvnNC1pi_Genie";
        }
        else if(vname.Contains("p2pi")){
          if(nuc_pdg == 2212)
            name = "genie_NonResRvbarpNC2pi_Genie";
          else if(nuc_pdg == 2112)
            name = "genie_NonResRvnNC2pi_Genie";
        }
      }
    }
  }
  if(name == ""){
    LOG("valor_sbn",pINFO) << " No alternative variation name was found for this event in which:" 
      << " \n  CC:            " << iscc
      << " \n  Nucleon PDG:   " << nuc_pdg
      << " \n  Original name: " << vname
      << " \n  --> Will set variation = 1 for all associated NonRes parameters ";
    vname = name;
  }
  else{
    LOG("valor_sbn",pINFO) << " Alternative variation name for: " << vname << " is: " << name << " where: "
      << " \n  CC:            " << iscc
      << " \n  Nucleon PDG:   " << nuc_pdg;

    vname = name;
  }
}
void GetPair(const TString &vname, std::pair<TString,TString> &new_vnames){
  std::pair<TString,TString> names;
  if(vname == "genie_NonResRvp1pi_Genie"){
    names = std::make_pair("genie_NonResRvpCC1pi_Genie","genie_NonResRvbarnCC1pi_Genie");
  }
  else if(vname == "genie_NonResRvbarp1pi_Genie"){
    names = std::make_pair("genie_NonResRvbarpCC1pi_Genie","genie_NonResRvnCC1pi_Genie");
  }
  else if(vname == "genie_NonResRvp2pi_Genie"){
    names = std::make_pair("genie_NonResRvpCC2pi_Genie","genie_NonResRvbarnCC2pi_Genie");
  }
  else if(vname == "genie_NonResRvbarp2pi_Genie"){
    names = std::make_pair("genie_NonResRvbarpCC2pi_Genie","genie_NonResRvnCC2pi_Genie");
  }
  else if(vname == "genie_NonResRvp1piAlt_Genie"){
    names = std::make_pair("genie_NonResRvpNC1pi_Genie","genie_NonResRvbarnNC1pi_Genie");
  }
  else if(vname == "genie_NonResRvbarp1piAlt_Genie"){
    names = std::make_pair("genie_NonResRvbarpNC1pi_Genie","genie_NonResRvnNC1pi_Genie");
  }
  else if(vname == "genie_NonResRvp2piAlt_Genie"){
    names = std::make_pair("genie_NonResRvpNC2pi_Genie","genie_NonResRvbarnNC2pi_Genie");
  }
  else if(vname == "genie_NonResRvbarp2piAlt_Genie"){
    names = std::make_pair("genie_NonResRvbarpNC2pi_Genie","genie_NonResRvnNC2pi_Genie");
  }
  else{
    LOG("valor_sbn",pERROR) << "No alternative parameter set found for input variation: " << vname;
    valor::util::Die(-1);
  }
  new_vnames = names;
}
