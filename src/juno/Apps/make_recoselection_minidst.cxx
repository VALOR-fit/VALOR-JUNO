//____________________________________________________________________________
/*!

  \app       valor_juno-make_recoselection_minidst

  \brief     Converts input recoselection analysis DSTs to the VALOR mini-DST format.

  \syntax    valor_juno-make_recoselection_minidst
                --input  input_analysis_file.root
                --output output_file_containing_valor_minidst.root
                --format format_of_the_input_file <hpgartpc, lartpc, fgt>
               [--disable-checks]

  \author   Yaoqi Cao <yaoqi.cao \at warwick.ac.uk>
            University of Warwick

  \created  Jan 2014

  \copyright Copyright (c) 2010-2024, The VALOR group. All rights reserved.

*/
//____________________________________________________________________________

// System
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

using std::string;
using std::ostringstream;

void GetCommandLineArgs (int /*argc*/, char ** /*argv*/);
static const string gTreeName = "valor_tree";
static const int gPrintEveryNEvents   = 50000;

typedef enum ScatterCode{

  kScatUnknownInteraction = -1,
  kNULL                   =  0,
  kScatQE                 =  1,
  kScatSingleKaon         =  2,
  kScatDIS                =  3,
  kScatRes                =  4,
  kScatCoh                =  5,
  kScatDiffractive        =  6,
  kScatElectronScattering =  7,
  kScatInverseMuonDecay   =  8,
  kScatAMNuGamma          =  9,
  kScatMEC                = 10,
  kScatCohElastic         = 11,
  kScatInverseBetaDecay   = 12,
  kScatGlashowResonance   = 13,
  kScatIMDAnnihilation    = 14,

} ScatterCode_t;

typedef enum NuMuCCExclusiveTopology{

  kUnknownTopology              = -1,
  kCC0Pi                        =  0,
  kCC1Pi                        =  1,
  kCC2Pi                        =  2,
  kCCPi0                        =  3,
  kCCOther                      =  4

} ExclusiveTopology_t;

int main(int argc, char ** argv)
{
  // Set GENIE print level
  genie::GHepRecord::SetPrintLevel(genie::RunOpt::Instance()->EventRecordPrintLevel());

  using namespace valor;

  CmdLnArgList arg;
  arg.Read(argc,argv, "$VALOR_EXPT/config/apps/make_recoselection_minidst.arg");
  LOG("valor_juno", pNOTICE) << arg;

  std::string thresholds = util::CompleteConfigPath(arg.GetString("mesg-thresholds"));
  Messenger::Instance()->Configure(thresholds);

  string opt_output_name       = arg.GetString("output");
  string opt_input_name        = arg.GetString("input");
  bool   opt_disable_checks    = arg.GetBool("disable-checks");

  int    valor_minidst_reac_mode      = -1;  // see Conventions/ReacMode.h
  int    valor_minidst_target         = -1;  // nuclear target -  GENIE PDG convention
  double valor_minidst_enu_true       = -1.; // true neutrino energy, in GeV
  double valor_minidst_Q2_true        = -1.; // true neutrino Q2
  double valor_minidst_enu_reco       = -1.; // reco neutrino energy, in GeV
  double valor_minidst_mu_momentum    = -1.; // reco muon momentum
  double valor_minidst_mu_cos_z       = -1.; // reco muon cos(angle to z)
  double valor_minidst_y_reco         = -1.; // inelasticity
  double valor_minidst_wght           = -1.; // misc weight
  double valor_minidst_wght_smear     = -1.; // energy smearing weight
  double valor_minidst_reco_exclusive = -1.; // exclusive reconstructed topology
  double valor_minidst_true_exclusive = -1.; // exclusive true topology

  TFile fout(opt_output_name.c_str(),"recreate");
  fout.cd();

  TTree * valor_minidst =
  new TTree("valor_juno_minidst","Mini-DST input for VALOR/JUNO fit");

  valor_minidst->Branch("reac_mode"     , &valor_minidst_reac_mode     ,"reac_mode/I"      );
  valor_minidst->Branch("target"        , &valor_minidst_target        ,"target/I"         );
  valor_minidst->Branch("enu_true"      , &valor_minidst_enu_true      ,"enu_true/D"       );
  valor_minidst->Branch("Q2_true"       , &valor_minidst_Q2_true       ,"q2_true/D"        );
  valor_minidst->Branch("enu_reco"      , &valor_minidst_enu_reco      ,"enu_reco/D"       );
  valor_minidst->Branch("enu_reco"      , &valor_minidst_enu_reco      ,"enu_reco/D"       );
  valor_minidst->Branch("mu_momentum"   , &valor_minidst_mu_momentum   ,"mu_momentum/D"    );
  valor_minidst->Branch("mu_cos_z"      , &valor_minidst_mu_cos_z      ,"mu_cos_z/D"       );
  valor_minidst->Branch("y_reco"        , &valor_minidst_y_reco        ,"y_reco/D"         );
  valor_minidst->Branch("wght"          , &valor_minidst_wght          ,"wght/D"           );
  valor_minidst->Branch("wght_smear"    , &valor_minidst_wght_smear    ,"wght_smear/D"     );
  valor_minidst->Branch("reco_exclusive", &valor_minidst_reco_exclusive,"reco_exclusive/D" );
  valor_minidst->Branch("true_exclusive", &valor_minidst_true_exclusive,"true_exclusive/D" );

   // open the ROOT file and get the TTree & its header
  TChain * analysis_tree = new TChain(gTreeName.c_str());
  analysis_tree->Add(opt_input_name.c_str());
  if(!analysis_tree || !analysis_tree->GetEntries()) {
    VLOG("valor_juno", pFATAL)
      << "Unable to read " << gTreeName
      <<" from "<< opt_input_name;
    exit(1);
  }

  Long64_t nev = analysis_tree->GetEntries();
  VLOG("valor_juno", pNOTICE)
    << "Input analysis tree has " << nev << ((nev==1) ? " entry." : " entries.");
  analysis_tree->GetEntry(0);


  // Variables to get from analysis tree
  bool iscc           = false;
  bool isnc           = false;
  double enu_true     = -999.;
  double enu_reco     = -999.;
  double qsqr         = -999.;
  double mu_momentum  = -999.;
  double mu_cos_z     = -999.;
  int nu_pdg          = -1;
  int nkaons          = -1;
  int npip            = -1;
  int npim            = -1;
  int npi0            = -1;
  int true_exclusive  = -1;
  int reco_exclusive  = -1;

  analysis_tree->SetBranchAddress("iscc",           &iscc            );
  analysis_tree->SetBranchAddress("isnc",           &isnc            );
  analysis_tree->SetBranchAddress("nu_pdg",         &nu_pdg          );
  analysis_tree->SetBranchAddress("enu_true",       &enu_true        );
  analysis_tree->SetBranchAddress("enu_reco",       &enu_reco        );
  analysis_tree->SetBranchAddress("mu_momentum",    &mu_momentum     );
  analysis_tree->SetBranchAddress("mu_cos_z",       &mu_cos_z        );
  analysis_tree->SetBranchAddress("qsqr",           &qsqr            );
  analysis_tree->SetBranchAddress("nkaons",         &nkaons          );
  analysis_tree->SetBranchAddress("npip",           &npip            );
  analysis_tree->SetBranchAddress("npim",           &npim            );
  analysis_tree->SetBranchAddress("npi0",           &npi0            );
  analysis_tree->SetBranchAddress("true_topology",  &true_exclusive  );
  analysis_tree->SetBranchAddress("reco_topology",  &reco_exclusive  );

  // event loop
  //

  std::vector<int> nev_accept_by_mode;
  for (int im = 0 ; im < valor::conv::kNReacModes; im++) {
    nev_accept_by_mode.push_back(0);
  }

  for(Long64_t iev = 0; iev < nev; iev++) {

    if (gPrintEveryNEvents && !(iev%gPrintEveryNEvents)) {
      VLOG("valor_juno", pNOTICE) << "** Reading event: " << iev;
    }
    else {
      VLOG("valor_juno", pDEBUG) << "** Reading event: " << iev;
    }

    analysis_tree->GetEntry(iev);
    conv::ReacType_t reac = conv::kmtUndefined; 
  if (iscc) { 
        if      (nkaons )                      {reac = conv::kmtCCoth;}
        else if (npip+npim == 1 && npi0 == 0)  {reac = conv::kmtCC1piC;}
        else if (npip+npim == 0 && npi0 == 1)  {reac = conv::kmtCC1pi0;}
        else if (npip+npim == 2 && npi0 == 0)  {reac = conv::kmtCC2piC;}
        else if (npip+npim == 0 && npi0 == 2)  {reac = conv::kmtCC2pi0;}
        else if (npip+npim == 1 && npi0 == 1)  {reac = conv::kmtCCpiCpi0;}
        else                                   {reac = conv::kmtCCoth;}
  }
    else if (isnc) { // start counting pions
        if      (nkaons)                       {reac = conv::kmtNCoth;}
        else if (npip+npim == 1 && npi0 == 0)  {reac = conv::kmtNC1piC;}
        else if (npip+npim == 0 && npi0 == 1)  {reac = conv::kmtNC1pi0;}
        else if (npip+npim == 2 && npi0 == 0)  {reac = conv::kmtNC2piC;}
        else if (npip+npim == 0 && npi0 == 2)  {reac = conv::kmtNC2pi0;}
        else if (npip+npim == 1 && npi0 == 1)  {reac = conv::kmtNCpiCpi0;}
        else                                   {reac = conv::kmtNCoth;}
    }
  else {
    reac = conv::kmtOther;
  }



    conv::NuFlavour_t final_flav = valor::util::FlavourFromPDG(nu_pdg);
    conv::NuFlavour_t init_flav = final_flav;

    // Get the real valor mode
    conv::ReacMode_t  reac_mode  = valor::util::ReacTypeToMode(reac,init_flav,final_flav);

    if(!opt_disable_checks) {
      //
      // ignore events that can not be used in the VALOR analysis
      //
      if(reac_mode == conv::kmUndefined) {
        VLOG("valor_juno", pWARN) << "- Rejecting event [bad mode]\n";
        continue;
      }
    }

    VLOG("valor_juno",pDEBUG) << "Read event: "
      << "Etrue         = "<< enu_true  << " GeV; ";

    valor_minidst_reac_mode      = reac_mode;
    valor_minidst_enu_true       = enu_true;
    valor_minidst_wght           = 1.;
    valor_minidst_Q2_true        = qsqr;
    valor_minidst_target         = 1000180400; // Argon40
    valor_minidst_y_reco         = 0.5; // FORCE y=0.5 for now
    valor_minidst_mu_momentum    = mu_momentum;
    valor_minidst_mu_cos_z       = mu_cos_z;
    valor_minidst_true_exclusive = true_exclusive;
    valor_minidst_reco_exclusive = reco_exclusive;

    VLOG("valor_juno", pINFO)
      << "\n Information copied to VALOR JUNO recoselection Mini-DST: "
      << "\n Reaction mode = " << valor_minidst_reac_mode
      << "\n Target        = " << valor_minidst_target
      << "\n Q2 (true)     = " << valor_minidst_Q2_true << " GeV^2"
      << "\n Ev (true)     = " << valor_minidst_enu_true << " GeV"
      << "\n Pmu (reco)    = " << valor_minidst_mu_momentum << "GeV/c^2"
      << "\n P cosZ (reco) = " << valor_minidst_mu_cos_z
      << "\n y (reco)      = " << valor_minidst_y_reco
      << "\n Smear weight  = " << valor_minidst_wght_smear
      << "\n Reco topology = " << valor_minidst_reco_exclusive
      << "\n True topology = " << valor_minidst_true_exclusive;

    VLOG("valor_juno",pDEBUG) << "Read event: "
      << "Ereco         = "<< enu_reco  << " GeV; ";
    if(!opt_disable_checks) {
        // ignore events that can not be used in the VALOR analysis
      //
      if(enu_reco <= 0) {
        VLOG("valor_juno", pINFO) << "- Rejecting event [Ev (reco) <= 0]: "<< enu_reco;
        continue;
      }

      VLOG("valor_juno", pINFO) << "+ Accepting event";
    }

    valor_minidst_enu_reco  = enu_reco;
    VLOG("valor_juno", pNOTICE)
      << "\n Ev (reco)     = " << valor_minidst_enu_reco << " GeV";

    valor_minidst->Fill();
  } // event loop

  std::cout << " Number in valor tree : " << valor_minidst->GetEntries() << std::endl;

  // Saving VALOR mini-FST to output file
  VLOG("valor_juno", pNOTICE)
    << "Saving VALOR/JUNO mini-DST to: " << opt_output_name;
  valor_minidst->Write();
  fout.Close();

  VLOG("valor_juno", pNOTICE)  << "Done!";

  return 0;
}
