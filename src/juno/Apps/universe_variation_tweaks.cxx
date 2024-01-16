//____________________________________________________________________________
/*!

  \app       valor_dune-universe_variation_tweaks

  \brief     Vary the MC spectra by throwing sets of systematics

  \options     --output              : Where to save tweaked histograms to. [def: tweaks_[PARAMS].root]
  --input               : Input DST file [required]
  --params              : Parameters to tweak - to name the output file [def: allsyst]
  --sample              : Topology [def: numulikeCChigh]
  --detector            : Detector [def: sbn_nd]
  --kinevar             : Kinematic paramter to parameterise in [def: ereco]
  --ntoys               : Number of toy-MC experiments [def: 1000]
  --config              : Analysis configuration file [required]
  --mesg-thresholds     : Message thresholds [def: "verbosity.xml"]


  \author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
  University of Liverpool & STFC Rutherford Appleton Laboratory
  Steve Dennis <s.r.dennis \at liverpool.ac.uk>
  University of Liverpool

  \created   October 2016

  \copyright (c) 2010-2016, The VALOR Fitting Group, All Rights Reserved

*/
//____________________________________________________________________________
// cpp
#include <stdlib.h>

// root
#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TString.h"
#include "TSystem.h"
// fmwk
#include "valor/Conventions/Opt.h"
#include "valor/CmdLine/CmdLnArgList.h"
#include "valor/FitCore/CorrelatedSystematics.h"
#include "valor/FitCore/PhysParamList.h"
#include "valor/FitCore/FitParamList.h"
#include "valor/FitCore/FitModelABC.h"
#include "valor/FitCore/FitDataABC.h"
#include "valor/FitCore/Likelihood.h"
#include "valor/FitCore/FitterABC.h"
#include "valor/FitCore/UtlFit.h"
#include "valor/FitCore/CompoundFitElements.h"
#include "valor/FitCore/SingleSampleFitModelABC.h"
#include "valor/FitCore/UtlFitParamList.h"
#include "valor/Utilities/OscParamValues.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Utilities/Env.h"
#include "valor/Utilities/VThreadPool.h"
#include "valor/Utilities/UtlParam.h"
#include "valor/Utilities/SingleSampleKineBinning.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Ntuple/StdFitNtp.h"
#include "valor/FitElements/GenericSingleSampleFitModel.h"
#include "valor/KineBinning/UtlBinning.h"

std::vector<std::string> modern_xsec = 
{ "genie_DISAth_Genie",
  "genie_DISBth_Genie",
  "genie_DISCv1u_Genie",
  "genie_DISCv2u_Genie",
  "genie_IntraNukeNabs_Genie",
  "genie_IntraNukeNcex_Genie",
  "genie_IntraNukeNinel_Genie",
  "genie_IntraNukeNmfp_Genie",
  "genie_IntraNukeNpi_Genie",
  "genie_IntraNukePIabs_Genie",
  "genie_IntraNukePIcex_Genie",
  "genie_IntraNukePIinel_Genie",
  "genie_IntraNukePImfp_Genie",
  "genie_IntraNukePIpi_Genie",
  "genie_ResDecayGamma_Genie",
  "genie_ccresVector_Genie",
  "genie_cohMA_Genie",
  "genie_cohR0_Genie",
  "genie_ncelAxial_Genie",
  "genie_ncelEta_Genie",
  "genie_ncresVector_Genie"
};
std::vector<std::string> proposal_xsec = 
{ "genie_ccresAxial_Genie",
  "genie_ncresAxial_Genie",
  "genie_qema_Genie",
  "genie_NC_Genie",
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
  "genie_NonResRvbarpNC2pi_Genie"
};
std::vector<std::string> OG_proposal_xsec = 
{ "genie_ccresAxial_Genie",
  "genie_ncresAxial_Genie",
  "genie_qema_Genie",
  "genie_NC_Genie",
  "genie_NonResRvbarp1piAlt_Genie",
  "genie_NonResRvbarp1pi_Genie",
  "genie_NonResRvbarp2piAlt_Genie",
  "genie_NonResRvbarp2pi_Genie",
  "genie_NonResRvp1piAlt_Genie",
  "genie_NonResRvp1pi_Genie",
  "genie_NonResRvp2piAlt_Genie",
  "genie_NonResRvp2pi_Genie",
};
std::vector<std::string> flux_uncorr = 
{ "expskin_FluxUnisim",
  "horncurrent_FluxUnisim",
  "nucleoninexsec_FluxUnisim",
  "nucleonqexsec_FluxUnisim",
  "nucleontotxsec_FluxUnisim",
  "pioninexsec_FluxUnisim",
  "pionqexsec_FluxUnisim",
  "piontotxsec_FluxUnisim"
};
std::vector<std::string> flux_corr = 
{ "kminus_PrimaryHadronNormalization",
  "piminus_PrimaryHadronSWCentralSplineVariation",
  "kplus_PrimaryHadronFeynmanScaling",
  "piplus_PrimaryHadronSWCentralSplineVariation",
  "kzero_PrimaryHadronSanfordWang"
};
std::vector<std::string> piplus = 
{  
  "piplus_PrimaryHadronSWCentralSplineVariation"
};

std::vector<std::string> modern_xsec_short = 
{ "DISAth",
  "DISBth",
  "DISCv1u",
  "DISCv2u",
  "IntraNukeNabs",
  "IntraNukeNcex",
  "IntraNukeNinel",
  "IntraNukeNmfp",
  "IntraNukeNpi",
  "IntraNukePIabs",
  "IntraNukePIcex",
  "IntraNukePIinel",
  "IntraNukePImfp",
  "IntraNukePIpi",
  "ResDecayGamma",
  "ccresVector",
  "cohMA",
  "cohR0",
  "ncelAxial",
  "ncelEta",
  "ncresVector"
};
std::vector<std::string> proposal_xsec_short = 
{ "ccresAxial",
  "ncresAxial",
  "qema",
  "NC",
  "NonResRvpCC1pi",
  "NonResRvnCC1pi",
  "NonResRvbarnCC1pi",
  "NonResRvbarpCC1pi",
  "NonResRvpCC2pi",
  "NonResRvnCC2pi",
  "NonResRvbarnCC2pi",
  "NonResRvbarpCC2pi",
  "NonResRvpNC1pi",
  "NonResRvnNC1pi",
  "NonResRvbarnNC1pi",
  "NonResRvbarpNC1pi",
  "NonResRvpNC2pi",
  "NonResRvnNC2pi",
  "NonResRvbarnNC2pi",
  "NonResRvbarpNC2pi"
};
std::vector<std::string> OG_proposal_xsec_short = 
{ "ccresAxial",
  "ncresAxial",
  "qema",
  "NC",
  "NonResRvbarp1piAlt",
  "NonResRvbarp1pi",
  "NonResRvbarp2piAlt",
  "NonResRvbarp2pi",
  "NonResRvp1piAlt",
  "NonResRvp1pi",
  "NonResRvp2piAlt",
  "NonResRvp2pi",
};

std::vector<std::string> flux_uncorr_short = 
{ "expskin",
  "horncurrent",
  "nucleoninexsec",
  "nucleonqexsec",
  "nucleontotxsec",
  "pioninexsec",
  "pionqexsec",
  "piontotxsec"
};
std::vector<std::string> flux_corr_short = 
{ "kminus",
  "piminus",
  "kplus",
  "piplus",
  "kzero"
};
std::vector<std::string> piplus_short = 
{  
  "piplus"
};

double enu_reco_max = 10.;

using namespace std;
using namespace valor;
using namespace valor::conv;

int main(int argc, char ** argv)
{
  time_t scriptstartrawtime;
  struct tm * scriptstarttimeinfo;
  time (&scriptstartrawtime);
  scriptstarttimeinfo = localtime (&scriptstartrawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Script start local time and date:  " << asctime(scriptstarttimeinfo);
  std::cout << "-----------------------------------------------------------" << std::endl;


  valor::util::PrintBanner(Env::Instance()->Banner(),500);
 
  CmdLnArgList arg;
  arg.Read(argc,argv,"$VALOR_EXPT/config/apps/universe_variation_tweaks.arg");
  LOG("universe_variation_tweaks", pNOTICE) << arg;
  
  // Read and iterpret the command line
  int ntoys               = arg.GetInt("ntoys");
  string config_file      = arg.GetString("config");
  string input_rootfile   = arg.GetString("input");
  string output_rootfile  = arg.GetString("output");
  string params_to_tweak  = arg.GetString("params");
  string mesg_thresholds  = arg.GetString("mesg-thresholds");
  string arg_sample       = arg.GetString("sample");
  string arg_detector     = arg.GetString("detector");
  string arg_kine_var     = arg.GetString("kine-var");
  string tree_name        = arg.GetString("tree-name");
  string param_name       = arg.GetString("param-name");

  Messenger::Instance()->Configure(util::CompleteConfigPath(mesg_thresholds));
  LOG("universe_variation_tweaks", pNOTICE) << *Rnd::Instance();
  
  Sample_t   sample       = util::SampleFromString(arg_sample);
  Detector_t detector     = util::DetectorFromString(arg_detector);

  //Define kine space
  conv::KineVar_t   kv(util::KineVarFromString(arg_kine_var));
  conv::KineSpace_t ks(util::KineSpaceFromKineVar(kv));
  assert(ks != conv::kspUndefined);
  if(util::IsTrue(kv))
    LOG("valor_sbn",pNOTICE) << " Kinematic variable " << util::AsString(kv) << " is a true parameter";
  if(util::IsReco(kv))
    LOG("valor_sbn",pNOTICE) << " Kinematic variable " << util::AsString(kv) << " is a reconstructed parameter";

  std::string elem_readopt("nomc nosplines");
  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(util::CompleteConfigPath(config_file),elem_readopt.c_str());
  FitModelABC * model = elements->FitModel();
  CompoundFitModel * cm = dynamic_cast<CompoundFitModel*>(model);
  assert(cm);
  
  GenericSingleSampleFitModel * gmodel = 0;
  bool got_appropriate_model = false;
  for(unsigned int imodel = 0; imodel < cm->GetNSubModels(); imodel++) {
    gmodel = dynamic_cast<GenericSingleSampleFitModel*>(cm->GetSubModel(imodel));
    if  (gmodel->GetDetector() == detector && 
         gmodel->GetSample()   == sample) 
      { got_appropriate_model = true; break; }
  }
  if (!got_appropriate_model) {
    LOG("universe_variation_tweaks", pFATAL) << "Failed to identify appropriate fit model to get binning.";
  }
  LOG("universe_variation_tweaks", pINFO) << " Successfully accessed the relevant model for binning information";
  
  std::vector<std::string> systematics, file_systematics;
  if (params_to_tweak == "proposal_xsec") {
    systematics      = proposal_xsec;
    file_systematics = proposal_xsec_short;
    LOG("universe_variation_tweaks", pNOTICE) << " Tweaking proposal_xsec systematics ";
    for(const std::string &s : file_systematics){
      LOG("universe_variation_tweaks", pINFO) << s;
    }
  }
  else if (params_to_tweak == "OG_proposal_xsec") {
    systematics      = OG_proposal_xsec;
    file_systematics = OG_proposal_xsec_short;
    LOG("universe_variation_tweaks", pNOTICE) << " Tweaking OG_proposal_xsec systematics ";
    for(const std::string &s : file_systematics){
      LOG("universe_variation_tweaks", pINFO) << s;
    }
  }
  else if (params_to_tweak == "modern_xsec") {
    systematics      = modern_xsec;
    file_systematics = modern_xsec_short;
    LOG("universe_variation_tweaks", pNOTICE) << " Tweaking modern_xsec systematics ";
    for(const std::string &s : file_systematics){
      LOG("universe_variation_tweaks", pINFO) << s;
    }
  }
  else if (params_to_tweak == "flux") {
    systematics      = flux_corr;
    systematics.insert(systematics.end(), flux_uncorr.begin(), flux_uncorr.end());
    file_systematics = flux_corr_short;
    file_systematics.insert(file_systematics.end(), flux_uncorr_short.begin(), flux_uncorr_short.end());
    LOG("universe_variation_tweaks", pNOTICE) << " Tweaking all flux systematics ";
    for(const std::string &s : file_systematics){
      LOG("universe_variation_tweaks", pINFO) << s;
    }
  }
  else if (params_to_tweak == "flux_corr") {
    systematics      = flux_corr;
    file_systematics = flux_corr_short;
    LOG("universe_variation_tweaks", pNOTICE) << " Tweaking correlated flux systematics ";
    for(const std::string &s : file_systematics){
      LOG("universe_variation_tweaks", pINFO) << s;
    }
  }
  else if (params_to_tweak == "flux_uncorr") {
    systematics      = flux_uncorr;
    file_systematics = flux_uncorr_short;
    LOG("universe_variation_tweaks", pNOTICE) << " Tweaking uncorrelated flux systematics ";
    for(const std::string &s : file_systematics){
      LOG("universe_variation_tweaks", pINFO) << s;
    }
  }
  else if (params_to_tweak == "pion") {
    systematics      = piplus;
    file_systematics = piplus_short;
    LOG("universe_variation_tweaks", pNOTICE) << " Tweaking piplus flux systematic ";
    for(const std::string &s : file_systematics){
      LOG("universe_variation_tweaks", pINFO) << s;
    }
  }
  else if (params_to_tweak == "all") {
    systematics      = flux_corr;
    file_systematics = flux_corr_short;
    systematics.insert(systematics.end(), flux_uncorr.begin(), flux_uncorr.end());
    systematics.insert(systematics.end(), proposal_xsec.begin(), proposal_xsec.end());
    systematics.insert(systematics.end(), modern_xsec.begin(), modern_xsec.end());
    file_systematics.insert(file_systematics.end(), flux_uncorr_short.begin(), flux_uncorr_short.end());
    file_systematics.insert(file_systematics.end(), proposal_xsec_short.begin(), proposal_xsec_short.end());
    file_systematics.insert(file_systematics.end(), modern_xsec_short.begin(), modern_xsec_short.end());
    LOG("universe_variation_tweaks", pNOTICE) << " Tweaking all systematics ";
    for(const std::string &s : file_systematics){
      LOG("universe_variation_tweaks", pINFO) << s;
    }
  }
  else if (params_to_tweak == "choose") {
    file_systematics.push_back(param_name);
    if (std::find(modern_xsec_short.begin(), modern_xsec_short.end(), param_name) != modern_xsec_short.end() ||
        std::find(proposal_xsec_short.begin(), proposal_xsec_short.end(), param_name) != proposal_xsec_short.end() ||
        std::find(OG_proposal_xsec_short.begin(), OG_proposal_xsec_short.end(), param_name) != OG_proposal_xsec_short.end()){
        systematics.push_back("genie_" + param_name + "_Genie");
    }
    else if (std::find(flux_uncorr_short.begin(), flux_uncorr_short.end(), param_name) != flux_uncorr_short.end()){
        systematics.push_back(param_name + "_FluxUnisim");
    }
    else if (std::find(flux_corr.begin(), flux_corr.end(), param_name) != flux_corr.end()){
        systematics.push_back(param_name);
    }
    else {
        LOG("universe_variation_tweaks",pFATAL) << "Unknown input systematic: "<< param_name; 
        util::Die(-1);
    }
    LOG("universe_variation_tweaks", pNOTICE) << " Tweaking single systematic " << param_name;
  }
  else {
    LOG("universe_variation_tweaks",pFATAL) << "Unknown systematics list: "<< params_to_tweak;
    util::Die(-1);
  }
  
  // Open input analysis file and get the DST
  TFile fin(input_rootfile.c_str(),"READ");
  if (!fin.IsOpen()) {
    LOG("universe_variation_tweaks", pFATAL) 
      << "Failed to open file "<<input_rootfile.c_str();
    util::Die(-1);
  }
  TTree * analysis_tree = dynamic_cast <TTree *> (fin.Get(tree_name.c_str()));
  if(!analysis_tree) {
    LOG("universe_variation_tweaks", pFATAL) 
      << "No event tree " << tree_name 
      << " in " << " input analysis file: "
      << input_rootfile;
    util::Die(-1);
  }
  LOG("universe_variation_tweaks", pINFO) << " Successfully accessed analysis TTree with the name " << tree_name;

  // Get the branches we're interested in
  //
  // Universes
  const size_t max_n_universes = 10000;
  std::map<std::string,int> branch_n_universes;
  std::map<std::string,std::vector<float>> branch_universes;
  for (const string &variation_name : systematics) {
    branch_universes.emplace(variation_name,std::vector<float>(max_n_universes));
    std::string branch_name   = Form("variation_%s",variation_name.c_str());
    std::string branch_name_n = Form("variation_n_%s",variation_name.c_str());
    analysis_tree->SetBranchStatus (branch_name_n.c_str(),1);
    analysis_tree->SetBranchStatus (branch_name  .c_str(),1);
    bool bstat_univ_n = analysis_tree->SetBranchAddress(branch_name_n.c_str(),&branch_n_universes[variation_name]);
    bool bstat_univ   = analysis_tree->SetBranchAddress(branch_name.c_str(),branch_universes[variation_name].data());
    if (bstat_univ_n || bstat_univ) {
      LOG("universe_variation_tweaks",pFATAL) << " Failed to load needed universe branch.";
      util::Die(-1);
    }
  }
  LOG("universe_variation_tweaks", pINFO) << " Successfully loaded universe branches";

  // Only apply simple binning for now
  const SingleSampleKineBinning & binning = *(gmodel->GetBinning());
  const VALORBinning * bins_sb            = util::IsTrue(kv) ? (binning.TrueSmpBinningSingle()) : (binning.RecoSmpBinning());
  if(bins_sb == nullptr && util::IsTrue(kv)){
    LOG("valor_sbn",pFATAL) << "Unable to determine the true binning, prehaps you have defined them per mode instead of using 'default'";
    util::Die(-1);
  }
  const int n_bins = bins_sb->NBins();
  LOG("valor_sbn",pINFO) << " Number of bins: " << n_bins;

  // Setup arrays to be used for constructing the histograms
  VALORArray *pdfnom = new VALORArray(n_bins);
  pdfnom->SetAll(0.);
  std::vector<VALORArray*> pdfuni;
  for(int itoy = 0; itoy < ntoys; ++itoy){
    pdfuni.push_back(new VALORArray(n_bins));
    pdfuni.at(itoy)->SetAll(0.);
  }
  LOG("universe_variation_tweaks", pINFO) << " Successfully setup PDFs";

  // Now get the map of universe id for each throw and parameter
  std::map< std::string, std::vector<int> > syst_toy_universe;
  analysis_tree->GetEntry(0);
  for(const std::string s : systematics){
    std::vector<int> universe_ids;
    LOG("universe_variation_tweaks",pDEBUG) << " For " << s << " with " << ntoys << " throws and 1000 possible universes";
    if(ntoys > branch_n_universes.at(s)){
      LOG("valor_sbn",pFATAL) << " Too many toys: " << ntoys << ", there are only " << branch_n_universes.at(s) << " entries for this parameter" << s;
      util::Die(-1);
    }
    for(int itoy = 0; itoy < ntoys; ++itoy){
      LOG("universe_variation_tweaks",pDEBUG) << " Current systematic has " << branch_n_universes.at(s) << " so will generate a random number from that range";
      int iuni = rand() % branch_n_universes.at(s);
      LOG("universe_variation_tweaks",pDEBUG) << " Throw " << itoy << " will access universe " << iuni;
      universe_ids.push_back(iuni);
    }
    syst_toy_universe.emplace(s,universe_ids);
  }

  // Ereco
  double enu_reco   = -999.;
  double enu_true   = -999.;
  double weight     = -999.;
  double pot_weight = -999.;
  analysis_tree->SetBranchAddress("enu_reco", &enu_reco);
  analysis_tree->SetBranchAddress("enu_true", &enu_true);
  analysis_tree->SetBranchAddress("wght",     &weight);
  analysis_tree->SetBranchAddress("wght_pot", &pot_weight);

  double data_pot = gmodel->GetPOT();
  int nev = analysis_tree->GetEntries();
  
  LOG("universe_variation_tweaks", pINFO) << " About to enter analysis tree loop with " << nev << " entries";


  unsigned int accepted_events = 0;
  unsigned int skipped_events  = 0;
  unsigned int one_varied = 0;

  // Vector to hold the IDs of any universe which contains events with no impact
  std::vector<int> skipped_universe_ids;

  // Histograms to assess the weights
  TH1D *h_weight_product         = new TH1D("h_weight_product","",200,-1,3);
  TH1D *h_weight_product_skipnom = new TH1D("h_weight_product_skipnom","",200,-1,3);
  for(int iev = 0; iev < nev; ++iev){
    analysis_tree->GetEntry(iev);

    // Exposure for assumed data - Weight MC templates to that exposure
    // Calculate POT weight (assumed data exposure) / (generated MC exposure)
    double pot_wght = pot_weight*data_pot;

    double this_entry = weight;
    this_entry *= pot_wght;

    // Get the current entry's bin
    std::vector<double> binsvec;
    if(util::IsReco(kv)){
      binsvec.push_back(enu_reco);
      LOG("valor_sbn",pDEBUG) << "Pushing a reconstructed energy of " << enu_reco << " onto pdf";
    }
    else{
      binsvec.push_back(enu_true);
      LOG("valor_sbn",pDEBUG) << "Pushing a true energy of " << enu_true << " onto pdf";
    }
    const int ibin = bins_sb->FindBin(binsvec);

    // Ensure the bin is refilled until everything has been accounted for
    pdfnom->AddToEntry(this_entry, ibin);

    for(int itoy = 0; itoy < ntoys; ++itoy){
      // Setup universe arrays
      double uni_wght = this_entry;
      double product_wght = 1;
      int iuni = itoy;

      unsigned int n_varied = 0;
      // Get the product of the 1 sigma weights for all sysematics
      for(const string &syst : systematics){
        // Get the total universe weight
        const std::vector<float> & universes = branch_universes.at(syst);
        const int &n = branch_n_universes.at(syst);
        // If we are NOT looking at a correlated flux parameter, extract the random universe variation
        if(std::find(flux_corr.begin(), flux_corr.end(), syst) == flux_corr.end())
          iuni = syst_toy_universe.at(syst).at(itoy);

        if(iuni > n){
          LOG("universe_variation_tweaks",pFATAL) << " The universe that is being accessed: " << iuni << " exceeds the number available for this parameter: " << n << ", " << syst;
          util::Die(-1);
        }

        // Get the weight for this parameter at that universe
        // Check the universe weight is a valid number and not 0
        if(isnan(universes.at(iuni)) || universes.at(iuni) < std::numeric_limits<float>::epsilon()){
          LOG("universe_variation_tweaks", pDEBUG) << " The weight for systematic " << syst << ", universe " << iuni << " and throw " << itoy << " is " << universes.at(iuni) << ", skipping";
          continue;
        }
        uni_wght *= universes.at(iuni);
        product_wght *= universes.at(iuni);
        LOG("universe_variation_tweaks", pDEBUG) << " The weight for systematic " << syst << ", universe " << iuni << " and throw " << itoy << " is " << universes.at(iuni);

        if(std::abs(universes.at(iuni)-1) > std::numeric_limits<double>::epsilon()){
          n_varied++;
        }
      }// Systematic
      
      h_weight_product->Fill(product_wght);
      // If we are only looking at correlated flux parameters, make sure at least one has been varied in this event and universe
      if(std::abs(uni_wght-this_entry) < std::numeric_limits<double>::epsilon() && (params_to_tweak == "flux_corr" || params_to_tweak == "pion")){
        // Check if the current universe ID is in the skipped vector, if not, add it
        if(std::find(skipped_universe_ids.begin(), skipped_universe_ids.end(), itoy) == skipped_universe_ids.end())
          skipped_universe_ids.push_back(itoy);
        skipped_events++;
        LOG("universe_variation_tweaks", pDEBUG) << " In event " << iev
                                                 << "\n The product of the universe weights is: " << uni_wght
                                                 << "\n The nominal event rate is: " << this_entry
                                                 << "\n These are too similar for any weight to have been applied, skipping";
      }
      LOG("universe_variation_tweaks", pDEBUG) << " Total universe weight for entry " << iev << " : " << uni_wght;
      h_weight_product_skipnom->Fill(product_wght);
      pdfuni[itoy]->AddToEntry(uni_wght,ibin);
      accepted_events++;
      if(n_varied == 1){
        one_varied++;
      }
    }// Throw 
  }// Events
  LOG("universe_variation_tweaks", pINFO) << " Finished analysis loop";

  // Now set the PDFs as histograms
  const char * nom_name = "h_nom_pdf";
  TH1D * nom_pdf = util::binning::Make1DHistogram(*pdfnom,binning,ks,true);
  nom_pdf->SetName(nom_name);

  LOG("universe_variation_tweaks", pINFO) << " Filled nominal histogram with " << nom_pdf->GetNbinsX() << " bins";

  std::vector<TH1D*> huni;
  for(int itoy = 0; itoy < ntoys; ++itoy){
    string name = "h_uni_var_pdf_"+itoy;
    huni.push_back(util::binning::Make1DHistogram(*pdfuni[itoy],binning,ks,true));
    huni.at(itoy)->SetName(name.c_str());
  }
  LOG("universe_variation_tweaks", pINFO) << " Translated PDFs to histograms";

  // Histograms to store the spectra for the default, +1 sigma and -1 sigma variations 
  // Setup the per-detector histograms for variations to determine sigma from nominal
  double max = -999.;
  TH1D * nom_pdf_werr = new TH1D(*nom_pdf);
  TH1D * err_m1 = new TH1D(*nom_pdf);
  TH1D * err_p1 = new TH1D(*nom_pdf);
  nom_pdf_werr->Reset();
  err_m1->Reset();
  err_p1->Reset();
  for(int ebin = 1; ebin <= nom_pdf->GetNbinsX(); ebin++) {
    double nom_val   = nom_pdf->GetBinContent(ebin);
    double one_sigma = 1.;
    double sumsq     = 0;
    double mean      = 0.;

    for(int itoy = 0; itoy < ntoys; itoy++) {
      mean += huni[itoy]->GetBinContent(ebin);
    }
    mean = mean/static_cast<double>(ntoys);

    for(int itoy = 0; itoy < ntoys; itoy++) {
      double curr_val = huni[itoy]->GetBinContent(ebin);
      sumsq += TMath::Power(curr_val-nom_val,2.);
  //    sumsq += TMath::Power(curr_val-mean,2.);
    }//itoy

    one_sigma = TMath::Sqrt(sumsq/static_cast<double>(ntoys));
    err_m1->SetBinContent(ebin, nom_val - one_sigma);
    err_p1->SetBinContent(ebin, nom_val + one_sigma);
    nom_pdf_werr->SetBinError(ebin,one_sigma);
    if( nom_pdf->GetBinCenter(ebin) < enu_reco_max ) {
      max = TMath::Max(max, nom_val + one_sigma);
    }
  }//ebin
  LOG("universe_variation_tweaks", pINFO) << " Filled histograms";
  
  // Draw and (if necessary), write to file
  TFile * outfile = NULL;
  if (output_rootfile!="-" && output_rootfile!="none") {
    output_rootfile = util::ReplaceAll(output_rootfile,"[PARAMS]",params_to_tweak);
    output_rootfile = util::ReplaceAll(output_rootfile,"[SAMPLE]",util::AsString(sample));
    outfile = new TFile(output_rootfile.c_str(),"RECREATE");
    LOG("universe_variation_tweaks",pNOTICE)<<"Writing to file: "<<output_rootfile;
  }
  
  BeamConf_t beam = gmodel->GetBeam();
  if (outfile && outfile->IsOpen()) {
    TString name = TString::Format("env_sb_%s_%s", util::AsString(beam).c_str(), util::AsString(sample).c_str() );
    nom_pdf->Write((name+"_nom").Data());
    nom_pdf_werr->Write((name+"_nom_err").Data());
    err_p1 ->Write((name+"_p1" ).Data());
    err_m1 ->Write((name+"_m1" ).Data());
    h_weight_product->Write("all_weights");
    h_weight_product_skipnom->Write("all_weights_skip_nominal");

    delete nom_pdf;
    delete nom_pdf_werr;
    delete err_p1;
    delete err_m1;
    delete h_weight_product;
    delete h_weight_product_skipnom;
  }
  
  if (outfile) {
    outfile->Close();
    delete outfile;
  }
  
  std::cout << "Done!"<<std::endl;
  std::cout << "Total number of events across all universes                   : " << accepted_events+skipped_events << std::endl;
  std::cout << "Total number of events with exactly one variation             : " << one_varied << std::endl;
  std::cout << "Total number of events with exactly no variations             : " << skipped_events << std::endl;
  std::cout << "Total number of universes containing events with no variations: " << skipped_universe_ids.size() << std::endl;
  std::cout << "Total number of universes                                     : " << ntoys << std::endl;
  
  
  time_t endrawtime;
  struct tm * endtimeinfo;
  time (&endrawtime);
  endtimeinfo = localtime (&endrawtime);
  time_t total_time = endrawtime-scriptstartrawtime;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End local time and date  :  " << asctime(endtimeinfo);
  std::cout << " Total time taken         :  " << total_time/static_cast<double>(60) << " minutes" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

}
