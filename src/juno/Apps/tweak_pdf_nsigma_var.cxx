//____________________________________________________________________________
/*!

\app       valor_sbn-tweak_pdf_nsigma_var

\brief     Vary the MC spectra by setting systematics to +/- N sigma

\options     --output              : Where to save tweaked histograms to. [def: sigma_var_[PARAMS]_[nsig].root]
             --nsigma              : Number of sigma to tweak the spectra by, default 1
             --physics-hypothesis  : Physics hypothesis [def: "3f_osc"]
             --fit-param-opt       : Options used when setting up the FitParamList
                                     [def: "NONE"]
             --config              : Analysis configuration file [required]
             --mesg-thresholds     : Message thresholds [def: "verbosity.xml"]
             --params              : Which parameters to tweak [def: allsyst]
             --kinevar             : The kine var to store 1D PDFs for

             --input-mixing12,etc  : Option to give input values to osc params (in .arg)

             --load-prior-from-file : Option for file to give constrained prior matrix [ use default prior ]
             --prior-fit-number     : Option to specify the entry in the tree for above file [ 0 ]
             
             [--nthreads ]         : Number of CPU threads to use [ std::thread::hardware_concurrency ]
             

\author    Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
           University of Liverpool & STFC Rutherford Appleton Laboratory
           Steve Dennis <s.r.dennis \at liverpool.ac.uk>
           University of Liverpool
           Francis Bench <f.bench \at liverpool.ac.uk>
           University of Liverpool

\created   June 2019

\copyright Copyright (c) 2010-2020, The VALOR group. All rights reserved.

*/
//____________________________________________________________________________

// root
#include "TFile.h"
#include "TObjArray.h"
#include "TString.h"
#include "TSystem.h"

// fmwk
#include "valor/Conventions/Opt.h"
#include "valor/Conventions/Names.h"
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
#include "valor/FitElements/GenericSingleSampleFitModel.h"
#include "valor/FitCore/UtlFitParamList.h"
#include "valor/Utilities/OscParamValues.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Utilities/Env.h"
#include "valor/Utilities/VThreadPool.h"
#include "valor/Utilities/NameRegistry.h"
#include "valor/Utilities/UtlParam.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Ntuple/StdFitNtp.h"

using namespace std;
using namespace valor;
using namespace valor::conv;


bool t2k = false;


void SpecialCaseTweak(PhysParamId_t p, FitParamList * fp, double& tweak_p, double& tweak_m, const int &nsigma, const double &err);
void GeneratePDFNominal(FitModelABC * model, CompoundFitModel * cm, FitParamList * fp, KineVar_t kv, TFile * outfile);
void GeneratePDFTweaked(FitModelABC * model, CompoundFitModel * cm, FitParamList * fp, KineVar_t kv,
                        PhysParamId_t p, double tweak_dial, TFile * outfile, std::string tag);

void PrintMCT(FitModelABC * model, CompoundFitModel * cm, FitParamList * fp, KineVar_t kv,
              PhysParamId_t p, double tweak_dial, int reco_bin, int true_bin);

int main(int argc, char ** argv)
{

  valor::util::PrintBanner(Env::Instance()->Banner(),500);
 
  CmdLnArgList arg;
  arg.Read(argc,argv,"$VALOR_EXPT/config/apps/tweak_pdf_nsigma_var.arg:$VALOR_EXPT/config/apps/seeds.arg:$VALOR_EXPT/config/apps/oscparams.arg");
  LOG("valor_sbn", pNOTICE) << arg;
  
  // Read the command line
  int n_threads           = arg.GetInt("nthreads");
  double nsigma           = arg.GetDouble("nsigma");
  string hypname          = arg.GetString("physics-hypothesis");
  string fplopt           = arg.GetString("fit-param-opt");
  string output_rootfile  = arg.GetString("output");
  string params_to_tweak  = arg.GetString("params");
  string kinevar_str      = arg.GetString("kinevar");
  string universe_file    = arg.GetString("universe-spectra");
  bool def_spline_debug   = arg.GetBool("def-spline-debug");
  int reco_bin            = arg.GetInt("reco-bin");
  int true_bin            = arg.GetInt("true-bin");

  TString universe_spectra(universe_file.data());
  bool sigma_from_universes = false;
  std::vector<double> nsigmas;
  const std::vector<int> &integer_sigmas = {-3,-2,-1,0,1,2,3};
  int nsigma_int = -999;

  // If a file containing the universe spectra has been specified, 
  // open it and fill the vector with the integer values of sigma
  // to construct the variations with
  if(!universe_spectra.Contains("NONE")){
    sigma_from_universes = true;
    TFile *f_univ = new TFile(universe_file.data(),"READ");
    if(!f_univ->IsOpen()){
      LOG("valor_sbn",pFATAL) << " Error: Couldn't open universe file at " << universe_file;
      util::Die(-1);
    }
    TVectorD knot_vals(*static_cast<TVectorD*>(f_univ->Get("actual_knot_values")));
    for(int i = 0; i < knot_vals.GetNrows(); ++i){
      nsigmas.push_back(knot_vals[i]);
    }
  }
  else{
    nsigma_int = static_cast<int>(round(nsigma));
    nsigmas.push_back(nsigma);
  }

  // Re-fill the nsigmas vector to map correctly to the order in integer_sigmas
  if(sigma_from_universes){
    std::vector<double> temp_nsigmas;
    for(const int &i : integer_sigmas){
      double sigma_diff = 999.;
      double closest_sigma = 999.;
      for(const double &u : nsigmas){
        if(abs(u-i) < sigma_diff){
          sigma_diff = abs(u-i);
          closest_sigma = u;
        }
      }
      if(sigma_diff > 0.5){
        std::cout << "Big difference (> 0.5) between the selecteded integer sigma and the closest universe sigma!" << std::endl;
        std::cout << "Integer sigma: " << i << "       Closest universe sigma: " << closest_sigma << std::endl;
        std::cout << "Manually setting closest universe sigma equal to the integer sigma (should probably ignore the plot from this sigma)" << std::endl;
        closest_sigma = i;
      }
      temp_nsigmas.push_back(closest_sigma);
    }
    nsigmas.clear();
    nsigmas = temp_nsigmas;
    assert(nsigmas.size() == integer_sigmas.size());
  }

  std::cout << "nsigmas vec: " << std::endl;
  for(auto i = nsigmas.begin(); i != nsigmas.end(); i++){
      std::cout << *i << std::endl;
  }
  KineVar_t kv = util::KineVarFromString(util::RemoveAllWhitespace(kinevar_str));
  if(kv == kvUndefined) {
    cerr << "ERROR: input kine var string [" << util::RemoveAllWhitespace(kinevar_str)
         << "] was not understood. Exiting." << std::endl;
    util::Die(-1);
  }

  string read_prior_from_file = arg.GetString("load-prior-from-file");
  
  std::string mesg_thresholds = util::CompleteConfigPath(arg.GetString("mesg-thresholds"));
  Messenger::Instance()->Configure(mesg_thresholds);
  
  // Set up the thread pool
  if(n_threads>0) {//without n-cores specified, use as many as the host has
    util::VThreadPool::Instance()->SetNThreads(n_threads);
  }
  
  // Set the random seeds
  Rnd::Instance()->SetSeeds(arg);
  
  // And OA parameters
  Opt::Instance()->gMassHierarchyAssumed = ( (arg.GetBool("inverted"))     ? conv::kNmhInverted      : conv::kNmhNormal);

  conv::NuOscDmsqScheme_t dmsq_scheme = util::DmsqSchemeFromString(arg.GetString("dmsq-scheme"));
  if(dmsq_scheme == kNuDmsqUndefined) {
    LOG("valor_sbn", pFATAL) << "Could not parse the input Dmsq scheme: [" << arg.GetString("dmsq-scheme") << "]";
    util::Die(-1);
  }
  Opt::Instance()->gOscDmsqScheme = dmsq_scheme;

  // TODO: Update this with the new input param search
  double mixing12 = arg.GetDouble("input-mixing12");
  double mixing13 = arg.GetDouble("input-mixing13");
  double mixing23 = arg.GetDouble("input-mixing23");
  double dm12sq   = arg.GetDouble("input-dm12sq");
  double dm23sq   = arg.GetDouble("input-dm23sq");
  double deltacp  = arg.GetDouble("input-deltaCP");
  
  // Config
  std::string config_path = util::CompleteConfigPath(arg.GetString("config"));

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config_path << endl;
  cout << "- Physics hypothesis : " << hypname << endl;
  
  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(config_path.c_str());
  FitModelABC * model = elements->FitModel();
  CompoundFitModel * cm = dynamic_cast<CompoundFitModel*>(model);
  assert(cm);

  // Get the parameters to vary
  std::vector<PhysParamId_t> param_set = util::GetParameterSet(params_to_tweak, true);    // Only get exact matches to any input param names (otherwise problems can be caused)
  assert(param_set.size());

  // Make the output file, if we are getting the number of sigma from the universes
  // replace [NSIG] with 'all-sigma'
  TFile * outfile = NULL;
  output_rootfile = util::ReplaceAll(output_rootfile,"[PARAMS]",params_to_tweak);
  if(sigma_from_universes)
    output_rootfile = util::ReplaceAll(output_rootfile,"[NSIG]", "all_sigma");
  else
    output_rootfile = util::ReplaceAll(output_rootfile,"[NSIG]", (std::to_string(nsigma_int)+"sigma").c_str());
  outfile = new TFile(output_rootfile.c_str(),"RECREATE");
  LOG("valor_sbn",pNOTICE)<<"Writing to file: "<<output_rootfile;

  const CorrelatedSystematics * syst = cm->GetSystematics();
  FitParamList * fp = util::GetFitParameterList(hypname.c_str(), fplopt.c_str());

  if(hypname == "eff_sterile_osc"){                                                                                                                                                   
    fp->Update(conv::kPhypDm41sq,   100);                                                                                                                                           
    fp->Update(conv::kPhypMixingMuMu, 0);                                                                                                                                           
    fp->Update(conv::kPhypMixingMuE,  1);                                                                                                                                           
    fp->Update(conv::kPhypMixingEE,   0);                                                                                                                                           
  }

  else{   
    if (fp->Exists(conv::kPhypMixing12)) { fp->Update(conv::kPhypMixing12,mixing12); }
    if (fp->Exists(conv::kPhypMixing13)) { fp->Update(conv::kPhypMixing13,mixing13); }
    if (fp->Exists(conv::kPhypMixing23)) { fp->Update(conv::kPhypMixing23,mixing23); }
    if (fp->Exists(conv::kPhypDm12sq  )) { fp->Update(conv::kPhypDm12sq  ,dm12sq  ); }
    if (fp->Exists(conv::kPhypDm23sq  )) { fp->Update(conv::kPhypDm23sq  ,dm23sq  ); }
    if (fp->Exists(conv::kPhypCPV     )) { fp->Update(conv::kPhypCPV     ,deltacp ); }
  }

  // Generate pdfs
  GeneratePDFNominal(model, cm, fp, kv, outfile);
  for (unsigned int ip = 0 ; ip < param_set.size() ; ip++) {
    PhysParamId_t p = param_set[ip];
    if(!util::IsSystematic(p)) { continue; }
    double bestfit = fp->InitValue(p);
    double sigma_err = syst->OneSigmaError(p);
    for(const double &nsig : nsigmas){
      double tweak_p = bestfit + (nsig * sigma_err);
      double tweak_m = bestfit - (nsig * sigma_err);

      // Apply any special cases here (e.g. to force the tweaks to a certain value)
      // Currently sets MEC, we don't want this right now
      //SpecialCaseTweak(p, fp, tweak_p, tweak_m, nsigma, 0.5);

      // Get the number of sigma as a rounded integer for naming 
      // (assuming we only ever vary by integer values of sigma but may want to 
      // select exactly the variation given by a universe)
      nsigma_int = static_cast<int>(round(nsig));

      LOG("valor_sbn", pNOTICE)
        << "Tweak systematic parameter: "  << util::NameRegistry::Instance()->ParamAsString(p)
        << "\n  Nominal value = "          << bestfit
        << "\n  sigma error = "            << sigma_err;
      if(sigma_from_universes){
        if(nsig < 0)
          LOG("valor_sbn", pNOTICE) << "\n  nominal " << nsig << " sigma tweak = " << tweak_p;
        else
          LOG("valor_sbn", pNOTICE) << "\n  nominal +" << nsig << " sigma tweak = " << tweak_p;
      }
      else{
        if(nsig < 0){
          LOG("valor_sbn", pNOTICE) 
            << "\n  nominal +" << nsig << " sigma tweak = " << tweak_m
            << "\n  nominal -" << nsig << " sigma tweak = " << tweak_p;
        }
        else{
          LOG("valor_sbn", pNOTICE) 
            << "\n  nominal +" << nsig << " sigma tweak = " << tweak_p
            << "\n  nominal -" << nsig << " sigma tweak = " << tweak_m;
        }
      }
      
      // If we are getting the values of sigma from the universes then only generate a single pdf
      // and determine whether it is positive or negative
      if(sigma_from_universes){
        if(nsig < 0 && nsigma_int != 0){
          GeneratePDFTweaked(model, cm, fp, kv, p, tweak_p, outfile, ("m"+std::to_string(abs(nsigma_int))).data());
          if(def_spline_debug){PrintMCT(model, cm, fp, kv, p, tweak_p, reco_bin, true_bin);}
        }
        else{
          GeneratePDFTweaked(model, cm, fp, kv, p, tweak_p, outfile, ("p"+std::to_string(abs(nsigma_int))).data());
          if(def_spline_debug){PrintMCT(model, cm, fp, kv, p, tweak_p, reco_bin, true_bin);}
        }
      }
      else{
        if(nsig < 0){
          GeneratePDFTweaked(model, cm, fp, kv, p, tweak_m, outfile, ("p"+std::to_string(abs(nsigma_int))).data());
          GeneratePDFTweaked(model, cm, fp, kv, p, tweak_p, outfile, ("m"+std::to_string(abs(nsigma_int))).data());
        }
        else{
          GeneratePDFTweaked(model, cm, fp, kv, p, tweak_p, outfile, ("p"+std::to_string(abs(nsigma_int))).data());
          GeneratePDFTweaked(model, cm, fp, kv, p, tweak_m, outfile, ("m"+std::to_string(abs(nsigma_int))).data());
        }
      }
    } // loop over nsigmas
    if(def_spline_debug){PrintMCT(model, cm, fp, kv, p, 0, reco_bin, true_bin);}
  } // loop over params

  outfile->Close();
  delete outfile;

  LOG("valor_sbn",pNOTICE)<<"Finished writing to file: "<<output_rootfile;
  std::cout << "Done!" << std::endl;
}
//________________________________________________________________________________________
void SpecialCaseTweak(conv::PhysParamId_t p, FitParamList * fp, double& tweak_p, double& tweak_m, const int &nsigma, const double &err = 0.5)
{
  if(t2k) {
    if(p == util::NameRegistry::Instance()->PhysParamFromString("2p2h_shape_O")) {
      double bestfit = fp->InitValue(p);
      tweak_p = bestfit + nsigma*err;
      tweak_m = bestfit - nsigma*err;
    }
  }
}
//________________________________________________________________________________________
void GeneratePDFNominal(FitModelABC * model, CompoundFitModel * cm, FitParamList * fp, KineVar_t kv, TFile * outfile)
{
  cout << "\nGenerating nominal PDF:" << std::endl;
  outfile->cd();
  model->ResetSystematics(*fp,"all");
  fp->FromPhysParamList(model->Tweaks(),"current-silent");

  // Generate predictions
  cout << "with fit param list:" << std::endl;
  fp->Print(cout);
  model->Generate(*fp,"+SB");

  // Array to store results
  TObjArray * array = new TObjArray();
  array->SetOwner(true);

  // Write stuff out
  Option_t * write_opt = "pdfts_rf+SB";
  if (util::IsTrue(kv)) { write_opt = "pdfcs_rf+SB"; }    // The true PDFs don't have totals in GSSFM, only the mode breakdown // TODO: util::binning::Make1DHistogram currently doesn't work with truth...
  model->Write(array,"nom",write_opt);

  // Pick results wanted
  for(unsigned int imodel = 0; imodel < cm->GetNSubModels(); imodel++) {
    GenericSingleSampleFitModel * this_model = dynamic_cast<GenericSingleSampleFitModel*>(cm->GetSubModel(imodel));
    Detector_t d = this_model->GetDetector();
    BeamConf_t b = this_model->GetBeam();
    Sample_t   s = this_model->GetSample();

    // Get the total PDF
    std::string model_name = util::nm::MdlName(d, b, s, util::KineSpaceFromKineVar(kv));
    if (util::IsReco(kv)) {
      TString name_tot     = TString::Format("nom_pdfts_sb_rf_%s_%s", model_name.c_str(), util::AsString(kv).c_str());
      TH1D * curr_pdf_tot = (TH1D*) array->FindObject(name_tot.Data());
      if (!curr_pdf_tot) {
        cerr<<"Could not find histogram "<<name_tot<<endl;
        exit(-1);
      }
      else {
        TString save_name = TString::Format("env_sb_%s_%s_nom", util::AsString(b).c_str(), util::AsString(s).c_str());
        curr_pdf_tot->Write(save_name.Data());
        LOG("valor_sbn",pDEBUG) << " Nominal integrated total: " << curr_pdf_tot->Integral();
      }
    }
    /*
    // Get the PDF broken down by mode contributions
    for(conv::ReacMode_t m : this_model->GetReactionModes()) {
      TString name_contrib = TString::Format("nom_pdfcs_sb_rf_%s_%s_%s", model_name.c_str(), util::AsString(kv).c_str(), util::AsString(m).c_str());
      TH1D * curr_pdf_contrib = (TH1D*) array->FindObject(name_contrib.Data());
      if (!curr_pdf_contrib) {
        cerr<<"Could not find histogram "<<name_contrib<<endl;
        exit(-1);
      }
      TString save_name = TString::Format("env_sb_%s_%s_%s", model_name.c_str(), util::AsString(m).c_str(), "nom");
      curr_pdf_contrib->Write(save_name.Data());
    } //imode
  */
  } //imodel
  array->Clear("C");
  delete array;
}
//________________________________________________________________________________________
void GeneratePDFTweaked(FitModelABC * model, CompoundFitModel * cm, FitParamList * fp, KineVar_t kv,
                        PhysParamId_t p, double tweak_dial, TFile * outfile, std::string tag)
{
  outfile->cd();
  std::string param_name = util::NameRegistry::Instance()->ParamAsString(p);
  cout << "\nGenerating " << tag << " PDF with tweak dial = " << tweak_dial << " for param " << param_name << std::endl;
  model->ResetSystematics(*fp,"all");

  model->TweakSystematic(p, tweak_dial, "");
  fp->FromPhysParamList(model->Tweaks(),"current-silent");

  // Generate predictions
  model->Generate(*fp,"+SB");

  // Array to store results
  TObjArray * array = new TObjArray();
  array->SetOwner(true);

  // Write stuff out
  Option_t * write_opt = "pdfts_rf+SB";
  if (util::IsTrue(kv)) { write_opt = "pdfcs_rf+SB"; }    // The true PDFs don't have totals in GSSFM, only the mode breakdown
  model->Write(array,"curr",write_opt);

  // Pick results wanted
  for(unsigned int imodel = 0; imodel < cm->GetNSubModels(); imodel++) {
    GenericSingleSampleFitModel * this_model = dynamic_cast<GenericSingleSampleFitModel*>(cm->GetSubModel(imodel));
    Detector_t d = this_model->GetDetector();
    BeamConf_t b = this_model->GetBeam();
    Sample_t   s = this_model->GetSample();

    // Get the total PDF
    std::string model_name = util::nm::MdlName(d, b, s, util::KineSpaceFromKineVar(kv));
    if (util::IsReco(kv)) {
      TString name = TString::Format("curr_pdfts_sb_rf_%s_%s", model_name.c_str(), util::AsString(kv).c_str());
      TH1D * curr_pdf = (TH1D*) array->FindObject(name.Data());
      if (!curr_pdf) {
        cerr<<"Could not find histogram "<<name<<endl;
        exit(-1);
      }
      TString save_name = TString::Format("env_sb_%s_%s_%s", util::AsString(b).c_str(), util::AsString(s).c_str(), tag.c_str());
      curr_pdf->Write(save_name.Data());
      LOG("valor_sbn",pDEBUG) << " Tweaked integrated total: " << curr_pdf->Integral();
    }

    /*
    for (conv::ReacMode_t m : this_model->GetReactionModes()) {
      TString name_contrib = TString::Format("curr_pdfcs_sb_rf_%s_%s_%s", model_name.c_str(), util::AsString(kv).c_str(), util::AsString(m).c_str());
      TH1D * curr_pdf_contrib = (TH1D*) array->FindObject(name_contrib.Data());
      if (!curr_pdf_contrib) {
        cerr<<"Could not find histogram "<<name_contrib<<endl;
        exit(-1);
      }
      TString save_name = TString::Format("env_sb_%s_%s_%s_%s", model_name.c_str(), util::AsString(m).c_str(),
                        param_name.c_str(), tag.c_str());
      curr_pdf_contrib->Write(save_name.Data());
    } //mode
    */
  } //imodel
  array->Clear("C");
  delete array;
}
//________________________________________________________________________________________
void PrintMCT(FitModelABC * model, CompoundFitModel * cm, FitParamList * fp, KineVar_t kv,
            PhysParamId_t p, double tweak_dial, int reco_bin, int true_bin)
{
  std::string param_name = util::NameRegistry::Instance()->ParamAsString(p);
  model->ResetSystematics(*fp,"all");
  model->TweakSystematic(p, tweak_dial, "");
  fp->FromPhysParamList(model->Tweaks(),"current-silent");

  // Generate predictions
  model->Generate(*fp,"+SB");

  // Array to store results
  TObjArray * array = new TObjArray();
  //array->SetOwner(true);

  // Write stuff out
  Option_t * write_opt = "pdfts_rf+SB";
  if (util::IsTrue(kv)) { write_opt = "pdfcs_rf+SB"; }    // The true PDFs don't have totals in GSSFM, only the mode breakdown // TODO: util::binning::Make1DHistogram currently doesn't work with truth...
  model->Write(array,"nom","mct_vf");

  // Pick results wanted
  for(unsigned int imodel = 0; imodel < cm->GetNSubModels(); imodel++) {
    GenericSingleSampleFitModel * this_model = dynamic_cast<GenericSingleSampleFitModel*>(cm->GetSubModel(imodel));
    Detector_t d = this_model->GetDetector();
    BeamConf_t b = this_model->GetBeam();
    Sample_t   s = this_model->GetSample();

    const unsigned int n_mcts = array->GetEntries();
    std::map<ReacMode_t, VALORArray> mct_modes;
    for(unsigned int imct = 0; imct < n_mcts; ++imct){
      const VALORArray mct = *static_cast<VALORArray*>((*array)[imct]);
      const std::string mct_name = mct.GetName();
      const TString t_mct_name = mct_name;
      const size_t underscore_index = mct_name.find_last_of("_");
      if(underscore_index == string::npos){
        LOG("valor_sbn",pERROR) << " No '_' found in mct_name: " << mct_name;
        util::Die(-1);
      }
      const std::string mct_mode_name = mct_name.substr(underscore_index+1);
      LOG("valor_sbn",pINFO) << " MCT name " << mct_name << ", mct_mode_name: " << mct_mode_name;
      if(!(t_mct_name.Contains("mctfin"))) continue;
      for(const ReacMode_t &mode : this_model->GetReactionModes()) {
        const TString mode_name = util::AsString(mode);
        if(mct_mode_name == mode_name){
          if(!mct_modes.emplace(mode,mct).second) {
            LOG("valor_sbn",pERROR) << " Couldn't place the template in the map entry due to repeated entry for mode " << util::AsString(mode);
            util::Die(-1);
          }
          break;
        }
      }
    }
    
    std::cout << "Looking at " << param_name << " with sigma variation = " << tweak_dial << std::endl;                                                                                              
    std::cout << "Printing MCT event rate breakdown for reco bin = " << reco_bin << " and true bin = " << true_bin << ": " << std::endl;   
    
    std::vector<double> tot_event_mode(mct_modes.size(), 0);   
    int n_true_bins_sb = 0;
    const SingleSampleKineBinning & binning = *(this_model->GetBinning());

    for(auto elem : mct_modes){
      std::cout << util::AsString(elem.first) << " " << elem.second.GetEntry(reco_bin,true_bin) <<  "\n";
      // Get the number of true bins
      const VALORBinning * true_sb = binning.TrueSmpBinning(elem.first);
      n_true_bins_sb = true_sb->NBins();

    }
    // sum up the event rates from all the true bins and store in vector 
    std::cout << "\nPrinting MCT event rate breakdown for the sum of reco bin " << reco_bin << ": " << std::endl;
    for(int i = 0; i < n_true_bins_sb; i++){
      int j = 0;

      for(auto elem : mct_modes){
        tot_event_mode[j] += elem.second.GetEntry(reco_bin,i);
        j++;
      }
    }

    int i = 0;
    double tot_bin_event_rate = 0;
    for(auto elem : mct_modes){
      std::cout << util::AsString(elem.first) << " " << tot_event_mode[i] <<  "\n";
      tot_bin_event_rate += tot_event_mode[i];
      i++;
    }
    std::cout << "\nTotal event rate for reco bin " << reco_bin << ": " << std::setprecision(10) << tot_bin_event_rate << "\n\n";

    /*
    // Get the PDF broken down by mode contributions
    for(conv::ReacMode_t m : this_model->GetReactionModes()) {
      TString name_contrib = TString::Format("nom_pdfcs_sb_rf_%s_%s_%s", model_name.c_str(), util::AsString(kv).c_str(), util::AsString(m).c_str());
      TH1D * curr_pdf_contrib = (TH1D*) array->FindObject(name_contrib.Data());
      if (!curr_pdf_contrib) {
        cerr<<"Could not find histogram "<<name_contrib<<endl;
        exit(-1);
      }
      TString save_name = TString::Format("env_sb_%s_%s_%s", model_name.c_str(), util::AsString(m).c_str(), "nom");
      curr_pdf_contrib->Write(save_name.Data());
    } //imode
  */
  } //imodel
  array->Clear("C");
  delete array;
}

