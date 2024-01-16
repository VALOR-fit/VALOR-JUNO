//____________________________________________________________________________                                                       
/*!                                                                                                                                  

  \app       valor_sbn-write_spectra_using_universes

  \brief     Write out spectra to make a covariance for SBN

  \options
  --config              : Analysis configuration file [required]
  --output-name         : Output file name [def:test]
  --mesg-thresholds     : Message thresholds [def: "verbosity.xml"] 
  --input               : Input values of the oscillation parameters

  \author    Steve
  \created   June 2019

  \copyright (c) 2010-2019, The VALOR Fitting Group, All Rights Reserved

*/
//____________________________________________________________________________ 

#include <iostream>
// root
#include "TFile.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TString.h"
#include "TVirtualFitter.h"
#include "TPRegexp.h"
// fmwk
#include "valor/Conventions/Opt.h"
#include "valor/CmdLine/CmdLnArgList.h"
#include "valor/FitCore/PhysParamList.h"
#include "valor/FitCore/FitParamList.h"
#include "valor/FitCore/FitModelABC.h"
#include "valor/FitCore/FitDataABC.h"
#include "valor/FitCore/Likelihood.h"
#include "valor/FitCore/FitterABC.h"
#include "valor/FitCore/UtlFit.h"
#include "valor/FitCore/CompoundFitElements.h"
#include "valor/FitCore/UtlFitParamList.h"
#include "valor/FitElements/GenericSingleSampleFitModel.h"
#include "valor/Utilities/OscParamValues.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Utilities/Env.h"
#include "valor/Utilities/NameRegistry.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Ntuple/StdFitNtp.h"

// Proposal list from here:
// https://cdcvs.fnal.gov/redmine/projects/sbncode/wiki/Fitter_Validation_NuMu_sample

std::vector<std::string> interaction_systematics =
{ "genie_ccresAxial_Genie",
  "genie_ncresAxial_Genie",
  "genie_qema_Genie",
  "genie_NC_Genie",
  "genie_NonResRvbarp1pi_Genie",
  "genie_NonResRvbarp2pi_Genie",
  "genie_NonResRvp1pi_Genie",
  "genie_NonResRvp2pi_Genie",
};

std::vector<std::string> flux_systematics =
{ "expskin_FluxUnisim",
  "horncurrent_FluxUnisim",
  "kminus_PrimaryHadronNormalization",
  "kplus_PrimaryHadronFeynmanScaling",
  "kzero_PrimaryHadronSanfordWang",
  "nucleoninexsec_FluxUnisim",
  "nucleonqexsec_FluxUnisim",
  "nucleontotxsec_FluxUnisim",
  "piminus_PrimaryHadronSWCentralSplineVariation",
  "pioninexsec_FluxUnisim",
  "pionqexsec_FluxUnisim",
  "piontotxsec_FluxUnisim",
  "piplus_PrimaryHadronSWCentralSplineVariation"
};

std::vector<std::string> hadronic_flux_systematics =
{ "kminus_PrimaryHadronNormalization",
  "kplus_PrimaryHadronFeynmanScaling",
  "kzero_PrimaryHadronSanfordWang",
  "piminus_PrimaryHadronSWCentralSplineVariation",
  "piplus_PrimaryHadronSWCentralSplineVariation"
};
std::vector<std::string> piplus_flux_systematics =
{
  "piplus_PrimaryHadronSWCentralSplineVariation"
};
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

  CmdLnArgList arg;
  arg.Read(argc,argv,"$VALOR_EXPT/config/apps/universes.arg:$VALOR_EXPT/config/apps/oscparams.arg:");
  LOG("valor_sbn", pNOTICE) << arg;

  string mesg_thresholds  = arg.GetString("mesg-thresholds");
  string tag              = arg.GetString("tag");
  string location         = arg.GetString("output-location");
  string config_file      = arg.GetString("config");

  size_t ntoys          = arg.GetInt("n-toys");
  size_t first_universe = arg.GetInt("first-universe");

  std::string syst_list = arg.GetString("syst-list");

  std::vector<std::string> param_list;
  if (syst_list == "interaction") {
    param_list = interaction_systematics;
  }
  else if (syst_list == "flux") {
    param_list = flux_systematics;
  }
  else if (syst_list == "hadron") {
    param_list = hadronic_flux_systematics;
  }
  else if (syst_list == "pion") {
    param_list = piplus_flux_systematics;
  }
  else {
    LOG("write_sbn_spectra_using_universes",pFATAL) << "Unknown systematics list: "<<syst_list;
  }

  Messenger::Instance()->Configure(util::CompleteConfigPath(mesg_thresholds));

  LOG("valor_sbn", pNOTICE) << *Rnd::Instance();

  int max_universes = 499;

  // Get model & configure
  // Make an elements without reading the MC
  std::string elem_readopt("nomc nosplines");
  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(util::CompleteConfigPath(config_file),elem_readopt.c_str());

  // Then tell it about the universes and read the MC
  FitModelABC * model = elements->FitModel();
  CompoundFitModel * cmodel = dynamic_cast<CompoundFitModel*>(model);
  assert(cmodel);
  for (size_t i = 0 ; i < cmodel->GetNSubModels() ; i++) {
    GenericSingleSampleFitModel * gmodel = dynamic_cast<GenericSingleSampleFitModel*>(cmodel->GetSubModel(i));
    assert(gmodel);
    gmodel->SetUniverseRangeToRead(first_universe,first_universe+max_universes);
    for (const std::string & param : param_list) {
      gmodel->SetUniverseParameterToRead(param);
    }
    // And finally read the MC properly
    std::string model_readopt("forcedst nosplines");
    gmodel->Read(gmodel->GetLastConfigFileRead().c_str(),model_readopt.c_str());
    
  }
  
  // Generate predictions
  // Setup some global information
  FitParamList *fp = 0;                                                                                                                                                                 
  fp = util::GetFitParameterList("eff_sterile_osc");  

  const util::NameRegistry * nr = util::NameRegistry::Instance();                                                                                                                       
  conv::PhysParamId_t param = nr->PhysParamFromString("dm41sq");                                                                                                                        
  fp->Update(param,100);                                                                                                                                                                
  fp->SetInitValue(param,100);                                                                                                                                                          
  conv::PhysParamId_t param2 = nr->PhysParamFromString("sinsq2thMuE");                                                                                                                  
  fp->Update(param2,1);                                                                                                                                                                 
  fp->SetInitValue(param2,1);                                                                                                                                                           
  conv::PhysParamId_t param3 = nr->PhysParamFromString("sinsq2thMuMu");                                                                                                                 
  fp->Update(param3,0.0);                                                                                                                                                               
  fp->SetInitValue(param3,0.0);                                                                                                                                                         
  conv::PhysParamId_t param4 = nr->PhysParamFromString("sinsq2thEE");                                                                                                                   
  fp->Update(param4,0.0);                                                                                                                                                               
  fp->SetInitValue(param4,0.0);

  PhysHypothesis ph(conv::kHypVoid);
  Option_t * generate_opt = ""; // also generate non-default (for fitting) predictions in extra binning schemes
  std::string write_opt = "pdfts_rf pdfcs_rf";
  
  for (size_t itoy = 0 ; itoy < ntoys ; itoy++) {

    // Set our variations
    for (size_t i = 0 ; i < cmodel->GetNSubModels() ; i++) {
      GenericSingleSampleFitModel * gmodel = dynamic_cast<GenericSingleSampleFitModel*>(cmodel->GetSubModel(i));
      assert(gmodel);
      gmodel->SetSystematicVariationIDToApply("product",itoy);
    }

    // Regenerate predictions with systematic variations
    //model->Generate(ph,generate_opt);
    model->Generate(*fp,generate_opt);

    // Write stuff out                              
    TObjArray myarr;
    myarr.SetOwner(true);
    model->Write(&myarr,"valor",write_opt.c_str());

    TString output_fn = TString::Format("%s/%s_univ_%lu.root",location.c_str(),tag.c_str(),itoy);
    TFile output_file(output_fn.Data(),"RECREATE");

    TH1* testhist = dynamic_cast<TH1*>(myarr[0]);
    assert(testhist);

    for (int ihist = 0 ; ihist < myarr.GetEntries() ; ihist++) {
      TH1* hist = dynamic_cast<TH1*>(myarr[ihist]);
      assert(hist);
      TString old_name = hist->GetName();
      output_file.cd();
      hist->Write();

      if (old_name.Contains("pdfts")) {
        LOG("valor_sbn", pNOTICE) << "Wrote histogram "<<hist->GetName()<<" with "<<hist->Integral()<<" events.";
      }
    } // histograms
  } // itoy

  std::cout << "Done. Have a good day!"<<std::endl << "GO VALOR!!!" << std::endl;

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
