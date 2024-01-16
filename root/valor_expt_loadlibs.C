int valor_expt_loadlibs()
{
  TString fmwk_dir = gSystem->Getenv("VALOR_FMWK");
  TString fmwk_loadlibs = "#include \""+fmwk_dir+"/root/valor_sdk_loadlibs.C\"";
  gROOT->ProcessLine(fmwk_loadlibs.Data());
  int valor_fmwk_libs_status = 0;
  gROOT->ProcessLine(TString::Format("*(int*)%p = valor_sdk_loadlibs();",&valor_fmwk_libs_status));
  
  TString baselibs = gSystem->GetDynamicPath();
  TString expt_libs = "$VALOR_EXPT/lib:";
  expt_libs += baselibs;
  
  gSystem->SetDynamicPath(expt_libs.Data());
  
  //gSystem->Load("libValorDuneConv.so");
  //gSystem->Load("libValorDuneUtil.so");
  
  return valor_fmwk_libs_status;
}

