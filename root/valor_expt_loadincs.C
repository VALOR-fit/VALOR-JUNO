int valor_expt_loadincs () {
  TString topdir = gSystem->Getenv("VALOR_EXPT");

  TString mp = gROOT->GetMacroPath();
  TString ip;

  TString inc = topdir + "/src/";

  const char* p = inc.Data();
  if (p) {
    mp += ":";
    mp += p;
    ip += " -I";
    ip += p;
  }

  mp += ":/usr/local/include/";
  ip += "  -I/usr/local/include/";

  gROOT->SetMacroPath(mp.Data());
  gSystem->SetIncludePath(ip);

  // additions to .include must be done individually or CINT will
  // try to quote all the spaces as a single path

  TString dip = ".include ";
  dip += inc.Data();
  gROOT->ProcessLine(dip.Data());

  dip = ".include /usr/local/include/";
  gROOT->ProcessLine(dip.Data());

  dip = ".include /usr/local/include/log4cpp";
  gROOT->ProcessLine(dip.Data());

  dip = ".include /opt/local/include/";
  gROOT->ProcessLine(dip.Data());

  dip = ".include /opt/local/include/log4cpp";
  gROOT->ProcessLine(dip.Data());
  
  return 0;
}

