void valor_expt()
{   
    TString curr_dir = gSystem->pwd();
    gSystem->cd(curr_dir.Data());

    TString top_dir = gSystem->Getenv("VALOR_EXPT");
    
    int valor_incs_status = -999;
    int valor_libs_status = -999;
    
    TString loadincs = "#include \""+top_dir+"/root/valor_expt_loadincs.C\"";
    TString execute_loadincs = TString::Format("*(int*)%p = valor_expt_loadincs();",&valor_incs_status);
    
    TString loadlibs = "#include \""+top_dir+"/root/valor_expt_loadlibs.C\"";
    TString execute_loadlibs = TString::Format("*(int*)%p = valor_expt_loadlibs();",&valor_libs_status);
    
    gROOT->ProcessLine(loadincs.Data());
    gROOT->ProcessLine(execute_loadincs.Data());
    gROOT->ProcessLine(loadlibs.Data());
    gROOT->ProcessLine(execute_loadlibs.Data());
    
    //~ std::cerr<<"Incs = "<<valor_incs_status<<" ; libs = "<<valor_libs_status<<std::endl;
    
    TString loadutils = ".x " + top_dir + "/root/valor_expt_utils.C";
    gROOT->ProcessLine(loadutils.Data());
}
