void define_valor_style()
{
  // VALOR style definition
  // Adopted from BaBar collaboration
  // Add the following lines to the start of your rootlogon.C file
  TStyle *valorStyle= new TStyle("VALOR","VALOR plots style");

  // use plain black on white colors
  valorStyle->SetFrameBorderMode(0);
  valorStyle->SetCanvasBorderMode(0);
  valorStyle->SetPadBorderMode(0);
  valorStyle->SetPadColor(0);
  valorStyle->SetCanvasColor(0);
  valorStyle->SetStatColor(0);
//valorStyle->SetFillColor(0); // this doesn't seem to be over-ridable
  valorStyle->SetLegendBorderSize(1); 

  // set the paper & margin sizes
  valorStyle->SetPaperSize(20,26);
  valorStyle->SetPadTopMargin(0.05);
  valorStyle->SetPadRightMargin(0.1);
  valorStyle->SetPadBottomMargin(0.16);
  valorStyle->SetPadLeftMargin(0.13);

  // use large Times-Roman fonts
  valorStyle->SetTextFont(132);
  valorStyle->SetTextSize(0.08);
  valorStyle->SetLabelFont(132,"x");
  valorStyle->SetLabelFont(132,"y");
  valorStyle->SetLabelFont(132,"z");
  valorStyle->SetLabelSize(0.05,"x");
  valorStyle->SetTitleSize(0.06,"x");
  valorStyle->SetLabelSize(0.05,"y");
  valorStyle->SetTitleSize(0.06,"y");
  valorStyle->SetLabelSize(0.05,"z");
  valorStyle->SetTitleSize(0.06,"z");
  valorStyle->SetLabelFont(132,"t");
  valorStyle->SetTitleFont(132,"x");
  valorStyle->SetTitleFont(132,"y");
  valorStyle->SetTitleFont(132,"z");
  valorStyle->SetTitleFont(132,"t");
  valorStyle->SetTitleFillColor(0);
  valorStyle->SetTitleX(0.25);
  valorStyle->SetTitleFontSize(0.08);
  valorStyle->SetTitleFont(132,"pad");

  // use bold lines and markers
  valorStyle->SetMarkerStyle(20);
  valorStyle->SetHistLineWidth(2);
  valorStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  valorStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  valorStyle->SetOptTitle(0);
  valorStyle->SetOptStat(0);
  valorStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  valorStyle->SetPadTickX(1);
  valorStyle->SetPadTickY(1);

  //add the stat box
  /*
  valorStyle->SetOptFit(1); // show fit results
  valorStyle->SetOptStat(1100); // show mean & RMS only
//valorStyle->SetStatStyle(0); //make stat boxes see through
  valorStyle->SetStatBorderSize(0);
  valorStyle->SetTextFont(42); //less bold

  // Add a greyscale palette for 2D plots
  
  int ncol=50;
  double dcol = 1./float(ncol);
  double gray = 1;
  TColor **theCols = new TColor*[ncol];
  for (int i=0;i<ncol;i++) theCols[i] = new TColor(999-i,0.0,0.7,0.7);
  for (int j = 0; j < ncol; j++) {
    theCols[j]->SetRGB(gray,gray,gray);
    gray -= dcol;
  }
  int ColJul[100];
  for  (int i=0; i<100; i++) ColJul[i]=999-i;
  valorStyle->SetPalette(ncol,ColJul);
  */

  //Define a nicer color palette (red->blue)
  //Uncomment these lines for a color palette (default is B&W)
  valorStyle->SetPalette(1,0);  // use the nice red->blue palette
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  valorStyle->SetNumberContours(NCont); 

}//End of definition of VALOR style

void print_banner()
{
  TString data_dir = valor::Env::Instance()->TopLevelDataDir();
  TString banner   = data_dir + "/banner/banner.dat";
  
  valor::util::PrintBanner(banner.Data(),500);
}

void config_mesg()
{
  TString conf_dir  = valor::Env::Instance()->TopLevelConfDir();
  TString verbosity = conf_dir + "/verbosity.xml";

  valor::util::Verbosity(verbosity.Data());
}

void valor_expt_utils()
{
    // Configure the messenger.
    config_mesg();
    
    std::string experiment_name = "untitled";
    std::string prompt = std::string("VALOR/")+experiment_name+": [%d] ";
    valor::util::ChangeInterpreterPrompt(prompt);
    
    define_valor_style();
    gROOT->SetStyle("VALOR");

    
    print_banner();
    
}
