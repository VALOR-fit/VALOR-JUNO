#
# Makefile for VALOR/JUNO apps
#
# Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(VALOR_EXPT)/make/Make.include

JUNO_LIBRARIES =

TGT = make_recoselection_minidst \
      make_caf_minidst \
      make_xsec_splines \
      make_flux_splines \
      universe_variation_tweaks \
      write_spectra_using_universes \
      tweak_pdf_nsigma_var 

all: $(TGT)

# Generate DSTs for JUNOCode fit
#
make_junocode_minidst: FORCE
ifeq ($(JUNOANALYSIS_ENABLED),1)
ifdef GENIE
	$(CXX) $(CXXFLAGS) -c make_recoselection_minidst.cxx $(GINCLUDES) $(JUNOANALYSIS_INC) $(LARSOFT_INCLUDES)
	$(LD) $(LDFLAGS) make_recoselection_minidst.o $(GLIBRARIES) $(VALOR_GLIBRARIES) $(JUNOANALYSIS_LIB) -o $(VALOR_EXPT)/bin/valor_juno-make_junocode_minidst
else
	@echo "Skipping make_recoselection_minidst due to no GENIE"
endif
else
	@echo "Skipping make_recoselection_minidst due to no JUNOanalysis"
endif

# Generate DSTs from CAF files
#
make_caf_minidst: FORCE
#ifeq ($(JUNOANALYSIS_ENABLED),1)
ifdef GENIE
	$(CXX) $(CXXFLAGS) -c make_caf_minidst.cxx $(GINCLUDES) $(JUNOANALYSIS_INC) $(LARSOFT_INCLUDES) 
	$(LD) $(LDFLAGS) make_caf_minidst.o $(GLIBRARIES) $(VALOR_GLIBRARIES) $(JUNOANALYSIS_LIB) -o $(VALOR_EXPT)/bin/valor_juno-make_caf_minidst
else
	@echo "Skipping make_caf_minidst due to no GENIE"
endif
#else
#	@echo "Skipping make_caf_minidst due to no JUNOanalysis"
#endif

# Generate DSTs for RecoSelection fit
#
make_recoselection_minidst: FORCE
ifeq ($(JUNOANALYSIS_ENABLED),1)
ifdef GENIE
	$(CXX) $(CXXFLAGS) -c make_recoselection_minidst.cxx $(GINCLUDES) $(JUNOANALYSIS_INC) $(LARSOFT_INCLUDES)
	$(LD) $(LDFLAGS) make_recoselection_minidst.o $(GLIBRARIES)  $(VALOR_GLIBRARIES) $(JUNOANALYSIS_LIB)  -o $(VALOR_EXPT)/bin/valor_juno-make_recoselection_minidst
else
	@echo "Skipping make_recoselection_minidst due to no GENIE"
endif
else
	@echo "Skipping make_recoselection_minidst due to no JUNOanalysis"
endif

# Generate JUNO flux splines
# Probably don't require the GENIE dependence but method isn't finalised yet so keeping just in case
#
make_flux_splines: FORCE
ifdef GENIE
	$(CXX) $(CXXFLAGS) -c make_flux_splines.cxx $(GINCLUDES)
	$(LD) $(LDFLAGS) make_flux_splines.o $(GLIBRARIES) $(VALOR_GLIBRARIES) -o $(VALOR_EXPT)/bin/valor_juno-make_flux_splines
else
	@echo "Skipping make_flux_splines due to no GENIE"
endif

# Generate JUNO interaction splines
# Probably don't require the GENIE dependence but method isn't finalised yet so keeping just in case
#
make_xsec_splines: FORCE
ifdef GENIE
	$(CXX) $(CXXFLAGS) -c make_xsec_splines.cxx $(GINCLUDES)
	$(LD) $(LDFLAGS) make_xsec_splines.o $(GLIBRARIES) $(VALOR_GLIBRARIES) -o $(VALOR_EXPT)/bin/valor_juno-make_xsec_splines
else
	@echo "Skipping make_xsec_splines due to no GENIE"
endif

# Tweaks taken directly from universe variations to validate VALOR systematic parameterisation
#
universe_variation_tweaks: FORCE
ifdef GENIE
	$(CXX) $(CXXFLAGS) -c universe_variation_tweaks.cxx $(GINCLUDES)
	$(LD) $(LDFLAGS) universe_variation_tweaks.o $(GLIBRARIES) $(VALOR_GLIBRARIES) -o $(VALOR_EXPT)/bin/valor_juno-universe_variation_tweaks
else
	@echo "Skipping universe_variation_tweaks due to no GENIE"
endif

write_spectra_using_universes: FORCE
ifdef GENIE
	$(CXX) $(CXXFLAGS) -c write_spectra_using_universes.cxx $(GINCLUDES)
	$(LD) $(LDFLAGS) write_spectra_using_universes.o $(GLIBRARIES) $(VALOR_GLIBRARIES) -o $(VALOR_EXPT)/bin/valor_juno-write_spectra_using_universes
else
	@echo "Skipping write_spectra_using_universes due to no GENIE"
endif

# Randomly throw the pdf to make envelopes
# 
tweak_pdf_nsigma_var: FORCE
ifdef GENIE
	$(CXX) $(CXXFLAGS) -c tweak_pdf_nsigma_var.cxx $(GINCLUDES)
	$(LD) $(LDFLAGS) tweak_pdf_nsigma_var.o $(GLIBRARIES) $(VALOR_GLIBRARIES) -o $(VALOR_EXPT)/bin/valor_juno-tweak_pdf_nsigma_var
else
	@echo "Skipping tweak_pdf_nsigma_var due to no GENIE"
endif

#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(VALOR_EXPT)/bin/valor_juno-make_junocode_minidst
	$(RM) $(VALOR_EXPT)/bin/valor_juno-make_caf_minidst
	$(RM) $(VALOR_EXPT)/bin/valor_juno-make_recoselection_minidst
	$(RM) $(VALOR_EXPT)/bin/valor_juno-make_flux_splines
	$(RM) $(VALOR_EXPT)/bin/valor_juno-make_xsec_splines
	$(RM) $(VALOR_EXPT)/bin/valor_juno-universe_variation_tweaks
	$(RM) $(VALOR_EXPT)/bin/valor_juno-write_spectra_using_universes
	$(RM) $(VALOR_EXPT)/bin/valor_juno-tweak_pdf_nsigma_var

distclean: FORCE
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_juno-make_junocode_minidst
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_juno-make_caf_minidst
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_juno-make_recoselection_minidst
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_juno-make_xsec_splines
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_juno-make_flux_splines
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_juno-universe_variation_tweaks
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_juno-write_spectra_using_universes
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_juno-tweak_pdf_nsigma_var
FORCE:

# DO NOT DELETE
