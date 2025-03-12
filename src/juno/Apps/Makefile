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

TGT = make_dst # \
      #make_response \


all: $(TGT)

# Generate DSTs
make_dst: FORCE
ifeq ($(JUNOANALYSIS_ENABLED),0)
ifdef GENIE
	$(CXX) $(CXXFLAGS) -c make_dst.cxx $(GINCLUDES)
	$(LD) $(LDFLAGS) make_dst.o $(GLIBRARIES) $(VALOR_GLIBRARIES) -o $(VALOR_EXPT)/bin/valor-juno_make_dst
else
	@echo "Skipping make_valortree due to no GENIE"
endif
else
	@echo "Skipping make_valortree due to no JUNOanalysis"
endif

# Generate response DSTs
make_response: FORCE
ifeq ($(JUNOANALYSIS_ENABLED),0)
ifdef GENIE
	$(CXX) $(CXXFLAGS) -c make_response.cxx $(GINCLUDES)
	$(LD) $(LDFLAGS) make_response.o $(GLIBRARIES) $(VALOR_GLIBRARIES) -o $(VALOR_EXPT)/bin/valor-make_response
else
	@echo "Skipping make_response due to no GENIE"
endif
else
	@echo "Skipping make_response due to no JUNOanalysis"
endif


#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(VALOR_EXPT)/bin/valor_juno-make_dst
	$(RM) $(VALOR_EXPT)/bin/valor_juno-make_response

distclean: FORCE
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_juno-make_dst
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_juno-make_response
FORCE:

# DO NOT DELETE
