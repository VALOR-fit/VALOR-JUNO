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

TGT = make_dst \
      make_dst_j24_flat
      #make_response \


all: $(TGT)

# Generate DSTs from GENIE/GHEP
make_dst: FORCE
ifeq ($(JUNOANALYSIS_ENABLED),0)
ifdef GENIE
	$(CXX) $(CXXFLAGS) -c make_dst.cxx $(GINCLUDES)
	$(LD) $(LDFLAGS) make_dst.o $(GLIBRARIES) $(VALOR_GLIBRARIES) -o $(VALOR_EXPT)/bin/valor_juno-make_dst
else
	@echo "Skipping make_dst due to no GENIE"
endif
else
	@echo "Skipping make_dst due to JUNOanalysis enabled"
endif

# Generate DSTs from J24 flat truth/reco + modeFriend
make_dst_j24_flat: FORCE
ifeq ($(JUNOANALYSIS_ENABLED),0)
	$(CXX) $(CXXFLAGS) -c make_dst_j24_flat.cxx $(GINCLUDES)
	$(LD) $(LDFLAGS) make_dst_j24_flat.o $(GLIBRARIES) $(VALOR_GLIBRARIES) -o $(VALOR_EXPT)/bin/valor_juno-make_dst_j24_flat
else
	@echo "Skipping make_dst_j24_flat due to JUNOanalysis enabled"
endif

# Generate response DSTs
make_response: FORCE
ifeq ($(JUNOANALYSIS_ENABLED),0)
ifdef GENIE
	$(CXX) $(CXXFLAGS) -c make_response.cxx $(GINCLUDES)
	$(LD) $(LDFLAGS) make_response.o $(GLIBRARIES) $(VALOR_GLIBRARIES) -o $(VALOR_EXPT)/bin/valor_juno-make_response
else
	@echo "Skipping make_response due to no GENIE"
endif
else
	@echo "Skipping make_response due to JUNOanalysis enabled"
endif


#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(VALOR_EXPT)/bin/valor_juno-make_dst
	$(RM) $(VALOR_EXPT)/bin/valor_juno-make_dst_j24_flat
	$(RM) $(VALOR_EXPT)/bin/valor_juno-make_response

distclean: FORCE
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_juno-make_dst
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_juno-make_dst_j24_flat
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_juno-make_response

FORCE:

# DO NOT DELETE