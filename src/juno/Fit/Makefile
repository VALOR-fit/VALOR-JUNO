#
# Makefile for VALOR/JUNO Fit library
#
# Ziou He <ziou.he \at warwick.ac.uk>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

include $(VALOR_EXPT)/make/Make.include

LIBNAME = libValorJunoFit.so
LIBRARY = $(VALOR_EXPT)/lib/$(LIBNAME)

SRCS = \
	JunoDetectorResponse.cxx \
	JunoSingleSampleFitModel.cxx \
	JunoFitDict.cxx

OBJS = $(SRCS:.cxx=.o)

DICT = JunoFitDict.cxx
DICT_PCM = JunoFitDict_rdict.pcm

all: $(LIBRARY)

$(DICT): JunoSingleSampleFitModel.h JunoDetectorResponse.h LinkDef.h
	@echo "** Generating JUNO Fit ROOT dictionary"
	rootcling -f $(DICT) -c $(GINCLUDES) \
		JunoSingleSampleFitModel.h \
		JunoDetectorResponse.h \
		LinkDef.h

$(LIBRARY): $(OBJS)
	@echo " "
	@echo "** Building JUNO Fit library: $(LIBRARY)"
	$(LD) $(LDFLAGS) $(OBJS) $(GLIBRARIES) $(VALOR_GLIBRARIES) -shared -o $(LIBRARY)
	@if [ -f $(DICT_PCM) ]; then cp $(DICT_PCM) $(VALOR_EXPT)/lib/; fi

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -fPIC -c $< $(GINCLUDES)

purge: FORCE
	$(RM) *.o *~ core $(DICT) $(DICT_PCM)

clean: FORCE
	$(RM) *.o *~ core $(DICT) $(DICT_PCM)
	$(RM) $(LIBRARY)
	$(RM) $(VALOR_EXPT)/lib/$(DICT_PCM)

distclean: clean

FORCE: