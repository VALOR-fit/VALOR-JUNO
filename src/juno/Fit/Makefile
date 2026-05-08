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
	JunoDetectorResponse.cxx

OBJS = $(SRCS:.cxx=.o)

all: $(LIBRARY)

$(LIBRARY): $(OBJS)
	@echo " "
	@echo "** Building JUNO Fit library: $(LIBRARY)"
	$(LD) $(LDFLAGS) $(OBJS) $(GLIBRARIES) $(VALOR_GLIBRARIES) -shared -o $(LIBRARY)

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -fPIC -c $< $(GINCLUDES)

purge: FORCE
	$(RM) *.o *~ core

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(LIBRARY)

distclean: clean

FORCE: