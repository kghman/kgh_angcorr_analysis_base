CPPC=g++
FC=gfortran
ROOTFLAGS=`root-config --cflags --glibs`
CFLAGS = -g -L/usr/lib -I/usr/local/include/gsl -I$(HOME)/vol/include -I$(HOME)/kgh_angcorr_analysis_base/common/include `gsl-config --cflags --libs` -DIKPLIBS_LOCAL -DHAVE_GNU_READLINE
GSLFLAGS = -lgsl -lgslcblas -lm

SRCDIR = $(HOME)/kgh_angcorr_analysis_base/common/src
INCLDIR = $(HOME)/kgh_angcorr_analysis_base/common/include
OBJDIR = $(HOME)/kgh_angcorr_analysis_base/common/objs

REQOBJS = $(OBJDIR)/Reaction.o $(OBJDIR)/Decay.o $(OBJDIR)/WikoParticleClass.o $(OBJDIR)/wiko_wiko5.o $(OBJDIR)/BeesHolder.o $(OBJDIR)/GrandBeesHolder.o $(OBJDIR)/FrameConverter.o $(OBJDIR)/SabreDetector.o $(OBJDIR)/Vec3.o $(OBJDIR)/Rotation.o
EXECS = pda

all: $(EXECS)

pda: $(REQOBJS) phi_decomp_analysis.cpp
	@echo Generating pda...
	@$(CPPC) $(CFLAGS) -o $@ $^ $(GSLFLAGS) $(ROOTFLAGS)

$(OBJDIR)/wiko_wiko5.o: $(SRCDIR)/wiko_wiko5.cpp $(INCLDIR)/wiko_wiko5.hpp $(INCLDIR)/wiko_types5.h
	@echo Compiling $@...
	@$(CPPC) $(CFLAGS) -o $@ -c $< $(GSLFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INCLDIR)/%.h
	@echo Compiling $@...
	@$(CPPC) $(CFLAGS) -o $@ -c $< $(ROOTFLAGS)

%.o: %.cpp
	@echo Compiling $@...
	@$(CPPC) $(CFLAGS) -c $^ $(ROOTFLAGS)

sorter: mreadrot_IW.f
	@echo Compiling $@...
	@$(FC) -o $@ $^

clean_bees:
	@rm -f B* sigma densmat msubs

clean:
	@rm -f *~ B* sigma densmat msubs *.o $(REQOBJS) $(EXECS)
