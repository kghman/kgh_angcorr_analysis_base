CPPC=g++
ROOTFLAGS=`root-config --cflags --glibs`
CFLAGS = -g -L/usr/lib -I/usr/local/include/gsl -I$(HOME)/vol/include -I$(HOME)/kgh_angcorr_analysis_base/common/include `gsl-config --cflags --libs` -DIKPLIBS_LOCAL -DHAVE_GNU_READLINE
GSLFLAGS = -lgsl -lgslcblas -lm

SRCDIR = $(HOME)/kgh_angcorr_analysis_base/common/src
INCLDIR = $(HOME)/kgh_angcorr_analysis_base/common/include
OBJDIR = $(HOME)/kgh_angcorr_analysis_base/common/objs

REQOBJS = $(OBJDIR)/Reaction.o $(OBJDIR)/Decay.o $(OBJDIR)/WikoParticleClass.o $(OBJDIR)/WikoGammaClass.o $(OBJDIR)/BeesHolder.o $(OBJDIR)/GrandBeesHolder.o $(OBJDIR)/FrameConverter.o $(OBJDIR)/SabreDetector.o $(OBJDIR)/Vec3.o $(OBJDIR)/Rotation.o $(OBJDIR)/wiko_wiko5.o
EXECS = bees geebees gambees

all: $(EXECS)

geebees: $(REQOBJS) geebees.cpp
	@echo Generating geebees...
	@$(CPPC) $(CFLAGS) -o $@ $^ $(ROOTFLAGS) $(GSLFLAGS)

bees: $(REQOBJS) bees.cpp
	@echo Generating bees...
	@$(CPPC) $(CFLAGS) -o $@ $^ $(ROOTFLAGS) $(GSLFLAGS)

gambees: $(REQOBJS) gambees.cpp
	@echo Generating gambees...
	@$(CPPC) $(CFLAGS) -o $@ $^ $(ROOTFLAGS) $(GSLFLAGS)

$(OBJDIR)/wiko_wiko5.o: $(SRCDIR)/wiko_wiko5.cpp $(INCLDIR)/wiko_wiko5.hpp $(INCLDIR)/wiko_types5.h
	@echo Compiling $@...
	@$(CPPC) $(CFLAGS) -o $@ -c $< $(GSLFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INCLDIR)/%.h
	@echo Compiling $@...
	@$(CPPC) $(CFLAGS) -o $@ -c $< $(ROOTFLAGS)

%.o: %.cpp
	@echo Compiling $@...
	@$(CPPC) $(CFLAGS) -c $^ $(ROOTFLAGS)

clean_bees:
	@rm -f B* sigma densmat msubs

clean:
	@rm -f *~ B* sigma densmat msubs *.o $(REQOBJS) $(EXECS)
