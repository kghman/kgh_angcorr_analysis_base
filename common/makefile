ROOTFLAGS = `root-config --cflags --glibs`
GSLFLAGS = `gsl-config --cflags --libs` -lgsl -lgslcblas -lm
CFLAGS = -g -I/home/kghman/gsl-2.7.1/specfunc -I/home/kghman/root/include -I/home/kghman/postgrad_work_2023/general_sim/include -DIKPLIBS_LOCAL -DHAVE_GNU_READLINE
CC=g++

SRCDIR = ./src
INCLDIR = ./include
OBJDIR = ./objs

REQOBJS = $(OBJDIR)/Reaction.o $(OBJDIR)/Decay.o $(OBJDIR)/WikoClass.o $(OBJDIR)/WikoParticleClass.o $(OBJDIR)/WikoGammaClass.o $(OBJDIR)/BeesHolder.o $(OBJDIR)/wiko_wiko5.o $(OBJDIR)/SabreDetector.o $(OBJDIR)/SelfFitter.o $(OBJDIR)/mathfuncs.o $(OBJDIR)/Vec3.o $(OBJDIR)/Rotation.o $(OBJDIR)/FrameConverter.o $(OBJDIR)/AnnularDetector.o $(OBJDIR)/StripDetector.o $(OBJDIR)/TargetLadder.o $(OBJDIR)/GrandBeesHolder.o

EXECS = Ecm purekin mykinace tls tsa trd

all: $(EXECS)

objects: $(REQOBJS)

mykinace: $(REQOBJS) mykinace_printer.cpp
	@echo Generating $@...
	@$(CC) $(CFLAGS) -o $@ $^ $(GSLFLAGS) $(ROOTFLAGS)
trd: $(REQOBJS) testing_ring_det.cpp
	@echo Generating $@...
	@$(CC) $(CFLAGS) -o $@ $^ $(GSLFLAGS) $(ROOTFLAGS)
tsa: $(REQOBJS) testing_solid_angle.cpp
	@echo Generating $@...
	@$(CC) $(CFLAGS) -o $@ $^ $(GSLFLAGS) $(ROOTFLAGS)
tls: $(REQOBJS) target_ladder_shadow.cpp
	@echo Generating $@...
	@$(CC) $(CFLAGS) -o $@ $^ $(GSLFLAGS) $(ROOTFLAGS)
	@$(CC) $(CFLAGS) -o $@ $^ $(GSLFLAGS) $(ROOTFLAGS)
Ecm: $(REQOBJS) Ecm.cpp
	@echo Generating $@...
	@$(CC) $(CFLAGS) -o $@ $^ $(GSLFLAGS) $(ROOTFLAGS)
purekin: $(REQOBJS) purekin.cpp
	@echo Generating $@...
	@$(CC) $(CFLAGS) -o $@ $^ $(GSLFLAGS) $(ROOTFLAGS)

$(OBJDIR)/SelfFitter.o: $(SRCDIR)/SelfFitter.cpp $(INCLDIR)/SelfFitter.h $(INCLDIR)/mathfuncs.h
	@echo Compiling $@...
	@$(CC) $(CFLAGS) -o $@ -c $< $(GSLFLAGS) $(ROOTFLAGS)

$(OBJDIR)/wiko_wiko5.o: $(SRCDIR)/wiko_wiko5.cpp $(INCLDIR)/wiko_wiko5.hpp $(INCLDIR)/wiko_types5.h
	@echo Compiling $@...
	@$(CC) $(CFLAGS) -o $@ -c $< $(GSLFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INCLDIR)/%.h
	@echo Compiling $@...
	@$(CC) $(CFLAGS) -o $@ -c $< $(GSLFLAGS) $(ROOTFLAGS)

clean:
	@rm -f $(EXECS) $(OBJDIR)/*.o *~ $(INCLDIR)/*.gch
