CPPC=g++
CC=gcc
ROOTFLAGS=`root-config --cflags --glibs`
CFLAGS = -g -I$(HOME)/vol/include -I/home2/stripsi/include `gsl-config --cflags --libs` -DIKPLIBS_LOCAL -DHAVE_GNU_READLINE

rootfile: Atester
	./$^

Atester: testing_Adeltas.o wiko_wiko5.o
	@echo Generating Atester...
	@$(CPPC) -o $@ $^ $(CFLAGS) $(ROOTFLAGS)

testing_Adeltas.o: testing_Adeltas.cpp
	@echo Compiling $@...
	@$(CPPC) -c $^ $(CFLAGS) $(ROOTFLAGS)

wiko_wiko5.o: wiko_wiko5.c wiko_wiko5.h wiko_types5.h
	@echo Compiling $@...
	@$(CC) -c $< $(CFLAGS)

clean:
	@rm -f *.gch *.o *~ B* sigma
