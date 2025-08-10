CC = gcc -g
CFLAGS = -I$(HOME)/vol/include -I/home2/stripsi/include -g -DIKPLIBS_LOCAL -DHAVE_GNU_READLINE

MAFIA_OBJS = 2d_fit_io.o
MAFIA_LIBS = -lvs -lhwc 

OBJS = corleone.o wiko_wiko.o menue.o
MATHLIBS = -L/usr/local/cernlib/99/lib -lmathlin -lkernlib -lfor
LIBS = -L$(HOME)/vol/lib/`uname` -L/home2/stripsi/lib -L$(HOME)/vol/lib/ -lmfile -lfitlib -lreadline -ltermcap -lgsl -lgslcblas -lm

corleone: $(OBJS)
	$(CC) $(CFLAGS) -o corleone $(OBJS) $(LIBS) 
capone: capone.o wiko_wiko4.o menue.o
	$(CC) $(CFLAGS) -o capone capone.o wiko_wiko4.o menue.o $(LIBS)

menue.o: menue.c menue.h 
2d_fit_io.o: 2d_fit_io.c 2d_fit_io.h
corleone.o: corleone.c
capone.o: capone.c
install: corleone
	cp corleone $(HOME)/bin/corleone
install_cap: capone
	cp capone $(HOME)/bin/capone
clean: 
	-rm corleone *.o core *~
wiko_wiko.o: wiko_wiko.c wiko_wiko.h
	$(CC) $(CFLAGS) -c wiko_wiko.c
wiko_wiko2.o: wiko_wiko2.c wiko_wiko2.h
	$(CC) $(CFLAGS) -c wiko_wiko2.c
wiko_wiko3.o: wiko_wiko3.c wiko_wiko3.h wiko_types3.h
	$(CC) $(CFLAGS) -c wiko_wiko3.c
wiko_wiko4.o: wiko_wiko4.c wiko_wiko4.h wiko_types4.h
	$(CC) $(CFLAGS) -c wiko_wiko4.c
wiko_test.o: wiko_test.c wiko_wiko.h
	$(CC) $(CFLAGS) -c wiko_test.c
wiko_alphatest4.o: wiko_alphatest4.c wiko_wiko4.h
	$(CC) $(CFLAGS) -c wiko_alphatest4.c
wiko_alphatest5.o: wiko_alphatest5.c wiko_wiko5.h
	$(CC) $(CFLAGS) -c wiko_alphatest5.c
msutest.o: msutest.c wiko_wiko2.h
	$(CC) $(CFLAGS) -c msutest.c
wiko_test: wiko_test.o wiko_wiko.o wiko_wiko.h 
	$(CC) $(CFLAGS) -o wiko_test wiko_test.o wiko_wiko.o -lm 
wiko_test_old: wiko_test.o wiko_wiko.o wiko_wiko.h 
	$(CC) $(CFLAGS) -o wiko_test_old wiko_test.o wiko_wiko.o -lm -lmathlib -lkernlib
wiko_test_new: wiko_test.o wiko_wiko2.o wiko_wiko2.h 
	$(CC) $(CFLAGS) -o wiko_test_new wiko_test.o wiko_wiko2.o -lm -lmathlib -lkernlib -lfor
wiko_alphatest: wiko_alphatest.o wiko_wiko2.o wiko_wiko2.h 
	$(CC) $(CFLAGS) $(LIBS) -o wiko_alphatest wiko_alphatest.o wiko_wiko2.o -lm -lmathlib -lkernlib -lfor
wiko_alphatest5: wiko_alphatest5.o wiko_wiko5.o wiko_wiko5.h 
	$(CC) $(CFLAGS) $(LIBS) -o wiko_alphatest5 wiko_alphatest5.o wiko_wiko5.o -lm 
wiko_alphatest4: wiko_alphatest4.o wiko_wiko4.o wiko_wiko4.h 
	$(CC) $(CFLAGS) $(LIBS) -o wiko_alphatest4 wiko_alphatest4.o wiko_wiko4.o -lm 
wiko_alphatest3: wiko_alphatest3.o wiko_wiko3.o wiko_wiko3.h wiko_types3.h 
	$(CC) $(CFLAGS) $(LIBS) -o wiko_alphatest3 wiko_alphatest3.o wiko_wiko3.o -lm -lmathlib -lkernlib -lfor
msutest: msutest.o wiko_wiko2.o wiko_wiko2.h 
	$(CC) $(CFLAGS) $(LIBS) -o msutest msutest.o wiko_wiko2.o -lm -lmathlib -lkernlib -lfor
rot_mat.o: rot_mat.c wiko_wiko.h
	$(CC) $(CFLAGS) -c rot_mat.c
rot_mat: rot_mat.o wiko_wiko.o wiko_wiko.h 
	$(CC) $(CFLAGS) -o rot_mat rot_mat.o wiko_wiko.o -lm 

