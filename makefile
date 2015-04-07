helix_s: helix.o readdata.o readseq.o readcmdseq.o readparam.o setcharge.o \
	findYFW.o printparam.o findRH.o findNcap.o findCcap.o \
	findNend.o findCend.o findlipo.o findHbond.o findcoul.o \
        calcprob.o clearmatrix.o mtplymatrix.o printprob.o
	gcc helix.o readdata.o readseq.o readcmdseq.o readparam.o setcharge.o \
	findYFW.o printparam.o findRH.o findNcap.o findCcap.o \
	findNend.o findCend.o findlipo.o findHbond.o findcoul.o \
	calcprob.o clearmatrix.o mtplymatrix.o printprob.o -lm -o helix

helix.o: helix.c
	gcc helix.c -c
readdata.o: readdata.c
	gcc readdata.c -c
readseq.o: readseq.c
	gcc readseq.c -c
readcmdseq.o: readcmdseq.c
	gcc readcmdseq.c -c	
readparam.o: readparam.c
	gcc readparam.c -c
setcharge.o: setcharge.c
	gcc setcharge.c -c
findYFW.o: findYFW.c
	gcc findYFW.c -c
printparam.o: printparam.c
	gcc printparam.c -c
findRH.o: findRH.c
	gcc findRH.c -c
findNcap.o: findNcap.c
	gcc findNcap.c -c
findCcap.o: findCcap.c
	gcc findCcap.c -c
findNend.o: findNend.c
	gcc findNend.c -c
findCend.o: findCend.c
	gcc findCend.c -c
findlipo.o: findlipo.c
	gcc findlipo.c -c
findHbond.o: findHbond.c
	gcc findHbond.c -c
findcoul.o: findcoul.c
	gcc findcoul.c -c
calcprob.o: calcprob.c
	gcc calcprob.c -c
clearmatrix.o: clearmatrix.c
	gcc clearmatrix.c -c
mtplymatrix.o: mtplymatrix.c
	gcc mtplymatrix.c -c
printprob.o: printprob.c
	gcc printprob.c -c
