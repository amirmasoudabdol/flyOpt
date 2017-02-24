# DO NOT CHANGE ANYTHING HERE!!! ##################################
# (unless you know *exactly* what you're doing...) 

# Utilites objects
FOBJ = zygotic.o fly_io.o maternal.o integrate.o translate.o solvers.o score.o \
         ../utils/error.o ../utils/distributions.o ../utils/random.o ../utils/ioTools.o ../utils/dSFMT.o ../utils/dSFMT_str_state.o

# Fly object
FSOBJ =  fly.o 

# Scatter Search objects
SSOBJ = ../ss/allocate.o ../ss/evaluate.o ../ss/init.o ../ss/local_search.o ../ss/recombine.o ../ss/refine.o ../ss/report.o ../ss/sort.o ../ss/ss.o ../ss/ssTools.o ../ss/stats.o ../ss/update.o 

# Enhanced Scatter Search objects
ESSOBJ = ../ess/ess.o ../ess/essAllocate.o ../ess/essEvaluate.o ../ess/essGoBeyond.o ../ess/essIO.o ../ess/essInit.o ../ess/essLocalSearch.o ../ess/essProblem.o ../ess/essRand.o ../ess/essRecombine.o ../ess/essSort.o ../ess/essStats.o ../ess/essTools.o



ifeq ($(METHOD), -DSS)
	METHODOBJ = $(SSOBJ)
else ifeq ($(METHOD), -DESS)
	METHODOBJ = $(ESSOBJ)
endif



#printscore objects
POBJ = zygotic.o fly_io.o maternal.o integrate.o translate.o solvers.o score.o printscore.o \
       ../utils/error.o ../utils/distributions.o ../utils/random.o ../utils/ioTools.o ../utils/dSFMT.o ../utils/dSFMT_str_state.o

#unfold objects
UOBJ = zygotic.o fly_io.o maternal.o integrate.o translate.o solvers.o score.o unfold.o \
	  ../utils/error.o ../utils/distributions.o ../utils/random.o ../utils/ioTools.o ../utils/dSFMT.o ../utils/dSFMT_str_state.o

#scramble objects
SOBJ = zygotic.o fly_io.o maternal.o integrate.o translate.o solvers.o score.o scramble.o \
	  ../utils/error.o ../utils/distributions.o ../utils/random.o ../utils/ioTools.o ../utils/dSFMT.o ../utils/dSFMT_str_state.o

SOURCES = `ls *.c`

#Below here are the rules for building things

all: $(FLYEXECS)

# special cases: dependencies and flags for individual .c files

fly.o: fly.c
	$(CC) -c $(CFLAGS) $(VFLAGS) fly.c

printscore.o: printscore.c
	$(CC) -c $(CFLAGS) $(VFLAGS) printscore.c

scramble.o: scramble.c
	$(CC) -c $(CFLAGS) $(VFLAGS) scramble.c

unfold.o: unfold.c
	$(CC) -c $(CFLAGS) $(VFLAGS) unfold.c

zygotic.o: zygotic.c
	$(CC) -c $(CFLAGS) $(KFLAGS) zygotic.c

# executable targets: 

$(execFile): $(FOBJ) $(FSOBJ)
	$(CC) -o $(execFile) $(CFLAGS) $(LDFLAGS) $(FOBJ) $(METHODOBJ) $(FLIBS) $(FSOBJ)

printscore: $(POBJ)
	$(CC) -o printscore $(CFLAGS) $(LDFLAGS) $(POBJ) $(LIBS) 

unfold: $(UOBJ)
	$(CC) -o unfold $(CFLAGS) $(LDFLAGS) $(UOBJ) $(LIBS) 

scramble: $(SOBJ)
	$(CC) -o scramble $(CFLAGS) $(LDFLAGS) $(SOBJ) $(LIBS) 

# ... and here are the cleanup and make deps rules

clean:
	rm -f *.o core*

Makefile: ${FRC}
	rm -f $@
	cp basic.mk $@
	echo "#Automatically generated dependencies list#" >> $@
	${CC} $(INCLUDES) -M ${SOURCES} >> $@
	chmod -w $@

