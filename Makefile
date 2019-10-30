FF	= gfortran
LD	= gfortran
CC	= g++

LDFLAGS	= `root-config --libs`  -lgfortran
FFFLAGS = -g -c
CCFLAGS = `root-config --cflags` -g -c

all: CF1D

CF1D: CF1D.o FsiTools.o FsiLednicky.o
	$(CC) $^ -o $@ $(LDFLAGS)



%.o: %.F
	$(FF) $^ -o $@ $(FFFLAGS)

%.o: %.cxx
	$(CC) $^ -o $@ $(CCFLAGS)

clean:
	rm -f *.o corr_test *.d
