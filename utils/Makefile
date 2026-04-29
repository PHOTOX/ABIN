# To compile Fortran utilities, type 'make utils' in the root folder
# Note that some of these utilities use subroutines on ABIN.

all: abin-randomint analyze_movie
	cp abin-randomint analyze_movie Exitabin.sh fetchabin.sh prum.awk checkenergy.awk ../bin

abin-randomint: abin-randomint.o
	${FC} abin-randomint.o $(LDFLAGS) ${LDLIBS} -o $@

analyze_movie: analyze_movie.o
	${FC} analyze_movie.o $(LDFLAGS) ${LDLIBS} -o $@

clean:
	/bin/rm -f *.o *.mod abin-randomint analyze_movie

.PHONY: clean

%.o: %.f90
	$(FC) $(FFLAGS) $(DFLAGS) -I../src/ -c $<
