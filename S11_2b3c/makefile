COMPILER = ifort
STDLIB_CFLAGS := $(shell pkg-config --cflags fortran_stdlib)

CFLAGS   = -O2 -qmkl $(STDLIB_CFLAGS)
CACFLAGS   = -coarray -O2 -qmkl -coarray-num-images=1 $(STDLIB_CFLAGS)
# CFLAGS = -qmkl -O0 -traceback -warn interfaces -fpe0 -check all -real-size 64 -check noarg_temp_created -real-size 64 -coarray-num-images=1 $(STDLIB_CFLAGS)
SINGLE_CACFLAGS   = -coarray -O2 -qmkl -real-size 64 -coarray-num-images=1 $(STDLIB_CFLAGS)

# Lapack and numFort libs
LALIBS = -lmkl_lapack95_lp64
NUMLIBS = -I$(NumFortPath) $(NumFortPath)kinds.o $(NumFortPath)lapack.o $(NumFortPath)NumFort.o $(NumFortPath)minf.o $(NumFortPath)quadpack.o
STDLIB_LIBS := $(shell pkg-config --libs fortran_stdlib)

LIBS = $(NUMLIBS) $(LALIBS) $(STDLIB_LIBS)

#  Objects to create
OBJS = heft.o \
       hamiltonian.o \
       SMatrix.o \
       bobyqa_module.o \
       co_random_numbers.o

FIN_OBJS = heft.o \
	   hamiltonian.o

INF_OBJS = heft.o \
	   SMatrix.o \
	   bobyqa_module.o

#  Pattern rule(s)
%.o : %.f90
	$(COMPILER) $(CFLAGS) -c $< -I$(NumFortPath)

#  targets
hamiltonian.o : ../heftCode/hamiltonian.f90
	$(COMPILER) $(CFLAGS) -c $< -I$(NumFortPath)

SMatrix.o : ../heftCode/SMatrix.f90
	$(COMPILER) $(CFLAGS) -c $< -I$(NumFortPath)

bobyqa_module.o : ../heftCode/bobyqa_module.f90
	$(COMPILER) $(CFLAGS) -c $< -I$(NumFortPath)

co_random_numbers.o : ../heftCode/co_random_numbers.f90
	$(COMPILER) $(CACFLAGS) -c $< -I$(NumFortPath)

heft.o : ../heftCode/heft.f90
	$(COMPILER) $(CACFLAGS) -c $< -I$(NumFortPath)

fin.x : ../heftCode/finiteVol.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(FIN_OBJS) $(LIBS)

mpiFin.x : ../heftCode/mpiFiniteVol.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(FIN_OBJS) $(LIBS)

lqcdFin.x : ../heftCode/lqcdFiniteVol.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(FIN_OBJS) $(LIBS)

multifitFin.x : ../heftCode/multifitFinite.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(FIN_OBJS) $(LIBS)

fitBare.x : ../heftCode/fitBareMassSlope.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(FIN_OBJS) bobyqa_module.o $(LIBS)

corr.x : ../heftCode/correlationFn.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(FIN_OBJS) $(LIBS)

inf.x : ../heftCode/infiniteVol.f90 $(OBJS)
	$(COMPILER) $(SINGLE_CACFLAGS) -o $@ $< $(INF_OBJS) $(LIBS)

fit.x : ../heftCode/fitScattering.f90 $(OBJS)
	$(COMPILER) $(SINGLE_CACFLAGS) -o $@ $< $(INF_OBJS) $(LIBS)

fitBQ.x : ../heftCode/fitScattering_bobyqa.f90 $(OBJS)
	$(COMPILER) $(SINGLE_CACFLAGS) -o $@ $< $(INF_OBJS) $(LIBS)

poles.x : ../heftCode/poleSearch.f90 $(OBJS)
	$(COMPILER) $(SINGLE_CACFLAGS) -o $@ $< $(INF_OBJS) $(LIBS)

multifitPoles.x : ../heftCode/multifitPoles.f90 $(OBJS)
	$(COMPILER) $(SINGLE_CACFLAGS) -o $@ $< $(INF_OBJS) $(LIBS)

TGrid.x : ../heftCode/TMatGrid.f90 $(OBJS)
	$(COMPILER) $(SINGLE_CACFLAGS) -o $@ $< $(INF_OBJS) $(LIBS)

trainingInf.x : ../heftCode/trainingInfiniteVol.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(INF_OBJS) co_random_numbers.o $(LIBS)

MCInf.x : ../heftCode/MCInfiniteVol.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(INF_OBJS) co_random_numbers.o $(LIBS)

test.x : test.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(OBJS) $(LIBS)


#  Dependencies

heft.o:

hamiltonian.o: heft.o

SMatrix.o: heft.o

bobyqa_module.o:

co_random_numbers.o:


# Clean up rule
clean:
	rm *.o *.mod *.x
