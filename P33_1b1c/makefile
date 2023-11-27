COMPILER = ifort

CFLAGS   = -O2 -qmkl -pg
CACFLAGS   = -coarray -O2 -qmkl -coarray-num-images=4
# CFLAGS = -qmkl -O0 -traceback -warn interfaces -fpe0 -check all -real-size 64 -check noarg_temp_created -real-size 64 -coarray-num-images=1
SINGLE_CACFLAGS   = -coarray -O2 -qmkl -coarray-num-images=1

# Lapack library
LALIBS = -lmkl_lapack95_lp64
LIBS = $(LALIBS)

#  Objects to create
OBJS = kinds.o \
       heft.o \
       hamiltonian.o \
       SMatrix.o \
       minf.o \
       quadpack.o \
       NumFort.o \
       bobyqa_module.o \
       co_random_numbers.o

FIN_OBJS = kinds.o \
           heft.o \
           NumFort.o \
	   hamiltonian.o

INF_OBJS = kinds.o \
           heft.o \
	   SMatrix.o \
           minf.o \
           quadpack.o \
           NumFort.o \
	   bobyqa_module.o


#  Pattern rule(s)
%.o : ../src/%.f90
	$(COMPILER) $(CFLAGS) -c $<


#  targets
# NumFort.o : ../src/NumFort.f90
# 	$(COMPILER) $(CFLAGS) -c $<

# minf.o : ../src/minf.f90
# 	$(COMPILER) $(CFLAGS) -c $<

# quadpack.o : ../src/quadpack.f90
# 	$(COMPILER) $(CFLAGS) -c $<

# hamiltonian.o : ../src/hamiltonian.f90
# 	$(COMPILER) $(CFLAGS) -c $<

# SMatrix.o : ../src/SMatrix.f90
# 	$(COMPILER) $(CFLAGS) -c $<

# bobyqa_module.o : ../src/bobyqa_module.f90
# 	$(COMPILER) $(CFLAGS) -c $<

co_random_numbers.o : ../src/co_random_numbers.f90
	$(COMPILER) $(CACFLAGS) -c $<

heft.o : ../src/heft.f90
	$(COMPILER) $(CACFLAGS) -c $<

fin.x : ../src/finiteVol.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(FIN_OBJS) $(LIBS)

mpiFin.x : ../src/mpiFiniteVol.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(FIN_OBJS) $(LIBS)

lqcdFin.x : ../src/lqcdFiniteVol.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(FIN_OBJS) $(LIBS)

multifitFin.x : ../src/multifitFinite.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(FIN_OBJS) $(LIBS)

fitBare.x : ../src/fitBareMassSlope.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(FIN_OBJS) bobyqa_module.o $(LIBS)

inf.x : ../src/infiniteVol.f90 $(OBJS)
	$(COMPILER) $(SINGLE_CACFLAGS) -o $@ $< $(INF_OBJS) $(LIBS)

fit.x : ../src/fitScattering.f90 $(OBJS)
	$(COMPILER) $(SINGLE_CACFLAGS) -o $@ $< $(INF_OBJS) $(LIBS)

fitBQ.x : ../src/fitScattering_bobyqa.f90 $(OBJS)
	$(COMPILER) $(SINGLE_CACFLAGS) -o $@ $< $(INF_OBJS) $(LIBS)

poles.x : ../src/poleSearch.f90 $(OBJS)
	$(COMPILER) $(SINGLE_CACFLAGS) -o $@ $< $(INF_OBJS) $(LIBS)

multifitPoles.x : ../src/multifitPoles.f90 $(OBJS)
	$(COMPILER) $(SINGLE_CACFLAGS) -o $@ $< $(INF_OBJS) $(LIBS)

TGrid.x : ../src/TMatGrid.f90 $(OBJS)
	$(COMPILER) $(SINGLE_CACFLAGS) -o $@ $< $(INF_OBJS) $(LIBS)

trainingInf.x : ../src/trainingInfiniteVol.f90 $(OBJS)
	$(COMPILER) $(CACFLAGS) -o $@ $< $(INF_OBJS) co_random_numbers.o $(LIBS)



#  Dependencies
kinds.o:

quadpack.o: kinds.o

minf.o: kinds.o

NumFort.o: kinds.o quadpack.o minf.o

heft.o: NumFort.o kinds.o

hamiltonian.o: kinds.o heft.o

SMatrix.o: kinds.o heft.o

bobyqa_module.o: kinds.o

co_random_numbers.o: kinds.o



# Clean up rule
clean:
	rm *.o *.mod *.x
