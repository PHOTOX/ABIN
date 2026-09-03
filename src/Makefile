F_OBJS := constants.o fortran_interfaces.o error.o modules.o mpi_wrapper.o files.o utils.o io.o random.o arrays.o qmmm.o fftw_interface.o \
          shake.o nosehoover.o gle.o transform.o potentials.o force_spline.o estimators.o ekin.o vinit.o plumed.o \
          remd.o force_bound.o water.o h2o_schwenke.o h2o_cvrqd.o force_h2o.o force_cp2k.o sh_integ.o surfacehop.o landau_zener.o potentials_sh.o\
          force_mm.o tera_mpi_api.o mace_mpi_api.o force_abin.o force_tcpb.o force_tera.o force_terash.o force_mace.o en_restraint.o analyze_ext_template.o geom_analysis.o analysis.o \
          minimizer.o mdstep.o forces.o cmdline.o init.o

C_OBJS := water_interface.o

# This is the default target
${BIN} : libabin.a
	${FC} ${LDFLAGS} abin.o ${LDLIBS} -o $@
	cp $(BIN) ../bin/

# Used for building Unit Tests
libabin.a: compile_info.o
	ar rcs libabin.a ${C_OBJS} $(F_OBJS) compile_info.o

# compile_info.F90 must be always recompiled to get the current date/time and git commit
compile_info.o: compile_info.F90 ${C_OBJS} ${F_OBJS} abin.o
	$(FC) $(FFLAGS) $(DFLAGS) -DGIT_COMMIT="'${COMMIT}'" -c $<

clean:
	/bin/rm -f *.o *.mod *.gcno *.gcda *.gcov libabin.a $(BIN)

.PHONY: clean

%.o: %.F90
	$(FC) $(FFLAGS) $(DFLAGS) -c $<

%.o: %.f
	$(FC) $(FFLAGS) $(DFLAGS) -c $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c $<
