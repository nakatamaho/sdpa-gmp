bin_PROGRAMS = sdpa_gmp

if SYSTEM_MPLAPACK
sdpa_gmp_SOURCES = sdpa_chordal.cpp sdpa_dataset.cpp sdpa_io.cpp \
sdpa_jordan.cpp sdpa_linear.cpp sdpa_main.cpp sdpa_newton.cpp sdpa_parts.cpp \
sdpa_struct.cpp sdpa_tool.cpp sdpa_chordal.h sdpa_dataset.h sdpa_algebra.h \
sdpa_include.h sdpa_io.h sdpa_jordan.h sdpa_linear.h sdpa_newton.h \
sdpa_parts.h sdpa_right.h sdpa_struct.h sdpa_tool.h
else
sdpa_gmp_SOURCES = sdpa_chordal.cpp sdpa_dataset.cpp sdpa_io.cpp \
sdpa_jordan.cpp sdpa_linear.cpp sdpa_main.cpp sdpa_newton.cpp sdpa_parts.cpp \
sdpa_struct.cpp sdpa_tool.cpp sdpa_chordal.h sdpa_dataset.h sdpa_algebra.h \
sdpa_include.h sdpa_io.h sdpa_jordan.h sdpa_linear.h sdpa_newton.h \
sdpa_parts.h sdpa_right.h sdpa_struct.h sdpa_tool.h \
mplapack/Rtrsm.cpp mplapack/Rsyrk.cpp\
mplapack/Rtrmm.cpp mplapack/Rtrsv.cpp mplapack/iMlaenv.cpp\
mplapack/Rlamch.cpp mplapack/Rlascl.cpp mplapack/Rsytrd.cpp mplapack/Rsterf.cpp\
mplapack/Rorgtr.cpp mplapack/Rlatrd.cpp mplapack/Rsyr2k.cpp mplapack/Rsytd2.cpp \
mplapack/Rlanst.cpp mplapack/Rlae2.cpp mplapack/Rlapy2.cpp mplapack/Rlasrt.cpp \
mplapack/Rorgql.cpp mplapack/Rorgqr.cpp mplapack/Rsymv.cpp mplapack/Rlarfg.cpp \
mplapack/Rsyr2.cpp mplapack/Rlassq.cpp mplapack/Rorg2l.cpp mplapack/Rlarft.cpp \
mplapack/Rlarfb.cpp mplapack/Rorg2r.cpp mplapack/Rnrm2.cpp mplapack/Rlarf.cpp \
mplapack/Rger.cpp mplapack/Rpotf2.cpp mplapack/Rpotrf.cpp mplapack/Rpotrf2.cpp \
mplapack/iMparmq.cpp\
mplapack/iMladlr.cpp\
mplapack/iMladlc.cpp\
mplapack/iMieeeck.cpp\
mplapack/Mxerbla.cpp\
mplapack/Rcombssq.cpp mplapack/Risnan.cpp\
mplapack/Mlsame.cpp mplapack/Rscal.cpp mplapack/Rgemv.cpp \
mplapack/Rtrmv.cpp mplapack/Rsteqr.cpp mplapack/Rlaset.cpp mplapack/Rlaev2.cpp \
mplapack/Rlasr.cpp mplapack/Rlartg.cpp mplapack/Rswap.cpp mplapack/Rsyev.cpp\
mplapack/Rlansy.cpp mplapack/Rlaisnan.cpp\
mplapack/Raxpy.cpp\
mplapack/Raxpy_omp.cpp\
mplapack/Rcopy.cpp\
mplapack/Rcopy_omp.cpp\
mplapack/Rdot.cpp\
mplapack/Rdot_serial.cpp\
mplapack/Rdot_omp.cpp\
mplapack/Rgemm_omp.cpp\
mplapack/Rgemm_NN_omp.cpp\
mplapack/Rgemm_NT_omp.cpp\
mplapack/Rgemm_TN_omp.cpp\
mplapack/Rgemm_TT_omp.cpp\
mplapack/mpblas_gmp.h\
mplapack/mplapack_config.h\
mplapack/mplapack_gmp.h\
mplapack/mplapack_utils_gmp.h
endif

AUTOMAKE_OPTIONS = subdir-objects

sdpa_gmp_LDFLAGS = -Wl,-rpath,$(GMP_LIBDIR)
sdpa_gmp_LDADD = $(pthread_lib) -L$(SPOOLES_LIBDIR) -lspooles -L$(GMP_LIBDIR) $(_OPENMP_CXXFLAGS) -L$(MPLAPACK_LIBDIR)
sdpa_gmp_CXXFLAGS= -g -O2 -march=native -fno-omit-frame-pointer $(pthread_cflags) -I$(SPOOLES_INCLUDEDIR) -DUseMETIS=0 -I$(GMP_INCLUDEDIR) $(_OPENMP_CXXFLAGS) -I$(MPLAPACK_INCLUDEDIR)
sdpa_gmp_LDADD +=-lgmpxx -lgmp

SUBDIRS = external

EXTRA_DIST = \
example1.dat \
example1.dat-s \
example1.ini \
example1.ini-s \
example2.dat \
param.sdpa \
external/gmp/download/*.tar.bz2 \
external/gmp/Makefile.am \
external/gmp/Makefile.in \
external/spooles/download/*.tgz \
external/spooles/patches/patch* \
external/spooles/Makefile.am \
external/spooles/Makefile.in
