nl
dnl Configure for SDPA-GMP 7.1.2
dnl

AC_INIT(sdpa-gmp, 7.1.2)
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_PROG_CXX
AC_CANONICAL_HOST

AC_ARG_WITH(system-spooles,
[  --with-system-spooles      Use spooles already on system],
WITH_SYSTEM_SPOOLES=$withval, WITH_SYSTEM_SPOOLES=no)
AM_CONDITIONAL(BUILD_SPOOLES, test x$WITH_SYSTEM_SPOOLES = xno)

AC_ARG_WITH(system-mpack,
[  --with-system-mpack      Use mpack already on system],
WITH_SYSTEM_MPACK=$withval, WITH_SYSTEM_MPACK=no)
AM_CONDITIONAL(BUILD_MPACK, test x$WITH_SYSTEM_MPACK = xno)

AC_ARG_WITH(gmp-includedir,
[  --with-gmp-includedir    Gives gmp include dir you use(e.g., /usr/local/include).],
GMP_INCLUDEDIR=$withval
)

AC_ARG_WITH(gmp-libdir,
[  --with-gmp-libdir    Gives gmp lib dir you use(e.g., /usr/local/lib).],
GMP_LIBDIR=$withval
)

AC_ARG_WITH(mpack-includedir,
[  --with-mpack-includedir    Gives mpack include dir you use(e.g., /usr/local/include).],
MPACK_INCLUDEDIR=$withval
)

AC_ARG_WITH(mpack-libdir,
[  --with-mpack-libdir    Gives mpack lib dir you use(e.g., /usr/local/lib).],
MPACK_LIBDIR=$withval
)

AC_ARG_WITH(spooles-includedir,
[  --with-spooles-includedir    Gives spooles include dir you use(e.g., /usr/local/include/spooles).],
SPOOLES_INCLUDEDIR=$withval
)

AC_ARG_WITH(spooles-libdir,
[  --with-spooles-libdir    Gives spooles lib dir you use(e.g., /usr/local/lib).],
SPOOLES_LIBDIR=$withval
)

AC_ARG_ENABLE(metis,
[  --enable-metis    Build with metis],
[case "${enableval}" in
  yes) with_metis=yes ;;
  no)  with_metis=no  ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-metis) ;;
esac],[with_metis=no])
AM_CONDITIONAL(ENABLE_METIS, test x$with_metis = xyes)

AC_ARG_WITH(metis-includedir,
[  --with-metis-includedir    Gives metis include dir you use(e.g., /usr/local/include/metis).],
METIS_INCLUDEDIR=$withval
)

AC_ARG_WITH(metis-libdir,
[  --with-metis-libdir    Gives metis lib dir you use(e.g., /usr/local/lib).],
METIS_LIBDIR=$withval
)

AC_ARG_ENABLE(openmp,
[  --enable-openmp       Enables openmp
],,)

AC_ARG_WITH(openmplib,
[  --with-openmplib    Gives openmplib],
OPENMPLIB=$withval
)

AC_ARG_WITH(openmpflag,
[  --with-openmpflag    Gives openmpflag],
OPENMPFLAG=$withval
)

AC_ARG_WITH(openmpthreads,
[  --with-openmpthread   Gives openmpthreads],
OPENMPTHREADS=$withval
)

AC_ARG_ENABLE(shared,
[  --enable-shared       Enables shared link to the libraries (slower)
],,)


dnl spooles and metis requires -lm when linking.
LIBS="-lm"

dnl How thread we link thread library esp. needed for spooles.
case "$host" in
        *-*-linux*|*-*-uclinux*)
            pthread_cflags="-D_REENTRANT"
            pthread_lib="-lpthread"
            ;;
        *-*-freebsd*|*-*-dragonfly*)
            pthread_cflags=""
            pthread_lib="-pthread"
            ;;
        *-apple-darwin*)
           enable_shared="yes"
        ;;
esac
LIBS="$LIBS $pthread_lib"
CXXFLAGS="$CXXFLAGS $pthread_cflags"
AC_SUBST(pthread_lib)

CFLAGS="$CFLAGS -O2 -funroll-all-loops"
CXXFLAGS="$CXXFLAGS -O2 -funroll-all-loops"
AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)

dnl some dummy settings
if test x$SPOOLES_INCLUDEDIR = x; then
SPOOLES_INCLUDEDIR="."
fi
if test x$SPOOLES_LIBDIR = x; then
SPOOLES_LIBDIR="."
fi

dnl

AC_MSG_CHECKING([whether to enable metis])
if test "$with_metis" = "no"; then
  AC_MSG_RESULT([no])
  ENABLE_METIS=NO
  AC_SUBST(ENABLE_METIS)
else
  AC_MSG_RESULT([yes])
  ENABLE_METIS=YES
  AC_SUBST(ENABLE_METIS)

  if test x$METIS_INCLUDEDIR = x; then
  METIS_INCLUDEDIR="."
  fi
  if test x$METIS_LIBDIR = x; then
  METIS_LIBDIR="."
  fi

  dnl Check for METIS include and library files
  HAVE_METIS_HEADER=""
  HAVE_METIS=""
  save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS -I$METIS_INCLUDEDIR"
  AC_CHECK_HEADER(metis.h, HAVE_METIS_HEADER=yes)
  CPPFLAGS=$save_CPPFLAGS
  if test x$HAVE_METIS_HEADER = xyes; then
    save_LDFLAGS="$LDFLAGS"
    LDFLAGS="$LDFLAGS -L$METIS_LIBDIR"
    AC_CHECK_LIB(metis, METIS_NodeND, [METIS_LIBS="-lmetis"; HAVE_METIS=yes])
    LDFLAGS=$save_LDFLAGS
  fi

  if test x$HAVE_METIS != xyes; then
    echo "Please install METIS"
    echo "METIS http://glaros.dtc.umn.edu/gkhome/views/metis/index.html"
    AC_MSG_ERROR([No METIS found])
  fi
  AC_SUBST(METIS_LIBDIR)
  AC_SUBST(METIS_INCLUDEDIR)
fi

if test x$WITH_SYSTEM_SPOOLES != xno; then
  dnl Check for SPOOLES include and library files
  HAVE_SPOOLES_HEADER=""
  HAVE_SPOOLES=""
  save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS -I$SPOOLES_INCLUDEDIR"
  AC_CHECK_HEADER(SPOOLES.h, HAVE_SPOOLES_HEADER=yes)
  CPPFLAGS=$save_CPPFLAGS

  if test x$HAVE_SPOOLES_HEADER = xyes; then
    save_LDFLAGS="$LDFLAGS"
    LDFLAGS="$LDFLAGS -L$SPOOLES_LIBDIR"
    AC_CHECK_LIB(spooles, InpMtx_init, [SPOOLES_LIB="-lspooles"; HAVE_SPOOLES=yes])
    LDFLAGS=$save_LDFLAGS
  fi

  if test x$HAVE_SPOOLES != xyes; then
    echo "Please install SPOOLES"
    echo "SPOOLES http://www.netlib.org/linalg/spooles/spooles.2.2.html"
    AC_MSG_ERROR([No SPOOLES found])
  fi
else
SPOOLES_LIBDIR=`pwd`/spooles/build
SPOOLES_INCLUDEDIR=`pwd`/spooles/build
fi
AC_SUBST(SPOOLES_LIBDIR)
AC_SUBST(SPOOLES_INCLUDEDIR)

dnl Check for GMP include and library files
dnl some dummy settings
if test x$GMP_INCLUDEDIR = x; then
GMP_INCLUDEDIR="."
fi
if test x$GMP_LIBDIR = x; then
GMP_LIBDIR="."
fi

HAVE_GMP_HEADER=""
HAVE_GMP=""
save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS -I$GMP_INCLUDEDIR"
AC_CHECK_HEADER(gmp.h, HAVE_GMP_HEADER=yes)
CPPFLAGS=$save_CPPFLAGS

if test x$HAVE_GMP_HEADER = xyes; then
  save_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS -L$GMP_LIBDIR"
  AC_CHECK_LIB(gmp, __gmpf_init_set, [GMP_LIBS="-lgmpxx -lgmp"; HAVE_GMP=yes])
  LDFLAGS=$save_LDFLAGS
fi

if test x$HAVE_GMP != xyes; then
  echo "Please install GMP"
  echo "GMP http://gmplib.org/"
  AC_MSG_ERROR([No GMP found])
fi
AC_SUBST(GMP_LIBDIR)
AC_SUBST(GMP_INCLUDEDIR)

dnl Check for MPACK include and library files
dnl some dummy settings
if test x$MPACK_INCLUDEDIR = x; then
MPACK_INCLUDEDIR="."
fi
if test x$MPACK_LIBDIR = x; then
MPACK_LIBDIR="."
fi

if test x$WITH_SYSTEM_MPACK != xno; then
HAVE_MPACK_HEADER=""
HAVE_MPACK=""
save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS -I$MPACK_INCLUDEDIR"
AC_LANG_PUSH([C++])
AC_CHECK_HEADER(mblas_gmp.h, HAVE_MPACK_HEADER=yes)
CPPFLAGS=$save_CPPFLAGS

if test x$HAVE_MPACK_HEADER = xyes; then
  save_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS -L$MPACK_LIBDIR"
  LIBS="$LDFLAGS -L$MPACK_LIBDIR -lmlapack_gmp -lmblas_gmp -lgmpxx -lgmp";
 # AC_CHECK_FUNC(Rsyev, [HAVE_MPACK=yes])
  HAVE_MPACK=yes
  LDFLAGS=$save_LDFLAGS
fi

if test x$HAVE_MPACK != xyes; then
  echo "Please install MPACK (MBLAS/MLAPACK)"
  echo "MPACK (MBLAS/MLAPACK) http://mplapack.cvs.sourceforge.net/"
  AC_MSG_ERROR([No MPACK found])
fi
else
MPACK_LIBDIR=`pwd`/mpack
MPACK_INCLUDEDIR=`pwd`/mpack
fi
AM_CONDITIONAL(SYSTEM_MPACK, test x$with_system_mpack = xyes)
AC_SUBST(MPACK_LIBDIR)
AC_SUBST(MPACK_INCLUDEDIR)

dnl openmp
WITH_OPENMP=NO
if test -n "$enable_openmp" && test "$enable_openmp" != "no"; then
OPENMPFLAGS="$OPENMPFLAG -DENABLE_OPENMP -DNUM_OF_THREADS=$OPENMPTHREADS"
WITH_OPENMP=YES
fi
AC_SUBST(OPENMPLIB)
AC_SUBST(OPENMPFLAGS)

AM_CONDITIONAL(ENABLE_SHARED, test x$enable_shared = xyes)

AC_OUTPUT(Makefile)
