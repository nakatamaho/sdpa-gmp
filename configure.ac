dnl
dnl Configure for SDPA-GMP 7.1.3
dnl

AC_INIT(sdpa-gmp, 7.1.3)
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE([subdir-objects])
AC_PROG_CC
AC_PROG_CXX
AC_CANONICAL_HOST
AC_OPENMP
AC_LANG([C++])
AC_OPENMP

AC_PROG_CC
IS_INTELCC=0
AC_COMPILE_IFELSE([
    AC_LANG_SOURCE([[
    #if !defined __INTEL_COMPILER
    #error
    #endif
    #include <stdio.h>
    int main() {
        return 0;
    }
   ]])
], [is_intelcc=1])
AM_CONDITIONAL(IS_INTELCC, test x$is_intelcc = x1)
if test x"$is_intelcc" = x"1"; then
    AC_MSG_RESULT([Intel One API (icc) detected])
fi

AC_PROG_CXX
IS_INTELCXX=0
AC_COMPILE_IFELSE([
    AC_LANG_SOURCE([[
    #if !defined __INTEL_COMPILER
    #error
    #endif
    #include <stdio.h>
    int main() {
        return 0;
    }
   ]])
], [is_intelcxx=1])
AM_CONDITIONAL(IS_INTELCXX, test x$is_intelcxx = x1)
if test x"$is_intelcxx" = x"1"; then
    AC_MSG_RESULT([Intel One API (icpc) detected])
fi

highest_abs_top_srcdir=`cd "$srcdir" && pwd`

AC_ARG_WITH(system-spooles,
[  --with-system-spooles      Use spooles already on system],
WITH_SYSTEM_SPOOLES=$withval, WITH_SYSTEM_SPOOLES=no)
AM_CONDITIONAL(BUILD_SPOOLES, test x$WITH_SYSTEM_SPOOLES = xno)

AC_ARG_WITH(spooles-includedir,
[  --with-spooles-includedir    Gives spooles include dir you use(e.g., /usr/local/include/spooles).],
SPOOLES_INCLUDEDIR=$withval
)

AC_ARG_WITH(spooles-libdir,
[  --with-spooles-libdir    Gives spooles lib dir you use(e.g., /usr/local/lib).],
SPOOLES_LIBDIR=$withval
)

AC_ARG_WITH(system-gmp,
[  --with-system-gmp      Use gmp already on system],
WITH_SYSTEM_GMP=$withval, WITH_SYSTEM_GMP=no)
AM_CONDITIONAL(BUILD_GMP, test x$WITH_SYSTEM_GMP = xno)

AC_ARG_WITH(gmp-includedir,
[  --with-gmp-includedir    Gives gmp include dir you use(e.g., /usr/local/include).],
GMP_INCLUDEDIR=$withval
)

AC_ARG_WITH(gmp-libdir,
[  --with-gmp-libdir    Gives gmp lib dir you use(e.g., /usr/local/lib).],
GMP_LIBDIR=$withval
)

AC_ARG_WITH(system-mplapack,
[  --with-system-mplapack      Use mplapack already on system],
WITH_SYSTEM_MPLAPACK=$withval, WITH_SYSTEM_MPLAPACK=no)
AM_CONDITIONAL(BUILD_MPLAPACK, test x$WITH_SYSTEM_MPLAPACK = xno)

AC_ARG_WITH(mplapack-includedir,
[  --with-mplapack-includedir    Gives mplapack include dir you use(e.g., /usr/local/include).],
MPLAPACK_INCLUDEDIR=$withval
)

AC_ARG_WITH(mplapack-libdir,
[  --with-mplapack-libdir    Gives mplapack lib dir you use(e.g., /usr/local/lib).],
MPLAPACK_LIBDIR=$withval
)

dnl spooles requires -lm when linking.
LIBS="-lm"

IS_LINUX=0
case "$host" in
    *linux*)
        IS_LINUX=1
        ;;
    *)
esac
AM_CONDITIONAL(IS_LINUX, test x$IS_LINUX = x1)

IS_MINGW=0
case "$host" in
        *-mingw32*) # mingw
          IS_MINGW=1
        ;;
        *)
esac
AM_CONDITIONAL(IS_MINGW, test x$IS_MINGW = x1)

IS_MACOS=0
case "$host" in
         *apple-darwin*) # MacOSX
          IS_MACOS=1
        ;;
        *)
esac
AM_CONDITIONAL(IS_MACOS, test x$IS_MACOS = x1)

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
SPOOLES_LIBDIR=$highest_abs_top_srcdir/external/i/SPOOLES/lib
SPOOLES_INCLUDEDIR=$highest_abs_top_srcdir/external/spooles/work/internal #just too messy
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

dnl Check for GMP include and library files
if test x$WITH_SYSTEM_GMP != xno; then
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
else
  GMP_LIBDIR=$highest_abs_top_srcdir/external/i/GMP/lib
  GMP_INCLUDEDIR=$highest_abs_top_srcdir/external/i/GMP/include
fi
AC_SUBST(GMP_LIBDIR)
AC_SUBST(GMP_INCLUDEDIR)

dnl Check for MPLAPACK include and library files
dnl some dummy settings
if test x$MPLAPACK_INCLUDEDIR = x; then
MPLAPACK_INCLUDEDIR="."
fi
if test x$MPLAPACK_LIBDIR = x; then
MPLAPACK_LIBDIR="."
fi

if test x$WITH_SYSTEM_MPLAPACK != xno; then
HAVE_MPLAPACK_HEADER=""
HAVE_MPLAPACK=""
save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS -I$MPLAPACK_INCLUDEDIR"
AC_LANG_PUSH([C++])
AC_CHECK_HEADER(mpblas_dd.h, HAVE_MPLAPACK_HEADER=yes)
CPPFLAGS=$save_CPPFLAGS

if test x$HAVE_MPLAPACK_HEADER = xyes; then
  save_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS -L$MPLAPACK_LIBDIR"
  LIBS="$LDFLAGS -L$MPLAPACK_LIBDIR -lmplapack_dd -lmpblas_dd_opt"
 # AC_CHECK_FUNC(Rsyev, [HAVE_MPLAPACK=yes])
  HAVE_MPLAPACK=yes
  if test x$IS_LINUX = x1; then
    LDFLAGS="$save_LDFLAGS -Wl,--rpath=$MPLAPACK_LIBDIR"
  fi
fi

if test x$HAVE_MPLAPACK != xyes; then
  echo "Please install MPLAPACK (MPBLAS/MPLAPACK)"
  echo "https://github.com/nakatamaho/mplapack/"
  AC_MSG_ERROR([No MPLAPACK found])
fi
else
MPLAPACK_LIBDIR=$highest_abs_top_srcdir/mplapack
MPLAPACK_INCLUDEDIR=$highest_abs_top_srcdir/mplapack
fi
AM_CONDITIONAL(SYSTEM_MPLAPACK, test x$with_system_mplapack = xyes)
AC_SUBST(MPLAPACK_LIBDIR)
AC_SUBST(MPLAPACK_INCLUDEDIR)

AC_CONFIG_FILES([
external/Makefile
external/spooles/Makefile
external/gmp/Makefile
])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
