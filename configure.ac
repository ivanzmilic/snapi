#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(INVERT,1.01,noort@mps.mpg.de)
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_LANG([C++])
AC_PROG_CXX
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_SUBST(RANLIB)
COMPILER=$CXX
AC_SUBST(COMPILER)
AR=ar
AC_SUBST(AR)

# Checks for libraries.
AC_CHECK_LIB([m], [exp])
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([z], [compress2])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdio.h stdlib.h string.h pthread.h])

AC_MSG_CHECKING([if the compiler requires -pthread])
AC_COMPILE_IFELSE(AC_LANG_SOURCE([[
#include <pthread.h>
int main(void){
  return 0;
}
]]),
       [AC_MSG_RESULT([no])],
       [PTHREADFLAG="-pthread"
       AC_MSG_RESULT([yes])
       ])
AC_SUBST(PTHREADFLAG)

AC_MSG_CHECKING([for int64_t definition])
AC_RUN_IFELSE(AC_LANG_SOURCE([[
#include <sys/types.h>
typedef long int int64_t;
int main(void){
  return sizeof(int64_t)!=8;
}
]]),
       [AC_DEFINE([LONG_INT64_T],1,[Define if int64_t may be defined as a long int.])
        AC_MSG_RESULT([long int])],
       [AC_MSG_RESULT([long long int])])

AC_MSG_CHECKING([for socklen_t in unistd.h or netdb.h])
AC_RUN_IFELSE(AC_LANG_SOURCE([[
#include <unistd.h>
#include <netdb.h>
int main(void){
  socklen_t s;
  return 0;
}
]]),
       [AC_DEFINE([HAVE_SOCKLEN_T],1,[define if socklen_t exists.])
       AC_MSG_RESULT([yes])],
       [AC_MSG_RESULT([no])])

#
#  check for CFITSIO support
#

AC_ARG_WITH([cfitsio],
  [AS_HELP_STRING([--with-cfitsio],
    [path to CFITSIO library @<:@required@:>@])],
  [],
  [with_cfitsio=yes])

#check for CFITSIO library and header files
AS_IF([test "x$with_cfitsio" != "xno"], [
  OCXXFLAGS=$CXXFLAGS; OLDFLAGS=$LDFLAGS; OLIBS=$LIBS
  FITSLIBS=" -lcfitsio"
  AS_IF([test "x$with_cfitsio" != "xyes"],[
    FITSINCFLAGS+=" -I$with_cfitsio -I$with_cfitsio/include"
    FITSLDFLAGS+=" -L$with_cfitsio -L$with_cfitsio/lib -L$with_cfitsio/lib64"
    echo "manual config: $FITSINCFLAGS"
  ],[
    PKG_CHECK_MODULES([CFITSIO],[cfitsio >= 2.0],[
      FITSINCFLAGS=$CFITSIO_CFLAGS
      FITSLIBS=$CFITSIO_LIBS],[
      FITSLDFLAGS+=" -L/usr/lib -L/usr/lib64"
      FITSINCFLAGS+=" -I/usr/include"
    ])
  ])

  CXXFLAGS+=" $FITSINCFLAGS"
  LDFLAGS+=" $FITSLDFLAGS"
  LIBS="$OLIBS $FITSLIBS"

  AC_MSG_CHECKING([if CFITSIO functional])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <fitsio.h>]],[[
  int status=0;
  fitsfile *fptr;
  fits_open_file(&fptr,"bla.fits",READONLY,&status)]])],[AC_MSG_RESULT([yes])],[AC_MSG_ERROR([no])])
  AC_SUBST(FITSINCFLAGS)
  AC_SUBST(FITSLDFLAGS)
  AC_SUBST(FITSLIBS)

  CXXFLAGS=$OCXXFLAGS; LDFLAGS=$OLDFLAGS; LIBS=$OLIBS],
  [echo no]
)

#
# Check for optimized FFT libraries
#
AC_ARG_WITH([fftw3],
  [AS_HELP_STRING([--with-fftw3],
    [support for FFTW3 library @<:@default=check@:>@])],
  [],
  [with_fftw3=yes])
  
FFT_SRC=fft

#check for FFTW3 library and header files
if test "x$with_fftw3" != "xno" ; then
  OCXXFLAGS=$CXXFLAGS; OLDFLAGS=$LDFLAGS; OLIBS=$LIBS
  FFTW3LIBS=" -lfftw3"
  if test "x$with_fftw3" != "xyes" ; then
    FFTW3INCFLAGS+=" -I$with_fftw3 -I$with_fftw3/include"
    FFTW3LDFLAGS+=" -L$with_fftw3 -L$with_fftw3/lib -L$with_fftw3/lib64"
    echo "manual config: $FFTW3INCFLAGS"
  else
    PKG_CHECK_MODULES([fftw3],[fftw3 >= 3.2],[
# this actually doesn't work!?
      FFTW3INCFLAGS=$fftw3_CFLAGS
      FFTW3LIBS=$fftw3_LIBS],[
      FFTW3LDFLAGS+=" -L/usr/lib -L/usr/lib64"
      FFTW3INCFLAGS+=" -I/usr/include"
    ])
  fi

  CXXFLAGS+=" $FFTW3INCFLAGS"
  LDFLAGS+=" $FFTW3LDFLAGS"
  LIBS="$OLIBS $FFTW3LIBS"

  AC_MSG_CHECKING([if FFTW3 functional])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <fftw3.h>]],[[
  fftw_plan plan;
  plan=fftw_plan_dft_2d(64,64,0,0,FFTW_FORWARD,FFTW_MEASURE)]])],[AC_MSG_RESULT([yes])
                                                                                                                                                                        FFT_SRC=fft_fftw],[AC_MSG_ERROR([no])])
  AC_SUBST(FFTW3INCFLAGS)
  AC_SUBST(FFTW3LDFLAGS)
  AC_SUBST(FFTW3LIBS)

  CXXFLAGS=$OCXXFLAGS; LDFLAGS=$OLDFLAGS; LIBS=$OLIBS
fi

AC_CONFIG_FILES([Makefile
                 cmlib/Makefile
                 grid/Makefile
                 atmos/Makefile
                 atmos/plpar/Makefile
                 atmos/atmol/Makefile
                 atmos/atmol/atom/H/Makefile
                 atmos/atmol/atom/He/Makefile
                 atmos/atmol/atom/Makefile
                 atmos/atmol/mol/Makefile
                 jdel/Makefile
                 jstat/Makefile
                 jsub/Makefile
                 master/Makefile
                 sdel/Makefile
                 slave/Makefile
                 slave/nr/Makefile])
AC_OUTPUT
