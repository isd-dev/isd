#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([ISD], [0.8dev], [no@mail.yet], [isd])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/isd.h])

# Checks for programs.
AC_PROG_SED

AC_PROG_CXXCPP
AC_PROG_CXX
AC_PROG_LIBTOOL

AC_LANG([C++])

# Checks for libraries.
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_rng_get])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_SIZEOF([unsigned long])

# Checks for library functions.

AM_INIT_AUTOMAKE([1.11 no-define])

AC_CANONICAL_HOST

# Cleanup, we don't use those MACROS (AC_DEFINE)
if test -f confdefs.h; then
   sed -i '/#define PACKAGE_/d' confdefs.h
   sed -i '/#define STDC_HEADERS/d' confdefs.h
   sed -i '/#define HAVE_DLFCN/d' confdefs.h
   sed -i '/#define HAVE_INTTYPES/d' confdefs.h
   sed -i '/#define HAVE_MEMORY/d' confdefs.h
   sed -i '/#define HAVE_STDLIB/d' confdefs.h
   sed -i '/#define HAVE_STDINT/d' confdefs.h
   sed -i '/#define HAVE_STRING/d' confdefs.h
   sed -i '/#define HAVE_SYS_/d' confdefs.h
   sed -i '/#define HAVE_UNISTD/d' confdefs.h
   sed -i '/#define HAVE_LIBM/d' confdefs.h
   sed -i '/#define HAVE_LIBGSLCBLAS/d' confdefs.h
   sed -i '/#define HAVE_LIBGSL/d' confdefs.h
fi

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile examples/Makefile])
AC_OUTPUT
