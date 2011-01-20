# This file is processed by autoconf to create a Makefile for the GARLI
# This file was written by hand (by mth using PAUP's configure.ac and http://autotoolset.sourceforge.net/tutorial.html as guides).

AC_PREREQ(2.59)
AC_INIT([Garli], [0.97], [zwickl@nescent.org])

AC_CONFIG_SRCDIR([src/garlimain.cpp])
AC_CONFIG_HEADER([config.h])

#	Directory that contains install-sh and other auxiliary files
AC_CONFIG_AUX_DIR([config])

################################################################################
# Checks for build-platform and target info
#	this defines the "target" variable that is used later in this file
################################################################################
AC_CANONICAL_TARGET

################################################################################
#	According to (http://www.mail-archive.com/autoconf@gnu.org/msg14232.html)
#		this macro should be after AC_INIT but before AM_INIT_AUTOMAKE
################################################################################
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.9 foreign dist-tarZ tar-ustar filename-length-max=299])

# Checks for programs.
AC_PROG_LN_S
AC_PROG_CC
AC_PROG_CPP
INP_CXXFLAGS="$CXXFLAGS"
AC_PROG_CXX

# Initialize CXXFLAGS to prevent it from defaulting to "-g -O2"
CXXFLAGS="$INP_CXXFLAGS -DUNIX -DNCL_CONST_FUNCS"

# required because we are building a library 
AC_PROG_RANLIB


AC_PROG_INSTALL
#AC_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h malloc.h stddef.h stdlib.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memmove memset pow sqrt strchr strdup strtol])



#--------------------------------------------------------------------------------------------------#
#  Set CXXFLAGS                                                                                    #
#--------------------------------------------------------------------------------------------------#


# Check for debugging mode.
AC_ARG_ENABLE(debugging,
	AC_HELP_STRING([--enable-debugging],[build for debugging]), , [enable_debugging=no])
if test "$enable_debugging" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: debugging is enabled; optimization is suppressed!
	   ])
fi


# we tell NCL that we will ignore asserts
NCL_ASSERT_FLAG="-DIGNORE_NXS_ASSERT"
# Check whether asserts should be allowed.
AC_ARG_ENABLE(asserts,
	AC_HELP_STRING([--enable-asserts],[build with asserts on (NDEBUG not defined)]), , [enable_asserts=no])
if test "$enable_asserts" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: compiling with assertions on (NDEBUG not defined)
	   ])
	if test "$enable_debugging" = yes; then
		NCL_ASSERT_FLAG=""
	fi
fi

# Compile openMP multithreaded version
AC_ARG_ENABLE(openmp,
	AC_HELP_STRING([--enable-openmp],[build OpenMP multithreaded version]), , [enable_openmp=no])
if test "$enable_openmp" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: compiling OpenMP multithreaded version
	   ])
fi

# Mpi run forking version
AC_ARG_ENABLE(mpi,
	AC_HELP_STRING([--enable-mpi],[build MPI run distributing version]), , [enable_mpi=no])
if test "$enable_mpi" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: compiling MPI run distributing version
	   ])
fi

# old Mpi wrapper version - run series of configs named run0.conf, run1.conf, etc
AC_ARG_ENABLE(oldmpi,
        AC_HELP_STRING([--enable-oldmpi],[build old MPI batch run version (YOU PROBABLY DON'T WANT THIS)]), , [enable_oldmpi=no])
if test "$enable_oldmpi" = yes; then
        AC_MSG_NOTICE([
           *** NOTE: compiling old MPI batch run version
           ])
fi

# single precision version (still being tested)
AC_ARG_ENABLE(single-prec,
        AC_HELP_STRING([--enable-single-prec],[use single precision floating point variables (EXPERIMENTAL)]), , [enable_single_prec=no])
if test "$enable_single_prec" = yes; then
        AC_MSG_NOTICE([
           *** NOTE: compiling single precision floating point version
           ])
fi


# Initialize optimization flag in case it doesn't get set below.
CXXFLAGS_OPTIM_SPEED="-O"

# "-g" may not work with some compilers, but end users shouldn't be
if test "$enable_debugging" = yes; then
	CXXFLAGS_OPTIM_SPEED="-O0"
	CXXFLAGS="$CXXFLAGS -Wall -g -Wreturn-type -Wunused -Wredundant-decls -Wcast-align -Wcomment -Wextra"
fi

#DJZ allow turning off NDEBUG only
if test  "$enable_asserts" = no; then
	CXXFLAGS="$CXXFLAGS -DNDEBUG"
fi

#DJZ OpenMP compilation
if test  "$enable_openmp" = yes; then
	if test "$CC" = "icc"; then
		CXXFLAGS="$CXXFLAGS -openmp"
	else
		AC_MSG_ERROR([the --enable-openmp option can only be used with the intel compiler (CC=icc)])
	fi
fi

if test "$CC" = "icc" -o "$CC" = "icc" ; then
#	Intel C compiler for Linux
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-O2 -ip -funroll-loops -fno-alias"
	fi
elif test "$CC" = "ccc"; then
#	Compaq C compiler for Linux
	if test "x$arch" = "x"; then
		arch="host"
	fi
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-fast -inline speed -arch $arch"
	fi
elif test "$CC" = "xlc"; then
#	IBM XL C compiler
	CCFLAGS="$CXXFLAGS -qsourcetype=c++ -qenablevmx -qchars=signed"
	if test "x$arch" = "x"; then
		arch="auto"
	fi
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-O3 -qarch=$arch -qtune=$arch -qalias=ansi -qunroll=yes"
	fi
elif test "x$GCC" = "xyes" ; then
	CXXFLAGS="$CXXFLAGS -Wno-uninitialized"
	if test "$enable_debugging" = yes; then
		CXXFLAGS_OPTIM_SPEED="-O0 -Wimplicit"
	else
		CXXFLAGS_OPTIM_SPEED="-O3 -ffast-math -funroll-loops -fstrict-aliasing"
	fi
	case "$build_os" in
	 	darwin*) CXXFLAGS="$CXXFLAGS" ;;
	 	*) CXXFLAGS="$CXXFLAGS -fsigned-char";;
	esac	
fi

#single precision
if test  "$enable_single_prec" = yes; then
	CXXFLAGS="$CXXFLAGS -DSINGLE_PRECISION_FLOATS"
	if test "$CC" = "gcc" ; then
		CXXFLAGS="$CXXFLAGS -fno-caller-saves"
	fi
fi

AC_LANG(C++)
ACX_MPI

#MPI run forker
if test  "$enable_mpi" = yes; then
	CXXFLAGS="$CXXFLAGS -DSUBROUTINE_GARLI"
	CC="$MPICC"
	CXX="$MPICXX"
	LIBS="$MPILIBS $LIBS"
fi

#old MPI versionr
if test  "$enable_oldmpi" = yes; then
	CXXFLAGS="$CXXFLAGS -DOLD_SUBROUTINE_GARLI"
	CC="$MPICC"
	CXX="$MPICXX"
	LIBS="$MPILIBS $LIBS"
fi


CXXFLAGS="$CXXFLAGS $CXXFLAGS_OPTIM_SPEED"


################################################################################
#	By default, use the bundled ncl code from source, 
#		but let the --with-ncl argument override this.
################################################################################
NCL_INC_DIR="/usr/local/include"
NCL_LIB_DIR="/usr/local/lib"
NCL_BIN_DIR="/usr/local/bin"

AC_ARG_WITH(
	[ncl],
	AC_HELP_STRING(
		[--with-ncl=DIR],
		[Specify the root directory for the ncl version 2.1 or greater library (parent of the include/ncl and lib directories).  Omit this flag and configure will look for a system level NCL installation.  You might use the build_garli.sh script to automate the build process.]
		),
	[
	if ! test "$withval" = "yes" -o "$withval" = "no" ; then
		NCL_INC_DIR="$withval/include"
		NCL_LIB_DIR="$withval/lib"
		NCL_BIN_DIR="$withval/bin"
	fi
	])

if ! test -d "$NCL_INC_DIR/ncl" ; then
	AC_MSG_ERROR([NCL 2.1 or higher is a prerequisite for building Garli. Expecting a directory called $NCL_INC_DIR/ncl, verify that the argument to --with-ncl is correct.])
else
	if ! test -f "$NCL_INC_DIR/ncl/nxsdefs.h" ; then
		AC_MSG_ERROR([NCL 2.1 or higher is a prerequisite for building Garli. Expecting a directory called $NCL_INC_DIR/ncl, verify that the argument to --with-ncl is correct.])
	fi
fi

if ! test -d "$NCL_LIB_DIR" ; then
	AC_MSG_ERROR([NCL 2.1 or higher is a prerequisite for building Garli. Expecting a directory called $NCL_LIB_DIR, verify that the argument to --with-ncl is correct.])
fi

CPPFLAGS="$CPPFLAGS -I$NCL_INC_DIR -I$NCL_INC_DIR/ncl $NCL_ASSERT_FLAG"
if test -d "$NCL_LIB_DIR/ncl" ; then
	LDFLAGS="$LDFLAGS -L$NCL_LIB_DIR/ncl"
else
	LDFLAGS="$LDFLAGS -L$NCL_LIB_DIR"
fi

AC_TRY_LINK(
[],
[int foo=2;],
[],
[AC_MSG_ERROR([Problem with link flags.  Most likely a library cannot be found.  
Check that any libraries listed below can be found in a system default location or a location specified in LDFLAGS with a -L. 
LIBS=$LIBS
LDFLAGS=$LDFLAGS
])])

LIBS="$LIBS -lncl"
AC_TRY_LINK( 
[#include "ncl/nxsdefs.h" ], 
[onlyDefinedInCouldBeConst ();], 
[], 
[AC_MSG_ERROR([NCL must be compiled with the NCL_CONST_FUNCS macro defined. Usually adding -DNCL_CONST_FUNCS to CPPFLAGS in your invocation of the NCL configure script is the solution.])]) 

AC_COMPILE_IFELSE( 
[#include "ncl/nxsdefs.h" 
# if ! defined(NCL_MAJOR_VERSION) ||  NCL_MAJOR_VERSION < 2
#	error "need at least v2"
# else
# 	 if ! defined(NCL_MINOR_VERSION) ||  NCL_MINOR_VERSION < 1
#		error "need at least v2.1"
#	 endif
# endif
], 
[], 
[AC_MSG_ERROR([Sorry. NCL 2.1 or greater is required.])]) 

GARLI_LIBS="-lncl"

AC_SUBST([NCL_BIN_DIR])

AC_SUBST([GARLI_LIBS])

AC_SUBST([CPPFLAGS])
AC_SUBST([LDFLAGS])


#currently we are bundling NCL in the svn, so we know the path

# Variables to substitute in Makefile:
AC_SUBST(CPPFLAGS)

AC_CONFIG_FILES([
	Makefile 			\
	src/Makefile		\
	tests/Makefile	\
	])

AC_OUTPUT
