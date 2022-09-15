# This file is processed by autoconf to create a Makefile for the GARLI
# This file was written by hand (by mth using PAUP's configure.ac and http://autotoolset.sourceforge.net/tutorial.html as guides).

AC_PREREQ(2.59)
AC_INIT([Garli], [2.2], [garli.support@gmail.com])

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

# A few miscelaneous features, not of general interest
AC_ARG_ENABLE(profiler, AC_HELP_STRING([--enable-profiler], [built in runtime profiler]),
[AC_DEFINE([ENABLE_CUSTOM_PROFILER], [1], [profiler for assessing execution times])],
[])

AC_ARG_ENABLE(leak_detection, AC_HELP_STRING([--enable-leak-detection], [memory leak detection is turned on]),
[AC_DEFINE([MONITORING_ALLOCATION], [1], [monitors all allocation/deallocation, writes report])],
[])


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

# Compile BOINC version
AC_ARG_ENABLE(boinc,
        AC_HELP_STRING([--enable-boinc],[build BOINC version (requires BOINC library and api installation.  Building may require much fiddling.)]), , [enable_boinc=no])
if test "$enable_boinc" = yes; then
        AC_MSG_NOTICE([
           *** NOTE: compiling BOINC version
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

# using hmsbeagle library 
AC_ARG_ENABLE(beagle,
        AC_HELP_STRING([--enable-beagle],[use hmsbeagle calculation library (EXPERIMENTAL)]), , [enable_beagle=no])
if test "$enable_beagle" = yes; then
        AC_MSG_NOTICE([
           *** NOTE: compiling version using hmsbeagle calculation library
           ])
fi
AM_CONDITIONAL([BEAGLE], [test "x$enable_beagle" = xyes] )

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
	elif test "$CC" = "gcc"; then
		CXXFLAGS="$CXXFLAGS -fopenmp"
	else
		AC_MSG_ERROR([the --enable-openmp option can only be used with the intel compiler and newer versions of gcc (CC=icc or CC=gcc)])
	fi
fi

#BOINC compilation
if test  "$enable_boinc" = yes; then
        CXXFLAGS="$CXXFLAGS -DBOINC"
		LIBS="$LIBS -lpthread -lboinc_api -lboinc"
fi

#Compiler choice
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

#beagle library support - this is a crappy temporary hack
AC_ARG_WITH(
	[beagle],
	AC_HELP_STRING(
		[--with-beagle=DIR],
		[Specify the root directory for the hmsbeagle library installation. This is required to compile the beagle version with the --enable-beagle flag.]
		),
	[
	if ! test "$withval" = "yes" -o "$withval" = "no" ; then
		BEAGLE_INSTALL_PATH="$withval"
	fi
	])
		
if test  "$enable_beagle" = yes; then
   BEAGLE_LDFLAGS="-L$BEAGLE_INSTALL_PATH/lib"
   BEAGLE_LIB="-lhmsbeagle"
   BEAGLE_INCLUDES="-I$BEAGLE_INSTALL_PATH/include/libhmsbeagle-1 -I$BEAGLE_INSTALL_PATH/include/libhmsbeagle-1/libhmsbeagle"
   LDFLAGS="$LDFLAGS $BEAGLE_LDFLAGS"
   LIBS="$LIBS $BEAGLE_LIB"
   CXXFLAGS="$CXXFLAGS $BEAGLE_INCLUDES -DNEW_MANAGEMENT -DUSE_BEAGLE -DBEAGLEPART"
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

#Location of BOINC library and include files
BOINC_INSTALL="/usr"
AC_ARG_WITH(
	[boinc],
	AC_HELP_STRING(
		[--with-boinc=DIR],
		[Specify the root directory for the BOINC library (parent of the include/boinc and lib directories). Omit this flag to use a system level install of BOINC.]
		),        
	[
	if ! test "$withval" = "yes" -o "$withval" = "no" ; then
		BOINC_INSTALL="$withval"
	fi
	])

if test  "$enable_boinc" = yes; then
	BOINC_INCL="$BOINC_INSTALL/include/boinc"
	BOINC_LIB_DIR="$BOINC_INSTALL/lib"
	if ! test -d "$BOINC_INSTALL" ; then
		AC_MSG_ERROR([A BOINC library and api installation is required for building the BOINC version of GARLI. Expecting a directory called $BOINC_INCL and libboinc library files in $BOINC_LIB_DIR. Use the --with-boinc flag with configure to specify a location other than /usr/local, or verify that the specified location is correct.])
	fi
	CXXFLAGS="$CXXFLAGS -I$BOINC_INCL"
	LDFLAGS="$LDFLAGS -L$BOINC_LIB_DIR"
fi

################################################################################
#	Require the builder to provide the --with-ncl argument,
#	otherwise default to looking for a system level installation.
#	ideally people should be using the build_garli.sh install script
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

####DJZ This was to look for an NCL install in other places, but 
#I gave up on it because the user should either use the script
#or know what they are doing.
#if ! test -n "$WNCL" ; then
#	echo WITH_NCL NOT SPECIFIED
#	if test -d "$srcdir/ncl/" ; then
#		echo FOUND NCL DIR
#		nl=`ls -ld $srcdir/ncl*/ | wc -l`
#		echo ###		
#		echo $nl
#		echo ###
#		if ! test $nl = 1 ; then
#			AC_MSG_ERROR([Multiple NCL directories found.])
#		fi
#		NCL_BASE=`ls -d $srcdir/ncl*/`
#		NCL_INC_DIR="$NCL_BASE/include"
#		NCL_LIB_DIR="$NCL_BASE/lib"
#		echo NCL_INC_DIR $NCL_INC_DIR
#	else
#		echo NO NCL DIR FOUND	
#	fi
#else
#	echo WITH_NCL SPECIFIED
#fi
###

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

LIBS="$LIBS -lncl"

AC_TRY_LINK(
[],
[int foo=2;],
[],
[AC_MSG_ERROR([Problem with link flags.  Most likely a library cannot be found.  
Check that any libraries listed below can be found in a system default location or a location specified in LDFLAGS with a -L. 
LIBS=$LIBS
LDFLAGS=$LDFLAGS
])])

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
