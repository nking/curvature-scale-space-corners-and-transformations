## adapted from example's given by http://www.bioinf.uni-freiburg.de/~mmann/HowTo/automake.html

# initial information about the project
AC_INIT([GEVSimilarityParameters],[2.0],[nichole@climbwithyourfeet.com])

# check if the source folder is correct
AC_CONFIG_SRCDIR([src/main/c++/algorithms/curves/GEVSimilarityParameters.cpp])

# Checks for programs

# check for C++ compiler
AC_PROG_CXXCPP
AC_PROG_CXX
AC_PROG_RANLIB

# automake initialisation (mandatory) and check for minimal automake API version 1.11
AM_INIT_AUTOMAKE([1.11])

# use the C++ compiler for the following checks
AC_LANG([C++])

AC_ARG_WITH(boostinc, AS_HELP_STRING([--with-boostinc=PATH],
	[specify directory containing boost headers, e.g. /usr/include/c++/4.2.1/tr1]))

## remove '/tr1' from boost path
## -n $var returns 'true' if $var is non-zero
if test -n "$with_boostinc"; then
    boostinc_path=`echo "$with_boostinc" | sed -e 's/\/tr1//'`
    #boostinc_path=${with_boostinc%/tr1} 
    boostinc_include="-I$boostinc_path"
else
    echo "ERROR: set --with-boostinc=PATH\n      example: ./configure --with-boostinc=/usr/include/c++/4.2.1/tr1"
    exit 1 
fi


# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T


# distribute additional compiler and linker flags
## NOTE: to make files as input for coverage reports, add to AM_CXXFLAGS: -fprofile-arcs -ftest-coverage -lgcov
AC_SUBST(AM_CXXFLAGS, "$AM_CXXFLAGS -stdlib=libstdc++")
AC_SUBST(CPPFLAGS, "-Wall -W -pipe -g -Wno-unused-parameter -Wno-unused-variable -Wno-reorder $CPPFLAGS $boostinc_include")
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])


# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([string])
AC_CHECK_HEADERS([iostream])
AC_CHECK_HEADERS(tr1/unordered_map, [], [], 
    [[#include <tr1/unordered_map>
          ]])


# files to generate via autotools (.am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/main/c++/algorithms/curves/Makefile])
AC_CONFIG_FILES([src/test/c++/algorithms/curves/Makefile])

# generate the final Makefile etc.
AC_OUTPUT
