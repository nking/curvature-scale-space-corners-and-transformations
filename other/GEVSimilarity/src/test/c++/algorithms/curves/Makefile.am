AM_CXXFLAGS = -I$(top_srcdir)/src @AM_CXXFLAGS@

bin_PROGRAMS = \
    GEVSimilarityParameters_test \
    ParametersKey_test

GEVSimilarityParameters_test_SOURCES = \
    GEVSimilarityParameters_test.cpp

ParametersKey_test_SOURCES = \
    ParametersKey_test.cpp

GEVSimilarityParameters_test_LDADD = \
    $(top_builddir)/src/main/c++/algorithms/curves/libGEVSimilarityParameters.a \
    $(AM_LDFLAGS)

ParametersKey_test_LDADD = \
    $(top_builddir)/src/main/c++/algorithms/curves/libGEVSimilarityParameters.a \
    $(AM_LDFLAGS)