AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I m4

# bin_SCRIPTS = mincdti_maps

bin_PROGRAMS = volperf

EXTRA_DIST = $(m4_files) \
	epm-header.in \
	example_MINC_Vector.mv

m4_files = m4/mni_REQUIRE_LIB.m4 m4/mni_REQUIRE_MNILIBS.m4	 \
	m4/mni_REQUIRE_OPENINVENTOR.m4				 \
	m4/mni_cxx_have_koenig_lookup.m4 m4/smr_CGAL_MAKEFILE.m4 \
	m4/smr_OPTIONAL_LIB.m4 m4/smr_REQUIRED_LIB.m4		 \
	m4/smr_WITH_BUILD_PATH.m4

volperf_SOURCES = volperf.c \
	perf_util.h perf_util.c \
	minc_vector_io.h minc_vector_io.c
