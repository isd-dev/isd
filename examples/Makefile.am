INCLUDES = -I$(top_srcdir)/src -I$(top_builddir)/src
AM_CPPFLAGS = -ansi -Wall -Wshadow

LDADD = $(top_builddir)/src/libisd.la

findlowweight_SOURCES = findlowweight.c cmdline.c cmdline.h

EXTRA_PROGRAMS = findlowweight

EXTRA_DIST = script.sh test_320_160 findlowweight.ggo

CLEANFILES = $(EXTRA_PROGRAMS)