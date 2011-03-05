AM_CPPFLAGS = -ansi -pedantic -Wall

LDADD = $(top_builddir)/src/libisd.la

check_HEADERS = test.h
check_PROGRAMS = systematic_io systematic_swap isd_run isd_c_api

systematic_io_SOURCES = systematic_io.cc
systematic_swap_SOURCES = systematic_swap.cc
isd_run_SOURCES = isd_run.cc

isd_c_api_SOURCES =  isd_c_api.c

EXTRA_PROGRAMS = isd_bench_1 isd_bench_2

isd_bench_1_SOURCES = isd_bench_1.cc
isd_bench_1_CPPFLAGS = -DISD_BENCH -ansi -Wall -Werror
isd_bench_2_SOURCES = isd_bench_2.cc
isd_bench_2_CPPFLAGS = -DISD_BENCH -ansi -Wall -Werror

CLEANFILES = $(EXTRA_PROGRAMS)

bench: test.h isd_bench_1$(EXEEXT) isd_bench_2$(EXEEXT)
	./isd_bench_1
	./isd_bench_2

TESTS = $(check_PROGRAMS)