AUTOMAKE_OPTIONS = gnu
ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src tests examples

info_TEXINFOS = isd.texi
isd_TEXINFOS = fdl-1.3.texi

bench: all check
	cd tests && $(MAKE) $(AM_MAKEFLAGS) bench
