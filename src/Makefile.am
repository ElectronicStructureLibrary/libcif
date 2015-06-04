#
# Copyright (C) 2015 Y. Pouillon
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 3 of the License, or 
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
#

#
# Makefile for the sources of LibCIF
#

                    # ------------------------------------ #

#
# Source files
#

# C source - keep this in alphabetical order
cif_core_srcs = \
  dummy.c

# Exported C headers - keep this in alphabetical order
cif_core_hdrs = \
  matrix3.h \
  vector3.h

# Internal C headers - keep this in alphabetical order
cif_hidden_hdrs =

                    # ------------------------------------ #

#
# Build targets
#

# Essential build parameters
AM_CPPFLAGS = @cif_core_incs@
AM_CFLAGS = @CFLAGS_COVERAGE@
AM_LDFLAGS = @LDFLAGS_COVERAGE@

# Libraries
lib_LTLIBRARIES = libcif.la

libcif_la_LDFLAGS = -version-info 0:0:0
libcif_la_SOURCES = $(cif_core_srcs)

# Headers
include_HEADERS = $(cif_core_hdrs)
noinst_HEADERS = $(cif_hidden_hdrs)

# Fortran
#if DO_BUILD_FORTRAN
#  SUBDIRS = . fortran
#endif

                    # ------------------------------------ #

#
# Test programs
#

cif_basic_tests = \
  test_error
#cif_xfail_tests = \
#  test_io
cif_gcov_pre_tests = \
  gcov_clean
cif_gcov_post_tests = \
  gcov_report
cif_memprof_tests = \
  test_memory

check_PROGRAMS = $(cif_basic_tests) $(cif_xfail_tests)
nodist_check_SCRIPTS = \
  $(cif_gcov_pre_tests) \
  $(cif_gcov_post_tests) \
  $(cif_memprof_tests)

test_error_SOURCES = test_error.c
test_error_LDADD = -lcif $(LIBS_COVERAGE)
test_error_CPPFLAGS = -I$(srcdir)/../src
test_error_DEPENDENCIES = libcif.la

TESTS = \
  $(cif_gcov_pre_tests) \
  $(cif_basic_tests) \
  $(cif_xfail_tests) \
  $(cif_gcov_post_tests) \
  $(cif_memprof_tests)
#XFAIL_TESTS = \
#  $(cif_xfail_tests)

                    # ------------------------------------ #

#
# Explicit dependencies
#

-include gcov-tests-deps.mf

                    # ------------------------------------ #

#
# Clean targets
#

# Temporary test outputs
MOSTLYCLEANFILES = *.gcda *.gcno *.gcov test_*.tmp
