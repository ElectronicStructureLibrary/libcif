#!/bin/sh
#
# Copyright (C) 2012-2015 Y. Pouillon
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

# Note: this script is temporary and will be removed upon release.

# Stop at first error and echo commands
set -ev

# Check that we are in the correct directory
test -s "configure.ac" -a -s "src/vector3.h" || exit 0

# Init build parameters
export CC="gcc"
export CFLAGS="-O0 -g3 -ggdb -Wall -Wextra -fbounds-check -fno-inline"
export FC="gfortran"
export FCFLAGS="-O0 -g3 -ggdb -Wall -Wextra -fbounds-check -fno-inline"

# Prepare source tree
./wipeout.sh
./autogen.sh

# Check default build
mkdir tmp-minimal
cd tmp-minimal
../configure
make dist
make
make clean && make -j4
make check
mkdir install-minimal
make install DESTDIR="${PWD}/install-minimal"
ls -lR install-minimal >install-minimal.log
cd ..

# Make distcheck
mkdir tmp-distcheck
cd tmp-distcheck
../configure
make distcheck -j4
make distcleancheck

# Clean-up the mess
rm -rf tmp-minimal tmp-distcheck
