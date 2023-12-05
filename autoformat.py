#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
#    This file is part of fprettify.
#    Copyright (C) 2016-2019 Patrick Seewald, CP2K developers group
#
#    fprettify is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    fprettify is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with fprettify. If not, see <http://www.gnu.org/licenses/>.
###############################################################################

# Vendored to ABIN by Daniel Hollas, with minor modifications.
# The source code in fprettify/ was taken from the 'abin' branch of the fork at
# https://github.com/danielhollas/fprettify

"""Wrapper script to run fprettify.
   When no arguments are provided, we run on all Fortran files in src/
"""

import sys
from fprettify import run
try:
    import configargparse
except ImportError:
    print("Could not find configargparse package")
    print("Please install it e.g. via 'pip install --user configargparse")
    sys.exit(1)

if len(sys.argv) == 1:
    print("Autoformatting files in src/")
    sys.argv.append('-r')
    sys.argv.append('src/')

if __name__ == '__main__':
    run()
