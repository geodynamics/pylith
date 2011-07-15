# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

from archimedes import use_merlin
use_merlin()

from merlin import setup, find_packages

setup(
    
    name = 'PyLith', 
    version = '1.6.2',

    zip_safe = False,
    packages = find_packages(),
    
    install_requires = [
    'spatialdata',
    'pythia[mpi] >= 0.8.1.6, < 0.8.2a',
    ],

    author = 'Brad Aagaard, Charles A. Williams, and Matt Knepley',
    author_email = 'cig-short@geodynamics.org',
    description = """A finite element code for the solution of crustal deformation problems associated with earthquakes, including kinematic ruptures, spontaneous ruptures, and pre- and post-seismic deformation.""",
    license = 'other',
    url = 'http://www.geodynamics.org/cig/software/packages/short/pylith/',

)
