
from archimedes import use_merlin
use_merlin()

from merlin import setup, find_packages

setup(
    
    name = 'PyLith', 
    version = '1.0',

    zip_safe = False,
    packages = find_packages(),
    
    install_requires = [
    'spatialdata',
    'pythia[mpi] >= 0.8.1.0, < 0.8.2a',
    ],

    author = 'Charles A. Williams, Brad Aagaard, and Matt Knepley',
    author_email = 'cig-short@geodynamics.org',
    description = """A finite element code for the solution of visco-elastic/plastic deformation that was designed for lithospheric modeling problems.""",
    license = 'other',
    url = 'http://www.geodynamics.org/cig/software/packages/short/pylith/',

)
