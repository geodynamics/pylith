
from archimedes import use_merlin
use_merlin()

from merlin import setup, find_packages

setup(
    
    name = 'pylith3d', 
    version = '0.8.3',

    zip_safe = False,
    package_dir = { "":"pylith3d" },
    packages = find_packages("pylith3d"),
    
    install_requires = [
    'pythia[mpi] >= 0.8.1.3, < 0.8.2a',
    ],

    author = 'Charles A. Williams, Brad Aagaard, and Matt Knepley',
    author_email = 'cig-short@geodynamics.org',
    description = """A finite element code for the solution of visco-elastic/plastic deformation that was designed for lithospheric modeling problems.""",
    license = 'other',
    url = 'http://www.geodynamics.org/cig/software/packages/short/pylith/',

)
