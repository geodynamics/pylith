#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

__all__ = ['PhysicalProperties',
           'Solution',
           'StateVariables',
           'Fault',
           ]


# ----------------------------------------------------------------------
def has_h5py():
  if not "flag" in dir(has_h5py):
    try:
      import h5py
      has_h5py.flag = True
    except ImportError:
      print "WARNING: Cannot find h5py Python modele."
      print "         Tests limited to running PyLith without errors."
      print "         Install h5py (available via the installer utility) "
      print "         in order to enable verification of output."
      has_h5py.flag = False
  return has_h5py.flag


# ----------------------------------------------------------------------
def run_pylith(appClass, dbClass=None, nprocs=1):
  """
  Helper function to generate spatial databases and run PyLith.
  """
  if not str(appClass) in dir(run_pylith):
    if not dbClass is None:
      # Generate spatial databases
      db = dbClass()
      db.run()

    # Run PyLith, limiting number of processes to number of local CPUs or maximum specified by environment.
    import multiprocessing
    numCPUs = multiprocessing.cpu_count()
    maxThreads = numCPUs
    import os
    if "MAX_PYLITH_PROCS" in os.environ:
      maxThreads = int(os.environ["MAX_PYLITH_PROCS"])
    if numCPUs < maxThreads:
      appNumProcs = min(numCPUs, nprocs)
      if appNumProcs < nprocs:
        print("WARNING: Detected %d CPUs. Reducing number of processes from %d to %d." % (maxNumProcs, nprocs, maxNumProcs))
    else:
        appNumProcs = min(maxThreads, nprocs)
        if appNumProcs < nprocs:
          print("WARNING: Detected environment with MAX_PYLITH_PROCS=%d. Reducing number of processes from %d to %d." % (appNumProcs, nprocs, appNumProcs))
      
    app = appClass()
    app.nodes = appNumProcs
    setattr(run_pylith, str(appClass), True)
    app.run()
  return


# End of file
