# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file pylith/utils/profiling.py

# ----------------------------------------------------------------------
def resourceUsage():
    """Get CPU time (hh:mm:ss) and memory use (MB).
    """

    try:
        import os
        import commands
        cmd = "ps -p %d -o cputime,rss" % os.getpid()
        info = commands.getoutput(cmd).split()
        cputime = info[2]
        memory = float(info[3])/1024.0
    except:
        cputime = "n/a"
        memory = 0
    return (cputime, memory)


# ----------------------------------------------------------------------
def resourceUsageString():
    """Get CPU time and memory usage as a string.
    """
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()
    (cputime, memory) = resourceUsage()
    return "[%d] CPU time: %s, Memory usage: %.2f MB" % \
        (comm.rank, cputime, memory)


# End of file
