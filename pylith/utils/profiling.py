#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ----------------------------------------------------------------------
#

## @file pylith/utils/profile.py

# ----------------------------------------------------------------------
def resourceUsage():
  """
  Get CPU time (hh:mm:ss) and memory use (MB).
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
  """
  Get CPU time and memory usage as a string.
  """
  return "CPU time: %s, Memory usage: %.2f MB" % resourceUsage()


# End of file
