#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ======================================================================
#

root="splitcube"

# ----------------------------------------------------------------------
def setupInput(nprocs):
  dupext = [".fuldat", ".prop", ".statevar", ".time"]
  sinext = [".coord", ".connect", ".bc", ".split"]

  print "Setting up symbolic links with prefix '%s_%d':" % (root, nprocs)
  import os

  dirFiles = os.listdir(os.getcwd())
  for ext in sinext:
    src = "%s%s" % (root, ext)
    dest = "%s_%s%s" % (root, nprocs, ext)
    if not dest in dirFiles:
      print "  %s -> %s... created" % (dest, src)
      os.symlink(src, dest)
    else:
      print "  %s -> %s... already exists" % (dest, src)

  for ext in dupext:
    src = "%s%s" % (root, ext)
    dest = "%s_%s%s" % (root, nprocs, ext)
    if not dest in dirFiles:
      print "  %s -> %s... created" % (dest, src)
      os.symlink(src, dest)
    else:
      print "  %s -> %s... already exists" % (dest, src)

  return

# ----------------------------------------------------------------------
if __name__ == "__main__":

  nprocs = 1
  setupInput(nprocs)
  nprocs = 2
  setupInput(nprocs)

# End of file
