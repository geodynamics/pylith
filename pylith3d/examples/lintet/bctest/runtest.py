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

root="bctest"

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
    for iproc in range(nprocs):
      dest = "%s_%s.%d%s" % (root, nprocs, iproc, ext)
      if not dest in dirFiles:
        print "  %s -> %s... created" % (dest, src)
        os.symlink(src, dest)
      else:
        print "  %s -> %s... already exists" % (dest, src)
  return


# ----------------------------------------------------------------------
def run(nprocs):
  print "Running PyLith..."

  # TODO: Replace the use of launching via 'system' with use
  # of Leif's architecture independent utility.

  cmd = "mpirun -np %d `which pylith3dapp.py` " \
        "--typos=relaxed " \
        "--scanner.fileRoot=%s_%d " \
        "--scanner.asciiOutput=full " \
        "--scanner.ucdOutput=ascii " \
        "-log_summary -pc_type bjacobi -sub_pc_type ilu " \
        "-ksp_monitor -ksp_view -ksp_rtol 1e-09" % (nprocs, root, nprocs)
  import os
  print cmd
  os.system(cmd)
  return


# ----------------------------------------------------------------------
if __name__ == "__main__":
  from optparse import OptionParser

  parser = OptionParser()
  parser.add_option("-n", "--numprocs", dest="nprocs",
                    type="int", metavar="NPROCS",
                    help="Set number of processors.")
  (options, args) = parser.parse_args()
  if len(args) != 0:
    parser.error("Incorrent number of arguments.")

  nprocs = 1
  if not options.nprocs is None:
    nprocs = options.nprocs

  setupInput(nprocs)
  run(nprocs)


# End of file
