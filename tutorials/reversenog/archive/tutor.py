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

archive_dir="../archive"
root="bmrsnog"

# ----------------------------------------------------------------------
def assembleStrings(fragments):
  """Assemble a list of strings from a list of list of fragments."""

  values = fragments[0]
  for fragset in fragments[1:]:
    newvalues = []
    for value in values:
      for fragment in fragset:
        newvalues.append(value + fragment)
    values = newvalues
  return values


# ----------------------------------------------------------------------
def getStepFiles(step):

  files = {'input': [],
           'output': []}
  if "mesh" == step:
    fragments = [ [root], [".geo"] ]
    files['input'] += assembleStrings(fragments)

    fragments = [ [root], [".netgen" ] ]
    files['output'] += assembleStrings(fragments)
    
  elif "setup" == step:
    fragments = [ [root], [".netgen", ".fault.par", ".par", ".aux"] ]
    files['input'] += assembleStrings(fragments)

    fragments = [ [root],
                  [".coord", ".connect", ".split", ".bc", ".w01.wink",
                   ".1.fcoord", ".1.fbc"] ]
    files['output'] += assembleStrings(fragments)

  elif "run1" == step:
    fragments = [ [root], [".coord", ".connect", ".split", ".bc",
                           ".fuldat", ".prop", ".statevar", ".time"] ]
    files['input'] += assembleStrings(fragments)
    files['input'] += ["runbm.py"]

    root1 = root + "_1"
    fragments = [ [root1], [".coord", ".connect", ".split", ".bc"] ]
    files['output'] += assembleStrings(fragments)

    fragments = [ [root1], [".0"],
                  [".coord", ".connect", ".split", ".bc",
                   ".ascii", ".vtk",
                   ".fuldat", ".prop", ".statevar", ".time"] ]
    files['output'] += assembleStrings(fragments)
    
    fragments = [ [root1], [".0"],
                  [".gmesh", ".mesh", ".mesh.split"],
                  ["", \
                   ".time.00000", ".time.00010", ".time.00050", ".time.00100"],
                  [".inp"] ]
    files['output'] += assembleStrings(fragments)
  elif "viz1" == step:
    root1 = root + "_1"
    fragments = [ [root1], [".0"],
                  [".mesh", ".mesh.time.00010"], [".inp"] ]
    files['input'] += assembleStrings(fragments)

    fragments = [ [root1], [".0"],
                  ["t00010"], [".inp"] ]
    files['output'] += assembleStrings(fragments)
  elif "run2" == step:
    fragments = [ [root], [".coord", ".connect", ".split", ".bc",
                           ".fuldat", ".prop", ".statevar", ".time"] ]
    files['input'] += assembleStrings(fragments)
    files['input'] += ["runbm.py"]

    root1 = root + "_2"
    fragments = [ [root1], [".coord", ".connect", ".split", ".bc"] ]
    files['output'] += assembleStrings(fragments)

    fragments = [ [root1], [".0", ".1"],
                  [".coord", ".connect", ".split", ".bc",
                   ".ascii", ".vtk",
                   ".fuldat", ".prop", ".statevar", ".time"] ]
    files['output'] += assembleStrings(fragments)
    
    fragments = [ [root1], [".0", ".1"],
                  [".gmesh", ".mesh", ".mesh.split"],
                  ["", \
                   ".time.00000", ".time.00010", ".time.00050", ".time.00100"],
                  [".inp"] ]
    files['output'] += assembleStrings(fragments)
  elif "viz2" == step:
    root1 = root + "_2"
    fragments = [ [root1], [".0", ".1"],
                  [".mesh", ".mesh.time.00010"], [".inp"] ]
    files['input'] += assembleStrings(fragments)

    fragments = [ [root1], [".0", ".1"],
                  [".t00010"], [".inp"] ]
    files['output'] += assembleStrings(fragments)
  else:
    raise ValueError("Unknown tutor step: %s" % step)
  return files


# ----------------------------------------------------------------------
def tutor(step, mode):
  files = getStepFiles(step)

  import os
  dirFiles = os.listdir(os.getcwd())

  if "check" == mode:
    print "Checking to make sure you have the following input files " \
          "for step '%s':" % step
    for filename in files['input']:
      if filename in dirFiles:
        print "  %s...found" % filename
      else:
        print "  %s...MISSING" % filename
  elif "clean" == mode:
    print "Removing any old output files from step '%s':" % step
    for filename in files['output']:
      if filename in dirFiles:
        os.remove(filename)
        print "  %s...removed" % filename
  elif "retrieve" == mode:
    import shutil
    print "Retrieving input files from archive for step '%s':" % step
    for filename in files['input']:
      if filename in dirFiles:
        print "  %s...already present" % filename
      else:
        shutil.copy("../archive/%s" % filename,
                    "./%s" % filename)
        print "  %s...retrieved from archive" % filename
  else:
    raise ValueError("Unrecognized mode: %s" % mode)
  return


# ----------------------------------------------------------------------
if __name__ == "__main__":
  from optparse import OptionParser

  usage = "usage: %prog -m MODE -s STEP\n\n" \
          "Modes:\n" \
          "  check     Check to make sure files necessary for step " \
          "exist.\n" \
          "  clean     Remove old files that will be overwritten in " \
          "this step.\n" \
          "  retrieve  Copy files (as needed) from archive needed in " \
          "this step.\n" \
          "Steps:\n" \
          "  mesh      Generating the finite-element mesh.\n" \
          "  setup     Setting up the PyLith input files.\n" \
          "  run1      Running the simulation on 1 processor.\n" \
          "  viz1      Visualizing the output from a simulation on " \
          "1 processor.\n" \
          "  run2      Running the simulation on 2 processors.\n" \
          "  viz2      Visualizing the output from a simulation on " \
          "2 processors.\n" \
          "  all       Apply mode to each step in succession."

  parser = OptionParser(usage=usage)
  parser.add_option("-m", "--mode", dest="mode",
                    type="string", metavar="MODE",
                    help="Set tutor to MODE.")
  parser.add_option("-s", "--step", dest="step",
                    type="string", metavar="STEP",
                    help="Perform mode operation for STEP.")
  (options, args) = parser.parse_args()
  if len(args) != 0:
    parser.error("Incorrent number of arguments.")

  if options.step is None or options.mode is None:
    parser.error("Both step and mode are required options.")

  allsteps = "mesh setup run1 viz1 run2 viz2"
  if "all" == options.step:
    for istep in allsteps.split():
      tutor(istep, options.mode)
  else:
    tutor(options.step, options.mode)


# End of file 
