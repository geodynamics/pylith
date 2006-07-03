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
def showUsage():
  print "usage: tutor.py MODE STEP"
  print "  Modes:"
  print "    check     Check to make sure files necessary for step exist."
  print "    clean     Remove old files that will be overwritten in this step."
  print "    retrieve  Copy files (as needed) from archive needed in " \
        "this step."
  print "  Steps:"
  print "    mesh      Generating the finite-element mesh."
  print "    setup     Setting up the PyLith input files."
  print "    run1      Running the simulation on 1 processor."
  print "    viz1      Visualizing the output from a simulation on " \
        "1 processor."
  print "    run2      Running the simulation on 2 processors."
  print "    viz2      Visualizing the output from a simulation on " \
        "2 processors."
  print "    all       Apply mode to each step in succession."
  return


# ----------------------------------------------------------------------
def assembleStrings(prefixes, suffixes):
  """Assemble list of strings from lists of prefixes and suffixes."""

  values = []
  for prefix in prefixes:
    for suffix in suffixes:
      values.append(prefix+suffix)
  return values


# ----------------------------------------------------------------------
def getStepFiles(curstep):

  files = {'input': [],
           'output': []}
  if "mesh" == curstep:
    files['input'] += assembleStrings(prefixes=[root], suffixes=[".geo"])
    files['output'] += assembleStrings(prefixes=[root], suffixes=[".netgen"])
    
  elif "setup" == curstep:
    suffixes = [".netgen", ".fault.par", ".par", ".aux"]
    files["input"] += assembleStrings(prefixes=[root], suffixes=suffixes)

    suffixes = [".coord", ".connect", ".bc", ".w01.wink", ".1.fcoord",
                ".1.fbc", ".split"]
    files['output'] += assembleStrings(prefixes=[root], suffixes)
    files['output'] += "runbm5.sh"

  elif "run1" == curstep:
    suffixes = [".coord", ".coonect", ".split", ".bc", ".fuldat", ".prop", \
                ".statevar", ".time"]
    files['input'] += assembleStrings(prefixes=[root], suffixes=suffixes)

    prefixes = [root + "_1"]
    suffixes = [".coord", ".connect", ".split", ".bc"]
    files['output'] += assembleStrings(prefixes=prefixes, suffixes=suffixes)

    prefixes = [root + "_1.[0-9]*"]
    suffixes = [".bc", ".connect", ".coord", ".split", ".ascii", ".vtk", \
                ".fuldat", ".prop", ".statevar", ".time"]
    files['output'] += assembleStrings(prefixes=[root1], suffixes=suffixes)

    prefixes = root + "_1.[0-9]*.gmesh"
    suffixes = ["", ".t10", ".time"]
    prefixesA = assembleStrings(prefixes=prefixes, suffixes=suffixes)
    prefixes = ["", ".t10", ".time"]
    suffixes = [".00000", ".00010", ".00050", ".00100"]
    suffixes = assembleStrings(prefixes=prefixes, suffixes=suffixes)
    suffixes = assembleString(prefixes=suffixes, suffixes=".inp")
    files['ouput'].append(assembleStrings(prefixes=prefixesA,
                                          suffixes=suffixes))
  else:
    error

  return files
                          

  elif [ "run1" == $curstep ]; then

    files_output="$files_output "`echo ${root}_1.[0-9]*.gmesh{,.t10,.time.{00000,00010,00050,00100}}.inp`
    files_output="$files_output "`echo ${root}_1.[0-9]*.mesh{,.t10,{.time.{00000,00010,00050,00100}}}.inp`
    files_output="$files_output "`echo ${root}_1.[0-9]*.mesh.split.time.{00000,00010,00050,00100}.inp`

  elif [ "viz1" == $curstep ]; then
    files_input=`echo ${root}_1.0.mesh{,.time.00010}.inp`
    files_output="${root}_1.0.mesh.t00010.inp"

  elif [ "run2" == $curstep ]; then
    files_input="runbm5.sh"
    files_input="$files_input "`echo $root.{coord,connect,split,bc}`
    files_input="$files_input "`echo $root.{fuldat,prop,statevar,time}`

    files_output=`echo ${root}_2.{coord,connect,split,bc}`
    files_output="$files_output "`echo ${root}_2.[0-9]*.{bc,connect,coord,split}`
    files_output="$files_output "`echo ${root}_2.[0-9]*.{ascii,vtk}`
    files_output="$files_output "`echo ${root}_2.[0-9]*.{fuldat,prop,statevar,time}`
    files_output="$files_output "`echo ${root}_2.[0-9]*.gmesh{,.t10,.time.{00000,00010,00050,00100}}.inp`
    files_output="$files_output "`echo ${root}_2.[0-9]*.mesh{,.t10,{.time.{00000,00010,00050,00100}}}.inp`
    files_output="$files_output "`echo ${root}_2.[0-9]*.mesh.split.time.{00000,00010,00050,00100}.inp`

  elif [ "viz2" == $curstep ]; then
    files_input=`echo ${root}_2.{0,1}.mesh{,.time.00010}.inp`
    files_output=`echo ${root}_2.{0,1}.mesh.t00010.inp`

  else
    echo "Unrecognized step: $step"
    show_usage
    exit 1
  fi
}


# version
__id__ = "$Id$"

# End of file 
