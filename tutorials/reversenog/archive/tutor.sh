#!/bin/bash
# ======================================================================
#
# Tutor shell script to help user run tutorial using SCEC benchmark 5.
#
# ======================================================================

archive="../archive"
root="bm5"

# ----------------------------------------------------------------------
function show_usage () {
  echo "usage: tutor.sh MODE STEP"
  echo "  Modes:"
  echo "    check     Check to make sure files necessary for step exist."
  echo "    clean     Remove old files that will be overwritten in this step."
  echo "    retrieve  Copy files (as needed) from archive needed in this step."
  echo "  Steps:"
  echo "    mesh      Generating the finite-element mesh."
  echo "    setup     Setting up the PyLith input files."
  echo "    run1      Running the simulation on 1 processor."
  echo "    viz1      Visualizing the output from a simulation on 1 processor."
  echo "    run2      Running the simulation on 2 processors."
  echo "    viz2      Visualizing the output from a simulation on 2 processors."
  echo "    all       Apply mode to each step in succession."
}

# ----------------------------------------------------------------------
function get_step_lists () {
  if [ "mesh" == $curstep ]; then
    files_input=$root.geo
    files_output=$root.netgen

  elif [ "setup" == $curstep ]; then
    files_input=`echo $root.{netgen,fault.par,par,aux}`
    files_output=`echo $root.{coord,connect,bc,w01.wink,1.fcoord,1.fbc,split}`

  elif [ "run1" == $curstep ]; then
    files_input="runbm5.sh"
    files_input="$files_input "`echo $root.{coord,connect,split,bc}`
    files_input="$files_input "`echo $root.{fuldat,prop,statevar,time}`

    files_output=`echo ${root}_1.{coord,connect,split,bc}`
    files_output="$files_output "`echo ${root}_1.[0-9]*.{bc,connect,coord,split}`
    files_output="$files_output "`echo ${root}_1.[0-9]*.{ascii,vtk}`
    files_output="$files_output "`echo ${root}_1.[0-9]*.{fuldat,prop,statevar,time}`
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

# ----------------------------------------------------------------------
function tutor () {
  get_step_lists

  if [ "check" == $mode ]; then
    echo "Checking to make sure you have the following files:"
    for file in $files_input; do
      echo -n "  $file..."
      if [ -r $file ]; then
        echo "found"
      else
        echo "MISSING"
      fi
    done
  elif [ "clean" == $mode ]; then
    echo "Removing any old files:"
    for file in $files_output; do
      if [ -r $file ] || [ -L $file ]; then
        echo "  $file...removed"
        rm $file
      fi
    done
  elif [ "retrieve" == $mode ]; then
    echo "Retrieving files from archive:"
    for file in $files_input; do
      echo -n "  $file..."
      if [ -r $file ]; then
        echo "already present"
      else
        echo "retrieving from archive"
        cp $archive/$file .
      fi
    done
  else
    echo "Unrecognized mode: $mode"
    show_usage
    exit 1
  fi
}

# ----------------------------------------------------------------------
# Check number of command line arguments
if [ $# != 2 ]; then
  show_usage
  exit 1
fi

# Get command line arguments
mode=$1
step=$2

allsteps="mesh setup run1 viz1 run2 viz2"

if [ "all" == $step ]; then
  for curstep in $allsteps; do
    echo "Running tutor in mode $mode for step $curstep"
    tutor
  done
else
  curstep=$step
  tutor
fi

# End of file
