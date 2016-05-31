#!/usr/bin/env python
#
# Python script to run initial static simulations, generate initial state
# variables spatial databases, run a restart simulation, and check the results.

import os
import sys
import subprocess
# import pdb
# pdb.set_trace()

materials = ["maxps", "genmaxps", "powerlawps", "dpps",
             "max3d", "genmax3d", "powerlaw3d", "dp3d"]

if len(sys.argv) != 1:
  raise ValueError("usage: run_all.py")

for d in ["output", "logs"]:
  if not os.path.isdir(d):
    os.mkdir(d)

# ----------------------------------------------------------------------
def runPyLith(args, logFilename):
  log = open("logs/" + logFilename, "w")
  subprocess.call("pylith " + args, stdout=log, stderr=log, shell=True)
  log.close()
  return

# ----------------------------------------------------------------------
for materialNum in range(len(materials)):
  material = materials[materialNum]
  print "Testing material %s:" % material
  meshCfg = "quad.cfg "
  if (materialNum > 3):
    meshCfg = "hex.cfg "
  staticArgs = meshCfg + "grav_static_" + material + ".cfg"
  staticLog = "grav_static_" + material + ".log"
  runPyLith(staticArgs, staticLog)
  subprocess.call("./gen_" + material + "_statedb.py", shell=True)
  restartArgs = meshCfg + "grav_restart_" + material + ".cfg"
  restartLog = "grav_restart_" + material + ".log"
  runPyLith(restartArgs, restartLog)
  subprocess.call("./check_" + material + "_restart.py", shell=True)
    
