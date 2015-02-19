#!/usr/bin/env python

# Run all of the examples. Execute this script from the examples
# directory.


import os
import glob
import sys

targets = ["all"]
if len(sys.argv) > 1:
    targets = sys.argv

examplesDir = os.getcwd()

# ======================================================================
def clean(dir):
    for pattern in ["*.vtk",
                    "output/*.vtk",
                    "output/*.h5",
                    "output/*.xmf",
                    ]:
        trash = glob.glob(pattern)
        for f in trash:
            os.remove(f)
    return


def run_pylith(dir, examples, nprocs=1):
    """
    Run pylith on examples in a given directory.

    @param dir Name of directory with input files.
    @param examples List or tuple of .cfg files for examples.
    @param nprocs Number of processes to use.
    """
    wkd = os.path.join(examplesDir, dir)
    os.chdir(wkd)
    clean(dir)
    for simfiles in examples:
        cmd = "pylith --nodes=%d %s" % (nprocs, simfiles)
        print "WORKDIR: %s, RUNNING %s" % (wkd, cmd)
        sys.stdout.flush()
        os.system(cmd)

# ======================================================================
if "twocells" in targets or "all" in targets:
    dir = "twocells/twotri3"
    examples = ("axialdisp.cfg", 
                "sheardisp.cfg", 
                "dislocation.cfg",
                )
    run_pylith(dir, examples)


    dir = "twocells/twoquad4"
    examples = ("axialdisp.cfg",
                "sheardisp.cfg",
                "axialtract.cfg",
                "dislocation.cfg",
                )
    run_pylith(dir, examples)


    dir = "twocells/twotet4"
    examples = ("axialdisp.cfg",
                "dislocation.cfg",
                )
    run_pylith(dir, examples)


    dir = "twocells/twohex8"
    examples = ("axialdisp.cfg",
                "sheardisp.cfg",
                "dislocation.cfg",
                )
    run_pylith(dir, examples)


    dir = "twocells/twotet4-geoproj"
    examples = ("dislocation.cfg",
                )
    run_pylith(dir, examples)
    

# ----------------------------------------------------------------------
if "3d/tet4" in targets or "all" in targets:
    dir = "3d/tet4"
    examples = ("step01.cfg",
                "step02.cfg",
                "step03.cfg",
                "step04.cfg",
                )
    run_pylith(dir, examples, nprocs=1)
    run_pylith(dir, examples, nprocs=2)
    run_pylith(dir, examples, nprocs=3)
    run_pylith(dir, examples, nprocs=4)
    run_pylith(dir, examples, nprocs=6)


# ----------------------------------------------------------------------
if "3d/hex8" in targets or "all" in targets:
    dir = "3d/hex8"
    examples = ("step01.cfg",
                "step02.cfg",
                "step03.cfg",
                "step04.cfg",
                "step05.cfg",
                "step06.cfg",
                "step07.cfg",
                "step08.cfg",
                "step09.cfg",
                "step10.cfg",
                "step11.cfg",
                "step12.cfg",
                "step13.cfg",
                "step14.cfg",
                "step15.cfg",
                "step16.cfg",
                "step17.cfg",
                "step18.cfg",
                "step19.cfg",
                "step20.cfg",
                "step21.cfg --problem=pylith.problems.GreensFns",
                )
    run_pylith(dir, examples, nprocs=1)

    examples = ("step01.cfg",
                "step03.cfg",
                "step06.cfg",
                "step15.cfg",
                "step19.cfg",
                "step20.cfg",
                "step21.cfg --problem=pylith.problems.GreensFns",
                )
    run_pylith(dir, examples, nprocs=2)
    run_pylith(dir, examples, nprocs=4)


# ----------------------------------------------------------------------
if "2d/subduction" in targets or "all" in targets:
    dir = "2d/subduction"
    examples = ("step01.cfg",
                "step02.cfg",
                "step03.cfg",
                "step04.cfg",
                )
    run_pylith(dir, examples, nprocs=1)
    run_pylith(dir, examples, nprocs=2)
    run_pylith(dir, examples, nprocs=4)
    run_pylith(dir, examples, nprocs=6)


# ----------------------------------------------------------------------
if "bar_shearwave" in targets or "all" in targets:
    dir = "bar_shearwave/tri3"
    examples = ("pylithapp.cfg",
                )
    run_pylith(dir, examples, nprocs=1)
    run_pylith(dir, examples, nprocs=3)


    dir = "bar_shearwave/quad4"
    examples = ("prescribedrup.cfg",
                "spontaneousrup.cfg spontaneousrup_staticfriction.cfg",
                "spontaneousrup.cfg spontaneousrup_slipweakening.cfg",
                "spontaneousrup.cfg spontaneousrup_ratestateageing.cfg",            
                )
    run_pylith(dir, examples, nprocs=1)
    run_pylith(dir, examples, nprocs=3)


    dir = "bar_shearwave/tet4"
    examples = ("pylithapp.cfg",
                )
    run_pylith(dir, examples, nprocs=1)
    run_pylith(dir, examples, nprocs=3)
    

    dir = "bar_shearwave/hex8"
    examples = ("pylithapp.cfg",
                )
    run_pylith(dir, examples, nprocs=1)
    run_pylith(dir, examples, nprocs=3)


# ----------------------------------------------------------------------
if "2d/greensfns" in targets or "all" in targets:
    dir = "2d/greensfns/strikeslip"
    examples = ("eqsim.cfg",
                "--problem=pylith.problems.GreensFns",
                )
    run_pylith(dir, examples, nprocs=1)
    run_pylith(dir, examples, nprocs=2)
    run_pylith(dir, examples, nprocs=4)
    run_pylith(dir, examples, nprocs=6)

    dir = "2d/greensfns/reverse"
    examples = ("eqsim.cfg",
                "--problem=pylith.problems.GreensFns",
                )
    run_pylith(dir, examples, nprocs=1)
    run_pylith(dir, examples, nprocs=2)
    run_pylith(dir, examples, nprocs=4)
    run_pylith(dir, examples, nprocs=6)
