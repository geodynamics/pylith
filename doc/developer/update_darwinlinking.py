#!/usr/bin/env python

import sys
import glob
import os
import os.path
import fnmatch
import subprocess

python_version = "%d.%d" % (sys.version_info.major, sys.version_info.minor)
pylith_version = "2.1.0"
spatialdata_version = "1.0.0"

binDir = "bin"
libDir = "lib"
modulesDir = "lib/python%s/site-packages" % python_version

def find_files(dir, pattern):
    matches = []
    for root, dirnames, filenames in os.walk(dir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches


def update_deplibs(file):
    # Ignore symbolic links and directories
    if os.path.islink(file) or os.path.isdir(file):
        return

    print("Updating %s..." % file)
    
    p = subprocess.Popen(["otool", "-L", file], stdout=subprocess.PIPE)
    output = p.communicate()[0]
    deplibs = []
    lines = output.split("\t")
    for line in lines[1:]:
        deplibs.append(line.split()[0])
    for deplibOrig in deplibs:
        # Ignore system files
        if deplibOrig.startswith("/usr") or deplibOrig.startswith("/System"):
            continue
        
        deplibName = os.path.split(deplibOrig)[1]
        deplibNew = "@executable_path/../lib/%s" % deplibName
        cmd = "install_name_tool -change %s %s %s" % (deplibOrig, deplibNew, file)
        subprocess.check_call(cmd, shell=True)

    return


# bin
for bin in glob.glob(os.path.join(binDir, "*")):
    update_deplibs(bin)

# modules
modules = find_files(modulesDir, "*.so")
for module in modules:
    update_deplibs(module)

dylibs = find_files(libDir, "*.dylib")
for lib in dylibs:
    # Can't change stub files
    if lib.startswith("lib/libgcc_ext"):
        continue

    update_deplibs(lib)
            

# Set version number in pylith and spatialdata libraries
#
# NOTE: Currently disabled because changing the version messes up
# loading libraries from modules/libraries that have already linked
# against the previous version number.

#cmd = "install_name_tool -id %s lib/libpylith.dylib" % pylith_version
#subprocess.check_call(cmd, shell=True)
        
#cmd = "install_name_tool -id %s lib/libspatialdata.dylib" % spatialdata_version
#subprocess.check_call(cmd, shell=True)
