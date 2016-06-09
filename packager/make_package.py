#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# Create PyLith binary package.
#
# Run from the top-level PyLith build dir.
#
# Usage: make_package.py

import os, sys, platform
from popen2 import Popen4
from glob import glob
from os.path import basename, dirname, isabs, isdir, isfile, join
import shutil
from distutils.sysconfig import parse_makefile
from ConfigParser import ConfigParser
import cygwin


class Packaging(object):
    _name = 'packaging'
    _attrs = ['bin_dirs', 
              'lib_dirs', 
              'misc_dirs',
              'files',
              'strip_list',
              'exclude',
              'scripts',
              'urls',
              ]


class Config(object):

    sections = [Packaging]

    #@classmethod
    def readPackingList(cls):

        parser = ConfigParser()
        curdir = dirname(__file__)
        parser.read(join(curdir, "packinglist.cfg"))

        config = Config()

        for sectionCls in cls.sections:
            sectionDct = parser._sections.get(sectionCls._name, None)
            section = None
            if sectionDct is not None:
                section = sectionCls()
                for attr in section._attrs:
                    value = sectionDct.get(attr, "")
                    value = value.split()
                    setattr(section, attr, value)
            setattr(config, sectionCls._name, section)
        
        return config
    readPackingList = classmethod(readPackingList)


def getMakeInfo():
    # NYI: This is Autotools-specific.
    makefile = parse_makefile("Makefile")
    info = {'package_name': makefile['PACKAGE_NAME'],
            'package': makefile['PACKAGE'],
            'version': makefile['VERSION'],
            'srcdir': makefile['top_srcdir'],
            'prefix': makefile['prefix'],
            'python_version': makefile['PYTHON_VERSION'],
            }
    return info



def walkDirTree(dirname):
    """
    Get list of files in tree under dirname.
    """
    fullNames = []
    for dirpath, dirs, files in os.walk(dirname):
        for f in files:
            fullNames.append(os.path.join(dirpath, f))
    return fullNames


def filterList(l, excludes):
    try:
        import re, functools, operator

        pats = [re.compile(ex) for ex in excludes]
        return [i for i in l if not functools.reduce(operator.or_, [bool(pat.match(i)) for pat in pats])]
    except (ImportError, AttributeError):
        import re, operator

        pats = [re.compile(ex) for ex in excludes]
        return [i for i in l if not reduce(operator.or_, [bool(pat.match(i)) for pat in pats])]


class PackingList(object):
    def __init__(self, config, opSys, python):
        
        binDirs = config.packaging.bin_dirs
        libDirs = config.packaging.lib_dirs
        miscDirs = config.packaging.misc_dirs
        if platform.machine() == "x86_64":
            libDirs.append("lib64")
        self.directories = binDirs + libDirs + miscDirs

        self.extraFiles = config.packaging.files
        self.exclude = config.packaging.exclude

        self.programs = []
        for d in binDirs:
            self.programs.extend(walkDirTree(d))
        self.programs = filterList(self.programs, self.exclude)

        self.libraries = []
        for d in libDirs:
            self.libraries.extend(walkDirTree(d))
        self.libraries = filterList(self.libraries, self.exclude)
            
        self.misc = []
        for d in miscDirs:
            self.misc.extend(walkDirTree(d))
        self.misc = filterList(self.misc, self.exclude)

        # KLUDGE Explicitly add include/pythonX.X/pyconfig.h if it exists.
        pythonConfigFile = "include/" + python + "/pyconfig.h"
        if isfile(pythonConfigFile):
            self.misc.append(pythonConfigFile)
            
        # Scripts are application specific.
        self.scripts = []
        for s in config.packaging.scripts:
            s = "bin/" + s
            self.scripts.append(s)

        # Suffix for libraries that will be stripped of symbols.
        if opSys == "linux":
            libSuffix = ".so"
        elif opSys == "darwin":
            libSuffix = ".dylib"
        elif opSys == "win":
            libSuffix = ".dll"
        else:
            sys.exit("Unknown OS: " + opSys)

        self.stripList = []
        if opSys == "win":
            for f in config.packaging.strip_list:
                fdll = f.replace("lib/lib", "bin/cyg")+"-0"+libSuffix
                if fdll in self.libraries:
                    self.stripList.append(fdll)
        else:
            for f in config.packaging.strip_list:
                flib = f+libSuffix
                if flib in self.libraries:
                    self.stripList.append(flib)

        cig = [("CIG", "cig", "http://www.geodynamics.org/")]
        self.urls = cig + tupleUp(config.packaging.urls, 3)

        return


    def addFile(self, f):
        self.misc.append(f)

    def addDirectory(self, d):
        self.directories.append(d)
        newFiles = walkDirTree(d)
        newFiles = filterList(newFiles, self.exclude)
        self.misc.extend(newFiles)


    def files(self):
        for f in self.programs:
            yield f
        for f in self.libraries:
            yield f
        for f in self.misc:
            yield f
        for f in self.extraFiles:
            yield f
        return


    def all(self):
        for d in self.directories:
            yield d
        for f in self.files():
            yield f
        return


def tupleUp(l, n):
    tuples = []
    i = iter(l)
    try:
        while True:
            t = []
            for count in xrange(n):
                t.append(i.next())
            tuples.append(tuple(t))
    except StopIteration:
        pass
    return tuples


def spawn(*argv):
    print ' '.join(argv)
    status = os.spawnvp(os.P_WAIT, argv[0], argv)
    if status != 0:
        statusMsg = "%s: %s: exit %d" % (sys.argv[0], argv[0], status)
        sys.exit(statusMsg)
    return


def ospawn(*argv):
    print ' '.join(argv)
    child = Popen4(argv)

    child.tochild.close()

    output = child.fromchild.readlines()
    status = child.wait()

    exitStatus = None
    if (os.WIFSIGNALED(status)):
        statusStr = "signal %d" % os.WTERMSIG(status)
    elif (os.WIFEXITED(status)):
        exitStatus = os.WEXITSTATUS(status)
        statusStr = "exit %d" % exitStatus
    else:
        statusStr = "status %d" % status
    if exitStatus != 0:
        sys.exit("%s: %s: %s" % (sys.argv[0], argv[0], statusStr))

    return output


def getGitInfo(srcdir):
    workdir = os.getcwd()
    os.chdir(srcdir)
    output = ospawn("git", "branch", "-v")
    os.chdir(workdir)
    revision = "unknown"
    for line in output:
        if line[0] == "*":
            values = line.split()
            revision = values[2]
    return revision


def stripBinaries(pl, opSys):
    strip = "strip"
    if opSys == "darwin":
        # Just plain "strip" renders our Python interpreter unusable
        # by extension modules.
        strip = "strip -S"
    if len(pl.stripList) > 0:
        status = os.system(strip + " " + " ".join(pl.stripList))
        if status != 0:
            sys.exit("strip: exit %d" % status)
    return


def rewriteScripts(pl, prefix, opSys):
    # Tweak the shebang line of scripts so they are relative instead
    # of absolute.
    
    if opSys == "win":
        # chrooted environment
        relative = "#!/bin/%s"
    else:
        relative = "#!/usr/bin/env %s"
    absolute = "#!" + prefix + "/bin/"
    
    for script in pl.scripts:
        s = open(script, "r")
        lines = s.readlines()
        s.close()
        shebang = lines[0]
        if shebang.startswith(absolute):
            interpreter = shebang[len(absolute):]
            shebang = relative % interpreter
            s = open(script, "w")
            s.write(shebang)
            for line in lines[1:]:
                s.write(line)
            s.close()
            
    return


def installSource(pl, taggedArchive):
    src = "src"
    shutil.rmtree(src, ignore_errors=True)
    if not isdir(src):
        os.makedirs(src)
    spawn("tar", "-C", src, "-xf", taggedArchive)
    pl.addDirectory(src)
    return


def mkpkg():
    # WithProperties simply doesn't work in this case.
    if False:
        revision = sys.argv[2]
        if revision:
            revision = "r" + revision
        else:
            revision = "forced"
    
    makeinfo = getMakeInfo()
    python = "python"+makeinfo['python_version']
    distdir = makeinfo['package'] + "-" + makeinfo['version']
    prefix = makeinfo['prefix']
    package = makeinfo['package']
    print "version:", makeinfo['version']
    print "python:", python
    print "prefix:", prefix
    print "distdir:", distdir
    
    revision = getGitInfo(makeinfo['srcdir'])
    opSys = platform.system().lower()
    if opSys.startswith("cygwin"):
        opSys = "win"
    if opSys=="darwin":
        arch = opSys + "-" + platform.mac_ver()[0]
    else:
        arch = opSys + "-" + (platform.processor() or platform.machine()) # Why does every last detail have to be painful?

    config = Config.readPackingList()
    if config.packaging is None:
        sys.exit("ERROR: No packaging configuation.")

    # Make source distribution (also extracted into binary package)
    archive = distdir + ".tar.gz"
    taggedArchive = revision + "-" + archive
    status = os.system("make dist")
    if status != 0:
        sys.exit("make: exit %d" % status)
    os.system("mv " + archive + " " + taggedArchive)

    workdir = os.getcwd()
    os.chdir(prefix)

    pl = PackingList(config, opSys, python)

    stripBinaries(pl, opSys)
    rewriteScripts(pl, prefix, opSys)
    installSource(pl, join(workdir, taggedArchive))

    distdir_arch = distdir + "-" + arch
    if opSys == "darwin" or opSys == "linux":
        # No .dmg on Mac for now.
        os.chdir(workdir)
        os.system("rm " + distdir_arch)
        os.system("ln -s " + prefix + " " + distdir_arch)

        # Write packing list to file
        plfile = open("tar_list", "w")
        plfile.write('\n'.join([os.path.join(distdir_arch, m) for m in pl.files()]))
        plfile.close()

        archive = distdir_arch + ".tar.gz"
        taggedArchive = revision + "-" + archive
        status = os.system("tar cvzf " + archive + " -T tar_list")
        if status != 0:
            sys.exit("tar: exit %d" % status)

        os.system("mv " + archive + " " + taggedArchive)
        os.system("rm tar_list " + distdir_arch)
        shutil.rmtree("src", ignore_errors=True)

    elif opSys == "win":
        name = "pylith"
        version = makeinfo['version']
        stage = cygwin.stageInstallation(prefix, workdir, python, pl)
        cygwin.createInstaller(python, stage, workdir, name, package, version, pl,
                               revision + "-" + distdir_arch)
    else:
        sys.exit("unknown OS: " + opSys)

    return


if __name__ == "__main__":
    mkpkg()
