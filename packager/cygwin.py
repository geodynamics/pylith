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
# Functions used in creating PyLith binary with cygwin for windows.
#
# Usage: make_package.py

import os, sys
from os.path import basename, dirname, isabs, isdir, isfile, join
from popen2 import Popen4
import shutil

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


def cygpath(*args):
    output = ospawn("cygpath", *args)
    return output[0].rstrip()


def itwindirs(l, sourceDir):
    for src in l:
        if isinstance(src, tuple):
            src, dest = src
        elif isabs(src):
            dest = src[1:]
        else:
            dest = src
        src = cygpath("-w", src)
        if src.startswith(sourceDir):
            src = src[len(sourceDir)+1:]
        dest = "{app}\\" + dest.replace("/", "\\")
        yield src, dest
    return


def itwinfiles(l, sourceDir):
    for src in l:
        if isinstance(src, tuple):
            src, dest = src
        elif isabs(src):
            dest = dirname(src)[1:]
        else:
            dest = dirname(src)
        src = cygpath("-w", src)
        if src.startswith(sourceDir):
            src = src[len(sourceDir)+1:]
        dest = "{app}\\" + dest.replace("/", "\\")
        yield src, dest
    return


def copyAll(srcList, prefix):
    for src in srcList:
        if isinstance(src, tuple):
            src, dest = src
        elif isabs(src):
            dest = dirname(src)[1:]
        else:
            dest = dirname(src)
        if not isabs(src):
            src = join(prefix, src)
        if not isdir(dest):
            os.makedirs(dest)

        if not isdir(src):
            shutil.copy(src, dest) # faster than os.system() for small files
        else:
            os.system("cp -r %s %s" % (src, dest))
    return


def generateBashrc(prefix, package, info):
    bashrc = "." + package.lower() + "rc"
    s = open(prefix + "/" + bashrc, "w")

    stuff = {
        "line": "-" * len(info["AppVerName"]),
        }
    stuff.update(info)

    s.write(
r"""

export PATH=/usr/bin:/bin:/lib:/lib/lapack:$PATH

echo %(line)s
echo %(AppVerName)s
echo %(line)s

PS1='\[\e]0;\w\a\]\n\[\e[32m\]\u@\h \[\e[33m\]\w\[\e[0m\]\n\$ '
""" % stuff
    )
    
    return bashrc


def installCigIcon(prefix):
    cigIco = "cig.ico"
    os.system("cp %s %s" % (cigIco, prefix))
    return


def minCygwin(python):
    # Minimal set of Cygwin binaries, as determined while preparing for CFEM 2006.
    l = [
        "[.exe",
        "bash.exe",
        "cat.exe",
        "cp.exe",
        "env.exe",
        "cygattr-1.dll",
        "cygiconv-2.dll",
        "cygintl-8.dll",
        "cygncurses-9.dll",
        "cygpath.exe",
        "cygcheck.exe",
        "cygreadline7.dll",
        "cyggfortran-3.dll",
        "cyggcc_s-1.dll",
        "cygstdc++-6.dll",
        "cygncursesw-10.dll",
        # This is given special treatment.
        "cygwin1.dll",
        "cygz.dll", # for HDF5 (used by Cigma)
        "echo.exe",
        "gunzip",
        "gzip.exe",
        "lib" + python + ".dll",
        "ln.exe",
        "ls.exe",
        "ldd.exe",
        "pwd.exe",
        "python",
        python + ".exe",
        "rm.exe",
        "sed.exe",
        "sh.exe",
        "tar.exe",
        "vim-nox.exe",
        "which.exe",
        ]
    l = ["/bin/" + b for b in l]

    # Additional Python stuff.
    l.extend([
        "/lib/" + python,
        "/usr/include/" + python,
        ])

    # Lapack libraries (used by numpy) now placed in own location.
    l.extend([
        "/lib/lapack",
        ])

    return l

def installMungedCygwinDLL(munged):
    # We don't require that the user has Cygwin installed.  At the
    # same time, we need to avoid conflicts if they *do* have it
    # installed.

    # Therefore, we install a munged copy of "cygwin1.dll", performing
    # the following transformations:
    #     Cygnus Solutions -> Cigwin Solutions
    #     cygwin1S4 -> cigwin1S4
    # These are the registry key name and the shared memory object
    # name, respectively.

    # N.B.: It is significant that the replacement string is the same
    # length as the original string.  Otherwise, the DLL would become
    # corrupted.

    # This gives us the freedom to create our own chrooted environment
    # using Cygwin's mounts, without messing-up the Cygwin registry on
    # systems that do have Cygwin installed.  Registry clean-up at
    # uninstall is also made easy.

    # Ugly, yes; but so much easier than building a custom Cygwin from
    # source...
    status = os.system("""sed -b """
                       """-e "s/Cygnus Solutions/Cigwin Solutions/g" """
                       """-e s/cygwin1S4/cigwin1S4/g """
                       """/bin/cygwin1.dll > %s""" % munged)
    if status != 0:
        sys.exit("sed: exit %d" % status)
    return


def installMinimalCygwin(python, prefix):
    cygwin = minCygwin(python)
    copyAll(cygwin, prefix)
    installMungedCygwinDLL("bin/cygwin1.dll")
    return


def rebaseDLLs(stage):
    # Work-around "unable to remap xyz.dll to same address as parent".
    # The base and offset were determined through trial-and-error.

    print "Rebasing DLLs."
    
    dlls = []

    for dirpath, dirnames, filenames in os.walk(stage):
        for filename in filenames:
            if filename.endswith(".dll"):
                dll = join(dirpath, filename)
                print dll
                dlls.append(dll)

    if not dlls:
        sys.exit("no DLLs found!")
    dlls = ' '.join(dlls)

    status = os.system("""rebase -v -d -b 0x70000000 -o 0x100000 %s""" % dlls)
    if status != 0:
        sys.exit("rebase: exit %d" % status)
    return


def stageInstallation(prefix, workdir, python, pl):
    
    stage = join(workdir, "buildbot-win-stage")
    shutil.rmtree(stage, ignore_errors=True)
    os.mkdir(stage)
    os.chdir(stage)

    print "Copying packing list to %s." % stage
    copyAll(pl.files(), prefix)
    
    print "Installing minimal cygwin."
    installMinimalCygwin(python, prefix)
    
    rebaseDLLs(stage)
    
    return stage


def generateISS(workdir, sourceDir, info, pl, bashrc, python, stage=None):
    s = open(workdir + "/buildbot.iss", "w")


    # [Setup]
    s.write(
"""; Inno Setup Script generated by BuildBot

[Setup]
AppName=%(AppName)s
AppVerName=%(AppVerName)s
AppPublisher=Computational Infrastructure for Geodynamics
AppPublisherURL=http://www.geodynamics.org/
AppSupportURL=http://www.geodynamics.org/
AppUpdatesURL=http://www.geodynamics.org/
DefaultDirName={pf}\\%(AppName)s
DefaultGroupName=%(AppName)s
SourceDir=%(SourceDir)s

""" % info
    )


    # [Tasks]
    s.write(
"""
[Tasks]
Name: "desktopicon"; Description: "Create a &desktop icon"; GroupDescription: "Additional icons:"

""")


    # [Files]
    #  NOTE: Don't use "Flags: ignoreversion" on any shared system files
    s.write(
"""
[Files]
"""
    )
    if stage:
        for pathname in os.listdir(stage):
            if isfile(pathname):
                s.write("""Source: "%s"; DestDir: "{app}"; Flags: ignoreversion\n""" % pathname)
            else:
                s.write("""Source: "%s\\*"; DestDir: "{app}\\%s"; Flags: ignoreversion recursesubdirs\n""" % (pathname, pathname))
    else:
        for src, dest in itwindirs(pl.directories, sourceDir):
            s.write("""Source: "%s\\*"; DestDir: "%s"; Flags: ignoreversion recursesubdirs\n""" % (src, dest))
        for src, dest in itwinfiles(pl.files(), sourceDir):
            s.write("""Source: "%s"; DestDir: "%s"; Flags: ignoreversion\n""" % (src, dest))
        for f in [bashrc, "cig.ico"]:
            s.write("""Source: "%s"; DestDir: "{app}"; Flags: ignoreversion\n""" % f)


    # [Dirs]
    s.write(
"""
[Dirs]
; Cygwin demands a /tmp directory
Name: "{app}\\tmp"
; for our /usr mounts
Name: "{app}\\usr"
Name: "{app}\\usr\\bin"
Name: "{app}\\usr\\lib"

"""
    )


    # [Registry]
    s.write(
"""
[Registry]
Root: HKLM; Subkey: "Software\\Cigwin Solutions"; Flags: uninsdeletekey noerror
Root: HKLM; Subkey: "Software\\Cigwin Solutions\\Cygwin"; Flags: uninsdeletekey noerror
Root: HKLM; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2"; Flags: uninsdeletekey noerror
Root: HKLM; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/"; Flags: uninsdeletekey noerror
Root: HKLM; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/"; Flags: noerror; ValueType: string; ValueName: "native"; ValueData: "{app}"
Root: HKLM; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/"; Flags: noerror; ValueType: dword; ValueName: "flags"; ValueData: "00000002"
Root: HKLM; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/usr/bin"; Flags: uninsdeletekey noerror
Root: HKLM; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/usr/bin"; Flags: noerror; ValueType: string; ValueName: "native"; ValueData: "{app}\\bin"
Root: HKLM; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/usr/bin"; Flags: noerror; ValueType: dword; ValueName: "flags"; ValueData: "00000002"
Root: HKLM; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/usr/lib"; Flags: uninsdeletekey noerror
Root: HKLM; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/usr/lib"; Flags: noerror; ValueType: string; ValueName: "native"; ValueData: "{app}\\lib"
Root: HKLM; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/usr/lib"; Flags: noerror; ValueType: dword; ValueName: "flags"; ValueData: "00000002"
; for non-admin install
Root: HKCU; Subkey: "Software\\Cigwin Solutions"; Flags: uninsdeletekey
Root: HKCU; Subkey: "Software\\Cigwin Solutions\\Cygwin"; Flags: uninsdeletekey
Root: HKCU; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2"; Flags: uninsdeletekey
Root: HKCU; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/"; Flags: uninsdeletekey
Root: HKCU; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/"; ValueType: string; ValueName: "native"; ValueData: "{app}"
Root: HKCU; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/"; ValueType: dword; ValueName: "flags"; ValueData: "00000002"
Root: HKCU; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/usr/bin"; Flags: uninsdeletekey
Root: HKCU; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/usr/bin"; ValueType: string; ValueName: "native"; ValueData: "{app}\\bin"
Root: HKCU; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/usr/bin"; ValueType: dword; ValueName: "flags"; ValueData: "00000002"
Root: HKCU; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/usr/lib"; Flags: uninsdeletekey
Root: HKCU; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/usr/lib"; ValueType: string; ValueName: "native"; ValueData: "{app}\\lib"
Root: HKCU; Subkey: "Software\\Cigwin Solutions\\Cygwin\\mounts v2\\/usr/lib"; ValueType: dword; ValueName: "flags"; ValueData: "00000002"

"""
    )


    # [INI]
    s.write("""[INI]\n""")
    for name, filename, url in pl.urls:
        s.write("""Filename: "{app}\\%(filename)s.url"; Section: "InternetShortcut"; Key: "URL"; String: "%(url)s"\n""" %
                {'filename': filename, 'url': url})


    # [Icons]
    s.write(
"""
[Icons]
Name: "{group}\\%(AppName)s"; Filename: "{app}\\bin\\bash.exe"; Parameters: "--init-file %(bashrc)s -i"; WorkingDir: "{app}"; IconFilename: "{app}\\cig.ico"
Name: "{userdesktop}\\%(AppName)s"; Filename: "{app}\\bin\\bash.exe"; Tasks: desktopicon; Parameters: "--init-file %(bashrc)s -i"; WorkingDir: "{app}"; IconFilename: "{app}\\cig.ico"
""" % {'AppName': info['AppName'], 'bashrc': bashrc}
    )
    for name, filename, url in pl.urls:
        s.write("""Name: "{group}\\%s"; Filename: "{app}\\%s.url"\n""" % (name, filename))
    s.write("""Name: "{group}\Uninstall %(AppName)s"; Filename: "{uninstallexe}"\n""" % info)


    # [Run]
    s.write(
"""
[Run]
Filename: "{app}\\bin\\bash.exe"; Description: "Launch %(AppName)s"; Parameters: "--init-file %(bashrc)s -i"; WorkingDir: "{app}"; Flags: nowait postinstall skipifsilent

""" % {'AppName': info['AppName'], 'bashrc': bashrc}
    )


    # [UninstallDelete]
    s.write("""[UninstallDelete]\n""")
    for name, filename, url in pl.urls:
        s.write("""Type: files; Name: "{app}\%s.url"\n""" % filename)
    for filename in [".bash_history"]:
        s.write("""Type: files; Name: "{app}\%s"\n""" % filename)

    s.write("\n; end of file\n")

    return


def createInstaller(python, stage, workdir, name, package, version, pl,
                           installer):

    sourceDir = cygpath("-w", stage)
    
    info = {
        "AppName": name,
        "AppVerName": name + " v" + version,
        "SourceDir": sourceDir,
        }

    installCigIcon(stage)
    bashrc = generateBashrc(stage, package, info)

    generateISS(workdir, sourceDir, info, pl, bashrc, python, stage)

    os.chdir(workdir)

    os.system("unix2dos buildbot.iss")
    
    status = os.system("iscc buildbot.iss")
    if status != 0:
        sys.exit("iscc: exit %d" % status)

    status = os.system("mv %s/Output/setup.exe %s.exe" % (stage, installer))
    if status != 0:
        sys.exit("mv: exit %d" % status)

    return


