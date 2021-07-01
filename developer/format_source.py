#!/usr/bin/env python3
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

"""Format C++ and Python to match coding style conventions used in the PyLith projectself.

Uses uncrustify to format C++ code and autopep8 to format Python code.
"""

import os
import subprocess
import glob
import fnmatch


def run_cmd(cmd, show_progress=False):
    """Run command as a subprocess.

    :param cmd: Command to run [program, options].
    """
    if show_progress:
        print(" ".join(cmd))

    subprocess.check_call(cmd)
    return


def verify_cmd(cmd):
    """Verify command runs without errors."""

    try:
        subprocess.check_call(cmd)
    except:
        raise IOError("Error running '{}'.".format(" ".join(cmd)))
    return


class Formatter(object):
    """Abstract base class for formatting objects."""

    def __init__(self, config_dir, show_progress=True):
        """Constructor.

        :config_dir: Directory containing uncrustify configuration file (uncrustify.cfg).

        :show_progress: Show progress if True, be silent if False.
        """
        self.config_dir = config_dir
        self.show_progress = show_progress

        self.config_filename = None
        return

    def _config_filename(self):
        raise NotImplementedError("Implement _config_filename() in subclass.")

    def _prog_test(self):
        raise NotImplementedError("Implement _prog_test() in subclass.")

    def _prog_format(self, filename):
        raise NotImplementedError("Implement _prog_format() in subclass.")

    def initialize(self):
        """Verify formatting program is installed and project configuration file exists."""
        self._verify_configuration()
        self._verify_installation()
        return

    def format_file(self, filename):
        """Format file.

        :param filename: Name of file to format.
        """
        run_cmd(self._prog_format(filename), show_progress=self.show_progress)
        return

    def _verify_configuration(self):
        """Verify project configuration file exists."""
        path_filename = os.path.join(self.config_dir, self._config_filename())
        if not os.path.isfile(path_filename):
            raise IOError(
                "Could not find configuration file '{}'.".format(path_filename))
        self.config_filename = path_filename
        return

    def _verify_installation(self):
        """Verify formatting program is installed."""
        verify_cmd(self._prog_test())
        return


class FormatterCPlusCplus(Formatter):
    """Formatter for C/C++ source files."""

    LANGUAGE = {
        ".h": "C",
        ".c": "C",
        ".hh": "CPP",
        ".cc": "CPP",
    }

    def _config_filename(self):
        return "uncrustify.cfg"

    def _prog_test(self):
        return ["uncrustify", "--version"]

    def _prog_format(self, filename):
        suffix = "." + filename.split(".")[-1]
        language = self.LANGUAGE[suffix]
        options = ["-l", language, "-c",
                   self.config_filename, "--no-backup", filename]
        return ["uncrustify"] + options


class FormatterPython(Formatter):
    """Formatter for Python source files."""

    def _config_filename(self):
        return "autopep8.cfg"

    def _prog_test(self):
        return ["autopep8", "--version"]

    def _prog_format(self, filename):
        options = [
            "--global-config={}".format(self.config_filename), "--in-place", filename]
        return ["autopep8"] + options


class App(object):
    """Application to reformat C++ and Python source code.
    """

    def __init__(self):
        """Constructor."""
        self.show_progress = False
        self.config_dir = None
        return

    def main(self):
        """Application driver."""
        args = self._parse_command_line()
        self.config_dir = args.config_dir
        self.show_progress = args.show_progress

        if args.cplusplus:
            self._reformat_cplusplus(args.cplusplus)

        if args.python:
            self._reformat_python(args.python)
        return

    def _reformat_cplusplus(self, flags):
        """Reformat C++ files.

        :param flags: Flags indicating what files to reformat ["all", FILE, PATTERN]
        """

        formatter = FormatterCPlusCplus(self.config_dir, self.show_progress)
        formatter.initialize()

        if flags == "all":
            targets = self._find_files(["*.cc", "*.hh", "*.c", "*.h"])
        else:
            targets = glob.glob(flags)

        for target in targets:
            formatter.format_file(target)
        return

    def _reformat_python(self, flags):
        """Reformat Python files.

        :param flags: Flags indicating what files to reformat ["all", FILE, PATTERN]
        """

        formatter = FormatterPython(self.config_dir, self.show_progress)
        formatter.initialize()

        if flags == "all":
            targets = self._find_files([".py"])
        else:
            targets = glob.glob(flags)

        for target in targets:
            formatter.format_file(target)
        return

    def _find_files(self, patterns):
        """Find files matching pattern.

        :param patterns: Patterns to search for.

        :returns: Files matching pattern.
        """
        files = []
        root, _, names = os.walk(os.getcwd())
        for pattern in patterns:
            files += [os.path.join(root, filename)
                      for filename in fnmatch.filter(names, pattern)]
        return files

    def _parse_command_line(self):
        """Parse command line arguments.
        """
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--cplusplus", action="store", dest="cplusplus",
                            default=None, help="[None, 'all', FILE, PATTERN]")
        parser.add_argument("--python", action="store", dest="python",
                            default=None, help="[None, 'all', FILE, PATTERN]")
        parser.add_argument("--config-dir", action="store", dest="config_dir",
                            default="doc/developer", help="Directory containing config files for formatters.")
        parser.add_argument("--quiet", action="store_false",
                            dest="show_progress", default=True)
        parser.add_argument("--debug", action="store_true",
                            dest="debug", default=True)
        return parser.parse_args()


if __name__ == "__main__":
    App().main()
