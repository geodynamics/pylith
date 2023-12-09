# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# Application for searching for and running PyLith .cfg files.

import sys
import argparse
import pathlib
import textwrap
import os

from pylith.utils.SimulationMetadata import fromFile
from pylith.apps.PyLithApp import PyLithApp

class RunnerApp():
    """Application for searching for and running PyLith .cfg files.
    """

    def main(self, **kwargs):
        """Main entry point.

        Keyword arguments:
            searchpath (str), default: "."
                Search path for .cfg files.
        """
        args = argparse.Namespace(**kwargs) if kwargs else self._parse_command_line()

        for filename in sorted(pathlib.Path(args.searchpath).glob("**/*.cfg")):
            metadata = fromFile(filename)
            if metadata:
                if metadata.arguments:
                    pylith_arguments = metadata.arguments
                    if args.nodes > 1:
                        pylith_arguments += [f"--nodes={args.nodes}"]
                    self._run_pylith(filename, pylith_arguments)
                elif args.verbose:
                    print(f"WARNING: File {filename} missing PyLith arguments.")
            elif args.verbose:
                print(f"INFO: File {filename} missing simulation metadata.")

    def _run_pylith(self, filename, arguments):
        """Run PyLith simulation using given arguments.

        Args:
            filename (str)
                Path to simulation parameter file.
            arguments (list of str)
                Command line arguments.
        """
        workdir = filename.parent
        cwd = os.getcwd()

        args = " ".join(arguments)
        if workdir.name:
            os.chdir(workdir)
            print(f"RUNNING {workdir} - pylith {args}...")
        else:
            print(f"RUNNING: pylith {args}...")

        app = PyLithApp()
        app.run(argv=["pylith"] + arguments)
        os.chdir(cwd)

    def _parse_command_line(self):
        """Parse command line arguments.

        Returns (argsparse.Namespace)
           Command line arguments.
        """
        DESCRIPTION = (
            "Application for searching for and running PyLith simulations."
        )

        parser = argparse.ArgumentParser(description=DESCRIPTION,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("--path", action="store",
                            dest="searchpath", default=".", help="Search path for .cfg files.")
        parser.add_argument("--verbose", action="store_true",
                            dest="verbose", help="Report missing metadata.")
        parser.add_argument("--nodes", action="store", default=1, type=int,
                            dest="nodes", help="Number of processes to use when running PyLith.")

        args = parser.parse_args()
        return args


# End of file
