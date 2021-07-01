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
#
# Application for searching PyLith .cfg files.

import sys
import argparse
import pathlib
import textwrap
import os

from pylith.utils.converters import string_to_list
from pylith.utils.SimulationMetadata import fromFile

class ConfigSearchApp():
    """Application for searching PyLith .cfg files.
    """

    def __init__(self):
        """Constructor.
        """
        self.filters = {}

    def main(self, **kwargs):
        """Main entry point.

        Keyword arguments:
            searchpath (str), default: "."
                Search path for .cfg files.
            display (str), default: "all"
                List of metadata to display in search results.
            keywords (str), default: None
                Comma delimited list of keywords for filtering search results.
            features (str), default: None
                Comma delimited list of features for filtering search results.
            authors (str), default: None
                Comma delimited list of authors for filtering search results.
            version (str), default: None
                PyLith version for filtering search results.
        """
        args = argparse.Namespace(
            **kwargs) if kwargs else self._parse_command_line()
        self._set_filters(args)

        for filename in sorted(pathlib.Path(args.searchpath).glob("**/*.cfg")):
            metadata = fromFile(filename)
            if metadata:
                if not len(metadata.arguments):
                    if args.verbose:
                        print(f"INFO: Skipping file {filename} with only base metadata.")
                    continue
                filter_fn = self._apply_filters_incompatible if args.incompatible else self._apply_filters
                if filter_fn(metadata):
                    self._display_metadata(filename, metadata, args.display)
                elif args.verbose:
                    print(f"MISMATCH: File {filename} did not pass metadata filter.")
            elif args.verbose:
                print(f"INFO: File {filename} missing simulation metadata.")

    def _set_filters(self, options):
        """Set filters for display from command line option.

        Args:
            options (argsparse.Namespace)
                Command line options.
        """
        if options.keywords:
            self.filters["keywords"] = string_to_list(options.keywords)
        if options.features:
            self.filters["features"] = string_to_list(options.features)
        if options.authors:
            self.filters["authors"] = string_to_list(options.authors)
        if options.version:
            self.filters["version"] = options.version

    def _apply_filters(self, metadata):
        """Apply filters to metadata.

        Args:
            metadata (pylith.utils.SimulationMetadata)
                Simulation metadata.
        
        Returns: (bool)
            True if metadata meets filter requirements, False otherwise.
        """
        if "keywords" in self.filters:
            if not metadata.keywords:
                return False
            if not all(keyword in metadata.keywords for keyword in self.filters["keywords"]):
                return False
        if "features" in self.filters:
            if not metadata.features:
                return False
            if not all(feature in metadata.features for feature in self.filters["features"]):
                return False
        if "authors" in self.filters:
            if not metadata.authors:
                return False
            if not all(author in metadata.authors for author in self.filters["authors"]):
                return False
        if "version" in self.filters:
            if not metadata.pylith_version:
                return False
            for verMeta in metadata.pylith_version:
                if not eval("{ver} {verMeta}".format(ver=self.filters["version"], verMeta=verMeta)):
                    return False
        return True

    def _apply_filters_incompatible(self, metadata):
        """Apply filters to metadata to find incompatible parameter files.

        Args:
            metadata (pylith.utils.SimulationMetadata)
                Simulation metadata.
        
        Returns: (bool)
            True if metadata is incompatible with filter requirements, False otherwise.
        """
        if "keywords" in self.filters:
            if not metadata.keywords:
                return True
        if "features" in self.filters:
            if not "features" in metadata:
                return True
        if "authors" in self.filters:
            if not "authors" in metadata:
                return True
        if "version" in self.filters:
            if not metadata.pylith_version:
                return True
            for verMeta in metadata.pylith_version:
                if not eval("{ver} {verMeta}".format(ver=self.filters["version"], verMeta=verMeta)):
                    return True
        return False

    def _display_metadata(self, filename, metadata, options):
        """Print metadata to stdout.

        Args:
            filename (str)
                Name of simulation .cfg file.
            metadata (pylith.utils.SimulationMetadata)
                Simulation metadata.
            options (list of str)
                List of metadata to display.
        """
        INDENT = " "*4

        show_all = "all" in options
        options = string_to_list(options)
        line0 = f"{filename}"
        if "version" in options or show_all:
            if metadata.version:
                line0 += f" v{metadata.version}"
            else:
                line0 += " missing 'version'"
        if "pylith_version" in options or show_all:
            if metadata.pylith_version:
                line0 += "; requires PyLith " + " and ".join(metadata.pylith_version)
            else:
                line0 += "; missing 'pylith_version'"

        lines = []
        if "description" in options or show_all:
            if metadata.description:
                lines += [metadata.description]
            else:
                lines += ["missing 'description'"]
        if "authors" in options or show_all:
            if metadata.authors:
                lines += ["Authors: " + ", ".join(metadata.authors)]
            else:
                lines += ["missing 'authors'"]
        if "keywords" in options or show_all:
            if metadata.keywords:
                lines += ["Keywords: " + ", ".join(metadata.keywords)]
            else:
                lines += ["missing 'keywords'"]
        if "features" in options or show_all:
            if metadata.features:
                features = textwrap.fill(", ".join(metadata.features), width=120)
                lines += ["Features:"] + textwrap.indent(features, INDENT).split("\n")
            else:
                lines += ["missing 'features'"]
        if "arguments" in options or show_all:
            if metadata.arguments:
                lines += ["pylith " + " ".join(metadata.arguments)]
            else:
                lines += ["missing 'arguments'"]
        print(line0)
        if len(lines):
            print(textwrap.indent("\n".join(lines), INDENT))

    def _parse_command_line(self):
        """Parse command line arguments.

        Returns (argsparse.Namespace)
           Command line arguments.
        """
        DESCRIPTION = (
            "Application for searching PyLith .cfg parameter files."
        )

        parser = argparse.ArgumentParser(description=DESCRIPTION,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("--path", action="store",
                            dest="searchpath", default=".", help="Search path for .cfg files.")
        parser.add_argument("--display", action="store",
                            dest="display", default="all", help="List of metadata to display in search results.")
        parser.add_argument("--verbose", action="store_true", dest="verbose",
                            help="Report missing metadata.")

        parser.add_argument("--keywords", action="store", dest="keywords",
                            help="Comma delimited list of keywords for filtering search results.")
        parser.add_argument("--features", action="store", dest="features",
                            help="Comma delimited list of features for filtering search results.")
        parser.add_argument("--authors", action="store", dest="authors",
                            help="Comma delimited list of authors for filtering search results.")
        parser.add_argument("--version", action="store", dest="version",
                            help="PyLith version for filtering search results.")
        parser.add_argument("--incompatible", action="store_true", dest="incompatible",
                            help="Filter search results to show incompatible parameter files.")

        args = parser.parse_args()

        return args


# End of file
