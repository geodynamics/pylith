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
# @file pylith/utils/DumpParametersJson.py
#
# @brief Python DumpParameters object for dumping PyLith parameter information to a JSON file.

from .DumpParameters import DumpParameters


class DumpParametersJson(DumpParameters):
    """Python DumpParameters object for dumping PyLith parameter information to a JSON file.

    FACTORY: dump_parameters
    """

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="pylith_parameters.json")
    filename.meta["tip"] = "Name of file written with parameters."

    style = pythia.pyre.inventory.str("style", default="normal", validator=pythia.pyre.inventory.choice(["normal", "compact"]))
    style.meta['tip'] = "Style of JSON file [compact, normal]."

    indent = pythia.pyre.inventory.int("indent", default=4)
    indent.meta['tip'] = "Nmber of spaces to indent, use a negative number for no newlines."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="dumpparametersjson"):
        """Constructor.
        """
        DumpParameters.__init__(self, name)
        return

    def write(self, app):
        """Write parameters to JSON file.
        """
        if self.info is None:
            self.collect(app)

        self._createPath(self.filename)
        with open(self.filename, "w") as fout:
            import json
            if self.style == "compact":
                indent = None
                separators = (",", ":")
            elif self.style == "normal":
                indent = self.indent
                separators = (",", ": ")
            else:
                raise ValueError("Unknown JSON style '%s'." % self.style)

            json.dump(self.info, fout, indent=indent, separators=separators)
        return

# FACTORIES ////////////////////////////////////////////////////////////

def dump_parameters():
    """Factory associated with DumpParametersJson.
    """
    return DumpParametersJson()


# End of file
