# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/utils/DumpParametersJson.py
#
# @brief Python DumpParameters object for dumping PyLith parameter information to a JSON file.

from .DumpParameters import DumpParameters


class DumpParametersJson(DumpParameters):
    """
    Python DumpParameters object for dumping PyLith parameter information to a JSON file.

    INVENTORY

    Properties
      - *filename* Name of file written with parameters.
      - *style* Style of JSON file [compact, normal].
      - *indent* Number of spaces to indent, use negative number for no newlines.

    Facilities
      - None
    """

    class Inventory(DumpParameters.Inventory):
        """
        Python object for managing DumpParametersJson facilities and properties.
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
        """
        Constructor.
        """
        DumpParameters.__init__(self, name)
        return

    def write(self, app):
        """
        Write parameters to JSON file.
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

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Configure object.
        """
        DumpParameters._configure(self)
        self.filename = self.inventory.filename
        self.indent = self.inventory.indent
        self.style = self.inventory.style
        return

# FACTORIES ////////////////////////////////////////////////////////////


def dump_parameters():
    """
    Factory associated with DumpParametersJson.
    """
    return DumpParametersJson()


# End of file
