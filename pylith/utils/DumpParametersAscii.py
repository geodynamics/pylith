# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

from .DumpParameters import DumpParameters


class DumpParametersAscii(DumpParameters):
    """
    Dump PyLith parameter information to an ASCII file.

    Implements `DumpParameters`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp]
            dump_parameters = pylith.utils.DumpParametersAscii

            [pylithapp.dump_parameters]
            filename = output/parameters.txt
            verbose = True
        """
    }

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="pylith_paramters.txt")
    filename.meta["tip"] = "Name of file written with parameters."

    indent = pythia.pyre.inventory.int("indent", default=4)
    indent.meta["tip"] = "Nmber of spaces to indent."

    verbose = pythia.pyre.inventory.bool("verbose", default=True)
    verbose.meta["tip"] = "Include description, location, and aliases."

    def __init__(self, name="dumpparametersascii"):
        """Constructor.
        """
        DumpParameters.__init__(self, name)

    def write(self, app):
        """Write parameters to ASCII file.
        """
        from pylith.mpi.Communicator import mpi_is_root
        if not mpi_is_root():
            return

        if self.info is None:
            self.collect(app)

        parameters = self.info["application"]
        self._createPath(self.filename)
        with open(self.filename, "w") as fout:
            from .CollectVersionInfo import CollectVersionInfo
            import datetime
            now = datetime.datetime.now()
            fout.write("Generated %s\n\n" % now.strftime("%A, %d %B %Y %I:%M%p"))
            fout.write(CollectVersionInfo.asString())
            fout.write("\nApplication: %(name)s %(class)s\n" % parameters)
            depth = 0
            self._writeComponent(fout, parameters, depth + 1)

    def _configure(self):
        """Configure object.
        """
        DumpParameters._configure(self)
        self.filename = self.inventory.filename
        self.indent = self.inventory.indent
        self.verbose = self.inventory.verbose
        self.tab = " " * self.indent

    def _writeComponent(self, fout, obj, depth):
        """Write component parameters to file.
        """
        indent = self.tab * depth
        for (key, item) in obj["properties"].items():
            if self.verbose:
                fout.write("\n")
            fout.write("%s%s (%s) = %s\n" % (indent, key, item["type"], item["value"]))
            if self.verbose:
                indent2 = indent + self.tab
                fout.write("%sDescription: %s\n" % (indent2, item["description"]))
                fout.write("%sSet from: %s\n" % (indent2, item["setFrom"]))

        for (key, item) in obj["components"].items():
            fout.write("\n%s%s = %s (%s)\n" % (indent, key, item["name"], item["class"]))
            if self.verbose:
                indent2 = indent + self.tab
                fout.write("%sConfigurable as: %s\n" % (indent2, ", ".join(item["aliases"])))
                fout.write("%sDescription: %s\n" % (indent2, item["description"]))
                fout.write("%sSet from: %s\n" % (indent2, item["setFrom"]))

            self._writeComponent(fout, item, depth + 1)

# FACTORIES ////////////////////////////////////////////////////////////


def dump_parameters():
    """Factory associated with DumpParametersAscii.
    """
    return DumpParametersAscii()

# End of file
