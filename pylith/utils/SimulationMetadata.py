# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------

import os

from .converters import string_to_list
from pythia.pyre.components.Component import Component

class SimulationMetadata(Component):
    """Python object for holding simulation metadata.

    FACTORY: simulation_metadata
    """

    import pythia.pyre.inventory

    description = pythia.pyre.inventory.str("description")
    description.meta["tip"] = "Description of simulation."

    authors = pythia.pyre.inventory.list("authors")
    authors.meta["tip"] = "Creator(s) of simulation."

    keywords = pythia.pyre.inventory.list("keywords")
    keywords.meta["tip"] = "Keywords describing simulation."

    features = pythia.pyre.inventory.list("features")
    features.meta["tip"] = "PyLith features used in simulation."

    arguments = pythia.pyre.inventory.list("arguments")
    arguments.meta["tip"] = "Command line arguments for running simulation."

    version = pythia.pyre.inventory.str("version")
    version.meta["tip"] = "Version number for simulation."

    pylith_version = pythia.pyre.inventory.list("pylith_version")
    pylith_version.meta["tip"] = "PyLith versions compatible with simulation input files."

    def __init__(self, name="metadata", facility="simulationmetadata"):
        """Constructor.
        """
        Component.__init__(self, name, facility)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def simulation_metadata():
    """Factory associated with SimulationMetadata.
    """
    return SimulationMetadata()


def fromFile(filename, codec="cfg"):
    """Set simulation metadata from file.

    Args:
        filename (str)
            Name of file
        codec (str)
            Codec for file
    """
    CONVERTERS = {
        "description": str,
        "keywords": string_to_list,
        "features": string_to_list,
        "authors": string_to_list,
        "arguments": string_to_list,
        "version": str,
        "pylith_version": string_to_list,
    }

    base,_ = os.path.splitext(str(filename))
    if codec == "cfg":
        from pythia.pyre.inventory.cfg.CodecConfig import CodecConfig
        reader = CodecConfig()
    else:
        raise NotImplementedError(f"Importing '{codec}' file not implemented.")
    registry = reader.open(base)
    metadata = None
    if registry and \
        "inventory" in registry.keys() and \
            "pylithapp" in registry["inventory"].facilities.keys():
        facilities = registry["inventory"].facilities["pylithapp"].facilities
        properties = facilities["metadata"].properties if "metadata" in facilities else None
        if properties:
            metadata = SimulationMetadata()
            for key, descriptor in properties.items():
                setattr(metadata, key, CONVERTERS[key](descriptor.value))
    return metadata


# End of file
