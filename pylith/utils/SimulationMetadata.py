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

import os

from .converters import string_to_list
from pylith.utils.validators import notEmptyString
from pythia.pyre.components.Component import Component


class SimulationMetadata(Component):
    """Python object for holding simulation metadata.

    FACTORY: simulation_metadata
    """

    import pythia.pyre.inventory

    description = pythia.pyre.inventory.str("description", validator=notEmptyString)
    description.meta["tip"] = "Description of simulation."

    authors = pythia.pyre.inventory.list("authors")
    authors.meta["tip"] = "Creator(s) of simulation."

    keywords = pythia.pyre.inventory.list("keywords")
    keywords.meta["tip"] = "Keywords describing simulation."

    features = pythia.pyre.inventory.list("features")
    features.meta["tip"] = "PyLith features used in simulation."

    arguments = pythia.pyre.inventory.list("arguments")
    arguments.meta["tip"] = "Command line arguments for running simulation."

    base = pythia.pyre.inventory.list("base")
    base.meta["tip"] = "Parameter files with metadata that complement this metadata."

    version = pythia.pyre.inventory.str("version")
    version.meta["tip"] = "Version number for simulation."

    pylith_version = pythia.pyre.inventory.list("pylith_version")
    pylith_version.meta["tip"] = "PyLith versions compatible with simulation input files."

    def __init__(self, name="metadata", facility="simulationmetadata"):
        """Constructor.
        """
        Component.__init__(self, name, facility)
        return

    def _validate(self, context):
        if not self.arguments:
            trait = self.inventory.getTrait("arguments")
            self._validationError(context, trait, "List of command line arguments required.")

        if not self.pylith_version:
            trait = self.inventory.getTrait("pylith_version")
            self._validationError(context, trait, "List of PyLith version constraints required.")

    def _validationError(self, context, trait, msg):
        from pythia.pyre.inventory.Item import Item
        error = ValueError(msg)
        descriptor = self.getTraitDescriptor(trait.name)
        context.error(error, items=[Item(trait, descriptor)])


# FACTORIES ////////////////////////////////////////////////////////////

def simulation_metadata():
    """Factory associated with SimulationMetadata.
    """
    return SimulationMetadata()


def fromFile(filename):
    """Set simulation metadata from file.

    Args:
        filename (str)
            Name of file
    """
    def _get_properties(filename):
        """Get metadata from file.
        """
        CONVERTERS = {
            "description": str,
            "keywords": string_to_list,
           "features": string_to_list,
            "authors": string_to_list,
            "arguments": string_to_list,
            "base": string_to_list,
            "version": str,
            "pylith_version": string_to_list,
        }
        if not os.path.isfile(filename):
            raise IOError(f"Could not open file '{filename}' to read metadata.")
        base, ext = os.path.splitext(str(filename))
        if ext == ".cfg":
            from pythia.pyre.inventory.cfg.CodecConfig import CodecConfig
            reader = CodecConfig()
        else:
            raise NotImplementedError(f"Reading '{ext}' file not implemented.")
        registry = reader.open(base)
        properties = None
        if registry and \
            "inventory" in registry.keys() and \
                "pylithapp" in registry["inventory"].facilities.keys():
            facilities = registry["inventory"].facilities["pylithapp"].facilities
            properties = facilities["metadata"].properties if "metadata" in facilities else None
            if properties:
                for key, descriptor in properties.items():
                    descriptor.value = CONVERTERS[key](descriptor.value)
        return properties

    def _merge(target, source):
        """Merge properties.

        Args:
            target (dict)
                Object to update.
            source (dict)
                Object with data to use for update.
        """
        if not source:
            return
        APPEND = ["features", "keywords"]
        for key, descriptor in source.items():
            if descriptor.value:
                if key in APPEND:
                    assert isinstance(descriptor.value, list)
                    target[key].value += descriptor.value
                else:
                    target[key] = source[key]

    metadata = None
    properties = _get_properties(filename)
    if properties:
        metadata = SimulationMetadata()
        if "base" in properties:
            base = properties["base"].value
            baseProperties = None
            for baseFilename in base:
                basePath = os.path.join(filename.parent, baseFilename)
                if not baseProperties:
                    baseProperties = _get_properties(basePath)
                else:
                    _merge(baseProperties, _get_properties(basePath))
            _merge(baseProperties, properties)
            properties = baseProperties

        for key, descriptor in properties.items():
            setattr(metadata, key, descriptor.value)
    return metadata


# End of file
