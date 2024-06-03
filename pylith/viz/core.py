# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
"""Module with core functionality for visualizing PyLith output using PyVista."""

from sys import float_info
from dataclasses import dataclass

import numpy
import pyvista

@dataclass
class Field:
    values: numpy.ndarray
    field_type: str
    basis_order: int
    #units: str # Reminder for future use

    COMPONENT = {
        "x": 0,
        "y": 1,
        "z": 2,
        "xx": 0,
        "yy": 1,
        "zz": 2,
        "xy": 3,
        "yz": 4,
        "xz": 5,
        "magnitude": None,
    }

    def as_scalar(self, t_index: int, component_name: str):
        scalar_field = None
        if component_name == "magnitude":
            if self.field_type == "scalar":
                scalar_field = numpy.abs(self.values[t_index,:,0])
            elif self.field_type == "vector":
                scalar_field = numpy.sqrt(self.values[t_index,:,0]**2 + self.values[t_index,:,1]**2 + self.values[t_index,:,2]**2)
            else:
                raise ValueError(f"Getting magnitude for vector field type {self.field_type} is not supported.")
        else:
            component = self.COMPONENT[component_name]
            scalar_field = self.values[t_index, :, component]
        return scalar_field

    def get_range(self, component_name: str):
        min_value = None
        max_value = None
        if component_name == "magnitude":
            min_value = 0
            if self.field_type == "scalar":
                max_value = numpy.abs(self.values).max()
            elif self.field_type == "vector":
                max_value = numpy.sqrt(self.values[:,:,0]**2 + self.values[:,:,1]**2 + self.values[:,:,2]**2).max()
            else:
                max_value = numpy.abs(self.values).max()
        else:
            component = self.COMPONENT[component_name]
            min_value = self.values[:,:,component].min()
            max_value = self.values[:,:,component].max()
        return (min_value, max_value)

@dataclass
class Mesh:
    grid: object
    time: numpy.ndarray
    fields: dict
    name: str

    def get_field(self, name):
        if not name in self.fields:
            raise ValueError(f"Could not find field {name} in data fields {list(self.fields.keys())}.")
        return self.fields[name]


def create_component_colormap(min_value: float, max_value: float, component_name: str):
    cmap = "plasma" if component_name == "magnitude" else "RdBu_r"
    lut = pyvista.LookupTable(cmap=cmap)
    if component_name == "magnitude":
        lut.scalar_range = (min_value, max_value)
    else:
        max_mag = max(abs(min_value), abs(max_value))
        lut.scalar_range = (-max_mag, +max_mag)
    return lut

def create_colormap(dataset, field_name, component_name):
    min_value = float_info.max
    max_value = -float_info.max
    for domain in dataset:
        data_field = domain.get_field(field_name)
        (field_min, field_max) = data_field.get_range(component_name)
        min_value = min(min_value, field_min)
        max_value = max(max_value, field_max)
    return create_component_colormap(min_value, max_value, component_name)
