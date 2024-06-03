# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
"""Module for plotting solution and derived fields from PyLith output."""

from pylith.viz import core

class PlotField:
    TIME_DEPENDENT = True

    def __init__(self, plotter, data, interactive=True):
        self.plotter = plotter
        self.data = data
        self.interactive = interactive

    def plot(self, field_name: str, component_name: str):
        self.field_name = field_name
        self.component_name = component_name
        self.lut = core.create_colormap(dataset=self.data, field_name=field_name, component_name=component_name)
        for domain in self.data:
            self._plot_data(domain)
        self.plotter.add_scalar_bar(title=f"{self.component_name} {self.field_name}")

    def update_time(self, t_index: int):
        for domain in self.data:
            self._update_data(domain=domain, t_index=t_index)

    def _plot_data(self, domain):
        self._set_scalar(domain=domain, t_index=0)
        self.plotter.add_mesh(domain.grid, name=domain.name, show_edges=True, show_scalar_bar=False, cmap=self.lut, scalars=self.field_name)

    def _update_data(self, domain, t_index: int):
        self._set_scalar(domain=domain, t_index=t_index)
        self.plotter.render()

    def _set_scalar(self, domain, t_index: int):
        grid_field = domain.get_field(self.field_name)
        scalar_field = grid_field.as_scalar(component_name=self.component_name, t_index=t_index)
        if grid_field.basis_order == 0:
            domain.grid.cell_data[self.field_name] = scalar_field
        elif grid_field.basis_order == 1:
            domain.grid.point_data[self.field_name] = scalar_field


def cli(plotter, data, args):
    figure = PlotField(plotter=plotter, data=data)
    figure.plot(field_name=args.field_name, component_name=args.component)
    return figure


def add_args(parser):
    parser.add_argument("--field", action="store", help="Field for shading. Currently, must be a vertex field.", dest="field_name", default="displacement")
    parser.add_argument("--component", action="store", help="Component of field for shading.", dest="component", default="magnitude", choices=core.Field.COMPONENT.keys())
    parser.set_defaults(func=cli)
