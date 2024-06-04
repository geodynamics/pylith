# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
"""Module for plotting warped grid from PyLith output."""

import pyvista

from pylith.viz import core

class WarpGrid:
    """Plot warped grid (deformed domain) colored by user-selected field along with outline of
    undeformed domain.
    """
    TIME_DEPENDENT = True

    def __init__(self, plotter: pyvista.Plotter, data: core.Mesh, interactive=True):
        self.plotter = plotter
        self.data = data
        self.interactive = interactive

    def plot(self, exaggeration: float, field_name: str, component_name: str, show_outline: bool=True, show_edges: bool=True):
        self.field_name = field_name
        self.component_name = component_name
        self.exaggeration = exaggeration
        self.show_edges = show_edges

        if show_outline:
            self._create_outlines()

        self.lut = core.create_colormap(dataset=self.data, field_name=field_name, component_name=component_name)
        self.cb_title = f"{component_name} {field_name}"
        for domain in self.data:
            self._create_warped_grid(domain)

        if self.interactive:
            self._add_widgets()

    def update_exaggeration(self, value: float):
        self.exaggeration = value
        for domain in self.data:
            self._warp_grid(domain)

    def update_time(self, t_index: int):
        for domain in self.data:
            self._update_data(domain=domain, t_index=t_index)
            self._warp_grid(domain)

    def _update_data(self, domain, t_index: int):
        self._set_scalar(domain=domain, t_index=t_index)
        domain.grid.point_data["displacement-vector"] = domain.get_field("displacement").values[t_index,:,:].squeeze()
        self.plotter.render()

    def _create_outlines(self):
        self.outlines = [self.plotter.add_mesh(domain.grid, style="wireframe", color="gray") for domain in self.data]

    def _create_warped_grid(self, domain):
        t_index = 0
        self._set_scalar(domain=domain, t_index=t_index)
        domain.grid.point_data["displacement-vector"] = domain.get_field("displacement").values[t_index,:,:].squeeze()
        
        warped_grid = domain.grid.warp_by_vector("displacement-vector", factor=self.exaggeration)
        cb_args = {"title": self.cb_title}
        self.plotter.add_mesh(warped_grid, name=domain.name, show_edges=self.show_edges, cmap=self.lut, scalars=self.field_name, scalar_bar_args=cb_args)

    def _warp_grid(self, domain):
        warped_grid = domain.grid.warp_by_vector("displacement-vector", factor=self.exaggeration)
        cb_args = {"title": self.cb_title}
        self.plotter.add_mesh(warped_grid, name=domain.name, show_edges=self.show_edges, cmap=self.lut, scalars=self.field_name, scalar_bar_args=cb_args)

    def _set_scalar(self, domain, t_index: int):
        grid_field = domain.get_field(self.field_name)
        scalar_field = grid_field.as_scalar(component_name=self.component_name, t_index=t_index)
        if grid_field.basis_order == 0:
            domain.grid.cell_data[self.field_name] = scalar_field
        elif grid_field.basis_order == 1:
            domain.grid.point_data[self.field_name] = scalar_field

    def _add_widgets(self):
        self.plotter.add_slider_widget(self.update_exaggeration, rng=[1, 5*self.exaggeration], pointa=(0.69, 0.92), pointb=(0.97, 0.92), value=self.exaggeration, title="Exaggeration")


def cli(plotter, data, args):
    figure = WarpGrid(plotter=plotter, data=data)
    figure.plot(exaggeration=args.exaggeration, field_name=args.field_name, component_name=args.component, show_outline=args.show_outline, show_edges=args.show_edges)
    return figure


def add_args(parser):
    parser.add_argument("--field", action="store", help="Field for shading. Currently, must be a vertex field.", dest="field_name", default="displacement")
    parser.add_argument("--component", action="store", help="Component of field for shading.", dest="component", default="magnitude", choices=core.Field.COMPONENT.keys())
    parser.add_argument("--exaggeration", action="store", help="Exaggeration factor for displacement.", dest="exaggeration", default=1000.0, type=float)
    parser.add_argument("--hide-outline", action="store_false", help="Hide outline of undeformed domain.", dest="show_outline")
    parser.add_argument("--hide-edges", action="store_false", help="Hide cell edges of deformed domain.", dest="show_edges")
    parser.set_defaults(func=cli)
