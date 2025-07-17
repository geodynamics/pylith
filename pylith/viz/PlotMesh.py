# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
"""Module for plotting mesh and mesh quality from PyLith output."""

class PlotMesh:
    TIME_DEPENDENT = False

    def __init__(self, plotter, data, interactive: bool=True):
        self.plotter = plotter
        self.data = data
        self.interactive = interactive

    def plot(self, quality_metric: str):
        for domain in self.data:
            self._plot_domain(domain, quality_metric)

    def _plot_domain(self, domain, quality_metric):
        mesh_args = {
            "show_edges": True,
            "cmap": "plasma_r",
        }
        if quality_metric:
            quality = domain.grid.compute_cell_quality(quality_measure=quality_metric)
            self.plotter.add_mesh(domain.grid, style="wireframe", color="gray")
            if self.interactive:
                self.plotter.add_mesh_threshold(quality, method="upper", **mesh_args)
            else:
                self.plotter.add_mesh(quality, **mesh_args)
        else:
            self.plotter.add_mesh(domain.grid, show_edges=True)


def cli(plotter, data, args):
    figure = PlotMesh(plotter=plotter, data=data)
    figure.plot(quality_metric=args.quality_metric)
    return figure


def add_args(parser):
    QUALITY_METRICS = (
        "area",
        "aspect_beta",
        "aspect_frobenius",
        "aspect_gamma",
        "aspect_ratio",
        "collapse_ratio",
        "condition",
        "diagonal",
        "dimension",
        "distortion",
        "jacobian",
        "max_angle",
        "max_aspect_frobenius",
        "max_edge_ratio",
        "med_aspect_frobenius",
        "min_angle",
        "oddy",
        "radius_ratio",
        "relative_size_squared",
        "scaled_jacobian",
        "shape",
        "shape_and_size",
        "shear",
        "shear_and_size",
        "skew",
        "stretch",
        "taper",
        "volume",
        "warpage",
    )
    parser.add_argument("--mesh-quality", action="store", help="Show mesh quality metric.", dest="quality_metric", choices=QUALITY_METRICS)
    parser.set_defaults(func=cli)
