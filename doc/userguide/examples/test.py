from examples_info import Features

from collections import OrderedDict
examples = OrderedDict({
    "Step01": {
        "General": {
            "Dimension": 2,
            "Coordinate system": Features.CS_CART,
            "Mesh generator": Features.MESH_CUBIT,
            "Cells": Features.CELL_TRI,
            "Problem type": Features.PROB_TIMEDEPENDENT,
            "Time dependence": Features.TIMEDEP_STATIC,
        },
        "Boundary Conditions": {
            "Dirichlet": 2,
            "Neumann": 2,
        },
        "Bulk Rheology": {
            "Linear elastic": 1,
        },
        "Solver": {
            "Solver": Features.SOLVER_LINEAR,
        },
        "Output": {
            "Format": Features.OUTPUT_VTK,
        },
        "Spatial Database": {
            "Simple": 3,
        },
    },
    "Step02": {
        "General": {
            "Dimension": 3,
            "Coordinate system": Features.CS_CART,
            "Mesh generator": Features.MESH_CUBIT,
            "Cells": Features.CELL_TET,
            "Problem type": Features.PROB_TIMEDEPENDENT,
            "Time dependence": Features.TIMEDEP_DYNAMIC,
        },
        "Boundary Conditions": {
            "Absorbing": 5,
        },
        "Bulk Rheology": {
            "Linear elastic": 1,
        },
        "Solver": {
            "Solver": Features.SOLVER_LINEAR,
        },
        "Output": {
            "Format": Features.OUTPUT_HDF5,
        },
        "Spatial Database": {
            "Uniform": 1,
        },
    },
})
