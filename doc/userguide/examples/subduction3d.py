from examples_info import Features

from collections import OrderedDict
examples = OrderedDict()

examples["3d/subduction/step01"] = {
    "General": {
        "Dimension": 3,
        "Coordinate system": Features.CS_PROJ,
        "Mesh generator": Features.MESH_CUBIT,
        "Cells": Features.CELL_TET,
        "Problem type": Features.PROB_TIMEDEPENDENT,
        "Time dependence": Features.TIMEDEP_STATIC,
        "Reordering": True,
    },
    "Boundary Condition": {
        "Dirichlet": 5,
    },
    "Bulk Rheology": {
        "Linear elastic": 4,
        "Stress/strain formulation": Features.STRAINFORM_INFINITESIMAL,
    },
    "Solver": {
        "Solver": Features.SOLVER_LINEAR,
        "Preconditioner": Features.PRECOND_ILU,
    },
    "Output": {
        "Format": Features.OUTPUT_HDF5,
        "ParaView": True,
        "Domain output": 1,
        "Surface output": 1,
        "State variable output": 4,
    },
    "Spatial Database": {
        "Uniform": 2,
        "Simple": 4,
    },
}

examples["3d/subduction/step02"] = {
    "General": {
        "Dimension": 3,
        "Coordinate system": Features.CS_PROJ,
        "Mesh generator": Features.MESH_CUBIT,
        "Cells": Features.CELL_TET,
        "Problem type": Features.PROB_TIMEDEPENDENT,
        "Time dependence": Features.TIMEDEP_QUASISTATIC,
        "Reordering": True,
    },
    "Boundary Condition": {
        "Dirichlet": 5,
    },
    "Fault": {
        "Prescribed slip": 1,
        "Slip time function": Features.SLIPFN_STEP,
    },
    "Bulk Rheology": {
        "Linear elastic": 2,
        "Linear Maxwell viscoelastic": 2,
        "Stress/strain formulation": Features.STRAINFORM_INFINITESIMAL,
    },
    "Solver": {
        "Solver": Features.SOLVER_LINEAR,
        "Preconditioner": Features.PRECOND_ML_CUSTOM,
        "Time stepping": Features.TS_BWDEULER,
    },
    "Output": {
    "Format": Features.OUTPUT_HDF5,
        "ParaView": True,
        "Domain output": 1,
        "Surface output": 1,
        "State variable output": 4,
    },
    "Spatial Database": {
        "Uniform": 2,
        "Simple": 3,
        "Simple grid": 2,
        "Composite": 2,
    },
}

examples["3d/subduction/step03"] = {
    "General": {
        "Dimension": 3,
        "Coordinate system": Features.CS_PROJ,
        "Mesh generator": Features.MESH_CUBIT,
        "Cells": Features.CELL_TET,
        "Problem type": Features.PROB_TIMEDEPENDENT,
        "Time dependence": Features.TIMEDEP_QUASISTATIC,
        "Reordering": True,
    },
    "Boundary Condition": {
        "Dirichlet": 5,
    },
    "Fault": {
        "Prescribed slip": 2,
        "Slip time function": Features.SLIPFN_RATE,
    },
    "Bulk Rheology": {
        "Linear elastic": 2,
        "Linear Maxwell viscoelastic": 2,
        "Stress/strain formulation": Features.STRAINFORM_INFINITESIMAL,
    },
    "Solver": {
        "Solver": Features.SOLVER_LINEAR,
        "Preconditioner": Features.PRECOND_ML_CUSTOM,
        "Time stepping": Features.TS_BWDEULER,
    },
    "Output": {
        "Format": Features.OUTPUT_HDF5,
        "ParaView": True,
        "Domain output": 1,
        "Surface output": 1,
        "State variable output": 4,
    },
    "Spatial Database": {
        "Uniform": 4,
        "Simple": 3,
        "Simple grid": 2,
        "Composite": 2,
    },
}

examples["3d/subduction/step04"] = {
    "General": {
        "Dimension": 3,
        "Coordinate system": Features.CS_PROJ,
        "Mesh generator": Features.MESH_CUBIT,
        "Cells": Features.CELL_TET,
        "Problem type": Features.PROB_TIMEDEPENDENT,
        "Time dependence": Features.TIMEDEP_QUASISTATIC,
        "Reordering": True,
    },
    "Boundary Condition": {
        "Dirichlet": 5,
    },
    "Fault": {
        "Prescribed slip": 3,
        "Slip time function": Features.SLIPFN_STEP,
    },
    "Bulk Rheology": {
        "Linear elastic": 2,
        "Linear Maxwell viscoelastic": 2,
    "Stress/strain formulation": Features.STRAINFORM_INFINITESIMAL,
    },
    "Solver": {
        "Solver": Features.SOLVER_LINEAR,
        "Preconditioner": Features.PRECOND_ML_CUSTOM,
        "Time stepping": Features.TS_BWDEULER,
    },
    "Output": {
        "Format": Features.OUTPUT_HDF5,
        "ParaView": True,
        "Domain output": 1,
        "Surface output": 1,
        "State variable output": 4,
    },
    "Spatial Database": {
        "Uniform": 7,
        "Simple": 3,
        "Simple grid": 5,
        "Composite": 2,
    },
}

examples["3d/subduction/step05"] = {
    "General": {
        "Dimension": 3,
        "Coordinate system": Features.CS_PROJ,
        "Mesh generator": Features.MESH_CUBIT,
        "Cells": Features.CELL_TET,
        "Problem type": Features.PROB_TIMEDEPENDENT,
        "Time dependence": Features.TIMEDEP_QUASISTATIC,
        "Reordering": True,
    },
    "Boundary Condition": {
        "Dirichlet": 5,
    },
    "Fault": {
        "Prescribed slip": 1,
        "Slip time function": Features.SLIPFN_RATE,
    "Constitutive model": 1,
        "Slip-weakening friction": True,
        "Traction perturbation": True,
    },
    "Bulk Rheology": {
        "Linear elastic": 2,
        "Linear Maxwell viscoelastic": 2,
        "Stress/strain formulation": Features.STRAINFORM_INFINITESIMAL,
    },
    "Solver": {
        "Solver": Features.SOLVER_NONLINEAR,
        "Preconditioner": Features.PRECOND_ML_CUSTOM,
        "Time stepping": Features.TS_BWDEULER,
    },
    "Output": {
        "Format": Features.OUTPUT_HDF5,
        "ParaView": True,
        "Domain output": 1,
        "Surface output": 1,
        "State variable output": 4,
    },
    "Spatial Database": {
        "Uniform": 7,
        "Simple": 3,
        "Simple grid": 5,
        "Composite": 2,
    },
}


examples["3d/subduction/step06"] = {
    "General": {
        "Dimension": 3,
        "Coordinate system": Features.CS_PROJ,
        "Mesh generator": Features.MESH_CUBIT,
        "Cells": Features.CELL_TET,
        "Problem type": Features.PROB_TIMEDEPENDENT,
        "Time dependence": Features.TIMEDEP_QUASISTATIC,
        "Reordering": True,
    },
    "Boundary Condition": {
        "Dirichlet": 5,
    },
    "Fault": {
        "Prescribed slip": 1,
        "Slip time function": Features.SLIPFN_TIMEHISTORY,
    },
    "Bulk Rheology": {
        "Linear elastic": 4,
        "Stress/strain formulation": Features.STRAINFORM_INFINITESIMAL,
    },
    "Solver": {
        "Solver": Features.SOLVER_NONLINEAR,
        "Preconditioner": Features.PRECOND_ML_CUSTOM,
        "Time stepping": Features.TS_BWDEULER,
    },
    "Output": {
        "Format": Features.OUTPUT_HDF5,
        "ParaView": True,
        "Domain output": 1,
        "Surface output": 1,
        "Point output": 1,
        "State variable output": 4,
    },
    "Spatial Database": {
        "Uniform": 1,
        "Simple": 4,
        "Simple grid": 1,
        "Time history": 1,
    },
}


# End of file
