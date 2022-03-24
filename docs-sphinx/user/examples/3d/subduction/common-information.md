(sec:example:subduction:3d:organization)=
# Organization of Simulation Parameters

PyLith automatically reads in `pylithapp.cfg` from the current directory, if it exists.
As a result, we generally put all parameters common to a set of examples in this file to avoid duplicating parameters across multiple files.
Because we often use a single mesh for multiple simulations in a directory, we place all parameters related to our mesh and identifying the materials in our mesh in `pylithapp.cfg`.
We assign the bulk constitutive model and its parameters to each material in other files, because we vary those across the simulations.
In general, we place roller boundary conditions (Dirichlet boundary conditions constraining the degrees of freedom perpendicular to the boundary) on the lateral and bottom boundaries, so we include those in `pylithapp.cfg`.
In some simulations we will overwrite the values for parameters will values specific to a given example.
This file is also a convenient place to put basic solver parameters and to turn on Pyre journals for displaying informational messages during a run and journalling debugging flags.

Hence the settings contained in `pylithapp.cfg` include:

 **pylithapp.journal.info** Settings that control the verbosity of the output written to stdout for the different components.

 **pylithapp.mesh_generator** Parameters for the type of mesh importer (generator), reordering of the mesh, and the mesh coordinate system.

 **pylithapp.problem.materials** Basic parameters for each of the four materials, including the label, block id in the mesh file, discretization, and output writer.

 **pylithapp.problem.bc** Parameters for Dirichlet boundary conditions on the lateral and bottom boundaries of the domain.

 **pylithapp.problem.formulation.output** Settings related output of the solution over the domain and subdomain (ground surface).

 **pylithapp.petsc** PETSc solver and logging settings.

### Coordinate system

We generated the mesh in a Cartesian coordinate system corresponding to a transverse Mercator projection.
We specify this geographic projection coordinate system in the `pylithapp.cfg` file, so that we can use other convenient georeferenced coordinate systems in the spatial databases.
PyLith will automatically transform points between compatible coordinate systems.
Our spatialdata library uses Proj4 for geographic projections, so we specify the projection using Proj4 syntax in the **proj_options** property:

```{code-block} cfg
---
caption: Excerpt from `pylithapp.cfg`
---
[pylithapp.mesh_generator.reader]
coordsys = spatialdata.geocoords.CSGeoProj
coordsys.space_dim = 3
coordsys.datum_horiz = WGS84
coordsys.datum_vert = mean sea level
coordsys.projector.projection = tmerc
coordsys.projector.proj_options = +lon_0=-122.6765 +lat_0=45.5231 +k=0.9996
```

## Materials

The finite-element mesh marks cells for each material and the type of cell determines the type of basis functions we use in the discretization.
This means we can specify this information in the `pylithapp.cfg` file and avoid duplicating it in each simulation parameter file.
To set up the materials, we first create an array of materials that defines the name for each material component.
For example, we create the array of four materials and then set the parameters for the slab:

```{code-block} cfg
---
caption: Excerpt from `pylithapp.cfg`
---
[pylithapp.problem]
materials = [slab, wedge, crust, mantle]

[pylithapp.problem.materials.slab]
label = Subducting slab ; Label for informative error messages
id = 1 ; Block id in ExodusII file from CUBIT/Trelis
quadrature.cell = pylith.feassemble.FIATSimplex ; Tetrahedral cells
quadrature.cell.dimension = 3

# Average cell output over quadrature points, yielding one point per cell
output.cell_filter = pylith.meshio.CellFilterAvg
output.writer = pylith.meshio.DataWriterHDF5 ; Output using HDF5
```

In this set of examples, we will consider cases in which all materials are linear, isotropic elastic and cases where the crust and wedge are linear, isotropic elastic but the slab and mantle are linear Maxwell viscoelastic.
As a result, we put the parameters for these two cases in separate `cfg` files with `mat_elastic.cfg` for the case with purely elastic models and `mat_viscoelastic.cfg` for the case with a mix of elastic and viscoelastic models.
Each of these files specifies the bulk constitutive model and spatial database to use for the properties for each material.
The values for the material properties are loosely based on a 3-D seismic velocity model for the Pacific Northwest {cite}`Stephenson:2007`.

### Boundary Conditions

For the Dirichlet boundary conditions, we specify the degree of freedom constrained, the name of the nodeset in the ExodusII file from CUBIT/Trelis that defines the boundary, and a label for the spatial database (required for informative error messages).
These settings constrain the y-displacement on the north (+y) boundary:

```{code-block} cfg
---
caption: Excerpt from `pylithapp.cfg`
---
[pylithapp.problem.bc.y_pos]
bc_dof = [1] ; Degree of freedoms are: x=0, y=1, and z=2
label = boundary_ypos ; nodeset in ExodusII file form CUBIT/Trelis
db_initial.label = Dirichlet BC on +y ; label for informative error messages
```

## Solver Parameters

We group solver parameters into a few different files to handle different cases.
The `pylithapp.cfg` contains tolerance values for the linear and nonlinear solvers and turns on some simple diagnostic information.
The file also directs PyLith to use a direct solver, which is suitable for debugging and test problems that do not include a fault; a direct solver is not well-suited for production runs because it does not scale well and uses a lot of memory.

```{code-block} cfg
---
caption: Excerpt from `pylithapp.cfg`
---
[pylithapp.petsc]
malloc_dump = ; Dump information about PETSc memory not deallocated.

# Use LU preconditioner (helpful for learning and debugging, not production simulations)
pc_type = lu

# Convergence parameters.
ksp_rtol = 1.0e-10 ; Converge if residual norm decreases by this amount
ksp_atol = 1.0e-11 ; Converge if residual norm drops below this value
ksp_max_it = 500 ; Maximum number of iterations in linear solve
ksp_gmres_restart = 50 ; Restart orthogonalization in GMRES after this number of iterations

# Linear solver monitoring options.
ksp_monitor = true ; Show residual norm at each iteration
#ksp_view = true ; Show solver parameters (commented out)
ksp_converged_reason = true ; Show reason linear solve converged
ksp_error_if_not_converged = true ; Generate an error if linear solve fails to converge

# Nonlinear solver monitoring options.
snes_rtol = 1.0e-10 ; Converge if nonlinear residual norm decreases by this amount
snes_atol = 1.0e-9 ; Converge if nonlinear residual norm drops below this value
snes_max_it = 100 ; Maximum number of iterations in nonlinear solve
snes_monitor = true ; Show nonlinear residual norm at each iteration
snes_linesearch_monitor = true ; Show nonlinear solver line search information
#snes_view = true ; Show nonlinear solver parameters (commented out)
snes_converged_reason = true ; Show reason nonlinear solve converged
snes_error_if_not_converged = true ; Generate an error if nonlinear solve fails to converge

# PETSc summary -- useful for performance information.
log_view = true
```

The `solver_algebraicmultigrid.cfg` provides more optimal settings for simulations without a fault by using an algebraic multigrid preconditioner.
Similarly, for simulations with a fault `solver_fieldsplit.cfg` provides settings for applying the algebraic multigrid preconditioner to the elasticity portion of the system Jacobian matrix and our custom fault preconditioner to the Lagrange multiplier portion.
