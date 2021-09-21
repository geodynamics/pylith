# Vertical Cross-Section of Subduction Zone (2D)

PyLith features discussed in this example:

* Static solution
* Quasi-static solution
* CUBIT/Trelis mesh generation w/APREPRO
* Nonplanar geometry
* Variable mesh resolution
*  Linear triangular cells
*  HDF5 output
*  Dirichlet displacement and velocity boundary conditions
*  ZeroDispDB spatial database
*  UniformDB spatial database
*  SimpleDB spatial database
*  SimpleGridDB
*  Multiple materials
*  Nonlinear solver
*  Plane strain linearly elastic material
*  Plane strain linear Maxwell viscoelastic material
*  Prescribed slip
*  Spontaneous rupture
*  Multiple faults
*  Spatially variable coseismic slip
*  Spatially variable aseismic creep
*  Afterslip via fault friction
*  Static friction
*  Slip-weakening friction
*  Rate-state friction

All of the files necessary to run the examples are contained in the directory `examples/2d/subduction`.

## Overview

This example examines quasistatic interseismic and coseismic deformation in 2D for a subduction zone (see {numref}`fig:example:subduction:2d:overview`).
It is based on the 2011 M9.0 Tohoku earthquake off the east coast of Japan.
{numref}`fig:example:subduction:2d:steps` shows the three steps of increasing complexity.
Step 1 focuses on the coseismic slip, Step 2 focuses on interseismic deformation, and Step 3 combines the two into a pseudo-earthquake cycle deformation simulation.
Step 4 focuses on using the change in tractions from Step 1 to construct a simulation with afterslip controlled by frictional sliding.
Steps 5 and 6 replace the prescribed aseismic slip on the subducting slab in Step 2 with a frictional interface, producing spontaneous earthquake ruptures and creep.

:::{figure-md} fig:example:subduction:2d:overview
<img src="figs/cartoon_general.*" alt="Cartoon of subduction zone example." width="100%"/>

Cartoon of subduction zone example.
:::

:::{figure-md} fig:example:subduction:2d:steps
<img src="figs/steps.*" alt="Diagram of fault slip and boundary conditions for each step in the subduction zone example." width="100%" />

Diagram of fault slip and boundary conditions for each step in the subduction zone.
:::

## Mesh Description

We construct the mesh in CUBIT by constructing the geometry, prescribing the discretization, running the mesher, and then grouping cells and vertices for boundary conditions and materials.
We use the APREPRO programming language within the journal files to enable use of units and to set variables for values used many times.
An appendix in the CUBIT documentation discusses the features available with APREPRO in CUBIT.
The CUBIT commands are in three separate journal files.
The main driver is in the journal file `mesh_tri3.jou`.
It calls the journal file `geometry.jou` to construct the geometry and `createbc.jou` to set up the groups associated with boundary conditions and materials.
The journal files are documented and describe the various steps outlined below.

1.  Create the geometry defining the domain.
    1.  Create points.
    2.  Connect points into spline curves.
    3.  Split curves to separate them into sections bounding surfaces.
    4.  Connect curves into surfaces.
    5.  Stitch surfaces together.
2.  Define meshing scheme and cell size variation.
    1.  Define cell size along curves near fault.
    2.  Increase cell size away from fault at a geometric rate (bias).
3.  Generate mesh.
4.  Create blocks for materials and nodesets for boundary conditions.
5.  Export mesh.

:::{figure-md} fig:example:subduction:2d:mesh
<img src="figs/tri3.*" alt="Variable resolution finite-element mesh with triangular cells. The nominal cell size increases at a geometric rate of 1.2 away from the region of coseismic slip." width="100%"/>

Variable resolution finite-element mesh with triangular cells. The nominal cell size increases at a geometric rate of 1.2 away from the region of coseismic slip.
:::

## Common Information

As in previous examples, we place parameters common to the three steps in the `pylithapp.cfg` file so that we do not have to duplicate them for each step.
The settings contained in `pylithapp.cfg` for this problem consist of:

**pylithapp.journal.info** Settings that control the verbosity of the output written to stdout for the different components.

**pylithapp.mesh_generator** Settings that control mesh importing, such as the importer type, filename, and the spatial dimension of the mesh.

**pylithapp.timedependent** Settings that control the problem, such as the total time, time-step size, and spatial dimension.

**pylithapp.timedependent.materials** Settings that control the material type, specify which material IDs are to be associated with a particular material type, and give the name of the spatial database containing the physical properties for the material. The quadrature information is also given.

**pylithapp.problem.formulation.output** Settings related output of the solution over the domain and subdomain (ground surface).

**pylithapp.timedependent.materials.*MATERIAL*.output** Settings related to output of the state variables for material *MATERIAL*.

**pylithapp.petsc** PETSc settings to use for the problem such as the preconditioner type.

The physical properties for each material are specified in spatial database files.
For example, the elastic properties for the continental crust are in `mat_concrust.spatialdb`.
The provided spatial database files all use just a single point to specify uniform physical properties within each material.
A good exercise is to alter the spatial database files with the physical properties to match PREM.

## Example Workflow

:::{toctree}
step01-coseismic-slip.md
step02-interseismic-deformation.md
step03-psuedo-earthquake-cycle-model.md
step04-frictional-afterslip.md
step05-spontaneous-slip-weakening.md
step06-spontaneous-rate-state.md
step07-twofaults-maxwell.md
step08-twofaults-powerlaw.md
:::

## Exercises

:::{admonition} TODO
:class: error

Geophysics.ou.edu link not found: geophysics.ou.edu/solid_earth/prem.html
:::

The list below includes some suggested modifications to these examples that will allow you to become more familiar with PyLith while examining some interesting physics.

*  Change the resolution of the mesh by editing the `mesh_tri3.jou` journal file.
Change the resolution and bias factor.
*  Add depth dependent viscosity to the mantle and crust.
This requires using the linear Maxwell plane strain bulk constitutive model in the crust as well and creating spatial databases that include viscosity for the crust.
Specifying a depth dependent variation in the parameters will require adding points, updating num-locs accordingly, and changing data-dim to 1.
*  Modify the spatial database files for the material properties to use depth-dependent elastic properties based on PREM (Dziewonski and Anderson, 1981, 10.1016/0031-9201(81)90046-7). See <geophysics.ou.edu/solid_earth/prem.html> for a simple table of values. Add points, update num-locs accordingly, and change data-dim to 1.
*  Modify the CUBIT journal files to use quad4 cells rather than tri3 cells.
This requires using the pave mesh scheme.
*  Modify Steps 5 and 6 to use a user-defined variable time step.
Experiment with longer time steps between earthquake ruptures and smaller time steps around the time of the earthquake ruptures.
Can you develop a simple algorithm for choosing the time step?
*  Adjust the parameters of the friction models and examine the effects on the deformation and the convergence of the nonlinear solve.
In which cases do you need to adjust the time step to retain reasonable convergence?
