(sec:example:subduction:3d:step01)=
# Step 1: Axial Compression

We start with a very simple example of axial compression in the east-west direction with purely elastic material properties, and no faults ({numref}`fig:example:subduction:3d:step01:diagram`).
We impose axial compression using Dirichlet boundary conditions on the east (+x) and west (-x) boundaries and confine the domain in the north-south direction via zero displacement Dirichlet boundary conditions on the north (+y) and south (-y) boundaries.
We constrain the vertical displacement by imposing zero displacement boundary conditions on the bottom (-z) boundary.

:::{figure-md} fig:example:subduction:3d:step01:diagram
 <img src="figs/subduction3d_step01_diagram.*" alt="Diagram of Step 1 - Axial compression. This static simulation uses Dirichlet boundary conditions with axial compression in the east-west (x-direction), roller boundary conditions on the north, south, and bottom boundaries, and purely elastic properties." width="100%"/>

Diagram of Step 1 - Axial compression. This static simulation uses Dirichlet boundary conditions with axial compression in the east-west (x-direction), roller boundary conditions on the north, south, and bottom boundaries, and purely elastic properties.
:::

The `pylithapp.cfg` file creates an array of five boundary conditions, which impose zero displacements by default.
We overwrite this behavior in the `step01.cfg` file for the -x and +x boundaries using spatial databases with a single uniform displacement value to create the axial compression:

```{code-block} cfg
---
caption: Excerpt from `step01.cfg`
---
# -x face
[pylithapp.problem.bc.x_neg]
db_initial = spatialdata.spatialdb.UniformDB
db_initial.label = Dirichlet BC on -x
db_initial.values = [displacement-x]
db_initial.data = [+2.0*m]

# +x face
[pylithapp.problem.bc.x_pos]
db_initial = spatialdata.spatialdb.UniformDB
db_initial.label = Dirichlet BC on +x
db_initial.values = [displacement-x]
db_initial.data = [-2.0*m]
```

As discussed in {ref}`sec:example:subduction:3d:organization`, we use `mat_elastic.cfg` to specify the parameters associated with linear, isotropic elastic bulk constitutive models for all of the materials for convenient reuse across several different simulations.

```{code-block} cfg
---
caption: Excerpt from `mat_elastic.cfg`
---
[pylithapp.problem.materials]
slab = pylith.materials.ElasticIsotropic3D
wedge = pylith.materials.ElasticIsotropic3D
crust = pylith.materials.ElasticIsotropic3D
mantle = pylith.materials.ElasticIsotropic3D

# Slab
[pylithapp.problem.materials.slab]
db_properties = spatialdata.spatialdb.SimpleDB
db_properties.label = Properties for subducting slab
db_properties.iohandler.filename = spatialdb/mat_slab_elastic.spatialdb

# Wedge
[pylithapp.problem.materials.wedge]
db_properties = spatialdata.spatialdb.SimpleDB
db_properties.label = Properties for accretionary wedge
db_properties.iohandler.filename = spatialdb/mat_wedge_elastic.spatialdb

# Mantle
[pylithapp.problem.materials.mantle]
db_properties = spatialdata.spatialdb.SimpleDB
db_properties.label = Properties for mantle
db_properties.iohandler.filename = spatialdb/mat_mantle_elastic.spatialdb

# Crust
[pylithapp.problem.materials.crust]
db_properties = spatialdata.spatialdb.SimpleDB
db_properties.label = Properties for continental crust
db_properties.iohandler.filename = spatialdb/mat_crust_elastic.spatialdb
```

We specify different elastic properties for each material (slab, wedge, mantle, and crust) using `SimpleDB` spatial databases with a single point to specify uniform properties within a material.
We choose `SimpleDB` rather than `UniformDB`, because we will reuse some of these spatial databases for the elastic properties when we use linear Maxwell viscoelastic constitutive model.

The remaining parameters in the `step01.cfg` file are mostly associated with setting filenames for all of the various output, including all of the parameters used and version information in a JSON file (`output/step01-parameters.json`), a file reporting the progress of the simulation and estimated time of completion (`output/step01-progress.txt`), and the filenames for the HDF5 files (the corresponding Xdmf files will use the same filename with the `xmf` suffix).

```{code-block} shell
---
caption: Run Step 1 simulation
---
$ pylith step01.cfg mat_elastic.cfg

```

The simulation will produce ten pairs of HDF5/Xdmf files in the `output` directory:

**step01-domain.h5[.xmf]**  Time series of the solution field over the domain.

**step01-groundsurf.h5[.xmf]**  Time series of the solution field over the ground surface.

**step01-slab_info.h5[.xmf]**  Properties for the slab material.

**step01-slab.h5[.xmf]**  Time series of the state variables (stress and strain) for the slab material.

**step01-wedge)info.h5[.xmf]**  Properties for the wedge material.

**step01-wedge.h5[.xmf]** Time series of the state variables (stress and strain) for the wedge material.

**step01-crust_info.h5[.xmf]**  Properties for the crust material.

**step01-crust.h5[.xmf]** Time series of the tate variables (stress and strain) for the crust material.

**step01-mantle_info.h5[.xmf]** Properties for the mantle material.

**step01-mantle.h5[.xmf]** Time series of the state variables (stress and strain) for the mantle material.

The HDF5 files contain the data and the Xdmf files contain the metadata required by ParaView and Visit (and other visualization tools that use Xdmf files) to access the mesh and data sets in the HDF5 files.

{numref}`fig:example:subduction:3d:step01`, which was created using the ParaView Python script `plot_dispvec.py` (see {ref}`sec:ParaView:Python:scripts` for how to run ParaView Python scripts), displays the magnitude of the displacement field arrows showing the direction and magnitude of the deformation.
Material properties with a positive Poisson's ratio result in vertical deformation along with the axial compression.
The variations in material properties among the properties result in local spatial variations that are most evident in the horizontal displacement components.

:::{figure-md} fig:example:subduction:3d:step01
<img src="figs/subduction3d_step01_soln.*" alt="Solution over the domain for Step 1. The colors indicate the magnitude of the displacement and the arrows indicate the direction with the length of each arrow equal to 10,000 times the magnitude of the displacement." width="100%"/>

Solution over the domain for Step 1. The colors indicate the magnitude of the displacement and the arrows indicate the direction with the length of each arrow equal to 10,000 times the magnitude of the displacement.
:::

## Exercises

* Run PyLith again and add `solver_algebraicmultigrid.cfg` as an argument on the command line to switch to the algebraic multigrid preconditioner.
    * Using the PETSc log summary to compare the runtime and memory use between the original LU preconditioner and the ML algebraic multigrid preconditioner. Hint: The algebraic multigrid preconditioner is faster.
    * Run the simulation again with the algebraic multigrid preconditioner using multiple cores via the `--nodes=NCORES` argument, replacing `NCORES` with 2 or up to the number of cores on your machine. Examine the PETSc log summary for the various runs to see how the time spent at varies stages changes with the number of cores. Make a plot of runtime versus the number of cores.
* Adjust the material properties in the spatial databases so that the slab is stiffer and the wedge is more compliant. What happens to the solution if you make the materials nearly incompressible? Does this also affect the rate of convergence of the linear solve?
* Change the Dirichlet boundary conditions to impose pure shear instead of axial compression. Hint: You will need to change the boundary conditions on the east, west, north, and south boundaries.
