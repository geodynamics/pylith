See <https://github.com/geodynamics/pylith/commits/main> for the complete log of changes made to PyLith.

## Version 3.1.0

* Removed support for importing meshes from LaGriT.
* Add output of change in fault tractions for prescribed slip.
* State variables are now included in the default `data_fields` for simulation output.
* By default use PETSc proper orthogonal decomposition (POD) methodology for initial guess of solutions to improve convergence.
* Add demonstration of `pylith_powerlaw_gendb` in Step 8 of `examples/reverse-2d`.
* Switched from CppUnit to Catch2 as the C++ testing framework.
* Update to PETSc 3.19.5
* Improve integration with VSCode for testing and debugging (see Developer Guide)
* Bug fixes
  * Fix errors in KinSrcTimeHistory.py
  * Fix creation of PETSc label for edges when importing Gmsh files. This fixes creation of faults with buried edges for 3D meshes imported from Gmsh.

## Version 3.0.3

This is a bug fix release with no new features or changes to the user interface.

* Fixed duplicate integration of fault terms if a fault had one material on one side and multiple materials on the other side.
* Fixed bugs related to running in parallel.
  * Creating constraints on buried fault edges failed for some mesh distribution cases.
  * Green's function problems did not manage fault impulses on multiple processes.
  * Creating a point mesh for `OutputSolnPoints` failed when running in parallel.
  * PetscSF inconsistencies generated errors at various times when running in parallel.
* Update to PETSc 3.18.0.

**Note**: We now use PETSc routines to write the HDF5 files. As a result, there is one change to the layout: `topology/cells` is now `viz/topology/cells`.
The corresponding Xdmf files reflect this change.

### Binary packages

* Update to Python 3.10.6.
* Use `gmforker` process manager with MPICH to avoid localhost name issues.

## Version 3.0.2

This is a bug fix release with no new features or changes to the user interface.

* Add check of PyLith version against version requirements specified in metadata of parameter files.
* Update defaults to better match most use cases.
  * Use nonlinear solver.
  * Basis order is 1 for solution fields.
  * Basis order is 0 for Cauchy stress and strain.
  * Use ML algebraic multigrid preconditioner (from Trilinos) instead of GAMG preconditioner for more robust solves. This is a temporary change until we find better GAMG settings.
* Update PETSc to v3.17.3.
* Remove obsolete LaTeX documentation.
* Bug fixes
  * Add `viz` directory missing from `examples/subduction-2d` in source distribution.
  * Project output fields using correct PETSc routine (`DMProjectFieldLabel()`). Fixes memory access bugs in both serial and parallel.
  * Fix build warnings.
  * Fix reordering that causes errors when importing Gmsh files.
* Documentation
  * Add discussion of translating boundary value problem information to parameter settings. Add more code blocks to manual.
  * Add discussion of `examples/troubleshooting-2d` to manual.

### Binary packages

* Added PyQT5 Python module for interactive plotting with matplotlib.
* Update PyLith Parameter Viewer to v2.0.1 (fix errors in packaging).

### Known issues

* The default PETSc options provide a computationally expensive preconditioner when solving incompressible elasticity problems in parallel. We expect to have a more optimal preconditioner in the next release.
* You may still encounter a few bugs when running in parallel; they appear to cases with specific partitioning of the mesh in relation to one or more faults.

## Version 3.0.1

This is a bug fix release with no new features or changes to the user interface.

* Bug fixes
  * Fix lots of small bugs related to running in parallel
  * Fix several discrepancies among the code, examples, and manual
* Examples
  * Added `examples/subduction-3d` steps 1-4 (included in the manual)
  * Added `examples/troubleshooting-2d` (included in the PyLith v3.0 tutorials but not yet added to the manual)
* Documentation
  * Added instructions for how to remove Apple quarantine attributes
  * Fix LaTeX build of documentation (now available at https://pylith.readthedocs.io)
  * Improved instructions on how to run ParaView Python scripts when starting ParaView from a shortcut
  * Added notes indicating steps of examples are not yet updated for v3.0
  * Fix lots of typos

### Binary packages

* Updated PyLith Parameter Viewer (v2.0.0) for Python 3.

### Known issues

* The default PETSc options provide a computationally expensive preconditioner when solving incompressible elasticity problems in parallel. We expect to have a more optimal preconditioner in the next release.
* You may still encounter a few bugs when running in parallel; they appear to cases with specific partitioning of the mesh in relation to one or more faults.

## Version 3.0.0

Version 3.0.0 includes major changes to the underlying finite-element formulation and implementation in order to support a more flexible specification of the governing equations and higher order basis functions.
These changes affect how simulations are defined.
Parameter files for previous versions will need to be updated; the changes are too complex for a simple translation table.
Some features present in v2.2.2, such as spontaneous rupture and finite strain, have not yet been implemented in the new formulation.

* Multiphysics
  * Elasticity for linear isotropic materials and linear Maxwell, generalized Maxwell, and power law viscoelastic models
  * Incompressible elasticity for linear isotropic materials
  * Prescribed slip for quasistatic and dynamic simulations
* Higher order basis functions
    Allow user to select order of basis functions independent of the mesh (which defines the geometry). This permits higher resolution for a given mesh.
* Switch to using PETSc time-stepping (TS) algorithms
  Replace simple Python-based time-stepping implementations with PETSc time-stepping algorithms that provide support for higher order discretization in time and real adaptive time stepping.
* Static Green's functions with user-specified discretization of fault slip impulses
* Import finite-element meshes from Cubit (Exodus II), Gmsh, and LaGriT
* Modular approach for initial conditions
* Output of subfields with user-defined basis order
* Simulation metadata with command line utility for searching metadata
* Convert to Python 3
* Convert LaTeX documentation to Sphinx + MyST
* Testing with the Method of Manufactured Solutions
* Automatically assign label value for fault cohesive cells (`id` setting is obsolete).
* Use `description` for descriptive labels and `label` and `label_value` for tagging entities. PyLith's use of`label` and `label_value` now corresponds to PETSc labels and label values.

### Deprecated features

* We plan to discontinue support for reading LaGriT mesh files in version 3.2. \
  Gmsh provides an open-source alternative with a graphical user interface.

### SpatialData settings

```{code-block} cfg
db = spatialdata.spatialdb.UniformDB

# Old
db.label = Slip spatial database

# New
db.description = Slip spatial database
```

### Material settings

```{code-block} cfg
material = pylith.materials.Elasticity

# Old
material.label = Elastic material
material.id = 2

# New
material.description = Elastic material
material.label_value = 2
```

### Known issues

* Running in parallel has a few minor bugs due to communication mismatches and over-aggressive error checking. We will be fixing these in the next week.
* We will be updating the 3D subduction zone example (examples/subduction-3d) to v3.0.0 in the next week, including providing the input mesh file; in the meantime do not attempt to run this example.
* We have included Gmsh in the binary packages.
For Linux there are additional libraries that must be installed for Gmsh to run; these are associated with the graphical user interface and included in most default installations.

### Contributors

* Brad Aagaard
* Matthew Knepley
* Charles Williams
* Robert Walker
* Chris Mills
* Shengduo Liu
* Thea Ragon
* Alex Berne
* Jed Brown
* Rey Koki
* Kali Allison
* Lorraine Hwang

## Version 2.2.2

* Bug fixes
  * Fix several typos in the manual.
  * Fix order of deallocation of data members in ElasticMaterial to prevent a segmentation fault, thereby allowing error messages to be shown when throwing an exception.
  * Fix tests for MPI and PETSc version info for more use cases.
  * Ensure all Python script are executable and use nemesis is used instead of python for correct paths to modules on Darwin.
* Added ability to write residual to an HDF5 file during solves. This is intended for debugging and is enabled at runtime using `--journal.debug.formulation=1`. The residual will be written to `residual.h5`. To generate the associated `.xdmf` file run `pylith_genxdmd -f residual.h5`.
* Updated to PETSc 3.10.2

## Version 2.2.1

* Added new examples.

  * **examples/3d/subduction**: New suite of examples for a 3-D subduction zone. This intermediate level suite of examples illustrates a wide range of PyLith features for quasi-static simulations.

  * **examples/2d/subduction**: Added quasi-static spontaneous rupture earthquake cycle examples (Steps 5 and 6) for slip-weakening and rate- and state-friction.

  * These new examples make use of ParaView Python scripts to facilitate using ParaView with PyLith.

* Improved the PyLith manual

  * Added diagram to guide users on which installation method best meets their needs.

  * Added instructions for how to use the Windows Subsystem for Linux to install the PyLith Linux binary on systems running Windows 10.

* Fixed bug in generating Xdmf files for 2-D vector output. Converted Xdmf generator from C++ to Python for more robust generation of Xdmf files from Python scripts.

* Updated spatialdata to v1.9.10. Improved error messages when reading SimpleDB and SimpleGridDB files.

* Updated PyLith parameter viewer to v1.1.0. Application and documentation are now available on line at <https://geodynamics.github.io/pylith_parameters>. Small fix to insure hierarchy path listed matches the one for PyLith.

* Updated PETSc to v3.7.6. See the PETSc documentation for a summary of all of the changes.

* Switched to using CentOS 6.9 for Linux binary builds to insure compatibility with glibc 2.12 and later.

## Version 2.2.0

* Added a browser-based parameter viewer for interactive viewing of all PyLith parameters and version information. See Section 4.10 PyLith Parameter Viewer of the PyLith user manual.

* Adjusted packaging of the binary distributions so that they can be used to extend PyLith and/or integrate other code with PyLith.

* Converted the user manual from Lyx to LaTeX and added syntax highlighting of parameter and spatial database files. Fixed several typos.

* Fixed bug that sometimes resulted in an inconsistent fault orientation when running in parallel. The bug appears to have been introduced in v2.0.

* Fixed two bugs in output of solution at points that sometimes
  happened in parallel simulations. The errors include:

  * The order of the station names does not match the order of the points. The point data is written in parallel by process order, so the points for process 0 are written first, then those for process 1, etc. This often results in reordering of the points. The station names were written in the original order.

  * The output values for some points are incorrect. The wrong cells were being used in the interpolation.

* Updated PETSc to v3.7.5.


## Version 2.1.4

* Added --version command line argument to display version information for PyLith and its dependencies.

* Improved information displayed with the --help command line argument.

* Added `--include-citations` command line argument to display publications to cite when publishing results from computations using PyLith. General PyLith references are also displayed with the `--version` command line argument.

* Allow use of NetCDF versions greater than 4.1.3. Switch from using C++ API to C API.

* Fixed bug in Pythia associated with validation of parameters being done before help could be displayed.

* Fixed typos in manual for gravity and point forces.

* Added integration with Travis for automated testing.


## Version 2.1.3

* Add `generate_statedb.py` and `postseismic.pvsm` files missing from `examples/2d/gravity`.

* Update handling of fault intersection when creating boundary
condition nodesets in `examples/meshing/surface_nurbs/subduction`.

* Fixes to Darwin binary package.

  * Fix linking of netCDF4 Python module.

  * Fix linking and executable mode permissions for Python scripts in binary by using nemesis so relative links are valid.

## Version 2.1.2

* Bugfixes for finite-strain formulation.

  * Added output of the Cauchy stresses (cauchy_stress). The second Piola-Kirchoff stresses are output via the stress field.

  * Material properties and state variables were not retrieved properly when updating state variables.

* Bugfixes for setting initial stress and state variables for viscoelastic materials. The deviatoric stress state is carried forward using the state variables, so the initial deviatoric stress should not be considered when computing the stresses.

* Created new examples showing how to use gravity, initial stress, and finite-strain in 2-D (examples/2d/gravity).

* Reintroduced check (that had been inadvertently removed in v2.x) for ambiguous description of fault surface based on groups of vertices defining faces of cells.

* Flush the output of the progress monitor so progress reports are updated promptly.

* Updates to the user manual.

  * Added section on the debugging examples covered in recent tutorials.

  * Added tables describing the spatial database values for each material.

  * Included a more complete discussion of the finite-strain formulation.

* PETSc

  * Updated to PETSc v3.7.2 (`knepley/pylith` branch).

  * Fixed Trilinos/ML configuration and code so that it can be built without a Fortran compiler.


## Version 2.1.0

* Station names are required for output at arbitrary points (`OutputSolnPoints`) and are included in a /stations dataset in HDF5 files.

* A progress monitor will update a text file with the progress of a simulation (time in the time stepping loop or the number of impulses completed) and given an estimate of when the simulation will be completed.

### Bug fixes

* A few bugs related to creating cohesive cells for fault intersections have been fixed. Faults can now meet at T intersections provided the buried edges of the faults are clamped. In other words, the fault ending at the T intersection has a clamped edge along the intersection. The fault ending at the intersection must also come AFTER the through-going fault in the list of fault interfaces.

* There have been two major bug fixes for Drucker-Prager plasticity, for both DruckerPrager3D and DruckerPragerPlaneStrain. The first fix was a missing initial pressure term for the plastic multiplier in the Drucker-Prager formulation. This affects plasticity calculations when initial stresses are used. The error has been corrected in the code, the manual, and the unit tests. The second bug was an incorrect test for tensile yield that could cause PyLith to exit with an error when plastic yield had not actually occurred. The error would only occur when the allow_tensile_yield flag was set to False. This bug has been fixed in the code, and the new test is also described in the manual. This should prevent problems that previously existed when allow_tensile_yield was set to False (as it should be for most quasi-static problems).

* Fixed bug in DataWriterHDF5Ext associated with multiple processes writing information to the HDF5 file. With external datasets the HDF5 file is limited to metadata and is maintained by process 0.

* A two-dimensional gravity example has been added, based on the tutorial from the June, 2014 workshop at Stanford University.  The tutorial itself is in examples/2d/gravity, and a new section has also been added to the manual describing the example.

* Fixed inconsistent fault orientation when running in parallel for 2-D domains.

### Migrating from v2.0.x to v2.1.x

The points file for `OutputSolnPoints` must now contain station names as the first column.


## Version 2.0.3

### Bug fixes

* Updated autotools files (Makefile.am, configure.ac) for compatibility with recent versions of automake (up to and including v1.14.1).

## Version 2.0.2

### Bug fixes

* Fixed linking issue in Darwin binary distribution, primarily affecting systems with OS X 10.7 and 10.8.

* Improved example journal files for CUBIT/Trelis to improve compatibility (examples/meshing/surface_nurbs/dem).

* Updated more journal in examples so that APREPRO lines have a leading '$' instead of a '#' to differentiate from comments.

* Added examples/debugging files from Crustal Deformation Modeling workshop debugging tutorial.

## Version 2.0.1

### Bug fixes

* Improved example journal files for CUBIT/Trelis to improve compatibility. All journal files should work with CUBIT 14.1 and Trelis 15.0.

* Created examples of IDless journal files in examples/2d/greensfns.  These files should work with all recent versions of CUBIT and Trelis.

* Switched journal APREPRO lines to have leading '$' instead of '#' to differentiate from comments.


## Version 2.0.0

* Replaced C++ Sieve implementation of finite-element data structures with C DMPlex implementation.

  DMPlex provides a simpler, more efficient implementation of the finite-element data structures that conforms to the PETSc data management (DM) interface. This provides tighter integration with the rest of PETSc. Additionally, this rewrite of the data structures results in a more efficient memory layout, resulting in better performance.

* Improved treatment of buried fault edges, so that the slip naturally tapers to zero along the buried edges.

  An additional nodeset/pset is used to designate the buried edges of a fault. This allows the cohesive cells to be inserted up to the edge of the fault without splitting the mesh at the fault edge. The slip will naturally taper to zero at along the buried edges as a result of how the cohesive cells are created.

* Switched from using Subversion to Git for version control.

  The source code repository changed from a CIG maintained Subversion repository to a Git repository at GitHub.com. The URL for the Git repository is <https://github.com/geodynamics/pylith>. The installer has been updated accordingly.

* Added ability to recursively refine a mesh.

  Global uniform refinement can now be done recursively. Each refinement reduces the vertex spacing by a factor of two. Using more than one level of refinement should be done carefully as the mesh quality will generally deteriorate with more levels of refinement.

* Directories for output are created as necessary.

  Directories where output files are written will be created if necessary. Previously, the directories would not be created, so that opening the output files in a nonexistent directory would generate an error.

* Improved error messages.

  Error messages originating in PETSc will include a stack trace that includes both PyLith and PETSc code. Previously, only the PETSc code was included. This provides significantly more information for debugging.

* Improved CUBIT example for mesh sizing functions.

  Based on experimentation with CUBIT 14.0, 14.1, and Trelis 15.0, we have improved the CUBIT mesh sizing examples (`examples/meshing/cubit_cellsize`). We were able to simplify the journal files and use fewer CUBIT commands. The new procedure also eliminates some CUBIT warnings.

* Several small improvements to various sections of the manual based
  on feedback and questions from users.

  * Added more information about the workflow involved in using PyLith.

  * Added a discussion of how to set scales for nondimensionalization.

  * Added a discussion of how the stable time step is computed for the various materials.

  * Updated and expanded the discussion of using initial state variables.

### Bug fixes

* Fixed two MPI related bugs in computing Green's functions in parallel. The number of impulses corresponded to only those on process 0.

* Corrected computation of fault tractions (Lagrange multipliers) on process boundaries for prescribed slip with explicit time stepping.

* Fixed bug when reading in list of output points with just one point.

* Adjusted autoconf Python setup macro to remove temporary sysconfig.pyc file.

* Added check to make sure degree of freedom specified in Dirichlet BC is consistent with spatial dimension of problem.

* Corrected two typos in the manual related to fault opening and tractions in `examples/3d/hex8/step20` and updating to the use of cell.dimension for the quadrature scheme with tets.

* Fixed stable time step computed for power-law viscoelastic rheology to match manual.


### Migrating from v1.9.x to 2.0.x

Changes to various C++ objects permitted simplifying the specification
of a number of components. The map below indicates the name changes.

```
CellFilterAvgMesh -> CellFilterAvg
CellFilterAvgSubMesh -> CellFilterAvg
DataWriterVTKMesh -> DataWriterVTK
DataWriterVTKSubMesh -> DataWriterVTK
DataWriterVTKSubSubMesh -> DataWriterVTK
DataWriterHDF5Mesh -> DataWriterHDF5
DataWriterHDF5SubMesh -> DataWriterHDF5
DataWriterHDF5SubSubMesh -> DataWriterHDF5
DataWriterHDF5ExtMesh -> DataWriterHDF5Ext
DataWriterHDF5ExtSubMesh -> DataWriterHDF5Ext
DataWriterHDF5ExtSubSubMesh -> DataWriterHDF5Ext
```

Running the script:

```
bash $PYLITH_DIR/doc/developer/update_1.9to2.0.sh
```

will update all .cfg files in the current directory and all subdirectories with the new names (you will need to replace $PYLITH_DIR with the directory containing the PyLith source code).


PyLith allows use of the Chaco and ParMetis/Metis partitioners. The name of the ParMetis/Metis partitioner was changed from "parmetis" to "metis".

```ini
[pylithapp.mesh_generator]
distributor.partitioner = metis
```

Buried edges of faults are handled differently in v2.0. A separate nodeset/pset should be created and contain the vertices on the buried edges of the fault. See the Section 6.4.2 of the PyLith manual for more information.

## Version 1.9.0

### New features

* Added Newton-Raphson algorithm for spontaneous rupture simulations with explicit-stepping.

* Enforcing the friction criterion in a spontaneous rupture simulation with explicit time-stepping now uses a Newton-Raphson algorithm to find the correct traction increment. This provides a more stable numerical solution and eliminates oscillatory behavior when using rate-state friction.

* Added SCEC spontaneous rupture benchmark TPV102 to the benchmark repository. PyLith produces results very similar to several other finite-element codes.

### Bug fixes

* Fixed two MPI related bugs in computing Green's functions in
parallel. The number of impulses corresponded to only those on
process 0 and the output of the impulses for vertices on processor
boundaries was inconsistent.

* Corrected computation of fault tractions (Lagrange multipliers) on
process boundaries for prescribed slip with explicit time stepping.

* Fixed bug when reading in list of output points with just one
point.

* Adjusted autoconf Python setup macro to remove temporary sysconfig.pyc file.

* Added check to make sure degree of freedom specified in Dirichlet BC is consistent with spatial dimension of problem.

* Corrected two typos in the manual related to fault opening and tractions in examples/3d/hex8/step20 and updating to the use of cell.dimension for the quadrature scheme with tetrahedral cells.

### Migrating from v1.8.x to v1.9.x

No changes are needed in .cfg files to switch from v1.8.0 to v1.9.0. Version 1.9.0 does includes some changes to the friction and material model interfaces, so extensions do require changes. See the templates for details.

## Version 1.8.0

### New features

* Additional flexibility in PETSc nonlinear solver parameters

  The default line search type for the PETSc nonlinear (SNES) solver is a customized backtrace method included in PyLith. The user may now select alternative line search types (basic, bt, l2, cp) available in PETSc.

* Post-processing utility `pylith_eqinfo` to compute slip information.

  This post-processing utility computes the moment magnitude, seismic moment, seismic potency, and average slip at user-specified snapshots in time from PyLith HDF5 output.  Information is given for each fault and across all faults. See the Post-processing section in the Running PyLith chapter of the manual for more information.

* Computation of the stable time step for explicit time-stepping.

  The stable time step for explicit time-stepping is computed based on the CFL condition and minimum edge lengths. For triangular and tetrahedral cells we also account for a reduction in the stable time step due to distorted cells (e.g., slivers and needles). See the Stable time step section in the Materials chapter of the manual for more information.

* Output the stable time step for each cell in a material.

  Output cell_info_fields "stable_dt_implicit" and "stable_dt_explicit" can be included in material output.

* Added netCDF Python module to binary distribution to provide Python interface to NetCDF files, including Exodus-II files. This is used in a new meshing example for setting the discretization size using an Exodus-II vertex field. Note that this required updating the NetCDF library.

### Bug fixes

* Fixed omission of synchronization of stable time step computation among processors. Minimum time step among all processors rather than local value should be used.

* Fixed density scale not being set in `NondimElasticQuasistatic`.  Density scale should be set based on shear modulus, length scale, and relaxation time.

* Added warning when initial state for a fault constitutive model is not set. If an initial state value is not given, for rate-state friction using a default value of L / reference slip rate. Other fault constitutive models use a default value of 0.0 for initial state variables.

* Separated tensor components in Xdmf files to avoid confusion. The corresponding HDF5 files remain unchanged.

* Removed explicit time-stepping formulation with non-lumped Jacobian. This formulation was not setup properly for spontaneous rupture models and is too computationally expensive for practical problems. The `ExplicitLumped` formulations are now simply `Explicit`.

* Fixed parallel bug that resulting in inconsistent orientation of fault slip directions. Flipping the fault orientation was not synchronized across processors. This bug would only appear when running in parallel with faults that change from dipping in one direction to dipping in the opposite direction.

* Fixed bug in setting name of field in `OutputSolnPoints` when output multiple fields. This bug caused the name of the first output field to be used and output data to overwrite each other.

### Migrating from v1.7.x to v1.8.x

Explicit time stepping with a non-lumped Jacobian has been eliminated
and ExplicitLumped is now Explicit.

```
Old setting
------------------------------------------------
formulation = pylith.problems.ExplicitLumped
formulation = pylith.problems.ExplicitLumpedTri3
formulation = pylith.problems.ExplicitLumpedTet4

New setting
------------------------------------------------
formulation = pylith.problems.Explicit
formulation = pylith.problems.ExplicitTri3
formulation = pylith.problems.ExplicitTet4
```

## Version 1.7.1

### Bug fixes

* Fixed a couple of bugs in the spontaneous earthquake rupture for quasi-static problems when running in parallel. These prevented the nonlinear solve from converging and erroneously generated fault-opening in a some cases.

* Minor updates to the documentation and manual. Added Green's function examples to the manual.


## Version 1.7.0

### New features

* User-friendly interface for Green's functions

  A new problem type provides a user-friendly interface for computing Green's functions associated with fault slip for complex spatial variation in elastic properties. See examples/2d/greensfns in the tutorials for examples.

* Output of solution field at user-specified locations

  Added a new output manager for interpolation of the solution field to user-specified point locations. This feature is useful for comparison of the solution with observations and in computing Green's functions. See examples/3d/hex8/step19 and examples/2d/greensfns in the tutorials for examples.

* Plane strain version of Drucker-Prager elastoplastic model

  Added a plane strain version of the Drucker-Prager elastoplastic model. Additionally, the user can now select whether to use an inscribed, intermediate, or circumscribed fit to the Mohr Coulomb criterion.

* Spatial and temporal variation in tractions for spontaneous earthquake rupture

  Switched from a simple constant spatial variation in initial fault tractions to the more flexible spatial and temporal variation consistent with the Dirichlet, Neumann, and point force boundary conditions. Also added a switch to turn on/off applying prescribed fault tractions when the fault opens; the default behavior is to stop applying prescribed fault tractions when the fault opens, but turning this off allows simulation of dike intrusions via prescribed fault tractions. See examples/3d/hex8/step20 in the tutorials for an example of how to specify fault tractions with the new implementation.
  
* Ability to use PETSc GPU solvers

  Added ability to build PyLith with either double (default) or single precision floating point values to facilitate use of GPUs. In order to use PETSc GPU solvers, CUDA and cusp must be installed and PETSc must be configured to use CUDA. See the PyLith manual and PETSc documentation for details.

* User-specified start time for simulations.

  Users can set the simulation start time to any desired value. This facilitates combining simulations to model the earthquake cycle.

* Elastic prestep in quasi-static simulations is optional.

  The elastic prestep in quasi-static simulations can be skipped (the default is to include the elastic prestep). This facilitates combining simulations to model the earthquake cycle.

### Bug fixes

* Fixed bug in the spontaneous earthquake rupture for quasi-static problems when running in parallel. 

### Migrating from v1.6.x to v1.7.x

Two changes are required when migrating from version 1.6 to 1.7.

1. The FIATSimplex object now has the same parameters as the `FIATLagrange` object.

```
Old setting                   New setting
------------------------      ------------------
cell.shape = line             cell.dimension = 1
cell.shape = triangle         cell.dimension = 2
cell.shape = tetrahedron      cell.dimension = 3
```

2. Prescribed fault tractions for spontaneous earthquake rupture use a new, more flexible implementation that follows the same functional form for spatial and temporal variation as that used in the Dirichlet and Neumann boundary conditions. Consequently, the output info fields are also different and follow the naming scheme used in the other time-dependent boundary conditions.

```
Old settings
--------------------------------------------------------------------
[pylithapp.timedependent.interfaces.fault]
db_initial_tractions = spatialdata.spatialdb.SimpleDB
db_initial_tractions.iohandler.filename = tractions.spatialdb
db_initial_tractions.label = Initial fault tractions

New settings
--------------------------------------------------------------------
traction_perturbation = pylith.faults.TractPerturbation
[pylithapp.timedependent.interfaces.fault.traction_perturbation]
db_initial = spatialdata.spatialdb.SimpleDB
db_initial.iohandler.filename = tractions.spatialdb
db_initial.label = Initial fault tractions
```


## Version 1.6.3

### Bug fixes

* Improved error messages for problems encountered during processing of parameters. A backtrace of the object hierarchy is now included to pinpoint in which object the error occurred.

* Added a line search to the inner friction solve in quasi-static simulations to increase the robustness of the nonlinear solve. Simulations using rate and state friction now converge under a much wider range of circumstances.

* Fixed bug in updating slip state variable in slip-weakening friction. This caused slight errors in the cumulative slip. We also added a parameter that forces healing to occur in a single time step. This is used to confine slip to a single time step in quasi-static simulations. See examples/3d/hex8/step13.cfg for an example.

* Tuned parameters in the slip-weakening friction and rate and state friction examples (step13.cfg and step14.cfg, respectively) in examples/3d/hex8 to give stick-slip behavior.

* Fixed communication issue associated with writing boundary condition information output in parallel.

* Changed info in Xdmf file for fields that are not scalars, vectors, or tensors so that the each component is extracted, facilitating visualization in ParaView. The corresponding HDF5 file remains the same.

* Added the ability to specify non-derived units (e.g., degree and radian). This is useful in specifying parameters for the Drucker-Prager elastoplastic rheologies. If no units are specified, radians are assumed.

### Internal changes

* Rate and state friction with ageing law

  The implementation of rate and state friction with ageing law was modified to work better with the iterative solver. We switched to the conventional, unregularized formulation but added a minimum cutoff for the slip rate. Below this cutoff friction has a linear rather than logarithmic dependence on slip rate. As long as this cutoff is close to the SNES solver tolerance, the difference in behavior is negligible while improving the ability of the solver to converge for very small deformations.

### Known Issues

* The rate and state friction with ageing law has not been tested for dynamic rupture simulations. We plan to run the SCEC Dynamic Rupture benchmarks for rate and state friction as soon as we add a spatial-temporal specification of initial fault tractions, which are required for the benchmark problems.

* Running simulations with more than a million cells and large faults in parallel can result in severe memory imbalances among processors. Some processors around the fault may use 10x more memory than processors away from the fault. We expect this problem to disappear in v1.7 when we switch to new, more efficient Sieve implementation.


## Version 1.6.2

### Bug fixes

* Fixed bug in writing tensor data for Xdmf files. Switched Tensor to Tensor6 to account for symmetry.

* Fixed bug in writing HDF5 files in parallel when one processor does not write any information (e.g., faults and boundary conditions).

* Added dimensioning of time dataset in HDF5 files. The units are now seconds rather than nondimensional time.

* Fixed memory allocation error (`std::bad_alloc`) when a processor did not contain cells for a boundary condition or output. This bug did not show up on all architectures.

* Increased robustness of spontaneous rupture (fault friction) implementation to broaden the range of conditions it can handle. The implementation now properly handles cases with fault opening and cases with zero shear or normal tractions.
    

### Internal changes

* Fault implementation

  Several changes have been made to the fault implementation, but none of these affect the user interface. The runtime performance is nearly identical with improved accuracy for spontaneous rupture (fault friction) simulations. These changes involved switching to using tractions (non-integrated quantities) for the Lagrange multipliers in the global coordinate system rather than integrated quantities in the fault coordinate system. Additionally, initial fault tractions are associated with the fault vertices and their interpolation uses the finite-element basis functions.

* Distribution of mesh among processors

  The data structures used to distribute the mesh among processors have been improved. This reduces memory use and runtime for this stage of the simulations.


### Known Issues

The custom line search used with the PETSc nonlinear solver (SNES)has difficulty handling some loading cases. In cases where the direction of the line search tends to be nearly orthogonal to the residual, the rate of convergence in the SNES iterations is extremely slow. In other cases the nonlinear solver gets stuck in a local minimum. We plan to improve the line search algorithm in a future release in order to resolve this issue and improve the rate of convergence in spontaneous rupture simulations.


## Version 1.6.1

* Validation of user input

  Added stricter requirements for descriptive labels of various objects, including spatial databases and friction models. The default labels are empty strings which do not result in useful error messages; the user is now required to specify a non-empty string for the labels. This makes errors related to spatial databases much easier to diagnose.


* Updates to manual

  * Updated description of cell_info_fields for Neumann boundary condition. The description had not been updated to reflect the time-dependence introduced in version 1.4.

  * Added steps 18 and 19 that discuss time-dependent Neumann boundary conditions to examples/3d/hex8.


### Bug fixes

* Fixed bug in writing rupture information to VTK and HDF5 files when using multiple earthquake sources. Field names did not include name of rupture. This caused loss of information in VTK output and a corrupted Xdmf metadata file for HDF5 output.

* Fixed error in use of initial stress tensor with generalized Maxwell models. The initial stress tensor was added to the current stress tensor twice.

* Fixed two bugs in the fault friction implementation. One bug pertained to accounting for roundoff errors and convergence tolerances in computing the slip rate. Slip rates less than 1.0e-12 (nondimensionalized) are set to zero. The friction implementation for quasi-static problems contained a bug that resulted in slip extending over all of the fault rather than the appropriate isolated patch.

* Cleaned up Green's function example (examples/greensfns/hex8) so that it runs without errors. Eliminated extraneous processing.

* Cleaned up meshing examples (examples/meshing), including elimination of superfluous pre-processing.

* Adjusted absolute tolerances for PETSc solves in examples/3d/hex8 so that solver terminates with desired convergence criterion.

* Updated examples/2d/subduction/geometry.jou to use APREPRO functions and variables to store id values.


## Version 1.6.0

### New features

* Parallel binary output via HDF5

  Provides much faster output by writing HDF5 files in parallel, which can be accessed directly from Matlab or indirectly from ParaView or Visit via automatically created Xdmf files. Temporal data is stored in 3-D arrays, permitting slicing in time and/or space. See examples/3d/hex8 Steps 6-9 and examples/2d/subduction in the tutorials for examples.

* 2-D generalized Maxwell viscoelastic bulk rheology

  Added a 2-D generalized Maxwell viscoelastic bulk rheology corresponding to the plane strain version of the 3-D generalized Maxwell viscoelastic model.

* Time-weakening fault constitutive model

  Added a linear time-weakening fault constitutive model. Some spontaneous rupture modelers prefer this model over linear slip-weakening because it is easier to maintain resolution of the cohesive zone.

* Global uniform parallel mesh refinement

  Permits running larger problems through uniform global refinement of the mesh by a factor of 2 (reduces the node spacing by a factor of 2) after the mesh is distributed among processors. This allows running problems that are 4x larger in 2-D and 8x larger in 3-D. See examples/3d/tet4 Steps 2 and 4 for examples.
  
* Custom algebraic multigrid preconditioner

  Adds a custom preconditioner for Lagrange multiplier degrees of freedom associated with fault slip via prescribed slip or spontaneous ruptures with algebraic multigrid preconditioning for quasi-static solutions. In most cases, this results in fewer iterations in the linear solve and the number of iterations increases much less with problem size. See examples/3d/tet4 Steps 2 and 4 for examples.

* PyLith installer utility

  This utility provides a much more robust method for building PyLith and all of its dependencies from source, including dependency checking, installation to a central location, and creation of a shell script to set environment variables.

### Bug fixes

* Fixed the fault friction implementation to correctly update Lagrange multiplier values when the slip is overestimated in an iteration. This primary fixes problems encountered with the use of the Dieterich-Ruina rate and state fault constitutive model.

* Corrected viscoelastic rheologies to properly account for a nonzero initial strain tensor.


### Migrating from v1.5.x to v1.6.x

No changes in parameters are required. Version 1.6.1 does require users to specify descriptive labels for spatial databases and friction models.


## Version 1.5.2

* PyLith 1.5.2 requires FIAT version 0.9.9 or later and an updated PETSc development version. It also requires users to update to the latest spatialdata version for compatibility of the SWIG generated files. These are included in the binary distribution, but users building PyLith from source will need to update FIAT, PETSc, and spatialdata.

### Bug fixes

* Fixed setting of elastic constants in DruckerPrager3D and computation of the yield function. Some off-diagonal elasticity constants were off by a factor of 2.0 and the yield function was missing a factor of 0.5 and sqrt().

* Fixed computation of stable time step when using initial stresses with PowerLaw3D. If effective stress is zero, then stable time step is infinite.

* Re-enabled check for compatibility of quadrature scheme and cells for bulk rheologies.

* Added check to configure for compatible version of FIAT.

* Fixed bug where buffer for output of initial stresses for dynamic (spontaneous) rupture.


## Version 1.5.1

### Bug fixes

* Fixed dimensioning of velocity and acceleration fields in output. The scales were set to just the length scale rather than the length scale divided by the time scale and length scale divided by the time scale squared.

* Fixed partitioning of cohesive cells. Cohesive cells were ignored during partitioning of the mesh, so they were randomly distributed among processors.


## Version 1.5.0

* Fault constitutive models

  Added fault friction interface conditions with static friction, linear slip-weakening friction, and rate- and state-friction with the ageing law. The implementation can be used in static, quasi-static, and dynamic problems.

* Drucker-Prager elastoplastic bulk rheology

  Added a Drucker-Prager elastoplastic bulk rheology. This is a perfect plasticity implementation (no hardening). This is a nonlinear constitutive model, so the nonlinear solver is required when this rheology is used. Refer to the 'Material Models' section of the manual.

* Plane strain Maxwell viscoelastic bulk rheology

  Linear Maxwell viscoelastic rheology for plane strain problems.

* Finite-deformation formulation

  Added a finite-deformation (rigid body motion and small strains) implementation of elasticity with stress calculated using the Second Piola Kirchhoff stress tensor and strains calculated using the Green-Lagrange strain tensor.

* Lumped Jacobian for explicit-time stepping

  Added the option to lump cell Jacobian matrices to form a diagonal system Jacobian matrix for explicit time stepping. This decouples all degrees of freedom and permits use of a fast, trivial, direct solver.

* Optimized elasticity objects

  Added optimized elasticity objects for the most popular cell types and basis functions (linear polynomials). For tri3 and tet4 cells with one quadrature point, the optimized implementations do not use reference (mapped) cells in order to reduce the number of operations.

* Scientific notation for ASCII VTK files

  Data values in ASCII data files are written in scientific notation with user-specified precision.

* Nodeset names in CUBIT Exodus files

  Use of nodeset names in CUBIT Exodus files for boundary conditions and faults. Users can specify to use nodeset names (default behavior) or ids.

* Velocity and slip rate as output fields

  Velocity (domain and subdomain) and slip rate (fault) fields are can be requested as output fields. The fields are computed using the time-stepping algorithm and alleviates the need to compute them via post-processing.

* Dimensionless values in spatial databases no longer need artificial dimensions. Values without dimensions are understood by the parser as dimensionless quantities.

* Bug fixes

* Updating state variables did not retrieve physical properties for cell. Last physical properties retrieved were used. Physical properties are now retrieved when updating state variables.

* Fixed incorrect dimensioning of physical properties and state variables for the power-law rheology in output.

* Fixed memory bug for a fault in a 1-D mesh when constructing the cohesive cells.

### Contributors

* Surendra Somala - fault friction implementation.

### Migrating from v1.4.x to 1.5.x

Three changes to the code require updating old parameters settings for use with version 1.5.

1. Recent releases of CUBIT include nodeset names in the Exodus file and PyLith now uses them to associate vertices with boundary conditions and faults. Use the NetCDF utility `ncdump` to examine the contents of the Exodus (.exo) file to see it it includes the variable ns_names. If it does, then use nodeset names rather than nodeset ids for boundary condition label properties. If your Exodus file does not contain nodeset names, then set the MeshIOCubit property use_nodeset_names to False to continue to use nodeset id values for boundary condition labels.

2. The power-law constitutive parameters have been changed so that the parameter units are no longer dependent on the power-law exponent. This is a more logical implementation and allows (among other things) users to vary power-law parameters using a spatial database. Previously, it was not possible to vary power-law parameters unless everything used the same power-law exponent. The new implementation uses reference-strain-rate, reference-stress, and power-law-exponent to describe the material. This is described in the 'Material Models' section of the manual.

3. The fault property `normal_dir` is obsolete. Only the property `up_dir` is required to enforce that positive slip is left-lateral, reverse, and fault-opening for dipping faults in 2-D and horizontal fault surfaces in 3-D. Previously, in 2-D positive slip was always left-lateral, but now the up-direction is used to enforce positive slip corresponds to reverse motion for dipping faults. For horizontal fault surfaces in 3-D a normal of (0,0,1) is assumed in determining the up-dip direction.

### Migrating from v1.3.x to 1.4.x

A number of changes to the code require updating old parameter settings for use with version 1.4.

1. The mesh "importer" is now called "reader".

2. The spatial database facility for a material, `db`, is separated into a `db_properties` and a `db_initial_state`. The initial stress and strain tensors are specified using the `db_initial_stress` and `db_initial_strain` facilities. The names of some of the spatial database values for physical properties for viscoelastic properties have changed.

3. The code is now intelligent enough to determine the dimensions of the quadrature required (e.g., Quadrature2D and Quadrature2Din3D, etc). Setting the quadrature to the object for a given spatial dimension and cell dimension is no longer allowed because it is done automatically.

4. The names of the output filters have changed and include suffixes Mesh or SubMesh to indicate that they operate on a mesh or submesh (e.g., CellFilterAvg is now CellFilterAvgMesh or CellFilterAvgSubMesh). This is related to the use of C++ templates.

5. The DirichletPoints boundary condition has been renamed to DirichletBC.

6. The procedure for enabling certain features no longer involves setting a "use" property to True. Instead, the features are enabled when the user sets the component to a facility. This applies to gravity, initial stresses, initial strains, and initial state variables, and time-dependent boundary conditions (Dirichlet, Neumann, and point force).

7. Nondimensionalization of the problem eliminates the need to condition the fault constraints. The "mat_db" facility was removed.

8. The Dirichlet and Neumann boundary conditions now follow a more general time dependence. The names of the facilities and the names of the values in the spatial databases are, in most cases, different.

9. The `FixedDOFDB` has been renamed to `ZeroDispDB` in order to better reflect the type of spatial database.


## Version 1.3.1

* Added stages to PETSc logging (`--petsc.log_summary`) to collect event logging into groups.

### Bug fixes

* Fixed partitioning options. Partitioning options were ignored in the 1.3.0 release.

* Fixed assembling of Jacobian, residual, and fault sections across processors. This bug caused errors in the computation of the change in tractions over the fault surface.


## Version 1.3.0

* New time stepping options

  In addition to a uniform, user-specified time step, which is the default, there are two new time-stepping options. The user may supply a file with nonuniform time steps or, for quasi-static simulations, the user can request the code to compute the time step automatically. For the current bulk constitutive models, the automatically determined time step is independent of the deformation rate, so it is uniform.

* Initial stresses

  Users may optionally supply an initial stress state for each material via a spatial database. The initial stress state can balance the gravitational body forces so that the model is in equilibrium without any deformation. This implementation of an initial stress state is a prelude to specifying an initial state for each material, which will be available in a future release.

### Bug fixes

* Fixed labeling of physical properties in output for the Maxwell viscoelastic and generalized Maxwell viscoelastic materials (mu and lambda were switched).

### Migrating from v1.2.x to v1.3.x

The implementation of different options for controlling the time step requires adjusting input parameters from those used with PyLith 1.2. The time stepping is specified under the time-stepping formulation rather than the problem (i.e., one level deeper).


## Version 1.2.0

* New Sieve implementation

  The previous implementation of Sieve provided a very generalized implementation of data structures and operations for finite-element meshes. Switching to a more rigid implementation in the new implementation streamlined the data structures, resulting in a significant reduction in the memory use for storing the mesh. This leads to an overall reduction in memory use of 25-30% in many cases.

* Multiple kinematic ruptures

  A single kinematic rupture on a fault has been replaced by a dynamic array of kinematic ruptures. This allows creation of an arbitrary number of kinematic ruptures on each fault surface. By using spatial databases to control the spatial and temporal extent of slip in each rupture independently, slip from different earthquake ruptures can overlap in space and/or time. Additionally, the rupture time at each location is specified with respect to the origin time of the corresponding earthquake rupture.

* New slip time functions

  * Step slip time function (now the default)

    This slip time function simplifies specifying instantaneous slip in a quasi-static simulation compared with using the Brune slip time function.

  * Constant slip rate slip time function

    This slip time function permits prescribing a constant slip rate on the fault surface.

* Gravitational body forces

  Gravitational body forces are implemented (they are turned off by default). The direction and acceleration of gravity may be specified.

* Fixed Makefile.am files to not delete source files during "make clean" when building in the source tree.


### Migrating from v1.1.x to 1.2.x

There are two new features in PyLith version 1.2 that require users to adjust input parameters from those used with PyLith 1.1. A dynamic array of kinematic rupture replaces a single kinematic rupture on a fault. Additionally, the default slip time function is now a step-function. This eliminates the need to specify a peak slip rate for quasi-static simulations. When using PyLith version 1.2 with a problem previously setup for PyLith 1.1, look for warnings about unknown components and settings in the screen output at the beginning of a run.


## Version 1.1.2

* Fixed bug in output of solution over sub-domain boundary surfaces in parallel.

* Fixed Makefile.am files to include documentation files in source distribution.


## Version 1.1.1

* Fixed Makefile.am files to include files missing from the source distribution.


## Version 1.1.0

* New boundary conditions

  * Neumann (traction) boundary conditions

  * Absorbing boundary conditions via simple, tuned dampers

  * Dirichlet boundary conditions with displacement and/or velocity values

* New bulk constitutive models

  * Generalized Maxwell viscoelastic model 

* New output implementation

  The output to VTK files has been completely rewritten. This new implementation includes output of physical properties and state variables associated with the bulk constitutive models, as well as output of fault information (earthquake rupture parameters and slip and traction time histories). Additionally, the VTK file with the solution no longer includes fault related values- it contains just the displacement field over the domain as one would expect. A user can now also request output of the solution over an arbitrary number of sub-domains of the domain boundary, e.g., the ground surface. For each of these different kinds of output, the frequency of output and the values included can be customized by the user. The names of the VTK files and the variable names have also been adjusted to permit animation of solutions within most VTK visualization tools.

* New spatial database implementations

  Spatialdata includes two new spatial database implementations. The `SCECCVMHDB` provides a seamless interface to the SCEC CVM-H seismic velocity model for elastic material properties. The UniformDB permits creating a spatial database for uniform values using only .cfg files or the command line; this eliminates the need to create a SimpleDB database file with one location.

* Dynamic arrays of components in Pyre

  Pyre now contains dynamic arrays of components, eliminating the need for containers for materials, boundary conditions, and faults.

* Better consistency checking of input parameters

  * Uniqueness of material identifiers for materials and faults is enforced.

  * The material identifier of each cell in the mesh is checked to make sure it matches a material model.

  * Each boundary, interface condition, and output group is checked to make sure it exists in the mesh.

### Bug fixes

* Fixed bug causing segmentation fault with multiple, non-overlapping Dirichlet boundary conditions applied to vertices.

* Fixed numerous bugs related to explicit time integration for dynamic problems.

* Eliminated several small memory errors.

* Fixed several bugs associated with writing VTK files in parallel.

### Known issues

* PyLith still uses much more memory that PyLith 0.8 due to the current general Sieve implementation. A much more efficient, albeit less general Sieve implementation is under development. Additionally, distribution of the mesh will also be improved in a future release.

* The preconditioner for explicit time stepping provides relatively poor overall performance compared to a direct solve with traditional mass lumping. An appropriate preconditioner and traditional mass lumping will be supported in a future release.



### Migrating from v1.0.x to 1.1.x

There are two new features in PyLith version 1.1 that require users to adjust input parameters from those used with PyLith 1.0. The elimination of containers in favor of the dynamic arrays of components present in the latest version of Pyre requires switching from setting the container to specifying the array of components on the command line or .cfg file. Additionally, the new implementation of output requires a completely new set of parameters. When using PyLith version 1.1 with a problem previously setup for PyLith 1.0, look for warnings about unknown components and settings in the output at the beginning of a run.


## Version 1.0.2

* Performance optimizations have significantly reduced runtime and memory use relative to version 1.0.1. The default quadrature order for tetrahedral cells is now 1, which is appropriate for the default basis functions.

* Added checks to verify the compatibility of quadrature scheme for solid and cohesive cells.

### Bug fixes

* In some cases, cohesive cells were not inserted into the finite-element mesh properly. The cells mixed together vertices from the different sides of the fault. A more efficient procedure for creating cohesive cells fixed this problem.

* Cell adjacency graph was created incorrectly which resulted in a poor quality of partitioning among processors.

* VTK output for meshes with N faults included cohesive cells for N-1 faults. Since VTK output does not understand cohesive cells, we now remove all cohesive cells from the VTK output.

* Using the SimpleDB in Spatialdata from Python limited interpolation to the "linear" scheme instead of allowing use of the "nearest" scheme. Setting the SimpleDB property to "nearest" and "linear" now works as expected.

* The reader for Spatialdata coordinate systems information did not correctly putback characters in the input stream, resulting in reading errors. The putback routines were fixed.

* Fault "up" and "normal" directions remained as string arrays when passed to the module, instead of being converted to float arrays.


## Version 1.0.1

### Bug fixes

* Cohesive cells lacked consistent orientation (inconsistent normals) in cases where cells were not ordered one side of the fault and then the other.

* Final slip of zero resulted in fault slip and slip increments of Nan.

* Parallel importing of meshes from LaGrit and CUBIT lacked guards against all processors reading the files.


## Version 1.0.0

* Code now includes both dynamic and quasi-static solutions.

* Completely rewritten in Python and C++, with bindings provided by Pyrex/Pyrexembed.

* Easier specification of simulations:

  * Parameters are all set using .cfg files (or .pml/command-line).

  * Mesh may be directly imported from CUBIT, LaGriT, or using PyLith mesh ASCII format.

  * Material properties, fault dislocations, and BC are all given using spatial databases, which are independent of mesh discretization.

* Faults are now implemented using cohesive elements:

  * Easy specification of kinematic fault slip using a spatial database.

  * Cohesive elements generate offsets in the mesh corresponding to fault slip, which increase the accuracy of displacement fields near faults and facilitate visualization of fault slip.

  * Usage of cohesive elements will facilitate the upcoming addition of fault constitutive relations, where fault slip occurs in response to prescribed physics.

* Improved implicit time-stepping eliminates need to perform more than one iteration for linear rheologies.

* Code is now completely modular and object-oriented, which allows much easier addition of new features.  Modules may be added without having to recompile the code.

* Features present in 0.8 that are not present in 1.0 that will be added in the near future.

  * Traction boundary conditions

  * Generalized Maxwell and Power-law Maxwell viscoelastic models
