# Examples: 3-D Subduction Zone

This suite of examples demonstrates use of a wide variety of features
and the general workflow often used in research simulations.  The
examples consistent of a step-by-step sequence of eight problems
involving a 3D subduction zone. They focus on modeling the
deformation associated with the the subducting slab, including
interseismic deformation with aseismic slip (creep) and viscoelastic
relaxation, coseismic slip on the slab interface and a splay fault,
and slow slip events on the subduction interface. We want to account
for the 3D material properties associated with different elastic
properties for the subducting slab, mantle, continental crust, and an
accretionary wedge. All of the examples use the same mesh, generated
using Cubit.

Directory structure:

* **input**: Large input files (Cubit and Gmsh mesh files) not included in the PyLith repository.
* **output**: Simulation output files, created automatically by PyLith.
* **viz**: Python scripts used in visualization.

Examples:

* `Step 1`: Static axial compression
* `Step 2`: Quasistatic coseismic and postseismic deformation from an earthquake rupture in the center of the subduction zone interface
* `Step 3`: Quasi-static interseisic deformation with creep on the top and bottom of the slab
* `Step 4`: Quasi-static earthquake cycle with prescribed earthquake rupture and creep
* `Step 5`: Spontaneous rupture driven by subducting slab (not yet updated for v3.0)
* `Step 6`: Prescribed slow-slip event
* `Step 7a,b`: Inversion of slow-slip event using 3D Green's functions
* `Step 8a,b,c`: Stress field due to gravitational body forces

Refer to the examples section of the PyLith documentation for more information.
