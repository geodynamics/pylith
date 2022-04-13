# Mesh Description

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
