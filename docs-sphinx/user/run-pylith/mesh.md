# Finite-Element Mesh

Geometrical and topological information for the finite element mesh may be provided by exporting an Exodus II format file from CUBIT, by exporting a GMV file and an accompanying Pset file from LaGriT, or by specifying the information in PyLith mesh ASCII format.

PyLith supports triangular and quadrilateral cells in 2D and tetrahedral and hexahedral cells in 3D.
The vertex ordering must follow the convention shown in {numref}`fig:2d:cells` and {numref}`fig:3d:cells`.
The cells define the geometry of the domain; the basis order and quadrature order used to discretize the solution subfields are specified separately.

The mesh information specifies the vertex coordinates and the vertices composing each cell in the mesh.
The mesh information must also define at least one set of vertices for which displacement (Dirichlet) boundary conditions will be provided.
In most realistic problems, there will be several vertex groups, each with a unique identifying label.
For example, one group might define a surface of the mesh where displacement (Dirichlet) boundary conditions will be applied, another might define a surface where traction (Neumann) boundary conditions will be applied, while a third might specify a surface that defines a fault.
Similarly, the mesh information contains cell labels that define the material type for each cell in the mesh.
For a mesh with a single material type, there will only be a single label for every cell in the mesh.
See {ref}`sec-user-physics-materials` and {ref}`src-user-physics-boundary-conditions` for more detailed discussions of setting the materials and boundary conditions.

:::{figure-md} fig:2d:cells
<img src="figs/cells2d.*" alt="2D cell types" width="400px">

Cells available for 2D problems are the triangle and the quadrilateral.
:::

:::{figure-md} fig:3d:cells
<img src="figs/cells3d.*" alt="3D cell types" width="400px">

Cells available for 3D problems are the tetrahedron and the hexahedron.
:::

## Mesh Importer

The default component for the PyLithApp `mesher` facility is `MeshImporter`, which provides the capabilities of reading the mesh from files.
The `MeshImporter` includes a facility for reordering the mesh.
Reordering the mesh so that vertices and cells connected topologically also reside close together in memory improves overall performance.

:::{admonition} User Interface
:class: seealso
See [`MeshImporter` component](../components/topology/MeshImporter.md)
:::

## ASCII Mesh Files - `MeshIOAscii`

The `MeshIOAscii` object is intended for reading small, simple ASCII files containing a mesh constructed by hand.
We use this file format extensively in small tests.
{ref}`src-user-formats-MeshIOAscii` describes the format of the files.

:::{seealso}
[`MeshIOAscii` Component](../components/meshio/MeshIOAscii.md)
:::

(sec-user-run-pylith-meshiocubit)=
## CUBIT (Exodus II) Mesh Files - `MeshIOCubit`

The `MeshIOCubit` object reads the NetCDF Exodus II files output from CUBIT/Trelis.
Beginning with CUBIT 11.0, the names of the nodesets are included in the Exodus II files and PyLith can use these nodeset names or revert to using the nodeset ids.

:::{seealso}
[`MeshIOCubit` Component](../components/meshio/MeshIOCubit.md)
:::

## LaGriT Mesh Files - `MeshIOLagrit`

The `MeshIOLagrit` object is used to read ASCII and binary GMV and PSET files output from [LaGriT](https://lanl.github.io/LaGriT/`).
PyLith will automatically detect whether the files are ASCII or binary.
We attempt to provide support for experimental 64-bit versions of LaGriT via flags indicating whether the FORTRAN code is using 32-bit or 64-bit integers.

:::{danger}
The PyLith developers have not used LaGriT since around 2008 and there have been a few releases since then so the interface may not be up to date.
:::

:::{seealso}
[`MeshIOLagrit` Component](../components/meshio/MeshIOLagrit.md)
:::

## Distribution among Processes - `Distributor`

The distributor uses a partitioner to compute which cells should be placed on each processor, computes the overlap among the processors, and then distributes the mesh among the processors.
The type of partitioner is set via PETSc settings.

:::{note}
METIS/ParMETIS are not included in the PyLith binaries due to licensing issues.
:::

:::{seealso}
[`Distributor` Component](../components/topology/Distributor.md)
:::

## Uniform Global Refinement - `Refiner`

The refiner is used to decrease node spacing by a power of two by recursively subdividing each cell by a factor of two.
In a 2D triangular mesh a node is inserted at the midpoint of each edge, splitting each cell into four cells (see {numref}`fig:uniform:refinement:2x`).
In a 2D quadrilateral mesh a node is inserted at the midpoint of each edge and at the centroid of the cell, splitting each cell into four cells.
In a 3D tetrahedral mesh a node is inserted at the midpoint of each edge, splitting each cell into eight cells.
In a 3D hexahedral mesh a node is inserted at the midpoint of each edge, the centroid of each face, and at the centroid of the cell, splitting each cell into eight cells.

:::{figure-md} fig:uniform:refinement:2x
<img src="figs/refinement2x.*" alt="Global refinement" width="100%"/>

Global uniform mesh refinement of 2D and 3D linear cells.
The blue lines and orange circles identify the edges and vertices in the original cells.
The purple lines and green circles identify the new edges and vertices added to the original cells to refine the mesh by a factor of two.
:::

Refinement occurs after distribution of the mesh among processors.
This allows one to run much larger simulations by (1) permitting the mesh generator to construct a mesh with a node spacing larger than that needed in the simulation and (2) operations performed in serial during the simulation setup phase, such as, adjusting the topology to insert cohesive cells and distribution of the mesh among processors uses this much smaller coarse mesh.
For 2D problems the global mesh refinement increases the maximum problem size by a factor of $4^{n}$, and for 3D problems it increases the maximum problem size by a factor of $8^{n}$, where $n$ is the number of recursive refinement levels.
For a tetrahedral mesh, the element quality decreases with refinement so $n$ should be limited to 1-2.

% End of file
