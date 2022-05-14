# Mesh Description

For this small, simple problem we constructed the finite-element mesh manually in a text editor using the [`MeshIOAscii` file format](../../file-formats/meshio-ascii.md).
{numref}`fig:example:box:2d:mesh` shows our mesh with quadrilateral cells and a discretization size of 4 km from the file `quad.mesh`; we also provide a mesh with triangular cells in `tri.mesh`.
The mesh file contains the coordinates of the vertices, the index of the vertices in each cell, the material identifier for each cell, and groups of vertices associated with boundaries.
The cells are defined by its vertices listed in a counterclockwise traversal of the cell boundary; the starting point is arbitrary.
For example, cell 0 is defined by the vertices (0, 5, 6, 1);
we also could have defined the cell using (5, 6, 1, 0), (6, 1, 0, 5), or (1, 0, 5, 6).
We use the groups of vertices to apply boundary conditions.

:::{figure-md} fig:example:box:2d:mesh
<img src="figs/mesh.*" alt="Uniform resolution finite-element mesh with quadrilateral." scale="75%"/>

Uniform resolution finite-element mesh with quadrilateral cells.
Vertices are numbered sequentially from zero with the labels adjacent to each vertex.
Cells are numbered sequentially from zero with the labels at the center of each cell.
:::
