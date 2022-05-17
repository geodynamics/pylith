# Cubit Mesh

## Geometry

We construct the geometry by taking a horizontal cross-section of a 3D block.
Alternatively, we could have constructed the geometry by building it up from points and curves like we did with Gmsh.
{numref}`fig:example:strikeslip:2d:geometry:cubit` shows the geometry and variables names of the vertices and curves.

:::{figure-md} fig:example:strikeslip:2d:geometry:cubit
<img src="figs/geometry-cubit.*" alt="Geometry created in Cubit for generating the mesh." scale="75%"/>

Geometry created in Cubit for generating the finite-element mesh.
The names of the verties and curves match the ones we use in the Cubit journal files.
:::

## Meshing using Journal Scripts

We use Cubit journal files `mesh_tri.jou`  and `mesh_quad.jou` to generate triangular and quadrilateral meshes, respectively.
Both of these journal files make use of the `geometry.jou`, `gradient.jou`, and `createbc.jou` files for creating the geometry, setting the discretization size, and tagging boundary conditions, faults, and materials, respectively.
We use the Cubit graphical user interface to play the Journal files.

:::{note}
Examine how we set the discretization size in Gmsh and Cubit.
In both cases the discretization size increases at a geometric rate with distance from the fault.
For this simple geometry, it required less than 10 lines of Python code in Gmsh and about the same number of lines in Cubit.
In Gmsh the code is very general and remains the same even as the domain geometry becomes more complex, whereas in Cubit the number of commands increases with the complexity of the geometry.
:::
