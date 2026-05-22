(sec-user-femesh)=
# Finite-Element Mesh

The finite-element mesh specifies the geometry and topology of the discretization of the domain.
It must be generated using external software before running PyLith.
PyLith supports triangular and quadrilateral cells in 2D and tetrahedral and hexahedral cells in 3D.
PyLith does not support meshes with mixtures of cell shapes (meshes with triangular and quadrilateral cells in 2D or meshes with tetrahedral, hexahedal and wedges in 3D).
The vertex ordering must follow the convention shown in {numref}`fig:2d:cells` and {numref}`fig:3d:cells`.
The cells define the geometry of the domain; the basis order and quadrature order used to discretize the solution subfields are specified separately.

The mesh information specifies the vertex coordinates and the vertices composing each cell in the mesh.
The mesh information must also define at least one boundary for which displacement (Dirichlet) boundary conditions will be provided.
In most simulations, there will be several groups of faces (or vertices) for boundary conditions, each with a unique identifying label.
For example, one group might define a surface of the mesh where displacement (Dirichlet) boundary conditions will be applied, another might define a surface where traction (Neumann) boundary conditions will be applied, while a third might specify a surface that defines a fault.
Similarly, the mesh information contains cell labels that define the material type for each cell in the mesh.
For a mesh with a single material type, there will only be a single label for every cell in the mesh.
See {ref}`sec-user-physics-materials` and {ref}`sec-user-physics-boundary-conditions` for more detailed discussions of setting the materials and boundary conditions.

:::{figure-md} fig:2d:cells
<img src="figs/cells2d.*" alt="2D cell types" width="400px">

Cells available for 2D problems are the triangle and the quadrilateral.
:::

:::{figure-md} fig:3d:cells
<img src="figs/cells3d.*" alt="3D cell types" width="400px">

Cells available for 3D problems are the tetrahedron and the hexahedron.
:::

:::{toctree}
initializer.md
phases.md
:::
