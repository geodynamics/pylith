(cha:glossary)=
# Glossary

## Pyre

component
: Basic building block of a Pyre application.
A component may be built-up from smaller building blocks, where simple data types are called properties and data structures and objects are called facilities.
In general a component is a specific implementation of the functionality of a facility.
For example, SimpleDB is a specific implementation of the spatial database facility.
A component is generally composed of a Python object and a C++ object, although either one may be missing.
We nearly always use the naming convention such that for an object called Foo the Python object is in file Foo.py, the C++ class definition is in Foo.hh, inline C++ functions are in foo.icc, the C++ class implementation is in Foo.cc, and the SWIG interface file that glues the C++ and Python code together is in Foo.i.

facility
: Complex data type (object or data structure) building block of a component.
See "component".

property
: Simple data type (string, integer, real number, or boolean value) parameter for a component.

## DMPlex

The plex construction is a representation of the topology of the finite-element mesh based upon a covering relation.
For example, segments are covered by their endpoints, faces by their bounding
edges, etc.
Geometry is absent from the plex, and is represented instead by a field with the coordinates of the vertices.
Meshes can also be understood as directed acyclic graphs, where we call the constituents points and arrows.
See {ref}`sec-developer-petsc-fem` for additional details.


cell
: The highest dimensional elements of a mesh, or mesh entities of codimension zero.

cone
: The set of entities which cover any entity in a mesh.
For example, the cone of a triangle is its three edges.

dimension
: The topological dimension of the mesh, meaning the cell dimension.
It can also mean the dimension of the space in which the mesh is embedded, but this is properly called the embedding dimension.

dual space
: For any vector space $V$ over a field $F$, the dual space $V^{*}$ is the set of all linear functions $\phi : V \rightarrow F$.
If $V$ is finite dimensional, then $V^{*}$ has the same dimension as $V$.
See <https://en.wikipedia.org/wiki/Dual_space> and <https://https://finite-element.github.io/L2_fespaces.html>.

face
: Mesh elements that separate cells, or mesh entities of codimension one.

field
: A parallel section which can be completed, or made consistent, across process boundaries.
These are used to represent continuum fields.

mesh
: A finite element mesh, used to partition space and provide support for the basis functions.

projection
: Interpolation of analytical or discretized functions into the finite-element space.
See {ref}`sec-developer-petsc-projection` for more information.

section
: These objects associate values in vectors to points (vertices, edges, faces, and cells) in a mesh.
The section describes the offset into the vector along with the number of values associated with each point.

support
: The set of mesh entities which are covered by any entity in a mesh.
For example, the support of a triangle is the two tetrahedra it separates.

vertex
: The zero dimensional mesh elements.

## PyLith

auxiliary field
: A field used to specify parameters and state variables for physics.
For an elastic materials the auxiliary field contains the elastic properties.
For a Dirichlet boundary condition the auxiliary field contains the parameters specifying the solution as a function of space and time.

basis order
: Order of the polynomial basis functions used to represent a field.

cohesive cell
: A zero volume cell inserted between two cells which share a fault face.
They are prisms with a fault face as the base.

derived field
: A field that can be computed from the solution field.
For example, the infinitesimal strain can be computed from the displacement field.

quadrature order
: Order of the quadrature scheme used in numerical integration.
The integrals in the weak form are integrated using numerical quadrature.
The location of the quadrature points are usually found to achieve a desired order of accuracy in integrated fields.

solution field
: A field with subfields for the unkown values that will be determined by solving a system of equations.

spatial database
: Data specifying values of one or more fields as a function of space.
The topology and resolution of the spatial database is set to match the data and is independent of the finite-element mesh.
The `spatialdata` library contains a few different implementations of spatial databases.
