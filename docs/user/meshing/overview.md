# Overview

Cubit and Gmsh can generate quadrilateral or triangular meshes in 2D and hexahedral or tetrahedral meshes in 3D.
In each case, you first create the geometry, specify the meshing algorithm and discretization size, and then generate the mesh.
You can build up the geometry from points, curves, surfaces, and volumes or use the geometry engines to construct the domain using simple shapes.

## 2D meshing

Constructing surfaces from points and curves for 2D meshing with Cubit and Gmsh is very similar.
Cubit provides more geometric operations than Gmsh, but many simple geometric operations in Gmsh can be implemented by the user when using the Python interface.
Generating unstructured quadrilateral meshes for complex geometry is often easier in Cubit, whereas generating meshes with complex specification of discretization size is often easier with Gmsh.

## 3D meshing

Cubit provides an extensive suite of tools for constructing complex 3D geometry.
This includes building surfaces and performing geometric operations on surfaces and volumes.
The suite of tools in the Gmsh geometry engine is more limited; the OpenCASCADE engine interface provides additional tools.
With either Cubit or Gmsh, you can use external CAD tools to generate the geometry.
As in the case with 2D meshing, generating unstructured hexahedral meshes is often easier in Cubit, whereas generating meshes with complex specification of discretization size is often easier in Gmsh.

## Marking surfaces

*New in v5.0.0.*

Beginning with v5.0.0, PyLith supports marking surfaces associated with boundary conditions, fault interfaces, and other surfaces using the cell faces rather than vertices.
In Cubit, marking faces rather than vertices corresponds to using sidesets rather than nodesets.
In Gmsh, marking faces rather than vertices means a physical group needs to contain only entities at the highest dimension (surfaces in 3D and curves in 2D).
Marking surfaces using vertices will be removed in v6.0.0.

## Reading meshes in parallel

*New in v5.0.0.*

When running PyLith with multiple processes (parallel processing), Exodus II files from Cubit and Gmsh files are read in serial and then distributed among the processes.
This means that one process must be able to hold the entire mesh in memory, and this approach does not scale to very large simulations.
The preferred approach for large simulations is to convert the finite-element mesh to an HDF5 file in the PETSc mesh format.
This format is designed for parallel I/O, and PyLith will read the mesh in parallel.
Refer to {ref}`sec-user-run-pylith-convertmesh` for more information about converting mesh files.
