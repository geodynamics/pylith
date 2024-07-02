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
As in the case with 2D meshing, generating unstructure hexahedral meshes is often easier in Cubit, whereas generating meshes with complex specification of discretization size is often eaiser in Gmsh.
