(sec-user-meshing-gmsh)=
# Gmsh

We only recently started using Gmsh and have only used the [Python interface](https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/api/gmsh.py).
Gmsh also offers a simple scripting language, similar to Cubit journal files.
The Gmsh Python interface integrates well with the rest of Python; it can be installed so that it is compatible with the Python interpreter used by PyLith.
This means one can leverage additional Python packages, such as geographic projection libraries. 
Gmsh includes its own geometry engine as well as an interface to the OpenCASCADE engine.

:::{warning}
Gmsh does not construct quadrilateral or hexahedral meshes directly; instead, it first constructs a triangular or tetrahedral mesh and then combines triangles or tetrahedra to form quadrilaterals or hexahedra.

Sometimes it will not be able to remove all triangular or hexahedral cells, resulting in meshes with multiple shapes, which PyLith does not support.
:::

## Gmsh tips

The function `gmsh.model.geo.synchronize()` copies entities from the geometry engine (built in or OpenCASCADE) into the model.
You need to call this function after you create geometry and _before_ you attempt to use it in the model.

### Built-in CAD engine

This is a list of some of the more useful functions in the built-in CAD engine.
Some functions include optional parameters that are not shown; refer to the Python API documentation for complete specification of the API.

`gmsh.model.geo.add_point(x: float, y: float, z: float) -> int` 
: Create a point at (x, y, z) and return the tag of the point.

`gmsh.model.geo.add_line(startTag: int, endTag: int) -> int` 
: Create a line from point `startTag` to point `endTag` and return the tag of the line.

`gmsh.model.geo.add_spline(pointTags: list) -> int` 
: Create a spline going through the points in `pointsTag` and return the tag of the spline.

`gmsh.model.geo.add_curve_loop(curveTags: list) -> int` 
: Create a curve loop formed by the curves `curveTags` and return the tag of the curve loop.

`gmsh.model.geo.add_plane_surface(wireTags: list) -> int` 
: Create a planar surface defined by one or more curve loops `wireTags` and return the tag of the surface. The first curve loop defines the exterior boundary, and any additional curve loops define holes. If the curve loops are not on a plane, no surface will be created.

`gmsh.model.geo.add_surface_filling(wireTags: list) -> int` 
: Create a surface, filling the curve loops in `wireTags` using transfinite interpolation, and return the tag of the surface. Only a single curve loop is supported and it must contain 3 or 4 curves.

`gmsh.model.geo.add_surface_loop(surfaceTags: list) -> int` 
: Create a surface loop formed by `surfaceTags` and return the tag of the surface.

`gmsh.model.geo.add_volume(shellTags: list) -> int` 
: Create a volume formed by `shellTags` and return the tag of the volume. The first surface loop defines the exterior boundary, and any additional surface loops define holes.

`gmsh.model.geo.extrude(dimTags: list, dx: float, dy: float, dz: float) -> list` 
: Extrude the entities `dimTags` (list of (dim, tag) tuples) by translating the entities along (dx, dy, dz) and returning a list of (dim, tag) tuples of the newly created entities.

`gmsh.model.geo.translate(dimTags: list, dx: float, dy: float, dz: float)` 
: Translate the entities `dimTags` (list of (dim, tag) tuples) along (dx, dy, dz).

`gmsh.model.geo.remove_all_duplicates()` 
: Remove all duplicate entities.

`gmsh.model.geo.split_curve(tag, pointTags) -> list` 
: Split the curve `tag` (line, spline, or b-spline) at the specified points `pointTags` and return a list of the newly created curves.

### OpenCASCADE geometry engine

`gmsh.model.occ.add_point(x: float, y: float, z: float) -> int` 
: Create a point at (x, y, z) and return the tag of the point.

`gmsh.model.occ.add_line(startTag: int, endTag: int) -> int` 
: Create a line from point `startTag` to point `endTag` and return the tag of the line.

`gmsh.model.occ.add_spline(pointTags: list) -> int` 
: Create a spline going through the points in `pointsTag` and return the tag of the spline.

`gmsh.model.occ.add_curve_loop(curveTags: list) -> int` 
: Create a curve loop formed by the curves `curveTags` and return the tag of the curve loop.

`gmsh.model.occ.add_rectangle(x, y, z, dx, dy) -> int` 
: Create a rectangle with lower left corner at `(x, y, z)` and upper right corner at `(x+dx, y+dy, z)` and return the tag of the rectangle.

`gmsh.model.occ.add_plane_surface(wireTags: list) -> int` 
: Create a planar surface defined by one or more curve loops `wireTags` and return the tag of the surface. The first curve loop defines the exterior boundary, and any additional curve loops define holes. If the curve loops are not on a plane, no surface will be created.

`gmsh.model.occ.add_surface_filling(wireTag: int, pointTags: list=[]) -> int` 
: Create a surface, filling the curve loop `wireTag` while going through points `pointTags` and return the tag of the surface.

`gmsh.model.occ.add_bspline_surface(pointTags: list, numPointsU: int) -> int` 
: Create a b-spline surface with control points `pointTags` (list with `numPointsU*numPointsV` points) and return the tag of the surface.

`gmsh.model.occ.add_trimmed_surface(surfaceTag: int, wireTags: list=[]) -> int` 
: Trim the surface `surfaceTag` with the wires `wireTags` and return the tag of the trimmed surface. The first wire defines the external contour, and any additional wires define holes.

`gmsh.model.occ.add_surface_loop(surfaceTags: list) -> int` 
: Create a surface loop formed by `surfaceTags` and return the tag of the surface.

`gmsh.model.occ.add_volume(shellTags: list) -> int` 
: Create a volume formed by `shellTags` and return the tag of the volume. The first surface loop defines the exterior boundary, and any additional surface loops define holes.

`gmsh.model.occ.add_sphere(x: float, y: float, z: float, radius: float) -> int` 
: Add a sphere centered at `(x, y, z)` with radius `radius` and return the tag of the sphere.

`gmsh.model.occ.add_box(x: float, y: float, z: float, dx: float, dy: float, dz: float) -> int` 
: Add a prism with a corner at `(x, y, z)` and extent `(dx, dy, dz)` and return the tag of the sphere.

`gmsh.model.occ.extrude(dimTags: list, dx: float, dy: float, dz: float) -> list` 
: Extrude the entities `dimTags` (list of (dim, tag) tuples) by translating the entities along (dx, dy, dz) and returning a list of (dim, tag) tuples of the newly created entities.

`gmsh.model.occ.translate(dimTags: list, dx: float, dy: float, dz: float)` 
: Translate the entities `dimTags` (list of (dim, tag) tuples) along (dx, dy, dz).

`gmsh.model.occ.fuse(objectDimTags: list, toolDimTags: list) -> tuple` 
: Compute the union of entities `objectDimTags` (list of (dim, tag) tuples) and `toolDimTags` (list of (dim, tag) tuples). The return value is a tuple with the a list of (dim, tag) tuples of the objects created and a mapping from the origin entities to the new entities.

`gmsh.model.occ.intersect(objectDimTags: list, toolDimTags: list) -> tuple` 
: Compute the intersection of entities `objectDimTags` (list of (dim, tag) tuples) and `toolDimTags` (list of (dim, tag) tuples). The return value is a tuple with the a list of (dim, tag) tuples of the objects created and a mapping from the origin entities to the new entities.

`gmsh.model.occ.cut(objectDimTags: list, toolDimTags: list) -> tuple` 
: Compute the difference of entities `objectDimTags` (list of (dim, tag) tuples) and `toolDimTags` (list of (dim, tag) tuples). The return value is a tuple with the a list of (dim, tag) tuples of the objects created and a mapping from the origin entities to the new entities.

`gmsh.model.occ.fragment(objectDimTags: list, toolDimTags: list) -> tuple` 
: Compute the boolean fragments resulting from the intersection of entities `objectDimTags` (list of (dim, tag) tuples) and `toolDimTags` (list of (dim, tag) tuples). The return value is a tuple with the a list of (dim, tag) tuples of the objects created and a mapping from the origin entities to the new entities.

`gmsh.model.occ.remove_all_duplicates()` 
: Remove all duplicate entities.

## Troubleshooting

### Using the Python debugger

You can use the Python debugger to have a Python script break at a specific point, so you can inspect the current state.
For example, insert `import pdb; pdb.set_trace()` at the desired breakpoint.

### Missing geometry

Did you forget to call `gmsh.model.geo.synchronize()` after creating the geometry?

### PETSc error when reading mesh

When a curve is not embedded in a surface or a surface is not embedded in a volume, Gmsh will generate independent meshes on each entity.
This will generate the following error when the mesh is read via the `MeshIOPetsc` mesh importer.

```{code-block} bash
[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: No support for this operation for this object type
[0]PETSC ERROR: Could not determine Plex facet for Gmsh element 23 (Plex cell 48)

-- lines omitted --

[0]PETSC ERROR: #1 DMPlexCreateGmsh() at petsc-dev/src/dm/impls/plex/plexgmsh.c:1766
[0]PETSC ERROR: #2 DMPlexCreateGmshFromFile() at petsc-dev/src/dm/impls/plex/plexgmsh.c:1470
[0]PETSC ERROR: #3 DMPlexCreateFromFile() at petsc-dev/src/dm/impls/plex/plexcreate.c:5889
[0]PETSC ERROR: #4 DMPlexCreateFromOptions_Internal() at petsc-dev/src/dm/impls/plex/plexcreate.c:4001
[0]PETSC ERROR: #5 DMSetFromOptions_Plex() at petsc-dev/src/dm/impls/plex/plexcreate.c:4418
[0]PETSC ERROR: #6 DMSetFromOptions() at petsc-dev/src/dm/interface/dm.c:914
[0]PETSC ERROR: #7 virtual void pylith::meshio::MeshIOPetsc::_read()() at src/cig/pylith/libsrc/pylith/meshio/MeshIOPetsc.cc:147
```
