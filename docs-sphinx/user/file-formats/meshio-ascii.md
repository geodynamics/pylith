(sec-user-file-formats-meshio-ascii)=
# PyLith Mesh ASCII File

PyLith mesh ASCII files allow quick specification of the mesh information for very small, simple meshes that are most easily written by hand.
We do not recommend using this format for anything other than these very small, simple meshes.

:::{figure-md} fig:meshioascii:diagram
<img src="figs/mesh_quad.*" alt="Diagram of mesh specified in the file."  width="400px"/>

Diagram of mesh specified in the file.
:::

```{code-block} c++
// This mesh file defines a finite-element mesh composed of two
// square cells of edge length 2.
//
// Comments can appear almost anywhere in these files and are
// delimited with two slashes (//) just like in C++. All text and
// whitespace after the delimiter on a given line is ignored.
mesh = { // begin specification of the mesh
  dimension = 2 // spatial dimension of the mesh
  // Begin vertex and cell labels with 0. This is the default so
  // this next line is optional
  use-index-zero = true

  vertices = { // vertices or nodes of the finite-element cells
      dimension = 2 // spatial dimension of the vertex coordinates
      count = 6 // number of vertices in the mesh
      coordinates = { // list of vertex index and coordinates
        // the coordinates must coincide with the coordinate
        // system specified in the Mesh object
        // exactly one vertex must appear on each line
        // (excluding whitespace)
        0 -2.0 -1.0
        1 -2.0 +1.0
        2 0.0 -1.0
        3 0.0 +1.0
        4 +2.0 -1.0
        5 +2.0 +1.0
      } // end of coordinates list
    } // end of vertices

    cells = { // finite-element cells
      count = 2 // number of cells in the mesh
      num-corners = 4 // number of vertices defining the cell
      simplices = { // list of vertices in each cell
        // see Section 4.2 for diagrams giving the order for each
        // type of cell supported in PyLith
        // index of cell precedes the list of vertices for the cell
        // exactly one cell must appear on each line
        // (excluding whitespace)
        0 0 2 3 1
        1 4 5 3 2
      } // end of simplices list

      material-ids = { // associated each cell with a material model
        // the material id is specified using the index of the cell
        // and then the corresponding material id
        0 0 // cell 0 has a material id of 0
        1 2 // cell 1 has a material id of 2
      } // end of material-ids list
    } // end of cells

    // This next section lists groups of vertices that can be used
    // in applying boundary conditions to portions of the domain
    group = { // start of a group
      // the name can have whitespace, so no comments are allowed
      // after the name
      name = face +y

      // Either groups of vertices or groups of cells can be
      // specified, but currently PyLith only makes use of groups
      // of vertices
      type = vertices // ’vertices’ or ’cells’
      count = 2 // number of vertices in the group
      indices = { // list of vertex indices in the group
        // multiple vertices may appear on a line
        0 4 // this group contains vertices 0 and 4
      } // end of list of vertices
    } // end of group

// additional groups can be listed here
```
