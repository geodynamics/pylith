(sec:format:SimpleGridDB)=
# *SimpleGridDB* Spatial Database Files

*SimpleGridDB* spatial database files contain a header describing the grid of points and then the data with each line listing the coordinates of a point followed by the values of the fields for that point.
The coordinates for each dimension of the grid do not need to be uniformly spaced.
The coordinate systems are specified the same way as they are in *SimpleDB* spatial database files as described in {ref}`sec:format:SimpleIOAscii`.

```{code-block} cfg
// This spatial database specifies the elastic properties on a
// 2-D grid in 3-D space.
//
// Comments can appear almost anywhere in these files and are
// delimited with two slashes (//) just like in C++. All text and
// whitespace after the delimiter on a given line is ignored.
// The next line is the magic header for spatial database files
// in ASCII format.
#SPATIAL_GRID.ascii 1
  SimpleGridDB { // start specifying the database parameters
    num-values = 3 // number of values in the database

    // Specify the names and the order of the values as they appear
    // in the data. The names of the values must correspond to the
    // names PyLith requests in querying the database.
    value-names = Vp Vs Density

    // Specify the units of the values in Python syntax.
    value-units = km/s km/s kg/m{*}{*}3
    num-x = 3 // Number of locations along x coordinate direction
    num-y = 1 // Number of locations along y coordinate direction
    num-z = 2 // Number of locations along z coordinate direction
    space-dim = 3 // Spatial dimension in which data resides

    // Specify the coordinate system associated with the
    // coordinates of the locations where data is given
    cs-data = cartesian \{ // Use a Cartesian coordinate system
      to-meters = 1.0e+3 // Coordinates are in km

      // Specify the spatial dimension of the coordinate system
      // This value must match the one associated with the database
      space-dim = 3
    } // cs-data // end of coordinate system specification
  } // end of SimpleGridDB specification
  // x coordinates
  -3.0 1.0 2.0
  // y coordinates
  8.0
  // z coordinates
  2.0 4.0
  // The locations and values are listed after the parameters.
  // Columns are coordinates of the points (1 column for each
  // spatial dimension) followed by the data values in the order
  // specified by the value-names field. The points can be in any order.
 -3.0 8.0 2.0 6.0 4.0 2500.0
  1.0 8.0 2.0 6.2 4.1 2600.0
  2.0 8.0 2.0 5.8 3.9 2400.0
 -3.0 8.0 4.0 6.1 4.1 2500.0
  1.0 8.0 4.0 5.9 3.8 2450.0
  2.0 8.0 4.0 5.7 3.7 2400.0
```
