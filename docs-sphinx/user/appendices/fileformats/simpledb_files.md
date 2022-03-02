(sec:format:SimpleIOAscii)=
# *SimpleDB* Spatial Database Files

*SimpleDB* spatial database files contain a header describing the set of points and then the data with each line listing the coordinates of a point followed by the values of the fields for that point.

```{code-block} cfg
// This spatial database specifies the distribution of slip on the
// fault surface. In this case we prescribe a piecewise linear,
// depth dependent distribution of slip. The slip is 2.0 m
// right-lateral with 0.25 m of reverse slip at the surface with
// a linear taper from 2.0 m to 0.0 m from -2 km to -4 km.
//
// Comments can appear almost anywhere in these files and are
// delimited with two slashes (//) just like in C++. All text and
// whitespace after the delimiter on a given line is ignored.
//
// The next line is the magic header for spatial database files
// in ASCII format.
#SPATIAL.ascii 1
SimpleDB { // start specifying the database parameters
  num-values = 3 // number of values in the database

  // Specify the names and the order of the values as they appear
  // in the data. The names of the values must correspond to the
  // names PyLith requests in querying the database.
  value-names = left-lateral-slip reverse-slip fault-opening

  // Specify the units of the values in Python syntax (e.g., kg/m**3).
  value-units = m m m
  num-locs = 3 // Number of locations where values are given
  data-dim = 1 // Locations of data points form a line
  space-dim = 3 // Spatial dimension in which data resides

  // Specify the coordinate system associated with the
  // coordinates of the locations where data is given
  cs-data = cartesian { // Use a Cartesian coordinate system
    to-meters = 1.0e+3 // Coordinates are in km

    // Specify the spatial dimension of the coordinate system
    // This value must match the one associated with the database
    space-dim = 3

  } // cs-data // end of coordinate system specification

} // end of SimpleDB parameters
// The locations and values are listed after the parameters.
// Columns are coordinates of the points (1 column for each
// spatial dimension) followed by the data values in the order
// specified by the value-names field.
0.0 0.0 0.0 -2.00 0.25 0.00
0.0 0.0 -2.0 -2.00 0.00 0.00
0.0 0.0 -4.0 0.00 0.00 0.00
```

## Spatial Database Coordinate Systems

The spatial database files support four different types of coordinate systems.
Conversions among the three geographic coordinate systems are supported in 3D.
Conversions among the geographic and geographic projected coordinate systems are supported in 2D.
In using the coordinate systems, we assume that you have installed the `proj` program in addition to the Proj.4 libraries, so that you can obtain a list of support projections, datums, and ellipsoids.
Alternatively, refer to the documentation for the Proj.4 Cartographic Projections library <https://proj.org/>.

### Cartesian

This is a conventional Cartesian coordinate system.
Conversions to other Cartesian coordinate systems are possible.

```{code-block} cfg
cs-data = cartesian {
  to-meters = 1.0e+3 // Locations in km
  space-dim = 2 // 1, 2, or 3 dimensions
}
```

### Geographic

This coordinate system is for geographic coordinates, such as longitude and latitude.
Specification of the location in three-dimensions is supported.
The vertical datum can be either the reference ellipsoid or mean sea level.
The vertical coordinate is positive upwards.

```{code-block} cfg
cs-data = geographic {
  // Conversion factor to get to meters (only applies to vertical
  // coordinate unless you are using a geocentric coordinate system).
  to-meters = 1.0e+3
  space-dim = 2 // 2 or 3 dimensions

  // Run ‘‘proj -le’’ to see a list of available reference ellipsoids.
  // Comments are not allowed at the end of the next line.
  ellipsoid = WGS84

  // Run ‘‘proj -ld’’ to see a list of available datums.
  // Comments are not allowed at the end of the next line.
  datum-horiz = WGS84

  // ‘‘ellipsoid’’ or ‘‘mean sea level’’
  // Comments are not allowed at the end of the next line.
  datum-vert = ellipsoid

  // Use a geocentric coordinate system?
  is-geocentric = false // true or false
}
```

### Geographic Projection

This coordinate system applies to geographic projections.
As in the geographic coordinate system, the vertical coordinate (if used) can be specified with respect to either the reference ellipsoid or mean sea level.
The coordinate system can use a local origin and be rotated about the vertical direction.

```{code-block} cfg
cs-data = geo-projected {
  to-meters = 1.0e+3 // Conversion factor to get to meters.
  space-dim = 2 // 2 or 3 dimensions

  // Run ‘‘proj -le’’ to see a list of available reference ellipsoids.
  // Comments are not allowed at the end of the next line.
  ellipsoid = WGS84

  // Run ‘‘proj -ld’’ to see a list of available datums.
  // Comments are not allowed at the end of the next line.
  datum-horiz = WGS84

  // ‘‘ellipsoid’’ or ‘‘mean sea level’’
  // Comments are not allowed at the end of the next line.
  datum-vert = ellipsoid

  // Longitude of local origin in WGS84.
  origin-lon = -120.0

  // Latitude of local origin in WGS84.
  origin-lat = 37.0

  // Rotation angle in degrees (CCW) or local x-axis from east.
  rotation-angle = 23.0

  // Run ‘‘proj -lp’’ to see a list of available geographic
  // projections.
    projector = projection {
    // Name of the projection. run ‘‘proj -lp’’ to see a list of
    // supported projections. Use the Universal Transverse Mercator
    // projection.
    projection = utm
    units = m // Units in the projection.

    // Provide a list of projection options; these are the command
    // line arguments you would use with the proj program. Refer to
    // the Proj.4 library documentation for complete details.
    // Comments are not allowed at the end of the next line.
    proj-options = +zone=10
}
```

### Geographic Local Cartesian

This coordinate system is a geographically referenced, local 3D Cartesian coordinate system.
This allows use of a conventional Cartesian coordinate system with accurate georeferencing.
For example, one can construct a finite-element model in this coordinate system and use spatial databases in geographic coordinates.
This coordinate system provides an alternative to using a geographic projection as the Cartesian grip.
The advantage of this coordinate system is that it retains the curvature of the Earth, while a geographic projection does not.

```{code-block} cfg
cs-data = geo-local-cartesian {
  // Conversion factor to get to meters (only applies to vertical
  // coordinate unless you are using a geocentric coordinate system).
  to-meters = 1.0 // use meters
  space-dim = 2 // 2 or 3 dimensions

  // Run ‘‘proj -le’’ to see a list of available reference ellipsoids.
  // Comments are not allowed at the end of the next line.
  ellipsoid = WGS84

  // Run ‘‘proj -ld’’ to see a list of available datums.
  // Comments are not allowed at the end of the next line.
  datum-horiz = WGS84

  // ‘‘ellipsoid’’ or ‘‘mean sea level’’
  // Comments are not allowed at the end of the next line.
  datum-vert = ellipsoid

  // Origin of the local Cartesian coordinate system. To avoid
  // round-off errors it is best to pick a location near the center of
  // the region of interest. An elevation on the surface of the Earth
  // in the middle of the region also works well (and makes the
  // vertical coordinate easy to interpret).
  origin-lon = -116.7094 // Longitude of the origin in decimal degrees
                        // (west is negative).

  origin-lat = 36.3874 // Latitude of the origin in decimal degrees
                      // (north is positive).

  // Elevation with respect to the vertical datum. Units are the
  // same as the Cartesian coordinate system (in this case meters).
  origin-elev = 3.5
}
```
