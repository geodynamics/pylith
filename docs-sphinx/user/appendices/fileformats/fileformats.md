## File Formats

### PyLith Mesh ASCII Files

PyLith mesh ASCII files allow quick specification of the mesh information for
very small, simple meshes that are most easily written by hand. We do not
recommend using this format for anything other than these very small, simple
meshes.

<figure>
<img src="fileformats/figs/meshquad4.*" id="fig:meshioascii:diagram" alt="Diagram of mesh specified in the file." /><figcaption aria-hidden="true">Diagram of mesh specified in the file.</figcaption>
</figure>

###  Spatial Database Files

spatial database files contain a header describing the set of points and then
the data with each line listing the coordinates of a point followed by the
values of the fields for that point.

<div class="SimpleIOAscii">

// This spatial database specifies the distribution of slip on the // fault
surface. In this case we prescribe a piecewise linear, // depth dependent
distribution of slip. The slip is 2.0 m // right-lateral with 0.25 m of
reverse slip at the surface with // a linear taper from 2.0 m to 0.0 m from -2
km to -4 km. // // Comments can appear almost anywhere in these files and are
// delimited with two slashes (//) just like in C++. All text and //
whitespace after the delimiter on a given line is ignored. // // The next line
is the magic header for spatial database files // in ASCII format.
#SPATIAL.ascii 1 SimpleDB

// start specifying the database parameters num-values = 3 // number of values
in the database

// Specify the names and the order of the values as they appear // in the
data. The names of the values must correspond to the // names PyLith requests
in querying the database. value-names = left-lateral-slip reverse-slip
fault-opening

// Specify the units of the values in Python syntax (e.g., kg/m\*\*3).
value-units = m m m num-locs = 3 // Number of locations where values are given
data-dim = 1 // Locations of data points form a line space-dim = 3 // Spatial
dimension in which data resides

// Specify the coordinate system associated with the // coordinates of the
locations where data is given cs-data = cartesian

// Use a Cartesian coordinate system to-meters = 1.0e+3 // Coordinates are in
km

// Specify the spatial dimension of the coordinate system // This value must
match the one associated with the database space-dim = 3

// cs-data // end of coordinate system specification

// end of SimpleDB parameters // The locations and values are listed after the
parameters. // Columns are coordinates of the points (1 column for each //
spatial dimension) followed by the data values in the order // specified by
the value-names field. 0.0 0.0 0.0 -2.00 0.25 0.00 0.0 0.0 -2.0 -2.00 0.00
0.00 0.0 0.0 -4.0 0.00 0.00 0.00

</div>

#### Spatial Database Coordinate Systems

The spatial database files support four different types of coordinate systems.
Conversions among the three geographic coordinate systems are supported in 3D.
Conversions among the geographic and geographic projected coordinate systems
are supported in 2D. In using the coordinate systems, we assume that you have
installed the program in addition to the Proj.4 libraries, so that you can
obtain a list of support projections, datums, and ellipsoids. Alternatively,
refer to the documentation for the Proj.4 Cartographic Projections library
[trac.osgeo.org/proj][].

##### Cartesian

This is a conventional Cartesian coordinate system. Conversions to other
Cartesian coordinate systems are possible.

<div class="SimpleIOAscii">

cs-data = cartesian to-meters = 1.0e+3 // Locations in km space-dim = 2 // 1,
2, or 3 dimensions

</div>

##### Geographic

This coordinate system is for geographic coordinates, such as longitude and
latitude. Specification of the location in three-dimensions is supported. The
vertical datum can be either the reference ellipsoid or mean sea level. The
vertical coordinate is positive upwards.

<div class="SimpleIOAscii">

cs-data = geographic

// Conversion factor to get to meters (only applies to vertical // coordinate
unless you are using a geocentric coordinate system). to-meters = 1.0e+3
space-dim = 2 // 2 or 3 dimensions

// Run &ldquo;proj -le&rdquo; to see a list of available reference ellipsoids.
// Comments are not allowed at the end of the next line. ellipsoid = WGS84

// Run &ldquo;proj -ld&rdquo; to see a list of available datums. // Comments
are not allowed at the end of the next line. datum-horiz = WGS84

// &ldquo;ellipsoid&rdquo; or &ldquo;mean sea level&rdquo; // Comments are not
allowed at the end of the next line. datum-vert = ellipsoid

// Use a geocentric coordinate system? is-geocentric = false // true or false

</div>

##### Geographic Projection

This coordinate system applies to geographic projections. As in the geographic
coordinate system, the vertical coordinate (if used) can be specified with
respect to either the reference ellipsoid or mean sea level. The coordinate
system can use a local origin and be rotated about the vertical direction.

##### Geographic Local Cartesian

This coordinate system is a geographically referenced, local 3D Cartesian
coordinate system. This allows use of a conventional Cartesian coordinate
system with accurate georeferencing. For example, one can construct a
finite-element model in this coordinate system and use spatial databases in
geographic coordinates. This coordinate system provides an alternative to
using a geographic projection as the Cartesian grip. The advantage of this
coordinate system is that it retains the curvature of the Earth, while a
geographic projection does not.

<div class="SimpleIOAscii">

cs-data = geo-local-cartesian

// Conversion factor to get to meters (only applies to vertical // coordinate
unless you are using a geocentric coordinate system). to-meters = 1.0 // use
meters space-dim = 2 // 2 or 3 dimensions

// Run &ldquo;proj -le&rdquo; to see a list of available reference ellipsoids.
// Comments are not allowed at the end of the next line. ellipsoid = WGS84

// Run &ldquo;proj -ld&rdquo; to see a list of available datums. // Comments
are not allowed at the end of the next line. datum-horiz = WGS84

// &ldquo;ellipsoid&rdquo; or &ldquo;mean sea level&rdquo; // Comments are not
allowed at the end of the next line. datum-vert = ellipsoid

// Origin of the local Cartesian coordinate system. To avoid // round-off
errors it is best to pick a location near the center of // the region of
interest. An elevation on the surface of the Earth // in the middle of the
region also works well (and makes the // vertical coordinate easy to
interpret). origin-lon = -116.7094 // Longitude of the origin in decimal
degrees // (west is negative).

origin-lat = 36.3874 // Latitude of the origin in decimal degrees // (north is
positive).

// Elevation with respect to the vertical datum. Units are the // same as the
Cartesian coordinate system (in this case meters). origin-elev = 3.5

</div>

###  Spatial Database Files

spatial database files contain a header describing the grid of points and then
the data with each line listing the coordinates of a point followed by the
values of the fields for that point. The coordinates for each dimension of the
grid do not need to be uniformly spaced. The coordinate systems are specified
the same way as they are in spatial database files as described in Section
[\[sec:format:SimpleIOAscii\]][1].

###  Database Files

database files contain a header describing the number of points in the time
history and the units for the time stamps followed by a list with pairs of
time stamps and amplitude values. The amplitude at an arbitrary point in time
is computed via interpolation of the values in the database. This means that
the time history database must span the range of time values of interest. The
points in the time history must also be ordered in time.

<div class="TimeHistoryIO">

// This time history database specifies temporal variation in // amplitude. In
this case we prescribe a triangular slip time // history. // // Comments can
appear almost anywhere in these files and are // delimited with two slashes
(//) just like in C++. All text and // whitespace after the delimiter on a
given line is ignored. // // The next line is the magic header for spatial
database files // in ASCII format. #TIME HISTORY ascii TimeHistory

// start specifying the database parameters num-points = 5 // number of points
in time history

// Specify the units used in the time stamps. time-units = year

// end of TimeHistory header // The time history values are listed after the
parameters. // Columns time and amplitude where the amplitude values are
unitless. 0.0 0.00 2.0 1.00 6.0 4.00 10.0 2.00 11.0 0.00

</div>

### User-Specified Time-Step File

This file lists the time-step sizes for nonuniform, user-specified time steps
associated with the object. The file&rsquo;s format is an ASCII file that
includes the units for the time-step sizes and then a list of the time steps.

<div class="TimeStepUser">

// This time step file specifies five time steps with the units in years. //
Comments can appear almost anywhere in these files and are // delimited with
two slashes (//) just like in C++. All text and // whitespace after the
delimiter on a given line is ignored. // // Units for the time steps units =
year 1.0 // Comment 2.0 3.0 2.5 3.0

</div>

###  File

This file lists the coordinates of the locations where output is requested for
the `OutputSolnPoints` component. The coordinate system is specified in the
`OutputSolnPoints` component.

<div class="PointsList">

# Comments are limited to complete lines. The default delimiter for comments #
is &rsquo;#&rsquo;, which can be changed via parameters. Additionally, the
delimiter # separating values can also be customized (default is whitespace).
# # The first column is the station name. The coordinates of the points are
given # in the subsequent columns. P0 1.0 -2.0 0.0 P1 2.0 -4.0 -0.1 P2 0.0
+2.0 0.0 P3 2.5 -0.2 -0.2 P4 0.0 2.0 +0.2

</div>

  [trac.osgeo.org/proj]: trac.osgeo.org/proj
  [1]: #sec:format:SimpleIOAscii
