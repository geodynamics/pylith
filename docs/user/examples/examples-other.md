# Additional Examples

## Cubit Meshing Examples

The directory `examples/meshing-cubit` contains several examples of using Cubit to construct finite-element meshes for complex geometry.
This includes features such as constructing nonplanar fault geometry from contours, constructing topography from a DEM, and merging sheet bodies (surfaces).
A separate examples discusses defining the discretization size using a vertex field in an Exodus-II file.
See the `README` files in the subdirectories for more detailed descriptions of these examples.

## Troubleshooting Examples

The directory `examples/troubleshooting-2d` contains a few examples to practice troubleshooting a variety of user errors in parameters files and problem setup.
The files with the errors corrected are in `examples/troubleshooting-2d/correct`.

## Code Verification Benchmarks

The CIG GitHub software repository <https://github.com/geodynamics/pylith_benchmarks> contains input files for a number of community benchmarks.
The benchmarks do not include the mesh files because they are so large; instead they include the Cubit journal files that can be used to generate the meshes.
Most, but not all, of the input files in the repository are updated for PyLith v2.0.0, so you will need to modify them if you use another version of PyLith.
