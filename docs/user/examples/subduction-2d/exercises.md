# Suggested Exercises

1. Change the resolution of the mesh by editing the `generate_gmsh.py` Gmsh Python script or the `mesh_tri.jou` Cubit journal files. Change the resolution and bias factor.
2. Add depth dependent viscosity to the mantle and crust.
This requires using the linear Maxwell plane strain bulk constitutive model in the crust as well and creating spatial databases that include viscosity for the crust.
Specifying a depth dependent variation in the parameters will require adding points, updating num-locs accordingly, and changing data-dim to 1.
3. Modify the spatial database files for the material properties to use depth-dependent elastic properties based on PREM (Dziewonski and Anderson, 1981, 10.1016/0031-9201(81)90046-7). See <http://ds.iris.edu/ds/products/emc-prem/> for a simple table of values. Add points, update num-locs accordingly, and change data-dim to 1.
4. Create a Cubit journal file `mesh_quad.jou` for generating a mesh with quadrilateral cells instead of triangular cells. This requires using the pave mesh scheme.
