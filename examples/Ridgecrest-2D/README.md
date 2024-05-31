## Made during the 2023 Pylith Hackathon ##
# Contributors 
#   Zechao Zhuo: Gmesh
#   Evan Marschall: Cubit

# Example: 2D Simple Ridgecrest Model

This example folder has two examples for a static uniform slip model on a simple Ridecgrest fault model taken from geographic coordinates (Lat-Lon) and converted to a planar UTM system.

# Going from Geogrphic to UTM

We manually picked a few points along the desired faults and then used proj to convert the lat lon coordinates into  UTM coordinates. Both of these files are in the folder "Lon_Lat.txt" i the lat lon coordinates and "xy.txt" has the utm conversions


## Meshing

We provide mesh files generated using Gmsh and Cubit.
We also include a Python script for generating the finite-element mesh with
triangular cells using Gmsh and Journal files for generating the
finite-element mesh with triangular cell using Cubit.

****** Important note about meshes **************
The faults made in Gmsh are constructed from linear segements while the cubit meshed fualts are smoothed splines. This leads to some slight differences in surface deformation but are pretty much consistant.



## Step 1 Gmsh: Static Uniform Coseismic Slip

This example involves a static simulation that solves for the deformation from prescribed uniform coseismic slip on the Ridgecrest fault geometry. The slip on the main fault is 4 meters of right lateral slip, the large branch has 3 meters of left lateral slip and the small branch has 1 m of left lateral slip.

To run the example:
```
pylith step01_slip.cfg 

```

## Step 1 Cubit: Static Uniform Coseismic Slip

Same as Gmsh but for the cubit mesh. The fault line will be smoothed in comparison to gmsh

This example involves a static simulation that solves for the deformation from prescribed uniform coseismic slip on the Ridgecrest fault geometry. The slip on the main fault is 4 meters of right lateral slip, the large branch has 3 meters of left lateral slip and the small branch has 1 m of left lateral slip.

To run the example:
```
pylith step01_slip_cubit.cfg 

```


