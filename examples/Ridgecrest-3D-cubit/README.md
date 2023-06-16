## Made during the 2023 Pylith Hackathon ##
# Contributors 
#   Evan Marschall
#   Zechao Zhuo

# Examples: 3D Ridgecrest 

This example is similiar to the 2D case, but done in 3D. Note the mesh is revlatively shallow (7.5 km) because we were limited to 50,000 elements for the cubit student license.



# Going from Geogrphic to UTM

We manually picked a few points along the desired faults and then used proj to convert the lat lon coordinates into  UTM coordinates. Both of these files are in the folder "Lon_Lat.txt" i the lat lon coordinates and "xy.txt" has the utm conversions


## Meshing

We provide mesh files generated using Cubit.
We also include  Journal files for generating the
finite-element mesh with triangular cell using Cubit.



## Step 1 Cubit: Static Uniform Coseismic Slip

This example involves a static simulation that solves for the deformation from prescribed uniform coseismic slip on the Ridgecrest fault geometry. The slip on the main fault is 4 meters of right lateral slip, the large branch has 3 meters of left lateral slip and the small branch has 1 m of left lateral slip.

To run the example:
```
pylith step01_slip.cfg 

```


