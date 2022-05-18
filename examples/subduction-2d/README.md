# Examples: 2D Subduction

This suite of examples simulates quasistatic interseismic and
coseismic deformation for a vertical cross-section of a subduction zone.
It is based on the 2011 M9.0 Tohoku earthquake off the east coast of Japan.

For each of the simulations, we recommend examining the displacements,
stress field, and fault slip.

## Mesh generation using Gmsh (optional)

**WARNING**: The result of this step will overwrite the included file
tri_mesh.msh. You may want to copy/rename this file so that you have a
backup copy in case you have difficulty running Gmsh.

Run the Gmsh Python script:
```
generate_gmsh.py --write --filename=mesh_tri.msh
```

We highly recommend that you study the contents of the Gmsh Python script
to understand the mesh generation process.

## Step 1: Coseismic slip simulation

This simulation involves coseismic slip between the continental crust
and top of the subducting oceanic crust. The slip also extends down
into the top of the mantle below the continental crust.

The parameters for the earthquake slip are defined in
`fault_slip_coesismic.spatialdb`.

To run the example:
```
pylith step01_coesismic.cfg
```

## Step 2: Interseismic deformation simulation

This simulation involves aseismic creep along the interfaces between
the subducting oceanic crust and the mantle. The slip rate is a
constant 8 cm/yr.

The parameters for the creep along the top of the slab are defined in 
`fault_slabtop_creep.spatialdb`.

To run the example:
```
pylith step02_interseismic.cfg
```

## Step 3: Pseudo-earthquake cycle model

This simulation combines the interseismic deformation from Step 2
with the coseismic slip from Step 1. The earthquake rupture occurs
midway (150 years) through the 300 year simulation.

To run the example:
```
pylith step03_eqcycle.cfg
```

## Step 4: Friction controlled afterslip

**WARNING**: This example will not work with PyLith v3.0. It uses
spontaneous fault rupture, which has not been implemented in v3.0. Use
the PyLith v2.2.2 release and the example files included with it to
simulate spontaneous fault rupture.

This simulation uses the stress changes associated with the the
coseismic deformation in Step 1 in a simulation of afterslip
governed by static friction. The afterslip occurs at the down-dip
extent of rupture where the coseismic slip increases the shear
tractions.

Run the simulation via the following command:
```
pylith step04_afterslip.cfg
```

# Step 5: Earthquake cycle with slip-weakening friction

**WARNING**: This example will not work with PyLith v3.0. It uses
spontaneous fault rupture, which has not been implemented in v3.0. Use
the PyLith v2.2.2 release and the example files included with it to
simulate spontaneous fault rupture.

This simulation uses drives sponaneous rupture on the subducting
interface using prescribed asesismic slip on the bottom of the
slab. We use the slip-weakening friction model on the subducting
interface.

Run the simulation and capture the output via the following commands:
```
pylith step05_eqcycleslipweakening.cfg >& step05.log &
tail -f step05.log
```

# Step 6: Earthquake cycle with rate-state friction

**WARNING**: This example will not work with PyLith v3.0. It uses
spontaneous fault rupture, which has not been implemented in v3.0. Use
the PyLith v2.2.2 release and the example files included with it to
simulate spontaneous fault rupture.

This simulation uses drives sponaneous rupture on the subducting
interface using prescribed asesismic slip on the bottom of the
slab. We use the rate-state friction model on the subducting
interface.

Run the simulation and capture the output via the following commands:
```
pylith step06_eqcycleratestate.cfg >& step06.log &
tail -f step06.log
```
