# Examples: 2-D Subduction

This suite of examples simulates quasi-static interseismic and
coseismic deformation for a subduction zone. It is based on the 2011
M9.0 Tohoku earthquake off the east coast of Japan.

The main features of this example are:

* Generating a finite-element mesh using CUBIT/Trelis
  * Nonplanar geometry
  * Variable mesh resolution
* Spatially variable coseismic slip and creep
* Maxwell viscoelastic relaxation

The example is broken up into three steps of increasing
complexity. Step 1 focuses on the coseismic slip, Step 2 focuses on
interseismic deformation, and Step 3 combines the two into a pseudo
earthquake cycle deformation simulation. In a similar fashion as other
examples, we place all of the general parameters in `pylithapp.cfg` and
the parameters specific to each simulation in `stepXX_XXXX.cfg`.

We model the crust with a linear elastic bulk constitutive model and
the mantle with a linear Maxwell viscoelastic model. Both constitutive
models use plane strain formulations in these 2-D models.

The parameters for the bulk constitutive models are defined in
  * `mat_concrust.spatialdb`
  * `mat_oceancrust.spatialdb`
  * `mat_conmantle.spatialdb`
  * `mat_oceanmantle.spatialdb`

The simulation will output the displacements on the ground surface at
every time step, the displacements over the entire domain every 10
years, the fault slip and tractions every time step, and the stresses
and strains for each material every 10 years.

For each of the simulations, we recommend examining the displacements,
stress field, and fault slip.

## Mesh generation using CUBIT/Trelis (optional)

**WARNING**: The result of this step will overwrite the included file
tri_mesh.exo. You may want to copy/rename this file so that you have a
backup copy in case you have difficulty running CUBIT/Trelis.

Start CUBIT/Trelis and play the journal file `mesh_tri.jou`. We highly
recommend that you study the contents of the journal files to
understand the mesh generation process.


## Step01: Coseismic slip simulation

This simulation involves coseismic slip between the continental crust
and top of the subducting oceanic crust. The slip also extends down
into the top of the mantle below the continental crust.

The parameters for the earthquake slip are defined in
`fault_slip_coesismic.spatialdb`.

To run the example:
```
pylith step01_coesismic.cfg
```


## Step02: Interseismic deformation simulation

This simulation involves aseismic creep along the interfaces between
the subducting oceanic crust and the mantle. The slip rate is a
constant 8 cm/yr.

The parameters for the creep along the top of the slab are defined in 
`fault_slabtop_creep.spatialdb`.

To run the example:
```
pylith step02_interseismic.cfg
```


## Step03: Pseudo-earthquake cycle model

This simulation combines the interseismic deformation from Step 2
with the coseismic slip from Step 1. The earthquake rupture occurs
midway (150 years) through the 300 year simulation.

To run the example:
```
pylith step03_eqcycle.cfg
```

## Step04: Friction controlled afterslip

**IMPORTANT**: This example will not work with PyLith v3.0. It uses
spontaneous fault rupture, which has not been implemented in v3.0. Use
the PyLith v2.2.1 release to simulate spontaneous fault rupture.

This simulation uses the stress changes associated with the the
coseismic deformation in Step 1 in a simulation of afterslip
governed by static friction. The afterslip occurs at the down-dip
extent of rupture where the coseismic slip increases the shear
tractions.

Run the simulation via the following command:
```
pylith step04_afterslip.cfg
```


# Step05: Earthquake cycle with slip-weakening friction

**IMPORTANT**: This example will not work with PyLith v3.0. It uses
spontaneous fault rupture, which has not been implemented in v3.0. Use
the PyLith v2.2.1 release to simulate spontaneous fault rupture.

This simulation uses drives sponaneous rupture on the subducting
interface using prescribed asesismic slip on the bottom of the
slab. We use the slip-weakening friction model on the subducting
interface.

Run the simulation and capture the output via the following commands:
```
pylith step05_eqcycleslipweakening.cfg >& step05.log &
tail -f step05.log
```


# Step06: Earthquake cycle with rate-state friction

**IMPORTANT**: This example will not work with PyLith v3.0. It uses
spontaneous fault rupture, which has not been implemented in v3.0. Use
the PyLith v2.2.1 release to simulate spontaneous fault rupture.

This simulation uses drives sponaneous rupture on the subducting
interface using prescribed asesismic slip on the bottom of the
slab. We use the rate-state friction model on the subducting
interface.

Run the simulation and capture the output via the following commands:
```
pylith step06_eqcycleratestate.cfg >& step06.log &
tail -f step06.log
```

