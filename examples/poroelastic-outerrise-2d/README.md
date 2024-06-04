# Examples: 2D Outer-Rise Hydration Using Poroelasticity

This example examines how bending stresses and permeability structure influence the hydration state of a slab of oceanic lithosphere in the outer-rise of a subduction zone.

## Meshing

We provide mesh files generated using Gmsh. We also include a Python script for generating
the finite-element mesh with triangular cells using Gmsh.

## Generate spatial databases

```bash
./generate_spatialdb_ypos.py
./generate_spatialdb_matprops.py
```

## Step 1: No faults, no flexure

The permeability field is depth dependent, decreasing with depth but does not vary laterally.
The lithosphere is not subject to any deformation, but a fluid pressure is applied to the top
boundary that is equivalent to the pressure exerted on the seafloor by the water column.

To run the example:

```bash
pylith step01_no_faults_no_flexure.cfg
```

## Step 2: No faults with flexure

The permeability field is depth dependent, decreasing with depth but does not vary laterally.
The lithosphere is now subject to deformation, over 300 kyr the slab bends to simulate extensional
stresses in the outer-rise of a subduction zone.

To run the example:

```bash
pylith step02_no_faults_flexure.cfg
```

## Step 1: Faults with flexure

The permeability field is depth dependent, decreasing with depth and also varies laterally,
simulating the enhanced permeability within normal faults in the outer-rise. The lithosphere is
subject to deformation, over 300 kyr the slab bends to simulate extensional stresses in the
outer-rise of a subduction zone.

To run the example:

```bash
step03_faults_flexure.cfg
```
