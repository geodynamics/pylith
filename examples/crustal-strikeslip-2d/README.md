# Examples: Crustal 2D strike-slip faults

This suite of examples simulates 2D (plane strain) quasistatic coseismic deformation
for a set of intersecting crustal strike-slip faults. It is based on the 2019 Ridgecrest
earthquake sequence.

## Mesh generation using Gmsh (optional)

**WARNING**: The result of this step will overwrite the included file
`mesh_tri.msh`. You may want to copy/rename this file so that you have a
backup copy in case you have difficulty running Gmsh.

Run the Gmsh Python script:

```bash
generate_gmsh.py --write
```

We highly recommend that you study the contents of the Gmsh Python script
to understand the mesh generation process.

## Mesh generation using Cubit (optional)

**WARNING**: The result of this step will overwrite the included file
`mesh_tri.exo`. You may want to copy/rename this file so that you have a
backup copy in case you have difficulty running Cubit.

Start the Cubit application and load the Python script `generate_cubit.py`
into the Cubit Journal editor. Make sure the current directory is set to this
directory and then run the Python script.

We highly recommend that you study the contents of the Cubit Python script
to understand the mesh generation process.

## Step 1: Coseismic slip simulation

This simulation involves uniform coseismic slip on the three faults. We provide
parameter files for finite-element meshes generated for both Gmsh (`step01_slip.cfg`)
and Cubit (`step01_slip_cubit.cfg`).

To run the example:

```bash
# Gmsh
pylith step01_slip.cfg

# Cubit
pylith step01_slip_cubit.cfg
```
