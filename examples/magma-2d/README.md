# Examples: 2D Magma Reservoir Using Poroelasticity

This example demonstrates the use of poroelasticity. It includes a magma reservoir with
poroelastic properties that differ from the surrounding domain.

## Meshing

We provide mesh files generated using Gmsh and Cubit. The parameter files are setup to
use the mesh from Gmsh.

## Step 1: Inflation

We impose flow into the conduit at the external boundary, which leads to inflation of the
magma reservoir.

To run the example:
```
pylith step01_inflation.cfg
```

## Step 2: Inflation with porosity state variable

Same as Step 1 with evolution of porosity state variable.

To run the example:
```
pylith step02_inflation_statevars.cfg
```
