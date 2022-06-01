# Examples: 2D Magma Reservoir Using Poroelasticity

This example demonstrates the use of poroelasticity. It includes a magma reservoir with
poroelastic properties that differ from the surrounding domain.

## Meshing

We provide mesh files generated using Cubit. We also include Cubit Journal files for
reproducing the provided finite-element meshes with triangular or quadrilateral cells.

## Step 1: Inflation

We impose flow into the conduit at the external boundary, which leads to inflation of the
magma reservoir.

To run the example:
```
pylith step01_inflation.cfg
```
