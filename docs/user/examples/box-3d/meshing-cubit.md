# Cubit Mesh

## Overview

We create the geometry, set the mesh parameters, generate the mesh, mark the materials and boundaries, and then write the mesh to a file.

## Meshing using Python Script

*New in v4.1.0*

We use the Python script `generate_cubit.py` to generate the mesh.
The Python script is setup so that it can be run from within Cubit or as a standalone Python script without the Cubit GUI interface.
In this example, we will run the script from within Cubit using the Journal editor.

Open the Python script `generate_cubit.py` in the Cubit journal editor.
Play the selected script or play the lines, making sure you play the first line so that Cubit uses the Python interpreter when running the script.
We specify the parameters controlling the geometry, mesh size, and cell shape near the top of the script.

:::{important}
The numbering of geometric entities, such as points, curves, and surfaces can depend upon the version of the ACIS geometry library used in Cubit.
Strategies to make scripts independent of the version of the ACIS geometry library being used include:
1. Use a variable to store the geometry `id`.
2. Use "idless" scripts that refer to entities by location rather than `id`. Refer to the Cubit documentation for more informaiton.
3. Name entities when they are created to minimize use of the entity `id` in scripts.
We often make use of all three of these methods.
:::

Once you have run the Python script to construct the geometry and generate the mesh, you will have a corresponding Exodus-II file (`mesh_tet.exo` or `mesh_hex.exo`).
These are NetCDF files, and they can be loaded into ParaView.

:::{figure-md} fig:example:box:3d:cubit:hex
<img src="figs/cubit-hex.*" alt="Finite-element mesh with hexahedral cells generated by Cubit." width="400px"/>

Finite-element mesh with hexahedral cells generated by Cubit.
We show the nodeset on the `boundary_zpos` boundary.
:::

:::{figure-md} fig:example:box:3d:cubit:tet
<img src="figs/cubit-tet.*" alt="Finite-element mesh with tetrahedral cells generated by Cubit." width="400px"/>

Finite-element mesh with tetrahedral cells generated by Cubit.
We show the nodeset on the `boundary_zpos` boundary.
:::
