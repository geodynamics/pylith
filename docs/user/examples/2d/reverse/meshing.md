# Generating the Finite-Element Mesh

We use CUBIT to generate the finite-element mesh.
Due to the small size of these 2D meshes, we include them in the PyLith source and binary distributions.
If you do not have CUBIT, you can use the provided meshes.

Mesh generation is controlled from either `mesh_tri.jou` (triangular meshes) or `mesh_quad.jou` (quadrilateral meshes).
In addition to creating the desired meshes, these scripts call the following additional journal files:

* **geometry.jou** Journal file to create the 2D geometry.
* **gradient.jou** Journal file to assign sizing information for the mesh.
* **createbc.jou** Journal file to define material blocks and nodesets for boundary conditions.

The first step is to create the geometry.
This consists of creating a brick, extracting a midsurface from it, and then splitting the remaining surface with an extended fault and a splay surface.
The surfaces, curves, and important vertices are then assigned names that are then used when setting up mesh sizing information and defining blocks and nodesets.

:::{important}
We use IDless journaling in CUBIT.
This allows us to reference objects in a manner that should be independent of the version of CUBIT that is being used.
In the journal files, the original command used is typically commented out, and the following command is the equivalent IDless command.
:::

Once the geometry has been generated, we then set sizing information using both a user-defined sizing function as well as the CUBIT curve and surface bias schemes.
Using this sizing information, an initial mesh is generated, and then one iteration of smoothing is performed to improve the cell quality.
Finally, blocks are defined for the three materials in the problem, and nodesets are also defined for the fault and splay surfaces.

:::{important}
In addition to providing nodesets for the fault and splay, it is also important to provide nodesets defining the buried edges of these two surfaces.
In 2D, this will consist of a single vertex for each surface.
This information is required by PyLith to form the corresponding cohesive cells defining fault surfaces.
:::

Once you have run either the `mesh_tri.jou` or `mesh_quad.jou` journal file to construct the geometry and generate the mesh, you will have a corresponding Exodus-II file (`mesh_tri.exo` or `mesh_quad.exo`).
These are NetCDF files, and they can be loaded into ParaView.
This can be done by either running ParaView and loading the file, or using the script provided in the viz directory.
For example, if ParaView is in your path, you can run the
following command:

```{code-block} console
paraview --script=viz/plot_mesh.py
```

This will open ParaView, load the mesh, and produce views of both the mesh (with fault nodeset) and the mesh quality.
