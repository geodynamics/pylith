(sec-user-physics-materials)=
# Materials

The material objects encapsulate the bulk behavior of the domain.
This includes both the governing equations and the associated bulk rheology.

## Specifying Material Properties

Associating material properties with a given cell involves several steps.

1. In the mesh generation process, assign a material identifier to each cell.

2. Define material property groups corresponding to each material identifier.
    In CUBIT/Trelis this is done by creating the blocks as part of the boundary conditions.

3. Provide the settings for each material group in the parameters, i.e., `cfg` file.

4. Specify the parameters for the material properties, e.g., linear variation in density with depth, using a spatial database file.
This allows variation of the material properties across cells with the same material identifier.

## Material Implementations

:::{toctree}
elasticity.md
incompressible-elasticity.md
poroelasticity.md
:::

:::{seealso}
See {ref}`sec-user-governing-eqns` for the derivation of the finite-element formulation for each of the materials.
:::
