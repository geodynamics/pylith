# Mesh Initialization

The `Initializer` object manages initializing the mesh, which involves reading it in from a file and performing some pre-processing operations.
The phases (steps) in the initialization process depend on whether the mesh can be read in parallel or it must be read in serial.

:::{admonition} Pyre User Interface
:class: seealso
See [`Initializer` component](../../components/initializers/Initializer.md)
:::

## Serial Phases

Phases for reading in a mesh in serial.

1. Read mesh
2. Reorder mesh
3. Distribute mesh
4. Insert fault interfaces
5. Refine mesh

:::{admonition} Pyre User Interface
:class: seealso
See [`Serial` component](../../components/initializers/Serial.md)
:::

## Parallel Phases

Phases for reading in a mesh in parallel.

1. Read mesh
2. Distribute mesh (fine tune distribution)
3. Insert fault interfaces
4. Refine mesh

:::{admonition} Pyre User Interface
:class: seealso
See [`Parallel` component](../../components/initializers/Parallel.md)
:::
