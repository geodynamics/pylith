(sec-user-physics-sources)=
# Point Sources

Point sources provide a mechanism for specifying internal sources within the domain, such as moment tensor sources for earthquake simulations or other force sources.
These sources are implemented at discrete points rather than on surfaces like boundary conditions or faults.

## Point Source Implementation

Point sources in PyLith are implemented using a moment tensor force representation.
The source is specified at discrete points in the domain, and the moment tensor components define the force couples applied at each point.
The temporal behavior of the source is controlled by a source time function, which can be an analytical wavelet or a user-specified time history.

## Specifying Point Sources

There are three basic steps in specifying point sources:

1. Create a file specifying the locations of the point sources.
2. Set the parameters for each source using `cfg` files or command line arguments.
3. Specify the spatial variation in source parameters (moment tensor, time delay) using a spatial database.
4. Select a source time function and specify its parameters.

## Arrays of Source Components

A dynamic array of source components associates a name (string) with each source.
The default source is `MomentTensorForce`.

```{code-block} cfg
---
caption: Array of sources in a `cfg` file.
---
[pylithapp.problem]
# Array of point sources
sources = [source1, source2]

# Default source is MomentTensorForce
# Both sources use the default type
```

## Source Implementations

:::{toctree}
moment-tensor-force.md
:::

