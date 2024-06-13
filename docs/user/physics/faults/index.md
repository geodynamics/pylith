(sec-user-physics-faults)=
# Fault Interface Conditions

Fault interfaces are used to create dislocations (jumps in the displacement field) in the model.
The dislocations arise from slip across a fault surface.
Both shear and tensile dislocations are supported.
Dislocations in 2D are specified in terms of left-lateral-slip and fault opening, and in 3D they are specified in terms of left-lateral-slip, reverse-slip, and fault opening.
PyLith supports kinematic (prescribed) slip and dynamic (spontaneous) rupture simulations.

:::{warning}
Spontaneous rupture is not available in PyLith v3.0; we plan to have it reimplemented in v3.1
:::

## Conventions

Slip corresponds to relative motion across a fault surface.
{numref}`fig:fault:orientation` shows the orientation of the slip vector in 3D with respect to the fault surface and coordinate axes.
PyLith automatically determines the local orientation of the fault surface.
This alleviates the user from having to compute the strike, dip, and rake angles over potentially complex, nonplanar fault surfaces.
Instead, the user specifies fault parameters in terms of lateral motion, reverse motion, and fault opening as shown in {numref}`fig:fault:slip:motions`.

:::{figure-md} fig:fault:orientation
<img src="figs/fault-orientation.*" alt="Orientation of a fault surface in 3D, where $\phi$ denotes the angle of the fault strike, $\delta$ denotes the angle of the fault dip, and $\lambda$ the rake angle." width="600px"/>

Orientation of a fault surface in 3D, where $\phi$ denotes the angle of the fault strike, $\delta$ denotes the angle of the fault dip, and $\lambda$ the rake angle.
:::

:::{figure-md} fig:fault:slip:motions
<img src="figs/slip-motions.*" alt="Sign conventions associated with fault slip. Positive values are associated with left-lateral, reverse, and fault opening motions." width="600px"/>

Sign conventions associated with fault slip.
Positive values are associated with left-lateral, reverse, and fault opening motions.
:::

## Fault Implementation

In order to create relative motion across the fault surface in the finite-element mesh, PyLith adjusts the topology of the mesh by inserting zero area (2D) or volume (3D) cohesive cells along the fault surface.
The cohesive cells allow control of the relative motion between vertices on the two sides of the fault.
{numref}`fig:fault:cohesive:cells` illustrates the results of inserting cohesive cells in a mesh consisting of triangular cells.
This example also shows the distinction between how buried fault edges are handled differently than fault edges that reach the edge of the domain, such as the ground surface.

:::{figure-md} fig:fault:cohesive:cells
<img src="figs/cohesive-cell.*" alt="Example of cohesive cells inserted into a mesh of triangular cells. The zero thickness cohesive cells control slip on the fault via the relative motion between the vertices on the positive and negative sides of the fault." width="800px" />

Example of cohesive cells inserted into a mesh of triangular cells.
The zero thickness cohesive cells control slip on the fault via the relative motion between the vertices on the positive and negative sides of the fault.
:::

:::{figure-md} fig:fault:fault_edge
<img src="figs/fault-edge.*" alt="Example of how faults with buried edges must be described with two sets of vertices. All of the vertices on the fault are included in the `fault` group; the subset of vertices along the buried edges are included in the `fault_edge` group. In 2D the fault edges are just a single vertex as shown in {numref}`fig:fault:cohesive:cells`." width="600px"/>

Example of how faults with buried edges must be described with two sets of vertices.
All of the vertices on the fault are included in the `fault` group; the subset of vertices along the buried edges are included in the `fault_edge` group.
In 2D the fault edges are just a single vertex as shown in {numref}`fig:fault:cohesive:cells`
:::

For faults that have buried edges, splitting the mesh apart and inserting the cohesive cells becomes complex at the buried edges due to the ambiguity of defining where the fault ends and how to insert the cohesive cell.
Starting in PyLith v2.0.0 we changed how the buried edges of the fault are managed.
An additional group of fault nodes is specified (for example, via a nodeset from Cubit) that marks the buried edges of the fault (see {numref}`fig:fault:fault_edge`).
This allows the cohesive cell insertion algorithm to adjust the topology so that cohseive cells are inserted up to the buried edge of the fault but no additional degrees of freedom are added on the fault edge.
This naturally forces slip to zero along the buried edges.

In 2D the default in-plane slip is left-lateral, so we use the reference directions to resolve the ambiguity in specifying reverse slip.
In 3D the reference directions are used to resolve the ambiguity in the along-strike and dip-dir directions.
If the fault plane is horizontal, then the up-dir corresponds to the reverse-motion on the +z side of the fault.
The only requirement for this direction is that it not be colinear with the fault normal direction.
The default value of [0, 0, 1] is appropriate for most 3D problems.

By default the output observers write both diagnostic information (for example, fault orientation directions) and the slip at each time step.
The fault coordinate system is shown in {numref}`fig:fault:slip:motions`.
The vectors in the fault coordinate system can be transformed to the global coordinate system using the direction vectors in the diagnostic output ({ref}`sec-user-output-observers`).

:::{important}
The normal direction is chosen based on how the cells are split to create cohesive cells.
If the normal direction contains a positive z component, then the directions conform to traditional seismologic conventions (along strike and up dip); however, if the normal direction contains a negative z component, then the directions correspond to along strike and down dip directions.
:::

:::{toctree}
prescribed-slip.md
slip-impulses.md
:::
