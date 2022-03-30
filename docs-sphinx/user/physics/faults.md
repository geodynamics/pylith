# Fault Interface Conditions

Fault interfaces are used to create dislocations (jumps in the displacement field) in the model.
The dislocations arise from slip across a fault surface.
Both shear and tensile dislocations are supported.
For fault interfaces, dislocations in 2D left-lateral-slip and fault opening, and in 3D left-lateral-slip, reverse-slip, and fault opening.
PyLith supports kinematic (prescribed) slip and dynamic (spontaneous) rupture simulations.

:::{warning}
Spontaneous rupture is not available in PyLith v3.0; we plan to have it reimplemented in v3.1
:::

## Conventions

Slip corresponds to relative motion across a fault surface.
{numref}`fig:fault:orientation` shows the orientation of the slip vector in 3D with respect to the fault surface and coordinate axes.
PyLith automatically determines the orientation of the fault surface.
This alleviates the user from having to compute the strike, dip, and rake angles over potentially complex, nonplanar fault surfaces.
Instead, the user specifies fault parameters in terms of lateral motion, reverse motion, and fault opening as shown in {numref}`fig:fault:slip:motions`.

:::{figure-md} fig:fault:orientation
<img src="figs/faultOrientation.*" alt="Orientation of a fault surface in 3D, where $\phi$ denotes the angle of the fault strike, $\delta$ denotes the angle of the fault dip, and $\lambda$ the rake angle." width="100px"/>


Orientation of a fault surface in 3D, where <span class="math inline"><em>&#x3D5;</em></span> denotes the angle of the fault strike, <span class="math inline"><em>&#x3B4;</em></span> denotes the angle of the fault dip, and <span class="math inline"><em>&#x3BB;</em></span> the rake angle.
:::

:::{figure-md} fig:fault:slip:motions
<img src="figs/slipmotions.*" alt="Sign conventions associated with fault slip. Positive values are associated with left-lateral, reverse, and fault opening motions." width="100px"/>

Sign conventions associated with fault slip. Positive values are associated with left-lateral, reverse, and fault opening motions.
:::

## Fault Implementation

In order to create relative motion across the fault surface in the finite-element mesh, additional degrees of freedom are added along with adjustment of the topology of the mesh.
These additional degrees of freedom are associated with cohesive cells.
These zero-volume cells allow control of the relative motion between vertices on the two sides of the fault.
PyLith automatically adds cohesive cells for each fault surface.
{numref}`fig:fault:cohesive:cells` illustrates the results of inserting cohesive cells in a mesh consisting of triangular cells.
This example also shows the distinction between how buried fault edges are handled differently than fault edges that reach the edge of the domain, such as the ground surface.

:::{figure-md} fig:fault:cohesive:cells
<img src="figs/cohesivecell.*" alt="Example of cohesive cells inserted into a mesh of triangular cells. The zero thickness cohesive cells control slip on the fault via the relative motion between the vertices on the positive and negative sides of the fault." width="100px" />

Example of cohesive cells inserted into a mesh of triangular cells. The zero thickness cohesive cells control slip on the fault via the relative motion between the vertices on the positive and negative sides of the fault.
:::

:::{figure-md} fig:fault:fault_edge
<img src="figs/faultEdge.*" alt="Example of how faults with buried edges must be described with two sets of vertices. All of the vertices on the fault are included in the `fault` group; the subset of vertices along the buried edges are included in the `fault_edge` group. In 2-D the fault edges are just a single vertex as shown in {numref}`fig:fault:cohesive:cells`." width="100px"/>

Example of how faults with buried edges must be described with two sets of vertices. All of the vertices on the fault are included in the `fault` group; the subset of vertices along the buried edges are included in the `fault_edge` group. In 2-D the fault edges are just a single vertex as shown in {numref}`fig:fault:cohesive:cells`
:::

For faults that have buried edges, splitting the mesh apart and inserting the cohesive cells becomes complex at the buried edges due to the ambiguity of defining where the fault ends and how to insert the cohesive cell.
In PyLith v2.0.0 we have changed how the buried edges of the fault are managed.
An additional group of fault nodes is specified (e.g., via a nodeset from CUBIT) that marks the buried edges of the fault (see {numref}`fig:fault:fault_edge`).
This allows the cohesive cell insertion algorithm to adjust the topology so that cohseive cells are inserted up to the buried edge of the fault but no additional degrees of freedom are added on the fault edge.
This naturally forces slip to zero along the buried edges.

## Fault Parameters

The principal parameters for fault interface conditions are:

:id: This is an integer identifier for the fault surface. It is used to specify the **material-id** of the cohesive cells in the mesh. Material identifiers must be unique across al materials and fault interfaces. Because PyLith creates the cohesive cells at runtime, there is no correspondence between the **id** property and information in the input mesh like there is for materials.
:label: Name of group of vertices associated with the fault surface. This label is also used in error and diagnostic reports (default="").
:edge: Name of group of vertices marking the buried edges of the fault (default="").
:ref_dir_1: First choice for reference direction to discriminate among tangential directions in 3-D (default = [0,0,1]);
:ref_dir_2: Second choice for reference direction to discriminate among tangential directions in 3-D (default=[0,1,0]);
:observers: Observers of boundary condition, e.g., output (default=*[PhysicsObserver]*):

In 2D the default in-plane slip is left-lateral, so we use the reference directions to resolve the ambiguity in specifying reverse slip.
In 3D the reference directions are used to resolve the ambiguity in the along-strike and dip-dir directions.
If the fault plane is horizontal, then the up-dir corresponds to the reverse-motion on the +z side of the fault.
The only requirement for this direction is that it not be collinear with the fault normal direction.
The default value of [0, 0, 1] is appropriate for most 3D problems.

By default the output observers write both diagnostic information (e.g., fault orientation directions) and the slip at each time step.
The fault coordinate system is shown in {numref}`fig:fault:slip:motions`.
The vectors in the fault coordinate system can be transformed to the global coordinate system using the direction vectors in the diagnostic output.

:::{warning}
Output of fault orientation information is not yet available in the beta release.
:::

:::{admonition} TODO
:class: error

Implement output of fault orientation information.
:::

## Kinematic Earthquake Rupture (*FaultCohesiveKin*)

Kinematic earthquake ruptures use the *FaultCohesiveKin* object to prescribe the slip as a function of time on the fault surface.
Slip may evolve simultaneously over the fault surface instantaneously in a single time step (as is usually done in quasistatic simulations) or propagate over the fault surface over hundreds and up to thousands of time steps (as is usually done in a dynamic simulation).

Multiple earthquake ruptures can be specified on a single fault surface.
This permits repeatedly rupturing the same portion of a fault or combining earthquake rupture on one subset of the fault surface with steady aseismic slip on another subset (the two subsets may overlap in both time and space).
A dynamic array of kinematic earthquake rupture components associates a name (string) with each kinematic rupture.
The default dynamic array contains a single earthquake rupture, &ldquo;rupture.&rdquo;

:::{admonition} TODO
:class: error

Add note about discretization of "slip" auxiliary subfield.
:::

### Kinematic Rupture Parameters (*KinSrc*)

The kinematic rupture parameters include the origin time and slip time function.
The slip initiation time in the slip time function is relative to the origin time (default is 0).
This means that slip initiates at a point at a time corresponding to the sum of the kinematic rupture&rsquo;s origin time and the slip initiation time for that point.

```{code-block} cfg
---
caption: FaultCohesiveKin parameters in a `cfg` file
---
[pylithapp.problem]
interfaces = [fault]
# FaultCohesiveKin is the default for a fault.

[pylithapp.problem.interfaces.fault]
id = 10
label = fault_x
edge = fault_x_buried_edge

# KinSrcStep is the default slip-time function.

[pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.label = Fault rupture auxiliary field spatial
database db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
db_auxiliary_field.data = [0.0*s, -2.0*m, 0.0*m]
```

### Slip Time Function

The current release of PyLith supports specification of the evolution of fault slip using analytical expressions for the slip time history at each point, where the parameters for the slip time function may vary over the fault surface.
Currently, two slip time functions are available: (1) a step-function for quasistatic modeling of earthquake rupture and (2) a constant slip rate time function for modeling steady aseismic slip.
Additional slip time functions will likely be available in future releases.
The default slip time function is the step-function slip function.

#### Step-Function Slip Time Function (*KinSrcStep*)

This slip function prescribes a step in slip at a given time at a point:
%
\begin{gather}
D(t)=\left\{ \begin{array}{cc}
0 & 0\leq t<t_{r}\\
D_{final} & t\ge t_{r}
\end{array}\right.\,,
\end{gather}
%
where $D(t)$ is slip at time $t$, $D_{final}$ is the final slip, and $t_{r}$ is the slip initiation time (time when rupture reaches the location).
The slip is specified independently for each of the components of slip, and the slip and slip starting time may vary over the fault surface.

```{table} Values in the auxiliary field spatial database for *KinSrcStep*.
:name: tab:slip:function:step
|      Subfield   |           Components           |
|:----------------|:-------------------------------|
| initiation_time |     --                         |
| final_slip      | opening, left_lateral, reverse |
```

#### Constant Slip Rate Slip Time Function (*KinSrcConstRate*)

This slip function prescribes a constant slip rate for the evolution of slip at a point:
%
\begin{gather}
  D(t)=\left\{ \begin{array}{cc}
0 & 0\leq t<t_{r}\\
V(t-t_{r}) & t\ge t_{r}
\end{array}\right.\,,
\end{gather}
%
where $D(t)$ is slip at time $t$, $V$ is the slip rate, and $t_{r}$ is the slip initiation time (time when rupture reaches the location).
The slip rate is specified independently for each of the components of slip, and the slip rate and slip starting time may vary over the fault surface.

```{table} Values in the auxiliary field spatial database for *KinSrcConstRate*.
:name: tab:slip:function:step
|  Subfield       |  Components                    |
|:----------------|:-------------------------------|
| initiation_time |            --                  |
| slip_rate       | opening, left_lateral, reverse |
```
