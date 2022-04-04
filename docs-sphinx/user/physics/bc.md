# Boundary Conditions

## Assigning Boundary Conditions

There are three basic steps in assigning a specific boundary condition to a portion of the domain.

1.  Create sets of vertices in the mesh generation process for each boundary condition.

2.  Set the parameters for each boundary condition group using `cfg` files and/or command line arguments.

3.  Specify the spatial variation in parameters for the boundary condition using a spatial database file.

## Creating Sets of Vertices

The procedure for creating sets of vertices differs depending on the mesh generator.
For meshes specified using the PyLith mesh ASCII format, the sets of vertices are specified using groups (see {ref}`sec-user-file-formats-meshio-ascii`.
In CUBIT/Trelis the groups of vertices are created using nodesets.
Similarly, in LaGriT, psets are used.
Note that we chose to associate boundary conditions with groups of vertices because nearly every mesh generation package supports associating a string or integer with groups of vertices.
Note also that we currently associate boundary conditions with string identifiers, so even if the mesh generator uses integers, the name is specified as the digits of the integer value.
Finally, note that every vertex set that ultimately is associated with a boundary condition on a cell face (e.g., Neumann boundary conditions and fault interface conditions) must correspond to a simply-connected surface.

## Arrays of Boundary Condition Components

A dynamic array of boundary condition components associates a name (string) with each boundary condition. The default boundary condition for each component in the array is the *DirichletTimeDependent* object.
Other boundary conditions can be bound to the named items in the array via a file or the command line.
The parameters for the boundary condition are set using the name of the boundary condition.

```{code-block} cfg
---
caption: Array of boundary conditions in a `cfg` file
---
[pylithapp.problem]
# Array of four boundary conditions
bc = [x_neg, x_pos, y_pos, z_neg]
# Default boundary condition is DirichletBC
# Keep default value for x_neg and x_pos
bc.y_pos = pylith.bc.AbsorbingDampers
bc.z_neg = pylith.bc.NeumannTimeDependent
```

# Time-Dependent Boundary Conditions

Several boundary conditions use a common formulation for the spatial and temporal variation of the boundary condition parameters,
%
\begin{equation}
f(\vec{x})=f_{0}(\vec{x})+\dot{f}_{1}(\vec{x})(t-t_{1}(\vec{x}))+f_{2}(\vec{x})a(t-t_{2}(\vec{x})),\
\end{equation}
%
where $f(\vec{x})$ may be a scalar or vector parameter, $f_{0}(\vec{x})$ is a constant value, $\dot{f}_{1}(\vec{x})$ is a constant rate of change in the value, $t_{1}(\vec{x})$ is the onset time for the constant rate of change, $f_{2}(\vec{x})$ is the amplitude for the temporal modulation, $a(t)$ is the variation in amplitude with time, $t_{2}(\vec{x})$ is the onset time for the temporal modulation, and $\vec{x}$ is the position of a location in space.
This common formulation permits easy specification of a scalar or vector with a constant value, constant rate of change of a value, and/or modulation of a value in time.
One can specify just the initial value, just the rate of change of the value (along with the corresponding onset time), or just the modulation in amplitude (along with the corresponding temporal variation and onset time), or any combination of the three.

## Time-Dependent Dirichlet Boundary Conditions (`DirichletTimeDependent`)

Dirichlet boundary conditions in PyLith prescribe the a solution subfield on a subset of the vertices of the finite-element mesh.
Currently, these constraints are required to be associated with vertices on a simply-connected boundary surface.

:::{seealso}
[`DirichletTimeDependent` Component](../components/bc/DirichletTimeDependent.md)
:::

## Neumann Boundary Conditions

Neumann boundary conditions are surface tractions applied over a boundary.
As with the *DirichletTimeDependent* condition, each Neumann boundary condition can only be applied to a simply-connected surface.
The surface over which the tractions are applied always has a spatial dimension that is one less than the dimension of the finite-element mesh.

:::{seealso}
[`Neumann` Component](../components/bc/Neumann.md)
:::

### Neumann Boundary Condition Spatial Database Files

The spatial database file the auxiliary subfields for the Neumann boundary condition specify the parameters for the time-dependent expressions.

```{table} Values in the auxiliary field spatial database used for Neumman boundary conditions.
| Dimension    |    Flag |                       Required Values                                              |
|:----|:----|:----------------------------------------------------------------------------------------------------------------------------|
| 2   |  **use_initial**   | initial_amplitude_normal, initial_amplitude_tangential                                                                      |
|     |  **use_rate**   | rate_start_time, rate_amplitude_normal, rate_amplitude_tangential                                                           |
|     |  **use_time_history**   | time_history_start, time_history_amplitude_normal, time_history_amplitude_tangential                                        |
| 3   |  **use_initial**   | initial_amplitude_normal, initial_amplitude_tangential_1, initial_amplitude_tangential_2                                    |
|     | **use_rate**     | rate_start_time, rate_amplitude_normal, rate_amplitude_tangential_1, rate_amplitude_tangential_2                            |
|     |  **use_time_history**    | time_history_start, time_history_amplitude_normal, time_history_amplitude_tangential_1, time_history_amplitude_tangential_2 |
```

# Absorbing Boundary Conditions (`AbsorbingDampers`)

This *AbsorbingDampers* boundary condition attempts to prevent seismic waves reflecting off of a boundary by placing simple dashpots on the boundary.
Normally incident dilatational and shear waves are perfectly absorbed.
Waves incident at other angles are only partially absorbed.
This boundary condition is simpler than a perfectly matched layer (PML) boundary condition but does not perform quite as well, especially for surface waves.
If the waves arriving at the absorbing boundary are relatively small in amplitude compared to the amplitudes of primary interest, this boundary condition gives reasonable results.

:::{seealso}
[`AbsorbingDampers` Component](../components/bc/AbsorbingDampers.md)
:::

## Finite-Element Implementation of Absorbing Boundary

:::{admonition} TODO
:class: error

Move this to the multiphysics implementation section.
:::

Consider a plane wave propagating at a velocity $c$.
We can write the displacement field as
%
\begin{equation}
\vec{u}(\vec{x},t)=\vec{u^{t}}(t-\frac{\vec{x}}{c}),
\end{equation}
%
where $\vec{x}$ is position, $t$ is time, and $\vec{u^{t}}$ is the shape of the propagating wave.
For an absorbing boundary we want the traction on the boundary to be equal to the traction associated with the wave propagating out of the domain.
Starting with the expression for the traction on a boundary, $T_{i}=\sigma_{ij}n_{j},$ and using the local coordinate system for the boundary $s_{h}s_{v}n,$ where $\vec{n}$ is the direction normal to the boundary, $\overrightarrow{s}_{h}$ is the horizontal direction tangent to the boundary, and $\overrightarrow{s}_{v}$ is the vertical direction tangent to the boundary, the tractions on the boundary are
%
\begin{gather}
T_{s_{h}}=\sigma_{s_{h}n}\\
T_{s_{v}}=\sigma_{s_{v}n}\\
T_{n}=\sigma_{nn}.
\end{gather}

In the case of a horizontal boundary, we can define an auxiliary direction in order to assign unique tangential directions.
For a linear elastic isotropic material, $\sigma_{ij}=\lambda\epsilon_{kk}\delta_{ij}+2\mu\epsilon_{ij},$ and we can write the tractions as
%
\begin{gather}
T_{s_{h}}=2\mu\epsilon_{s_{h}n}\\
T_{s_{v}}=2\epsilon_{s_{v}n}\\
T_{n}=(\lambda+2\mu)\epsilon_{nn}+\lambda(\epsilon_{s_{h}s_{h}}+\epsilon_{s_{v}s_{v}}).
\end{gather}
%
For infinitesimal strains, $\epsilon_{ij}=\frac{1}{2}(u_{i,j}+u_{j,i})$ and we have
%
\begin{gather}
\epsilon_{s_{h}n}=\frac{1}{2}(u_{s_{h},n}+u_{n,s_{h}})\\
\epsilon_{s_{v}n}=\frac{1}{2}(u_{s_{v},n}+u_{n,s_{v}})\\
\epsilon_{nn}=u_{n,n}.
\end{gather}
%
For our propagating plane wave, we recognize that
%
\begin{equation}
\frac{\partial\vec{u^{t}}(t-\frac{\vec{x}}{c})}{\partial x_{i}}=-\frac{1}{c}\frac{\partial\vec{u^{t}}(t-\frac{\vec{x}}{c})}{\partial t},
\end{equation}
%
so that our expressions for the tractions become
%
\begin{gather}
T_{s_{h}}=-\frac{\mu}{c}\left(\frac{\partial u_{s_{h}}^{t}(t-\frac{\vec{x}}{c})}{\partial t}+\frac{\partial u_{n}^{t}(t-\frac{\vec{x}}{c})}{\partial t}\right),\\
T_{s_{v}}=-\frac{\mu}{c}\left(\frac{\partial u_{s_{v}}^{t}(t-\frac{\vec{x}}{c})}{\partial t}+\frac{\partial u_{n}^{t}(t-\frac{\vec{x}}{c})}{\partial t}\right).
\end{gather}
%
For the normal traction, consider a dilatational wave propagating normal to the boundary at speed $v_{p}$; in this case $u_{s_{h}}=u_{s_{v}}=0$ and $c=v_{p}$.
For the shear tractions, consider a shear wave propagating normal to the boundary at speed $v_{s}$; we can decompose this into one case where $u_{n}=u_{s_{v}}=0$ and another case where $u_{n}=u_{s_{h}}=0$, with $c=v_{s}$ in both cases.
We also recognize that $\mu=\rho v_{s}^{2}$ and $\lambda+2\mu=\rho v_{p}^{2}$.
This leads to the following expressions for the tractions:
%
\begin{gather}
T_{s_{h}}=-\rho v_{s}\frac{\partial u_{s_{h}}^{t}(t-\frac{\vec{x}}{c})}{\partial t}\\
T_{s_{v}}=-\rho v_{s}\frac{\partial u_{v}^{t}(t-\frac{\vec{x}}{c})}{\partial t}\\
T_{n}=-\rho v_{p}\frac{\partial u_{n}^{t}(t-\frac{\vec{x}}{c})}{\partial t}
\end{gather}
%
We write the weak form of the boundary condition as
%
\begin{equation}
\int_{S_{T}}T_{i}\phi_{i}\, dS=\int_{S_{T}}-\rho c_{i}\frac{\partial u_{i}}{\partial t}\phi_{i}\, dS,
\end{equation}
%
where $c_{i}$ equals $v_{p}$ for the normal traction and $v_{s}$ for the shear tractions, and $\phi_{i}$ is our weighting function.
We express the trial solution and weighting function as linear combinations of basis functions,
%
\begin{gather}
u_{i}=\sum_{m}a_{i}^{m}N^{m},\\
\phi_{i}=\sum_{n}c_{i}^{n}N^{n}.
\end{gather}
%
Substituting into our integral over the absorbing boundaries yields
%
\begin{equation}
\int_{S_{T}}T_{i}\phi_{i}\, dS=\int_{S_{T}}-\rho c_{i}\sum_{m}\dot{a}_{i}^{m}N^{m}\sum_{n}c_{i}^{n}N^{n}\, dS.
\end{equation}
%
In the derivation of the governing equations, we recognized that the weighting function is arbitrary, so we form the residual by setting the terms associated with the coefficients $c_{i}^{n}$ to zero,
%
\begin{equation}
r_{i}^{n}=\sum_{\text{tract cells}}\sum_{\text{quad pts}}-\rho(x_{q})c_{i}(x_{q})\sum_{m}\dot{a}_{i}^{m}N^{m}(x_{q})N^{n}(x_{q})w_{q}|J_{cell}(x_{q})|,
\end{equation}
%
where $x_{q}$ are the coordinates of the quadrature points, $w_{q}$ are the weights of the quadrature points, and $|J_{cell}(x_{q})|$ is the determinant of the Jacobian matrix evaluated at the quadrature points associated with mapping the reference cell to the actual cell.

The appearance of velocity in the expression for the residual means that the absorbing dampers also contribute to the system Jacobian matrix.
Using the central difference method, the velocity is written in terms of the displacements,
%
\begin{equation}
\dot{u}_{i}(t)=\frac{1}{2\Delta t}(u_{i}(t+\Delta t)-u_{i}(t-\Delta t)).
\end{equation}
%
Expressing the displacement at time $t+\Delta t$ in terms of the displacement at time $t$ ($u_{i}(t)$) and the increment in the displacement at time $t$ ($du_{i}(t)$) leads to
%
\begin{equation}
\dot{u}_{i}(t)=\frac{1}{2\Delta t}(du_{i}(t)+u_{i}(t)-u_{i}(t-\Delta t))
\end{equation}
%
The terms contributing to the system Jacobian are associated with the increment in the displacement at time $t$.
Substituting into the governing equations and isolating the term associated with the increment in the displacement at time t yields
%
\begin{equation}
A_{ij}^{nm}=\sum_{\text{tract cells}}\sum_{\text{quad pts}}\delta_{ij}\frac{1}{2\Delta t}\rho(x_{q})v_{i}(x_{q})N^{m}(x_{q})N^{n}(x_{q})w_{q}|J_{cells}(x_{q})|,
\end{equation}
%
where $A_{ij}^{mn}$ is an $nd$ by $md$ matrix ($d$ is the dimension of the vector space), $m$ and $n$ refer to the basis functions and $i$ and $j$ are vector space components.
