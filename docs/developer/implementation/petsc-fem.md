(sec-developer-petsc-fem)=
# PETSc Finite-Element Implementation

Formulating the weak form of the governing equation in terms of pointwise functions allows the PyLith implementation of the equations to be done at a rather high level.
Most of the finite-element details are encapsulated in PETSc routines that compute the integrals and solve the system of equations.
In fact, adding materials and boundary conditions requires calling only a few PETSc finite-element functions to register the point-wise functions.
A new material may need to add to the library of solution subfields and auxiliary subfields, but adding these fields is also done at a high-level.

The remainder of this section discusses three aspects of the finite-element implementation handled by PyLith to give you a peek of what is going on under the hood.
Refer to the {ref}`cha:glossary` for definitions of the terminology.

## DMPlex

The finite-element mesh is stored as a `DMPlex` object.
This is a particular implementation of the PETSc Data Management (`DM`, called `PetscDM` within PyLith) object.
Within a `DMPlex` object, vertices, edges, faces, and cells are all called points.
The points are numbered sequentially, beginning with cells, followed by vertices, edges, and then faces.
Treating all topological pieces of the mesh the same way, as points in an abstract graph, allows us to write algorithms which apply to very different meshes without change.
For example, we can write a finite element assembly loop that applies to meshes of any dimension, with any cell shape, with (almost) any finite element.

### Point Depth and Height

In general, vertices are at a _depth_ of 0 and cells are at the maximum depth.
Similarly, cells are at a _height_ of 0 and vertices are at the maximum height.
Notice that depth and height correspond to the usual dimension and codimension.
{numref}`tab-developer-plex-depth-height` shows the heights and depths of the vertices edges, faces, and cells in a 3D mesh.


```{table} Depth and height of various topological pieces in DMPlex.
:name: tab-developer-plex-depth-height
| Point type | Depth | Height |
| :--------- | :---: | :----: |
| Vertices   |   0   |   3    |
| Edges      |   1   |   2    |
| Faces      |   2   |   1    |
| Cells      |   3   |   0    |
```

For a boundary mesh, we currently store the full set of points (vertices, edges, faces, and cells).
Obviously for 2-D meshes, the boundary mesh doesn't contain "volume" cells, but just vertices, edges, and faces.
This means the "boundary" cells are at a height of 1 and a depth equal to the maximum depth of 1.

:::{figure-md} fig-developer-mesh-topology
<img src="figs/meshtopology.*" alt="Mesh topology" width="100%"/>

Conventional numbering (left) with vertices and cells numbered independently and `DMPlex` numbering (right) with cells, vertices, edges, (and faces), numbered sequentially.
:::

(sec-developer-pylith-field)=
## `pylith::topology::Field`, `PetscSection`, and `PetscVec`

A PyLith `Field` (`pylith::topology::Field`) is composed of a `PetscSection`, which associates the points (vertices, edges, faces, and cells) with the finite-element coefficients, and a `PetscVec`, which is a vector storing the finite-element coefficients.
A `Field` may hold a single subfield, such as displacement, or it may hold several subfields, such as the density, shear modulus, and bulk modulus for an isotropic, linear elastic material.

The chart of the `PetscSection` defines the range of points (minimum and maximum) over which the section is defined.
For each point in the chart, the section holds the number of degrees of freedom a point has and the offset in the `PetscVec` for its first degree of freedom.
The section also holds the number of degrees of freedom and offset for each individual subfield within the section for each point in the chart.

Because the `PetscSection` knows only about points, and not about topology or geometry, it allows us to express the mathematical concepts without getting lost in the details.
For example, a $P_1$ 3D vector field assigns 3 degrees of freedom to each vertex.
This same section could be used to layout a field over a circle, surface of a cylinder, Mobius strip, sphere, cube, ball, or PyLith mesh.
The section separates the layout of data over a mesh from the actual storage of values, which is done by the `PetscVec`.
The section tells you how to look up the values, which are associated with a piece of the mesh, inside the vector.
For example, you can ask for the offset into the vector for the values associated with an edge in the mesh, and how many degrees of freedom there are on the edge.

The `PetscSection` also includes the information about the constrained degrees of freedom.
We refer to a `PetscVec` that includes the values for the constrained degrees of freedom as a _local_ vector and a `PetscVec` with the values for the constrained degrees of freedom removed as a _global_ vector.
Constraints often arise from Dirichlet boundary conditions, which change the basis functions of the approximating space, but can also arise from algebraic relations between unknowns.
A local vector is used for assembly of the residual vector and Jacobian matrix, because we need the boundary values in order to compute those integrals.
Global vectors are used for the algebraic solver because we do not want solution values fixed by constraints to participate in the solve.

The `Field::view()` method will print the `Field` metadata and the `PetscSection` in a more useful format compared to the corresponding PETSc viewer functions.

% ./mmstest_linearelasticity_faults2d OneFaultShearNoSlip::QuadQ1::testResidual --journal.debug=timedependent
```{code-block} bash
---
caption: Output from `pylith::topology::Field::view()` for a solution field with two subfields showing only the sections (layout of the field). The displacement subfield has degrees of freedom on the vertices of the bulk cells. The lagrange_multiplier_fault subfield has degrees of freedom on the edges of the cohesive cells. The order of the values (offsets in the PetscVec) follows the ordering of the points (cells, vertices, edges, and faces). The local section includes the constrained degrees of freedom, whereas the global section does not. The sections list the point and then the degree of freedom associated with the components of the fields.
---
Viewing field 'DM_0x84000000_1 solution' Solution field.
  Subfields:
    Subfield displacement, index: 0, components: displacement_x displacement_y, scale: 1000, basisOrder: 1, quadOrder: 1, dimension: -1, cellBasis: tensor
    Subfield lagrange_multiplier_fault, index: 1, components: lagrange_multiplier_fault_x lagrange_multiplier_fault_y, scale: 2.25e+10, basisOrder: 1, quadOrder: 1, dimension: 1, cellBasis: tensor
        displacement  lagrange_multiplier_fault
Local Section
Processor 0
     0    --    --      --    -- # cells
     1    --    --      --    --
     2    --    --      --    --
     3    --    --      --    --
     4    --    --      --    --
     5    --    --      --    --
     6    --    --      --    --
     7    --    --      --    --
     8    --    --      --    --
     9    --    --      --    --
    10    --    --      --    --
    11    --    --      --    --
    12    24    25      --    --    (constrained: 24 25) # vertices
    13    26    27      --    --    (constrained: 26 27)
    14    28    29      --    --    (constrained: 28 29)
    15    30    31      --    --    (constrained: 30 31)
    16    32    33      --    --
    17    34    35      --    --
    18    36    37      --    --
    19    38    39      --    --
    20    40    41      --    --    (constrained: 40 41)
    21    42    43      --    --    (constrained: 42 43)
    22    44    45      --    --    (constrained: 44 45)
    23    46    47      --    --    (constrained: 46 47)
    24     0     1      --    --
    25     2     3      --    --
    26     6     7      --    --
    27     8     9      --    --
    28    12    13      --    --
    29    14    15      --    --
    30    18    19      --    --
    31    20    21      --    --
    32    --    --      --    -- # edges
    33    --    --      --    --
    34    --    --      --    --
    35    --    --      --    --
    36    --    --      --    --
    37    --    --      --    --
    38    --    --      --    --
    39    --    --      --    --
    40    --    --      --    --
    41    --    --      --    --
    42    --    --      --    --
    43    --    --      --    --
    44    --    --      --    --
    45    --    --      --    --
    46    --    --      --    --
    47    --    --      --    --
    48    --    --      --    --
    49    --    --      --    --
    50    --    --      --    --
    51    --    --      --    --
    52    --    --      --    --
    53    --    --      --    --
    54    --    --      --    --
    55    --    --      --    --
    56    --    --      --    --
    57    --    --      --    --
    58    --    --      --    --
    59    --    --       4     5 # cohesive edges
    60    --    --      10    11
    61    --    --      16    17
    62    --    --      22    23
Global Section
Processor 0
     0    --    --      --    -- # cells
     1    --    --      --    --
     2    --    --      --    --
     3    --    --      --    --
     4    --    --      --    --
     5    --    --      --    --
     6    --    --      --    --
     7    --    --      --    --
     8    --    --      --    --
     9    --    --      --    --
    10    --    --      --    --
    11    --    --      --    --
    12    --    --      --    -- # vertices
    13    --    --      --    --
    14    --    --      --    --
    15    --    --      --    --
    16    24    25      --    --
    17    26    27      --    --
    18    28    29      --    --
    19    30    31      --    --
    20    --    --      --    --
    21    --    --      --    --
    22    --    --      --    --
    23    --    --      --    --
    24     0     1      --    --
    25     2     3      --    --
    26     6     7      --    --
    27     8     9      --    --
    28    12    13      --    --
    29    14    15      --    --
    30    18    19      --    --
    31    20    21      --    --
    32    --    --      --    -- # edges
    33    --    --      --    --
    34    --    --      --    --
    35    --    --      --    --
    36    --    --      --    --
    37    --    --      --    --
    38    --    --      --    --
    39    --    --      --    --
    40    --    --      --    --
    41    --    --      --    --
    42    --    --      --    --
    43    --    --      --    --
    44    --    --      --    --
    45    --    --      --    --
    46    --    --      --    --
    47    --    --      --    --
    48    --    --      --    --
    49    --    --      --    --
    50    --    --      --    --
    51    --    --      --    --
    52    --    --      --    --
    53    --    --      --    --
    54    --    --      --    --
    55    --    --      --    --
    56    --    --      --    --
    57    --    --      --    --
    58    --    --      --    --
    59    --    --       4     5 # cohesive edges
    60    --    --      10    11
    61    --    --      16    17
    62    --    --      22    23
```

```{code-block} bash
---
caption: Corresponding output of PetscSectionVew() for the local section of the solution field.
---
PetscSection Object: 1 MPI process
  type not yet set
2 fields
  field 0 "displacement" with 2 components
Process 0:
  (   0) dof  0 offset  24 # cells
  (   1) dof  0 offset  24
  (   2) dof  0 offset  24
  (   3) dof  0 offset  24
  (   4) dof  0 offset  24
  (   5) dof  0 offset  24
  (   6) dof  0 offset  24
  (   7) dof  0 offset  24
  (   8) dof  0 offset  24
  (   9) dof  0 offset  12
  (  10) dof  0 offset  18
  (  11) dof  0 offset  24
  (  12) dof  2 offset  24 constrained 0 1 # vertices
  (  13) dof  2 offset  26 constrained 0 1
  (  14) dof  2 offset  28 constrained 0 1
  (  15) dof  2 offset  30 constrained 0 1
  (  16) dof  2 offset  32
  (  17) dof  2 offset  34
  (  18) dof  2 offset  36
  (  19) dof  2 offset  38
  (  20) dof  2 offset   0
  (  21) dof  2 offset   6
  (  22) dof  2 offset  12
  (  23) dof  2 offset  18
  (  24) dof  2 offset  40 constrained 0 1
  (  25) dof  2 offset  42 constrained 0 1
  (  26) dof  2 offset  44 constrained 0 1
  (  27) dof  2 offset  46 constrained 0 1
  (  28) dof  2 offset   2
  (  29) dof  2 offset   8
  (  30) dof  2 offset  14
  (  31) dof  2 offset  20
  (  32) dof  0 offset  48 # edges
  (  33) dof  0 offset  48
  (  34) dof  0 offset  48
  (  35) dof  0 offset  48
  (  36) dof  0 offset  48
  (  37) dof  0 offset  48
  (  38) dof  0 offset  48
  (  39) dof  0 offset  48
  (  40) dof  0 offset  48
  (  41) dof  0 offset  48
  (  42) dof  0 offset  48
  (  43) dof  0 offset   0
  (  44) dof  0 offset  48
  (  45) dof  0 offset  12
  (  46) dof  0 offset  48
  (  47) dof  0 offset  18
  (  48) dof  0 offset  48
  (  49) dof  0 offset  48
  (  50) dof  0 offset  48
  (  51) dof  0 offset  48
  (  52) dof  0 offset  48
  (  53) dof  0 offset  48
  (  54) dof  0 offset  48
  (  55) dof  0 offset  48
  (  56) dof  0 offset   0
  (  57) dof  0 offset  12
  (  58) dof  0 offset  18
  (  59) dof  0 offset   4 # cohesive edges
  (  60) dof  0 offset  10
  (  61) dof  0 offset  16
  (  62) dof  0 offset  22
  field 1 "lagrange_multiplier_fault" with 2 components
Process 0:
  (   0) dof  0 offset  24 # cells
  (   1) dof  0 offset  24
  (   2) dof  0 offset  24
  (   3) dof  0 offset  24
  (   4) dof  0 offset  24
  (   5) dof  0 offset  24
  (   6) dof  0 offset  24
  (   7) dof  0 offset  24
  (   8) dof  0 offset  24
  (   9) dof  0 offset  12
  (  10) dof  0 offset  18
  (  11) dof  0 offset  24
  (  12) dof  0 offset  26 # vertices
  (  13) dof  0 offset  28
  (  14) dof  0 offset  30
  (  15) dof  0 offset  32
  (  16) dof  0 offset  34
  (  17) dof  0 offset  36
  (  18) dof  0 offset  38
  (  19) dof  0 offset  40
  (  20) dof  0 offset   2
  (  21) dof  0 offset   8
  (  22) dof  0 offset  14
  (  23) dof  0 offset  20
  (  24) dof  0 offset  42
  (  25) dof  0 offset  44
  (  26) dof  0 offset  46
  (  27) dof  0 offset  48
  (  28) dof  0 offset   4
  (  29) dof  0 offset  10
  (  30) dof  0 offset  16
  (  31) dof  0 offset  22
  (  32) dof  0 offset  48 # edges
  (  33) dof  0 offset  48
  (  34) dof  0 offset  48
  (  35) dof  0 offset  48
  (  36) dof  0 offset  48
  (  37) dof  0 offset  48
  (  38) dof  0 offset  48
  (  39) dof  0 offset  48
  (  40) dof  0 offset  48
  (  41) dof  0 offset  48
  (  42) dof  0 offset  48
  (  43) dof  0 offset   0
  (  44) dof  0 offset  48
  (  45) dof  0 offset  12
  (  46) dof  0 offset  48
  (  47) dof  0 offset  18
  (  48) dof  0 offset  48
  (  49) dof  0 offset  48
  (  50) dof  0 offset  48
  (  51) dof  0 offset  48
  (  52) dof  0 offset  48
  (  53) dof  0 offset  48
  (  54) dof  0 offset  48
  (  55) dof  0 offset  48
  (  56) dof  0 offset   0
  (  57) dof  0 offset  12
  (  58) dof  0 offset  18
  (  59) dof  2 offset   4 # cohesive edges
  (  60) dof  2 offset  10
  (  61) dof  2 offset  16
  (  62) dof  2 offset  22
```

## Integration

Integration involves integrals over the domain (materials), over the boundary of the domain (boundary conditions), or over interior interfaces (fault surface).
These three operations are done by different PETSc functions and integrator objects.

`DMPlexComputeResidual_Internal`
: Compute the contribution to the LHS or RHS residual for a single material.
This function and the functions it calls handle looping over the cells in the material, integrating the weak form for each of the fields, and adding them to the residual.
A more appropriate name would be `DMPlexComputeResidualSingle`, and that may be used in the future.
This function is called from `IntegratorDomain`.

`DMPlexComputeBdResidualSingle`
: Compute the contribution to the LHS or RHS residual for a single boundary condition.
This function and the functions it calls handle looping over faces (3D) or edges (2D) on the boundary, integrating the weak form for each of the fields, and adding them to the residual.
This function is called from `IntegratorBoundary`.

`DMPlexComputeResidual_Hybrid_Internal`
: Compute the contribution to the LHS or RHS residual for a single fault surface.
This function and the functions it calls handle looping over the cohesive cells , integrating the weak form for each of the fields, and adding them to the residual.
This function is called from `IntegratorInterface`.

(sec-developer-petsc-projection)=
## Projection

Input and output often involve projecting fields to/from the finite-element space.
PETSc provides a family of functions for this.
We generally use two of these, one for analytical functions and one for discretized fields.
Projection may be a misleading term here, since we are not referring to the common $L_2$ projection, but rather interpolation of the function by functions in our finite-element space.

Let's start with the simple example of Fourier analysis, which most people have experience with.
If we want the Fourier interpolant $\tilde f$ for a given function $f$, then we need to determine its Fourier coefficients, $f_k$, where
\begin{equation}
  \tilde f = \sum_k f_k e^{i k x}.
\end{equation}

This is straightforward because the basis functions in the Fourier representation are orthogonal,
\begin{equation}
  \int^{2\pi}_0 e^{-i m x} e^{i k x} \, dx = 2\pi \delta_{km}.
\end{equation}
To find the coefficient $f_m$, we just multiply by the conjugate of the basis function and integrate,
\begin{align}
  \int^{2\pi}_0 e^{-i m x} \tilde f \, dx &= \int^{2\pi}_0 e^{-i m x} \sum_k f_k e^{i k x} \, dx, \\
                                  &= \sum_k f_k \int^{2\pi}_0 e^{-i m x} e^{i k x} \, dx, \\
                                  &= \sum_k f_k 2\pi \delta_{km}, \\
  &= 2\pi f_m,
\end{align}
and we have our coefficient $f_m$,
\begin{equation}
  f_m = \frac{1}{2 \pi}\int^{2\pi}_0 e^{-i m x} \tilde f \, dx.
\end{equation}

The finite element basis $\phi_i$ is not orthogonal, so we have an extra step.
We could take the inner product of $f$ with all the basis functions, and then sort out the dependencies by solving a linear system (the mass matrix), which is what happens in $L_2$
projection.
However, suppose we have another basis $\psi_i$ of linear functionals which is _biorthogonal_ to $\phi_i$, meaning
\begin{equation}
  \psi_i(\phi_j) = \delta_{ij}.
\end{equation}
We can easily pick out the coefficient of $\tilde f$ by acting with the corresponding basis functional.
Our interpolant is
\begin{equation}
  \tilde f = \sum_k f_k \phi_k(x).
\end{equation}
Acting on the interpolant with the biorthogonal basis results in
\begin{equation}
  \psi_i(\tilde f) = \psi_i\left( \sum_k f_k \phi_k(x) \right).
\end{equation}
Making use of the fact that the finite-element basis $\phi_i$ is linear, yields
\begin{align}
  \psi_i(\tilde f) &= \sum_k f_k \psi_i\left( \phi_k(x) \right),\\
  \psi_i(\tilde f) &= \sum_k f_k \delta_{ik}, \\
  \psi_i(\tilde f) &= f_i.
\end{align}
We call $\phi_i$ the _primal_ basis, and $\psi_i$ the _dual_ basis.
We note that if $f$ does not lie in our approximation space spanned by $\phi_i$, then interpolation is not equivalent to $L_2$ projection.
This will not usually be important for our purposes.

`DMProjectFunctionLocal`
: Project an analytical function into the given finite-element space.

`DMProjectFieldLocal`
: Project a discretized field in one finite-element space into another finite-element space.

## Pointwise functions (kernels)

The following code blocks show the function prototypes for pointwise functions for the residual and the Jacobian.
We use the same prototype for the boundary and fault interfaces residuals.

```{code-block} c++
---
caption: Function prototype for pointwise functions for integrating the residual over the domain (materials).
---
/** PetscPointFunc prototype.
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field.
 * @param[in] numA Number of registered subfields in auxiliary field.
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] numConstants Number of registered constants.
 * @param[in] constants Array of registered constants.
 * @param[out] f [dim].
 */
 void
 func(const PylithInt dim,
      const PylithInt numS,
      const PylithInt numA,
      const PylithInt sOff[],
      const PylithInt sOff_x[],
      const PylithScalar s[],
      const PylithScalar s_t[],
      const PylithScalar s_x[],
      const PylithInt aOff[],
      const PylithInt aOff_x[],
      const PylithScalar a[],
      const PylithScalar a_t[],
      const PylithScalar a_x[],
      const PylithReal t,
      const PylithScalar x[],
      const PylithInt numConstants,
      const PylithScalar constants[],
      PylithScalar f[]);
```

```{code-block} c++
---
caption: Function prototype for pointwise functions for integrating the Jacobian over the domain (materials).
---
/** PetscPointJac prototype.
 *
 * This is identical to the PetscPointFunc with the addition of the
 * s_tshift argument.
 *
 * @param[in] s_tshift The multiplier for dF/dS_t.
 */
 void
 func(const PylithInt dim,
      const PylithInt numS,
      const PylithInt numA,
      const PylithInt sOff[],
      const PylithInt sOff_x[],
      const PylithScalar s[],
      const PylithScalar s_t[],
      const PylithScalar s_x[],
      const PylithInt aOff[],
      const PylithInt aOff_x[],
      const PylithScalar a[],
      const PylithScalar a_t[],
      const PylithScalar a_x[],
      const PylithReal t,
      const PylithReal s_tshift,
      const PylithScalar x[],
      const PylithInt numConstants,
      const PylithScalar constants[],
      PylithScalar J[]);
```

```{code-block} c++
---
caption: Function prototype for pointwise functions for integrating the residual over domain boundaries and interior interfaces.
---
/** PetscBdPointFunc prototype.
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field.
 * @param[in] numA Number of registered subfields in auxiliary field.
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] numConstants Number of registered constants.
 * @param[in] constants Array of registered constants.
 * @param[out] f0 [dim].
 */
void
func(const PylithInt dim,
     const PylithInt numS,
     const PylithInt numA,
     const PylithInt sOff[],
     const PylithInt sOff_x[],
     const PylithScalar s[],
     const PylithScalar s_t[],
     const PylithScalar s_x[],
     const PylithInt aOff[],
     const PylithInt aOff_x[],
     const PylithScalar a[],
     const PylithScalar a_t[],
     const PylithScalar a_x[],
     const PylithReal t,
     const PylithScalar x[],
     const PylithReal n[],
     const PylithInt numConstants,
     const PylithScalar constants[],
     PylithScalar f0[]);
```

```{code-block} c++
---
caption: Function prototype for pointwise functions for integrating the residual over domain boundaries.
---
/** PetscBdPointJac prototype.
 *
 * This is identical to the PetscBdPointFunc with the addition of the
 * s_tshift argument.
 *
 * @param[in] s_tshift The multiplier for dF/dS_t.
 */
void
func(const PylithInt dim,
     const PylithInt numS,
     const PylithInt numA,
     const PylithInt sOff[],
     const PylithInt sOff_x[],
     const PylithScalar s[],
     const PylithScalar s_t[],
     const PylithScalar s_x[],
     const PylithInt aOff[],
     const PylithInt aOff_x[],
     const PylithScalar a[],
     const PylithScalar a_t[],
     const PylithScalar a_x[],
     const PylithReal t,
     const PylithReal s_tshift,
     const PylithScalar x[],
     const PylithReal n[],
     const PylithInt numConstants,
     const PylithScalar constants[],
     PylithScalar J[]);
```

The integrators handle calling the appropriate functions for setting the kernels for the integration.
The `Physics` objects (materials, boundary conditions, and faults) tell the integrators what kernels to use.
Pointwise functions not used should be set to `NULL`.

```{code-block} c++
---
caption: Data structure used to specify kernels for integration of the residual over the domain (materials).
---
struct ResidualKernels {
    std::string subfield; ///< Name of subfield
    PetscPointFunc r0; ///< f0 (RHS) or g0 (LHS) function.
    PetscPointFunc r1; ///< f1 (RHS) or g1 (LHS) function.
};
```

```{code-block} c++
---
caption: Data structure used to specify kernels for integration of the Jacobian over the domain (materials).
---
struct JacobianKernels {
    std::string subfieldTrial; ///< Name of subfield associated with trial function (row in Jacobian).
    std::string subfieldBasis; ///< Name of subfield associated with basis function (column in Jacobian).
    PetscPointJac j0; ///< J0 function.
    PetscPointJac j1; ///< J1 function.
    PetscPointJac j2; ///< J2 function.
    PetscPointJac j3; ///< J3 function.
};
```

```{code-block} c++
---
caption: Data structure used to specify kernels for integration of the residual over the domain (materials).
---
struct ResidualKernels {
    std::string subfield; ///< Name of subfield
    PetscBdPointFunc r0; ///< f0 (RHS) or g0 (LHS) function.
    PetscBdPointFunc r1; ///< f1 (RHS) or g1 (LHS) function.
};
```

```{code-block} c++
---
caption: Data structure used to specify kernels for integration of the Jacobian over the domain (materials).
---
struct JacobianKernels {
    std::string subfieldTrial; ///< Name of subfield associated with trial function (row in Jacobian).
    std::string subfieldBasis; ///< Name of subfield associated with basis function (column in Jacobian).
    PetscBdPointJac j0; ///< J0 function.
    PetscBdPointJac j1; ///< J1 function.
    PetscBdPointJac j2; ///< J2 function.
    PetscBdPointJac j3; ///< J3 function.
};
```
