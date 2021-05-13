# Method of Manufactured Solutions

We use the Method of Manufactured Solutions (MMS) to verify the order of accuracy of our solution to the governing equations.
See [Code Verification by the Method of Manufactured Solutions](https://www.osti.gov/servlets/purl/759450) for a more detailed discussion of MMS applied to code verification for simulations of partial differential equations.
The general formulation is to assume a solution for the governing equations.
If the solution is not an exact solution, then plugging it into the governing equation will result in a residual that is a body force.
This body force can be included as an additional term in the governing equation so that the solution becomes exact.

In our application of the MMS, we test four features:

1. Representation of the solution in the finite-element space;
2. Residual for the governing equation is zero within some tolerance;
3. A Taylor series expansion of the Jacobian provides the expected order of accuracy; and
4. A finite-difference Jacobian matches the computed Jacobian within some tolerance.

The C++ class `pylith::testing::MMSTest` implements these four tests.
In implementing an MMS test case, we need to provide:

* Analytical functions for each subfield of the solution;
* Analytical functions for each subfield in the auxiliary fields;
* Dirichlet boundary conditions;
* Finite-element mesh;
* Setup of the problem, including the solution, physics, boundary conditions, and the PETSc TS.

:::{warning}
In setting up a MMS test, remember that the solution should be able to be represented in the finite-element space.
If you use a polynomial of order 2 for a subfield of the solution, then you should use a basis order of at least 2 for that subfield of the solution.
:::

## Example

In `tests/mmstests/elasticity` we create a suite of MMS tests that all use a single material (`Elasticity`), a single Dirichlet boundary condition, and a data object to hold test-specific parameters.
The `pylith::mmstests::TestElasticity` object holds these parameters with the `pylith::mmstests::TestIsotropicLinearElasticity` object adding a rheology.
The individual test cases, such as `pylith::mmstests::TestIsotropicLinearElasticity2D_UniformStrain` provides the problem-specific parameters, such as analytical functions for the solution and auxiliary fields and sets up the problem.
Child classes specify different cell shapes and orders for the basis functions and numerical quadrature.

:::{warning}
The quadrature order must be the same across all solution subfields and auxiliary subfields.
:::
