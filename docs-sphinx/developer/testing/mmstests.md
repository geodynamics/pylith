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

In mathematical terms, the residual test is that the residual computed for a known solution is below some tolerance, $\epsilon$,
\begin{equation}
  || F(\vec{s}) - G(\vec{s}) || \le \epsilon,
\end{equation}
where $F(\vec{s})$ is the LHS residual and $G(\vec{s})$ is the RHS residual.
In the Taylor series Jacobian test, we verify that
\begin{equation}
  || F(\vec{s} + \epsilon \vec{\delta s}) - F(\vec{s}) - \epsilon J \vec{v} || < \epsilon^2,
\end{equation}
where $\vec{\delta s}$ is a perturbation in the soluton and $J$ is the Jacobian matrix.

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

## Debugging residual errors

When the residual test fails, we generally use the following procedure to diagnose the problem.

1. Verify that the discretization check passes indicating accurate representation of the solution in the finite-element space.
2. Run just the residual test for a single discretization and turn on the debug journal corresponding to the name of the MMS test as set by `GenericComponent::setName(JOURNAL_NAME)` in `setUp()`; this is done via the `--journal.debug=JOURNAL_NAME` command line argument to the MMS test driver.
3. Verify that the residual kernels show up correctly in the view of the PETSc discretization.
4. Analyze the residual vector to see which degrees of freedom have nonzero terms. Look at the solution section to see what solution subfield and point are associated with those degrees for freedom.
5. Check the pointwise functions for the residual and solution associated with the subfield with nonzero terms.

## Debugging Jacobian errors

When one of the Jacobian tests fails, we focus on the finite-difference Jacobian test.

1. Verify that the residual check passes.
2. Run just the Jacobian finite-difference test for a single discretization and turn on the debug journal corresponding to the name of the MMS test as set by `GenericComponent::setName(JOURNAL_NAME)` in `setUp()`; this is done via the `--journal.debug=JOURNAL_NAME` command line argument to the MMS test driver.
3. Examine any differences between the hand-coded Jacobian and the finite-difference Jacobian.
4. Check the differences against pointwise functions listed in the view of the PETSc discretization. Are any expected functions missing from the list?
5. Isolate the error by dropping corresponding terms from the residual, Jacobian, and solution by editing the pointwise functions until the hand-coded and finite-difference Jacobians agree; that is, drop terms from the governing equation until the MMS tests pass. Then, add the terms back into the governing equation (residual, Jacobian, and solution pointwise functions) one by one, fixing errors as they are detected.

## Example

In `tests/mmstests/elasticity`, we create a suite of MMS tests that all use a single material (`Elasticity`), a single Dirichlet boundary condition, and a data object to hold test-specific parameters.
The `pylith::mmstests::TestElasticity` object holds these parameters with the `pylith::mmstests::TestIsotropicLinearElasticity` object adding a rheology.
The individual test cases, such as `pylith::mmstests::TestIsotropicLinearElasticity2D_UniformStrain` provides the problem-specific parameters, such as analytical functions for the solution and auxiliary fields and sets up the problem.
Child classes specify different cell shapes and orders for the basis functions and numerical quadrature.

:::{warning}
The quadrature order must be the same across all solution subfields and auxiliary subfields.
:::

:::{admonition} TODO
:class: error

Add a simple test case that introduces known errors into the discretization, residual, and Jacobian to illustrate these steps.
:::
