# PETSc TS Notes

## Explicit Time Stepping

Explicit time stepping with the PETSc TS requires {math}`F(t,s,\dot{s}) = \dot{s}`.

* We do not specify the functions {math}`\vec{f}_0(t,s,\dot{s})` and {math}`\mathbf{f}_1(t,s,\dot{s})` because the PETSc TS will assume {math}`F(t,s,\dot{s}) = \dot{s}` if no LHS (or I) function is given.
* The PETSc TS will verify that the LHS (or I) function is null.
* We also do not specify {math}`J_F` or {math}`J_G`.
* This leaves us with only needing to specify {math}`\vec{g}_0(t,s)` and {math}`\mathbf{g}_1(t,s)`.

For explicit time stepping with the PETSc TS, we need {math}`F(t,s,\dot{s}) = \dot{s}`.
Using a finite-element formulation for elastodynamics, {math}`F(t,s,\dot{s})` generally involves integrals of the inertia over the domain.
It is tempting to simply move these terms to the RHS, but that introduces inertial terms into the boundary conditions, which makes them less intuitive.
Instead, we transform our equation into the form {math}`\dot{s} = G^*(t,s)` where {math}`G^*(t,s) = M^{-1} G(t,s)`.
We take {math}`M` to be a lumped (diagonal) matrix, so that {math}`M^{-1}` is trivial to compute.
In computing the RHS function, {math}`G^*(t,s)`, we compute {math}`G(t,s)`, then compute {math}`M` and {math}`M^{-1}`, and then {math}`M^{-1}G(t,s)`.
For the elasticity equation with inertial terms, {math}`M` contains the mass matrix.

## Implicit Time Stepping

The LHS (or I) function is associated with implicit time-stepping. When using implicit time-stepping, we place all of the terms on the LHS. Even though placing all of the terms on the LHS sometimes requires different pointwise functions for implicit and explicit time stepping, it minimizes the number of pointwise functions needed for implicit time stepping. If no RHS function is given, then the PETSc TS assumes {math}`G(t,s) = 0`, so we only need to specify {math}`F(t,s,\dot{s})` and {math}`J_f`.

## Implicit-Explicit Time Stepping

For implicit-explicit time stepping algorithms, the equations integrated with explicit time stepping have {math}`\dot{s}` as the LHS function, and the equations integrated with implicit time stepping have 0 as the RHS function.
