# Examples: Shear wave propagating down a bar

This suite of exampes demonstrates some basic concepts of using PyLith
to solve the dynamic elasticity equation. We discretize the bar with either
triangular or quadrilateral cells.

Concepts discussed include:

* Solution of the dynamic elasticity equation
* Mesh from CUBIT/Trelis
* Dirichlet boundary condition
* Absorbing boundary condition
* Isotropic, linear elasticity with uniform material properties

## Step01: Dilatational wave using a time-dependent Dirichlet boundary condition

We generate a dilatational (P) wave using time-dependent Dirichlet
boundary condition on one of of a bar. We place a Dirichlet BC on the
other end so that the wave continues to propagate back and forth along
the bar.

The simulation parameters are in the `pylithapp.cfg` and
`step01_pwave_reflected.cfg` files.

To run the example:
```
pylith step01_pwave_reflected.cfg
```


## Step02: Shear wave using a time-dependent Neumann boundary condition

We generate a shear (S) wave using time-dependent Neumann boundary
condition on one of of a bar. We place a Dirichlet BC on the other end
so that the wave continues to propagate back and forth along the bar.

The simulation parameters are in the `pylithapp.cfg` and
`step02_swave_reflected.cfg` files.

To run the example:
```
pylith step02_swave_reflected.cfg
```


## Step03: Shear wave and absorbing boundary

We replace the Dirichlet BC in Step02 with an absorbing boundary so that the
shear wave does not reflect off of the boundary.

The simulation parameters are in the `pylithapp.cfg` and
`step03_swave_absorbed.cfg` files.

To run the example:
```
pylith step03_swave_absorbed.cfg
```

## Step04: Shear wave from prescribed slip



