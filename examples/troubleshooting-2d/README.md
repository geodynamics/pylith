# Troubleshooting: 2D Reverse Fault

This suite of examples demonstrates how to use a variety of techniques
to troubleshoot errors in running PyLith simulations.  We introduce 
several common and a few uncommon mistakes into input files
from `exmaples/reverse-2d`.

## Step 1: Gravitational body forces and linear isotropic elasticity

We use zero-displacement Dirichlet boundary conditions on the +x, -x, and -y
boundaries, and gravitational body forces.

To run the example:
```
pylith step01a_gravity.cfg
```

## Step 6: Earthquake rupture on two faults and linear isotropic linear elasticity

This problem is similar to Step 5, but includes earthquake rupture on both
the main thrust fault and a splay fault.

To run the example:
```
pylith step06_twofaults_elastic.cfg
```

