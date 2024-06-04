# Examples: 2-D Reverse Fault and Gravitational Body Forces

This suite of examples demonstrates use of gravitational body forces, surface
loadging, slip on multiple, intersecting reverse faults, and elastic and viscoelastic
bulk rheologies in 2D.

## Step 1: Gravitational body forces and linear isotropic elasticity

We use zero-displacement Dirichlet boundary conditions on the +x, -x, and -y
boundaries, and gravitational body forces.

To run the example:

```bash
# Coarse mesh
pylith step01a_gravity.cfg

# Refined mesh
pylith step01b_gravity.cfg

# Higher order discretization
pylith step01c_gravity.cfg
```

## Step 2: Gravitational body forces and linear isotropic elasticity with a reference stress state

Same as Step 1 but with a reference stress state to balance the gravitational body forces.

To run the example:

```bash
pylith step02_gravity_refstate.cfg
```

## Step 3: Gravitational body forces and linear isotropic incompressible elasticity

Same as Step 1 but with incompressible elastasticity.

To run the example:

```bash
pylith step03_gravity_incompressible.cfg
```

## Step 4: Surface tractions and linear isotropic linear elasticity

This problem imposes normal tractions on the ground surface, which could
be a proxy for loading by water or ice.

To run the example:

```bash
# Coarse mesh
pylith step04a_surfload.cfg

# Refined mesh
pylith step04b_surfload.cfg

# Higher order discretization
pylith step04c_surfload.cfg
```

## Step 5: Earthquake rupture on one fault and linear isotropic linear elasticity

This problem has the same zero-displacement BC used in Steps 1-3, but
includes earthquake rupture of the main thrust fault.

To run the example:

```bash
# Coarse mesh
pylith step05a_onefault.cfg

# Refined mesh
pylith step05b_onefault.cfg

# Higher order discretization
pylith step05c_onefault.cfg
```

## Step 6: Earthquake rupture on two faults and linear isotropic linear elasticity

This problem is similar to Step 5, but includes earthquake rupture on both
the main thrust fault and a splay fault.

To run the example:

```bash
pylith step06_twofaults_elastic.cfg
```

## Step 7: Earthquake rupture on two faults and linear isotropic Maxwell viscoelastic rheology

We replace the linear, isotropic elastic bulk rheology for the slab in Step 6 with a linear, isotropic Maxwell viscoelastic bulk rheology. We also run the simulation for a duration of 100 years.

To run the example:

```bash
pylith step07_twofaults_maxwell.cfg
```

## Step 8: Earthquake rupture on two faults and linear isotropic powerlaw viscoelastic rheology

We replace the linear, isotropic Maxwell viscoelastic bulk rheology for the slab in Step 7 with an isotropic powerlaw viscoelastic bulk rheology.

To run the example:

```bash
pylith step08_twofaults_powerlaw.cfg
```
