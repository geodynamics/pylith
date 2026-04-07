# 2-D Linear Elasticity with Prescribed Slip

We evaluate solver performance for 2-D plane strain linear elasticity with prescribed slip using Step 2 (variable slip)  from `examples/crustal-strikeslip-2d`.
We create a suite of simulations in which we use uniform global refinement to reduce the discretization size by a factor of 2 with up to 5 levels of refinement.
We assess the performance based on the number of iterations in the linear solver, the time required to solve, and the total runtime.

:::{note}
We generate the meshes with a gradient in the mesh size in which the mesh size increases at a geometric rate with distance from the faults.
Uniform global refinement of the meshes does not preserve this geometric rate, because the refinement within a cell is uniform.
Generating new meshes for each desired discretization size on the fault allows preserving the geometric rate and will result in a smaller number of unknowns relative to the discretization size on the fault compared with uniform global refinement.
:::

## Running the Benchmark

```{code-block} console
# Generate the finite-element mesh using Gmsh
./generate_gmsh.py --write --filename=mesh_tri.msh

# Run the benchmark simulations
pylith step02_fieldsplit_selfp.cfg >& logs/step02a_fieldsplit_selfp.stdout
pylith step02_fieldsplit_selfp.cfg step02b_fieldsplit_selfp.cfg >& logs/step02b_fieldsplit_selfp.stdout

pylith step03_vpbjacobi.cfg >& step03a_vpbjacobi.stdout
pylith step03_vpbjacobi.cfg step03b_vpbjacobi.cfg >& step03b_vpbjacobi.stdout
pylith step03_vpbjacobi.cfg step03c_vpbjacobi.cfg >& step03c_vpbjacobi.stdout
pylith step03_vpbjacobi.cfg step03d_vpbjacobi.cfg >& step03d_vpbjacobi.stdout
pylith step03_vpbjacobi.cfg step03e_vpbjacobi.cfg >& step03e_vpbjacobi.stdout
pylith step03_vpbjacobi.cfg step03f_vpbjacobi.cfg >& step03f_vpbjacobi.stdout

pylith step04_vpbjacobi_tunefine.cfg >& step04a_vpbjacobi_tunefine.stdout
pylith step04_vpbjacobi_tunefine.cfg step04b_vpbjacobi_tunefine.cfg >& step04b_vpbjacobi_tunefine.stdout
pylith step04_vpbjacobi_tunefine.cfg step04c_vpbjacobi_tunefine.cfg >& step04c_vpbjacobi_tunefine.stdout
pylith step04_vpbjacobi_tunefine.cfg step04d_vpbjacobi_tunefine.cfg >& step04d_vpbjacobi_tunefine.stdout
pylith step04_vpbjacobi_tunefine.cfg step04e_vpbjacobi_tunefine.cfg >& step04e_vpbjacobi_tunefine.stdout
pylith step04_vpbjacobi_tunefine.cfg step04f_vpbjacobi_tunefine.cfg >& step04f_vpbjacobi_tunefine.stdout

pylith step06_vpbjacobi_dispscale.cfg >& step06a_vpbjacobi_dispscale.stdout
pylith step06_vpbjacobi_dispscale.cfg step06b_vpbjacobi_dispscale.cfg >& step06b_vpbjacobi_dispscale.stdout
pylith step06_vpbjacobi_dispscale.cfg step06c_vpbjacobi_dispscale.cfg >& step06c_vpbjacobi_dispscale.stdout
pylith step06_vpbjacobi_dispscale.cfg step06d_vpbjacobi_dispscale.cfg >& step06d_vpbjacobi_dispscale.stdout
pylith step06_vpbjacobi_dispscale.cfg step06e_vpbjacobi_dispscale.cfg >& step06e_vpbjacobi_dispscale.stdout
pylith step06_vpbjacobi_dispscale.cfg step06f_vpbjacobi_dispscale.cfg >& step06f_vpbjacobi_dispscale.stdout
```

## Results

For small problems (about 10,000) unknowns we find little difference in the runtimes across the different solver settings.
As the number of unknowns increases towards 100,000, the number of iterations required by the Schur complement field split preconditioner increases dramatically, resulting in dramatic increases in runtime.
The geometric algegraic multigrid preconditioner with variable point block Jacobi preconditioner at the fine scale with default parameters works well up to about 1 million unknowns.
Fine-tuning the fine-scale multigrid preconditioner parameters gives good performance across the entire range of problem sizes considered.

:::{note}
Using independent displacement amd length scales does not improve solver performance, but it allows use of solver tolerances that are relative to the displacement scale and independent of discretization size.
:::

:::{figure-md} fig:benchmarks:performance:linear:elasticity:solver:2d
<img src="figs/solver-performance-2d.*" alt="" width="80em"/>

Solver performance for 2-D plane strain elasticity with prescribed slip.
The left panel shows the number of iterations used by the linear solver as a function of the problem size (number of unknowns).
The right panels shows the total (squares) and solve (triangles) runtime as a function of hte problem size (number of unknowns).
The dashed lines indicate the solves required more than 10,000 iterations at the next level of refinement.
:::
