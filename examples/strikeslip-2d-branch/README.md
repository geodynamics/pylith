## Made during the 2023 Pylith Hackathon ##
# Contributors 
#   Zechao Zhuo: Gmesh
#   Evan Marschall: Cubit

# Examples: 2D Strike-Slip Branch Fault

This suite of examples demonstrates some basic concepts of using
PyLith with a simple through going branch fault in 2D. The angle between the main and branch fault is held at 30 degrees, but this can be adjusted in the meshing files. 

## Meshing

We provide mesh files generated using Gmsh and Cubit.
We also include a Python script for generating the finite-element mesh with
triangular cells using Gmsh and Journal files for generating the
finite-element mesh with triangular cell using Cubit.

## Notes on Running

For these pyltih examples must be run with several flags following the intial inputs. This may be fixed in future releases

example of running : pylith file.cfg --petsc.snes_view --petsc.ksp_monitor_true_residual --petsc.pc_fieldsplit_schur_factorization_type=full --petsc.pc_fieldsplit_schur_precondition=full --petsc.fieldsplit_displacement_pc_type=lu --petsc.fieldsplit_lagrange_multiplier_fault_pc_type=lu

## Step 1 Gmsh: Static Uniform Coseismic Slip

This example involves a static simulation that solves for the deformation from prescribed coseismic slip on the branch fault system for the Gmsh mesh. We specify 4 meters of right-lateral slip on the main fault and 2 meters of right-lateral slip .

To run the example:
```
pylith step01_slip.cfg --petsc.snes_view --petsc.ksp_monitor_true_residual --petsc.pc_fieldsplit_schur_factorization_type=full --petsc.pc_fieldsplit_schur_precondition=full --petsc.fieldsplit_displacement_pc_type=lu --petsc.fieldsplit_lagrange_multiplier_fault_pc_type=lu
```
## Step 1 Cubit: Static Uniform Coseismic Slip

The same as step 1 Gmsh, but for the cubit mesh

This example involves a static simulation that solves for the deformation from prescribed coseismic slip on the branch fault system for the Gmsh mesh. We specify 4 meters of right-lateral slip on the main fault and 2 meters of right-lateral slip .

To run the example:
```
pylith step01_slip_cubit.cfg --petsc.snes_view --petsc.ksp_monitor_true_residual --petsc.pc_fieldsplit_schur_factorization_type=full --petsc.pc_fieldsplit_schur_precondition=full --petsc.fieldsplit_displacement_pc_type=lu --petsc.fieldsplit_lagrange_multiplier_fault_pc_type=lu
```

## Step 2 Gmsh: Single Earthquake Rupture and Velocity Boundary Conditions

This example involves a quasistatic simulation that solves for the deformation from velocity boundary conditions and prescribed coseismic slip on the branch fault system.

We let strain accumulate due to the motion of the boundaries (1 cm/year) and then release the strain by prescribing slip at t=100 years.

Slip values are 4 meters of right-lateral slip on the main fault and 2 meters of right-lateral slip on the branch segment.

To run the example:
```
pylith step02_slip_velbc.cfg --petsc.snes_view --petsc.ksp_monitor_true_residual --petsc.pc_fieldsplit_schur_factorization_type=full --petsc.pc_fieldsplit_schur_precondition=full --petsc.fieldsplit_displacement_pc_type=lu --petsc.fieldsplit_lagrange_multiplier_fault_pc_type=lu
```

## Step 3 Cubit: Single Earthquake with Variable slip 

This example involves a static simulation that solves for the deformation from prescribed coseismic slip on the branch fault system for the Gmsh mesh. We specify variable slip along both the main and branch faults which are called in from external .db files.

To run the example:
```
pylith step03_slip_cubit.cfg --petsc.snes_view --petsc.ksp_monitor_true_residual --petsc.pc_fieldsplit_schur_factorization_type=full --petsc.pc_fieldsplit_schur_precondition=full --petsc.fieldsplit_displacement_pc_type=lu --petsc.fieldsplit_lagrange_multiplier_fault_pc_type=lu
```
## Step 4 Cubit: Creeping Main Fault with Coseismic Slip on Branch 

This example has left lateral creep of 2 mm/yr on the main fault and a single rupture of 2 m on the branch fault after 50 years.

To run the example:
```
pylith step04_slip_cubit.cfg --petsc.snes_view --petsc.ksp_monitor_true_residual --petsc.pc_fieldsplit_schur_factorization_type=full --petsc.pc_fieldsplit_schur_precondition=full --petsc.fieldsplit_displacement_pc_type=lu --petsc.fieldsplit_lagrange_multiplier_fault_pc_type=lu
```
