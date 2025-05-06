# Component Implementations

## Problems

[`TimeDependent`](problems/TimeDependent.md)
: Time-dependent boundary value problems, including problems with a single time step (static).

[`GreensFns`](problems/GreensFns.md)
: Computation of Green's functions for static fault slip.

## Materials

[`Elasticity`](materials/Elasticity.md)
: Elasticity with and without faults.

[`IncompressibleElasticity`](materials/IncompressibleElasticity.md)
: Elasticity formulation for nearly incompressible materials.

[`Poroelasticity`](materials/Poroelasticity.md)
: Poroelasticity

## Bulk rheologies

[`IsotropicLinearElasticity`](materials/IsotropicLinearElasticity.md)
: Isotropic, linear elastic bulk rheology for `Elasticity`.

[`IsotropicLinearMaxwell`](materials/IsotropicLinearMaxwell.md)
: Isotropic, linear Maxwell visoelastic bulk rheology for `Elasticity`.

[`IsotropicLinearGenMaxwell`](materials/IsotropicLinearGenMaxwell.md)
: Isotropic, linear generalized Maxwell visoelastic bulk rheology for `Elasticity`.

[`IsotropicPowerLaw`](materials/IsotropicPowerLaw.md)
: Isotropic power-law visoelastic bulk rheology for `Elasticity`.

[`IsotropicLinearIncompElasticity`](materials/IsotropicLinearIncompElasticity.md)
: Isotropic, linear elastic bulk rheology for `IncompressibleElasticity`.

[`IsotropicLinearPoroelasticity`](materials/IsotropicLinearPoroelasticity.md)
: Isotropic, linear elastic bulk rheology for `Poroelasticity`.

## Boundary conditions

[`DirichletTimeDependent`](bc/DirichletTimeDependent.md)
: Time-dependent Dirichlet boundary condition.

[`NeumannTimeDependent`](bc/NeumannTimeDependent.md)
: Time-dependent Neumann boundary condition.

[`AbsorbingDampers`](bc/AbsorbingDampers.md)
: Absorbing boundary condition for dynamic simulations.

## Interface conditions

[`FaultCohesiveKin`](faults/FaultCohesiveKin.md)
: Earthquake ruptures using prescribed slip (kinematic rupture model) in `TimeDependent` problems.

[`FaultCohesiveImpulses`](faults/FaultCohesiveImpulses.md)
: Slip impulses for computing Green's functions in `GreensFns` problems.

## Kinematic earthquake ruptures

[`KinSrcStep`](faults/KinSrcStep.md)
: Slip time function using a step function for use in `FaultCohesiveKin`.

[`KinSrcRamp`](faults/KinSrcRamp.md)
: Slip time function using a ramp function for use in `FaultCohesiveKin`.

[`KinSrcConstRate`](faults/KinSrcConstRate.md)
: Slip time function using a constant rate function for use in `FaultCohesiveKin`.

[`KinSrcBrune`](faults/KinSrcBrune.md)
: Slip time function using the Brune far-field time function for use in `FaultCohesiveKin`.

[`KinSrcLiuCos`](faults/KinSrcLiuCos.md)
: Slip time function using three cosine functions for use in `FaultCohesiveKin`.

[`KinSrcTimeHistory`](faults/KinSrcTimeHistory.md)
: User-specified slip time function for use in `FaultCohesiveKin`.

## Initial conditions

[`InitialConditionDomain`](problems/InitialConditionDomain.md)
: Initial conditions (solution values) over the domain for use in `TimeDependent`.

[`InitialConditionPatch`](problems/InitialConditionPatch.md)
: Initial conditions (solution values) over a patch of the domain for use in `TimeDependent`.

## Solution fields

[`SolnDisp`](problems/SolnDisp.md)
: Solution field array with `displacement` subfield.

[`SolnDispLagrange`](problems/SolnDispLagrange.md)
: Solution field with `displacement` and `lagrange_multiplier_fault` subfields.

[`SolnDispPres`](problems/SolnDispPres.md)
: Solution field with `displacement` and `pressure` subfields.

[`SolnDispPresLagrange`](problems/SolnDispPresLagrange.md)
: Solution field with `displacement`, `pressure`, and `lagrange_multiplier_fault` subfields.

[`SolnDispPresTracStrain`](problems/SolnDispPresTracStrain.md)
: Solution field with `displacement`, `pressure`, and `trace_strain` subfields.

[`SolnDispPresTracStrainVelPdotTdot`](problems/SolnDispPresTracStrainVelPdotTdot.md)
: Solution field with `displacement`, `pressure`, `trace_strain`, `velocity`, `pressure_dot`, `trace_strain_dot` subfields.

[`SolnDispPresVel`](problems/SolnDispPresVel.md)
: Solution field with `displacement`, `pressure`, and `velocity` subfields.

[`SolnDispVel`](problems/SolnDispVel.md)
: Solution field with `displacement` and `velocity` subfields.

[`SolnDispVelLagrange`](problems/SolnDispVelLagrange.md)
: Solution field with `displacement`, `velocity`, `lagrange_multiplier_fault` subfields.

## Mesh readers

[`MeshIOAscii`](meshio/MeshIOAscii.md)
: Reader for PyLith ASCII finite-element mesh files.

[`MeshIOCubit`](meshio/MeshIOCubit.md)
: Reader for Cubit Exodus II finite-element mesh files.

[`MeshIOPetsc`](meshio/MeshIOPetsc.md)
: Reader for Gmsh finite-element mesh files and meshes generated with PETSc command line arguments.

## Mesh distributors

[`Distributor`](topology/Distributor.md)
: Distributor for distributing mesh information among processes.

## Mesh refiners

[`RefineUniform`](topology/RefineUniform.md)
: Uniform refinement for decreasing discretization size by a power of 2.

## Solution output

[`OutputSolnDomain`](meshio/OutputSolnDomain.md)
: Output of solution subfields over the entire domain.

[`OutputSolnBoundary`](meshio/OutputSolnBoundary.md)
: Output of solution subfields over an external boundary of the domain.

[`OutputSolnPoints`](meshio/OutputSolnPoints.md)
: Output of solution subfields at specific points in the domain.

## Output writers

[`DataWriterHDF5`](meshio/DataWriterHDF5.md)
: Write output to HDF5 files.

[`DataWriterHDF5Ext`](meshio/DataWriterHDF5Ext.md)
: Write output to HDF5 files with datasets in external binary files written using MPI-IO.

[`DataWriterVTK`](meshio/DataWriterVTK.md)
: Write output to VTK files.

## Output triggers

[`OutputTriggerStep`](meshio/OutputTriggerStep.md)
: Trigger output by skipping set number of time steps.

[`OutputTriggerTime`](meshio/OutputTriggerTime.md)
: Trigger output by skipping a set time interval.

## Progress monitors

[`ProgressMonitorStep`](problems/ProgressMonitorStep.md)
: Trigger output of simulation progress based on number of time steps.

[`ProgressMonitorTime`](problems/ProgressMonitorTime.md)
: Trigger output of simulation progress based on a time interval.

## SpatialData implementations

<https://spatialdata.readthedocs.io/en/latest/user/components/index.html>
