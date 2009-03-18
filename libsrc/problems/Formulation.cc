// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "Formulation.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Quadrature<Mesh>
#include "pylith/topology/SubMesh.hh" // USES Quadrature<SubMesh>
#include "pylith/feassemble/Quadrature.hh" // USES Integrator<Quadrature>
#include "pylith/feassemble/Integrator.hh" // USES Integrator
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Formulation::Formulation(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::Formulation::~Formulation(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Generic C interface for reformResidual for integration with
// PETSc SNES solvers.
void
pylith::problems::Formulation::reformResidual(PetscSNES snes,
					      PetscVec solutionVec,
					      PetscVec residualVec,
					      void* context)
{ // reformResidual
  assert(0 != context);
  ArgsResidual* args = (ArgsResidual*) context;
  assert(0 != args);
  assert(0 != args->object);
  assert(0 != args->fields);
  assert(0 != args->residual);

  // Copy solution information from PETSc vector into field
  args->fields->solution().scatterVectorToSection();

  // Reform residual
  args->object->reformResidual(args->residual, args->fields, 
			       args->t, args->dt);  

  // Copy residual information from field into PETSc vector
  args->residual->scatterSectionToVector();
} // reformResidual

// ----------------------------------------------------------------------
// Generic C interface for reformJacobian for integration with
// PETSc SNES solvers.
void
pylith::problems::Formulation::reformJacobian(PetscSNES snes,
					      PetscVec solutionVec,
					      PetscMat jacobianMat,
					      PetscMat preconditionerMat,
					      int* preconditionerLayout,
					      void* context)
{ // reformJacobian
  assert(0 != context);
  ArgsJacobian* args = (ArgsJacobian*) context;
  assert(0 != args);
  assert(0 != args->object);
  args->object->reformJacobian(args->jacobian, args->fields, args->t, args->dt);
} // reformJacobian

// ----------------------------------------------------------------------
// Set integrators over the mesh.
void
pylith::problems::Formulation::meshIntegrators(IntegratorMesh** integrators,
					       const int numIntegrators)
{ // meshIntegrators
  assert( (0 == integrators && 0 == numIntegrators) ||
	  (0 != integrators && 0 < numIntegrators) );
  _meshIntegrators.resize(numIntegrators);
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i] = integrators[i];
} // meshIntegrators
  
// ----------------------------------------------------------------------
// Set integrators over lower-dimension meshes.
void
pylith::problems::Formulation::submeshIntegrators(IntegratorSubMesh** integrators,
						  const int numIntegrators)
{ // submeshIntegrators
  assert( (0 == integrators && 0 == numIntegrators) ||
	  (0 != integrators && 0 < numIntegrators) );
  _submeshIntegrators.resize(numIntegrators);
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i] = integrators[i];
} // submeshIntegrators

// ----------------------------------------------------------------------
// Reform system residual.
void
pylith::problems::Formulation::reformResidual(
			     topology::Field<topology::Mesh>* const residual,
			     topology::SolutionFields* const fields,
			     const double t,
			     const double dt)
{ // reformResidual
  assert(0 != residual);
  assert(0 != fields);

  // Set residual to zero.
  residual->zero();

  // Add in contributions that require assembly.
  int numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i) {
    _meshIntegrators[i]->timeStep(dt);
    _meshIntegrators[i]->integrateResidual(*residual, t, fields);
  } // for
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i) {
    _submeshIntegrators[i]->timeStep(dt);
    _submeshIntegrators[i]->integrateResidual(*residual, t, fields);
  } // for

  // Assemble residual.
  residual->complete();

  // Add in contributions that do not require assembly.
  numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i]->integrateResidualAssembled(*residual, t, fields);
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i]->integrateResidual(*residual, t, fields);
} // reformResidual

// ----------------------------------------------------------------------
// Reform system Jacobian.
void
pylith::problems::Formulation::reformJacobian(
				     topology::Jacobian* jacobian,
				     topology::SolutionFields* const fields,
				     const double t,
				     const double dt)
{ // reformJacobian
  assert(0 != jacobian);
  assert(0 != fields);

  // Set residual to zero.
  jacobian->zero();

  // Add in contributions that do not require assembly.
  int numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i]->integrateJacobianAssembled(jacobian, t, fields);
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i]->integrateJacobianAssembled(jacobian, t, fields);

  // Assemble residual.
  jacobian->assemble("flush_assembly");

  // Add in contributions that require assembly.
  numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i]->integrateJacobian(jacobian, t, fields);
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i]->integrateJacobian(jacobian, t, fields);
  
  // Assemble residual.
  jacobian->assemble("final_assembly");
} // reformJacobian


// End of file 
