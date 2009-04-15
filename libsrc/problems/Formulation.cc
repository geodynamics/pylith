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
pylith::problems::Formulation::Formulation(void) :
  _t(0.0),
  _dt(0.0),
  _jacobian(0),
  _fields(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::Formulation::~Formulation(void)
{ // destructor
  _jacobian = 0; // Handle only, we do not manage the memory.
  _fields = 0; // Handle only, we do not manage the memory.
} // destructor

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
// Update handles and parameters for reforming the Jacobian and
// residual.
void
pylith::problems::Formulation::updateSettings(topology::Jacobian* jacobian,
					      topology::SolutionFields* fields,
					      const double t,
					      const double dt)
{ // updateSettings
  assert(0 != jacobian);
  assert(0 != fields);
  assert(dt > 0.0);

  _jacobian = jacobian;
  _fields = fields;
  _t = t;
  _dt = dt;
} // updateSettings


// ----------------------------------------------------------------------
// Reform system residual.
void
pylith::problems::Formulation::reformResidual(Vec solutionVec, Vec residualVec)
{ // reformResidual
  assert(0 != _fields);

  // Need to pass these Vecs for updating

  // Update section view of field.
  topology::Field<topology::Mesh>& solution = _fields->get("solution");
  solution.scatterVectorToSection(solutionVec);

  // Set residual to zero.
  topology::Field<topology::Mesh>& residual = _fields->get("residual");
  residual.zero();

  // Add in contributions that require assembly.
  int numIntegrators = _meshIntegrators.size();
  assert(numIntegrators > 0); // must have at least 1 bulk integrator
  for (int i=0; i < numIntegrators; ++i) {
    _meshIntegrators[i]->timeStep(_dt);
    _meshIntegrators[i]->integrateResidual(residual, _t, _fields);
  } // for
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i) {
    _submeshIntegrators[i]->timeStep(_dt);
    _submeshIntegrators[i]->integrateResidual(residual, _t, _fields);
  } // for

  // Assemble residual.
  residual.complete();

  // Add in contributions that do not require assembly.
  numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i]->integrateResidualAssembled(residual, _t, _fields);
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i]->integrateResidual(residual, _t, _fields);

  // Update PETSc view of residual
  residual.scatterSectionToVector(residualVec);

  // TODO: Move this to SolverLinear
  VecScale(residualVec, -1.0);
} // reformResidual

// ----------------------------------------------------------------------
// Reform system Jacobian.
void
pylith::problems::Formulation::reformJacobian(void)
{ // reformJacobian
  reformJacobian(NULL);
} // reformJacobian
void
pylith::problems::Formulation::reformJacobian(Vec solutionVec)
{ // reformJacobian
  assert(0 != _jacobian);
  assert(0 != _fields);

  // Need to pass these Vecs for updating

  // Update section view of field.
  if (solutionVec != NULL) {
    topology::Field<topology::Mesh>& solution = _fields->get("solution");
    solution.scatterVectorToSection(solutionVec);
  }

  // Set residual to zero.
  _jacobian->zero();

  // Add in contributions that do not require assembly.
  int numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i]->integrateJacobianAssembled(_jacobian, _t, _fields);
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i]->integrateJacobianAssembled(_jacobian, _t, _fields);

  // Assemble residual.
  _jacobian->assemble("flush_assembly");

  // Add in contributions that require assembly.
  numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i]->integrateJacobian(_jacobian, _t, _fields);
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i]->integrateJacobian(_jacobian, _t, _fields);
  
  // Assemble residual.
  _jacobian->assemble("final_assembly");
} // reformJacobian


// End of file 
