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
  _jacobianLumped(0),
  _fields(0),
  _needVelocity(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::Formulation::~Formulation(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::Formulation::deallocate(void)
{ // deallocate
  _jacobian = 0; // :TODO: Use shared pointer.
  _jacobianLumped = 0; // :TODO: Use shared pointer.
  _fields = 0; // :TODO: Use shared pointer.
} // deallocate
  
// ----------------------------------------------------------------------
// Return the fields
const pylith::topology::SolutionFields&
pylith::problems::Formulation::fields(void) const
{ // fields
  return *this->_fields;
} // fields

// ----------------------------------------------------------------------
// Get flag indicating whether we need to compute velocity at time t.
bool
pylith::problems::Formulation::needVelocity(void) const
{ // needVelocity
  return _needVelocity;
} // needVelocity
  
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
// Initialize formulation.
void
pylith::problems::Formulation::initialize(void)
{ // initialize

  // Determine whether we need to compute velocity
  _needVelocity = false;
  int numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i) {
    assert(0 != _meshIntegrators[i]);
    if (_meshIntegrators[i]->needVelocity())
      _needVelocity = true;
  } // fpr
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i) {
    assert(0 != _submeshIntegrators[i]);
    if (_submeshIntegrators[i]->needVelocity())
      _needVelocity = true;
  } // for
} // initialize

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
// Update handles and parameters for reforming the Jacobian and
// residual.
void
pylith::problems::Formulation::updateSettings(topology::Field<topology::Mesh>* jacobian,
					      topology::SolutionFields* fields,
					      const double t,
					      const double dt)
{ // updateSettings
  assert(0 != jacobian);
  assert(0 != fields);
  assert(dt > 0.0);

  _jacobianLumped = jacobian;
  _fields = fields;
  _t = t;
  _dt = dt;
} // updateSettings

// ----------------------------------------------------------------------
// Reform system residual.
void
pylith::problems::Formulation::reformResidual(const PetscVec* tmpResidualVec,
					      const PetscVec* tmpSolutionVec)
{ // reformResidual
  assert(0 != _fields);

  // Update section view of field.
  if (0 != tmpSolutionVec) {
    topology::Field<topology::Mesh>& solution = _fields->solution();
    solution.scatterVectorToSection(*tmpSolutionVec);
  } // if

  // Compute velocity if necessary.
  if (_needVelocity)
    _calcVelocity();

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
    _submeshIntegrators[i]->integrateResidualAssembled(residual, _t, _fields);

  // Update PETSc view of residual
  if (0 != tmpResidualVec)
    residual.scatterSectionToVector(*tmpResidualVec);
  else
    residual.scatterSectionToVector();

  // TODO: Move this to SolverLinear 
  if (0 != tmpResidualVec)
    VecScale(*tmpResidualVec, -1.0);
} // reformResidual

// ----------------------------------------------------------------------
// Reform system residual.
void
pylith::problems::Formulation::reformResidualLumped(const PetscVec* tmpResidualVec,
                const PetscVec* tmpSolutionVec)
{ // reformResidualLumped
  assert(0 != _fields);

  // Update section view of field.
  if (0 != tmpSolutionVec) {
    topology::Field<topology::Mesh>& solution = _fields->solution();
    solution.scatterVectorToSection(*tmpSolutionVec);
  } // if

  // Compute velocity if necessary.
  if (_needVelocity)
    _calcVelocity();

  // Set residual to zero.
  topology::Field<topology::Mesh>& residual = _fields->get("residual");
  residual.zero();

  // Add in contributions that require assembly.
  int numIntegrators = _meshIntegrators.size();
  assert(numIntegrators > 0); // must have at least 1 bulk integrator
  for (int i=0; i < numIntegrators; ++i) {
    _meshIntegrators[i]->timeStep(_dt);
    _meshIntegrators[i]->integrateResidualLumped(residual, _t, _fields);
  } // for
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i) {
    _submeshIntegrators[i]->timeStep(_dt);
    _submeshIntegrators[i]->integrateResidualLumped(residual, _t, _fields);
  } // for

  // Assemble residual.
  residual.complete();

  // Add in contributions that do not require assembly.
  numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i]->integrateResidualLumpedAssembled(residual, _t, _fields);
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i]->integrateResidualLumpedAssembled(residual, _t, _fields);

  // Update PETSc view of residual
  if (0 != tmpResidualVec)
    residual.scatterSectionToVector(*tmpResidualVec);
  else
    residual.scatterSectionToVector();

  // TODO: Move this to SolverLinear
  if (0 != tmpResidualVec)
    VecScale(*tmpResidualVec, -1.0);
} // reformResidualLumped

// ----------------------------------------------------------------------
// Reform system Jacobian.
void
pylith::problems::Formulation::reformJacobian(const PetscVec* tmpSolutionVec)
{ // reformJacobian
  assert(0 != _jacobian);
  assert(0 != _fields);

  // Update section view of field.
  if (0 != tmpSolutionVec) {
    topology::Field<topology::Mesh>& solution = _fields->solution();
    solution.scatterVectorToSection(*tmpSolutionVec);
  } // if

  // Set jacobian to zero.
  _jacobian->zero();

  // Add in contributions that do not require assembly.
  int numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i]->integrateJacobianAssembled(_jacobian, _t, _fields);
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i]->integrateJacobianAssembled(_jacobian, _t, _fields);

  // Flush assembled portion.
  _jacobian->assemble("flush_assembly");

  // Add in contributions that require assembly.
  numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i]->integrateJacobian(_jacobian, _t, _fields);
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i]->integrateJacobian(_jacobian, _t, _fields);
  
  // Assemble jacobian.
  _jacobian->assemble("final_assembly");

} // reformJacobian

// ----------------------------------------------------------------------
// Reform system Jacobian.
void
pylith::problems::Formulation::reformJacobianLumped(void)
{ // reformJacobian
  assert(0 != _jacobianLumped);
  assert(0 != _fields);

  // Set jacobian to zero.
  _jacobianLumped->zero();

  // Add in contributions that require assembly.
  int numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i]->integrateJacobian(_jacobianLumped, _t, _fields);
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i]->integrateJacobian(_jacobianLumped, _t, _fields);
  
  // Assemble jacbian.
  _jacobianLumped->complete();

  // Add in contributions that do not require assembly.
  numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i]->integrateJacobianAssembled(_jacobianLumped, 
						    _t, _fields);
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i]->integrateJacobianAssembled(_jacobianLumped,
						       _t, _fields);

} // reformJacobian

// ----------------------------------------------------------------------
// Constrain solution space.
void
pylith::problems::Formulation::constrainSolnSpace(const PetscVec* tmpSolutionVec)
{ // constrainSolnSpace
  assert(0 != tmpSolutionVec);
  assert(0 != _fields);

  // Update section view of field.
  if (0 != tmpSolutionVec) {
    topology::Field<topology::Mesh>& solution = _fields->solution();
    solution.scatterVectorToSection(*tmpSolutionVec);
  } // if

  int numIntegrators = _meshIntegrators.size();
  assert(numIntegrators > 0); // must have at least 1 bulk integrator
  for (int i=0; i < numIntegrators; ++i) {
    _meshIntegrators[i]->timeStep(_dt);
    _meshIntegrators[i]->constrainSolnSpace(_fields, _t, *_jacobian);
  } // for
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i) {
    _submeshIntegrators[i]->timeStep(_dt);
    _submeshIntegrators[i]->constrainSolnSpace(_fields, _t, *_jacobian);
  } // for

  // Update PETScVec of solution for changes to Lagrange multipliers.
  if (0 != tmpSolutionVec) {
    topology::Field<topology::Mesh>& solution = _fields->solution();
    solution.scatterSectionToVector(*tmpSolutionVec);
  } // if
} // constrainSolnSpace

// ----------------------------------------------------------------------
// Adjust solution from solver with lumped Jacobian to match Lagrange
//  multiplier constraints.
void
pylith::problems::Formulation::adjustSolnLumped(void)
{ // adjustSolnLumped
  int numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i]->adjustSolnLumped(_fields, *_jacobianLumped);

  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i]->adjustSolnLumped(_fields, *_jacobianLumped);
} // adjustSolnLumped


// End of file 
