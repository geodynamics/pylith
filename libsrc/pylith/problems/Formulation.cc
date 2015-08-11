// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "Formulation.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/feassemble/Integrator.hh" // USES Integrator
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Formulation::Formulation(void) :
  _t(0.0),
  _dt(0.0),
  _jacobian(0),
  _customConstraintPCMat(0),
  _jacobianLumped(0),
  _fields(0),
  _isJacobianSymmetric(false),
  _splitFields(false)
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
  PYLITH_METHOD_BEGIN;

  _jacobian = 0; // :TODO: Use shared pointer.
  _jacobianLumped = 0; // :TODO: Use shared pointer.
  _fields = 0; // :TODO: Use shared pointer.

#if 0   // :KLUDGE: Assume Solver deallocates matrix.
  PetscErrorCode err = 0;
  if (_customConstraintPCMat) {
    err = PetscObjectDereference((PetscObject) _customConstraintPCMat);PYLITH_CHECK_ERROR(err);
    _customConstraintPCMat = 0;
  } // if
#else
  _customConstraintPCMat = 0;
#endif

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Set flag for splitting fields.
void
pylith::problems::Formulation::splitFields(const bool flag)
{ // splitFields
  _splitFields = flag;
} // splitFields

// ----------------------------------------------------------------------
// Get flag for splitting fields.
bool
pylith::problems::Formulation::splitFields(void) const
{ // splitFields
  return _splitFields;
} // splitFields

// ----------------------------------------------------------------------
// Set flag for using custom preconditioner for Lagrange constraints.
void
pylith::problems::Formulation::useCustomConstraintPC(const bool flag)
{ // useCustomConstraintPC
  _useCustomConstraintPC = flag;
} // useCustomConstraintPC

// ----------------------------------------------------------------------
// Get flag indicating use of custom conditioner for Lagrange constraints.
bool
pylith::problems::Formulation::useCustomConstraintPC(void) const
{ // useCustomConstraintPC
  return _useCustomConstraintPC;
} // useCustomConstraintPC

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
pylith::problems::Formulation::isJacobianSymmetric(void) const
{ // isJacobianSymmetric
  return _isJacobianSymmetric;
} // isJacobianSymmetric
  
// ----------------------------------------------------------------------
// Set integrators over the mesh.
void
pylith::problems::Formulation::integrators(feassemble::Integrator* integratorArray[],
					   const int numIntegrators)
{ // integrators
  assert( (!integratorArray && 0 == numIntegrators) ||
	  (integratorArray && 0 < numIntegrators) );
  _integrators.resize(numIntegrators);
  for (int i=0; i < numIntegrators; ++i)
    _integrators[i] = integratorArray[i];
} // integrators
  
// ----------------------------------------------------------------------
// Set handle to preconditioner.
void
pylith::problems::Formulation::customPCMatrix(PetscMat& mat)
{ // preconditioner
  _customConstraintPCMat = mat;

#if 0 // :KLUDGE: Assume solver deallocates matrix
  PetscErrorCode err = 0;
  err = PetscObjectReference((PetscObject) mat); PYLITH_CHECK_ERROR(err);
#endif
} // preconditioner

// ----------------------------------------------------------------------
// Update handles and parameters for reforming the Jacobian and
// residual.
void
pylith::problems::Formulation::updateSettings(topology::Jacobian* jacobian,
					      topology::SolutionFields* fields,
					      const PylithScalar t,
					      const PylithScalar dt)
{ // updateSettings
  assert(jacobian);
  assert(fields);
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
pylith::problems::Formulation::updateSettings(topology::Field* jacobian,
					      topology::SolutionFields* fields,
					      const PylithScalar t,
					      const PylithScalar dt)
{ // updateSettings
  assert(jacobian);
  assert(fields);
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
  PYLITH_METHOD_BEGIN;

  assert(_fields);

  // Update section view of field.
  if (tmpSolutionVec) {
    topology::Field& solution = _fields->solution();
    solution.scatterGlobalToLocal(*tmpSolutionVec);
  } // if

  // Update rate fields (must be consistent with current solution).
  calcRateFields();  

  // Set residual to zero.
  topology::Field& residual = _fields->get("residual");
  residual.zeroAll();

  // Add in contributions that require assembly.
  const int numIntegrators = _integrators.size();
  assert(numIntegrators > 0); // must have at least 1 integrator
  for (int i=0; i < numIntegrators; ++i) {
    _integrators[i]->timeStep(_dt);
    _integrators[i]->integrateResidual(residual, _t, _fields);
  } // for

  // Assemble residual.
  residual.complete();

  // Update PETSc view of residual
  if (tmpResidualVec)
    residual.scatterLocalToGlobal(*tmpResidualVec);
  else
    residual.scatterLocalToGlobal();

  // TODO: Move this to SolverLinear 
  if (tmpResidualVec)
    VecScale(*tmpResidualVec, -1.0);

  PYLITH_METHOD_END;
} // reformResidual

// ----------------------------------------------------------------------
// Reform system Jacobian.
void
pylith::problems::Formulation::reformJacobian(const PetscVec* tmpSolutionVec)
{ // reformJacobian
  PYLITH_METHOD_BEGIN;

  assert(_jacobian);
  assert(_fields);

  // Update section view of field.
  if (tmpSolutionVec) {
    topology::Field& solution = _fields->solution();
    solution.scatterGlobalToLocal(*tmpSolutionVec);
  } // if

  // Set jacobian to zero.
  _jacobian->zero();

  // Add in contributions that require assembly.
  const int numIntegrators = _integrators.size();
  for (int i=0; i < numIntegrators; ++i) {
    _integrators[i]->integrateJacobian(_jacobian, _t, _fields);
  } // for
  
  // Assemble jacobian.
  _jacobian->assemble("final_assembly");

  if (_customConstraintPCMat) {
    // Recalculate preconditioner.
    for (int i=0; i < numIntegrators; ++i) {
      _integrators[i]->calcPreconditioner(&_customConstraintPCMat, _jacobian, _fields);
    } // for

    MatAssemblyBegin(_customConstraintPCMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_customConstraintPCMat, MAT_FINAL_ASSEMBLY);

#if 0 // debugging
    std::cout << "Preconditioner Matrix" << std::endl;
    MatView(_customConstraintPCMat, PETSC_VIEWER_STDOUT_WORLD);
#endif
  } // if

  PYLITH_METHOD_END;
} // reformJacobian

// ----------------------------------------------------------------------
// Reform system Jacobian.
void
pylith::problems::Formulation::reformJacobianLumped(void)
{ // reformJacobianLumped
  PYLITH_METHOD_BEGIN;

  assert(_jacobianLumped);
  assert(_fields);

  // Set jacobian to zero.
  _jacobianLumped->zeroAll();

  // Add in contributions that require assembly.
  const int numIntegrators = _integrators.size();
  for (int i=0; i < numIntegrators; ++i) {
    _integrators[i]->integrateJacobian(_jacobianLumped, _t, _fields);
  } // for
  
  // Assemble jacbian.
  _jacobianLumped->complete();

  PYLITH_METHOD_END;
} // reformJacobianLumped

// ----------------------------------------------------------------------
// Constrain solution space.
void
pylith::problems::Formulation::constrainSolnSpace(const PetscVec* tmpSolutionVec)
{ // constrainSolnSpace
  PYLITH_METHOD_BEGIN;

  assert(tmpSolutionVec);
  assert(_fields);

  topology::Field& solution = _fields->solution();

  if (!_fields->hasField("dispIncr adjust")) {
    _fields->add("dispIncr adjust", "dispIncr_adjust");
    topology::Field& adjust = _fields->get("dispIncr adjust");
    adjust.cloneSection(solution);
  } // for
  topology::Field& adjust = _fields->get("dispIncr adjust");
  adjust.zeroAll();

  // Update section view of field.
  if (tmpSolutionVec) {
    solution.scatterGlobalToLocal(*tmpSolutionVec);
  } // if

  const int numIntegrators = _integrators.size();
  assert(numIntegrators > 0); // must have at least 1 bulk integrator
  for (int i=0; i < numIntegrators; ++i) {
    _integrators[i]->timeStep(_dt);
    _integrators[i]->constrainSolnSpace(_fields, _t, *_jacobian);
  } // for

  adjust.complete();
  solution += adjust;  

  // Update PETScVec of solution for changes to Lagrange multipliers.
  if (tmpSolutionVec) {
    solution.scatterLocalToGlobal(*tmpSolutionVec);
  } // if

  PYLITH_METHOD_END;
} // constrainSolnSpace

// ----------------------------------------------------------------------
// Adjust solution from solver with lumped Jacobian to match Lagrange
//  multiplier constraints.
void
pylith::problems::Formulation::adjustSolnLumped(void)
{ // adjustSolnLumped
  PYLITH_METHOD_BEGIN;

  topology::Field& solution = _fields->solution();

  if (!_fields->hasField("dispIncr adjust")) {
    _fields->add("dispIncr adjust", "dispIncr_adjust");
    topology::Field& adjust = _fields->get("dispIncr adjust");
    adjust.cloneSection(solution);
  } // for
  topology::Field& adjust = _fields->get("dispIncr adjust");
  adjust.zeroAll();

  const int numIntegrators = _integrators.size();
  for (int i=0; i < numIntegrators; ++i) {
    _integrators[i]->adjustSolnLumped(_fields, _t, *_jacobianLumped);
  } // for

  adjust.complete();
  solution += adjust;

  PYLITH_METHOD_END;
} // adjustSolnLumped

#include "pylith/meshio/DataWriterHDF5.hh"
// ----------------------------------------------------------------------
void
pylith::problems::Formulation::printState(PetscVec* solutionVec,
					  PetscVec* residualVec,
					  PetscVec* solution0Vec,
					  PetscVec* searchDirVec)
{ // printState
  assert(solutionVec);
  assert(residualVec);
  assert(searchDirVec);

  meshio::DataWriterHDF5 writer;

  const topology::Mesh& mesh = _fields->mesh();

  writer.filename("state.h5");
  const int numTimeSteps = 1;
  writer.open(mesh, numTimeSteps);
   
  topology::Field& solution = _fields->solution();
  solution.scatterGlobalToLocal(*solutionVec);
  writer.writeVertexField(0.0, solution, mesh);
  solution.view("DIVERGED_SOLUTION");
  const char* label = solution.label();

  solution.label("solution_0");
  solution.scatterGlobalToLocal(*solution0Vec);
  writer.writeVertexField(0.0, solution, mesh);
  solution.view("DIVERGED_SOLUTION0");
  solution.label(label);

  topology::Field& residual = _fields->get("residual");
  residual.scatterGlobalToLocal(*residualVec);
  writer.writeVertexField(0.0, residual, mesh);
  residual.view("DIVERGED_RESIDUAL");

  residual.label("search_dir");
  residual.scatterGlobalToLocal(*searchDirVec);
  writer.writeVertexField(0.0, residual, mesh);
  residual.view("DIVERGED_SEARCHDIR");

  writer.close();
} // printState



// End of file 
