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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
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
  _customConstraintPCMat(0),
  _isJacobianSymmetric(false),
  _splitFields(false),
  _splitFieldComponents(false)
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

#if 0   // :KLUDGE: Assume Solver deallocates matrix.
  PetscErrorCode err = 0;
  if (_customConstraintPCMat) {
    err = PetscObjectDereference((PetscObject) _customConstraintPCMat);CHECK_PETSC_ERROR(err);
    _customConstraintPCMat = 0;
  } // if
#else
  _customConstraintPCMat = 0;
#endif
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
// Set flag for splitting field components.
void
pylith::problems::Formulation::splitFieldComponents(const bool flag)
{ // splitFieldComponents
  _splitFieldComponents = flag;
} // splitFieldComponents

// ----------------------------------------------------------------------
// Get flag for splitting field components.
bool
pylith::problems::Formulation::splitFieldComponents(void) const
{ // splitFieldComponents
  return _splitFieldComponents;
} // splitFieldComponents

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
pylith::problems::Formulation::meshIntegrators(IntegratorMesh** integrators,
					       const int numIntegrators)
{ // meshIntegrators
  assert( (!integrators && 0 == numIntegrators) ||
	  (integrators && 0 < numIntegrators) );
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
  assert( (!integrators && 0 == numIntegrators) ||
	  (integrators && 0 < numIntegrators) );
  _submeshIntegrators.resize(numIntegrators);
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i] = integrators[i];
} // submeshIntegrators

// ----------------------------------------------------------------------
// Set handle to preconditioner.
void
pylith::problems::Formulation::customPCMatrix(PetscMat& mat)
{ // preconditioner
  _customConstraintPCMat = mat;

#if 0 // :KLUDGE: Assume solver deallocates matrix
  PetscErrorCode err = 0;
  err = PetscObjectReference((PetscObject) mat); CHECK_PETSC_ERROR(err);
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
pylith::problems::Formulation::updateSettings(topology::Field<topology::Mesh>* jacobian,
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
  assert(_fields);

  // Update section view of field.
  if (tmpSolutionVec) {
    topology::Field<topology::Mesh>& solution = _fields->solution();
    solution.scatterVectorToSection(*tmpSolutionVec);
  } // if

  // Update rate fields (must be consistent with current solution).
  calcRateFields();  

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

  // Update PETSc view of residual
  if (tmpResidualVec)
    residual.scatterSectionToVector(*tmpResidualVec);
  else
    residual.scatterSectionToVector();

  // TODO: Move this to SolverLinear 
  if (tmpResidualVec)
    VecScale(*tmpResidualVec, -1.0);
} // reformResidual

// ----------------------------------------------------------------------
// Reform system Jacobian.
void
pylith::problems::Formulation::reformJacobian(const PetscVec* tmpSolutionVec)
{ // reformJacobian
  assert(0 != _jacobian);
  assert(0 != _fields);

  // Update section view of field.
  if (tmpSolutionVec) {
    topology::Field<topology::Mesh>& solution = _fields->solution();
    solution.scatterVectorToSection(*tmpSolutionVec);
  } // if

  // Set jacobian to zero.
  _jacobian->zero();

  // Add in contributions that require assembly.
  int numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i]->integrateJacobian(_jacobian, _t, _fields);
  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i]->integrateJacobian(_jacobian, _t, _fields);
  
  // Assemble jacobian.
  _jacobian->assemble("final_assembly");

  if (_customConstraintPCMat) {
    // Recalculate preconditioner.
    numIntegrators = _meshIntegrators.size();
    for (int i=0; i < numIntegrators; ++i)
      _meshIntegrators[i]->calcPreconditioner(&_customConstraintPCMat,
					      _jacobian, _fields);
    numIntegrators = _submeshIntegrators.size();
    for (int i=0; i < numIntegrators; ++i)
      _submeshIntegrators[i]->calcPreconditioner(&_customConstraintPCMat,
						 _jacobian, _fields);

    MatAssemblyBegin(_customConstraintPCMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_customConstraintPCMat, MAT_FINAL_ASSEMBLY);

#if 0 // debugging
    std::cout << "Preconditioner Matrix" << std::endl;
    MatView(_customConstraintPCMat, PETSC_VIEWER_STDOUT_WORLD);
#endif
  } // if
} // reformJacobian

// ----------------------------------------------------------------------
// Reform system Jacobian.
void
pylith::problems::Formulation::reformJacobianLumped(void)
{ // reformJacobianLumped
  assert(_jacobianLumped);
  assert(_fields);

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

} // reformJacobianLumped

// ----------------------------------------------------------------------
// Constrain solution space.
void
pylith::problems::Formulation::constrainSolnSpace(const PetscVec* tmpSolutionVec)
{ // constrainSolnSpace
  assert(tmpSolutionVec);
  assert(_fields);

  topology::Field<topology::Mesh>& solution = _fields->solution();

  if (!_fields->hasField("dispIncr adjust")) {
    _fields->add("dispIncr adjust", "dispIncr_adjust");
    topology::Field<topology::Mesh>& adjust = _fields->get("dispIncr adjust");
    adjust.cloneSection(solution);
  } // for
  topology::Field<topology::Mesh>& adjust = _fields->get("dispIncr adjust");
  adjust.zero();

  // Update section view of field.
  if (tmpSolutionVec) {
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

  adjust.complete();
  solution += adjust;  

  // Update PETScVec of solution for changes to Lagrange multipliers.
  if (tmpSolutionVec) {
    solution.scatterSectionToVector(*tmpSolutionVec);
  } // if
} // constrainSolnSpace

// ----------------------------------------------------------------------
// Adjust solution from solver with lumped Jacobian to match Lagrange
//  multiplier constraints.
void
pylith::problems::Formulation::adjustSolnLumped(void)
{ // adjustSolnLumped
  topology::Field<topology::Mesh>& solution = _fields->solution();

  if (!_fields->hasField("dispIncr adjust")) {
    _fields->add("dispIncr adjust", "dispIncr_adjust");
    topology::Field<topology::Mesh>& adjust = _fields->get("dispIncr adjust");
    adjust.cloneSection(solution);
  } // for
  topology::Field<topology::Mesh>& adjust = _fields->get("dispIncr adjust");
  adjust.zero();

  int numIntegrators = _meshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _meshIntegrators[i]->adjustSolnLumped(_fields, _t, *_jacobianLumped);

  numIntegrators = _submeshIntegrators.size();
  for (int i=0; i < numIntegrators; ++i)
    _submeshIntegrators[i]->adjustSolnLumped(_fields, _t, *_jacobianLumped);

  adjust.complete();
  solution += adjust;
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

  meshio::DataWriterHDF5<topology::Mesh,topology::Field<topology::Mesh> > writer;

  const topology::Mesh& mesh = _fields->mesh();

  writer.filename("state.h5");
  const int numTimeSteps = 1;
  writer.open(mesh, numTimeSteps);
   
  topology::Field<topology::Mesh>& solution = _fields->solution();
  solution.scatterVectorToSection(*solutionVec);
  writer.writeVertexField(0.0, solution, mesh);
  solution.view("DIVERGED_SOLUTION");
  const char* label = solution.label();

  solution.label("solution_0");
  solution.scatterVectorToSection(*solution0Vec);
  writer.writeVertexField(0.0, solution, mesh);
  solution.view("DIVERGED_SOLUTION0");
  solution.label(label);

  topology::Field<topology::Mesh>& residual = _fields->get("residual");
  residual.scatterVectorToSection(*residualVec);
  writer.writeVertexField(0.0, residual, mesh);
  residual.view("DIVERGED_RESIDUAL");

  residual.label("search_dir");
  residual.scatterVectorToSection(*searchDirVec);
  writer.writeVertexField(0.0, residual, mesh);
  residual.view("DIVERGED_SEARCHDIR");

  writer.close();
} // printState



// End of file 
