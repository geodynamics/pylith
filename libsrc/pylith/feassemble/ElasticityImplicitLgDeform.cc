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

#include "ElasticityImplicitLgDeform.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature
#include "CellGeometry.hh" // USES CellGeometry

#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN
#include "pylith/utils/lapack.h" // USES LAPACKdgesvd

#include "petscmat.h" // USES PetscMat
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimendional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ElasticityImplicitLgDeform::ElasticityImplicitLgDeform(void) :
  _dtm1(-1.0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ElasticityImplicitLgDeform::~ElasticityImplicitLgDeform(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::ElasticityImplicitLgDeform::deallocate(void)
{ // deallocate
  IntegratorElasticityLgDeform::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::ElasticityImplicitLgDeform::timeStep(const PylithScalar dt)
{ // timeStep
  if (_dt != -1.0)
    _dtm1 = _dt;
  else
    _dtm1 = dt;
  _dt = dt;
  if (0 != _material)
    _material->timeStep(_dt);
} // timeStep

// ----------------------------------------------------------------------
// Get stable time step for advancing from time t to time t+dt.
PylithScalar
pylith::feassemble::ElasticityImplicitLgDeform::stableTimeStep(const topology::Mesh& mesh) const
{ // stableTimeStep
  assert(0 != _material);
  return _material->stableTimeStepImplicit(mesh);
} // stableTimeStep

// ----------------------------------------------------------------------
// Integrate constributions to residual term (r) for operator.
void
pylith::feassemble::ElasticityImplicitLgDeform::integrateResidual(
			  const topology::Field<topology::Mesh>& residual,
			  const PylithScalar t,
			  topology::SolutionFields* const fields)
{ // integrateResidual
  /// Member prototype for _elasticityResidualXD()
  typedef void (pylith::feassemble::ElasticityImplicitLgDeform::*elasticityResidual_fn_type)
    (const scalar_array&, const scalar_array&);
  
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(0 != _logger);
  assert(0 != fields);

  const int setupEvent = _logger->eventId("ElIR setup");
  const int geometryEvent = _logger->eventId("ElIR geometry");
  const int computeEvent = _logger->eventId("ElIR compute");
  const int restrictEvent = _logger->eventId("ElIR restrict");
  const int stateVarsEvent = _logger->eventId("ElIR stateVars");
  const int stressEvent = _logger->eventId("ElIR stress");
  const int updateEvent = _logger->eventId("ElIR update");

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const int tensorSize = _material->tensorSize();
  if (cellDim != spaceDim)
    throw std::logic_error("Integration for cells with spatial dimensions "
			   "different than the spatial dimension of the "
			   "domain not implemented yet.");

  // Set variables dependent on dimension of cell
  totalStrain_fn_type calcTotalStrainFn;
  elasticityResidual_fn_type elasticityResidualFn;
  if (1 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicitLgDeform::_elasticityResidual1D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicitLgDeform::_elasticityResidual2D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicitLgDeform::_elasticityResidual3D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain3D;
  } else
    assert(0);

  // Allocate vectors for cell values.
  scalar_array dispTCell(numBasis*spaceDim);
  scalar_array dispTIncrCell(numBasis*spaceDim);
  scalar_array dispTpdtCell(numBasis*spaceDim);
  scalar_array deformCell(numQuadPts*spaceDim*spaceDim);
  scalar_array strainCell(numQuadPts*tensorSize);
  scalar_array gravVec(spaceDim);
  scalar_array quadPtsGlobal(numQuadPts*spaceDim);

  // Get cell information
  DM              dmMesh = fields->mesh().dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMComplexGetStratumIS(dmMesh, "material-id", _material->id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

  // Get sections
  topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
  PetscSection dispTSection = dispT.petscSection();
  Vec          dispTVec     = dispT.localVector();
  assert(dispTSection);assert(dispTVec);

  topology::Field<topology::Mesh>& dispTIncr = fields->get("dispIncr(t->t+dt)");
  PetscSection dispTIncrSection = dispTIncr.petscSection();
  Vec          dispTIncrVec     = dispTIncr.localVector();
  assert(dispTIncrSection);assert(dispTIncrVec);

  PetscSection residualSection = residual.petscSection();
  Vec          residualVec     = residual.localVector();

  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);

  assert(0 != _normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar gravityScale = 
    _normalizer->pressureScale() / (_normalizer->lengthScale() *
				    _normalizer->densityScale());

  _logger->eventEnd(setupEvent);
  _logger->eventBegin(computeEvent);

  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(cell);
#else
    const PetscScalar *coords;
    PetscInt           coordsSize;
    err = DMComplexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    _quadrature->computeGeometry(coordinatesCell, cell);
    err = DMComplexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
#endif

    // Get state variables for cell.
    _material->retrievePropsAndVars(cell);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    const PetscScalar *dispTArray, *dispTIncrArray;
    PetscInt           dispTSize,   dispTIncrSize;
    err = DMComplexVecGetClosure(dmMesh, dispTSection, dispTVec, cell, &dispTSize, &dispTArray);CHECK_PETSC_ERROR(err);
    err = DMComplexVecGetClosure(dmMesh, dispTIncrSection, dispTIncrVec, cell, &dispTIncrSize, &dispTIncrArray);CHECK_PETSC_ERROR(err);
    assert(dispTSize == dispTIncrSize);

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& basisDeriv = _quadrature->basisDeriv();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();
    const scalar_array& quadPtsNondim = _quadrature->quadPts();

    // Compute current estimate of displacement at time t+dt using
    // solution increment.
    for(PetscInt i = 0; i < dispTSize; ++i) {dispTpdtCell[i] = dispTArray[i] + dispTIncrArray[i];}
    err = DMComplexVecRestoreClosure(dmMesh, dispTSection, dispTVec, cell, &dispTSize, &dispTArray);CHECK_PETSC_ERROR(err);
    err = DMComplexVecRestoreClosure(dmMesh, dispTIncrSection, dispTIncrVec, cell, &dispTIncrSize, &dispTIncrArray);CHECK_PETSC_ERROR(err);

    // Compute body force vector if gravity is being used.
    if (0 != _gravityField) {
      const spatialdata::geocoords::CoordSys* cs = fields->mesh().coordsys();
      assert(0 != cs);
      
      // Get density at quadrature points for this cell
      const scalar_array& density = _material->calcDensity();

      quadPtsGlobal = quadPtsNondim;
      _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
				  lengthScale);

      // Compute action for element body forces
      spatialdata::spatialdb::SpatialDB* db = _gravityField;
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const int err = db->query(&gravVec[0], gravVec.size(),
					     &quadPtsGlobal[0], spaceDim, cs);
	if (err)
	  throw std::runtime_error("Unable to get gravity vector for point.");
	_normalizer->nondimensionalize(&gravVec[0], gravVec.size(), 
				       gravityScale);
	const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad];
	for (int iBasis=0, iQ=iQuad*numBasis;
	     iBasis < numBasis; ++iBasis) {
	  const PylithScalar valI = wt*basis[iQ+iBasis];
	  for (int iDim=0; iDim < spaceDim; ++iDim) {
	    _cellVector[iBasis*spaceDim+iDim] += valI*gravVec[iDim];
	  } // for
	} // for
      } // for
      PetscLogFlops(numQuadPts*(2+numBasis*(1+2*spaceDim)));
    } // if

    // Compute B(transpose) * sigma, first computing deformation
    // tensor and strains
    _calcDeformation(&deformCell, basisDeriv, coordinatesCell, dispTpdtCell,
		     numBasis, numQuadPts, spaceDim);
    calcTotalStrainFn(&strainCell, deformCell, numQuadPts);
    const scalar_array& stressCell = _material->calcStress(strainCell, true);

    CALL_MEMBER_FN(*this, elasticityResidualFn)(stressCell, dispTpdtCell);

#if 0 // DEBUGGING
    std::cout << "Updating residual for cell " << *c_iter << std::endl;
    for(int i = 0; i < _quadrature->spaceDim() * _quadrature->numBasis(); ++i) {
      std::cout << "  v["<<i<<"]: " << _cellVector[i] << std::endl;
    }
#endif
    // Assemble cell contribution into field
    err = DMComplexVecSetClosure(dmMesh, residualSection, residualVec, cell, &_cellVector[0], ADD_VALUES);CHECK_PETSC_ERROR(err);
  } // for
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);
  
  _logger->eventEnd(computeEvent);
} // integrateResidual

// ----------------------------------------------------------------------
// Compute stiffness matrix.
void
pylith::feassemble::ElasticityImplicitLgDeform::integrateJacobian(
					topology::Jacobian* jacobian,
					const PylithScalar t,
					topology::SolutionFields* fields)
{ // integrateJacobian
  /// Member prototype for _elasticityJacobianXD()
  typedef void (pylith::feassemble::ElasticityImplicitLgDeform::*elasticityJacobian_fn_type)
    (const scalar_array&, const scalar_array&, const scalar_array&);

  assert(0 != _quadrature);
  assert(0 != _material);
  assert(0 != _logger);
  assert(0 != jacobian);
  assert(0 != fields);

  const int setupEvent = _logger->eventId("ElIJ setup");
  const int geometryEvent = _logger->eventId("ElIJ geometry");
  const int computeEvent = _logger->eventId("ElIJ compute");
  const int restrictEvent = _logger->eventId("ElIJ restrict");
  const int stateVarsEvent = _logger->eventId("ElIJ stateVars");
  const int updateEvent = _logger->eventId("ElIJ update");

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const int tensorSize = _material->tensorSize();
  if (cellDim != spaceDim)
    throw std::logic_error("Don't know how to integrate elasticity " \
			   "contribution to Jacobian matrix for cells with " \
			   "different dimensions than the spatial dimension.");

  // Set variables dependent on dimension of cell
  totalStrain_fn_type calcTotalStrainFn;
  elasticityJacobian_fn_type elasticityJacobianFn;
  if (1 == cellDim) {
    elasticityJacobianFn = 
      &pylith::feassemble::ElasticityImplicitLgDeform::_elasticityJacobian1D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    elasticityJacobianFn = 
      &pylith::feassemble::ElasticityImplicitLgDeform::_elasticityJacobian2D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    elasticityJacobianFn = 
      &pylith::feassemble::ElasticityImplicitLgDeform::_elasticityJacobian3D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticityLgDeform::_calcTotalStrain3D;
  } else
    assert(0);

  // Allocate vector for total strain
  scalar_array dispTCell(numBasis*spaceDim);
  scalar_array dispTIncrCell(numBasis*spaceDim);
  scalar_array dispTpdtCell(numBasis*spaceDim);
  scalar_array deformCell(numQuadPts*spaceDim*spaceDim);
  scalar_array strainCell(numQuadPts*tensorSize);
  scalar_array stressCell(numQuadPts*tensorSize);

  // Get cell information
  DM              dmMesh = fields->mesh().dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  const int materialId = _material->id();
  err = DMComplexGetStratumIS(dmMesh, "material-id", materialId, &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

  // Get sections
  topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
  PetscSection dispTSection = dispT.petscSection();
  Vec          dispTVec     = dispT.localVector();
  assert(dispTSection);assert(dispTVec);

  topology::Field<topology::Mesh>& dispTIncr = fields->get("dispIncr(t->t+dt)");
  PetscSection dispTIncrSection = dispTIncr.petscSection();
  Vec          dispTIncrVec     = dispTIncr.localVector();
  assert(dispTIncrSection);assert(dispTIncrVec);

  // Get sparse matrix
  const PetscMat jacobianMat = jacobian->matrix();
  assert(0 != jacobianMat);

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  assert(dt > 0);

  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);

  _logger->eventEnd(setupEvent);
  _logger->eventBegin(computeEvent);

  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(cell);
#else
    const PetscScalar *coords;
    PetscInt           coordsSize;
    err = DMComplexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    _quadrature->computeGeometry(coordinatesCell, cell);
    err = DMComplexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
#endif

    // Get state variables for cell.
    _material->retrievePropsAndVars(cell);

    // Reset element matrix to zero
    _resetCellMatrix();

    // Restrict input fields to cell
    const PetscScalar *dispTArray, *dispTIncrArray;
    PetscInt           dispTSize,   dispTIncrSize;
    err = DMComplexVecGetClosure(dmMesh, dispTSection, dispTVec, cell, &dispTSize, &dispTArray);CHECK_PETSC_ERROR(err);
    err = DMComplexVecGetClosure(dmMesh, dispTIncrSection, dispTIncrVec, cell, &dispTIncrSize, &dispTIncrArray);CHECK_PETSC_ERROR(err);
    assert(dispTSize == dispTIncrSize);

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& basisDeriv = _quadrature->basisDeriv();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();

    // Compute current estimate of displacement at time t+dt using
    // solution increment.
    for(PetscInt i = 0; i < dispTSize; ++i) {dispTpdtCell[i] = dispTArray[i] + dispTIncrArray[i];}
    err = DMComplexVecRestoreClosure(dmMesh, dispTSection, dispTVec, cell, &dispTSize, &dispTArray);CHECK_PETSC_ERROR(err);
    err = DMComplexVecRestoreClosure(dmMesh, dispTIncrSection, dispTIncrVec, cell, &dispTIncrSize, &dispTIncrArray);CHECK_PETSC_ERROR(err);
      
    // Compute deformation tensor, strains, and stresses
    _calcDeformation(&deformCell, basisDeriv, coordinatesCell, dispTpdtCell,
		     numBasis, numQuadPts, spaceDim);
    calcTotalStrainFn(&strainCell, deformCell, numQuadPts);

    // Get "elasticity" matrix at quadrature points for this cell
    const scalar_array& elasticConsts = 
      _material->calcDerivElastic(strainCell);

    // Get Second Priola-Kirchoff stress tensor
    const scalar_array& stressCell = _material->calcStress(strainCell, true);

    CALL_MEMBER_FN(*this, elasticityJacobianFn)(elasticConsts, stressCell,
						dispTpdtCell);

    if (_quadrature->checkConditioning()) {
      int n = numBasis*spaceDim;
      int lwork = 5*n;
      int idummy = 0;
      int lierr = 0;
      PylithScalar *elemMat = new PylithScalar[n*n];
      PylithScalar *svalues = new PylithScalar[n];
      PylithScalar *work    = new PylithScalar[lwork];
      PylithScalar minSV = 0;
      PylithScalar maxSV = 0;
      PylithScalar sdummy = 0;

      const int n2 = n*n;
      for (int i = 0; i < n2; ++i)
	elemMat[i] = _cellMatrix[i];
      lapack_dgesvd("N", "N", &n, &n, elemMat, &n, svalues, 
		    &sdummy, &idummy, &sdummy, &idummy, work,
		    &lwork, &lierr);
      if (lierr)
	throw std::runtime_error("Lapack SVD failed");
      minSV = svalues[n-7];
      maxSV = svalues[0];
      std::cout << "Element " << cell << std::endl;
      for(int i = 0; i < n; ++i)
	std::cout << "    sV["<<i<<"] = " << svalues[i] << std::endl;
      std::cout << "  kappa(elemMat) = " << maxSV/minSV << std::endl;
      delete [] elemMat;
      delete [] svalues;
      delete [] work;
    } // if

    // Assemble cell contribution into PETSc matrix.
    err = DMComplexMatSetClosure(dmMesh, dispTSection, PETSC_NULL, jacobianMat, cell, &_cellMatrix[0], ADD_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");
  } // for

  _logger->eventEnd(computeEvent);

  _needNewJacobian = false;
  _material->resetNeedNewJacobian();
} // integrateJacobian


// End of file 
