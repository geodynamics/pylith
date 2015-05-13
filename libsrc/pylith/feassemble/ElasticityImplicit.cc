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

#include "ElasticityImplicit.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature
#include "CellGeometry.hh" // USES CellGeometry

#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN
#include "pylith/utils/lapack.h" // USES LAPACKdgesvd

#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimendional

#include "petscmat.h" // USES PetscMat

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ElasticityImplicit::ElasticityImplicit(void) :
  _dtm1(-1.0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ElasticityImplicit::~ElasticityImplicit(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::ElasticityImplicit::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  IntegratorElasticity::deallocate();

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::ElasticityImplicit::timeStep(const PylithScalar dt)
{ // timeStep
  PYLITH_METHOD_BEGIN;

  if (_dt != -1.0)
    _dtm1 = _dt;
  else
    _dtm1 = dt;
  _dt = dt;
  if (_material)
    _material->timeStep(_dt);

  PYLITH_METHOD_END;
} // timeStep

// ----------------------------------------------------------------------
// Get stable time step for advancing from time t to time t+dt.
PylithScalar
pylith::feassemble::ElasticityImplicit::stableTimeStep(const topology::Mesh& mesh) const
{ // stableTimeStep
  PYLITH_METHOD_BEGIN;

  assert(_material);
  PYLITH_METHOD_RETURN(_material->stableTimeStepImplicit(mesh));
} // stableTimeStep

// ----------------------------------------------------------------------
void
pylith::feassemble::ElasticityImplicit::integrateResidual(const topology::Field& residual,
							  const PylithScalar t,
							  topology::SolutionFields* const fields)
{ // integrateResidual
  PYLITH_METHOD_BEGIN;

  /// Member prototype for _elasticityResidualXD()
  typedef void (pylith::feassemble::ElasticityImplicit::*elasticityResidual_fn_type)
    (const scalar_array&);
  
  assert(_quadrature);
  assert(_material);
  assert(_logger);
  assert(fields);

  const int setupEvent = _logger->eventId("ElIR setup");
  const int computeEvent = _logger->eventId("ElIR compute");
#if defined(DETAILED_EVENT_LOGGING)
  const int geometryEvent = _logger->eventId("ElIR geometry");
  const int restrictEvent = _logger->eventId("ElIR restrict");
  const int stateVarsEvent = _logger->eventId("ElIR stateVars");
  const int stressEvent = _logger->eventId("ElIR stress");
  const int updateEvent = _logger->eventId("ElIR update");
#endif

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == size_t(numQuadPts));
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
  if (2 == cellDim) {
    elasticityResidualFn = &pylith::feassemble::ElasticityImplicit::_elasticityResidual2D;
    calcTotalStrainFn = &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    elasticityResidualFn = &pylith::feassemble::ElasticityImplicit::_elasticityResidual3D;
    calcTotalStrainFn = &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else {
    assert(false);
    throw std::logic_error("Unsupported cell dimension in ElasticityImplicit::integrateResidual().");
  } // if/else		   

  // Allocate vectors for cell values.
  scalar_array dispTpdtCell(numBasis*spaceDim);
  scalar_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;
  scalar_array gravVec(spaceDim);
  scalar_array quadPtsGlobal(numQuadPts*spaceDim);

  // Get cell information
  PetscDM dmMesh = fields->mesh().dmMesh();assert(dmMesh);
  assert(_materialIS);
  const PetscInt* cells = _materialIS->points();
  const PetscInt numCells = _materialIS->size();

  // Setup field visitors.
  scalar_array dispCell(numBasis*spaceDim);
  topology::VecVisitorMesh dispVisitor(fields->get("disp(t)"), "displacement");
  dispVisitor.optimizeClosure();

  scalar_array dispIncrCell(numBasis*spaceDim);
  topology::VecVisitorMesh dispIncrVisitor(fields->get("dispIncr(t->t+dt)"), "displacement");
  dispIncrVisitor.optimizeClosure();

  topology::VecVisitorMesh residualVisitor(residual, "displacement");
  residualVisitor.optimizeClosure();

  scalar_array coordsCell(numBasis*spaceDim); // :KLUDGE: numBasis to numCorners after switching to higher order
  topology::CoordsVisitor coordsVisitor(dmMesh);

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar gravityScale = _normalizer->pressureScale() / (_normalizer->lengthScale() * _normalizer->densityScale());

  _logger->eventEnd(setupEvent);
  _logger->eventBegin(computeEvent);

  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, cell);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), cell);

    // Get state variables for cell.
    _material->retrievePropsAndVars(cell);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    dispVisitor.getClosure(&dispCell, cell);
    dispIncrVisitor.getClosure(&dispIncrCell, cell);

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& basisDeriv = _quadrature->basisDeriv();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();
    const scalar_array& quadPtsNondim = _quadrature->quadPts();

    // Compute current estimate of displacement at time t+dt using solution increment.
    for(PetscInt i = 0, dispSize = dispCell.size(); i < dispSize; ++i) {
      dispTpdtCell[i] = dispCell[i] + dispIncrCell[i];
    } // for

    // Compute body force vector if gravity is being used.
    if (_gravityField) {
      const spatialdata::geocoords::CoordSys* cs = fields->mesh().coordsys();assert(cs);

      // Get density at quadrature points for this cell
      const scalar_array& density = _material->calcDensity();

      quadPtsGlobal = quadPtsNondim;
      _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(), lengthScale);

      // Compute action for element body forces
      spatialdata::spatialdb::SpatialDB* db = _gravityField;
      for (int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
        const int err = db->query(&gravVec[0], gravVec.size(), &quadPtsGlobal[0], spaceDim, cs);
        if (err) {
	  throw std::runtime_error("Unable to get gravity vector for point.");
	} // if
        _normalizer->nondimensionalize(&gravVec[0], gravVec.size(), gravityScale);
        const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad];
        for (int iBasis = 0, iQ = iQuad * numBasis; iBasis < numBasis; ++iBasis) {
          const PylithScalar valI = wt * basis[iQ + iBasis];
          for (int iDim = 0; iDim < spaceDim; ++iDim) {
            _cellVector[iBasis * spaceDim + iDim] += valI * gravVec[iDim];
          } // for
        } // for
      } // for
      PetscLogFlops(numQuadPts * (2 + numBasis * (1 + 2 * spaceDim)));
    } // if

    // residualSection->view("After gravity contribution");
    // Compute B(transpose) * sigma, first computing strains
    calcTotalStrainFn(&strainCell, basisDeriv, &dispTpdtCell[0], numBasis, spaceDim, numQuadPts);
    const scalar_array& stressCell = _material->calcStress(strainCell, true);

    CALL_MEMBER_FN(*this, elasticityResidualFn)(stressCell);

#if 0 // DEBUGGING
    std::cout << "Updating residual for cell " << cell << std::endl;
    for(PetscInt i = 0; i < spaceDim*numBasis; ++i) {
      std::cout << "  v["<<i<<"]: " << _cellVector[i] << std::endl;
    } // for
#endif
    // Assemble cell contribution into field
    residualVisitor.setClosure(&_cellVector[0], _cellVector.size(), cell, ADD_VALUES);
  } // for

  _logger->eventEnd(computeEvent);

  PYLITH_METHOD_END;
} // integrateResidual

// ----------------------------------------------------------------------
// Compute stiffness matrix.
void
pylith::feassemble::ElasticityImplicit::integrateJacobian(topology::Jacobian* jacobian,
							  const PylithScalar t,
							  topology::SolutionFields* fields)
{ // integrateJacobian
  PYLITH_METHOD_BEGIN;

  /// Member prototype for _elasticityJacobianXD()
  typedef void (pylith::feassemble::ElasticityImplicit::*elasticityJacobian_fn_type)
    (const scalar_array&);

  assert(_quadrature);
  assert(_material);
  assert(_logger);
  assert(jacobian);
  assert(fields);

  const int setupEvent = _logger->eventId("ElIJ setup");
  const int computeEvent = _logger->eventId("ElIJ compute");

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == size_t(numQuadPts));
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
  if (2 == cellDim) {
    elasticityJacobianFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityJacobian2D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    elasticityJacobianFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityJacobian3D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else {
    assert(false);
    throw std::logic_error("Unsupported cell dimension in ElasticityImplicit::integrateJacobian().");
  } // if/else

  // Allocate vector for total strain
  scalar_array dispTpdtCell(numBasis*spaceDim);
  scalar_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;

  // Get cell information
  PetscDM dmMesh = fields->mesh().dmMesh();assert(dmMesh);
  assert(_materialIS);
  const PetscInt* cells = _materialIS->points();
  const PetscInt numCells = _materialIS->size();

  // Setup field visitors.
  scalar_array dispCell(numBasis*spaceDim);
  topology::VecVisitorMesh dispVisitor(fields->get("disp(t)"), "displacement");
  dispVisitor.optimizeClosure();

  scalar_array dispIncrCell(numBasis*spaceDim);
  topology::VecVisitorMesh dispIncrVisitor(fields->get("dispIncr(t->t+dt)"), "displacement");
  dispIncrVisitor.optimizeClosure();

  scalar_array coordsCell(numBasis*spaceDim); // :KLUDGE: numBasis to numCorners after switching to higher order
  topology::CoordsVisitor coordsVisitor(dmMesh);

  // Get sparse matrix
  const PetscMat jacobianMat = jacobian->matrix();assert(jacobianMat);
  topology::MatVisitorMesh jacobianVisitor(jacobianMat, fields->get("disp(t)"));

  // Get parameters used in integration.
  const PylithScalar dt = _dt;
  assert(dt > 0);

  _logger->eventEnd(setupEvent);
  _logger->eventBegin(computeEvent);

  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];

    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, cell);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), cell);

    // Get physical properties and state variables for cell.
    _material->retrievePropsAndVars(cell);

    // Reset element matrix to zero
    _resetCellMatrix();

    // Restrict input fields to cell
    dispVisitor.getClosure(&dispCell, cell);
    dispIncrVisitor.getClosure(&dispIncrCell, cell);

    // Get cell geometry information that depends on cell
    const scalar_array& basisDeriv = _quadrature->basisDeriv();

    // Compute current estimate of displacement at time t+dt using solution increment.
    for(PetscInt i = 0, dispSize = dispCell.size(); i < dispSize; ++i) {
      dispTpdtCell[i] = dispCell[i] + dispIncrCell[i];
    } // for
      
    // Compute strains
    calcTotalStrainFn(&strainCell, basisDeriv, &dispTpdtCell[0], numBasis, spaceDim, numQuadPts);
      
    // Get "elasticity" matrix at quadrature points for this cell
    const scalar_array& elasticConsts = _material->calcDerivElastic(strainCell);

    CALL_MEMBER_FN(*this, elasticityJacobianFn)(elasticConsts);

    if (_quadrature->checkConditioning()) {
      int n = numBasis*spaceDim;
      int lwork = 5*n;
      int idummy = 0;
      int lierr = 0;
      PylithScalar *elemMat = new PylithScalar[n*n];
      PylithScalar *svalues = new PylithScalar[n];
      PylithScalar *work    = new PylithScalar[lwork];
#if 0
      PylithScalar minSV = 0;
      PylithScalar maxSV = 0;
#endif
      PylithScalar sdummy = 0;

      const int n2 = n*n;
      for (int i = 0; i < n2; ++i)
        elemMat[i] = _cellMatrix[i];
      lapack_dgesvd("N", "N", &n, &n, elemMat, &n, svalues, 
		    &sdummy, &idummy, &sdummy, &idummy, work,
		    &lwork, &lierr);
      if (lierr)
        throw std::runtime_error("Lapack SVD failed");
#if 0
      minSV = svalues[n-7];
      maxSV = svalues[0];
      std::cout << "Element " << cell << std::endl;
      for(int i = 0; i < n; ++i)
        std::cout << "    sV["<<i<<"] = " << svalues[i] << std::endl;
      std::cout << "  kappa(elemMat) = " << maxSV/minSV << std::endl;
#endif
      delete [] elemMat;
      delete [] svalues;
      delete [] work;
    } // if

    // Assemble cell contribution into PETSc matrix.
    //   Notice that we are using the default sections
    jacobianVisitor.setClosure(&_cellMatrix[0], _cellMatrix.size(), cell, ADD_VALUES);
  } // for

  _needNewJacobian = false;
  _material->resetNeedNewJacobian();

  _logger->eventEnd(computeEvent);

  PYLITH_METHOD_END;
} // integrateJacobian


// End of file 
