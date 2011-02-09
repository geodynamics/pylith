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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "ElasticityExplicitTri3.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature

#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "pylith/utils/array.hh" // USES double_array
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN
#include "pylith/utils/lapack.h" // USES LAPACKdgesvd

#include "petscmat.h" // USES PetscMat
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimendional

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

//#define DETAILED_EVENT_LOGGING

// ----------------------------------------------------------------------
const int pylith::feassemble::ElasticityExplicitTri3::_spaceDim = 2;
const int pylith::feassemble::ElasticityExplicitTri3::_cellDim = 2;
const int pylith::feassemble::ElasticityExplicitTri3::_tensorSize = 3;
const int pylith::feassemble::ElasticityExplicitTri3::_numBasis = 3;
const int pylith::feassemble::ElasticityExplicitTri3::_numQuadPts = 1;

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ElasticityExplicitTri3::ElasticityExplicitTri3(void) :
  _dtm1(-1.0)
{ // constructor
  _basisDerivArray.resize(_numQuadPts*_numBasis*_spaceDim);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ElasticityExplicitTri3::~ElasticityExplicitTri3(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::ElasticityExplicitTri3::deallocate(void)
{ // deallocate
  IntegratorElasticity::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::ElasticityExplicitTri3::timeStep(const double dt)
{ // timeStep
  if (_dt != -1.0)
    _dtm1 = _dt;
  else
    _dtm1 = dt;
  _dt = dt;
  assert(_dt == _dtm1); // For now, don't allow variable time step
  if (0 != _material)
    _material->timeStep(_dt);
} // timeStep

// ----------------------------------------------------------------------
// Set flag for setting constraints for total field solution or
// incremental field solution.
void
pylith::feassemble::ElasticityExplicitTri3::useSolnIncr(const bool flag)
{ // useSolnIncr
  _material->useElasticBehavior(false);
  if (!flag)
    throw std::logic_error("Non-incremental solution not supported for "
			   "explicit time integration of elasticity "
			   "equation.");
} // useSolnIncr

// ----------------------------------------------------------------------
// Integrate constributions to residual term (r) for operator.
void
pylith::feassemble::ElasticityExplicitTri3::integrateResidual(
			  const topology::Field<topology::Mesh>& residual,
			  const double t,
			  topology::SolutionFields* const fields)
{ // integrateResidual
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
  assert(_quadrature->numQuadPts() == _numQuadPts);
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == _numQuadPts);
  assert(_quadrature->numBasis() == _numBasis);
  assert(_quadrature->spaceDim() == _spaceDim);
  assert(_quadrature->cellDim() == _cellDim);
  assert(_material->tensorSize() == _tensorSize);
  const int spaceDim = _spaceDim;
  const int cellDim = _cellDim;
  const int tensorSize = _tensorSize;
  const int numBasis = _numBasis;
  const int numQuadPts = _numQuadPts;
  /** :TODO:
   *
   * If cellDim and spaceDim are different, we need to transform
   * displacements into cellDim, compute action, and transform result
   * back into spaceDim. We get this information from the Jacobian and
   * inverse of the Jacobian.
   */
  if (cellDim != spaceDim)
    throw std::logic_error("Integration for cells with spatial dimensions "
			   "different than the spatial dimension of the "
			   "domain not implemented yet.");

  // Allocate vectors for cell values.
  double_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;
  double_array gravVec(spaceDim);
  double_array quadPtsGlobal(numQuadPts*spaceDim);

  // Get cell information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const int materialId = _material->id();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  double_array accCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& accSection = 
    fields->get("acceleration(t)").section();
  assert(!accSection.isNull());
  RestrictVisitor accVisitor(*accSection, accCell.size(), &accCell[0]);

  double_array dispCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& dispSection = 
    fields->get("disp(t)").section();
  assert(!dispSection.isNull());
  RestrictVisitor dispVisitor(*dispSection, dispCell.size(), &dispCell[0]);

  const ALE::Obj<RealSection>& residualSection = residual.section();
  UpdateAddVisitor residualVisitor(*residualSection, &_cellVector[0]);
  
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
  const double gravityScale = 
    _normalizer->pressureScale() / (_normalizer->lengthScale() *
				    _normalizer->densityScale());

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif

    // Compute geometry information for current cell
    coordsVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);
    const double area = _area(coordinatesCell);
    assert(area > 0.0);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(stateVarsEvent);
#endif

    // Get state variables for cell.
    _material->retrievePropsAndVars(*c_iter);

    #if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(stateVarsEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    accVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, accVisitor);
    
    dispVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispVisitor);
    
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    const double_array& density = _material->calcDensity();
    assert(density.size() == 1);

    // Compute body force vector if gravity is being used.
    if (0 != _gravityField) {
      const spatialdata::geocoords::CoordSys* cs = fields->mesh().coordsys();
      assert(0 != cs);

      quadPtsGlobal = 0.0;
      for (int iBasis=0; iBasis < numBasis; ++iBasis)
        for (int iDim=0; iDim < spaceDim; ++iDim)
          quadPtsGlobal[iDim] += 
	    coordinatesCell[iBasis*spaceDim+iDim] / numBasis;
      _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
          lengthScale);

      // Compute action for element body forces
      const int err = _gravityField->query(&gravVec[0], gravVec.size(),
        &quadPtsGlobal[0], spaceDim, cs);
      if (err)
        throw std::runtime_error("Unable to get gravity vector for point.");
      _normalizer->nondimensionalize(&gravVec[0], gravVec.size(),
          gravityScale);
      const double wtVertex = density[0] * area / 3.0;
      for (int iBasis=0; iBasis < numBasis; ++iBasis)
        for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis * spaceDim + iDim] += wtVertex * gravVec[iDim];
      PetscLogFlops(numBasis*spaceDim*2 + numBasis*spaceDim*2);
    } // if

    // Compute action for inertial terms
    const double wtVertex = density[0] * area / 9.0;
    for (int iBasis = 0; iBasis < numBasis; ++iBasis)
      for (int jBasis = 0; jBasis < numBasis; ++jBasis)
        for (int iDim = 0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis*spaceDim+iDim] -= 
	      wtVertex * accCell[jBasis*spaceDim+iDim];

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(3 + numBasis*numBasis*spaceDim*2);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(stressEvent);
#endif

    // Compute B(transpose) * sigma, first computing strains
    const double x0 = coordinatesCell[0];
    const double y0 = coordinatesCell[1];

    const double x1 = coordinatesCell[2];
    const double y1 = coordinatesCell[3];

    const double x2 = coordinatesCell[4];
    const double y2 = coordinatesCell[5];

    const double scaleB = 2.0 * area;
    const double b0 = (y1 - y2) / scaleB;
    const double c0 = (x2 - x1) / scaleB;

    const double b1 = (y2 - y0) / scaleB;
    const double c1 = (x0 - x2) / scaleB;

    const double b2 = (y0 - y1) / scaleB;
    const double c2 = (x1 - x0) / scaleB;

    assert(strainCell.size() == 3);
    strainCell[0] = b2*dispCell[4] + b1*dispCell[2] + b0*dispCell[0];
    
    strainCell[1] = c2*dispCell[5] + c1*dispCell[3] + c0*dispCell[1];

    strainCell[2] = (b2*dispCell[5] + c2*dispCell[4] + b1*dispCell[3] + 
		     c1*dispCell[2] + b0*dispCell[1] + c0*dispCell[0]) / 2.0;


    const double_array& stressCell = _material->calcStress(strainCell, true);

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(34);
    _logger->eventEnd(stressEvent);
    _logger->eventBegin(computeEvent);
#endif

    assert(_cellVector.size() == 6);
    assert(stressCell.size() == 3);
    _cellVector[0] -= (c0*stressCell[2] + b0*stressCell[0]) * area;
    _cellVector[1] -= (b0*stressCell[2] + c0*stressCell[1]) * area;
    _cellVector[2] -= (c1*stressCell[2] + b1*stressCell[0]) * area;
    _cellVector[3] -= (b1*stressCell[2] + c1*stressCell[1]) * area;
    _cellVector[4] -= (c2*stressCell[2] + b2*stressCell[0]) * area;
    _cellVector[5] -= (b2*stressCell[2] + c2*stressCell[1]) * area;

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(30);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Assemble cell contribution into field
    residualVisitor.clear();
    sieveMesh->updateClosure(*c_iter, residualVisitor);
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(cells->size()*(3 + numBasis*numBasis*spaceDim*2 + 34+30));
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate constributions to residual term (r) for operator.
void
pylith::feassemble::ElasticityExplicitTri3::integrateResidualLumped(
        const topology::Field<topology::Mesh>& residual,
        const double t,
        topology::SolutionFields* const fields)
{ // integrateResidualLumped
  /// Member prototype for _elasticityResidualXD()
  typedef void (pylith::feassemble::ElasticityExplicitTri3::*elasticityResidual_fn_type)
    (const double_array&);

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
  assert(_quadrature->numQuadPts() == _numQuadPts);
  assert(_quadrature->numBasis() == _numBasis);
  assert(_quadrature->spaceDim() == _spaceDim);
  assert(_quadrature->cellDim() == _cellDim);
  assert(_material->tensorSize() == _tensorSize);
  const int spaceDim = _spaceDim;
  const int cellDim = _cellDim;
  const int tensorSize = _tensorSize;
  const int numBasis = _numBasis;
  const int numQuadPts = _numQuadPts;
  /** :TODO:
   *
   * If cellDim and spaceDim are different, we need to transform
   * displacements into cellDim, compute action, and transform result
   * back into spaceDim. We get this information from the Jacobian and
   * inverse of the Jacobian.
   */
  if (cellDim != spaceDim)
    throw std::logic_error("Integration for cells with spatial dimensions "
         "different than the spatial dimension of the "
         "domain not implemented yet.");

  // Allocate vectors for cell values.
  double_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;
  double_array gravVec(spaceDim);
  double_array quadPtsGlobal(numQuadPts*spaceDim);

  // Get cell information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const int materialId = _material->id();
  const ALE::Obj<SieveMesh::label_sequence>& cells =
    sieveMesh->getLabelStratum("material-id", materialId);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());

  // Get sections
  double_array accCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& accSection = 
    fields->get("acceleration(t)").section();
  assert(!accSection.isNull());
  RestrictVisitor accVisitor(*accSection, accCell.size(), &accCell[0]);

  double_array dispCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& dispSection = 
    fields->get("disp(t)").section();
  assert(!dispSection.isNull());
  RestrictVisitor dispVisitor(*dispSection, dispCell.size(), &dispCell[0]);

  const ALE::Obj<RealSection>& residualSection = residual.section();
  UpdateAddVisitor residualVisitor(*residualSection, &_cellVector[0]);

  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates =
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
  const double gravityScale =
    _normalizer->pressureScale() / (_normalizer->lengthScale() *
            _normalizer->densityScale());

  // Get parameters used in integration.
  double_array valuesIJ(numBasis);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Restrict input fields to cell
#if 1
    coordsVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);

    accVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, accVisitor);

    dispVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispVisitor);
#else
    coordsVisitor.clear();
    sieve->orientedConeOpt(*c_iter, coordsVisitor, numBasis, spaceDim);

    accVisitor.clear();
    sieve->orientedConeOpt(*c_iter, accVisitor, numBasis, spaceDim);

    dispVisitor.clear();
    sieve->orientedConeOpt(*c_iter, dispVisitor, numBasis, spaceDim);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(geometryEvent);
#endif

    // Compute geometry information for current cell
    const double area = _area(coordinatesCell);
    assert(area > 0.0);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(stateVarsEvent);
#endif

    // Get state variables for cell.
    _material->retrievePropsAndVars(*c_iter);

    // Get density at quadrature points for this cell
    const double_array& density = _material->calcDensity();

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(stateVarsEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Compute body force vector if gravity is being used.
    if (0 != _gravityField) {
      const spatialdata::geocoords::CoordSys* cs = fields->mesh().coordsys();
      assert(0 != cs);

      quadPtsGlobal = 0.0;
      for (int iBasis=0; iBasis < numBasis; ++iBasis)
        for (int iDim=0; iDim < spaceDim; ++iDim)
          quadPtsGlobal[iDim] += 
	    coordinatesCell[iBasis*spaceDim+iDim] / numBasis;
      _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
          lengthScale);

      // Compute action for element body forces
      const int err = _gravityField->query(&gravVec[0], gravVec.size(),
        &quadPtsGlobal[0], spaceDim, cs);
      if (err)
        throw std::runtime_error("Unable to get gravity vector for point.");
      _normalizer->nondimensionalize(&gravVec[0], gravVec.size(),
          gravityScale);
      const double wtVertex = density[0] * area / 3.0;
      for (int iBasis=0; iBasis < numBasis; ++iBasis)
        for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis * spaceDim + iDim] += wtVertex * gravVec[iDim];
      PetscLogFlops(numBasis*spaceDim*2 + numBasis*spaceDim*2);
    } // if

    // Compute action for inertial terms
    const double wtVertex = density[0] * area / 3.0;
    _cellVector -= wtVertex * accCell;

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(2 + numBasis*spaceDim*2);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(stressEvent);
#endif

    // Compute B(transpose) * sigma, first computing strains
    const double x0 = coordinatesCell[0];
    const double y0 = coordinatesCell[1];

    const double x1 = coordinatesCell[2];
    const double y1 = coordinatesCell[3];

    const double x2 = coordinatesCell[4];
    const double y2 = coordinatesCell[5];

    const double scaleB = 2.0 * area;
    const double b0 = (y1 - y2) / scaleB;
    const double c0 = (x2 - x1) / scaleB;

    const double b1 = (y2 - y0) / scaleB;
    const double c1 = (x0 - x2) / scaleB;

    const double b2 = (y0 - y1) / scaleB;
    const double c2 = (x1 - x0) / scaleB;

    assert(strainCell.size() == 3);
    strainCell[0] = b2*dispCell[4] + b1*dispCell[2] + b0*dispCell[0];
    
    strainCell[1] = c2*dispCell[5] + c1*dispCell[3] + c0*dispCell[1];

    strainCell[2] = (b2*dispCell[5] + c2*dispCell[4] + b1*dispCell[3] + 
		     c1*dispCell[2] + b0*dispCell[1] + c0*dispCell[0]) / 2.0;

    const double_array& stressCell = _material->calcStress(strainCell, true);

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(34);
    _logger->eventEnd(stressEvent);
    _logger->eventBegin(computeEvent);
#endif

    assert(_cellVector.size() == 6);
    assert(stressCell.size() == 3);
    _cellVector[0] -= (c0*stressCell[2] + b0*stressCell[0]) * area;
    _cellVector[1] -= (b0*stressCell[2] + c0*stressCell[1]) * area;
    _cellVector[2] -= (c1*stressCell[2] + b1*stressCell[0]) * area;
    _cellVector[3] -= (b1*stressCell[2] + c1*stressCell[1]) * area;
    _cellVector[4] -= (c2*stressCell[2] + b2*stressCell[0]) * area;
    _cellVector[5] -= (b2*stressCell[2] + c2*stressCell[1]) * area;

#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(30);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Assemble cell contribution into field
    residualVisitor.clear();
    sieveMesh->updateClosure(*c_iter, residualVisitor);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(cells->size()*(2 + numBasis*spaceDim*2 + 34+30));
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidualLumped

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ElasticityExplicitTri3::integrateJacobian(
					topology::Jacobian* jacobian,
					const double t,
					topology::SolutionFields* fields)
{ // integrateJacobian
  assert(0 != _quadrature);
  assert(0 != _material);
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
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const int tensorSize = _material->tensorSize();
  if (cellDim != spaceDim)
    throw std::logic_error("Don't know how to integrate elasticity " \
			   "contribution to Jacobian matrix for cells with " \
			   "different dimensions than the spatial dimension.");

  // Get cell information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const int materialId = _material->id();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  double_array dispCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& dispSection = 
    fields->get("disp(t)").section();
  assert(!dispSection.isNull());

  // Get sparse matrix
  const PetscMat jacobianMat = jacobian->matrix();
  assert(0 != jacobianMat);

  // Get parameters used in integration.
  const double dt = _dt;
  const double dt2 = dt*dt;
  assert(dt > 0);

  const ALE::Obj<SieveMesh::order_type>& globalOrder = 
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", dispSection);
  assert(!globalOrder.isNull());
  // We would need to request unique points here if we had an interpolated mesh
  IndicesVisitor jacobianVisitor(*dispSection, *globalOrder,
				 (int) pow(sieveMesh->getSieve()->getMaxConeSize(),
					   sieveMesh->depth())*spaceDim);

  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif

    // Compute geometry information for current cell
    coordsVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);
    const double area = _area(coordinatesCell);
    assert(area > 0.0);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(stateVarsEvent);
#endif

    // Get state variables for cell.
    _material->retrievePropsAndVars(*c_iter);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(stateVarsEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Reset element matrix to zero
    _resetCellMatrix();

    // Get material physical properties at quadrature points for this cell
    const double_array& density = _material->calcDensity();
    assert(density.size() == 1);

    // Compute Jacobian for inertial terms
    const double wtVertex = density[0] * area / (9.0 * dt2);
    for (int iBasis = 0; iBasis < numBasis; ++iBasis)
      for (int jBasis = 0; jBasis < numBasis; ++jBasis)
	for (int iDim=0; iDim < spaceDim; ++iDim) {
	  const int iBlock = (iBasis*spaceDim + iDim) * (numBasis*spaceDim);
	  const int jBlock = (jBasis*spaceDim + iDim);
	  _cellMatrix[iBlock+jBlock] += wtVertex;
	} // for
    
#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(numQuadPts*(3+numBasis*numBasis*spaceDim*1));
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif
    
    // Assemble cell contribution into PETSc matrix.
    jacobianVisitor.clear();
    PetscErrorCode err = updateOperator(jacobianMat, *sieveMesh->getSieve(),
					jacobianVisitor, *c_iter,
					&_cellMatrix[0], ADD_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(cells->size()*(3+numBasis*numBasis*spaceDim*1));
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;
  _material->resetNeedNewJacobian();
} // integrateJacobian

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ElasticityExplicitTri3::integrateJacobian(
			    topology::Field<topology::Mesh>* jacobian,
			    const double t,
			    topology::SolutionFields* fields)
{ // integrateJacobian
  assert(0 != _quadrature);
  assert(0 != _material);
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
  assert(_quadrature->numBasis() == _numBasis);
  assert(_quadrature->spaceDim() == _spaceDim);
  assert(_quadrature->cellDim() == _cellDim);
  assert(_material->tensorSize() == _tensorSize);
  const int spaceDim = _spaceDim;
  const int cellDim = _cellDim;
  const int numBasis = _numBasis;
  if (cellDim != spaceDim)
    throw std::logic_error("Don't know how to integrate elasticity " \
			   "contribution to Jacobian matrix for cells with " \
			   "different dimensions than the spatial dimension.");

  // Get cell information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const int materialId = _material->id();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get parameters used in integration.
  const double dt = _dt;
  const double dt2 = dt*dt;
  assert(dt > 0);

  // Get sections
  const ALE::Obj<RealSection>& jacobianSection = jacobian->section();
  assert(!jacobianSection.isNull());
  UpdateAddVisitor jacobianVisitor(*jacobianSection, &_cellVector[0]);

  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif
  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
    coordsVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);
    const double area = _area(coordinatesCell);
    assert(area > 0.0);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(stateVarsEvent);
#endif

    // Get state variables for cell.
    _material->retrievePropsAndVars(*c_iter);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(stateVarsEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Compute Jacobian for inertial terms
    const double_array& density = _material->calcDensity();
    _cellVector = density[0] * area / (3.0 * dt2);
    
#if defined(DETAILED_EVENT_LOGGING)
    PetscLogFlops(3);
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif
    
    // Assemble cell contribution into lumped matrix.
    jacobianVisitor.clear();
    sieveMesh->updateClosure(*c_iter, jacobianVisitor);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  PetscLogFlops(cells->size()*3);
  _logger->eventEnd(computeEvent);
#endif

  _needNewJacobian = false;
  _material->resetNeedNewJacobian();
} // integrateJacobian

// ----------------------------------------------------------------------
// Compute area of triangular cell.
double
pylith::feassemble::ElasticityExplicitTri3::_area(
			     const double_array& coordinatesCell) const
{ // __area
  assert(6 == coordinatesCell.size());

  const double x0 = coordinatesCell[0];
  const double y0 = coordinatesCell[1];

  const double x1 = coordinatesCell[2];
  const double y1 = coordinatesCell[3];

  const double x2 = coordinatesCell[4];
  const double y2 = coordinatesCell[5];

  const double area = 0.5*((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0));
  PetscLogFlops(8);

  return area;  
} // _area


// End of file 
