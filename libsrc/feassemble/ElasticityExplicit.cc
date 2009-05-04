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

#include "ElasticityExplicit.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature
#include "CellGeometry.hh" // USES CellGeometry

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

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ElasticityExplicit::ElasticityExplicit(void) :
  _dtm1(-1.0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ElasticityExplicit::~ElasticityExplicit(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::ElasticityExplicit::timeStep(const double dt)
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
pylith::feassemble::ElasticityExplicit::useSolnIncr(const bool flag)
{ // useSolnIncr
  if (flag)
    throw std::logic_error("Incremental solution not supported for "
			   "explicit time integration of elasticity "
			   "equation.");
} // useSolnIncr

// ----------------------------------------------------------------------
// Integrate constributions to residual term (r) for operator.
void
pylith::feassemble::ElasticityExplicit::integrateResidual(
			  const topology::Field<topology::Mesh>& residual,
			  const double t,
			  topology::SolutionFields* const fields)
{ // integrateResidual
  /// Member prototype for _elasticityResidualXD()
  typedef void (pylith::feassemble::ElasticityExplicit::*elasticityResidual_fn_type)
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
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const int tensorSize = _material->tensorSize();
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

  // Set variables dependent on dimension of cell
  totalStrain_fn_type calcTotalStrainFn;
  elasticityResidual_fn_type elasticityResidualFn;
  if (1 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityExplicit::_elasticityResidual1D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityExplicit::_elasticityResidual2D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityExplicit::_elasticityResidual3D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else
    assert(0);

  // Allocate vectors for cell values.
  double_array dispTCell(numBasis*spaceDim);
  double_array dispTmdtCell(numBasis*spaceDim);
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
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();
  assert(!dispTSection.isNull());
  topology::Mesh::RestrictVisitor dispTVisitor(*dispTSection,
					       numBasis*spaceDim, 
					       &dispTCell[0]);
  const ALE::Obj<RealSection>& dispTmdtSection = 
    fields->get("disp(t-dt)").section();
  assert(!dispTmdtSection.isNull());
  topology::Mesh::RestrictVisitor dispTmdtVisitor(*dispTmdtSection,
					       numBasis*spaceDim, 
					       &dispTmdtCell[0]);
  const ALE::Obj<RealSection>& residualSection = residual.section();
  topology::Mesh::UpdateAddVisitor residualVisitor(*residualSection,
						   &_cellVector[0]);

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
  const double gravityScale = 
    _normalizer->pressureScale() / (_normalizer->lengthScale() *
				    _normalizer->densityScale());

  // Get parameters used in integration.
  const double dt = _dt;
  const double dt2 = dt*dt;
  assert(dt > 0);

  _logger->eventEnd(setupEvent);

  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    _logger->eventBegin(geometryEvent);
    _quadrature->retrieveGeometry(*c_iter);
    _logger->eventEnd(geometryEvent);

    // Get state variables for cell.
    _logger->eventBegin(stateVarsEvent);
    _material->retrievePropsAndVars(*c_iter);
    _logger->eventEnd(stateVarsEvent);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    _logger->eventBegin(restrictEvent);
    dispTVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispTVisitor);
    dispTmdtVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispTmdtVisitor);
    _logger->eventEnd(restrictEvent);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute action for inertial terms
    _logger->eventBegin(computeEvent);
    const double_array& density = _material->calcDensity();
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad] / dt2;
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQuad*numBasis+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQuad*numBasis+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis*spaceDim+iDim] += 
	      valIJ * (+ 2.0 * dispTCell[jBasis*spaceDim+iDim]
		       - dispTmdtCell[jBasis*spaceDim+iDim]);
        } // for
      } // for
    } // for
    PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(6*spaceDim))));
    _logger->eventEnd(computeEvent);

    // Compute B(transpose) * sigma, first computing strains
    _logger->eventBegin(stressEvent);
    calcTotalStrainFn(&strainCell, basisDeriv, dispTCell, 
		      numBasis, numQuadPts);
    const double_array& stressCell = _material->calcStress(strainCell, true);
    _logger->eventEnd(stressEvent);

    _logger->eventBegin(computeEvent);
    CALL_MEMBER_FN(*this, elasticityResidualFn)(stressCell);
    _logger->eventEnd(computeEvent);

    // Assemble cell contribution into field
    _logger->eventBegin(updateEvent);
    residualVisitor.clear();
    sieveMesh->updateAdd(*c_iter, residualVisitor);
    _logger->eventEnd(updateEvent);
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ElasticityExplicit::integrateJacobian(
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

  // Allocate vectors for cell data.
  double_array dispTCell(numBasis*spaceDim);

  // Get cell information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const int materialId = _material->id();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<RealSection>& dispTSection = 
    fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  // Get sparse matrix
  const PetscMat jacobianMat = jacobian->matrix();
  assert(0 != jacobianMat);

  // Get parameters used in integration.
  const double dt = _dt;
  const double dt2 = dt*dt;
  assert(dt > 0);

  const ALE::Obj<SieveMesh::order_type>& globalOrder = 
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", dispTSection);
  assert(!globalOrder.isNull());
  // We would need to request unique points here if we had an interpolated mesh
  topology::Mesh::IndicesVisitor jacobianVisitor(*dispTSection, *globalOrder,
		  (int) pow(sieveMesh->getSieve()->getMaxConeSize(),
			    sieveMesh->depth())*spaceDim);

  _logger->eventEnd(setupEvent);

  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    _logger->eventBegin(geometryEvent);
    _quadrature->retrieveGeometry(*c_iter);
    _logger->eventEnd(geometryEvent);

    // Get state variables for cell.
    _logger->eventBegin(stateVarsEvent);
    _material->retrievePropsAndVars(*c_iter);
    _logger->eventEnd(stateVarsEvent);

    // Reset element matrix to zero
    _resetCellMatrix();

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Get material physical properties at quadrature points for this cell
    const double_array& density = _material->calcDensity();

    // Compute Jacobian for inertial terms
    _logger->eventBegin(computeEvent);
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad] / dt2;
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQ+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQ+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim) {
            const int iBlock = (iBasis*spaceDim + iDim) * (numBasis*spaceDim);
            const int jBlock = (jBasis*spaceDim + iDim);
            _cellMatrix[iBlock+jBlock] += valIJ;
          } // for
        } // for
      } // for
    } // for
    PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(1+spaceDim))));
    _logger->eventEnd(computeEvent);
    
    // Assemble cell contribution into PETSc matrix.
    _logger->eventBegin(updateEvent);
    jacobianVisitor.clear();
    PetscErrorCode err = updateOperator(jacobianMat, *sieveMesh->getSieve(),
					jacobianVisitor, *c_iter,
					&_cellMatrix[0], ADD_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");
    _logger->eventEnd(updateEvent);
  } // for

  _needNewJacobian = false;
  _material->resetNeedNewJacobian();
} // integrateJacobian


// End of file 
