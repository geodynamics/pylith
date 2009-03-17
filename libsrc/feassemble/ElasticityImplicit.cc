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

#include "ElasticityImplicit.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature
#include "CellGeometry.hh" // USES CellGeometry

#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/array.hh" // USES double_array
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN
#include "pylith/utils/lapack.h" // USES LAPACKdgesvd

#include "petscmat.h" // USES PetscMat
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimendional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

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
} // destructor
  
// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::ElasticityImplicit::timeStep(const double dt)
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
double
pylith::feassemble::ElasticityImplicit::stableTimeStep(void) const
{ // stableTimeStep
  assert(0 != _material);
  return _material->stableTimeStepImplicit();
} // stableTimeStep

// ----------------------------------------------------------------------
// Set flag for setting constraints for total field solution or
// incremental field solution.
void
pylith::feassemble::ElasticityImplicit::useSolnIncr(const bool flag)
{ // useSolnIncr
  assert(0 != _material);
  _useSolnIncr = flag;
  _material->useElasticBehavior(!_useSolnIncr);
} // useSolnIncr

// ----------------------------------------------------------------------
// Integrate constributions to residual term (r) for operator.
void
pylith::feassemble::ElasticityImplicit::integrateResidual(
			  const topology::Field<topology::Mesh>& residual,
			  const double t,
			  topology::SolutionFields* const fields)
{ // integrateResidual
  /// Member prototype for _elasticityResidualXD()
  typedef void (pylith::feassemble::ElasticityImplicit::*elasticityResidual_fn_type)
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
  if (cellDim != spaceDim)
    throw std::logic_error("Integration for cells with spatial dimensions "
			   "different than the spatial dimension of the "
			   "domain not implemented yet.");

  // Set variables dependent on dimension of cell
  totalStrain_fn_type calcTotalStrainFn;
  elasticityResidual_fn_type elasticityResidualFn;
  if (1 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityResidual1D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityResidual2D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityResidual3D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else
    assert(0);

  // Allocate vectors for cell values.
  double_array dispTBctpdtCell(numBasis*spaceDim);
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
  const ALE::Obj<RealSection>& dispTBctpdtSection = 
    fields->get("disp(t), bc(t+dt)").section();
  assert(!dispTBctpdtSection.isNull());
  topology::Mesh::RestrictVisitor dispTBctpdtVisitor(*dispTBctpdtSection,
						     numBasis*spaceDim, 
						     &dispTBctpdtCell[0]);
  const ALE::Obj<RealSection>& residualSection = residual.section();
  topology::Mesh::UpdateAddVisitor residualVisitor(*residualSection,
						   &_cellVector[0]);

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
  const double gravityScale = 
    _normalizer->pressureScale() / (_normalizer->lengthScale() *
				    _normalizer->densityScale());

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
    dispTBctpdtVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispTBctpdtVisitor);
    _logger->eventBegin(restrictEvent);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();
    const double_array& quadPtsNondim = _quadrature->quadPts();

    // Compute body force vector if gravity is being used.
    if (0 != _gravityField) {
      _logger->eventBegin(computeEvent);
      const spatialdata::geocoords::CoordSys* cs = fields->mesh().coordsys();
      assert(0 != cs);
      
      // Get density at quadrature points for this cell
      const double_array& density = _material->calcDensity();

      quadPtsGlobal = quadPtsNondim;
      _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
				  lengthScale);

      // Compute action for element body forces
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const int err = _gravityField->query(&gravVec[0], gravVec.size(),
					     &quadPtsGlobal[0], spaceDim, cs);
	if (err)
	  throw std::runtime_error("Unable to get gravity vector for point.");
	_normalizer->nondimensionalize(&gravVec[0], gravVec.size(), 
				       gravityScale);
	const double wt = quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad];
	for (int iBasis=0, iQ=iQuad*numBasis;
	     iBasis < numBasis; ++iBasis) {
	  const double valI = wt*basis[iQ+iBasis];
	  for (int iDim=0; iDim < spaceDim; ++iDim) {
	    _cellVector[iBasis*spaceDim+iDim] += valI*gravVec[iDim];
	  } // for
	} // for
      } // for
      PetscLogFlops(numQuadPts*(2+numBasis*(1+2*spaceDim)));
      _logger->eventEnd(computeEvent);      
    } // if

    // Compute B(transpose) * sigma, first computing strains
    _logger->eventBegin(stressEvent);
    calcTotalStrainFn(&strainCell, basisDeriv, dispTBctpdtCell, 
		      numBasis, numQuadPts);
    const double_array& stressCell = _material->calcStress(strainCell, true);
    _logger->eventEnd(stressEvent);

    _logger->eventBegin(computeEvent);
    CALL_MEMBER_FN(*this, elasticityResidualFn)(stressCell);
    _logger->eventEnd(computeEvent);

#if 0
    std::cout << "Updating residual for cell " << *c_iter << std::endl;
    for(int i = 0; i < _quadrature->spaceDim() * _quadrature->numBasis(); ++i) {
      std::cout << "  v["<<i<<"]: " << _cellVector[i] << std::endl;
    }
#endif
    // Assemble cell contribution into field
    _logger->eventBegin(updateEvent);
    residualVisitor.clear();
    sieveMesh->updateAdd(*c_iter, residualVisitor);
    _logger->eventEnd(updateEvent);
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Compute stiffness matrix.
void
pylith::feassemble::ElasticityImplicit::integrateJacobian(
					topology::Jacobian* jacobian,
					const double t,
					topology::SolutionFields* fields)
{ // integrateJacobian
  /// Member prototype for _elasticityJacobianXD()
  typedef void (pylith::feassemble::ElasticityImplicit::*elasticityJacobian_fn_type)
    (const double_array&);

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

  // Set variables dependent on dimension of cell
  totalStrain_fn_type calcTotalStrainFn;
  elasticityJacobian_fn_type elasticityJacobianFn;
  if (1 == cellDim) {
    elasticityJacobianFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityJacobian1D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    elasticityJacobianFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityJacobian2D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    elasticityJacobianFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityJacobian3D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else
    assert(0);

  // Allocate vector for total strain
  double_array dispTBctpdtCell(numBasis*spaceDim);
  double_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;

  // Get cell information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const int materialId = _material->id();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<RealSection>& dispTBctpdtSection = 
    fields->get("disp(t), bc(t+dt)").section();
  assert(!dispTBctpdtSection.isNull());
  topology::Mesh::RestrictVisitor dispTBctpdtVisitor(*dispTBctpdtSection,
						     numBasis*spaceDim, 
						     &dispTBctpdtCell[0]);

  // Get sparse matrix
  const PetscMat jacobianMat = jacobian->matrix();
  assert(0 != jacobianMat);

  // Get parameters used in integration.
  const double dt = _dt;
  assert(dt > 0);

  const ALE::Obj<SieveMesh::order_type>& globalOrder = 
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default",
					    dispTBctpdtSection);
  assert(!globalOrder.isNull());
  // We would need to request unique points here if we had an interpolated mesh
  topology::Mesh::IndicesVisitor jacobianVisitor(*dispTBctpdtSection,
						 *globalOrder,
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

    // Restrict input fields to cell
    _logger->eventBegin(restrictEvent);
    dispTBctpdtVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispTBctpdtVisitor);
    _logger->eventBegin(restrictEvent);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    _logger->eventBegin(computeEvent);
    // Compute strains
    calcTotalStrainFn(&strainCell, basisDeriv, dispTBctpdtCell, 
		      numBasis, numQuadPts);
      
    // Get "elasticity" matrix at quadrature points for this cell
    const double_array& elasticConsts = 
      _material->calcDerivElastic(strainCell);

    CALL_MEMBER_FN(*this, elasticityJacobianFn)(elasticConsts);
    _logger->eventEnd(computeEvent);

    if (_quadrature->checkConditioning()) {
      int n = numBasis*spaceDim;
      int lwork = 5*n;
      int idummy = 0;
      int lierr = 0;
      double *elemMat = new double[n*n];
      double *svalues = new double[n];
      double *work    = new double[lwork];
      double minSV = 0;
      double maxSV = 0;
      double sdummy = 0;

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
      std::cout << "Element " << *c_iter << std::endl;
      for(int i = 0; i < n; ++i)
	std::cout << "    sV["<<i<<"] = " << svalues[i] << std::endl;
      std::cout << "  kappa(elemMat) = " << maxSV/minSV << std::endl;
      delete [] elemMat;
      delete [] svalues;
      delete [] work;
    } // if

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
