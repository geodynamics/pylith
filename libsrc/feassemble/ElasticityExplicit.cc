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
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
#include "pylith/utils/array.hh" //   USES double_array
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN

#include "petscmat.h" // USES PetscMat
#include "spatialdata/spatialdb/SpatialDB.hh"

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error

#define FASTER

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
			      const ALE::Obj<real_section_type>& residual,
			      const double t,
			      topology::FieldsManager* const fields,
			      const ALE::Obj<Mesh>& mesh)
{ // integrateResidual
  /// Member prototype for _elasticityResidualXD()
  typedef void (pylith::feassemble::ElasticityExplicit::*elasticityResidual_fn_type)
    (const double_array&);

  assert(0 != _quadrature);
  assert(0 != _material);
  assert(!residual.isNull());
  assert(0 != fields);
  assert(!mesh.isNull());

  PetscErrorCode err = 0;

  // Set variables dependent on dimension of cell
  const int cellDim = _quadrature->cellDim();
  int tensorSize = 0;
  totalStrain_fn_type calcTotalStrainFn;
  elasticityResidual_fn_type elasticityResidualFn;
  if (1 == cellDim) {
    tensorSize = 1;
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityExplicit::_elasticityResidual1D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    tensorSize = 3;
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityExplicit::_elasticityResidual2D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    tensorSize = 6;
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityExplicit::_elasticityResidual3D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else
    assert(0);

  // Get cell information
  const int materialId = _material->id();
  const ALE::Obj<ALE::Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", materialId);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator  cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const ALE::Obj<real_section_type>& dispT = fields->getFieldByHistory(1);
  assert(!dispT.isNull());
  const ALE::Obj<real_section_type>& dispTmdt = fields->getFieldByHistory(2);
  assert(!dispTmdt.isNull());

  // Get parameters used in integration.
  const double dt = _dt;
  const double dt2 = dt*dt;
  assert(dt > 0);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  /** :TODO:
   *
   * If cellDim and spaceDim are different, we need to transform
   * displacements into cellDim, compute action, and transform result
   * back into spaceDim. We get this information from the Jacobian and
   * inverse of the Jacobian.
   */
  if (cellDim != spaceDim)
    throw std::logic_error("Not implemented yet.");

  // Precompute the geometric and function space information
  _quadrature->precomputeGeometry(mesh, coordinates, cells);

  // Allocate vectors for cell values.
  _initCellVector();
  const int cellVecSize = numBasis*spaceDim;
  double_array dispTCell(cellVecSize);
  double_array dispTmdtCell(cellVecSize);

  // Allocate vector for total strain
  double_array totalStrain(numQuadPts*tensorSize);
  totalStrain = 0.0;

#ifdef FASTER
  fields->createCustomAtlas("material-id", materialId);
  const int dispTAtlasTag = 
    fields->getFieldAtlasTagByHistory(1, materialId);
  const int dispTmdtAtlasTag = 
    fields->getFieldAtlasTagByHistory(2, materialId);
  const int residualAtlasTag = 
    fields->getFieldAtlasTag("residual", materialId);
#endif

  int c_index = 0;
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++c_index) {
    // Compute geometry information for current cell
    _quadrature->retrieveGeometry(mesh, coordinates, *c_iter, c_index);

    // Get state variables for cell.
    _material->getStateVarsCell(*c_iter, numQuadPts);

    // Reset element vector to zero
    _resetCellVector();

#ifdef FASTER
    mesh->restrict(dispT, dispTAtlasTag, c_index, &dispTCell[0], 
		   cellVecSize);
    mesh->restrict(dispTmdt, dispTmdtAtlasTag, c_index, &dispTmdtCell[0], 
		   cellVecSize);
#else
    mesh->restrict(dispT, *c_iter, &dispTCell[0], cellVecSize);
    mesh->restrict(dispTmdt, *c_iter, &dispTmdtCell[0], cellVecSize);
#endif

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute action for inertial terms
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
    PetscLogFlopsNoCheck(numQuadPts*(3+numBasis*(1+numBasis*(6*spaceDim))));

    // Compute B(transpose) * sigma, first computing strains
    calcTotalStrainFn(&totalStrain, basisDeriv, dispTCell, 
		      numBasis, numQuadPts);
    const double_array& stress = _material->calcStress(totalStrain);
    CALL_MEMBER_FN(*this, elasticityResidualFn)(stress);

    // Assemble cell contribution into field
#ifdef FASTER
    mesh->updateAdd(residual, residualAtlasTag, c_index, _cellVector);
#else
    mesh->updateAdd(residual, *c_iter, _cellVector);
#endif
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ElasticityExplicit::integrateJacobian(
					PetscMat* jacobian,
					const double t,
					topology::FieldsManager* fields,
					const ALE::Obj<Mesh>& mesh)
{ // integrateJacobian
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(0 != jacobian);
  assert(0 != fields);
  assert(!mesh.isNull());

  // Get cell information
  const ALE::Obj<ALE::Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", _material->id());
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator  cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const ALE::Obj<real_section_type>& dispT = fields->getFieldByHistory(1);
  assert(!dispT.isNull());

  // Get parameters used in integration.
  const double dt = _dt;
  const double dt2 = dt*dt;
  assert(dt > 0);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();

  // Precompute the geometric and function space information
  _quadrature->precomputeGeometry(mesh, coordinates, cells);

  // Allocate vector for cell values (if necessary)
  _initCellMatrix();

  int c_index = 0;
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++c_index) {
    // Compute geometry information for current cell
    _quadrature->retrieveGeometry(mesh, coordinates, *c_iter, c_index);

    // Get state variables for cell.
    _material->getStateVarsCell(*c_iter, numQuadPts);

    // Reset element vector to zero
    _resetCellMatrix();

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Get material physical properties at quadrature points for this cell
    const double_array& density = _material->calcDensity();

    // Compute Jacobian for inertial terms
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
    PetscLogFlopsNoCheck(numQuadPts*(3+numBasis*(1+numBasis*(1+spaceDim))));
    
    // Assemble cell contribution into PETSc Matrix
    const ALE::Obj<Mesh::order_type>& globalOrder = 
      mesh->getFactory()->getGlobalOrder(mesh, "default", dispT);
    assert(!globalOrder.isNull());

    PetscErrorCode err = updateOperator(*jacobian, mesh, dispT, globalOrder,
					*c_iter, _cellMatrix, ADD_VALUES);
    if (err)
      throw std::runtime_error("Update to PETSc Mat failed.");
  } // for

  _needNewJacobian = false;
  _material->resetNeedNewJacobian();
} // integrateJacobian


// End of file 
