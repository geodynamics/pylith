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
#include "Elasticity.hh" // USES Elasticity

#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
#include "pylith/utils/array.hh" // USES double_array

#include "petscmat.h" // USES PetscMat
#include "spatialdata/spatialdb/SpatialDB.hh"

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ElasticityExplicit::ElasticityExplicit(void) :
  _dtm1(-1.0),
  _material(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ElasticityExplicit::~ElasticityExplicit(void)
{ // destructor
  _material = 0; // Don't manage memory for material
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
// Set material.
void
pylith::feassemble::ElasticityExplicit::material(materials::ElasticMaterial* m)
{ // material
  _material = m;
  if (0 != _material)
    _material->timeStep(_dt);  
} // material

// ----------------------------------------------------------------------
// Determine whether we need to recompute the Jacobian.
bool
pylith::feassemble::ElasticityExplicit::needNewJacobian(void)
{ // needNewJacobian
  assert(0 != _material);
  if (!_needNewJacobian)
    _needNewJacobian = _material->needNewJacobian();
  return _needNewJacobian;
} // needNewJacobian

// ----------------------------------------------------------------------
// Integrate constributions to residual term (r) for operator.
void
pylith::feassemble::ElasticityExplicit::integrateResidual(
			      const ALE::Obj<real_section_type>& residual,
			      const double t,
			      topology::FieldsManager* const fields,
			      const ALE::Obj<Mesh>& mesh)
{ // integrateResidual
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(!residual.isNull());
  assert(0 != fields);
  assert(!mesh.isNull());

  PetscErrorCode err = 0;

  // Get cell information
  const ALE::Obj<ALE::Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", _material->id());
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator  cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const ALE::Obj<real_section_type>& dispTpdt = fields->getHistoryItem(0);
  const ALE::Obj<real_section_type>& dispT = fields->getHistoryItem(1);
  const ALE::Obj<real_section_type>& dispTmdt = fields->getHistoryItem(2);
  assert(!dispTpdt.isNull());
  assert(!dispT.isNull());
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
  const int cellDim = _quadrature->cellDim();

  /** :TODO:
   *
   * If cellDim and spaceDim are different, we need to transform
   * displacements into cellDim, compute action, and transform result
   * back into spaceDim. We get this information from the Jacobian and
   * inverse of the Jacobian.
   */
  if (cellDim != spaceDim)
    throw std::logic_error("Not implemented yet.");

  // Allocate vectors for cell values.
  _initCellVector();
  const int cellVecSize = numBasis*spaceDim;
  double_array dispTpdtCell(cellVecSize);
  double_array dispTCell(cellVecSize);
  double_array dispTmdtCell(cellVecSize);

  // Allocate vector for total strain
  int tensorSize = 0;
  if (1 == cellDim)
    tensorSize = 1;
  else if (2 == cellDim)
    tensorSize = 3;
  else if (3 == cellDim)
    tensorSize = 6;
  else {
    std::cerr << "Unknown case for cellDim '" << cellDim << "'." << std::endl;
    assert(0);
  } // else
  std::vector<double_array> totalStrain(numQuadPts);
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    totalStrain[iQuad].resize(tensorSize);
    totalStrain[iQuad] = 0.0;
  } // for

  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(mesh, coordinates, *c_iter);

    // Get state variables for cell.
    _material->getStateVarsCell(*c_iter, numQuadPts);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    mesh->restrict(dispTpdt, *c_iter, &dispTpdtCell[0], cellVecSize);
    mesh->restrict(dispT, *c_iter, &dispTCell[0], cellVecSize);
    mesh->restrict(dispTmdt, *c_iter, &dispTmdtCell[0], cellVecSize);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute action for inertial terms
    const std::vector<double_array>& density = 
      _material->calcDensity();
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad][0] / dt2;
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQuad*numBasis+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQuad*numBasis+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis*spaceDim+iDim] += 
	      valIJ * (- dispTpdtCell[jBasis*spaceDim+iDim] 
		       + 2.0 * dispTCell[jBasis*spaceDim+iDim]
		       - dispTmdtCell[jBasis*spaceDim+iDim]);
        } // for
      } // for
    } // for
    err = PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(6*spaceDim))));
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");

    // Compute action for elastic terms
    if (1 == cellDim) {
      // Compute stresses
      Elasticity::calcTotalStrain1D(&totalStrain, basisDeriv,
				    dispTCell, numBasis);
      const std::vector<double_array>& stress = 
	_material->calcStress(totalStrain);

      // Compute elastic action
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad];
	const double s11 = stress[iQuad][0];
	for (int iBasis=0; iBasis < numBasis; ++iBasis) {
	  const double N1 = wt*basisDeriv[iQuad*numBasis+iBasis  ];
	  _cellVector[iBasis*spaceDim  ] -= N1*s11;
	} // for
      } // for
      err = PetscLogFlops(numQuadPts*(1+numBasis*5));
      if (err)
	throw std::runtime_error("Logging PETSc flops failed.");

    } else if (2 == cellDim) {
      // Compute stresses
      Elasticity::calcTotalStrain2D(&totalStrain, basisDeriv,
				    dispTCell, numBasis);
      const std::vector<double_array>& stress = 
	_material->calcStress(totalStrain);
      
      // Compute elastic action
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad];
	const double s11 = stress[iQuad][0];
	const double s22 = stress[iQuad][1];
	const double s12 = stress[iQuad][2];
	for (int iBasis=0, iQ=iQuad*numBasis*cellDim;
	     iBasis < numBasis;
	     ++iBasis) {
	  const double N1 = wt*basisDeriv[iQ+iBasis*cellDim  ];
	  const double N2 = wt*basisDeriv[iQ+iBasis*cellDim+1];
	  _cellVector[iBasis*spaceDim  ] -= N1*s11 + N2*s12;
	  _cellVector[iBasis*spaceDim+1] -= N1*s12 + N2*s22;
	} // for
      } // for
      err = PetscLogFlops(numQuadPts*(1+numBasis*(8+2+9)));
      if (err)
	throw std::runtime_error("Logging PETSc flops failed.");
      
    } else if (3 == cellDim) {
      // Compute stresses
      Elasticity::calcTotalStrain3D(&totalStrain, basisDeriv, 
				    dispTCell, numBasis);
      const std::vector<double_array>& stress = 
	_material->calcStress(totalStrain);

      // Compute elastic action
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad];
	const double s11 = stress[iQuad][0];
	const double s22 = stress[iQuad][1];
	const double s33 = stress[iQuad][2];
	const double s12 = stress[iQuad][3];
	const double s23 = stress[iQuad][4];
	const double s13 = stress[iQuad][5];

	for (int iBasis=0, iQ=iQuad*numBasis*cellDim;
	     iBasis < numBasis;
	     ++iBasis) {
	  const int iBlock = iBasis*spaceDim;
	  const double N1 = wt*basisDeriv[iQ+iBasis*cellDim+0];
	  const double N2 = wt*basisDeriv[iQ+iBasis*cellDim+1];
	  const double N3 = wt*basisDeriv[iQ+iBasis*cellDim+2];
	  _cellVector[iBlock  ] -= N1*s11 + N2*s12 + N3*s13;
	  _cellVector[iBlock+1] -= N1*s12 + N2*s22 + N3*s23;
	  _cellVector[iBlock+2] -= N1*s13 + N2*s23 + N3*s33;
	} // for
      } // for
      err = PetscLogFlops(numQuadPts*(1+numBasis*(3+12)));
      if (err)
	throw std::runtime_error("Logging PETSc flops failed.");
    } else {
      std::cerr << "Unknown case for cellDim '" << cellDim << "'."
		<< std::endl;
      assert(0);
    } // if/else
    // Assemble cell contribution into field
    mesh->updateAdd(residual, *c_iter, _cellVector);
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
  const ALE::Obj<real_section_type>& dispT = fields->getHistoryItem(1);
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

  // Allocate vector for cell values (if necessary)
  _initCellMatrix();

  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(mesh, coordinates, *c_iter);

    // Get state variables for cell.
    _material->getStateVarsCell(*c_iter, numQuadPts);

    // Reset element vector to zero
    _resetCellMatrix();

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Get material physical properties at quadrature points for this cell
    const std::vector<double_array>& density = _material->calcDensity();

    // Compute Jacobian for inertial terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad][0] / dt2;
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
    PetscErrorCode err = 
      PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(1+spaceDim))));
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");
    
    // Assemble cell contribution into PETSc Matrix
    const ALE::Obj<Mesh::order_type>& globalOrder = 
      mesh->getFactory()->getGlobalOrder(mesh, "default", dispT);
    assert(!globalOrder.isNull());

    err = updateOperator(*jacobian, mesh, dispT, globalOrder,
			 *c_iter, _cellMatrix, ADD_VALUES);
    if (err)
      throw std::runtime_error("Update to PETSc Mat failed.");
  } // for

  _needNewJacobian = false;
  _material->resetNeedNewJacobian();
} // integrateJacobian

// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::ElasticityExplicit::updateState(
				   const double t,
				   const ALE::Obj<real_section_type>& disp,
				   const ALE::Obj<Mesh>& mesh)
{ // updateState
  assert(0 != _material);
  assert(!disp.isNull());

  // No need to update state if using elastic behavior
  std::cout << "USES_UPDATE_STATE: " << _material->usesUpdateState() << std::endl;
  if (!_material->usesUpdateState())
    return;

  // Get cell information
  const ALE::Obj<ALE::Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", _material->id());
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();

  const int cellVecSize = numBasis*spaceDim;
  double_array dispCell(cellVecSize);

  // Allocate vector for total strain
  int tensorSize = 0;
  if (1 == cellDim)
    tensorSize = 1;
  else if (2 == cellDim)
    tensorSize = 3;
  else if (3 == cellDim)
    tensorSize = 6;
  else {
    std::cerr << "Unknown case for cellDim '" << cellDim << "'." << std::endl;
    assert(0);
  } // else
  std::vector<double_array> totalStrain(numQuadPts);
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    totalStrain[iQuad].resize(tensorSize);
    totalStrain[iQuad] = 0.0;
  } // for

  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(mesh, coordinates, *c_iter);

    // Restrict input fields to cell
    mesh->restrict(disp, *c_iter, &dispCell[0], cellVecSize);

    // Get cell geometry information that depends on cell
    const double_array& basisDeriv = _quadrature->basisDeriv();
  
    // Compute action for elastic terms
    if (1 == cellDim)
      Elasticity::calcTotalStrain1D(&totalStrain, basisDeriv,
				    dispCell, numBasis);
    else if (2 == cellDim)
      Elasticity::calcTotalStrain2D(&totalStrain, basisDeriv,
				    dispCell, numBasis);
    else if (3 == cellDim)
      Elasticity::calcTotalStrain3D(&totalStrain, basisDeriv, 
				    dispCell, numBasis);
    else {
      std::cerr << "Unknown case for cellDim '" << cellDim << "'."
		<< std::endl;
      assert(0);
    } // else
    _material->updateState(totalStrain, *c_iter);
  } // for
} // updateState


// End of file 
