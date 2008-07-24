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

#include "IntegratorElasticity.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature
#include "CellGeometry.hh" // USES CellGeometry

#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
#include "pylith/utils/array.hh" // USES double_array

#include <string.h> // USES memcpy()
#include <strings.h> // USES strcasecmp()
#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorElasticity::IntegratorElasticity(void) :
  _material(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorElasticity::~IntegratorElasticity(void)
{ // destructor
  _material = 0; // Don't manage memory for material
} // destructor
  
// ----------------------------------------------------------------------
// Set material.
void
pylith::feassemble::IntegratorElasticity::material(materials::ElasticMaterial* m)
{ // material
  _material = m;
  if (0 != _material)
    _material->timeStep(_dt);  
} // material

// ----------------------------------------------------------------------
// Determine whether we need to recompute the Jacobian.
bool
pylith::feassemble::IntegratorElasticity::needNewJacobian(void)
{ // needNewJacobian
  assert(0 != _material);
  if (!_needNewJacobian)
    _needNewJacobian = _material->needNewJacobian();
  return _needNewJacobian;
} // needNewJacobian

// ----------------------------------------------------------------------
// Set flag for setting constraints for total field solution or
void
pylith::feassemble::IntegratorElasticity::useSolnIncr(const bool flag)
{ // useSolnIncr
  Integrator::useSolnIncr(flag);
  
  assert(0 != _material);
  _material->useElasticBehavior(!flag);
} // useSolnIncr

// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::IntegratorElasticity::updateState(
				   const double t,
				   topology::FieldsManager* const fields,
				   const ALE::Obj<Mesh>& mesh)
{ // updateState
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(0 != fields);

  // No need to update state if using elastic behavior
  if (!_material->usesUpdateProperties())
    return;

  // Set variables dependent on dimension of cell
  const int cellDim = _quadrature->cellDim();
  int tensorSize = 0;
  totalStrain_fn_type calcTotalStrainFn;
  if (1 == cellDim) {
    tensorSize = 1;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    tensorSize = 3;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    tensorSize = 6;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else
    assert(0);

  // Get cell information
  const int materialId = _material->id();
  const ALE::Obj<Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", materialId);
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

  const int cellVecSize = numBasis*spaceDim;
  double_array dispCell(cellVecSize);

  // Allocate vector for total strain
  double_array totalStrain(numQuadPts*tensorSize);
  totalStrain = 0.0;

  const ALE::Obj<real_section_type>& disp = fields->getSolution();
  /// const int dispAtlasTag = fields->getSolutionAtlasTag(materialId);
  
  // Loop over cells
  int c_index = 0;
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++c_index) {
    // Compute geometry information for current cell
    _quadrature->retrieveGeometry(mesh, coordinates, *c_iter, c_index);

    // Restrict input fields to cell
    mesh->restrictClosure(disp, *c_iter, &dispCell[0], cellVecSize);

    // Get cell geometry information that depends on cell
    const double_array& basisDeriv = _quadrature->basisDeriv();
  
    // Compute strains
    calcTotalStrainFn(&totalStrain, basisDeriv, dispCell, numBasis, numQuadPts);

    // Update material state
    _material->updateProperties(totalStrain, *c_iter);
  } // for
} // updateState

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::feassemble::IntegratorElasticity::verifyConfiguration(
						 const ALE::Obj<Mesh>& mesh) const
{ // verifyConfiguration
  assert(0 != _quadrature);
  assert(0 != _material);

  const int dimension = mesh->getDimension();

  // check compatibility of mesh and material
  if (_material->dimension() != dimension) {
    std::ostringstream msg;
    msg << "Material '" << _material->label()
	<< "' is incompatible with mesh.\n"
	<< "Dimension of mesh: " << dimension
	<< ", dimension of material: " << _material->dimension()
	<< ".";
    throw std::runtime_error(msg.str());
  } // if

  // check compatibility of mesh and quadrature scheme
  if (_quadrature->cellDim() != dimension) {
    std::ostringstream msg;
    msg << "Quadrature is incompatible with cells for material '"
	<< _material->label() << "'.\n"
	<< "Dimension of mesh: " << dimension
	<< ", dimension of quadrature: " << _quadrature->cellDim()
	<< ".";
    throw std::runtime_error(msg.str());
  } // if
  const int numCorners = _quadrature->refGeometry().numCorners();
  const ALE::Obj<Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", _material->id());
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    const int cellNumCorners = mesh->getNumCellCorners(*c_iter);
    if (numCorners != cellNumCorners) {
      std::ostringstream msg;
      msg << "Quadrature is incompatible with cell in material '"
	  << _material->label() << "'.\n"
	  << "Cell " << *c_iter << " has " << cellNumCorners
	  << " vertices but quadrature reference cell has "
	  << numCorners << " vertices.";
      throw std::runtime_error(msg.str());
    } // if
  } // for
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get cell field associated with integrator.
const ALE::Obj<pylith::real_section_type>&
pylith::feassemble::IntegratorElasticity::cellField(
				 VectorFieldEnum* fieldType,
				 const char* name,
				 const ALE::Obj<Mesh>& mesh,
				 topology::FieldsManager* const fields)
{ // cellField
  assert(0 != _material);

  // We assume the material stores the total-strain field if
  // usesUpdateProperties() is TRUE.

  const int numQuadPts = _quadrature->numQuadPts();

  if (!_material->usesUpdateProperties() &&
      (0 == strcasecmp(name, "total_strain") ||
       0 == strcasecmp(name, "stress")) ) {
    assert(0 != fields);
    _calcStrainStressField(&_bufferCellOther, name, mesh, fields);
    return _bufferCellOther;

  } else if (0 == strcasecmp(name, "stress")) {
    _material->propertyField(&_bufferCellOther,
			     "total_strain", mesh, numQuadPts);
    _calcStressFromStrain(&_bufferCellOther, mesh);
    return _bufferCellOther;

  } else {
    const VectorFieldEnum fieldType = _material->propertyFieldType(name);
    switch (fieldType)
      { // switch
      case SCALAR_FIELD :
	_material->propertyField(&_bufferCellScalar, name, mesh, numQuadPts);
	return _bufferCellScalar;
	break;
      case VECTOR_FIELD :
	_material->propertyField(&_bufferCellVector, name, mesh, numQuadPts);
	return _bufferCellVector;
	break;
      case TENSOR_FIELD :
	_material->propertyField(&_bufferCellTensor, name, mesh, numQuadPts);
	return _bufferCellTensor;
	break;
      case OTHER_FIELD :
	_material->propertyField(&_bufferCellOther, name, mesh, numQuadPts);
	return _bufferCellOther;
	break;
      default:
	assert(0);
      } // switch
  } // else

  // Return scalar section to satisfy member function definition.
  return _bufferCellScalar;
} // cellField

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticity::_calcStrainStressField(
				 ALE::Obj<real_section_type>* field,
				 const char* name,
				 const ALE::Obj<Mesh>& mesh,
				 topology::FieldsManager* const fields)
{ // _calcStrainStressField
  assert(0 != _quadrature);
  assert(0 != _material);

  const bool calcStress = (0 == strcasecmp(name, "stress")) ? true : false;
    
  const int cellDim = _quadrature->cellDim();
  int tensorSize = 0;
  totalStrain_fn_type calcTotalStrainFn;
  if (1 == cellDim) {
    tensorSize = 1;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    tensorSize = 3;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    tensorSize = 6;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else
    assert(0);
  
  // Get cell information
  const int materialId = _material->id();
  const ALE::Obj<Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", materialId);
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
  
  const int cellVecSize = numBasis*spaceDim;
  double_array dispCell(cellVecSize);
  
  // Allocate array for total strain
  const int totalFiberDim = numQuadPts * tensorSize;
  double_array totalStrain(totalFiberDim);
  totalStrain = 0.0;
  
  // Allocate buffer for property field.
  if (field->isNull() || 
      totalFiberDim != (*field)->getFiberDimension(*cells->begin())) {
    *field = new real_section_type(mesh->comm(), mesh->debug());
    (*field)->setChart(real_section_type::chart_type(*std::min_element(cells->begin(), cells->end()),
                                                     *std::max_element(cells->begin(), cells->end())+1));
    (*field)->setFiberDimension(cells, totalFiberDim);
    mesh->allocate(*field);
  } // if
  
  const ALE::Obj<real_section_type>& disp = fields->getSolution();
  /// const int dispAtlasTag = fields->getSolutionAtlasTag(materialId);
  
  // Loop over cells
  int c_index = 0;
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++c_index) {
    // Compute geometry information for current cell
    _quadrature->retrieveGeometry(mesh, coordinates, *c_iter, c_index);
    
    // Restrict input fields to cell
    mesh->restrictClosure(disp, *c_iter, &dispCell[0], cellVecSize);
    
    // Get cell geometry information that depends on cell
    const double_array& basisDeriv = _quadrature->basisDeriv();
    
    // Compute strains
    calcTotalStrainFn(&totalStrain, basisDeriv, dispCell, numBasis, 
		      numQuadPts);
    
    if (!calcStress) {
      (*field)->updatePoint(*c_iter, &totalStrain[0]);
    } else {
      _material->getPropertiesCell(*c_iter, numQuadPts);
      const double_array& stress = _material->calcStress(totalStrain);
      (*field)->updatePoint(*c_iter, &stress[0]);	
    } // else
  } // for
} // _calcStrainStressField

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticity::_calcStressFromStrain(
				 ALE::Obj<real_section_type>* field,
				 const ALE::Obj<Mesh>& mesh)
{ // _calcStressFromStrain
  assert(0 != _quadrature);
  assert(0 != _material);

  const int cellDim = _quadrature->cellDim();
  int tensorSize = 0;
  if (1 == cellDim) {
    tensorSize = 1;
  } else if (2 == cellDim) {
    tensorSize = 3;
  } else if (3 == cellDim) {
    tensorSize = 6;
  } else
    assert(0);
  
  // Get cell information
  const int materialId = _material->id();
  const ALE::Obj<Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", materialId);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  
  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  
  // Allocate vector for total strain
  double_array totalStrain(numQuadPts*tensorSize);
  totalStrain = 0.0;
  
  // Allocate buffer for tensor field.
  if (field->isNull()) {
    const int fiberDim = numQuadPts * tensorSize;
    *field = new real_section_type(mesh->comm(), mesh->debug());
    (*field)->setChart(real_section_type::chart_type(*std::min_element(cells->begin(), cells->end()),
                                                     *std::max_element(cells->begin(), cells->end())+1));
    (*field)->setFiberDimension(cells, fiberDim);
    mesh->allocate(*field);
  } // if
  
  // Loop over cells
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    const real_section_type::value_type* strainVals = 
      (*field)->restrictPoint(*c_iter);
    memcpy(&totalStrain[0], strainVals, sizeof(double)*totalStrain.size());

    _material->getPropertiesCell(*c_iter, numQuadPts);
    const double_array& stress = _material->calcStress(totalStrain);
    (*field)->updatePoint(*c_iter, &stress[0]);	
  } // for
} // _calcStressFromStrain

// ----------------------------------------------------------------------
// Integrate elasticity term in residual for 1-D cells.
void
pylith::feassemble::IntegratorElasticity::_elasticityResidual1D(
				     const double_array& stress)
{ // _elasticityResidual1D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const double_array& quadWts = _quadrature->quadWts();
  const double_array& jacobianDet = _quadrature->jacobianDet();
  const double_array& basisDeriv = _quadrature->basisDeriv();

  assert(1 == cellDim);
  assert(quadWts.size() == numQuadPts);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    const double s11 = stress[iQuad];
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const double N1 = wt*basisDeriv[iQuad*numBasis+iBasis  ];
      _cellVector[iBasis*spaceDim  ] -= N1*s11;
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*5));
} // _elasticityResidual1D

// ----------------------------------------------------------------------
// Integrate elasticity term in residual for 2-D cells.
void
pylith::feassemble::IntegratorElasticity::_elasticityResidual2D(
				     const double_array& stress)
{ // _elasticityResidual2D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const double_array& quadWts = _quadrature->quadWts();
  const double_array& jacobianDet = _quadrature->jacobianDet();
  const double_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(2 == cellDim);
  assert(quadWts.size() == numQuadPts);
  const int stressSize = 3;

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    const double s11 = stress[iQuad*stressSize+0];
    const double s22 = stress[iQuad*stressSize+1];
    const double s12 = stress[iQuad*stressSize+2];
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim;
	 iBasis < numBasis;
	 ++iBasis) {
      const double N1 = wt*basisDeriv[iQ+iBasis*spaceDim  ];
      const double N2 = wt*basisDeriv[iQ+iBasis*spaceDim+1];
      _cellVector[iBasis*spaceDim  ] -= N1*s11 + N2*s12;
      _cellVector[iBasis*spaceDim+1] -= N1*s12 + N2*s22;
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(8+2+9)));
} // _elasticityResidual2D

// ----------------------------------------------------------------------
// Integrate elasticity term in residual for 3-D cells.
void
pylith::feassemble::IntegratorElasticity::_elasticityResidual3D(
				     const double_array& stress)
{ // _elasticityResidual3D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const double_array& quadWts = _quadrature->quadWts();
  const double_array& jacobianDet = _quadrature->jacobianDet();
  const double_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(3 == cellDim);
  assert(quadWts.size() == numQuadPts);
  const int stressSize = 6;
  
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    const double s11 = stress[iQuad*stressSize+0];
    const double s22 = stress[iQuad*stressSize+1];
    const double s33 = stress[iQuad*stressSize+2];
    const double s12 = stress[iQuad*stressSize+3];
    const double s23 = stress[iQuad*stressSize+4];
    const double s13 = stress[iQuad*stressSize+5];
    
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim;
	 iBasis < numBasis;
	 ++iBasis) {
      const int iBlock = iBasis*spaceDim;
      const double N1 = wt*basisDeriv[iQ+iBasis*spaceDim+0];
      const double N2 = wt*basisDeriv[iQ+iBasis*spaceDim+1];
      const double N3 = wt*basisDeriv[iQ+iBasis*spaceDim+2];

      _cellVector[iBlock  ] -= N1*s11 + N2*s12 + N3*s13;
      _cellVector[iBlock+1] -= N1*s12 + N2*s22 + N3*s23;
      _cellVector[iBlock+2] -= N1*s13 + N2*s23 + N3*s33;
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(3+12)));
} // _elasticityResidual3D

// ----------------------------------------------------------------------
// Integrate elasticity term in Jacobian for 1-D cells.
void
pylith::feassemble::IntegratorElasticity::_elasticityJacobian1D(
			       const double_array& elasticConsts)
{ // _elasticityJacobian1D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const double_array& quadWts = _quadrature->quadWts();
  const double_array& jacobianDet = _quadrature->jacobianDet();
  const double_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(1 == cellDim);
  assert(quadWts.size() == numQuadPts);
  
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    const double C1111 = elasticConsts[iQuad];
    for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
      const double valI = wt*basisDeriv[iQ+iBasis]*C1111;
      for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	const double valIJ = valI * basisDeriv[iQ+jBasis];
	const int iBlock = iBasis*spaceDim * (numBasis*spaceDim);
	const int jBlock = jBasis*spaceDim;
	_cellMatrix[iBlock+jBlock] += valIJ;
      } // for
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(2+numBasis*3)));
} // _elasticityJacobian1D

// ----------------------------------------------------------------------
// Integrate elasticity term in Jacobian for 2-D cells.
void
pylith::feassemble::IntegratorElasticity::_elasticityJacobian2D(
			       const double_array& elasticConsts)
{ // _elasticityJacobian2D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const double_array& quadWts = _quadrature->quadWts();
  const double_array& jacobianDet = _quadrature->jacobianDet();
  const double_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(2 == cellDim);
  assert(quadWts.size() == numQuadPts);
  const int numConsts = 6;

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    // tau_ij = C_ijkl * e_kl
    //        = C_ijlk * 0.5 (u_k,l + u_l,k)
    //        = 0.5 * C_ijkl * (u_k,l + u_l,k)
    // divide C_ijkl by 2 if k != l
    const double C1111 = elasticConsts[iQuad*numConsts+0];
    const double C1122 = elasticConsts[iQuad*numConsts+1];
    const double C1112 = elasticConsts[iQuad*numConsts+2]/2.0;
    const double C2222 = elasticConsts[iQuad*numConsts+3];
    const double C2212 = elasticConsts[iQuad*numConsts+4]/2.0;
    const double C1212 = elasticConsts[iQuad*numConsts+5]/2.0;
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim;
	 iBasis < numBasis;
	 ++iBasis) {
      const double Nip = wt*basisDeriv[iQ+iBasis*spaceDim  ];
      const double Niq = wt*basisDeriv[iQ+iBasis*spaceDim+1];
      const int iBlock = (iBasis*spaceDim  ) * (numBasis*spaceDim);
      const int iBlock1 = (iBasis*spaceDim+1) * (numBasis*spaceDim);
      for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	const double Njp = basisDeriv[iQ+jBasis*spaceDim  ];
	const double Njq = basisDeriv[iQ+jBasis*spaceDim+1];
	const double ki0j0 = 
	  C1111 * Nip * Njp + C1112 * Niq * Njp +
	  C1112 * Nip * Njq + C1212 * Niq * Njq;
	const double ki0j1 =
	  C1122 * Nip * Njq + C2212 * Niq * Njq +
	  C1112 * Nip * Njp + C1212 * Niq * Njp;
	const double ki1j0 =
	  C1122 * Niq * Njp + C2212 * Niq * Njq +
	  C1112 * Nip * Njp + C1212 * Nip * Njq;
	const double ki1j1 =
	  C2222 * Niq * Njq + C2212 * Nip * Njq +
	  C2212 * Niq * Njp + C1212 * Nip * Njp;
	const int jBlock = (jBasis*spaceDim  );
	const int jBlock1 = (jBasis*spaceDim+1);
	_cellMatrix[iBlock +jBlock ] += ki0j0;
	_cellMatrix[iBlock +jBlock1] += ki0j1;
	_cellMatrix[iBlock1+jBlock ] += ki1j0;
	_cellMatrix[iBlock1+jBlock1] += ki1j1;
      } // for
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(2+numBasis*(3*11+4))));
} // _elasticityJacobian2D

// ----------------------------------------------------------------------
// Integrate elasticity term in Jacobian for 3-D cells.
void
pylith::feassemble::IntegratorElasticity::_elasticityJacobian3D(
			       const double_array& elasticConsts)
{ // _elasticityJacobian3D
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const double_array& quadWts = _quadrature->quadWts();
  const double_array& jacobianDet = _quadrature->jacobianDet();
  const double_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(3 == cellDim);
  assert(quadWts.size() == numQuadPts);
  const int numConsts = 21;

  // Compute Jacobian for consistent tangent matrix
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    // tau_ij = C_ijkl * e_kl
    //        = C_ijlk * 0.5 (u_k,l + u_l,k)
    //        = 0.5 * C_ijkl * (u_k,l + u_l,k)
    // divide C_ijkl by 2 if k != l
    const double C1111 = elasticConsts[iQuad*numConsts+ 0];
    const double C1122 = elasticConsts[iQuad*numConsts+ 1];
    const double C1133 = elasticConsts[iQuad*numConsts+ 2];
    const double C1112 = elasticConsts[iQuad*numConsts+ 3] / 2.0;
    const double C1123 = elasticConsts[iQuad*numConsts+ 4] / 2.0;
    const double C1113 = elasticConsts[iQuad*numConsts+ 5] / 2.0;
    const double C2222 = elasticConsts[iQuad*numConsts+ 6];
    const double C2233 = elasticConsts[iQuad*numConsts+ 7];
    const double C2212 = elasticConsts[iQuad*numConsts+ 8] / 2.0;
    const double C2223 = elasticConsts[iQuad*numConsts+ 9] / 2.0;
    const double C2213 = elasticConsts[iQuad*numConsts+10] / 2.0;
    const double C3333 = elasticConsts[iQuad*numConsts+11];
    const double C3312 = elasticConsts[iQuad*numConsts+12] / 2.0;
    const double C3323 = elasticConsts[iQuad*numConsts+13] / 2.0;
    const double C3313 = elasticConsts[iQuad*numConsts+14] / 2.0;
    const double C1212 = elasticConsts[iQuad*numConsts+15] / 2.0;
    const double C1223 = elasticConsts[iQuad*numConsts+16] / 2.0;
    const double C1213 = elasticConsts[iQuad*numConsts+17] / 2.0;
    const double C2323 = elasticConsts[iQuad*numConsts+18] / 2.0;
    const double C2313 = elasticConsts[iQuad*numConsts+19] / 2.0;
    const double C1313 = elasticConsts[iQuad*numConsts+20] / 2.0;
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim;
	 iBasis < numBasis;
	 ++iBasis) {
      const double Nip = wt*basisDeriv[iQ+iBasis*spaceDim+0];
      const double Niq = wt*basisDeriv[iQ+iBasis*spaceDim+1];
      const double Nir = wt*basisDeriv[iQ+iBasis*spaceDim+2];
      for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	const double Njp = basisDeriv[iQ+jBasis*spaceDim+0];
	const double Njq = basisDeriv[iQ+jBasis*spaceDim+1];
	const double Njr = basisDeriv[iQ+jBasis*spaceDim+2];
	const double ki0j0 = 
	  C1111 * Nip * Njp + C1112 * Niq * Njp + C1113 * Nir * Njp +
	  C1112 * Nip * Njq + C1212 * Niq * Njq + C1213 * Nir * Njq +
	  C1113 * Nip * Njr + C1213 * Niq * Njr + C1313 * Nir * Njr;
	const double ki0j1 =
	  C1122 * Nip * Njq + C2212 * Niq * Njq + C2213 * Nir * Njq +
	  C1112 * Nip * Njp + C1212 * Niq * Njp + C1213 * Nir * Njp +
	  C1123 * Nip * Njr + C1223 * Niq * Njr + C2313 * Nir * Njr;
	const double ki0j2 =
	  C1133 * Nip * Njr + C3312 * Niq * Njr + C3313 * Nir * Njr +
	  C1123 * Nip * Njq + C1223 * Niq * Njq + C2313 * Nir * Njq +
	  C1113 * Nip * Njp + C1213 * Niq * Njp + C1313 * Nir * Njp;
	const double ki1j0 =
	  C1122 * Niq * Njp + C1112 * Nip * Njp + C1123 * Nir * Njp +
	  C2212 * Niq * Njq + C1212 * Nip * Njq + C1223 * Nir * Njq +
	  C2213 * Niq * Njr + C1213 * Nip * Njr + C2313 * Nir * Njr;
	const double ki1j1 =
	  C2222 * Niq * Njq + C2212 * Nip * Njq + C2223 * Nir * Njq +
	  C2212 * Niq * Njp + C1212 * Nip * Njp + C1223 * Nir * Njp +
	  C2223 * Niq * Njr + C1223 * Nip * Njr + C2323 * Nir * Njr;
	const double ki1j2 =
	  C2233 * Niq * Njr + C3312 * Nip * Njr + C3323 * Nir * Njr +
	  C2223 * Niq * Njq + C1223 * Nip * Njq + C2323 * Nir * Njq +
	  C2213 * Niq * Njp + C1213 * Nip * Njp + C2313 * Nir * Njp;
	const double ki2j0 =
	  C1133 * Nir * Njp + C1123 * Niq * Njp + C1113 * Nip * Njp +
	  C3312 * Nir * Njq + C1223 * Niq * Njq + C1213 * Nip * Njq +
	  C3313 * Nir * Njr + C2313 * Niq * Njr + C1313 * Nip * Njr; 
	const double ki2j1 =
	  C2233 * Nir * Njq + C2223 * Niq * Njq + C2213 * Nip * Njq +
	  C3312 * Nir * Njp + C1223 * Niq * Njp + C1213 * Nip * Njp +
	  C3323 * Nir * Njr + C2323 * Niq * Njr + C2313 * Nip * Njr; 
	const double ki2j2 =
	  C3333 * Nir * Njr + C3323 * Niq * Njr + C3313 * Nip * Njr +
	  C3323 * Nir * Njq + C2323 * Niq * Njq + C2313 * Nip * Njq +
	  C3313 * Nir * Njp + C2313 * Niq * Njp + C1313 * Nip * Njp;
	const int iBlock = iBasis*spaceDim * (numBasis*spaceDim);
	const int iBlock1 = (iBasis*spaceDim+1) * (numBasis*spaceDim);
	const int iBlock2 = (iBasis*spaceDim+2) * (numBasis*spaceDim);
	const int jBlock = jBasis*spaceDim;
	const int jBlock1 = jBasis*spaceDim+1;
	const int jBlock2 = jBasis*spaceDim+2;
	_cellMatrix[iBlock +jBlock ] += ki0j0;
	_cellMatrix[iBlock +jBlock1] += ki0j1;
	_cellMatrix[iBlock +jBlock2] += ki0j2;
	_cellMatrix[iBlock1+jBlock ] += ki1j0;
	_cellMatrix[iBlock1+jBlock1] += ki1j1;
	_cellMatrix[iBlock1+jBlock2] += ki1j2;
	_cellMatrix[iBlock2+jBlock ] += ki2j0;
	_cellMatrix[iBlock2+jBlock1] += ki2j1;
	_cellMatrix[iBlock2+jBlock2] += ki2j2;
      } // for
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(3+numBasis*(6*26+9))));
} // _elasticityJacobian3D

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D(
					    double_array* strain,
					    const double_array& basisDeriv,
					    const double_array& disp,
					    const int numBasis,
					    const int numQuadPts)
{ // calcTotalStrain1D
  assert(0 != strain);

  const int dim = 1;
  
  assert(basisDeriv.size() == numQuadPts*numBasis*dim);
  assert(disp.size() == numBasis*dim);

  (*strain) = 0.0;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    for (int iBasis=0; iBasis < numBasis; ++iBasis)
      (*strain)[iQuad] += basisDeriv[iQuad*numBasis+iBasis] * disp[iBasis];
  } // for
} // calcTotalStrain1D

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D(
					    double_array* strain,
					    const double_array& basisDeriv,
					    const double_array& disp,
					    const int numBasis,
					    const int numQuadPts)
{ // calcTotalStrain2D
  assert(0 != strain);
  
  const int dim = 2;

  assert(basisDeriv.size() == numQuadPts*numBasis*dim);
  assert(disp.size() == numBasis*dim);

  (*strain) = 0.0;
  const int strainSize = 3;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    for (int iBasis=0, iQ=iQuad*numBasis*dim; iBasis < numBasis; ++iBasis) {
      (*strain)[iQuad*strainSize+0] += 
	basisDeriv[iQ+iBasis*dim  ] * disp[iBasis*dim  ];
      (*strain)[iQuad*strainSize+1] += 
	basisDeriv[iQ+iBasis*dim+1] * disp[iBasis*dim+1];
      (*strain)[iQuad*strainSize+2] += 
	0.5 * (basisDeriv[iQ+iBasis*dim+1] * disp[iBasis*dim  ] +
	       basisDeriv[iQ+iBasis*dim  ] * disp[iBasis*dim+1]);
    } // for
} // calcTotalStrain2D

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D(
					    double_array* strain,
					    const double_array& basisDeriv,
					    const double_array& disp,
					    const int numBasis,
					    const int numQuadPts)
{ // calcTotalStrain3D
  assert(0 != strain);

  const int dim = 3;

  assert(basisDeriv.size() == numQuadPts*numBasis*dim);
  assert(disp.size() == numBasis*dim);

  (*strain) = 0.0;
  const int strainSize = 6;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    for (int iBasis=0, iQ=iQuad*numBasis*dim; iBasis < numBasis; ++iBasis) {
      (*strain)[iQuad*strainSize+0] += 
	basisDeriv[iQ+iBasis*dim  ] * disp[iBasis*dim  ];
      (*strain)[iQuad*strainSize+1] += 
	basisDeriv[iQ+iBasis*dim+1] * disp[iBasis*dim+1];
      (*strain)[iQuad*strainSize+2] += 
	basisDeriv[iQ+iBasis*dim+2] * disp[iBasis*dim+2];
      (*strain)[iQuad*strainSize+3] += 
	0.5 * (basisDeriv[iQ+iBasis*dim+1] * disp[iBasis*dim  ] +
	       basisDeriv[iQ+iBasis*dim  ] * disp[iBasis*dim+1]);
      (*strain)[iQuad*strainSize+4] += 
	0.5 * (basisDeriv[iQ+iBasis*dim+2] * disp[iBasis*dim+1] +
	       basisDeriv[iQ+iBasis*dim+1] * disp[iBasis*dim+2]);
      (*strain)[iQuad*strainSize+5] += 
	0.5 * (basisDeriv[iQ+iBasis*dim+2] * disp[iBasis*dim  ] +
	       basisDeriv[iQ+iBasis*dim  ] * disp[iBasis*dim+2]);
    } // for
} // calcTotalStrain3D


// End of file 
