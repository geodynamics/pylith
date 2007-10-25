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
#include "Elasticity.hh" // USES Elasticity

#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
#include "pylith/utils/array.hh" // USES double_array
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN

#include "petscmat.h" // USES PetscMat
#include "spatialdata/spatialdb/SpatialDB.hh"

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error

#define FASTER

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ElasticityImplicit::ElasticityImplicit(void) :
  _dtm1(-1.0),
  _material(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ElasticityImplicit::~ElasticityImplicit(void)
{ // destructor
  _material = 0; // Don't manage memory for material
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
  assert(_dt == _dtm1); // For now, don't allow variable time step
  if (0 != _material)
    _material->timeStep(_dt);
} // timeStep

// ----------------------------------------------------------------------
// Set material.
void
pylith::feassemble::ElasticityImplicit::material(materials::ElasticMaterial* m)
{ // material
  _material = m;
  if (0 != _material)
    _material->timeStep(_dt);  
} // material

// ----------------------------------------------------------------------
// Determine whether we need to recompute the Jacobian.
bool
pylith::feassemble::ElasticityImplicit::needNewJacobian(void)
{ // needNewJacobian
  assert(0 != _material);
  if (!_needNewJacobian)
    _needNewJacobian = _material->needNewJacobian();
  return _needNewJacobian;
} // needNewJacobian

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
			      const ALE::Obj<real_section_type>& residual,
			      const double t,
			      topology::FieldsManager* const fields,
			      const ALE::Obj<Mesh>& mesh)
{ // integrateResidual
  /// Member prototype for _elasticityResidualXD()
  typedef void (pylith::feassemble::ElasticityImplicit::*elasticityResidual_fn_type)
    (const std::vector<double_array>&);
  
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(!residual.isNull());
  assert(0 != fields);
  assert(!mesh.isNull());

  static PetscEvent setupEvent = 0, cellGeomEvent = 0, stateVarsEvent = 0, restrictEvent = 0, computeEvent = 0, updateEvent = 0, stressEvent;

  if (!setupEvent)
    PetscLogEventRegister(&setupEvent, "IRSetup", 0);
  if (!cellGeomEvent)
    PetscLogEventRegister(&cellGeomEvent, "IRCellGeom", 0);
  if (!stateVarsEvent)
    PetscLogEventRegister(&stateVarsEvent, "IRStateVars", 0);
  if (!restrictEvent)
    PetscLogEventRegister(&restrictEvent, "IRRestrict", 0);
  if (!computeEvent)
    PetscLogEventRegister(&computeEvent, "IRCompute", 0);
  if (!updateEvent)
    PetscLogEventRegister(&updateEvent, "IRUpdate", 0);
  if (!stressEvent)
    PetscLogEventRegister(&stressEvent, "IRMaterialStress", 0);

  const Obj<sieve_type>& sieve = mesh->getSieve();

  PetscLogEventBegin(setupEvent,0,0,0,0);
  // Set variables dependent on dimension of cell
  const int cellDim = _quadrature->cellDim();
  int tensorSize = 0;
  Elasticity::totalStrain_fn_type calcTotalStrainFn;
  elasticityResidual_fn_type elasticityResidualFn;
  if (1 == cellDim) {
    tensorSize = 1;
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityResidual1D;
    calcTotalStrainFn = &pylith::feassemble::Elasticity::calcTotalStrain1D;
  } else if (2 == cellDim) {
    tensorSize = 3;
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityResidual2D;
    calcTotalStrainFn = &pylith::feassemble::Elasticity::calcTotalStrain2D;
  } else if (3 == cellDim) {
    tensorSize = 6;
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityResidual3D;
    calcTotalStrainFn = &pylith::feassemble::Elasticity::calcTotalStrain3D;
  } else
    assert(0);

  // Get cell information
  const ALE::Obj<ALE::Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", _material->id());
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator  cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const ALE::Obj<real_section_type>& dispTBctpdt = 
    fields->getReal("dispTBctpdt");
  assert(!dispTBctpdt.isNull());

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

#ifdef FASTER
  static std::map<int, int> tags;
  int                       c = 0;

  if (tags.find(_material->id()) == tags.end()) {
    tags[_material->id()] = mesh->calculateCustomAtlas(dispTBctpdt, cells);
    residual->copyCustomAtlas(dispTBctpdt, tags[_material->id()]);
  }
  const int tag = tags[_material->id()];
#endif
  // Precompute the geometric and function space information
  _quadrature->precomputeGeometry(mesh, coordinates, cells);

  // Allocate vector for cell values.
  _initCellVector();
  const int cellVecSize = numBasis*spaceDim;
  double_array dispTBctpdtCell(cellVecSize);
  //double_array gravCell(cellVecSize);

  // Allocate vector for total strain
  std::vector<double_array> totalStrain(numQuadPts);
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    totalStrain[iQuad].resize(tensorSize);
    totalStrain[iQuad] = 0.0;
  } // for
  PetscLogEventEnd(setupEvent,0,0,0,0);

  // Loop over cells
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    PetscLogEventBegin(cellGeomEvent,0,0,0,0);
    _quadrature->retrieveGeometry(mesh, coordinates, *c_iter);
    PetscLogEventEnd(cellGeomEvent,0,0,0,0);

    // Get state variables for cell.
    PetscLogEventBegin(stateVarsEvent,0,0,0,0);
    _material->getStateVarsCell(*c_iter, numQuadPts);
    PetscLogEventEnd(stateVarsEvent,0,0,0,0);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    PetscLogEventBegin(restrictEvent,0,0,0,0);
#ifdef FASTER
    mesh->restrict(dispTBctpdt, tag, c, &dispTBctpdtCell[0], cellVecSize);
#else
    mesh->restrict(dispTBctpdt, *c_iter, &dispTBctpdtCell[0], cellVecSize);
#endif
    PetscLogEventEnd(restrictEvent,0,0,0,0);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    if (cellDim != spaceDim)
      throw std::logic_error("Not implemented yet.");

#if 0
    // Comment out gravity section for now, until we figure out how to deal
    // with gravity vector.

    // Get density at quadrature points for this cell
    const std::vector<double_array>& density = _material->calcDensity();

    // Compute action for element body forces
    if (!grav.isNull()) {
      mesh->restrict(grav, *c_iter, &gravCell[0], cellVecSize);
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	const double wt = quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad];
	for (int iBasis=0, iQ=iQuad*numBasis*cellDim;
	     iBasis < numBasis; ++iBasis) {
	  const double valI = wt*basis[iQ+iBasis];
	  for (int iDim=0; iDim < spaceDim; ++iDim) {
	    _cellVector[iBasis*spaceDim+iDim] += valI*gravCell[iDim];
	  } // for
	} // for
      } // for
      PetscLogFlopsNoCheck(numQuadPts*(2+numBasis*(2+2*spaceDim)));
    } // if
#endif

    // Compute B(transpose) * sigma, first computing strains
    PetscLogEventBegin(stressEvent,0,0,0,0);
    calcTotalStrainFn(&totalStrain, basisDeriv, dispTBctpdtCell, numBasis);
    const std::vector<double_array>& stress = 
      _material->calcStress(totalStrain);
    PetscLogEventEnd(stressEvent,0,0,0,0);

    PetscLogEventBegin(computeEvent,0,0,0,0);
    CALL_MEMBER_FN(*this, elasticityResidualFn)(stress);
    PetscLogEventEnd(computeEvent,0,0,0,0);

#if 0
    std::cout << "Updating residual for cell " << *c_iter << std::endl;
    for(int i = 0; i < _quadrature->spaceDim() * _quadrature->numBasis(); ++i) {
      std::cout << "  v["<<i<<"]: " << _cellVector[i] << std::endl;
    }
#endif
    // Assemble cell contribution into field
    PetscLogEventBegin(updateEvent,0,0,0,0);
#ifdef FASTER
    mesh->updateAdd(residual, tag, c++, _cellVector);
#else
    mesh->updateAdd(residual, *c_iter, _cellVector);
#endif
    PetscLogEventEnd(updateEvent,0,0,0,0);
  } // for
} // integrateResidual


// ----------------------------------------------------------------------
// Compute stiffness matrix.
void
pylith::feassemble::ElasticityImplicit::integrateJacobian(
					PetscMat* mat,
					const double t,
					topology::FieldsManager* fields,
					const ALE::Obj<Mesh>& mesh)
{ // integrateJacobian
  /// Member prototype for _elasticityJacobianXD()
  typedef void (pylith::feassemble::ElasticityImplicit::*elasticityJacobian_fn_type)
    (const std::vector<double_array>&);

  assert(0 != _quadrature);
  assert(0 != _material);
  assert(0 != mat);
  assert(0 != fields);
  assert(!mesh.isNull());

  // Set variables dependent on dimension of cell
  const int cellDim = _quadrature->cellDim();
  int tensorSize = 0;
  Elasticity::totalStrain_fn_type calcTotalStrainFn;
  elasticityJacobian_fn_type elasticityJacobianFn;
  if (1 == cellDim) {
    tensorSize = 1;
    elasticityJacobianFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityJacobian1D;
    calcTotalStrainFn = &pylith::feassemble::Elasticity::calcTotalStrain1D;
  } else if (2 == cellDim) {
    tensorSize = 3;
    elasticityJacobianFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityJacobian2D;
    calcTotalStrainFn = &pylith::feassemble::Elasticity::calcTotalStrain2D;
  } else if (3 == cellDim) {
    tensorSize = 6;
    elasticityJacobianFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityJacobian3D;
    calcTotalStrainFn = &pylith::feassemble::Elasticity::calcTotalStrain3D;
  } else
    assert(0);

  // Get cell information
  const ALE::Obj<ALE::Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", _material->id());
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator  cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const ALE::Obj<real_section_type>& dispTBctpdt = 
    fields->getReal("dispTBctpdt");
  assert(!dispTBctpdt.isNull());

  // Get parameters used in integration.
  const double dt = _dt;
  assert(dt > 0);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  
  if (cellDim != spaceDim)
    throw std::logic_error("Don't know how to integrate elasticity " \
			   "contribution to Jacobian matrix for cells with " \
			   "different dimensions than the spatial dimension.");

  // Precompute the geometric and function space information
  _quadrature->precomputeGeometry(mesh, coordinates, cells);

  // Allocate matrix and vectors for cell values.
  _initCellMatrix();
  const int cellVecSize = numBasis*spaceDim;
  double_array dispTBctpdtCell(cellVecSize);

  // Allocate vector for total strain
  std::vector<double_array> totalStrain(numQuadPts);
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    totalStrain[iQuad].resize(tensorSize);
    totalStrain[iQuad] = 0.0;
  } // for

  // Loop over cells
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    _quadrature->retrieveGeometry(mesh, coordinates, *c_iter);

    // Get state variables for cell.
    _material->getStateVarsCell(*c_iter, numQuadPts);

    // Reset element vector to zero
    _resetCellMatrix();

    // Restrict input fields to cell
    mesh->restrict(dispTBctpdt, *c_iter, &dispTBctpdtCell[0], cellVecSize);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute strains
    calcTotalStrainFn(&totalStrain, basisDeriv, dispTBctpdtCell, numBasis);
      
    // Get "elasticity" matrix at quadrature points for this cell
    const std::vector<double_array>& elasticConsts = 
      _material->calcDerivElastic(totalStrain);

    CALL_MEMBER_FN(*this, elasticityJacobianFn)(elasticConsts);

    // Assemble cell contribution into field.  Not sure if this is correct for
    // global stiffness matrix.
    const ALE::Obj<Mesh::order_type>& globalOrder = 
      mesh->getFactory()->getGlobalOrder(mesh, "default", dispTBctpdt);
    PetscErrorCode err = updateOperator(*mat, mesh, dispTBctpdt, globalOrder,
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
pylith::feassemble::ElasticityImplicit::updateState(
				   const double t,
				   const ALE::Obj<real_section_type>& disp,
				   const ALE::Obj<Mesh>& mesh)
{ // updateState
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(!disp.isNull());

  // No need to update state if using elastic behavior
  if (!_material->usesUpdateState())
    return;

  // Set variables dependent on dimension of cell
  const int cellDim = _quadrature->cellDim();
  int tensorSize = 0;
  Elasticity::totalStrain_fn_type calcTotalStrainFn;
  if (1 == cellDim) {
    tensorSize = 1;
    calcTotalStrainFn = &pylith::feassemble::Elasticity::calcTotalStrain1D;
  } else if (2 == cellDim) {
    tensorSize = 3;
    calcTotalStrainFn = &pylith::feassemble::Elasticity::calcTotalStrain2D;
  } else if (3 == cellDim) {
    tensorSize = 6;
    calcTotalStrainFn = &pylith::feassemble::Elasticity::calcTotalStrain3D;
  } else
    assert(0);

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

  // Precompute the geometric and function space information
  _quadrature->precomputeGeometry(mesh, coordinates, cells);

  const int cellVecSize = numBasis*spaceDim;
  double_array dispCell(cellVecSize);

  // Allocate vector for total strain
  std::vector<double_array> totalStrain(numQuadPts);
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    totalStrain[iQuad].resize(tensorSize);
    totalStrain[iQuad] = 0.0;
  } // for

  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    _quadrature->retrieveGeometry(mesh, coordinates, *c_iter);

    // Restrict input fields to cell
    mesh->restrict(disp, *c_iter, &dispCell[0], cellVecSize);

    // Get cell geometry information that depends on cell
    const double_array& basisDeriv = _quadrature->basisDeriv();
  
    // Compute strains
    calcTotalStrainFn(&totalStrain, basisDeriv, dispCell, numBasis);

    // Update material state
    _material->updateState(totalStrain, *c_iter);
  } // for

  _material->useElasticBehavior(false);
} // updateState

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::feassemble::ElasticityImplicit::verifyConfiguration(
						 const ALE::Obj<Mesh>& mesh)
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
  const int numCorners = _quadrature->numBasis();
  const ALE::Obj<ALE::Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", _material->id());
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  assert(!sieve.isNull());
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    const int cellNumCorners = sieve->nCone(*c_iter, mesh->depth())->size();
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
// Integrate elasticity term in residual for 1-D cells.
void
pylith::feassemble::ElasticityImplicit::_elasticityResidual1D(
				     const std::vector<double_array>& stress)
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
    const double s11 = stress[iQuad][0];
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const double N1 = wt*basisDeriv[iQuad*numBasis+iBasis  ];
      _cellVector[iBasis*spaceDim  ] -= N1*s11;
    } // for
  } // for
  PetscLogFlopsNoCheck(numQuadPts*(1+numBasis*5));
} // _elasticityResidual1D

// ----------------------------------------------------------------------
// Integrate elasticity term in residual for 2-D cells.
void
pylith::feassemble::ElasticityImplicit::_elasticityResidual2D(
				     const std::vector<double_array>& stress)
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
  
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    const double s11 = stress[iQuad][0];
    const double s22 = stress[iQuad][1];
    const double s12 = stress[iQuad][2];
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim;
	 iBasis < numBasis;
	 ++iBasis) {
      const double N1 = wt*basisDeriv[iQ+iBasis*spaceDim  ];
      const double N2 = wt*basisDeriv[iQ+iBasis*spaceDim+1];
      _cellVector[iBasis*spaceDim  ] -= N1*s11 + N2*s12;
      _cellVector[iBasis*spaceDim+1] -= N1*s12 + N2*s22;
    } // for
  } // for
  PetscLogFlopsNoCheck(numQuadPts*(1+numBasis*(8+2+9)));
} // _elasticityResidual2D

// ----------------------------------------------------------------------
// Integrate elasticity term in residual for 3-D cells.
void
pylith::feassemble::ElasticityImplicit::_elasticityResidual3D(
				     const std::vector<double_array>& stress)
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
  
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    const double s11 = stress[iQuad][0];
    const double s22 = stress[iQuad][1];
    const double s33 = stress[iQuad][2];
    const double s12 = stress[iQuad][3];
    const double s23 = stress[iQuad][4];
    const double s13 = stress[iQuad][5];
    
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
  PetscLogFlopsNoCheck(numQuadPts*(1+numBasis*(3+12)));
} // _elasticityResidual3D

// ----------------------------------------------------------------------
// Integrate elasticity term in Jacobian for 1-D cells.
void
pylith::feassemble::ElasticityImplicit::_elasticityJacobian1D(
			       const std::vector<double_array>& elasticConsts)
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
    const double C1111 = elasticConsts[iQuad][0];
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
  PetscLogFlopsNoCheck(numQuadPts*(1+numBasis*(2+numBasis*3)));
} // _elasticityJacobian1D

// ----------------------------------------------------------------------
// Integrate elasticity term in Jacobian for 2-D cells.
void
pylith::feassemble::ElasticityImplicit::_elasticityJacobian2D(
			       const std::vector<double_array>& elasticConsts)
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
  
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    // tau_ij = C_ijkl * e_kl
    //        = C_ijlk * 0.5 (u_k,l + u_l,k)
    //        = 0.5 * C_ijkl * (u_k,l + u_l,k)
    // divide C_ijkl by 2 if k != l
    const double C1111 = elasticConsts[iQuad][0];
    const double C1122 = elasticConsts[iQuad][1];
    const double C1112 = elasticConsts[iQuad][2]/2.0;
    const double C2222 = elasticConsts[iQuad][3];
    const double C2212 = elasticConsts[iQuad][4]/2.0;
    const double C1212 = elasticConsts[iQuad][5]/2.0;
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
  PetscLogFlopsNoCheck(numQuadPts*(1+numBasis*(2+numBasis*(3*11+4))));
} // _elasticityJacobian2D

// ----------------------------------------------------------------------
// Integrate elasticity term in Jacobian for 3-D cells.
void
pylith::feassemble::ElasticityImplicit::_elasticityJacobian3D(
			       const std::vector<double_array>& elasticConsts)
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
  
  // Compute Jacobian for consistent tangent matrix
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double wt = quadWts[iQuad] * jacobianDet[iQuad];
    // tau_ij = C_ijkl * e_kl
    //        = C_ijlk * 0.5 (u_k,l + u_l,k)
    //        = 0.5 * C_ijkl * (u_k,l + u_l,k)
    // divide C_ijkl by 2 if k != l
    const double C1111 = elasticConsts[iQuad][ 0];
    const double C1122 = elasticConsts[iQuad][ 1];
    const double C1133 = elasticConsts[iQuad][ 2];
    const double C1112 = elasticConsts[iQuad][ 3]/2.0;
    const double C1123 = elasticConsts[iQuad][ 4]/2.0;
    const double C1113 = elasticConsts[iQuad][ 5]/2.0;
    const double C2222 = elasticConsts[iQuad][ 6];
    const double C2233 = elasticConsts[iQuad][ 7];
    const double C2212 = elasticConsts[iQuad][ 8]/2.0;
    const double C2223 = elasticConsts[iQuad][ 9]/2.0;
    const double C2213 = elasticConsts[iQuad][10]/2.0;
    const double C3333 = elasticConsts[iQuad][11];
    const double C3312 = elasticConsts[iQuad][12]/2.0;
    const double C3323 = elasticConsts[iQuad][13]/2.0;
    const double C3313 = elasticConsts[iQuad][14]/2.0;
    const double C1212 = elasticConsts[iQuad][15]/2.0;
    const double C1223 = elasticConsts[iQuad][16]/2.0;
    const double C1213 = elasticConsts[iQuad][17]/2.0;
    const double C2323 = elasticConsts[iQuad][18]/2.0;
    const double C2313 = elasticConsts[iQuad][19]/2.0;
    const double C1313 = elasticConsts[iQuad][20]/2.0;
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
  PetscLogFlopsNoCheck(numQuadPts*(1+numBasis*(3+numBasis*(6*26+9))));
} // _elasticityJacobian3D


// End of file 
