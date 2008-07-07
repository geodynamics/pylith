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
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
#include "pylith/utils/array.hh" // USES double_array
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN

#include "petscmat.h" // USES PetscMat
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include <assert.h> // USES assert()
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
			      const ALE::Obj<real_section_type>& residual,
			      const double t,
			      topology::FieldsManager* const fields,
			      const ALE::Obj<Mesh>& mesh,
			      const spatialdata::geocoords::CoordSys* cs)
{ // integrateResidual
  /// Member prototype for _elasticityResidualXD()
  typedef void (pylith::feassemble::ElasticityImplicit::*elasticityResidual_fn_type)
    (const double_array&);
  
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(!residual.isNull());
  assert(0 != fields);
  assert(!mesh.isNull());

  static PetscLogEvent setupEvent = 0, cellGeomEvent = 0, stateVarsEvent = 0, restrictEvent = 0, computeEvent = 0, updateEvent = 0, stressEvent;

  if (!setupEvent)
    PetscLogEventRegister("IRSetup", 0, &setupEvent);
  if (!cellGeomEvent)
    PetscLogEventRegister("IRCellGeom", 0, &cellGeomEvent);
  if (!stateVarsEvent)
    PetscLogEventRegister("IRProperties", 0, &stateVarsEvent);
  if (!restrictEvent)
    PetscLogEventRegister("IRRestrict", 0, &restrictEvent);
  if (!computeEvent)
    PetscLogEventRegister("IRCompute", 0, &computeEvent);
  if (!updateEvent)
    PetscLogEventRegister("IRUpdate", 0, &updateEvent);
  if (!stressEvent)
    PetscLogEventRegister("IRMaterialStress", 0, &stressEvent);

  const Obj<sieve_type>& sieve = mesh->getSieve();

  PetscLogEventBegin(setupEvent,0,0,0,0);
  // Set variables dependent on dimension of cell
  const int cellDim = _quadrature->cellDim();
  int tensorSize = 0;
  totalStrain_fn_type calcTotalStrainFn;
  elasticityResidual_fn_type elasticityResidualFn;
  if (1 == cellDim) {
    tensorSize = 1;
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityResidual1D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    tensorSize = 3;
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityResidual2D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    tensorSize = 6;
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityResidual3D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else
    assert(0);

  // Get cell information
  const int materialId = _material->id();
  const ALE::Obj<Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", materialId);
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

  // Precompute the geometric and function space information
  _quadrature->precomputeGeometry(mesh, coordinates, cells);

  // Allocate vector for cell values.
  _initCellVector();
  const int cellVecSize = numBasis*spaceDim;
  double_array dispTBctpdtCell(cellVecSize);

  // Allocate vector for total strain
  double_array totalStrain(numQuadPts*tensorSize);
  totalStrain = 0.0;
  PetscLogEventEnd(setupEvent,0,0,0,0);

  // Loop over cells
  int c_index = 0;
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++c_index) {
    // Compute geometry information for current cell
    PetscLogEventBegin(cellGeomEvent,0,0,0,0);
    _quadrature->retrieveGeometry(mesh, coordinates, *c_iter, c_index);
    PetscLogEventEnd(cellGeomEvent,0,0,0,0);

    // Get state variables for cell.
    PetscLogEventBegin(stateVarsEvent,0,0,0,0);
    _material->getPropertiesCell(*c_iter, numQuadPts);
    PetscLogEventEnd(stateVarsEvent,0,0,0,0);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    PetscLogEventBegin(restrictEvent,0,0,0,0);
    mesh->restrictClosure(dispTBctpdt, *c_iter, &dispTBctpdtCell[0], cellVecSize);
    PetscLogEventEnd(restrictEvent,0,0,0,0);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();
    const double_array& quadPts = _quadrature->quadPts();

    if (cellDim != spaceDim)
      throw std::logic_error("Not implemented yet.");


    // Compute body force vector if gravity is being used.
    if (0 != _gravityField) {

      // Make sure coordinate names exist in gravity field.
      _gravityField->open();
      if (1 == spaceDim){
	const char* queryNames[] = { "x"};
	_gravityField->queryVals(queryNames, spaceDim);
      } else if (2 == spaceDim){
	const char* queryNames[] = { "x", "y"};
	_gravityField->queryVals(queryNames, spaceDim);
      } else if (3 == spaceDim){
        const char* queryNames[] = { "x", "y", "z"};
	_gravityField->queryVals(queryNames, spaceDim);
      } else
	assert(0);

      // Get density at quadrature points for this cell
      const double_array& density = _material->calcDensity();

      // Compute action for element body forces
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	double gravVec[spaceDim];
	double coords[spaceDim];
	memcpy(coords, &quadPts[iQuad * spaceDim], sizeof(double)*spaceDim);

	const int err = _gravityField->query(gravVec, spaceDim,
				  coords, spaceDim, cs);
	if (err)
	  throw std::runtime_error("Unable to get gravity vector for point.");
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
      _gravityField->close();
    } // if

    // Compute B(transpose) * sigma, first computing strains
    PetscLogEventBegin(stressEvent,0,0,0,0);
    calcTotalStrainFn(&totalStrain, basisDeriv, dispTBctpdtCell, 
		      numBasis, numQuadPts);
    const double_array& stress = _material->calcStress(totalStrain, true);
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
    mesh->updateAdd(residual, *c_iter, _cellVector);
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
    (const double_array&);

  assert(0 != _quadrature);
  assert(0 != _material);
  assert(0 != mat);
  assert(0 != fields);
  assert(!mesh.isNull());
  typedef ALE::ISieveVisitor::IndicesVisitor<Mesh::real_section_type,Mesh::order_type,PetscInt> visitor_type;

  // Set variables dependent on dimension of cell
  const int cellDim = _quadrature->cellDim();
  int tensorSize = 0;
  totalStrain_fn_type calcTotalStrainFn;
  elasticityJacobian_fn_type elasticityJacobianFn;
  if (1 == cellDim) {
    tensorSize = 1;
    elasticityJacobianFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityJacobian1D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    tensorSize = 3;
    elasticityJacobianFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityJacobian2D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    tensorSize = 6;
    elasticityJacobianFn = 
      &pylith::feassemble::ElasticityImplicit::_elasticityJacobian3D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else
    assert(0);

  // Get cell information
  const int materialId = _material->id();
  const ALE::Obj<Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", materialId);
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
  double_array totalStrain(numQuadPts*tensorSize);
  totalStrain = 0.0;

  const ALE::Obj<Mesh::order_type>& globalOrder = mesh->getFactory()->getGlobalOrder(mesh, "default", dispTBctpdt);
  assert(!globalOrder.isNull());
  visitor_type iV(*dispTBctpdt, *globalOrder, (int) pow(mesh->getSieve()->getMaxConeSize(), mesh->depth())*spaceDim);

  // Loop over cells
  int c_index = 0;
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++c_index) {
    // Compute geometry information for current cell
    _quadrature->retrieveGeometry(mesh, coordinates, *c_iter, c_index);

    // Get state variables for cell.
    _material->getPropertiesCell(*c_iter, numQuadPts);

    // Reset element vector to zero
    _resetCellMatrix();

    // Restrict input fields to cell
    mesh->restrictClosure(dispTBctpdt, *c_iter, &dispTBctpdtCell[0], cellVecSize);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute strains
    calcTotalStrainFn(&totalStrain, basisDeriv, dispTBctpdtCell, 
		      numBasis, numQuadPts);
      
    // Get "elasticity" matrix at quadrature points for this cell
    const double_array& elasticConsts = 
      _material->calcDerivElastic(totalStrain);

    CALL_MEMBER_FN(*this, elasticityJacobianFn)(elasticConsts);

    // Assemble cell contribution into field.  Not sure if this is correct for
    // global stiffness matrix.
    PetscErrorCode err = updateOperator(*mat, *mesh->getSieve(), iV, *c_iter, _cellMatrix, ADD_VALUES);
    if (err)
      throw std::runtime_error("Update to PETSc Mat failed.");
    iV.clear();
  } // for
  _needNewJacobian = false;
  _material->resetNeedNewJacobian();
} // integrateJacobian


// End of file 
