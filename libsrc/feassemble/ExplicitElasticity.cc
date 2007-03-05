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

#include "ExplicitElasticity.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature
#include "ParameterManager.hh" // USES ParameterManager

#include "petscmat.h" // USES PetscMat
#include "spatialdata/spatialdb/SpatialDB.hh"

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ExplicitElasticity::ExplicitElasticity(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ExplicitElasticity::~ExplicitElasticity(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::ExplicitElasticity::ExplicitElasticity(const ExplicitElasticity& i) :
  IntegratorExplicit(i)
{ // copy constructor
  std::cout << "In ExplicitElasticity copy constructor." << std::endl;
} // copy constructor

// ----------------------------------------------------------------------
// Integrate residual term (b) for dynamic elasticity term for 3-D
// finite elements.
void
pylith::feassemble::ExplicitElasticity::integrateResidual(
			      const ALE::Obj<real_section_type>& residual,
			      const ALE::Obj<real_section_type>& dispT,
			      const ALE::Obj<real_section_type>& dispTmdt,
			      const ALE::Obj<real_section_type>& coordinates)
{ // integrateResidual
  assert(0 != _quadrature);

  // Get information about section
  const topology_type::patch_type patch = 0;
  const ALE::Obj<topology_type>& topology = dispT->getTopology();
  const ALE::Obj<topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator cellsEnd = cells->end();


  // Get parameters used in integration.
  const double dt = _dt;
  const ALE::Obj<real_section_type>& density = _parameters->getReal("density");

  // Allocate vector for cell values (if necessary)
  _initCellVector();

  for (topology_type::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(coordinates, *cellIter);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    const real_section_type::value_type* dispTCell = 
      dispT->restrict(patch, *cellIter);
    const real_section_type::value_type* dispTmdtCell = 
      dispTmdt->restrict(patch, *cellIter);

    // Get cell geometry information
    const int numQuadPts = _quadrature->numQuadPts();
    const double* basis = _quadrature->basis();
    const double* quadPts = _quadrature->quadPts();
    const double* quadWts = _quadrature->quadWts();
    const double* jacobianDet = _quadrature->jacobianDet();
    const int numBasis = _quadrature->numCorners();
    const int spaceDim = _quadrature->spaceDim();

    // Restrict material properties material database to quadrature 
    // points for this cell
    const real_section_type::value_type* densityCell = 
      density->restrict(patch, *cellIter);

    // Compute action for cell

    // Compute action for inertial terms
    const double dt2 = dt*dt;
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] * densityCell[iQuad] / dt2;
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const int iBlock = iBasis * spaceDim;
        const double valI = wt*basis[iQ+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const int jBlock = jBasis * spaceDim;
          const double valIJ = valI * basis[iQ+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBlock+iDim] += 
	      valIJ * 2.0 * (dispTCell[jBlock+iDim] - 
			     dispTmdtCell[jBlock+iDim]);
        } // for
      } // for
    } // for

    // Compute action for elastic terms
    // ADD STUFF HERE

    PetscErrorCode err = 
      PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(1+4*spaceDim))));
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");
    
    // Assemble cell contribution into field
    residual->updateAdd(patch, *cellIter, _cellVector);
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::ExplicitElasticity::integrateJacobian(
			     PetscMat* mat,
			     const ALE::Obj<real_section_type>& dispT,
			     const ALE::Obj<real_section_type>& coordinates)
{ // integrateJacobian
  assert(0 != mat);
  assert(0 != _quadrature);
  PetscErrorCode err;

  // Get information about section
  const topology_type::patch_type patch = 0;
  const ALE::Obj<topology_type>& topology = coordinates->getTopology();
  const ALE::Obj<topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<ALE::Mesh::order_type>& globalOrder = 
    ALE::New::NumberingFactory<topology_type>::singleton(
       topology->debug())->getGlobalOrder(topology, patch, 
					  dispT->getName(), 
					  dispT->getAtlas());

  // Setup symmetric, sparse matrix
  // :TODO: This needs to be moved outside Integrator object, because
  // integrator object will be specific to cell type and material type
  int localSize  = globalOrder->getLocalSize();
  int globalSize = globalOrder->getGlobalSize();
  err = MatCreate(topology->comm(), mat);
  err = MatSetSizes(*mat, localSize, localSize, globalSize, globalSize);
  err = MatSetFromOptions(*mat);
  err = preallocateMatrix(topology, dispT->getAtlas(), globalOrder, *mat);

  // Get parameters used in integration.
  const ALE::Obj<real_section_type>& density = _parameters->getReal("density");

  // Allocate matrix for cell values (if necessary)
  _initCellMatrix();

  for (topology_type::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(coordinates, *cellIter);

    // Reset element matrix to zero
    _resetCellMatrix();

    // Get cell geometry information
    const int numQuadPts = _quadrature->numQuadPts();
    const double* basis = _quadrature->basis();
    const double* quadPts = _quadrature->quadPts();
    const double* quadWts = _quadrature->quadWts();
    const double* jacobianDet = _quadrature->jacobianDet();
    const int numBasis = _quadrature->numCorners();
    const int spaceDim = _quadrature->spaceDim();

    // Restrict material properties material database to quadrature 
    // points for this cell
    const real_section_type::value_type* densityCell = 
      density->restrict(patch, *cellIter);

    // Integrate cell
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] * densityCell[iQuad];
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	const int iBlock = iBasis * spaceDim;
	const double valI = wt*basis[iQ+iBasis];
	for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	  const int jBlock = jBasis * spaceDim;
	  const double val = valI * basis[iQ+jBasis];
	  for (int iDim=0; iDim < spaceDim; ++iDim)
	    _cellMatrix[(iBlock+iDim)*(numBasis*spaceDim)+jBlock+iDim] += val;
	} // for
      } // for
    } // for
    err = PetscLogFlops(numQuadPts*(2+numBasis*(1+numBasis*(1+spaceDim))));
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");
    
    // Assemble cell contribution into sparse matrix
    err = updateOperator(*mat, dispT, globalOrder, *cellIter, _cellMatrix, 
			 ADD_VALUES);
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Compute lumped matrix associated with operator.
void
pylith::feassemble::ExplicitElasticity::integrateJacobian(
			     const ALE::Obj<real_section_type>& fieldOut,
			     const ALE::Obj<real_section_type>& dispT,
			     const ALE::Obj<real_section_type>& coordinates)
{ // integrateJacobian
  assert(0 != _quadrature);
  assert(0 != _parameters);

  // Get information about section
  const topology_type::patch_type patch = 0;
  const ALE::Obj<topology_type>& topology = coordinates->getTopology();
  const ALE::Obj<topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator cellsEnd = cells->end();

  const ALE::Obj<real_section_type>& density = _parameters->getReal("density");

  // Allocate matrix for cell values (if necessary)
  _initCellVector();

  for (topology_type::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(coordinates, *cellIter);

    // Reset element matrix to zero
    _resetCellVector();

    // Get cell geometry information
    const int numQuadPts = _quadrature->numQuadPts();
    const double* basis = _quadrature->basis();
    const double* quadPts = _quadrature->quadPts();
    const double* quadWts = _quadrature->quadWts();
    const double* jacobianDet = _quadrature->jacobianDet();
    const int numBasis = _quadrature->numCorners();
    const int spaceDim = _quadrature->spaceDim();

    // Restrict material properties material database to quadrature 
    // points for this cell
    const real_section_type::value_type* densityCell = 
      density->restrict(patch, *cellIter);

    // Compute lumped mass matrix for cell
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] * densityCell[iQuad];
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	const int iBlock = iBasis * spaceDim;
	const double valI = wt*basis[iQ+iBasis];
	for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	  const int jBlock = jBasis * spaceDim;
	  const double val = valI*basis[iQ+jBasis];
	  for (int iDim=0; iDim < spaceDim; ++iDim)
	    _cellVector[iBlock+iDim] += val;
	} // for
      } // for
    } // for

    PetscErrorCode err = 
      PetscLogFlops(numQuadPts*(2+numBasis*(1+numBasis*(1+spaceDim))));
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");
    
    // Assemble cell contribution into field
    fieldOut->updateAdd(patch, *cellIter, _cellVector);
  } // for
} // integrateLumped

// ----------------------------------------------------------------------
// Setup material property parameters by querying database.
void
pylith::feassemble::ExplicitElasticity::setupMatProp(ALE::Obj<ALE::Mesh>& mesh,
						     spatialdata::geocoords::CoordSys* cs,
						     spatialdata::spatialdb::SpatialDB* db)
{ // setupMatProp
  assert(0 != cs);
  assert(0 != db);
  assert(0 != _parameters);

  typedef ALE::Mesh::real_section_type real_section_type;
  typedef ALE::Mesh::topology_type topology_type;

  _parameters->addReal("density");
  const ALE::Obj<real_section_type>& density = _parameters->getReal("density");

  const int numQuadPts = _quadrature->numQuadPts();
  const ALE::Mesh::int_section_type::patch_type patch = 0;
  const int fiberDim = numQuadPts; // number of values in field per cell
  density->setFiberDimensionByDepth(patch, 0, fiberDim);
  density->allocate();

  // Open database
  db->open();
  const int numVals = 1;
  const char* names[numVals];
  names[0] = "density";
  db->queryVals(names, numVals);
  
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<topology_type>& topology = coordinates->getTopology();
  const ALE::Obj<topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator cellsEnd = cells->end();

  // Loop over cells
  double* cellDensity = (numQuadPts > 0) ? new double[numQuadPts] : 0;
  for (topology_type::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(coordinates, *cellIter);

    const double* quadPts = _quadrature->quadPts();
    const int spaceDim = _quadrature->spaceDim();

    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, index=0; 
	 iQuadPt < numQuadPts; 
	 ++iQuadPt, index+=spaceDim)
      // account for differences in spaceDim
      const int err = db->query(&cellDensity[iQuadPt], numVals, 
				&quadPts[index], spaceDim, cs);
    // Assemble cell contribution into field
    density->updateAdd(patch, *cellIter, cellDensity);
  } // for
  delete[] cellDensity; cellDensity = 0;

  // Close database
  db->close();
} // setupMatProp


// End of file 
