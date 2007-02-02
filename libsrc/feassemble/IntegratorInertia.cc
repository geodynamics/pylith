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

#include "IntegratorInertia.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature

#include "petscmat.h" // USES PetscMat
#include "spatialdata/spatialdb/SpatialDB.hh"

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorInertia::IntegratorInertia(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorInertia::~IntegratorInertia(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::IntegratorInertia::IntegratorInertia(const IntegratorInertia& i) :
  Integrator(i)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Integrate inertial term for 3-D finite elements.
void
pylith::feassemble::IntegratorInertia::integrateAction(
			      const ALE::Obj<real_section_type>& fieldOut,
			      const ALE::Obj<real_section_type>& fieldIn,
			      const ALE::Obj<real_section_type>& coordinates)
{ // integrateAction
  assert(0 != _quadrature);
  assert(0 != _density.isNull());

  // Get information about section
  const topology_type::patch_type patch = 0;
  const ALE::Obj<topology_type>& topology = fieldIn->getTopology();
  const ALE::Obj<topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator cellsEnd = cells->end();

  // Allocate vector for cell values (if necessary)
  _initCellVector();

  for (topology_type::label_sequence::iterator cellIter=cells->begin();
       cellIter != cellsEnd;
       ++cellIter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(coordinates, *cellIter);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input field to cell
    const real_section_type::value_type* fieldInCell = 
      fieldIn->restrict(patch, *cellIter);

    // Get cell geometry information
    const int numQuadPts = _quadrature->numQuadPts();
    const double* basis = _quadrature->basis();
    const double* quadPts = _quadrature->quadPts();
    const double* quadWts = _quadrature->quadWts();
    const double* jacobianDet = _quadrature->jacobianDet();
    const int numBasis = _quadrature->numCorners();
    const int spaceDim = _quadrature->spaceDim();

    // Restrict density from material database to quadrature points for this cell
    const real_section_type::value_type* density = 
      _density->restrict(patch, *cellIter);

    // Compute action for cell
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad];
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const int iBlock = iBasis * spaceDim;
        const double valI = wt*basis[iQ+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const int jBlock = jBasis * spaceDim;
          const double val = valI * basis[iQ+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBlock+iDim] += val * fieldInCell[jBlock+iDim];
        } // for
      } // for
    } // for
    PetscErrorCode err = 
      PetscLogFlops(numQuadPts*(2+numBasis*(1+numBasis*(1+2*spaceDim))));
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");
    
    // Assemble cell contribution into field
    fieldOut->updateAdd(patch, *cellIter, _cellVector);
  } // for
} // integrateAction

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::IntegratorInertia::integrate(
			     PetscMat* mat,
			     const ALE::Obj<real_section_type>& fieldIn,
			     const ALE::Obj<real_section_type>& coordinates)
{ // integrate
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
					  fieldIn->getName(), 
					  fieldIn->getAtlas());

  // Setup symmetric, sparse matrix
  // :TODO: This needs to be moved outside Integrator object, because
  // integrator object will be specific to cell type and material type
  int localSize  = globalOrder->getLocalSize();
  int globalSize = globalOrder->getGlobalSize();
  err = MatCreate(topology->comm(), mat);
  err = MatSetSizes(*mat, localSize, localSize, globalSize, globalSize);
  err = MatSetFromOptions(*mat);
  err = preallocateMatrix(topology, fieldIn->getAtlas(), globalOrder, *mat);

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

    // :TODO: Get mass density at quadrature points from material database
    // For now, hardwire mass density
    const double density = 1.0;

    // Integrate cell
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad] * density;
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
    err = updateOperator(*mat, fieldIn, globalOrder, *cellIter, _cellMatrix, 
			 ADD_VALUES);
  } // for
} // integrate

// ----------------------------------------------------------------------
// Compute lumped matrix associated with operator.
void
pylith::feassemble::IntegratorInertia::integrateLumped(
			     const ALE::Obj<real_section_type>& fieldOut,
			     const ALE::Obj<real_section_type>& coordinates)
{ // integrateLumped
  assert(0 != _quadrature);

  // Get information about section
  const topology_type::patch_type patch = 0;
  const ALE::Obj<topology_type>& topology = coordinates->getTopology();
  const ALE::Obj<topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator cellsEnd = cells->end();

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

    // :TODO: Get mass density at quadrature points from material database
    // For now, hardwire mass density
    const double density = 1.0;

    // Compute lumped mass matrix for cell
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad] * density;
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
// Initialize, get material property parameters from database.
void
pylith::feassemble::IntegratorInertia::initialize(
				     ALE::Obj<ALE::Mesh>& mesh,
				     spatialdata::geocoords::CoordSys* cs,
				     spatialdata::spatialdb::SpatialDB* db)
{ // initialize
  assert(0 != cs);
  assert(0 != db);

  typedef ALE::Mesh::real_section_type real_section_type;
  typedef ALE::Mesh::topology_type topology_type;

  // Create density section
  const int numQuadPts = _quadrature->numQuadPts();
  const ALE::Mesh::int_section_type::patch_type patch = 0;
  _density = mesh->getRealSection("density");
  const int fiberDim = numQuadPts; // number of values in field per cell
  _density->setName("density");
  _density->setFiberDimensionByDepth(patch, 0, fiberDim);
  _density->allocate();

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
    _density->updateAdd(patch, *cellIter, cellDensity);
  } // for
  delete[] cellDensity; cellDensity = 0;

  // Close database
  db->close();
} // initialize


// End of file 
