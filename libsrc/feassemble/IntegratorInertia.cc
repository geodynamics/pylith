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

    // FIX THIS
    // Hardwire mass density
    const double density = 1.0;

    // Compute action for cell
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad] * density;
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	const int iBlock = iBasis * spaceDim;
	double val = wt*basis[iQ+iBasis];
	for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	  val *= basis[iQ+jBasis];
	  for (int iDim=0; iDim < spaceDim; ++iDim)
	    _cellVector[iBlock+iDim] += val * fieldInCell[iBlock+iDim];
	} // for
      } // for
    } // for
    PetscErrorCode err = 
      PetscLogFlops(numQuadPts*(2+numBasis*(1+numBasis*(1+2*spaceDim))));
    if (err)
      throw std::runtime_error("Couldn't log PETSc flops.");
    
    // Assemble cell contribution into field
    fieldOut->updateAdd(patch, *cellIter, _cellVector);
  } // for
} // integrateAction

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::IntegratorInertia::integrate(
			     PetscMat* mat,
			     const ALE::Obj<real_section_type>& coordinates)
{ // integrate
  assert(0 != mat);
  assert(0 != _quadrature);

  // Get information about section
  const topology_type::patch_type patch = 0;
  const ALE::Obj<topology_type>& topology = coordinates->getTopology();
  const ALE::Obj<topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator cellsEnd = cells->end();

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

    // FIX THIS
    // Hardwire mass density
    const double density = 1.0;

    // Integrate cell
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad] * density;
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	const int iBlock = iBasis * spaceDim;
	double val = wt*basis[iQ+iBasis];
	for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	  const int jBlock = jBasis * spaceDim;
	  val *= basis[iQ+jBasis];
	  for (int iDim=0; iDim < spaceDim; ++iDim)
	    _cellMatrix[(iBlock+iDim)*(numBasis*spaceDim)+jBlock+iDim] += val;
	} // for
      } // for
    } // for
    PetscErrorCode err = 
      PetscLogFlops(numQuadPts*(2+numBasis*(1+numBasis*(1+spaceDim))));
    if (err)
      throw std::runtime_error("Couldn't log PETSc flops.");
    
    // Assemble cell contribution into sparse matrix
    // STUFF GOES HERE
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

    // FIX THIS
    // Hardwire mass density
    const double density = 1.0;

    // Compute lumped mass matrix for cell
    double sumdiag = 0;
    double diagScale = 0;
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad] * density;
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	const int iBlock = iBasis * spaceDim;
	const double val = wt*basis[iQ+iBasis]*basis[iQ+iBasis];
	for (int iDim=0; iDim < spaceDim; ++iDim)
	  _cellVector[iBlock+iDim] += val;
	sumdiag += val;
      } // for
      diagScale += numBasis*wt;
    } // for
    diagScale /= sumdiag*numBasis;

    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const int iBlock = iBasis * spaceDim;
      for (int iDim=0; iDim < numBasis; ++iDim)
	_cellVector[iBlock+iDim] *= diagScale;
    } // for
    PetscErrorCode err = 
      PetscLogFlops(numQuadPts*(2+numBasis*2) + numBasis);
    if (err)
      throw std::runtime_error("Couldn't log PETSc flops.");
    
    // Assemble cell contribution into field
    fieldOut->updateAdd(patch, *cellIter, _cellVector);
  } // for
} // integrateLumped


// End of file 
