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

#include "IntegratorInertia3D.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature

#include "petscmat.h" // USES PetscMat

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorInertia3D::IntegratorInertia3D(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorInertia3D::~IntegratorInertia3D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::IntegratorInertia3D::IntegratorInertia3D(const IntegratorInertia3D& i) :
  Integrator(i)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Integrate inertial term for 3-D finite elements.
void
pylith::feassemble::IntegratorInertia3D::integrateAction(
			      const ALE::Obj<real_section_type>& fieldOut,
			      const ALE::Obj<real_section_type>& fieldIn,
			      const ALE::Obj<real_section_type>& coordinates)
{ // integrateAction
  assert(0 != _quadrature);

  const topology_type::patch_type patch = 0;
  const ALE::Obj<topology_type>& topology = fieldIn->getTopology();
  const ALE::Obj<topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator cellsEnd = cells->end();

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

    const int numQuadPts = _quadrature->numQuadPts();
    const double* basis = _quadrature->basis();
    const double* quadPts = _quadrature->quadPts();
    const double* quadWts = _quadrature->quadWts();
    const double* jacobianDet = _quadrature->jacobianDet();
    const int numBasis = _quadrature->numCorners();
    const int spaceDim = _quadrature->spaceDim();

    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad]*jacobianDet[iQuad];
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	const int iC = iBasis*spaceDim;
	double val = wt*basis[iQ+iBasis];
	for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	  val *= basis[iQ+jBasis];
	  for (int iDim=0; iDim < spaceDim; ++iDim)
	    _cellVector[iC+iDim] += val*fieldInCell[iC+iDim];
	} // for
      } // for
    } // for
    
    // Assemble cell contribution into field
    fieldOut->updateAdd(patch, *cellIter, _cellVector);
  } // for
} // integrateAction

// ----------------------------------------------------------------------
// Compute matrix associated with operator.
void
pylith::feassemble::IntegratorInertia3D::integrate(PetscMat* mat,
			      const ALE::Obj<real_section_type>& coordinates)
{ // integrate
} // integrate


// End of file 
