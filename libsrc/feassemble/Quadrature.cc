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

#include "Quadrature.hh" // implementation of class methods

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature::Quadrature(void) :
  _pBasisFns(0),
  _pBasisFnsDeriv(0),
  _pQuadPts(0),
  _pQuadWts(0),
  _numDims(0),
  _numCorners(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature::~Quadrature(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Set basis functions and their derivatives and coordinates and
//   weights of the quadrature points.
void
pylith::feassemble::Quadrature::initialize(const double* pBasisFns,
					   const double* pBasisFnsDeriv,
					   const double* pQuadPts,
					   const double* pQuadWts,
					   const int numDims,
					   const int numCorners,
					   const int numQuadPts)
{ // initialize
  if (0 == pBasisFns ||
      0 == pBasisFnsDeriv ||
      0 == pQuadPts ||
      0 == pQuadWts ||
      numDims < 1 || numDims > 3)
    throw std::runtime_error("Basis functions and their derivatives "
			     "and quadrature points must all be specified.");

  int size = numCorners * numQuadPts;
  delete[] _pBasisFns;
  _pBasisFns = (size > 0) ? new double[size] : 0;

  size = numCorners * numQuadPts * numDims;
  delete[] _pBasisFnsDeriv;
  _pBasisFnsDeriv = (size > 0) ? new double[size] : 0;

  size = numQuadPts * numDims;
  delete[] _pQuadPts;
  _pQuadPts = (size > 0) ? new double[size] : 0;

  size = numQuadPts;
  delete[] _pQuadWts;
  _pQuadWts = (size > 0) ? new double[size] : 0;

  _numDims = numDims;
  _numCorners = numCorners;
  _numQuadPts = numQuadPts;
} // initialize

// End of file 
