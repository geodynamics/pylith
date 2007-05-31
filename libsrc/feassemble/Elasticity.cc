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

#include "Elasticity.hh" // implementation of class methods

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
void
pylith::feassemble::Elasticity::calcTotalStrain1D(
					    std::vector<double_array>* strain,
					    const double_array& basisDeriv,
					    const double_array& disp,
					    const int numBasis)
{ // calcTotalStrain1D
  assert(0 != strain);
  
  const int dim = 1;
  const int numQuadPts = strain->size();

  assert(basisDeriv.size() == numQuadPts*numBasis*dim);
  assert(disp.size() == numBasis*dim);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    assert(1 == (*strain)[iQuad].size());
    (*strain)[iQuad] *= 0.0;
    for (int iBasis=0; iBasis < numBasis; ++iBasis)
      (*strain)[iQuad][0] += 
	basisDeriv[iQuad*numBasis+iBasis] * disp[iBasis];
  } // for
} // calcTotalStrain1D

// ----------------------------------------------------------------------
void
pylith::feassemble::Elasticity::calcTotalStrain2D(
					    std::vector<double_array>* strain,
					    const double_array& basisDeriv,
					    const double_array& disp,
					    const int numBasis)
{ // calcTotalStrain2D
  assert(0 != strain);
  
  const int dim = 2;
  const int numQuadPts = strain->size();

  assert(basisDeriv.size() == numQuadPts*numBasis*dim);
  assert(disp.size() == numBasis*dim);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    assert(3 == (*strain)[iQuad].size());
    (*strain)[iQuad] *= 0.0;
    for (int iBasis=0, iQ=iQuad*numBasis*dim; iBasis < numBasis; ++iBasis) {
      (*strain)[iQuad][0] += basisDeriv[iQ+iBasis*dim  ] * disp[iBasis*dim  ];
      (*strain)[iQuad][1] += basisDeriv[iQ+iBasis*dim+1] * disp[iBasis*dim+1];
      (*strain)[iQuad][2] += 
	0.5 * (basisDeriv[iQ+iBasis*dim+1] * disp[iBasis*dim  ] +
	       basisDeriv[iQ+iBasis*dim  ] * disp[iBasis*dim+1]);
    } // for
  } // for
} // calcTotalStrain2D

// ----------------------------------------------------------------------
void
pylith::feassemble::Elasticity::calcTotalStrain3D(
					    std::vector<double_array>* strain,
					    const double_array& basisDeriv,
					    const double_array& disp,
					    const int numBasis)
{ // calcTotalStrain3D
  assert(0 != strain);

  const int dim = 3;
  const int numQuadPts = strain->size();

  assert(basisDeriv.size() == numQuadPts*numBasis*dim);
  assert(disp.size() == numBasis*dim);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    assert(6 == (*strain)[iQuad].size());
    (*strain)[iQuad] *= 0.0;
    for (int iBasis=0, iQ=iQuad*numBasis*dim; iBasis < numBasis; ++iBasis) {
      (*strain)[iQuad][0] += basisDeriv[iQ+iBasis*dim  ] * disp[iBasis*dim  ];
      (*strain)[iQuad][1] += basisDeriv[iQ+iBasis*dim+1] * disp[iBasis*dim+1];
      (*strain)[iQuad][2] += basisDeriv[iQ+iBasis*dim+2] * disp[iBasis*dim+2];
      (*strain)[iQuad][3] += 
	0.5 * (basisDeriv[iQ+iBasis*dim+1] * disp[iBasis*dim  ] +
	       basisDeriv[iQ+iBasis*dim  ] * disp[iBasis*dim+1]);
      (*strain)[iQuad][4] += 
	0.5 * (basisDeriv[iQ+iBasis*dim+2] * disp[iBasis*dim+1] +
	       basisDeriv[iQ+iBasis*dim+1] * disp[iBasis*dim+2]);
      (*strain)[iQuad][5] += 
	0.5 * (basisDeriv[iQ+iBasis*dim+2] * disp[iBasis*dim  ] +
	       basisDeriv[iQ+iBasis*dim  ] * disp[iBasis*dim+2]);
    } // for
  } // for
} // calcTotalStrain3D


// End of file 
