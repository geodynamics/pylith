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

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::logic_error

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticity::calcTotalStrain1D(
					    std::vector<double_array>* strain,
					    const double_array& basisDeriv,
					    const double* disp,
					    const int numBasis)
{ // calcTotalStrain1D
  assert(0 != strain);
  assert(0 != disp);
  
  const int dimension = 1;
  const int numQuadPts = strain->size();

  assert(basisDeriv.size() == numQuadPts*numBasis*dimension);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    assert(1 == (*strain)[iQuad].size());
    for (int iBasis=0; iBasis < numBasis; ++iBasis)
      (*strain)[iQuad][0] += 
	basisDeriv[iQuad*numBasis+iBasis] * disp[iBasis];
  } // for
} // calcTotalStrain1D

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticity::calcTotalStrain2D(
					    std::vector<double_array>* strain,
					    const double_array& basisDeriv,
					    const double* disp,
					    const int numBasis)
{ // calcTotalStrain2D
  assert(0 != strain);
  assert(0 != disp);
  
  const int dimension = 2;
  const int numQuadPts = strain->size();

  assert(basisDeriv.size() == numQuadPts*numBasis*dimension);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    assert(3 == (*strain)[iQuad].size());
    for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
      (*strain)[iQuad][0] += 
	basisDeriv[iQ+iBasis  ] * disp[iBasis  ];
      (*strain)[iQuad][1] += 
	basisDeriv[iQ+iBasis+1] * disp[iBasis+1];
      (*strain)[iQuad][2] += 
	0.5 * (basisDeriv[iQ+iBasis+1] * disp[iBasis  ] +
	       basisDeriv[iQ+iBasis  ] * disp[iBasis+1]);
    } // for
  } // for
} // calcTotalStrain2D

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticity::calcTotalStrain3D(
					    std::vector<double_array>* strain,
					    const double_array& basisDeriv,
					    const double* disp,
					    const int numBasis)
{ // calcTotalStrain3D
  assert(0 != strain);
  assert(0 != disp);

  const int dimension = 3;
  const int numQuadPts = strain->size();

  assert(basisDeriv.size() == numQuadPts*numBasis*dimension);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    assert(6 == (*strain)[iQuad].size());
    for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
      (*strain)[iQuad][0] += 
	basisDeriv[iQ+iBasis  ] * disp[iBasis  ];
      (*strain)[iQuad][1] += 
	basisDeriv[iQ+iBasis+1] * disp[iBasis+1];
      (*strain)[iQuad][2] += 
	basisDeriv[iQ+iBasis+2] * disp[iBasis+2];
      (*strain)[iQuad][3] += 
	0.5 * (basisDeriv[iQ+iBasis+1] * disp[iBasis  ] +
	       basisDeriv[iQ+iBasis  ] * disp[iBasis+1]);
      (*strain)[iQuad][4] += 
	0.5 * (basisDeriv[iQ+iBasis+2] * disp[iBasis+1] +
	       basisDeriv[iQ+iBasis+1] * disp[iBasis+2]);
      (*strain)[iQuad][5] += 
	0.5 * (basisDeriv[iQ+iBasis+2] * disp[iBasis  ] +
	       basisDeriv[iQ+iBasis  ] * disp[iBasis+2]);
    } // for
  } // for
} // calcTotalStrain3D


// End of file 
