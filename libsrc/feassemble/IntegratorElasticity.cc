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
pylith::feassemble::IntegratorElasticity::calcTotalStrain(
					  std::vector<double_array>* strain,
					  const double_array& basisDeriv,
					  const double* disp,
					  const int dimension,
					  const int numBasis)
{ // calcTotalStrain
  assert(0 != strain);

  const int numQuadPts = strain->size();

  assert(basisDeriv.size() == numQuadPts*numBasis*dimension);

  if (3 == dimension) {
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      assert(6 == (*strain)[iQuad].size());
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	strain[iQuad][0] += 
	  basisDeriv[iQ+iBasis  ] * disp[iBasis  ];
	strain[iQuad][1] += 
	  basisDeriv[iQ+iBasis+1] * disp[iBasis+1];
	strain[iQuad][2] += 
	  basisDeriv[iQ+iBasis+2] * disp[iBasis+2];
	strain[iQuad][3] += 
	  0.5 * (basisDeriv[iQ+iBasis+1] * disp[iBasis  ] +
		 basisDeriv[iQ+iBasis  ] * disp[iBasis+1]);
	strain[iQuad][4] += 
	  0.5 * (basisDeriv[iQ+iBasis+2] * disp[iBasis+1] +
		 basisDeriv[iQ+iBasis+1] * disp[iBasis+2]);
	strain[iQuad][5] += 
	  0.5 * (basisDeriv[iQ+iBasis+2] * disp[iBasis  ] +
		 basisDeriv[iQ+iBasis  ] * disp[iBasis+2]);
      } // for
    } // for
  } else if (2 == dimension) {
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      assert(3 == (*strain)[iQuad].size());
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
	strain[iQuad][0] += 
	  basisDeriv[iQ+iBasis  ] * disp[iBasis  ];
	strain[iQuad][1] += 
	  basisDeriv[iQ+iBasis+1] * disp[iBasis+1];
	strain[iQuad][2] += 
	  0.5 * (basisDeriv[iQ+iBasis+1] * disp[iBasis  ] +
		 basisDeriv[iQ+iBasis  ] * disp[iBasis+1]);
      } // for
    } // for
  } else if (1 == dimension) {
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      assert(1 == (*strain)[iQuad].size());
      for (int iBasis=0; iBasis < numBasis; ++iBasis)
	(*strain)[iQuad][0] += 
	  basisDeriv[iQuad*numBasis+iBasis] * disp[iBasis];
    } // for
  } else {
    throw std::logic_error("Dimension out of range.");
  } // else
} // calcTotalStrain


// End of file 
