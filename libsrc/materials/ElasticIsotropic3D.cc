// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "ElasticMaterial.hh" // ISA ElasticMaterial
#include "ElasticIsotropic3D.hh" // implementation of object methods

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    class _ElasticIsotropic3D;
  } // materials
} // pylith
class pylith::materials::_ElasticIsotropic3D {
public :
  static const char* NAMESPARAMS[];
  static const int NUMPARAMS;

  /// Indices of paramters in queries
  enum parameters {
    IVP = 0, ///< Index for Vp
    IVS = 1, ///< Index for Vs
    IDENSITY = 2 ///< Index for density
  };
}; // _ElasticIsotropic3D
const char* pylith::materials::_ElasticIsotropic3D::NAMESPARAMS[] =
  { "Vp", "Vs", "Density" };
const int pylith::materials::_ElasticIsotropic3D::NUMPARAMS = 3;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticIsotropic3D::ElasticIsotropic3D(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticIsotropic3D::~ElasticIsotropic3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::materials::ElasticIsotropic3D::ElasticIsotropic3D(
						const ElasticIsotropic3D& m) :
  ElasticMaterial(m)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Get names of parameters for material.
const char**
pylith::materials::ElasticIsotropic3D::namesParams(void) const
{ // namesParams
  return _ElasticIsotropic3D::NAMESPARAMS;
} // namesParams
  
// ----------------------------------------------------------------------
// Get number of parameters for material.
int
pylith::materials::ElasticIsotropic3D::numParams(void) const
{ // numParams
  return _ElasticIsotropic3D::NUMPARAMS;
} // numParams
  
// ----------------------------------------------------------------------
// Compute inertia at points using material parameters.
void
pylith::materials::ElasticIsotropic3D::calcInertia(double** ppInertia,
						   int* pSize,
						   const double* pParams,
						   const int npts,
						   const void* pState) const
{ // calcInertia
  assert(0 != ppInertia);
  assert(0 != pSize);
  assert( (0 == pParams && 0 == npts) ||
	  (0 != pParams && 0 < npts) );

  delete[] *ppInertia; *ppInertia = (npts > 0) ? new double[npts] : 0;
  *pSize = npts;

  for (int ipt=0; ipt < npts; ++ipt) {
    const double* pParamsPt = pParams + ipt*_ElasticIsotropic3D::NUMPARAMS;
    (*ppInertia)[ipt] = pParamsPt[_ElasticIsotropic3D::IDENSITY];
  } // for
} // calcInertia

// ----------------------------------------------------------------------
// Compute elasticity constants at points using material parameters.
void
pylith::materials::ElasticIsotropic3D::calcElasticityConsts(
						    double** ppElasticityC,
						    int* pSize,
						    const double* pParams,
						    const int npts,
						    const void* pState) const
{ // calcElasticityConsts
  assert(0 != ppElasticityC);
  assert(0 != pSize);
  assert( (0 == pParams && 0 == npts) ||
	  (0 != pParams && 0 < npts) );

  delete[] *ppElasticityC; 
  *ppElasticityC = (npts > 0) ? new double[NUMELASTCONSTS*npts] : 0;
  *pSize = NUMELASTCONSTS*npts;

  for (int ipt=0; ipt < npts; ++ipt) {
    const double* pParamsPt = pParams + ipt*_ElasticIsotropic3D::NUMPARAMS;
    const double density = pParamsPt[_ElasticIsotropic3D::IDENSITY];
    const double vp = pParamsPt[_ElasticIsotropic3D::IVP];
    const double vs = pParamsPt[_ElasticIsotropic3D::IVS];

    const double mu = density*vs*vs;
    const double lambda = density*vp*vp - 2.0*mu;
    
    const int iElasticityC = NUMELASTCONSTS*ipt;
    (*ppElasticityC)[iElasticityC   ] = lambda + 2.0*mu; // C1111
    (*ppElasticityC)[iElasticityC+ 1] = lambda; // C1122
    (*ppElasticityC)[iElasticityC+ 2] = lambda; // C1133
    (*ppElasticityC)[iElasticityC+ 3] = 0; // C1112
    (*ppElasticityC)[iElasticityC+ 4] = 0; // C1123
    (*ppElasticityC)[iElasticityC+ 5] = 0; // C1113
    (*ppElasticityC)[iElasticityC+ 6] = lambda + 2.0*mu; // C2222
    (*ppElasticityC)[iElasticityC+ 7] = lambda; // C2233
    (*ppElasticityC)[iElasticityC+ 8] = 0; // C2212
    (*ppElasticityC)[iElasticityC+ 9] = 0; // C2223
    (*ppElasticityC)[iElasticityC+10] = 0; // C2212
    (*ppElasticityC)[iElasticityC+11] = lambda + 2.0*mu; // C3333
    (*ppElasticityC)[iElasticityC+12] = 0; // C3312
    (*ppElasticityC)[iElasticityC+13] = 0; // C3323
    (*ppElasticityC)[iElasticityC+14] = 0; // C3312
    (*ppElasticityC)[iElasticityC+15] = mu; // C1212
    (*ppElasticityC)[iElasticityC+16] = 0; // C1223
    (*ppElasticityC)[iElasticityC+17] = 0; // C1213
    (*ppElasticityC)[iElasticityC+18] = mu; // C2323
    (*ppElasticityC)[iElasticityC+19] = 0; // C2313
    (*ppElasticityC)[iElasticityC+20] = mu; // C1313
  } // for
} // elasticityConsts

// version
// $Id$

// End of file 
