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

#include "ElasticMaterial.hh" // implementation of object methods

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
const int pylith::materials::ElasticMaterial::NUMELASTCONSTS = 21;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticMaterial::ElasticMaterial(void) :
  _pDB(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticMaterial::~ElasticMaterial(void)
{ // destructor
  // Do not own database, so just set pointer to null.
  _pDB = 0;
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::materials::ElasticMaterial::ElasticMaterial(const ElasticMaterial& m) :
  _pDB(m._pDB)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Set database for material parameters.
void
pylith::materials::ElasticMaterial::parametersDB(
                            spatialdata::spatialdb::SpatialDB* pDB)
{ // parametersDB
  _pDB = pDB;
} // parametersDB

// ----------------------------------------------------------------------
// Get intertia at points.
void
pylith::materials::ElasticMaterial::inertia(
			     double** ppInertia,
			     int* pSize,
			     const double* pPts,
			     const int npts,
			     const spatialdata::geocoords::CoordSys* pCS,
			     const void* pState) const
{ // inertia
  double* pParams = 0;
  int size = 0;
  _getParameters(&pParams, &size, pPts, npts, pCS);
  calcInertia(ppInertia, pSize, pParams, size, pState);
  delete[] pParams; pParams = 0;
} // inertia

// ----------------------------------------------------------------------
/* Get elasticity constants at quadrature points. */
void
pylith::materials::ElasticMaterial::elasticityConsts(
			 double** ppElasticityC,
			 int* pSize,
			 const double* pPts,
			 const int npts,
			 const spatialdata::geocoords::CoordSys* pCS,
			 const void* pState) const
{ // elasticityConsts
  double* pParams = 0;
  int size = 0;
  _getParameters(&pParams, &size, pPts, npts, pCS);
  calcElasticityConsts(ppElasticityC, pSize, pParams, size, pState);
  delete[] pParams; pParams = 0;
} // elasticityConsts

// ----------------------------------------------------------------------
// Get material parameters at a set of points.
void
pylith::materials::ElasticMaterial::_getParameters(double** ppParams,
						   int* pSize,
						   const double* pPts,
						   const int npts,
			 const spatialdata::geocoords::CoordSys* pCS) const
{ // _getParameters
  assert(0 != ppParams);
  assert(0 != pSize);
  assert( (0 == pPts && 0 == npts) ||
	  (0 != pPts && 0 < npts) );
  assert(0 != _pDB);

  const int size = npts*numParams();
  delete[] *ppParams; *ppParams = (size > 0) ? new double[size] : 0;
  *pSize = size;

  const int numCoords = 3;
  for (int ipt=0; ipt < npts; ++ipt) {
    double* pVals = 0;
    int nvals = numParams();
    const double* pCoords = pPts + numCoords*ipt;
    _pDB->query(&pVals, nvals, pCoords[0], pCoords[1], pCoords[2], pCS);
    for (int ival=0; ival < nvals; ++ival)
      (*ppParams)[ipt*nvals+ival] = pVals[ival];
  } // for
} // _getParameters

// version
// $Id$

// End of file 
