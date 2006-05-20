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

#include "spatialdata/spatialdb/SpatialDB.h" // USES SpatialDB

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
  delete _pDB; _pDB = 0;
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::materials::ElasticMaterial::ElasticMaterial(const ElasticMaterial& m) :
  _pDB(0)
{ // copy constructor
#if 0
  if (0 != m._pDB)
    _pDB = _pDB->clone();
#endif
} // copy constructor

// ----------------------------------------------------------------------
// Set database for material parameters.
void
pylith::materials::ElasticMaterial::parametersDB(
                            const spatialdata::spatialdb::SpatialDB* pDB)
{ // parametersDB
  delete _pDB;
#if 0
  _pDB = (0 != pDB) ? pDB->clone : 0;
#endif
} // parametersDB

// ----------------------------------------------------------------------
// Initialize material.
void
pylith::materials::ElasticMaterial::initialize(void)
{ // initialize
#if 0
  if (0 != _pDB)
    _pDB->initialize();
#endif
} // initialize

// ----------------------------------------------------------------------
// Open database.
void
pylith::materials::ElasticMaterial::openDB(void)
{ // openDB
  if (0 != _pDB)
    _pDB->open();
} // openDB

// ----------------------------------------------------------------------
// Close database.
void
pylith::materials::ElasticMaterial::closeDB(void)
{ // closeDB
  if (0 != _pDB)
    _pDB->close();
} // closeDB

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
