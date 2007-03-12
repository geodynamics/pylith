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

#include "pylith/feassemble/ParameterManager.hh" // USES ParameterManager

#include <petscmesh.h> // USES Mesh

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticMaterial::ElasticMaterial(void) :
  _density(0),
  _elasticConsts(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticMaterial::~ElasticMaterial(void)
{ // destructor
  delete[] _density; _density = 0;
  delete[] _elasticConsts; _elasticConsts = 0;
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::materials::ElasticMaterial::ElasticMaterial(const ElasticMaterial& m) :
  Material(m),
  _density(0),
  _elasticConsts(0)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute physical properties of cell at quadrature points.
void
pylith::materials::ElasticMaterial::calcProperties(
				     const topology_type::point_type& cell,
				     const topology_type::patch_type& patch,
				     const int numQuadPts)
{ // calcProperties
  _initCellData(numQuadPts);

  const int numParams = _numParameters();
  const char** paramNames = _parameterNames();

  int size = numParams * numQuadPts;
  double* paramsCell = (size > 0) ? new double[size] : 0;
  for (int iParam=0; iParam < numParams; ++iParam) {
    const ALE::Obj<real_section_type> parameter = 
      _parameters->getReal(paramNames[iParam]);
    assert(numQuadPts == parameter->getFiberDimension(patch, cell));
    const real_section_type::value_type* parameterCell =
      parameter->restrict(patch, cell);
    for (int iQuadPt=0; iQuadPt < numQuadPts; ++iQuadPt)
      paramsCell[iQuadPt*numParams+iParam] = parameterCell[iQuadPt];
  } // for
  _calcDensity(paramsCell, numParams, numQuadPts);
  _calcElasticConsts(paramsCell, numParams, numQuadPts);
} // calcProperties

// ----------------------------------------------------------------------
// Initialize arrays holding cell data.
void
pylith::materials::ElasticMaterial::_initCellData(const int numQuadPts)
{ // _initCellData
  int size = numQuadPts;
  delete[] _density; _density = (size > 0) ? new double[size] : 0;
  size = numQuadPts * numElasticConsts();
  delete[] _elasticConsts; _elasticConsts = (size > 0) ? new double[size] : 0;
} // _initCellData


// End of file 
