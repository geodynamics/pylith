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
  _stress(0),
  _elasticConsts(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticMaterial::~ElasticMaterial(void)
{ // destructor
  delete[] _density; _density = 0;
  delete[] _stress; _stress = 0;
  delete[] _elasticConsts; _elasticConsts = 0;
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::materials::ElasticMaterial::ElasticMaterial(const ElasticMaterial& m) :
  Material(m),
  _density(0),
  _stress(0),
  _elasticConsts(0)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute density for cell at quadrature points.
const double*
pylith::materials::ElasticMaterial::calcDensity(
				     const Mesh::point_type& cell,
				     const int numQuadPts)
{ // calcDensity
  assert(0 != _density);
  int size = numQuadPts;
  memset(_density, 0, size*sizeof(double));

  double* paramsCell = 0;
  _getParameters(&paramsCell, cell, numQuadPts);
  const int numParams = _numParameters();
  _calcDensity(paramsCell, numParams, numQuadPts);
  delete[] paramsCell; paramsCell = 0;

  return _density;
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor for cell at quadrature points.
const double*
pylith::materials::ElasticMaterial::calcStress(
				     const Mesh::point_type& cell,
				     const double* totalStrain,
				     const int numQuadPts,
				     const int spaceDim)
{ // calcStress
  assert(0 != totalStrain);
  assert(0 != _stress);

  int size = numQuadPts * stressSize();
  memset(_stress, 0, size*sizeof(double));

  double* paramsCell = 0;
  _getParameters(&paramsCell, cell, numQuadPts);
  const int numParams = _numParameters();
  _calcStress(paramsCell, numParams, totalStrain, numQuadPts, spaceDim);
  delete[] paramsCell; paramsCell = 0;

  return _stress;
} // calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix for cell at quadrature points.
const double*
pylith::materials::ElasticMaterial::calcDerivElastic(
				     const Mesh::point_type& cell,
				     const double* totalStrain,
				     const int numQuadPts,
				     const int spaceDim)
{ // calcDerivElastic
  assert(0 != totalStrain);
  assert(0 != _elasticConsts);

  int size = numQuadPts * numElasticConsts();
  memset(_elasticConsts, 0, size*sizeof(double));

  double* paramsCell = 0;
  _getParameters(&paramsCell, cell, numQuadPts);
  const int numParams = _numParameters();
  _calcElasticConsts(paramsCell, numParams, totalStrain, numQuadPts, spaceDim);
  delete[] paramsCell; paramsCell = 0;

  return _elasticConsts;
} // calcDerivElastic

// ----------------------------------------------------------------------
// Initialize arrays holding cell data.
void
pylith::materials::ElasticMaterial::_initCellData(const int numQuadPts)
{ // _initCellData
  int size = numQuadPts;
  delete[] _density; _density = (size > 0) ? new double[size] : 0;

  size = numQuadPts * stressSize();
  delete[] _stress; _stress = (size > 0) ? new double[size] : 0;

  size = numQuadPts * numElasticConsts();
  delete[] _elasticConsts; _elasticConsts = (size > 0) ? new double[size] : 0;
} // _initCellData


// ----------------------------------------------------------------------
// Get parameters for cell.
void
pylith::materials::ElasticMaterial::_getParameters(double** paramsCell,
				       const Mesh::point_type& cell,
				       const int numQuadPts)
{ // _getParameters
  const int numParams = _numParameters();
  const char** paramNames = _parameterNames();

  const int size = numParams * numQuadPts;
  delete[] *paramsCell; *paramsCell = (size > 0) ? new double[size] : 0;
  for (int iParam=0; iParam < numParams; ++iParam) {
    const ALE::Obj<real_section_type> parameter = 
      _parameters->getReal(paramNames[iParam]);
    
    assert(numQuadPts == parameter->getFiberDimension(cell));
    const real_section_type::value_type* parameterCell =
      parameter->restrictPoint(cell);
    for (int iQuadPt=0; iQuadPt < numQuadPts; ++iQuadPt)
      (*paramsCell)[iQuadPt*numParams+iParam] = parameterCell[iQuadPt];
  } // for
} // _getParameters


// End of file 
