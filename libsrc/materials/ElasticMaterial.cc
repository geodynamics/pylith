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
#include "pylith/utils/array.hh" // USES double_array

#include <petscmesh.h> // USES Mesh

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticMaterial::ElasticMaterial(void) :
  _numQuadPts(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticMaterial::~ElasticMaterial(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::materials::ElasticMaterial::ElasticMaterial(const ElasticMaterial& m) :
  Material(m),
  _numQuadPts(0)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute density for cell at quadrature points.
const std::vector<pylith::double_array>&
pylith::materials::ElasticMaterial::calcDensity(void)
{ // calcDensity
  const int numQuadPts = _numQuadPts;
  assert(_paramsCell.size() == numQuadPts);
  assert(_density.size() == numQuadPts);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcDensity(&_density[iQuad], _paramsCell[iQuad]);

  return _density;
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor for cell at quadrature points.
const std::vector<pylith::double_array>&
pylith::materials::ElasticMaterial::calcStress(
			         const std::vector<double_array>& totalStrain)
{ // calcStress
  const int numQuadPts = _numQuadPts;
  assert(_paramsCell.size() == numQuadPts);
  assert(totalStrain.size() == numQuadPts);
  assert(_stress.size() == numQuadPts);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcStress(&_stress[iQuad], _paramsCell[iQuad], totalStrain[iQuad]);

  return _stress;
} // calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix for cell at quadrature points.
const std::vector<pylith::double_array>&
pylith::materials::ElasticMaterial::calcDerivElastic(
				const std::vector<double_array>& totalStrain)
{ // calcDerivElastic
  const int numQuadPts = _numQuadPts;
  assert(_paramsCell.size() == numQuadPts);
  assert(totalStrain.size() == numQuadPts);
  assert(_elasticConsts.size() == numQuadPts);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcElasticConsts(&_elasticConsts[iQuad], _paramsCell[iQuad], 
		       totalStrain[iQuad]);

  return _elasticConsts;
} // calcDerivElastic

// ----------------------------------------------------------------------
// Initialize arrays holding cell data.
void
pylith::materials::ElasticMaterial::initCellData(const Mesh::point_type& cell,
						 const int numQuadPts)
{ // initCellData
  if (_numQuadPts != numQuadPts) {
    _numQuadPts = numQuadPts;

    _paramsCell.resize(numQuadPts);
    _density.resize(numQuadPts);
    _stress.resize(numQuadPts);
    _elasticConsts.resize(numQuadPts);
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      _paramsCell[iQuad].resize(_numParameters());
      _density[iQuad].resize(1);
      _stress[iQuad].resize(_tensorSize());
      _elasticConsts[iQuad].resize(_numElasticConsts());
    } // for
  } // if

  _getParameters(cell);
} // initCellData


// ----------------------------------------------------------------------
// Get parameters for cell.
void
pylith::materials::ElasticMaterial::_getParameters(const Mesh::point_type& cell)
{ // _getParameters
  const int numQuadPts = _numQuadPts;
  assert(_paramsCell.size() == numQuadPts);
  
  const int numParams = _numParameters();
  const char** paramNames = _parameterNames();
  
  for (int iParam=0; iParam < numParams; ++iParam) {
    const ALE::Obj<real_section_type> parameter = 
      _parameters->getReal(paramNames[iParam]);
    
    assert(parameter->getFiberDimension(cell) == numQuadPts);
    const real_section_type::value_type* parameterCell =
      parameter->restrictPoint(cell);
    for (int iQuadPt=0; iQuadPt < numQuadPts; ++iQuadPt)
      _paramsCell[iQuadPt][iParam] = parameterCell[iQuadPt];
  } // for
} // _getParameters


// End of file 
