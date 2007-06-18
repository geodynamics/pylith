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

#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
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
// Get cell's state variable information from material's sections.
void
pylith::materials::ElasticMaterial::getStateVarsCell(const Mesh::point_type& cell,
						 const int numQuadPts)
{ // getStateVarsCell
  if (_numQuadPts != numQuadPts) {
    _numQuadPts = numQuadPts;

    _paramsCell.resize(numQuadPts);
    _density.resize(numQuadPts);
    _stress.resize(numQuadPts);
    _elasticConsts.resize(numQuadPts);
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      _density[iQuad].resize(1);
      _stress[iQuad].resize(_tensorSize());
      _elasticConsts[iQuad].resize(_numElasticConsts());
      int_array numParamValues;
      _numParamValues(&numParamValues);
      const int numParams = numParamValues.size();
      _paramsCell[iQuad].resize(numParams);
      for (int iParam=0; iParam < numParams; ++iParam)
	_paramsCell[iQuad][iParam].resize(numParamValues[iParam]);
    } // for
  } // if

  _getParameters(cell);
} // getStateVarsCell

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::materials::ElasticMaterial::updateState(
				const std::vector<double_array>& totalStrain,
				const Mesh::point_type& cell)
{ // updateState
  const int numQuadPts = _numQuadPts;
  getStateVarsCell(cell, numQuadPts);

  assert(_paramsCell.size() == numQuadPts);
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _updateState(&_paramsCell[iQuad], totalStrain[iQuad]);

  int_array numParamValues;
  _numParamValues(&numParamValues);
  const int numParams = numParamValues.size();
  const char** paramNames = _parameterNames();
  
  for (int iParam=0; iParam < numParams; ++iParam) {
    const ALE::Obj<real_section_type> parameter = 
      _parameters->getReal(paramNames[iParam]);
    assert(!parameter.isNull());

    const int numValues = numParamValues[iParam];
    double_array parameterCell(numQuadPts*numValues);
    
    assert(parameter->getFiberDimension(cell) == numQuadPts*numValues);
    for (int iQuadPt=0; iQuadPt < numQuadPts; ++iQuadPt)
      for (int iValue=0; iValue < numValues; ++iValue)
	parameterCell[iQuadPt*numValues+iValue] = 
	  _paramsCell[iQuadPt][iParam][iValue];
    parameter->updatePoint(cell, &parameterCell[0]);
  } // for
} // updateState

// ----------------------------------------------------------------------
// Get parameters for cell.
void
pylith::materials::ElasticMaterial::_getParameters(const Mesh::point_type& cell)
{ // _getParameters
  assert(0 != _parameters);

  const int numQuadPts = _numQuadPts;
  assert(_paramsCell.size() == numQuadPts);
  
  int_array numParamValues;
  _numParamValues(&numParamValues);
  const int numParams = numParamValues.size();
  const char** paramNames = _parameterNames();
  
  for (int iParam=0; iParam < numParams; ++iParam) {
    const ALE::Obj<real_section_type> parameter = 
      _parameters->getReal(paramNames[iParam]);
    assert(!parameter.isNull());

    const int numValues = numParamValues[iParam];
    
    assert(parameter->getFiberDimension(cell) == numQuadPts*numValues);
    const real_section_type::value_type* parameterCell =
      parameter->restrictPoint(cell);
    for (int iQuadPt=0; iQuadPt < numQuadPts; ++iQuadPt)
      for (int iValue=0; iValue < numValues; ++iValue)
	_paramsCell[iQuadPt][iParam][iValue] = 
	  parameterCell[iQuadPt*numValues+iValue];
  } // for
} // _getParameters

// ----------------------------------------------------------------------
// Update parameters (for next time step).
void
pylith::materials::ElasticMaterial::_updateState(
				 std::vector<double_array>* const parameters,
				 const double_array& totalStrain)
{ // _updateState
} // _updateState


// End of file 
