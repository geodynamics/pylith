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

#include "pylith/utils/sievetypes.hh" // USES Mesh

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticMaterial::ElasticMaterial(const int* numParamValues,
						    const int size) :
  Material(numParamValues, size),
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
const pylith::double_array&
pylith::materials::ElasticMaterial::calcDensity(void)
{ // calcDensity
  const int numQuadPts = _numQuadPts;
  const int numParamsQuadPt = _numParamsQuadPt;
  assert(_paramsCell.size() == numQuadPts*numParamsQuadPt);
  assert(_density.size() == numQuadPts);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcDensity(&_density[iQuad], 
		 &_paramsCell[iQuad*numParamsQuadPt], numParamsQuadPt);

  return _density;
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor for cell at quadrature points.
const pylith::double_array&
pylith::materials::ElasticMaterial::calcStress(const double_array& totalStrain)
{ // calcStress
  const int numQuadPts = _numQuadPts;
  const int numParamsQuadPt = _numParamsQuadPt;
  const int tensorSize = _tensorSize();
  assert(_paramsCell.size() == numQuadPts*numParamsQuadPt);
  assert(totalStrain.size() == numQuadPts*tensorSize);
  assert(_stress.size() == numQuadPts*tensorSize);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcStress(&_stress[iQuad*tensorSize], tensorSize,
		&_paramsCell[iQuad*numParamsQuadPt], numParamsQuadPt,
		&totalStrain[iQuad*tensorSize], tensorSize);

  return _stress;
} // calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix for cell at quadrature points.
const pylith::double_array&
pylith::materials::ElasticMaterial::calcDerivElastic(
					       const double_array& totalStrain)
{ // calcDerivElastic
  const int numQuadPts = _numQuadPts;
  const int numParamsQuadPt = _numParamsQuadPt;
  const int tensorSize = _tensorSize();
  const int numElasticConsts = _numElasticConsts();
  assert(_paramsCell.size() == numQuadPts*numParamsQuadPt);
  assert(totalStrain.size() == numQuadPts*tensorSize);
  assert(_elasticConsts.size() == numQuadPts*numElasticConsts);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcElasticConsts(&_elasticConsts[iQuad*numElasticConsts], numElasticConsts,
		       &_paramsCell[iQuad*numParamsQuadPt], numParamsQuadPt, 
		       &totalStrain[iQuad*tensorSize], tensorSize);

  return _elasticConsts;
} // calcDerivElastic

// ----------------------------------------------------------------------
// Get cell's state variable information from material's sections.
void
pylith::materials::ElasticMaterial::getStateVarsCell(
						   const Mesh::point_type& cell,
						   const int numQuadPts)
{ // getStateVarsCell
  if (_numQuadPts != numQuadPts) {
    const int numParamsQuadPt = _numParamsQuadPt;
    const int tensorSize = _tensorSize();
    const int numElasticConsts = _numElasticConsts();

    _numQuadPts = numQuadPts;

    _paramsCell.resize(numQuadPts*numParamsQuadPt);
    _density.resize(numQuadPts*1);
    _stress.resize(numQuadPts*tensorSize);
    _elasticConsts.resize(numQuadPts*numElasticConsts);
  } // if

  _getParameters(cell);
} // getStateVarsCell

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::materials::ElasticMaterial::updateState(const double_array& totalStrain,
						const Mesh::point_type& cell)
{ // updateState
  const int numQuadPts = _numQuadPts;
  const int numParamsQuadPt = _numParamsQuadPt;
  const int tensorSize = _tensorSize();
  getStateVarsCell(cell, numQuadPts);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _updateState(&_paramsCell[iQuad*numParamsQuadPt], numParamsQuadPt,
		 &totalStrain[iQuad*tensorSize], tensorSize);
  
  _parameters->updatePoint(cell, &_paramsCell[0]);
} // updateState

// ----------------------------------------------------------------------
// Get parameters for cell.
void
pylith::materials::ElasticMaterial::_getParameters(const Mesh::point_type& cell)
{ // _getParameters
  assert(!_parameters.isNull());

  const int numQuadPts = _numQuadPts;
  const int numParamsQuadPt = _numParamsQuadPt;
  assert(_paramsCell.size() == numQuadPts*numParamsQuadPt);  
  assert(_parameters->getFiberDimension(cell) == numQuadPts*numParamsQuadPt);
  const real_section_type::value_type* parameterCell =
    _parameters->restrictPoint(cell);
  memcpy(&_paramsCell[0], parameterCell, 
		      numQuadPts*numParamsQuadPt*sizeof(double));
} // _getParameters

// ----------------------------------------------------------------------
// Update parameters (for next time step).
void
pylith::materials::ElasticMaterial::_updateState(double* const parameters,
						 const int numParamsQuadPt,
						 const double* totalStrain,
						 const int strainSize)
{ // _updateState
} // _updateState


// End of file 
