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

#include "pylith/utils/array.hh" // USES double_array
#include "pylith/utils/constdefs.h" // USES MAXDOUBLE
#include "pylith/utils/sievetypes.hh" // USES Mesh

#include <string.h> // USES memcpy()
#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticMaterial::ElasticMaterial(const int tensorSize,
						    const int numElasticConsts,
						    const char** dbValues,
						    const char** initialStateDBValues,
						    const int numDBValues,
						    const PropMetaData* properties,
						    const int numProperties) :
  Material(tensorSize, dbValues, initialStateDBValues, numDBValues, properties, numProperties),
  _numQuadPts(0),
  _numElasticConsts(numElasticConsts)
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
  const int totalPropsQuadPt = _totalPropsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*totalPropsQuadPt);
  assert(_density.size() == numQuadPts);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcDensity(&_density[iQuad], 
		 &_propertiesCell[iQuad*totalPropsQuadPt], totalPropsQuadPt);

  return _density;
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor for cell at quadrature points.
const pylith::double_array&
pylith::materials::ElasticMaterial::calcStress(const double_array& totalStrain,
					       const bool computeStateVars)
{ // calcStress
  const int numQuadPts = _numQuadPts;
  const int totalPropsQuadPt = _totalPropsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*totalPropsQuadPt);
  assert(totalStrain.size() == numQuadPts*_tensorSize);
  assert(_stress.size() == numQuadPts*_tensorSize);
  assert(_initialStateCell.size() == numQuadPts*_initialStateSize);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcStress(&_stress[iQuad*_tensorSize], _tensorSize,
		&_propertiesCell[iQuad*totalPropsQuadPt], totalPropsQuadPt,
		&totalStrain[iQuad*_tensorSize], _tensorSize, 
		&_initialStateCell[iQuad*_initialStateSize],
		_initialStateSize, computeStateVars);

  return _stress;
} // calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix for cell at quadrature points.
const pylith::double_array&
pylith::materials::ElasticMaterial::calcDerivElastic(
					       const double_array& totalStrain)
{ // calcDerivElastic
  const int numQuadPts = _numQuadPts;
  const int totalPropsQuadPt = _totalPropsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*totalPropsQuadPt);
  assert(totalStrain.size() == numQuadPts*_tensorSize);
  assert(_elasticConsts.size() == numQuadPts*_numElasticConsts);
  assert(_initialStateCell.size() == numQuadPts*_initialStateSize);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcElasticConsts(&_elasticConsts[iQuad*_numElasticConsts], 
		       _numElasticConsts,
		       &_propertiesCell[iQuad*totalPropsQuadPt], 
		       totalPropsQuadPt, 
		       &totalStrain[iQuad*_tensorSize], _tensorSize,
		       &_initialStateCell[iQuad*_initialStateSize],
		       _initialStateSize);

  return _elasticConsts;
} // calcDerivElastic

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::ElasticMaterial::stableTimeStepImplicit(void) const
{ // stableTimeStepImplicit
  const int numQuadPts = _numQuadPts;
  const int totalPropsQuadPt = _totalPropsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*totalPropsQuadPt);

  double dtStable = pylith::PYLITH_MAXDOUBLE;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const double dt = 
      _stableTimeStepImplicit(&_propertiesCell[iQuad*totalPropsQuadPt],
			      totalPropsQuadPt);
    if (dt < dtStable)
      dtStable = dt;
  } // for

  return dtStable;
} // stableTimeStepImplicit

// ----------------------------------------------------------------------
// Get cell's property information from material's sections.
void
pylith::materials::ElasticMaterial::getPropertiesCell(const Mesh::point_type& cell,
						      const int numQuadPts)
{ // getPropertiesCell
  if (_numQuadPts != numQuadPts) {
    const int totalPropsQuadPt = _totalPropsQuadPt;

    _numQuadPts = numQuadPts;

    _propertiesCell.resize(numQuadPts*totalPropsQuadPt);
    _density.resize(numQuadPts*1);
    _stress.resize(numQuadPts*_tensorSize);
    _elasticConsts.resize(numQuadPts*_numElasticConsts);
    _initialStateCell.resize(numQuadPts*_tensorSize);
  } // if

  _getProperties(cell);
} // getPropertiesCell

// ----------------------------------------------------------------------
// Update properties (state variables for next time step).
void
pylith::materials::ElasticMaterial::updateProperties(
					      const double_array& totalStrain,
					      const Mesh::point_type& cell)
{ // updateProperties
  const int numQuadPts = _numQuadPts;
  const int totalPropsQuadPt = _totalPropsQuadPt;
  getPropertiesCell(cell, numQuadPts);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _updateProperties(&_propertiesCell[iQuad*totalPropsQuadPt], 
		      totalPropsQuadPt,
		      &totalStrain[iQuad*_tensorSize], _tensorSize,
		      &_initialStateCell[iQuad*_initialStateSize],
		      _initialStateSize);
  
  _properties->updatePoint(cell, &_propertiesCell[0]);
} // updateProperties

// ----------------------------------------------------------------------
// Get properties for cell.
void
pylith::materials::ElasticMaterial::_getProperties(const Mesh::point_type& cell)
{ // _getProperties
  assert(!_properties.isNull());

  const int numQuadPts = _numQuadPts;
  const int totalPropsQuadPt = _totalPropsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*totalPropsQuadPt);  
  assert(_properties->getFiberDimension(cell) == numQuadPts*totalPropsQuadPt);
  assert(_initialStateCell.size() == numQuadPts*_initialStateSize);
  const real_section_type::value_type* parameterCell =
    _properties->restrictPoint(cell);
  memcpy(&_propertiesCell[0], parameterCell, 
		      numQuadPts*totalPropsQuadPt*sizeof(double));
  if (_initialState.isNull())
    _initialStateCell = 0.0;
  else {
    assert(_initialState->getFiberDimension(cell) == numQuadPts*_initialStateSize);
    const real_section_type::value_type* initialStateValuesCell =
      _initialState->restrictPoint(cell);
    memcpy(&_initialStateCell[0], initialStateValuesCell, 
	   numQuadPts*_initialStateSize*sizeof(double));
  } // else

} // _getProperties

// ----------------------------------------------------------------------
// Update properties (for next time step).
void
pylith::materials::ElasticMaterial::_updateProperties(double* const properties,
						      const int totalPropsQuadPt,
						      const double* totalStrain,
						      const int strainSize,
						      const double* initialState,
						      const int initialStateSize)
{ // _updateProperties
} // _updateProperties


// End of file 
