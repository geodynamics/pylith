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
#include "pylith/utils/sievetypes.hh" // USES Mesh

#include <string.h> // USES memcpy()
#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticMaterial::ElasticMaterial(const int tensorSize,
						    const int numElasticConsts,
						    const char** dbValues,
						    const int numDBValues,
						    const PropMetaData* properties,
						    const int numProperties) :
  Material(dbValues, numDBValues, properties, numProperties),
  _numQuadPts(0),
  _tensorSize(tensorSize),
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

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcStress(&_stress[iQuad*_tensorSize], _tensorSize,
		&_propertiesCell[iQuad*totalPropsQuadPt], totalPropsQuadPt,
		&totalStrain[iQuad*_tensorSize], _tensorSize, 
		computeStateVars);

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

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcElasticConsts(&_elasticConsts[iQuad*_numElasticConsts], 
		       _numElasticConsts,
		       &_propertiesCell[iQuad*totalPropsQuadPt], 
		       totalPropsQuadPt, 
		       &totalStrain[iQuad*_tensorSize], _tensorSize);

  return _elasticConsts;
} // calcDerivElastic

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
		      &totalStrain[iQuad*_tensorSize], _tensorSize);
  
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
  const real_section_type::value_type* parameterCell =
    _properties->restrictPoint(cell);
  memcpy(&_propertiesCell[0], parameterCell, 
		      numQuadPts*totalPropsQuadPt*sizeof(double));
} // _getProperties

// ----------------------------------------------------------------------
// Update properties (for next time step).
void
pylith::materials::ElasticMaterial::_updateProperties(double* const properties,
						      const int totalPropsQuadPt,
						      const double* totalStrain,
						      const int strainSize)
{ // _updateProperties
} // _updateProperties


// End of file 
