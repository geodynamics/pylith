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

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/utils/array.hh" // USES double_array, std::vector
#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES memcpy()
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticMaterial::ElasticMaterial(const int dimension,
						    const int tensorSize,
						    const int numElasticConsts,
						    const Metadata& metadata) :
  Material(dimension, tensorSize, metadata),
  _dbInitialStress(0),
  _dbInitialStrain(0),
  _initialStress(0),
  _initialStrain(0),
  _numQuadPts(0),
  _numElasticConsts(numElasticConsts)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticMaterial::~ElasticMaterial(void)
{ // destructor
  delete _initialStress; _initialStress = 0;
  delete _initialStrain; _initialStrain = 0;

  // Python db object owns databases, so just set point to null
  // :KLUDGE: Should use shared pointer
  _dbInitialStress = 0;
  _dbInitialStrain = 0;
} // destructor

// ----------------------------------------------------------------------
// Initialize material by getting physical property parameters from
// database.
void
initialize(const topology::Mesh& mesh,
	   feassemble::Quadrature<topology::Mesh>* quadrature)
{ // initialize
  Material::initialize(mesh, quadrature);

  delete _initialStress; _initialStress = 0;
  delete _initialStrain; _initialStrain = 0;

  if (0 != _dbInitialStress)
    _initializeInitialStress(mesh, quadrature);

  if (0 != _dbInitialStrain)
    _initializeInitialStrain(mesh, quadrature);
} // initialize

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();

  // Create arrays for querying
  const int tensorSize = tensorSize();
  double_array quadPtsGlobal(numQuadPts*spaceDim);
  double_array valuesQuery(tensorSize);
  double_array valuesCell(numQuadPts*tensorSize);

  // Get cells associated with material
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();

  // Create field to hold initial stress state.
  delete _initialStress; _initialStress = 0;
  if (0 != _dbInitialStress) {
    _initialStress = new topology::Field<topology::Mesh>(mesh);
    assert(0 != _initialStress);
    fiberDim = numQuadPts * tensorSize;
    _initialStress->newSection(cells, fiberDim);
    _initialStress->allocate();
    _initialStress->zero();
  } // if
  const ALE::Obj<RealSection>& initialStressSection = 
    (0 != _initialStress) ? _initialStress->section() : 0;

  // Create field to hold initial strain state.
  delete _initialStrain; _initialStrain = 0;
  if (0 != _dbInitialStrain) {
    _initialStrain = new topology::Field<topology::Mesh>(mesh);
    assert(0 != _initialStrain);
    fiberDim = numQuadPts * tensorSize;
    if (0 == _initialStress)
      _initialStrain->newSection(cells, fiberDim);
    else { // reuse chart from inital stress
      const ALE::Obj<RealSection::chart_type>& chart = 
	initialStressSection->getChart();
      assert(!chart.isNull());
      _initialStrain->newSection(*chart, fiberDim);
    } // else
    _initialStrain->allocate();
    _initialStrain->zero();
  } // if
  const ALE::Obj<RealSection>& initialStrainSection = 
    (0 != _initialStrain) ? _initialStrain->section() : 0;

  // Setup databases as necessary.
  if (0 != _dbInitialStress) {
    _dbInitialStress->open();
    _dbInitialStress->queryVals(stressDBValues, tensorSize);
  } // if
  if (0 != _dbInitialStrain) {
    _dbInitialStrain->open();
    _dbInitialStrain->queryVals(strainDBValues, tensorSize);
  } // if

  assert(0 != _normalizer);
  const double pressureScale = _normalizer->pressureScale();
  const double lengthScale = _normalizer->pressureScale();
    
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    quadrature->retrieveGeometry(*c_iter);

    const double_array& quadPtsNonDim = quadrature->quadPts();
    quadPtsGlobal = quadPtsNonDim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
				lengthScale);

    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, index=0; 
	 iQuadPt < numQuadPts; 
	 ++iQuadPt, index+=spaceDim) {
      int err = _dbInitalStress->query(&valuesQuery[0], tensorSize,
				       &quadPtsGlobal[index], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find initial stress state at \n"
	    << "(";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPtsGlobal[index+i];
	msg << ") in material " << _label << "\n"
	    << "using spatial database '" << _dbInitialStress->label() << "'.";
	throw std::runtime_error(msg.str());
      } // if
      _dbToProperties(&propertiesCell[iQuadPt*_numPropsQuadPt], 
		      propertiesQuery);
      _nondimProperties(&propertiesCell[iQuadPt*_numPropsQuadPt],
			_numPropsQuadPt);

      if (0 != _dbInitialState) {
	err = _dbInitialState->query(&stateVarsQuery[0], numDBStateVars,
				     &quadPtsGlobal[index], spaceDim, cs);
	if (err) {
	  std::ostringstream msg;
	  msg << "Could not find initial state values at \n" << "(";
	  for (int i=0; i < spaceDim; ++i)
	    msg << "  " << quadPtsGlobal[index+i];
	  msg << ") in material " << _label << "\n"
	      << "using spatial database '" << _dbInitialState->label()
	      << "'.";
	  throw std::runtime_error(msg.str());
	} // if
	_dbToStateVars(&stateVarsCell[iQuadPt*_numVarsQuadPt], 
		       stateVarsQuery);
	_nondimStateVars(&stateVarsCell[iQuadPt*_numVarsQuadPt],
			 _numVarsQuadPt);
      } // if

    } // for
    // Insert cell contribution into fields
    propertiesSection->updatePoint(*c_iter, &propertiesCell[0]);
    if (0 != _dbInitialState)
      stateVarsSection->updatePoint(*c_iter, &stateVarsCell[0]);
  } // for

  // Close databases
  _dbProperties->close();
  if (0 != _dbInitialState)
    _dbInitialState->close();
  
} // initialize
  
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
