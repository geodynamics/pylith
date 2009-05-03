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
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

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
pylith::materials::ElasticMaterial::initialize(
			    const topology::Mesh& mesh,
			    feassemble::Quadrature<topology::Mesh>* quadrature)
{ // initialize
  Material::initialize(mesh, quadrature);

  assert(0 != quadrature);
  _numQuadPts = quadrature->numQuadPts();

  _initializeInitialStress(mesh, quadrature);
  _initializeInitialStrain(mesh, quadrature);
  _allocateCellArrays();
} // initialize

// ----------------------------------------------------------------------
// Retrieve parameters for physical properties and state variables for cell.
void
pylith::materials::ElasticMaterial::retrievePropsAndVars(const int cell)
{ // retrievePropsAndVars
  assert(0 != _properties);
  assert(0 != _stateVars);

  const ALE::Obj<RealSection>& propertiesSection = _properties->section();
  assert(!propertiesSection.isNull());
  propertiesSection->restrictPoint(cell, &_propertiesCell[0],
				   _propertiesCell.size());

  if (hasStateVars()) {
    const ALE::Obj<RealSection>& stateVarsSection = _stateVars->section();
    assert(!stateVarsSection.isNull());
    stateVarsSection->restrictPoint(cell, &_stateVarsCell[0],
				    _stateVarsCell.size());
  } // if

  if (0 == _initialStress)
    _initialStressCell = 0.0;
  else {
    const ALE::Obj<RealSection>& stressSection = _initialStress->section();
    assert(!stressSection.isNull());
    stressSection->restrictPoint(cell, &_initialStressCell[0],
				 _initialStressCell.size());
  } // if/else
  if (0 == _initialStrain)
    _initialStrainCell = 0.0;
  else {
    const ALE::Obj<RealSection>& strainSection = _initialStrain->section();
    assert(!strainSection.isNull());
    strainSection->restrictPoint(cell, &_initialStrainCell[0],
				 _initialStrainCell.size());
  } // if/else
} // retrievePropsAndVars

// ----------------------------------------------------------------------
// Compute density for cell at quadrature points.
const pylith::double_array&
pylith::materials::ElasticMaterial::calcDensity(void)
{ // calcDensity
  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*numPropsQuadPt);
  assert(_stateVarsCell.size() == numQuadPts*numVarsQuadPt);
  assert(_densityCell.size() == numQuadPts*1);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcDensity(&_densityCell[iQuad], 
		 &_propertiesCell[iQuad*numPropsQuadPt], numPropsQuadPt,
		 &_stateVarsCell[iQuad*numVarsQuadPt], numVarsQuadPt);

  return _densityCell;
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor for cell at quadrature points.
const pylith::double_array&
pylith::materials::ElasticMaterial::calcStress(const double_array& totalStrain,
					       const bool computeStateVars)
{ // calcStress
  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int tensorSize = _tensorSize;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*numPropsQuadPt);
  assert(_stateVarsCell.size() == numQuadPts*numVarsQuadPt);
  assert(_stressCell.size() == numQuadPts*_tensorSize);
  assert(_initialStressCell.size() == numQuadPts*_tensorSize);
  assert(_initialStrainCell.size() == numQuadPts*_tensorSize);
  assert(totalStrain.size() == numQuadPts*_tensorSize);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcStress(&_stressCell[iQuad*_tensorSize], _tensorSize,
		&_propertiesCell[iQuad*numPropsQuadPt], numPropsQuadPt,
		&_stateVarsCell[iQuad*numVarsQuadPt], numVarsQuadPt,
		&totalStrain[iQuad*_tensorSize], _tensorSize, 
		&_initialStressCell[iQuad*_tensorSize], _tensorSize,
		&_initialStrainCell[iQuad*_tensorSize], _tensorSize,
		computeStateVars);

  return _stressCell;
} // calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix for cell at quadrature points.
const pylith::double_array&
pylith::materials::ElasticMaterial::calcDerivElastic(
					       const double_array& totalStrain)
{ // calcDerivElastic
  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int tensorSize = _tensorSize;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*numPropsQuadPt);
  assert(_stateVarsCell.size() == numQuadPts*numVarsQuadPt);
  assert(_elasticConstsCell.size() == numQuadPts*_numElasticConsts);
  assert(_initialStressCell.size() == numQuadPts*_tensorSize);
  assert(_initialStrainCell.size() == numQuadPts*_tensorSize);
  assert(totalStrain.size() == numQuadPts*_tensorSize);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcElasticConsts(&_elasticConstsCell[iQuad*_numElasticConsts], 
		       _numElasticConsts,
		       &_propertiesCell[iQuad*numPropsQuadPt], 
		       numPropsQuadPt, 
		       &_stateVarsCell[iQuad*numVarsQuadPt], numVarsQuadPt,
		       &totalStrain[iQuad*_tensorSize], _tensorSize,
		       &_initialStressCell[iQuad*_tensorSize], _tensorSize,
		       &_initialStrainCell[iQuad*_tensorSize], _tensorSize);

  return _elasticConstsCell;
} // calcDerivElastic

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::materials::ElasticMaterial::updateStateVars(
					      const double_array& totalStrain,
					      const int cell)
{ // updateStateVars
  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int tensorSize = _tensorSize;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*numPropsQuadPt);
  assert(_stateVarsCell.size() == numQuadPts*numVarsQuadPt);
  assert(_initialStressCell.size() == numQuadPts*_tensorSize);
  assert(_initialStrainCell.size() == numQuadPts*_tensorSize);
  assert(totalStrain.size() == numQuadPts*_tensorSize);

  const ALE::Obj<RealSection>& stateVarsSection = _stateVars->section();
  assert(!stateVarsSection.isNull());
  stateVarsSection->restrictPoint(cell, &_stateVarsCell[0],
				  _stateVarsCell.size());
  
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _updateStateVars(&_stateVarsCell[iQuad*numVarsQuadPt], numVarsQuadPt,
		     &_propertiesCell[iQuad*numPropsQuadPt], 
		     numPropsQuadPt,
		     &totalStrain[iQuad*_tensorSize], _tensorSize,
		     &_initialStressCell[iQuad*_tensorSize], _tensorSize,
		     &_initialStrainCell[iQuad*_tensorSize], _tensorSize);
  
  stateVarsSection->updatePoint(cell, &_stateVarsCell[0]);
} // updateStateVars

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::ElasticMaterial::stableTimeStepImplicit(const topology::Mesh& mesh)
{ // stableTimeStepImplicit
  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int tensorSize = _tensorSize;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*numPropsQuadPt);
  assert(_stateVarsCell.size() == numQuadPts*numVarsQuadPt);
  assert(_elasticConstsCell.size() == numQuadPts*_numElasticConsts);
  assert(_initialStressCell.size() == numQuadPts*_tensorSize);
  assert(_initialStrainCell.size() == numQuadPts*_tensorSize);

  double dtStable = pylith::PYLITH_MAXDOUBLE;

  // Get cells associated with material
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    retrievePropsAndVars(*c_iter);

    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double dt = 
	_stableTimeStepImplicit(&_propertiesCell[iQuad*numPropsQuadPt],
				numPropsQuadPt,
				&_stateVarsCell[iQuad*numVarsQuadPt],
				numVarsQuadPt);
      if (dt < dtStable)
	dtStable = dt;
    } // for
  } // for
  
  return dtStable;
} // stableTimeStepImplicit

// ----------------------------------------------------------------------
// Get initial stress field.
const pylith::topology::Field<pylith::topology::Mesh>&
pylith::materials::ElasticMaterial::initialStressField(void) const
{ // initialStressField
  if (0 == _initialStress)
    throw std::runtime_error("Request for initial stress field, but field "
			     "does not exist.");

  return *_initialStress;
} // initialStressField

// ----------------------------------------------------------------------
// Get initial strain field.
const pylith::topology::Field<pylith::topology::Mesh>&
pylith::materials::ElasticMaterial::initialStrainField(void) const
{ // initialStrainField
  if (0 == _initialStrain)
    throw std::runtime_error("Request for initial strain field, but field "
			     "does not exist.");

  return *_initialStrain;
} // initialStrainField

// ----------------------------------------------------------------------
// Allocate cell arrays.
void
pylith::materials::ElasticMaterial::_allocateCellArrays(void)
{ // _allocateCellArrays
  const int numQuadPts = _numQuadPts;
  const int tensorSize = _tensorSize;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int numVarsQuadPt = _numVarsQuadPt;
  const int numElasticConsts = _numElasticConsts;

  _propertiesCell.resize(numQuadPts * numPropsQuadPt);
  _stateVarsCell.resize(numQuadPts * numVarsQuadPt);
  _initialStressCell.resize(numQuadPts * tensorSize);
  _initialStrainCell.resize(numQuadPts * tensorSize);
  _densityCell.resize(numQuadPts);
  _stressCell.resize(numQuadPts * tensorSize);
  _elasticConstsCell.resize(numQuadPts * numElasticConsts);
} // _allocateCellArrays

// ----------------------------------------------------------------------
// Initialize initial stress field.
void
pylith::materials::ElasticMaterial::_initializeInitialStress(
			 const topology::Mesh& mesh,
			 feassemble::Quadrature<topology::Mesh>* quadrature)
{ // _initializeInitialStress
  delete _initialStress; _initialStress = 0;
  if (0 == _dbInitialStress)
    return;

  assert(0 != _dbInitialStress);
  assert(0 != quadrature);

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();

  // Create arrays for querying
  const int tensorSize = _tensorSize;
  double_array quadPtsGlobal(numQuadPts*spaceDim);
  double_array stressCell(numQuadPts*tensorSize);

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
  _initialStress = new topology::Field<topology::Mesh>(mesh);
  assert(0 != _initialStress);
  const int fiberDim = numQuadPts * tensorSize;
  assert(fiberDim > 0);
  _initialStress->newSection(cells, fiberDim);
  _initialStress->allocate();
  _initialStress->zero();
  const ALE::Obj<RealSection>& initialStressSection = 
    _initialStress->section();

  // Setup databases for querying
  _dbInitialStress->open();
  switch (dimension())
    { // switch
    case 1: {
      const char* stressDBValues[] = { "stress" };
      _dbInitialStress->queryVals(stressDBValues, tensorSize);
      break;
    } // case 1
    case 2 : {
      const char* stressDBValues[] = { 
	"stress-xx",
	"stress-yy",
	"stress-xy",
      };
      _dbInitialStress->queryVals(stressDBValues, tensorSize);
      break;
    } // case 2
    case 3 : {
      const char* stressDBValues[] = { 
	"stress-xx",
	"stress-yy",
	"stress-zz",
	"stress-xy",
	"stress-yz",
	"stress-xz",
      };
      _dbInitialStress->queryVals(stressDBValues, tensorSize);
      break;
    } // case 3
    default :
      assert(0);
    } // switch
  
  assert(0 != _normalizer);
  const double lengthScale = _normalizer->pressureScale();
  const double pressureScale = _normalizer->pressureScale();

  const ALE::Obj<RealSection>& stressSection = _initialStress->section();
  assert(!stressSection.isNull());
    
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    quadrature->retrieveGeometry(*c_iter);

    // Dimensionalize coordinates for querying
    const double_array& quadPtsNonDim = quadrature->quadPts();
    quadPtsGlobal = quadPtsNonDim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
				lengthScale);
    
    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, iCoord=0, iStress=0; 
	 iQuadPt < numQuadPts; 
	 ++iQuadPt, iCoord+=spaceDim, iStress+=tensorSize) {
      int err = _dbInitialStress->query(&stressCell[iStress], tensorSize,
					&quadPtsGlobal[iCoord], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find initial stress at (";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPtsGlobal[iCoord+i];
	msg << ") in material " << label() << "\n"
	    << "using spatial database '" << _dbInitialStress->label() << "'.";
	throw std::runtime_error(msg.str());
      } // if
    } // for

    // Nondimensionalize stress
    _normalizer->nondimensionalize(&stressCell[0], stressCell.size(), 
				   pressureScale);

    stressSection->updatePoint(*c_iter, &stressCell[0]);
  } // for

  // Close databases
  _dbInitialStress->close();
} // _initializeInitialStress

// ----------------------------------------------------------------------
// Initialize initial strain field.
void
pylith::materials::ElasticMaterial::_initializeInitialStrain(
			 const topology::Mesh& mesh,
			 feassemble::Quadrature<topology::Mesh>* quadrature)
{ // _initializeInitialStrain
  delete _initialStrain; _initialStrain = 0;
  if (0 == _dbInitialStrain)
    return;

  assert(0 != _dbInitialStrain);
  assert(0 != quadrature);

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();

  // Create arrays for querying
  const int tensorSize = _tensorSize;
  double_array quadPtsGlobal(numQuadPts*spaceDim);
  double_array strainCell(numQuadPts*tensorSize);

  // Get cells associated with material
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();

  // Create field to hold initial strain state.
  delete _initialStrain; _initialStrain = 0;
  _initialStrain = new topology::Field<topology::Mesh>(mesh);
  assert(0 != _initialStrain);
  const int fiberDim = numQuadPts * tensorSize;
  assert(fiberDim > 0);
  _initialStrain->newSection(cells, fiberDim);
  _initialStrain->allocate();
  _initialStrain->zero();
  const ALE::Obj<RealSection>& initialStrainSection = 
    _initialStrain->section();

  // Setup databases for querying
  _dbInitialStrain->open();
  switch (dimension())
    { // switch
    case 1: {
      const char* strainDBValues[] = { "strain" };
      _dbInitialStrain->queryVals(strainDBValues, tensorSize);
      break;
    } // case 1
    case 2 : {
      const char* strainDBValues[] = { 
	"strain-xx",
	"strain-yy",
	"strain-xy",
      };
      _dbInitialStrain->queryVals(strainDBValues, tensorSize);
      break;
    } // case 2
    case 3 : {
      const char* strainDBValues[] = { 
	"strain-xx",
	"strain-yy",
	"strain-zz",
	"strain-xy",
	"strain-yz",
	"strain-xz",
      };
      _dbInitialStrain->queryVals(strainDBValues, tensorSize);
      break;
    } // case 3
    default :
      assert(0);
    } // switch
  
  assert(0 != _normalizer);
  const double lengthScale = _normalizer->pressureScale();
  const double pressureScale = _normalizer->pressureScale();
    
  const ALE::Obj<RealSection>& strainSection = _initialStrain->section();
  assert(!strainSection.isNull());
    
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    quadrature->retrieveGeometry(*c_iter);

    // Dimensionalize coordinates for querying
    const double_array& quadPtsNonDim = quadrature->quadPts();
    quadPtsGlobal = quadPtsNonDim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
				lengthScale);
    
    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, iCoord=0, iStrain=0; 
	 iQuadPt < numQuadPts; 
	 ++iQuadPt, iCoord+=spaceDim, iStrain+=tensorSize) {
      int err = _dbInitialStrain->query(&strainCell[iStrain], tensorSize,
					&quadPtsGlobal[iCoord], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find initial strain at (";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPtsGlobal[iCoord+i];
	msg << ") in material " << label() << "\n"
	    << "using spatial database '" << _dbInitialStrain->label() << "'.";
	throw std::runtime_error(msg.str());
      } // if
    } // for

    // Nondimensionalize strain
    _normalizer->nondimensionalize(&strainCell[0], strainCell.size(), 
				   pressureScale);

    strainSection->updatePoint(*c_iter, &strainCell[0]);
  } // for

  // Close databases
  _dbInitialStrain->close();
} // _initializeInitialStrain

// ----------------------------------------------------------------------
// Update stateVars (for next time step).
void
pylith::materials::ElasticMaterial::_updateStateVars(
					    double* const stateVars,
					    const int numStateVars,
					    const double* properties,
					    const int numProperties,
					    const double* totalStrain,
					    const int strainSize,
					    const double* initialStress,
					    const int initialStressSize,
					    const double* initialStrain,
					    const int initialStrainSize)
{ // _updateStateVars
} // _updateStateVars


// End of file 
