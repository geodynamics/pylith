// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "ElasticMaterial.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/utils/array.hh" // USES scalar_array, std::vector
#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES memcpy()
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

typedef pylith::topology::Field<pylith::topology::Mesh>::RestrictVisitor RestrictVisitor;
typedef pylith::topology::Field<pylith::topology::Mesh>::UpdateAddVisitor UpdateAddVisitor;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticMaterial::ElasticMaterial(const int dimension,
						    const int tensorSize,
						    const int numElasticConsts,
						    const Metadata& metadata) :
  Material(dimension, tensorSize, metadata),
  _dbInitialStress(0),
  _dbInitialStrain(0),
  _initialFields(0),
  _numQuadPts(0),
  _numElasticConsts(numElasticConsts)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticMaterial::~ElasticMaterial(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::ElasticMaterial::deallocate(void)
{ // deallocate
  delete _initialFields; _initialFields = 0;

  _dbInitialStress = 0; // :TODO: Use shared pointer.
  _dbInitialStrain = 0; // :TODO: Use shared pointer.
} // deallocate
  
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

  if (0 != _dbInitialStress || 0 != _dbInitialStrain) {
    delete _initialFields; 
    _initialFields =
      new topology::Fields<topology::Field<topology::Mesh> >(mesh);
  } // if
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

  _initialStressCell = 0.0;
  _initialStrainCell = 0.0;
  if (0 != _initialFields) {
    if (_initialFields->hasField("initial stress")) {
      const ALE::Obj<RealSection>& stressSection = 
	_initialFields->get("initial stress").section();
      assert(!stressSection.isNull());
      stressSection->restrictPoint(cell, &_initialStressCell[0],
				   _initialStressCell.size());
    } // if
    if (_initialFields->hasField("initial strain")) {
      const ALE::Obj<RealSection>& strainSection =
	_initialFields->get("initial strain").section();
      assert(!strainSection.isNull());
      strainSection->restrictPoint(cell, &_initialStrainCell[0],
				   _initialStrainCell.size());
    } // if
  } // if
} // retrievePropsAndVars

// ----------------------------------------------------------------------
// Compute stress tensor for cell at quadrature points.
const pylith::scalar_array&
pylith::materials::ElasticMaterial::calcStress(const scalar_array& totalStrain,
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
const pylith::scalar_array&
pylith::materials::ElasticMaterial::calcDerivElastic(
					       const scalar_array& totalStrain)
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
					      const scalar_array& totalStrain,
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

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _updateStateVars(&_stateVarsCell[iQuad*numVarsQuadPt], numVarsQuadPt,
		     &_propertiesCell[iQuad*numPropsQuadPt], 
		     numPropsQuadPt,
		     &totalStrain[iQuad*_tensorSize], _tensorSize,
		     &_initialStressCell[iQuad*_tensorSize], _tensorSize,
		     &_initialStrainCell[iQuad*_tensorSize], _tensorSize);
  
  const ALE::Obj<RealSection>& stateVarsSection = _stateVars->section();
  assert(!stateVarsSection.isNull());  
  stateVarsSection->updatePoint(cell, &_stateVarsCell[0]);
} // updateStateVars

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
PylithScalar
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

  PylithScalar dtStable = pylith::PYLITH_MAXDOUBLE;

  // Get cells associated with material
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    retrievePropsAndVars(*c_iter);

    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar dt = 
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
  if (0 == _dbInitialStress)
    return;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Materials");

  assert(0 != _initialFields);
  _initialFields->add("initial stress", "initial_stress");
  topology::Field<topology::Mesh>& initialStress = 
    _initialFields->get("initial stress");

  assert(0 != _dbInitialStress);
  assert(0 != quadrature);

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();
  const int numBasis = quadrature->numBasis();

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  // Create arrays for querying
  const int tensorSize = _tensorSize;
  scalar_array quadPtsGlobal(numQuadPts*spaceDim);
  scalar_array stressCell(numQuadPts*tensorSize);

  // Get cells associated with material
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();

  // Create field to hold initial stress state.
  const int fiberDim = numQuadPts * tensorSize;
  assert(fiberDim > 0);
  initialStress.newSection(cells, fiberDim);
  initialStress.allocate();
  initialStress.zero();
  const ALE::Obj<RealSection>& initialStressSection = initialStress.section();

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
      std::cerr << "Bad dimension '" << dimension() << "'." << std::endl;
      assert(0);
      throw std::logic_error("Unknown dimension in elastic material.");
    } // switch
  
  assert(0 != _normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);
    quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    // Dimensionalize coordinates for querying
    const scalar_array& quadPtsNonDim = quadrature->quadPts();
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

    initialStressSection->updatePoint(*c_iter, &stressCell[0]);
  } // for

  // Close databases
  _dbInitialStress->close();

  logger.stagePop();
} // _initializeInitialStress

// ----------------------------------------------------------------------
// Initialize initial strain field.
void
pylith::materials::ElasticMaterial::_initializeInitialStrain(
			 const topology::Mesh& mesh,
			 feassemble::Quadrature<topology::Mesh>* quadrature)
{ // _initializeInitialStrain
  if (0 == _dbInitialStrain)
    return;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Materials");

  assert(0 != _initialFields);
  _initialFields->add("initial strain", "initial_strain");
  topology::Field<topology::Mesh>& initialStrain = 
    _initialFields->get("initial strain");

  assert(0 != _dbInitialStrain);
  assert(0 != quadrature);

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();
  const int numBasis = quadrature->numBasis();

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  // Create arrays for querying
  const int tensorSize = _tensorSize;
  scalar_array quadPtsGlobal(numQuadPts*spaceDim);
  scalar_array strainCell(numQuadPts*tensorSize);

  // Get cells associated with material
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();

  // Create field to hold initial strain state.
  const int fiberDim = numQuadPts * tensorSize;
  assert(fiberDim > 0);
  initialStrain.newSection(cells, fiberDim);
  initialStrain.allocate();
  initialStrain.zero();
  const ALE::Obj<RealSection>& initialStrainSection = initialStrain.section();

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
      std::cerr << "Bad dimension '" << dimension() << "'." << std::endl;
      assert(0);
      throw std::logic_error("Unknown dimension in elastic material.");
    } // switch
  
  assert(0 != _normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
    
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);
    quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    // Dimensionalize coordinates for querying
    const scalar_array& quadPtsNonDim = quadrature->quadPts();
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

    initialStrainSection->updatePoint(*c_iter, &strainCell[0]);
  } // for

  // Close databases
  _dbInitialStrain->close();

  logger.stagePop();
} // _initializeInitialStrain

// ----------------------------------------------------------------------
// Update stateVars (for next time step).
void
pylith::materials::ElasticMaterial::_updateStateVars(
					    PylithScalar* const stateVars,
					    const int numStateVars,
					    const PylithScalar* properties,
					    const int numProperties,
					    const PylithScalar* totalStrain,
					    const int strainSize,
					    const PylithScalar* initialStress,
					    const int initialStressSize,
					    const PylithScalar* initialStrain,
					    const int initialStrainSize)
{ // _updateStateVars
} // _updateStateVars


// End of file 
