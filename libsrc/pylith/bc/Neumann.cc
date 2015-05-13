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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "Neumann.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Fields.hh" // HOLDSA Fields
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/VisitorSubMesh.hh" // USES VecVisitorSubMesh
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES strcpy()
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::Neumann::Neumann(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::Neumann::~Neumann(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::Neumann::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  BCIntegratorSubMesh::deallocate();
  TimeDependent::deallocate();

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
pylith::bc::Neumann::initialize(const topology::Mesh& mesh,
				const PylithScalar upDir[3])
{ // initialize
  PYLITH_METHOD_BEGIN;

  _queryDatabases();
  _paramsLocalToGlobal(upDir);

  // Optimize coordinate retrieval in closure
  assert(_boundaryMesh);
  const PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::CoordsVisitor::optimizeClosure(dmSubMesh);

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::Neumann::integrateResidual(const topology::Field& residual,
				       const PylithScalar t,
				       topology::SolutionFields* const fields)
{ // integrateResidual
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);
  assert(_boundaryMesh);
  assert(_parameters);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == size_t(numQuadPts));
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Allocate vectors for cell values.
  _initCellVector();
  scalar_array tractionsCell(numQuadPts*spaceDim);

  // Get cell information
  PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::Stratum cellsStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  // Get sections
  _calculateValue(t);
  topology::Field& valueField = _parameters->get("value");
  topology::VecVisitorMesh valueVisitor(valueField);
  PetscScalar* valueArray = valueVisitor.localArray();

  // Get subsections
  topology::SubMeshIS submeshIS(*_boundaryMesh);
  topology::VecVisitorSubMesh residualVisitor(residual, submeshIS);
  submeshIS.deallocate();

  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmSubMesh);

  // Loop over faces and integrate contribution from each face
  for(PetscInt c = cStart; c < cEnd; ++c) {
    coordsVisitor.getClosure(&coordsCell, c);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), c);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict tractions to cell
    const PetscInt voff = valueVisitor.sectionOffset(c);
    assert(numQuadPts*spaceDim == valueVisitor.sectionDof(c));

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();

    // Compute action for traction bc terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const PylithScalar valI = wt*basis[iQuad*numBasis+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const PylithScalar valIJ = valI * basis[iQuad*numBasis+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis*spaceDim+iDim] += valueArray[voff+iQuad*spaceDim+iDim] * valIJ;
        } // for
      } // for
    } // for

    residualVisitor.setClosure(&_cellVector[0], _cellVector.size(), c, ADD_VALUES);

    PetscLogFlops(numQuadPts*(1+numBasis*(1+numBasis*(1+2*spaceDim))));
  } // for

  PYLITH_METHOD_END;
} // integrateResidual

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::Neumann::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  PYLITH_METHOD_BEGIN;

  BCIntegratorSubMesh::verifyConfiguration(mesh);

  PYLITH_METHOD_END;
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get cell field for tractions.
const pylith::topology::Field&
pylith::bc::Neumann::cellField(const char* name,
			       topology::SolutionFields* const fields)
{ // cellField
  PYLITH_METHOD_BEGIN;

  assert(_parameters);
  assert(name);

  if (0 == strcasecmp(name, "initial_value"))
    PYLITH_METHOD_RETURN(_parameters->get("initial"));

  else if (0 == strcasecmp(name, "rate_of_change"))
    PYLITH_METHOD_RETURN(_parameters->get("rate"));

  else if (0 == strcasecmp(name, "change_in_value"))
    PYLITH_METHOD_RETURN(_parameters->get("change"));

  else if (0 == strcasecmp(name, "rate_start_time"))
    PYLITH_METHOD_RETURN(_parameters->get("rate time"));

  else if (0 == strcasecmp(name, "change_start_time"))
    PYLITH_METHOD_RETURN(_parameters->get("change time"));

  else {
    std::ostringstream msg;
    msg << "Unknown field '" << name << "' requested for Neumann BC '" 
	<< _label << "'.";
    throw std::runtime_error(msg.str());
  } // else

  // Satisfy method definition
  PYLITH_METHOD_RETURN(_parameters->get("traction"));
} // cellField

// ----------------------------------------------------------------------
// Query databases for parameters.
void
pylith::bc::Neumann::_queryDatabases(void)
{ // _queryDatabases
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);
  assert(_boundaryMesh);
  
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar rateScale = pressureScale / timeScale;

  const int spaceDim = _quadrature->spaceDim();
  const int numQuadPts = _quadrature->numQuadPts();

  delete _parameters; _parameters = new topology::Fields(*_boundaryMesh);assert(_parameters);

  // Create section to hold time dependent values
  _parameters->add("value", "traction", topology::FieldBase::FACES_FIELD, numQuadPts*spaceDim);
  topology::Field& value = _parameters->get("value");
  value.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
  value.scale(pressureScale);
  value.allocate();
  if (_dbInitial) {
    _parameters->add("initial", "initial_traction", topology::FieldBase::FACES_FIELD, numQuadPts*spaceDim);
    topology::Field& initial = _parameters->get("initial");
    initial.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
    initial.scale(pressureScale);
    initial.allocate();
  }
  if (_dbRate) {
    _parameters->add("rate", "traction_rate", topology::FieldBase::FACES_FIELD, numQuadPts*spaceDim);
    topology::Field& rate = _parameters->get("rate");
    rate.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
    rate.scale(rateScale);
    rate.allocate();
    _parameters->add("rate time", "traction_rate_time", topology::FieldBase::FACES_FIELD, numQuadPts);
    topology::Field& rateTime = _parameters->get("rate time");
    rateTime.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
    rateTime.scale(timeScale);
    rateTime.allocate();
  } // if
  if (_dbChange) {
    _parameters->add("change", "change_traction", topology::FieldBase::FACES_FIELD, numQuadPts*spaceDim);
    topology::Field& change = _parameters->get("change");
    change.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
    change.scale(pressureScale);
    change.allocate();
    _parameters->add("change time", "change_traction_time", topology::FieldBase::FACES_FIELD, numQuadPts);
    topology::Field& changeTime = _parameters->get("change time");
    changeTime.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
    changeTime.scale(timeScale);
    changeTime.allocate();
  } // if

  if (_dbInitial) { // Setup initial values, if provided.
    _dbInitial->open();
    switch (spaceDim)
      { // switch
      case 1 : {
	const char* valueNames[] = {"traction-normal"};
	_dbInitial->queryVals(valueNames, 1);
	break;
      } // case 1
      case 2 : {
	const char* valueNames[] = {"traction-shear", "traction-normal"};
	_dbInitial->queryVals(valueNames, 2);
	break;
      } // case 2
      case 3 : {
	const char* valueNames[] = {"traction-shear-horiz",
				    "traction-shear-vert",
				    "traction-normal"};
	_dbInitial->queryVals(valueNames, 3);
	break;
      } // case 3
      default :
        std::ostringstream msg;
        msg << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
        throw std::logic_error(msg.str());
      } // switch
    _queryDB("initial", _dbInitial, spaceDim, pressureScale);
    _dbInitial->close();
  } // if

  if (_dbRate) { // Setup rate of change of values, if provided.
    _dbRate->open();
    switch (spaceDim)
      { // switch
      case 1 : {
	const char* valueNames[] = {"traction-rate-normal"};
	_dbRate->queryVals(valueNames, 1);
	break;
      } // case 1
      case 2 : {
	const char* valueNames[] = {"traction-rate-shear", 
				    "traction-rate-normal"};
	_dbRate->queryVals(valueNames, 2);
	break;
      } // case 2
      case 3 : {
	const char* valueNames[] = {"traction-rate-shear-horiz",
				    "traction-rate-shear-vert",
				    "traction-rate-normal"};
	_dbRate->queryVals(valueNames, 3);
	break;
      } // case 3
      default :
        std::ostringstream msg;
        msg << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
        throw std::logic_error(msg.str());
      } // switch
    _queryDB("rate", _dbRate, spaceDim, rateScale);

    const char* timeNames[1] = { "rate-start-time" };
    _dbRate->queryVals(timeNames, 1);
    _queryDB("rate time", _dbRate, 1, timeScale);
    _dbRate->close();
  } // if

  if (_dbChange) { // Setup change of values, if provided.
    _dbChange->open();
    switch (spaceDim)
      { // switch
      case 1 : {
	const char* valueNames[] = {"traction-normal"};
	_dbChange->queryVals(valueNames, 1);
	break;
      } // case 1
      case 2 : {
	const char* valueNames[] = {"traction-shear", "traction-normal"};
	_dbChange->queryVals(valueNames, 2);
	break;
      } // case 2
      case 3 : {
	const char* valueNames[] = {"traction-shear-horiz",
				    "traction-shear-vert",
				    "traction-normal"};
	_dbChange->queryVals(valueNames, 3);
	break;
      } // case 3
      default :
        std::ostringstream msg;
        msg << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
        throw std::logic_error(msg.str());
      } // switch
    _queryDB("change", _dbChange, spaceDim, pressureScale);

    const char* timeNames[1] = { "change-start-time" };
    _dbChange->queryVals(timeNames, 1);
    _queryDB("change time", _dbChange, 1, timeScale);
    _dbChange->close();

    if (_dbTimeHistory)
      _dbTimeHistory->open();
  } // if

  PYLITH_METHOD_END;
} // _queryDatabases

// ----------------------------------------------------------------------
// Query database for values.
void
pylith::bc::Neumann::_queryDB(const char* name,
			      spatialdata::spatialdb::SpatialDB* const db,
			      const int querySize,
			      const PylithScalar scale)
{ // _queryDB
  PYLITH_METHOD_BEGIN;

  assert(name);
  assert(db);
  assert(_boundaryMesh);
  assert(_quadrature);
  assert(_parameters);

  // Get 'surface' cells (1 dimension lower than top-level cells)
  PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::Stratum cellsStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = _quadrature->spaceDim();
  
  // Containers for database query results and quadrature coordinates in
  // reference geometry.
  scalar_array valuesCell(numQuadPts*querySize);
  scalar_array quadPtsGlobal(numQuadPts*spaceDim);

  // Get sections.
  topology::Field& valueField = _parameters->get(name);
  topology::VecVisitorMesh valueVisitor(valueField);
  PetscScalar* valueArray = valueVisitor.localArray();

  // Get coordinates
  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmSubMesh);

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();
  assert(cs);

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();

  // Compute quadrature information
  _quadrature->initializeGeometry();

  // Loop over cells in boundary mesh and perform queries.
  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, c);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), c);

    const scalar_array& quadPtsNondim = _quadrature->quadPts();
    quadPtsGlobal = quadPtsNondim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(), lengthScale);
    
    valuesCell = 0.0;
    for (int iQuad=0, iSpace=0; iQuad < numQuadPts; ++iQuad, iSpace+=spaceDim) {
      const int err = db->query(&valuesCell[iQuad*querySize], querySize,
                                &quadPtsGlobal[iSpace], spaceDim, cs);
      
      if (err) {
        std::ostringstream msg;
        msg << "Could not find values at (";
        for (int i=0; i < spaceDim; ++i)
          msg << " " << quadPtsGlobal[i+iSpace];
        msg << ") for traction boundary condition '" << _label
            << "' using spatial database '" << db->label() << "'.";
        throw std::runtime_error(msg.str());
      } // if
    } // for
    _normalizer->nondimensionalize(&valuesCell[0], valuesCell.size(), scale);

    // Update section
    const PetscInt voff = valueVisitor.sectionOffset(c);
    const PetscInt vdof = valueVisitor.sectionDof(c);
    assert(numQuadPts*querySize == vdof);
    for(PetscInt d = 0; d < vdof; ++d)
      valueArray[voff+d] = valuesCell[d];
  } // for

  PYLITH_METHOD_END;
} // _queryDB

// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
  pylith::bc::Neumann::_paramsLocalToGlobal(const PylithScalar upDir[3])
{ // _paramsLocalToGlobal
  PYLITH_METHOD_BEGIN;

  assert(_boundaryMesh);
  assert(_parameters);
  assert(_quadrature);

  scalar_array up(3);
  for (int i=0; i < 3; ++i)
    up[i] = upDir[i];

  // Get 'surface' cells (1 dimension lower than top-level cells)
  PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::Stratum cellsStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  // Quadrature related values.
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const int cellDim = _quadrature->cellDim() > 0 ? _quadrature->cellDim() : 1;
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = cellGeometry.spaceDim();
  const scalar_array& quadPtsRef = _quadrature->quadPtsRef();
  
  // Containers for orientation information
  const int orientationSize = spaceDim * spaceDim;
  const int jacobianSize = spaceDim * cellDim;
  scalar_array jacobian(jacobianSize);
  PylithScalar jacobianDet = 0;
  scalar_array orientation(orientationSize);

  // Get coordinates.
  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmSubMesh);

  // Get sections
  scalar_array tmpLocal(spaceDim);
  
  topology::Field* initialField = (_dbInitial) ? &_parameters->get("initial") : 0;
  topology::VecVisitorMesh* initialVisitor = (initialField) ? new topology::VecVisitorMesh(*initialField) : 0;
  PetscScalar* initialArray = (initialVisitor) ? initialVisitor->localArray() : NULL;

  topology::Field* rateField = (_dbRate) ? &_parameters->get("rate") : 0;
  topology::VecVisitorMesh* rateVisitor = (rateField) ? new topology::VecVisitorMesh(*rateField) : 0;
  PetscScalar* rateArray = (rateVisitor) ? rateVisitor->localArray() : NULL;

  topology::Field* changeField = (_dbChange) ? &_parameters->get("change") : 0;
  topology::VecVisitorMesh* changeVisitor = (changeField) ? new topology::VecVisitorMesh(*changeField) : 0;
  PetscScalar* changeArray = (changeVisitor) ? changeVisitor->localArray() : NULL;

  // Loop over cells in boundary mesh, compute orientations, and then
  // rotate corresponding traction vector from local to global coordinates.
  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, c);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), c);

    for(int iQuad=0, iRef=0, iSpace=0; iQuad < numQuadPts; ++iQuad, iRef+=cellDim, iSpace+=spaceDim) {
      // Compute Jacobian and determinant at quadrature point, then get orientation.
      cellGeometry.jacobian(&jacobian, &jacobianDet, &coordsCell[0], numBasis, spaceDim, &quadPtsRef[iRef], cellDim);
      cellGeometry.orientation(&orientation, jacobian, jacobianDet, up);
      assert(jacobianDet > 0.0);
      orientation /= jacobianDet;

      if (_dbInitial) {
        // Rotate traction vector from local coordinate system to global coordinate system
	assert(initialVisitor);
	const PetscInt ioff = initialVisitor->sectionOffset(c);
	assert(numQuadPts*spaceDim == initialVisitor->sectionDof(c));
        for(int iDim = 0; iDim < spaceDim; ++iDim) {
          tmpLocal[iDim] = initialArray[ioff+iSpace+iDim];
        } // for
        for(int iDim = 0; iDim < spaceDim; ++iDim) {
          initialArray[ioff+iSpace+iDim] = 0.0;
          for(int jDim = 0; jDim < spaceDim; ++jDim) {
            initialArray[ioff+iSpace+iDim] += orientation[jDim*spaceDim+iDim] * tmpLocal[jDim];
	  } // for
        } // for
      } // if

      if (_dbRate) {
        // Rotate traction vector from local coordinate system to global coordinate system
	assert(rateVisitor);
	const PetscInt roff = rateVisitor->sectionOffset(c);
	assert(numQuadPts*spaceDim == rateVisitor->sectionDof(c));
        for(int iDim = 0; iDim < spaceDim; ++iDim) {
          tmpLocal[iDim] = rateArray[roff+iSpace+iDim];
        } // for
        for(int iDim = 0; iDim < spaceDim; ++iDim) {
          rateArray[roff+iSpace+iDim] = 0.0;
          for(int jDim = 0; jDim < spaceDim; ++jDim) {
            rateArray[roff+iSpace+iDim] += orientation[jDim*spaceDim+iDim] * tmpLocal[jDim];
	  } // for
        } // for
      } // if

      if (_dbChange) {
        // Rotate traction vector from local coordinate system to global coordinate system
	assert(changeVisitor);
	const PetscInt coff = changeVisitor->sectionOffset(c);
	assert(numQuadPts*spaceDim == changeVisitor->sectionDof(c));
        for (int iDim = 0; iDim < spaceDim; ++iDim) {
          tmpLocal[iDim] = changeArray[coff+iSpace+iDim];
        } // for
        for(int iDim = 0; iDim < spaceDim; ++iDim) {
          changeArray[coff+iSpace+iDim] = 0.0;
          for(int jDim = 0; jDim < spaceDim; ++jDim) {
            changeArray[coff+iSpace+iDim] += orientation[jDim*spaceDim+iDim] * tmpLocal[jDim];
	  } // for
        } // for
      } // if
    } // for
  } // for

  delete initialVisitor; initialVisitor = 0;
  delete rateVisitor; rateVisitor = 0;
  delete changeVisitor; changeVisitor = 0;

  PYLITH_METHOD_END;
} // paramsLocalToGlobal

// ----------------------------------------------------------------------
// Calculate temporal and spatial variation of value over the list of Submesh.
void
pylith::bc::Neumann::_calculateValue(const PylithScalar t)
{ // _calculateValue
  PYLITH_METHOD_BEGIN;

  assert(_parameters);
  assert(_boundaryMesh);
  assert(_quadrature);

  const PylithScalar timeScale = _getNormalizer().timeScale();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::Stratum cellsStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  const int spaceDim = _quadrature->spaceDim();
  const int numQuadPts = _quadrature->numQuadPts();

  // Get sections
  topology::Field& valueField = _parameters->get("value");
  topology::VecVisitorMesh valueVisitor(valueField);
  PetscScalar* valueArray = valueVisitor.localArray();
  
  topology::Field* initialField = (_dbInitial) ? &_parameters->get("initial") : 0;
  topology::VecVisitorMesh* initialVisitor = (initialField) ? new topology::VecVisitorMesh(*initialField) : 0;
  PetscScalar* initialArray = (initialVisitor) ? initialVisitor->localArray() : NULL;

  topology::Field* rateField = (_dbRate) ? &_parameters->get("rate") : 0;
  topology::VecVisitorMesh* rateVisitor = (rateField) ? new topology::VecVisitorMesh(*rateField) : 0;
  PetscScalar* rateArray = (rateVisitor) ? rateVisitor->localArray() : NULL;

  topology::Field* rateTimeField = (_dbRate) ? &_parameters->get("rate time") : 0;
  topology::VecVisitorMesh* rateTimeVisitor = (rateTimeField) ? new topology::VecVisitorMesh(*rateTimeField) : 0;
  PetscScalar* rateTimeArray = (rateTimeVisitor) ? rateTimeVisitor->localArray() : NULL;

  topology::Field* changeField = (_dbChange) ? &_parameters->get("change") : 0;
  topology::VecVisitorMesh* changeVisitor = (changeField) ? new topology::VecVisitorMesh(*changeField) : 0;
  PetscScalar* changeArray = (changeVisitor) ? changeVisitor->localArray() : NULL;

  topology::Field* changeTimeField = (_dbChange) ? &_parameters->get("change time") : 0;
  topology::VecVisitorMesh* changeTimeVisitor = (changeTimeField) ? new topology::VecVisitorMesh(*changeTimeField) : 0;
  PetscScalar* changeTimeArray = (changeTimeVisitor) ? changeTimeVisitor->localArray() : NULL;

  for(PetscInt c = cStart; c < cEnd; ++c) {
    const PetscInt voff = valueVisitor.sectionOffset(c);
    const PetscInt vdof = valueVisitor.sectionDof(c);
    assert(numQuadPts*spaceDim == vdof);
    for (PetscInt d = 0; d < vdof; ++d) {
      valueArray[voff+d] = 0.0;
    } // for

    // Contribution from initial value
    if (_dbInitial) {
      assert(initialVisitor);
      const PetscInt ioff = initialVisitor->sectionOffset(c);
      const PetscInt idof = initialVisitor->sectionDof(c);
      assert(numQuadPts*spaceDim == idof);
      for (PetscInt d = 0; d < idof; ++d) {
        valueArray[voff+d] += initialArray[ioff+d];
      } // for
    } // if
    
    // Contribution from rate of change of value
    if (_dbRate) {
      assert(rateVisitor);
      const PetscInt roff = rateVisitor->sectionOffset(c);
      const PetscInt rdof = rateVisitor->sectionDof(c);
      assert(numQuadPts*spaceDim == rdof);
      assert(rateTimeVisitor);
      const PetscInt rtoff = rateTimeVisitor->sectionOffset(c);
      assert(numQuadPts == rateTimeVisitor->sectionDof(c));

      for(int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
        const PylithScalar tRel = t - rateTimeArray[rtoff+iQuad];
        if (tRel > 0.0)  // rate of change integrated over time
          for(int iDim = 0; iDim < spaceDim; ++iDim) {
            valueArray[voff+iQuad*spaceDim+iDim] += rateArray[roff+iQuad*spaceDim+iDim] * tRel;
	  } // for
      } // for
    } // if
    
    // Contribution from change of value
    if (_dbChange) {
      assert(changeVisitor);
      const PetscInt coff = changeVisitor->sectionOffset(c);
      assert(changeTimeField);
      const PetscInt ctoff = changeTimeVisitor->sectionOffset(c);
      assert(numQuadPts == changeTimeVisitor->sectionDof(c));

      for(int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
        const PylithScalar tRel = t - changeTimeArray[ctoff+iQuad];
        if (tRel >= 0) { // change in value over time
          PylithScalar scale = 1.0;
          if (_dbTimeHistory) {
            PylithScalar tDim = tRel;
            _getNormalizer().dimensionalize(&tDim, 1, timeScale);
            const int err = _dbTimeHistory->query(&scale, tDim);
            if (err) {
              std::ostringstream msg;
              msg << "Error querying for time '" << tDim 
                  << "' in time history database "
                  << _dbTimeHistory->label() << ".";
              throw std::runtime_error(msg.str());
            } // if
          } // if
          for (int iDim = 0; iDim < spaceDim; ++iDim) {
            valueArray[voff+iQuad*spaceDim+iDim] += changeArray[coff+iQuad*spaceDim+iDim]*scale;
	  } // for
        } // if
      } // for
    } // if
  } // for

  delete initialVisitor; initialVisitor = 0;
  delete rateVisitor; rateVisitor = 0;
  delete rateTimeVisitor; rateTimeVisitor = 0;
  delete changeVisitor; changeVisitor = 0;
  delete changeTimeVisitor; changeTimeVisitor = 0;

  PYLITH_METHOD_END;
}  // _calculateValue


// End of file 
