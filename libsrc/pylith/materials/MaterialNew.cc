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

#include "MaterialNew.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Stratum.hh" // USES StratumIS
#include "pylith/utils/array.hh" // USES scalar_array, std::vector

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <petscds.h> // USES PetscDS

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::Material::Material(const int dimension) :
  _auxFields(0),
  _normalizer(new spatialdata::units::Nondimensional),
  _materialIS(0),
  _dimension(dimension),
  _needNewJacobian(false),
  _isJacobianSymmetric(true),
  _dbAuxFields(0),
  _id(0),
  _label("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::Material::~Material(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::Material::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  delete _normalizer; _normalizer = 0;
  delete _materialIS; _materialIS = 0;
  delete _properties; _properties = 0;
  delete _stateVars; _stateVars = 0;

  _dbAuxFields = 0; // :TODO: Use shared pointer.

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Set scales used to nondimensionalize physical properties.
void
pylith::materials::Material::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
  PYLITH_METHOD_BEGIN;

  if (!_normalizer)
    _normalizer = new spatialdata::units::Nondimensional(dim);
  else
    *_normalizer = dim;

  PYLITH_METHOD_END;
} // normalizer

// ----------------------------------------------------------------------
// Get physical property parameters and initial state (if used) from database.
void
pylith::materials::Material::initialize(const topology::Mesh& mesh)
{ // initialize
  PYLITH_METHOD_BEGIN;

  // Get quadrature information
#if 0
  const int numQuadPts = quadrature->numQuadPts();
  const int numBasis = quadrature->numBasis();
  const int numCorners = quadrature->refGeometry().numCorners();
  const int spaceDim = quadrature->spaceDim();
#endif

  // Get cells associated with material
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  const bool includeOnlyCells = true;
  delete _materialIS; _materialIS = new topology::StratumIS(dmMesh, "material-id", _id, includeOnlyCells);assert(_materialIS);
  const PetscInt numCells = _materialIS->size();
  const PetscInt* cells = _materialIS->points();

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();assert(cs);

  delete _auxFields; _auxFields = new topology::Field(mesh);assert(_auxFields);
  _auxFields->label("auxiliary fields");
  const int fieldsSize = 0;
  int_array cellsTmp(cells, numCells);

  _auxFields->newSection(cellsTmp, propsFiberDim);
  _auxFields->allocate();
  _auxFields->zeroAll();

  // TODO Need to decide how to manage PetscDS
  PetscDS        prob;
  PetscFE        fe;
  PetscInt       dim, closureSize, vStart, vEnd, numVertices = 0;
  PetscInt      *closure   = NULL, c;
  PetscBool      isSimplex = PETSC_FALSE;
  PetscErrorCode err;

  err = DMGetDimension(dmMesh, &dim);PYLITH_CHECK_ERROR(err);

  // BEGIN REFACTOR THIS
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetTransitiveClosure(dmMesh, 0, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
  for (c = 0; c < closureSize*2; ++c) if ((closure[c] >= vStart) && (closure[c] < vEnd)) ++numVertices;
  if (numVertices == dim+1) isSimplex = PETSC_TRUE;
  err = DMPlexRestoreTransitiveClosure(dmMesh, 0, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
  // END REFACTOR THIS

  // :MATT: Need some comments here to explain what this is doing. 
  // :MATT: Why do we need fiberDim?
  err = PetscFECreateDefault(dmMesh, dim, fiberDim, isSimplex, NULL, -1, &fe);PYLITH_CHECK_ERROR(err);
  err = PetscDSCreate(mesh.comm(), &prob);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetDiscretization(prob, 0, (PetscObject) fe);PYLITH_CHECK_ERROR(err);
  err = PetscFEDestroy(&fe);PYLITH_CHECK_ERROR(err);
  err = PetscDSSetUp(prob);PYLITH_CHECK_ERROR(err);
  err = DMSetDS(_auxFields->dmMesh(), prob);PYLITH_CHECK_ERROR(err);
  err = PetscDSDestroy(&prob);PYLITH_CHECK_ERROR(err);

  if (_dbAuxFields) {
    _inittializeAuxFieldsDB(mesh);
  } else { // else
    assert(0);
    throw std::logic_error("Unknown case for setting up auxiliary fields.");
  } // if/else

  PYLITH_METHOD_END;
} // initialize


// ----------------------------------------------------------------------
// Check whether material has a given auxilirary field.
bool
pylith::materials::Material::hasAuxField(const char* name)
{ // hasAuxField
  PYLITH_METHOD_BEGIN;

  assert(_auxFields);

  PYLITH_METHOD_RETURN(_auxFields->hasSubfield(name));
} // hasAuxField


// ----------------------------------------------------------------------
// Get auxiliary field.
void
pylith::materials::Material::getAuxField(topology::Field *field,
					 const char* name) const
{ // getAuxField
  PYLITH_METHOD_BEGIN;

  assert(field);
  assert(_auxFields);

  _auxFields->copySubfield(*field, name);

  PYLITH_METHOD_END;
} // getAuxField
  
// ----------------------------------------------------------------------
// Initialize auxiliary fields using spatial database.
void
pylith::materials::_initializeAuxFieldsDB(void)
{ // _initializeAuxFields
  PYLITH_METHOD_BEGIN;

  assert(_dbAuxFields);

  topology::VecVisitorMesh auxFieldsVisitor(*_auxFields);
  PetscScalar* auxFieldsArray = auxFieldsVisitor.localArray();

  scalar_array coordsCell(numCorners*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmMesh);

  // Optimize coordinate retrieval in closure  
  topology::CoordsVisitor::optimizeClosure(dmMesh);

  // Create arrays for querying.
  const int numDBProperties = _metadata.numDBProperties();
  scalar_array coordsQuadPtsGlobal(numQuadPts*spaceDim);
  scalar_array valuesQuery(numDBProperties);
  scalar_array dbValuesCell(propsFiberDim);

  // Setup database for quering for physical properties
  assert(_dbAuxFields);
  _dbAuxFields->open();
  _dbAuxFields->queryVals(_metadata.dbProperties(), _metadata.numDBProperties());

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();

  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];

    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, cell);
#if 0
    quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), cell);

    const scalar_array& coordsQuadPtsNonDim = quadrature->quadPts();
    coordsQuadPtsGlobal = coordsQuadPtsNonDim;
#endif
    _normalizer->dimensionalize(&coordsQuadPtsGlobal[0], coordsQuadPtsGlobal.size(), lengthScale);

    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, index=0; iQuadPt < numQuadPts; ++iQuadPt, index+=spaceDim) {
      int err = _dbAuxFields->query(&dbValues[0], numDBValues, &coordsQuadPtsGlobal[index], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find parameters for auxilirary fieldsat " << "(";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPtsGlobal[index+i];
	msg << ") in material '" << label() << "' using spatial database '" << _dbAuxFields->label() << "'.";
	throw std::runtime_error(msg.str());
      } // if
      _dbToAuxFields(&auxFieldsCell[iQuadPt*numAuxValuesQuadPt], numAuxValuesQuadPt, dbValues);
      _nondimProperties(&auxFieldsCell[iQuadPt*_numAuxValuesQuadPt], _numAuxValuesQuadPt);

    } // for
    // Insert cell contribution into fields
    const PetscInt off = auxFieldsVisitor.sectionOffset(cell);
    assert(auxFieldsFiberDim == auxFieldsVisitor.sectionDof(cell));
    for(PetscInt d = 0; d < auxFieldsFiberDim; ++d) {
      auxFieldsArray[off+d] = auxFieldsCell[d];
    } // for
  } // for
  delete stateVarsVisitor; stateVarsVisitor = 0;

  // Close databases
  _dbAuxFields->close();


  PYLITH_METHOD_END;
} // _initializeAuxFields


// End of file 
