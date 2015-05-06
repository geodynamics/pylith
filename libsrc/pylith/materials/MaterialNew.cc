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
pylith::materials::MaterialNew::MaterialNew(const int dimension) :
  _materialIS(0),
  _dimension(dimension),
  _dbAuxFields(0),
  _id(0),
  _label("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::MaterialNew::~MaterialNew(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::MaterialNew::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  delete _materialIS; _materialIS = 0;

  _dbAuxFields = 0; // :TODO: Use shared pointer.

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Get physical property parameters and initial state (if used) from database.
void
pylith::materials::MaterialNew::initialize(const topology::Mesh& mesh)
{ // initialize
  PYLITH_METHOD_BEGIN;

  // Get cells associated with material
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  const bool includeOnlyCells = true;
  delete _materialIS; _materialIS = new topology::StratumIS(dmMesh, "material-id", _id, includeOnlyCells);assert(_materialIS);
  const PetscInt numCells = _materialIS->size();
  const PetscInt* cells = _materialIS->points();

  delete _auxFields; _auxFields = new topology::Field(mesh);assert(_auxFields);
  _auxFields->label("auxiliary fields");

  // :TODO: Update setup of auxiliary field
#if 0
  int_array cellsTmp(cells, numCells);
  _auxFields->newSection(cellsTmp, propsFiberDim);
  _auxFields->allocate();
  _auxFields->zeroAll();
#endif

  if (_dbAuxFields) {
    _initializeAuxFieldsFromDB();
  } else { // else
    assert(0);
    throw std::logic_error("Unknown case for setting up auxiliary fields.");
  } // if/else

  PYLITH_METHOD_END;
} // initialize


// ----------------------------------------------------------------------
// Initialize auxiliary fields using spatial database.
void
pylith::materials::MaterialNew::_initializeAuxFieldsFromDB(void)
{ // _initializeAuxFields
  PYLITH_METHOD_BEGIN;

  assert(_dbAuxFields);
  
  /* :TODO: Redo this function to look like:

  void (*matFuncs[1])(const PetscReal x[], PetscScalar *u, void *ctx) = {nu_2d};
  Vec            nu;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMCreateLocalVector(dmAux, &nu);CHKERRQ(ierr);
  ierr = DMPlexProjectFunctionLocal(dmAux, matFuncs, NULL, INSERT_ALL_VALUES, nu);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject) dm, "A", (PetscObject) nu);CHKERRQ(ierr);
  ierr = VecDestroy(&nu);CHKERRQ(ierr);

void nu_2d(const PetscReal x[], PetscScalar *u, void *ctx)
{
  *u = x[0] + x[1];
}

  */
  
#if 0 // OBSOLETE
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();assert(cs);

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

#endif


  PYLITH_METHOD_END;
} // _initializeAuxFields


// End of file 
