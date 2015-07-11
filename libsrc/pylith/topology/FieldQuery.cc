// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "FieldQuery.hh" // implementation of class methods

#include "Field.hh" // USES Field
#include "Mesh.hh" // USES Mesh

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::FieldQuery::FieldQuery(const Field& field) :
  _field(field),
  _functions(NULL),
  _contexts(NULL),
  _contextPtrs(NULL)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::FieldQuery::~FieldQuery(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::FieldQuery::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  delete[] _functions; _functions = NULL;
  delete[] _contexts; _contexts = NULL;
  delete[] _contextPtrs; _contextPtrs = NULL;
  
  _queryFns.clear();

  PYLITH_METHOD_END;
} // deallocate

 
// ----------------------------------------------------------------------
// Set query function for subfield.
void
pylith::topology::FieldQuery::queryFn(const char* subfield,
				      const queryfn_type fn)
{ // setQuery
  PYLITH_METHOD_BEGIN;

  assert(subfield);
  assert(fn);

  _queryFns[subfield] = fn;

  PYLITH_METHOD_END;
} // setQuery


// ----------------------------------------------------------------------
// Set query function for subfield.
const pylith::topology::FieldQuery::queryfn_type
pylith::topology::FieldQuery::queryFn(const char* subfield) const
{ // setQuery
  PYLITH_METHOD_BEGIN;

  assert(subfield);
  const queryfn_map_type::const_iterator& iter = _queryFns.find(subfield);
  if (iter == _queryFns.end()) {
    std::ostringstream msg;
    msg << "Could not find query function for subfield '" << subfield << "'." << std::endl;
    throw std::logic_error(msg.str());
  } // if

  PYLITH_METHOD_RETURN(*iter->second);
} // setQuery


// ----------------------------------------------------------------------
// Get array of query functions.
pylith::topology::FieldQuery::queryfn_type*
pylith::topology::FieldQuery::functions(void) const
{ // functions
  return _functions;
} // functions
 

// ----------------------------------------------------------------------
// Get array of pointers to contexts.
const pylith::topology::FieldQuery::DBQueryContext* const*
pylith::topology::FieldQuery::contextPtrs(void) const
{ // contextPtrs
  return _contextPtrs;
} // contextPtrs


// ----------------------------------------------------------------------
// Query spatial database to set values in field.
void
pylith::topology::FieldQuery::openDB(spatialdata::spatialdb::SpatialDB* db,
				     const PylithReal lengthScale)
{ // openDB
  PYLITH_METHOD_BEGIN;
  
  assert(db);

  // Create contexts and funcs. Need to put contexts into an array of
  // pointers, since Petsc function doesn't know the size of the
  // context.
  const Field::subfields_type& subfields = _field._subfields;
  const unsigned size = subfields.size();
  assert(_queryFns.size() == size);
  delete[] _functions; _functions = (size > 0) ? new queryfn_type[size] : NULL;
  delete[] _contexts; _contexts = (size > 0) ? new DBQueryContext[size] : NULL;
  delete[] _contextPtrs; _contextPtrs = (size > 0) ? new DBQueryContext*[size] : NULL;

  int i=0;
  for (Field::subfields_type::const_iterator iter=subfields.begin(); iter != subfields.end(); ++iter, ++i) {
    const std::string& name = iter->first;
    if (_queryFns.find(name) == _queryFns.end()) { // if
      std::ostringstream msg;
      msg << "FieldQuery for field '" << _field.label() << "' missing query function for subfield '" << name << "'";
      throw std::logic_error(msg.str());
    } // if/else
    _functions[i] = _queryFns[name];

    _contexts[i].db = db;
    _contexts[i].cs = _field.mesh().coordsys();
    _contexts[i].lengthScale = lengthScale;
    _contexts[i].valueScale = iter->second.metadata.scale;

    _contextPtrs[i] = &_contexts[i];
  } // for
  
    // Open spatial database.
  db->open();
  
  PYLITH_METHOD_END;
} // openDB


// ----------------------------------------------------------------------
// Query spatial database to set values in field.
void
pylith::topology::FieldQuery::queryDB(void)
{ // queryDB
  PYLITH_METHOD_BEGIN;
  
  PetscErrorCode err = 0;
  const PetscDM dm = _field.dmMesh();
  const PetscVec fieldVec = _field.localVector();
  err = DMPlexProjectFunctionLocal(dm, _functions, (void**)_contextPtrs, INSERT_ALL_VALUES, fieldVec);PYLITH_CHECK_ERROR(err);
  //err = PetscObjectCompose((PetscObject) dm, "A", (PetscObject) fieldVec);CHKERRQ(ierr); // :MATT: Which dm is this? Do we need this?
  
  PYLITH_METHOD_END;
} // queryDB


// ----------------------------------------------------------------------
// Query spatial database to set values in field.
void
pylith::topology::FieldQuery::closeDB(spatialdata::spatialdb::SpatialDB* db)
{ // queryDB
  PYLITH_METHOD_BEGIN;
  
  assert(db);

  delete[] _functions; _functions = NULL;
  delete[] _contexts; _contexts = NULL;
  delete[] _contextPtrs; _contextPtrs = NULL;
  
  // Close spatial database.
  db->close();

  PYLITH_METHOD_END;
} // queryDB


// End of file 
