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
pylith::topology::FieldQuery::FieldQuery(void)
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

  _queryFns.clear();

  PYLITH_METHOD_END;
} // deallocate

 
// ----------------------------------------------------------------------
// Set query function for subfield.
void
pylith::topology::FieldQuery::setQuery(const char* subfield,
				       const queryfn_type fn)
{ // setQuery
  PYLITH_METHOD_BEGIN;

  assert(subfield);
  assert(fn);

  _queryFns[subfield] = fn;

  PYLITH_METHOD_END;
} // setQuery


// ----------------------------------------------------------------------
// Query spatial database to set values in field.
void
pylith::topology::FieldQuery::queryDB(Field* field,
				      spatialdata::spatialdb::SpatialDB* db,
				      const PylithReal lengthScale)
{ // queryDB
  PYLITH_METHOD_BEGIN;
  
  assert(field);
  assert(db);

  // Create contexts and funcs. Need to put contexts into an array of
  // pointers, since Petsc function doesn't know the size of the
  // context.
  const Field::subfields_type& subfields = field->_subfields;
  const unsigned size = subfields.size();
  assert(_queryFns.size() == size);
  queryfn_type* fns = (size > 0) ? new queryfn_type[size] : NULL;
  DBQueryContext* contexts = (size > 0) ? new DBQueryContext[size] : NULL;
  DBQueryContext** contextPtrs = (size > 0) ? new DBQueryContext*[size] : NULL;

  int i=0;
  for (Field::subfields_type::const_iterator iter=subfields.begin(); iter != subfields.end(); ++iter, ++i) {
    const std::string& name = iter->first;
    if (_queryFns.find(name) == _queryFns.end()) { // if
      std::ostringstream msg;
      msg << "FieldQuery for field '" << field->label() << "' missing query function for subfield '" << name << "'";
      throw std::logic_error(msg.str());
    } // if/else
    fns[i] = _queryFns[name];
    contexts[i].db = db;
    contexts[i].cs = field->mesh().coordsys();
    contexts[i].lengthScale = lengthScale;
    contexts[i].valueScale = iter->second.metadata.scale;
    contextPtrs[i] = &contexts[i];
  } // for
  
    // Open spatial database.
  db->open();
  
  PetscErrorCode err = 0;
  const PetscDM dm = field->dmMesh();
  const PetscVec fieldVec = field->localVector();
  err = DMPlexProjectFunctionLocal(dm, fns, (void**)contextPtrs, INSERT_ALL_VALUES, fieldVec);PYLITH_CHECK_ERROR(err);
  //err = PetscObjectCompose((PetscObject) dm, "A", (PetscObject) fieldVec);CHKERRQ(ierr); // :MATT: Which dm is this? Do we need this?
  
  delete[] fns; fns = NULL;
  delete[] contexts; contexts = NULL;
  delete[] contextPtrs; contextPtrs = NULL;
  
  // Close spatial database.
  db->close();

  PYLITH_METHOD_END;
} // queryDB


// End of file 
