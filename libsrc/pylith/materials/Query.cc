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

#include "pylith/materials/Query.hh"

#include "pylith/topology/FieldQuery.hh" // USES DBQueryContext
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
// Query spatial database for shear modulus, $\mu$, in 2-D.
void
pylith::materials::Query::dbQueryMu2D(const PetscReal x[],
				      PetscScalar* values,
				      void* context)
{ // dbQueryMu
  PYLITH_METHOD_BEGIN;

  const int _dim = 2;
  const int _numComps = 1;

  const int dim = _dim; // :KLUDGE: hardwiring size
  const int numComps = 1; // :KLUDGE: hardwiring size

  assert(x);
  assert(values);
  assert(context);
  assert(_numComps == numComps);
  assert(_dim == dim);

  const pylith::topology::FieldQuery::DBQueryContext* queryctx = (pylith::topology::FieldQuery::DBQueryContext*)context;assert(queryctx);

  // Tell database which values we want.
  const int numDBValues = 2;
  PylithReal dbValues[numDBValues];
  const int i_density = 0;
  const int i_vs = 1;
  const char* dbValueNames[numDBValues] = {"density", "vs"};
  queryctx->db->queryVals(dbValueNames, numDBValues);

  // Dimensionalize query location coordinates.
  assert(queryctx->lengthScale > 0);
  double xDim[_dim];
  for (int i=0; i < _dim; ++i) {
    xDim[i] = x[i] * queryctx->lengthScale;
  } // for

  assert(queryctx->cs);
  const int err = queryctx->db->query(dbValues, numDBValues, xDim, _dim, queryctx->cs);
  if (err) {
    std::ostringstream msg;
    msg << "Could not find density and vs at (";
    for (int i=0; i < _dim; ++i)
      msg << "  " << xDim[i];
    msg << ") in using spatial database '" << queryctx->db->label() << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  const PylithReal density = dbValues[i_density];
  const PylithReal vs = dbValues[i_vs];

  if (density <= 0) {
    std::ostringstream msg;
    msg << "Found negative density (" << density << ") at location (";
    for (int i=0; i < _dim; ++i)
      msg << "  " << xDim[i];
    msg << ") in spatial database '" << queryctx->db->label() << "'.";
    throw std::runtime_error(msg.str());
  } // if
  if (vs <= 0) {
    std::ostringstream msg;
    msg << "Found negative shear wave speed (" << vs << ") at location (";
    for (int i=0; i < _dim; ++i)
      msg << "  " << xDim[i];
    msg << ") in spatial database '" << queryctx->db->label() << "'.";
    throw std::runtime_error(msg.str());
  } // if

  const PylithReal mu = density * vs * vs;assert(mu > 0);
  assert(queryctx->valueScale > 0);
  values[0] = mu / queryctx->valueScale;

  PYLITH_METHOD_END;
} // dbQueryMu


// End of file
