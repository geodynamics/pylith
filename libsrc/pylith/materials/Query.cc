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

#include <iostream>
// ----------------------------------------------------------------------
// Query spatial database for density in 2-D.
void
pylith::materials::Query::dbQueryDensity2D(const PetscReal x[],
					   PetscScalar* values,
					   void* context)
{ // dbQueryDensity2D
  PYLITH_METHOD_BEGIN;

  const int _dim = 2;
  const int _numComps = 1;

  const int dim = _dim; // :KLUDGE: TEMPORARY hardwiring size
  const int numComps = _numComps; // :KLUDGE: TEMPORARY hardwiring size

  assert(x);
  assert(values);
  assert(context);
  assert(_numComps == numComps);
  assert(_dim == dim);

  const pylith::topology::FieldQuery::DBQueryContext* queryctx = (pylith::topology::FieldQuery::DBQueryContext*)context;assert(queryctx);

  // Tell database which values we want.
  const int numDBValues = 1;
  PylithReal dbValues[numDBValues];
  const int i_density = 0;
  const char* dbValueNames[numDBValues] = {"density"};
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
    msg << "Could not find density at (";
    for (int i=0; i < _dim; ++i)
      msg << "  " << xDim[i];
    msg << ") in using spatial database '" << queryctx->db->label() << "'.";
    std::cout << "MSG: "<< msg.str() << std::endl;
    throw std::runtime_error(msg.str());
  } // if
  
  const PylithReal density = dbValues[i_density];

  if (density <= 0) {
    std::ostringstream msg;
    msg << "Found negative density (" << density << ") at location (";
    for (int i=0; i < _dim; ++i)
      msg << "  " << xDim[i];
    msg << ") in spatial database '" << queryctx->db->label() << "'.";
    throw std::runtime_error(msg.str());
  } // if

  assert(queryctx->valueScale > 0);
  values[0] = density / queryctx->valueScale;

  PYLITH_METHOD_END;
} // dbQueryDensity2D


// ----------------------------------------------------------------------
// Query spatial database for shear modulus, $\mu$, in 2-D.
void
pylith::materials::Query::dbQueryMu2D(const PetscReal x[],
				      PetscScalar* values,
				      void* context)
{ // dbQueryMu2D
  PYLITH_METHOD_BEGIN;

  const int _dim = 2;
  const int _numComps = 1;

  const int dim = _dim; // :KLUDGE: TEMPORARY hardwiring size
  const int numComps = _numComps; // :KLUDGE: TEMPORARY hardwiring size

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
    msg << "Could not find density and Vs at (";
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
} // dbQueryMu2D


// ----------------------------------------------------------------------
// Query spatial database for Lame's constant, $\lambda$, in 2-D.
void
pylith::materials::Query::dbQueryLambda2D(const PetscReal x[],
					  PetscScalar* values,
					  void* context)
{ // dbQueryLambda2D
  PYLITH_METHOD_BEGIN;

  const int _dim = 2;
  const int _numComps = 1;

  const int dim = _dim; // :KLUDGE: TEMPORARY hardwiring size
  const int numComps = _numComps; // :KLUDGE: TEMPORARY hardwiring size

  assert(x);
  assert(values);
  assert(context);
  assert(_numComps == numComps);
  assert(_dim == dim);

  const pylith::topology::FieldQuery::DBQueryContext* queryctx = (pylith::topology::FieldQuery::DBQueryContext*)context;assert(queryctx);

  // Tell database which values we want.
  const int numDBValues = 3;
  PylithReal dbValues[numDBValues];
  const int i_density = 0;
  const int i_vs = 1;
  const int i_vp = 2;
  const char* dbValueNames[numDBValues] = {"density", "vs", "vp"};
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
    msg << "Could not find density, Vs, and Vp at (";
    for (int i=0; i < _dim; ++i)
      msg << "  " << xDim[i];
    msg << ") in using spatial database '" << queryctx->db->label() << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  const PylithReal density = dbValues[i_density];
  const PylithReal vs = dbValues[i_vs];
  const PylithReal vp = dbValues[i_vp];

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
  if (vp <= 0) {
    std::ostringstream msg;
    msg << "Found negative dilatational wave speed (" << vp << ") at location (";
    for (int i=0; i < _dim; ++i)
      msg << "  " << xDim[i];
    msg << ") in spatial database '" << queryctx->db->label() << "'.";
    throw std::runtime_error(msg.str());
  } // if

  const PylithReal mu = density * vs * vs;assert(mu > 0);
  const PylithReal lambda = density * vs * vs - 2.0*mu;
  assert(queryctx->valueScale > 0);
  values[0] = lambda / queryctx->valueScale;

  PYLITH_METHOD_END;
} // dbQueryLambda2D


// ----------------------------------------------------------------------
// Query spatial database for body force vector in 2-D.
void
pylith::materials::Query::dbQueryBodyForce2D(const PetscReal x[],
					     PetscScalar* values,
					     void* context)
{ // dbQueryMu
  PYLITH_METHOD_BEGIN;

  const int _dim = 2;
  const int _numComps = _dim;

  const int dim = _dim; // :KLUDGE: TEMPORARY hardwiring size
  const int numComps = _numComps; // :KLUDGE: TEMPORARY hardwiring size

  assert(x);
  assert(values);
  assert(context);
  assert(_numComps == numComps);
  assert(_dim == dim);

  const pylith::topology::FieldQuery::DBQueryContext* queryctx = (pylith::topology::FieldQuery::DBQueryContext*)context;assert(queryctx);

  // Tell database which values we want.
  const int numDBValues = 2;
  PylithReal dbValues[numDBValues];
  const int i_fx = 0;
  const int i_fy = 1;
  const char* dbValueNames[numDBValues] = {"body-force-x", "body-force-y"};
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
    msg << "Could not find body-force-x and body-force-y at (";
    for (int i=0; i < _dim; ++i)
      msg << "  " << xDim[i];
    msg << ") in using spatial database '" << queryctx->db->label() << "'.";
    throw std::runtime_error(msg.str());
  } // if
  
  const PylithReal fx = dbValues[i_fx];
  const PylithReal fy = dbValues[i_fy];

  assert(queryctx->valueScale > 0);
  values[0] = fx / queryctx->valueScale;
  values[1] = fy / queryctx->valueScale;

  PYLITH_METHOD_END;
} // dbQueryBodyForce2D


// End of file
