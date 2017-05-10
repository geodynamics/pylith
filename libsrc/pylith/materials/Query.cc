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
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END, PYLITH_SET_ERROR

// ----------------------------------------------------------------------
// Query for Vs, Vp, and density to determine shear modulus.
PetscErrorCode
pylith::materials::Query::dbQueryShearModulus(PylithInt dim,
                                              PylithReal t,
                                              const PylithReal x[],
                                              PylithInt nvalues,
                                              PylithScalar* values,
                                              void* context)
{ // dbQueryShearModulus
    PYLITH_METHOD_BEGIN;

    const int _nvalues = 1;

    assert(x);
    assert(values);
    assert(context);
    assert(_nvalues == nvalues);
    assert(2 == dim || 3 == dim);

    const pylith::topology::FieldQuery::DBQueryContext* queryctx = (pylith::topology::FieldQuery::DBQueryContext*)context; assert(queryctx);

    // Tell database which values we want.
    const int numDBValues = 2;
    PylithReal dbValues[numDBValues];
    const int i_density = 0;
    const int i_vs = 1;
    const char* dbValueNames[numDBValues] = {"density", "vs"};
    try {
        queryctx->db->queryVals(dbValueNames, numDBValues);
    } catch (const std::exception& err) {
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, err.what());
    } // try/catch


    // Dimensionalize query location coordinates.
    assert(queryctx->lengthScale > 0);
    double xDim[3];
    for (int i=0; i < dim; ++i) {
        xDim[i] = x[i] * queryctx->lengthScale;
    } // for

    assert(queryctx->cs);
    const int err = queryctx->db->query(dbValues, numDBValues, xDim, dim, queryctx->cs);
    if (err) {
        std::ostringstream msg;
        msg << "Could not find density and Vs at (";
        for (int i=0; i < dim; ++i)
            msg << "  " << xDim[i];
        msg << ") using spatial database '" << queryctx->db->label() << "'.";
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
    } // if

    const PylithReal density = dbValues[i_density];
    const PylithReal vs = dbValues[i_vs];

    if (density <= 0) {
        std::ostringstream msg;
        msg << "Found negative density (" << density << ") at location (";
        for (int i=0; i < dim; ++i)
            msg << "  " << xDim[i];
        msg << ") in spatial database '" << queryctx->db->label() << "'.";
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
    } // if
    if (vs <= 0) {
        std::ostringstream msg;
        msg << "Found negative shear wave speed (" << vs << ") at location (";
        for (int i=0; i < dim; ++i)
            msg << "  " << xDim[i];
        msg << ") in spatial database '" << queryctx->db->label() << "'.";
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
    } // if

    const PylithReal shearModulus = density * vs * vs; assert(shearModulus > 0);
    assert(queryctx->valueScale > 0);
    values[0] = shearModulus / queryctx->valueScale;

    PYLITH_METHOD_RETURN(0);
} // dbQueryShearModulus


// ----------------------------------------------------------------------
// Query for Vs, Vp, and density to determine bulk modulus.
PetscErrorCode
pylith::materials::Query::dbQueryBulkModulus(PylithInt dim,
                                             PylithReal t,
                                             const PylithReal x[],
                                             PylithInt nvalues,
                                             PylithScalar* values,
                                             void* context)
{ // dbQueryBulkModulus
    PYLITH_METHOD_BEGIN;

    const int _nvalues = 1;

    assert(x);
    assert(values);
    assert(context);
    assert(_nvalues == nvalues);
    assert(2 == dim || 3 == dim);

    const pylith::topology::FieldQuery::DBQueryContext* queryctx = (pylith::topology::FieldQuery::DBQueryContext*)context; assert(queryctx);

    // Tell database which values we want.
    const int numDBValues = 3;
    PylithReal dbValues[numDBValues];
    const int i_density = 0;
    const int i_vs = 1;
    const int i_vp = 2;
    const char* dbValueNames[numDBValues] = {"density", "vs", "vp"};
    try {
        queryctx->db->queryVals(dbValueNames, numDBValues);
    } catch (const std::exception& err) {
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, err.what());
    } // try/catch


    // Dimensionalize query location coordinates.
    assert(queryctx->lengthScale > 0);
    double xDim[3];
    for (int i=0; i < dim; ++i) {
        xDim[i] = x[i] * queryctx->lengthScale;
    } // for

    assert(queryctx->cs);
    const int err = queryctx->db->query(dbValues, numDBValues, xDim, dim, queryctx->cs);
    if (err) {
        std::ostringstream msg;
        msg << "Could not find density, Vs, and Vp at (";
        for (int i=0; i < dim; ++i)
            msg << "  " << xDim[i];
        msg << ") using spatial database '" << queryctx->db->label() << "'.";
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
    } // if

    const PylithReal density = dbValues[i_density];
    const PylithReal vs = dbValues[i_vs];
    const PylithReal vp = dbValues[i_vp];

    if (density <= 0) {
        std::ostringstream msg;
        msg << "Found negative density (" << density << ") at location (";
        for (int i=0; i < dim; ++i)
            msg << "  " << xDim[i];
        msg << ") in spatial database '" << queryctx->db->label() << "'.";
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
    } // if
    if (vs <= 0) {
        std::ostringstream msg;
        msg << "Found negative shear wave speed (" << vs << ") at location (";
        for (int i=0; i < dim; ++i)
            msg << "  " << xDim[i];
        msg << ") in spatial database '" << queryctx->db->label() << "'.";
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
    } // if
    if (vp <= 0) {
        std::ostringstream msg;
        msg << "Found negative dilatational wave speed (" << vp << ") at location (";
        for (int i=0; i < dim; ++i)
            msg << "  " << xDim[i];
        msg << ") in spatial database '" << queryctx->db->label() << "'.";
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
    } // if

    const PylithReal bulkModulus = density*(vp*vp - 4.0/3.0*vs*vs); assert(bulkModulus > 0.0);
    assert(queryctx->valueScale > 0);
    values[0] = bulkModulus / queryctx->valueScale;

    PYLITH_METHOD_RETURN(0);
} // dbQueryBulkModulus


// ----------------------------------------------------------------------
// Query for Vs, density, and viscosity to determine Maxwell time.
PetscErrorCode
pylith::materials::Query::dbQueryMaxwellTime(PylithInt dim,
					     PylithReal t,
					     const PylithReal x[],
					     PylithInt nvalues,
					     PylithScalar* values,
					     void* context)
{ // dbQueryMaxwellTime
    PYLITH_METHOD_BEGIN;

    const int _nvalues = 1;

    assert(x);
    assert(values);
    assert(context);
    assert(_nvalues == nvalues);
    assert(2 == dim || 3 == dim);

    const pylith::topology::FieldQuery::DBQueryContext* queryctx = (pylith::topology::FieldQuery::DBQueryContext*)context; assert(queryctx);

    // Tell database which values we want.
    const int numDBValues = 3;
    PylithReal dbValues[numDBValues];
    const int i_density = 0;
    const int i_vs = 1;
    const int i_viscosity = 2;
    const char* dbValueNames[numDBValues] = {"density", "vs", "viscosity"};
    try {
        queryctx->db->queryVals(dbValueNames, numDBValues);
    } catch (const std::exception& err) {
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, err.what());
    } // try/catch


    // Dimensionalize query location coordinates.
    assert(queryctx->lengthScale > 0);
    double xDim[3];
    for (int i=0; i < dim; ++i) {
        xDim[i] = x[i] * queryctx->lengthScale;
    } // for

    assert(queryctx->cs);
    const int err = queryctx->db->query(dbValues, numDBValues, xDim, dim, queryctx->cs);
    if (err) {
        std::ostringstream msg;
        msg << "Could not find density, Vs, and viscosity at (";
        for (int i=0; i < dim; ++i)
            msg << "  " << xDim[i];
        msg << ") using spatial database '" << queryctx->db->label() << "'.";
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
    } // if

    const PylithReal density = dbValues[i_density];
    const PylithReal vs = dbValues[i_vs];
    const PylithReal viscosity = dbValues[i_viscosity];

    if (density <= 0) {
        std::ostringstream msg;
        msg << "Found negative density (" << density << ") at location (";
        for (int i=0; i < dim; ++i)
            msg << "  " << xDim[i];
        msg << ") in spatial database '" << queryctx->db->label() << "'.";
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
    } // if
    if (vs <= 0) {
        std::ostringstream msg;
        msg << "Found negative shear wave speed (" << vs << ") at location (";
        for (int i=0; i < dim; ++i)
            msg << "  " << xDim[i];
        msg << ") in spatial database '" << queryctx->db->label() << "'.";
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
    } // if

    const PylithReal shearModulus = density * vs * vs; assert(shearModulus > 0);
    const PylithReal maxwellTime = viscosity/shearModulus; assert(maxwellTime > 0);
    assert(queryctx->valueScale > 0);
    values[0] = maxwellTime / queryctx->valueScale;

    PYLITH_METHOD_RETURN(0);
} // dbQueryMaxwellTime


// ----------------------------------------------------------------------
// Query for components of gravity field.
PetscErrorCode
pylith::materials::Query::dbQueryGravityField(PylithInt dim,
                                              PylithReal t,
                                              const PylithReal x[],
                                              PylithInt nvalues,
                                              PylithScalar* values,
                                              void* context)
{ // dbQueryGravityField
    PYLITH_METHOD_BEGIN;

    const int _nvalues = dim;

    assert(x);
    assert(values);
    assert(context);
    assert(_nvalues == nvalues);
    assert(2 == dim || 3 == dim);

    const pylith::topology::FieldQuery::DBQueryContext* queryctx = (pylith::topology::FieldQuery::DBQueryContext*)context; assert(queryctx);

    // Tell database which values we want.
    const int numDBValues = dim;
    PylithReal dbValues[3];
    const int i_x = 0;
    const int i_y = 1;
    const int i_z = 2;
    const char* dbValueNames[3] = {"x", "y", "z"};
    try {
        queryctx->db->queryVals(dbValueNames, numDBValues);
    } catch (const std::exception& err) {
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, err.what());
    } // try/catch


    // Dimensionalize query location coordinates.
    assert(queryctx->lengthScale > 0);
    double xDim[3];
    for (int i=0; i < dim; ++i) {
        xDim[i] = x[i] * queryctx->lengthScale;
    } // for

    assert(queryctx->cs);
    const int err = queryctx->db->query(dbValues, numDBValues, xDim, dim, queryctx->cs);
    if (err) {
        std::ostringstream msg;
        msg << "Could not find gravity field at (";
        for (int i=0; i < dim; ++i)
            msg << "  " << xDim[i];
        msg << ") using spatial database '" << queryctx->db->label() << "'.";
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
    } // if

    PylithReal mag = 0.0;
    for (int i=0; i < dim; ++i) {
        mag += values[i]*values[i];
    } // for
    const PylithReal tolerance = 1.0e-6;
    if (mag < tolerance) {
        std::ostringstream msg;
        msg << "Found near zero magnitude (" << mag << ") for gravity field vector at location (";
        for (int i=0; i < dim; ++i)
            msg << "  " << xDim[i];
        msg << ") in spatial database '" << queryctx->db->label() << "'.";
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
    } // if

    assert(queryctx->valueScale > 0);
    if (2 == dim) {
        values[0] = dbValues[i_x] / queryctx->valueScale;
        values[1] = dbValues[i_y] / queryctx->valueScale;
    } else {
        values[0] = dbValues[i_x] / queryctx->valueScale;
        values[1] = dbValues[i_y] / queryctx->valueScale;
        values[2] = dbValues[i_z] / queryctx->valueScale;
    } // if/else

    PYLITH_METHOD_RETURN(0);
} // dbQueryGravityField


// End of file
