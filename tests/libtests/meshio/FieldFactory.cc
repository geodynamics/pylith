// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FieldFactory.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include <cppunit/extensions/HelperMacros.h>

// ----------------------------------------------------------------------
// Default constructor.
pylith::meshio::FieldFactory::FieldFactory(pylith::topology::Field& field) :
    _field(field) {}


// ----------------------------------------------------------------------
// Destructor.
pylith::meshio::FieldFactory::~FieldFactory(void) {}


// ----------------------------------------------------------------------
// Add scalar field.
void
pylith::meshio::FieldFactory::addScalar(const pylith::topology::FieldBase::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;

    const char* fieldName = "scalar";

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = "scalar";
    description.scale = 1.0;
    description.validator = NULL;

    _field.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addScalar


// ----------------------------------------------------------------------
// Add vector field.
void
pylith::meshio::FieldFactory::addVector(const pylith::topology::FieldBase::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;

    const char* fieldName = "vector";
    const char* components[3] = { "vector_x", "vector_y", "vector_z" };

    const int spaceDim = _field.getSpaceDim();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = spaceDim;
    description.componentNames.resize(spaceDim);
    for (int i = 0; i < spaceDim; ++i) {
        description.componentNames[i] = components[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addVector


// ----------------------------------------------------------------------
// Add tensor field.
void
pylith::meshio::FieldFactory::addTensor(const pylith::topology::FieldBase::Discretization& discretization) {
    const char* fieldName = "tensor";

    const int spaceDim = _field.getSpaceDim();
    CPPUNIT_ASSERT(2 == spaceDim || 3 == spaceDim);

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::TENSOR;
    if (2 == spaceDim) {
        const int tensorSize = 3;
        const char* componentNames[tensorSize] = { "tensor_xx", "tensor_yy", "tensor_xy" };
        description.numComponents = tensorSize;
        description.componentNames.resize(tensorSize);
        for (int i = 0; i < tensorSize; ++i) {
            description.componentNames[i] = componentNames[i];
        } // for
    } else if (3 == spaceDim) {
        const int tensorSize = 6;
        const char* componentNames[6] = { "tensor_xx", "tensor_yy", "tensor_zz", "tensor_xy", "tensor_yz", "tensor_xz" };
        description.numComponents = tensorSize;
        description.componentNames.resize(tensorSize);
        for (int i = 0; i < tensorSize; ++i) {
            description.componentNames[i] = componentNames[i];
        } // for
    } else {
        throw std::logic_error("Unknown spatial dimension.");
    } // if/else
    description.scale = 1.0;
    description.validator = NULL;

    _field.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addTensor


// ----------------------------------------------------------------------
// Add other field.
void
pylith::meshio::FieldFactory::addOther(const pylith::topology::FieldBase::Discretization& discretization) {
    const char* fieldName = "other";
    const int otherSize = 2;
    const char* componentNames[otherSize] = { "other_1", "other_2" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = otherSize;
    description.componentNames.resize(otherSize);
    for (int i = 0; i < otherSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addOther


// ----------------------------------------------------------------------
void
pylith::meshio::FieldFactory::setValues(const PylithScalar* values,
                                        const PylithInt numPoints,
                                        const PylithInt numDOF) {
    PYLITH_METHOD_BEGIN;

    pylith::topology::VecVisitorMesh fieldVisitor(_field);
    PylithScalar* fieldArray = fieldVisitor.localArray();CPPUNIT_ASSERT(fieldArray);
    const PylithInt fieldSize = numPoints * numDOF;
    CPPUNIT_ASSERT_EQUAL(fieldSize, _field.getStorageSize());
    for (PylithInt i = 0; i < fieldSize; ++i) {
        fieldArray[i] = values[i];
    } // for

    PYLITH_METHOD_END;
} // setField


// End of file
