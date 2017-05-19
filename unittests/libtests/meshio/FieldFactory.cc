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

#include "FieldFactory.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include <cppunit/extensions/HelperMacros.h>

// ----------------------------------------------------------------------
// Default constructor.
pylith::meshio::FieldFactory::FieldFactory(pylith::topology::Fields& fields) :
    _fields(fields)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::meshio::FieldFactory::~FieldFactory(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Add scalar field.
void
pylith::meshio::FieldFactory::scalar(const pylith::topology::FieldBase::Discretization& discretization,
                                     const PylithScalar* values,
                                     const PylithInt numPoints,
                                     const PylithInt numComponents)
{ // scalar
    PYLITH_METHOD_BEGIN;

    const char* fieldName = "scalar";

    _fields.add(fieldName, fieldName);
    pylith::topology::Field& field = _fields.get(fieldName);

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = "scalar";
    description.scale = 1.0;
    description.validator = NULL;

    field.subfieldAdd(description, discretization);
    field.subfieldsSetup();
    field.allocate();

    this->setField(&field, values, numPoints, numComponents);

    PYLITH_METHOD_END;
} // scalar

// ----------------------------------------------------------------------
// Add vector field.
void
pylith::meshio::FieldFactory::vector(const pylith::topology::FieldBase::Discretization& discretization,
                                     const PylithScalar* values,
                                     const PylithInt numPoints,
                                     const PylithInt numComponents)
{ // vector
    PYLITH_METHOD_BEGIN;

    const char* fieldName = "vector";
    const char* components[3] = { "vector_x", "vector_y", "vector_z" };

    _fields.add(fieldName, fieldName);
    pylith::topology::Field& field = _fields.get(fieldName);
    const int spaceDim = field.spaceDim();

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

    field.subfieldAdd(description, discretization);
    field.subfieldsSetup();
    field.allocate();

    this->setField(&field, values, numPoints, numComponents);

    PYLITH_METHOD_END;
} // vector

// ----------------------------------------------------------------------
// Add tensor field.
void
pylith::meshio::FieldFactory::tensor(const pylith::topology::FieldBase::Discretization& discretization,
                                     const PylithScalar* values,
                                     const PylithInt numPoints,
                                     const PylithInt numComponents)
{ // tensor

    const char* fieldName = "tensor";

    _fields.add(fieldName, fieldName);
    pylith::topology::Field& field = _fields.get(fieldName);
    const int spaceDim = field.spaceDim();
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

    field.subfieldAdd(description, discretization);
    field.subfieldsSetup();
    field.allocate();

    this->setField(&field, values, numPoints, numComponents);

    PYLITH_METHOD_END;
} // tensor

// ----------------------------------------------------------------------
// Add other field.
void
pylith::meshio::FieldFactory::other(const pylith::topology::FieldBase::Discretization& discretization,
                                    const PylithScalar* values,
                                    const PylithInt numPoints,
                                    const PylithInt numComponents)
{ // other
    const char* fieldName = "other";
    const int otherSize = 2;
    const char* componentNames[otherSize] = { "other_1", "other_2" };

    _fields.add(fieldName, fieldName);
    pylith::topology::Field& field = _fields.get(fieldName);

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

    field.subfieldAdd(description, discretization);
    field.subfieldsSetup();
    field.allocate();

    this->setField(&field, values, numPoints, numComponents);

    PYLITH_METHOD_END;
} // other


// ----------------------------------------------------------------------
void
pylith::meshio::FieldFactory::setField(pylith::topology::Field* field,
                                       const PylithScalar* values,
                                       const PylithInt numPoints,
                                       const PylithInt numComponents)
{ // setField
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(field);

    pylith::topology::VecVisitorMesh fieldVisitor(*field);
    PylithScalar* fieldArray = fieldVisitor.localArray();CPPUNIT_ASSERT(fieldArray);
    const PylithInt fieldSize = numPoints * numComponents;
    CPPUNIT_ASSERT_EQUAL(fieldSize, field->sectionSize());
    for (PylithInt i = 0; i < fieldSize; ++i) {
        fieldArray[i] = values[i];
    } // for
    field->view("values set");

    PYLITH_METHOD_END;
} // setField


// End of file
