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

#include "TestOutputSolnPoints.hh" // Implementation of class methods

#include "pylith/meshio/OutputSolnPoints.hh"

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/meshio/MeshIOCubit.hh" // USES MeshIOCubit

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "data/OutputSolnPointsDataTri3.hh"
#include "data/OutputSolnPointsDataQuad4.hh"
#include "data/OutputSolnPointsDataTet4.hh"
#include "data/OutputSolnPointsDataHex8.hh"

#include <string.h> // USES strcmp()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::meshio::TestOutputSolnPoints);

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestOutputSolnPoints::testConstructor(void) { // testConstructor
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
// Test setupInterpolator for tri3 mesh.
void
pylith::meshio::TestOutputSolnPoints::testSetupInterpolatorTri3(void) { // testSetupInterpolatorTri3
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataTri3 data;

    _testSetupInterpolator(data);

    PYLITH_METHOD_END;
} // testSetupInterpolatorTri3


// ----------------------------------------------------------------------
// Test interpolation for tri3 mesh.
void
pylith::meshio::TestOutputSolnPoints::testInterpolateTri3(void) { // testInterpolateTri3
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataTri3 data;

    _testInterpolate(data);

    PYLITH_METHOD_END;
} // testInterpolateTri3


// ----------------------------------------------------------------------
// Test setupInterpolator for quad4 mesh.
void
pylith::meshio::TestOutputSolnPoints::testSetupInterpolatorQuad4(void) { // testSetupInterpolatorQuad4
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataQuad4 data;

    _testSetupInterpolator(data);

    PYLITH_METHOD_END;
} // testSetupInterpolatorQuad4


// ----------------------------------------------------------------------
// Test interpolation for quad4 mesh.
void
pylith::meshio::TestOutputSolnPoints::testInterpolateQuad4(void) { // testInterpolateQuad4
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataQuad4 data;

    _testInterpolate(data);

    PYLITH_METHOD_END;
} // testInterpolateQuad4


// ----------------------------------------------------------------------
// Test setupInterpolator for tet4 mesh.
void
pylith::meshio::TestOutputSolnPoints::testSetupInterpolatorTet4(void) { // testSetupInterpolatorTet4
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataTet4 data;

    _testSetupInterpolator(data);

    PYLITH_METHOD_END;
} // testSetupInterpolatorTet4


// ----------------------------------------------------------------------
// Test interpolation for tet4 mesh.
void
pylith::meshio::TestOutputSolnPoints::testInterpolateTet4(void) { // testInterpolateTet4
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataTet4 data;

    _testInterpolate(data);

    PYLITH_METHOD_END;
} // testInterpolateTet4


// ----------------------------------------------------------------------
// Test setupInterpolator for hex8 mesh.
void
pylith::meshio::TestOutputSolnPoints::testSetupInterpolatorHex8(void) { // testSetupInterpolatorHex8
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataHex8 data;

    _testSetupInterpolator(data);

    PYLITH_METHOD_END;
} // testSetupInterpolatorHex8


// ----------------------------------------------------------------------
// Test interpolation for hex8 mesh.
void
pylith::meshio::TestOutputSolnPoints::testInterpolateHex8(void) { // testInterpolateHex8
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataHex8 data;

    _testInterpolate(data);

    PYLITH_METHOD_END;
} // testInterpolateHex8


// ----------------------------------------------------------------------
// Test setupInterpolator().
void
pylith::meshio::TestOutputSolnPoints::_testSetupInterpolator(const OutputSolnPointsData& data) { // _testSetupInterpolator
    PYLITH_METHOD_BEGIN;

    const int numPoints = data.numPoints;
    const int spaceDim = data.spaceDim;

    const int numVerticesE = numPoints;
    const int numCellsE = numPoints;
    const int numCornersE = 1;

    topology::Mesh mesh;
    spatialdata::geocoords::CSCart cs;
    spatialdata::units::Nondimensional normalizer;

    cs.setSpaceDim(spaceDim);
    mesh.setCoordSys(&cs);
    MeshIOCubit iohandler;
    iohandler.filename(data.meshFilename);
    iohandler.read(&mesh);

    OutputSolnPoints output;
    CPPUNIT_ASSERT(data.points);
    output.setupInterpolator(&mesh, data.points, numPoints, spaceDim, data.names, numPoints, normalizer);

    PetscDM dmMesh = output.pointsMesh().dmMesh();CPPUNIT_ASSERT(dmMesh);

    // Check vertices
    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();
    CPPUNIT_ASSERT_EQUAL(numVerticesE, verticesStratum.size());
    for (PetscInt v = vStart, index = 0; v < vEnd; ++v, ++index) {
        const int vertexE = numCellsE + index;
        CPPUNIT_ASSERT_EQUAL(vertexE, v);
    } // for

    // Check cells
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();
    PetscErrorCode err = 0;
    CPPUNIT_ASSERT_EQUAL(numCellsE, cellsStratum.size());
    for (PetscInt c = cStart, index = 0; c < cEnd; ++c) {
        const PetscInt *cone = NULL;
        PetscInt coneSize = 0;

        err = DMPlexGetConeSize(dmMesh, c, &coneSize);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetCone(dmMesh, c, &cone);PYLITH_CHECK_ERROR(err);

        CPPUNIT_ASSERT_EQUAL(numCornersE, coneSize);

        for (PetscInt p = 0; p < coneSize; ++p, ++index) {
            const int coneE = numCellsE+index;
            CPPUNIT_ASSERT_EQUAL(coneE, cone[p]);
        } // for
    } // for

    PYLITH_METHOD_END;
} // _testSetupInterpolator


// ----------------------------------------------------------------------
// Test interpolation.
void
pylith::meshio::TestOutputSolnPoints::_testInterpolate(const OutputSolnPointsData& data) { // _testInterpolate
    PYLITH_METHOD_BEGIN;

    const int numPoints = data.numPoints;
    const int spaceDim = data.spaceDim;

    topology::Mesh mesh;
    spatialdata::geocoords::CSCart cs;
    spatialdata::units::Nondimensional normalizer;

    cs.setSpaceDim(spaceDim);
    mesh.setCoordSys(&cs);
    MeshIOCubit iohandler;
    iohandler.filename(data.meshFilename);
    iohandler.read(&mesh);

    OutputSolnPoints output;
    CPPUNIT_ASSERT(data.points);
    output.setupInterpolator(&mesh, data.points, numPoints, spaceDim, data.names, numPoints, normalizer);

    // Create field with data.
    const char* fieldName = "data_field";
    const int fiberDim = data.fiberDim;
    pylith::topology::Field field(mesh);
    field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
    field.allocate();
    field.setLabel(fieldName);
    this->_calcField(&field, data);

    // Create field for interpolated data.
    pylith::topology::Field fieldInterp(output.pointsMesh());
    fieldInterp.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
    fieldInterp.allocate();
    fieldInterp.setLabel(field.getLabel());
    fieldInterp.vectorFieldType(field.vectorFieldType());
    fieldInterp.scale(field.scale());

    // Create field to populate with expected data.
    pylith::topology::Field fieldInterpE(fieldInterp.mesh());
    fieldInterpE.cloneSection(fieldInterp);
    fieldInterpE.allocate();
    this->_calcField(&fieldInterpE, data);

    PetscDM pointsMeshDM = output.pointsMesh().dmMesh();CPPUNIT_ASSERT(pointsMeshDM);
    PetscErrorCode err;
    err = DMInterpolationSetDof(output._interpolator, fiberDim);PYLITH_CHECK_ERROR(err);
    err = DMInterpolationEvaluate(output._interpolator, field.dmMesh(), field.localVector(), fieldInterp.localVector());PYLITH_CHECK_ERROR(err);

    // Check interpolated field
    topology::Stratum verticesStratum(pointsMeshDM, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();

    topology::VecVisitorMesh fieldInterpVisitor(fieldInterp);
    PetscScalar* fieldInterpArray = fieldInterpVisitor.localArray();CPPUNIT_ASSERT(fieldInterpArray);

    topology::VecVisitorMesh fieldInterpEVisitor(fieldInterpE);
    PetscScalar* fieldInterpEArray = fieldInterpEVisitor.localArray();CPPUNIT_ASSERT(fieldInterpEArray);

    const double tolerance = 1.0e-6;
    for (PetscInt v = vStart; v < vEnd; ++v) {
        const PetscInt off = fieldInterpVisitor.sectionOffset(v);
        CPPUNIT_ASSERT_EQUAL(fiberDim, fieldInterpVisitor.sectionDof(v));

        const PetscInt offE = fieldInterpEVisitor.sectionOffset(v);
        CPPUNIT_ASSERT_EQUAL(fiberDim, fieldInterpEVisitor.sectionDof(v));

        for (PetscInt d = 0; d < fiberDim; ++d) {
            const PylithScalar valueE = fieldInterpEArray[offE+d];
            const PylithScalar value = fieldInterpArray[off+d];
            if (fabs(valueE) > 1.0) {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0,  value / valueE, tolerance);
            } else {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(valueE, value, tolerance);
            } // else
        } // for
    } // for

    PYLITH_METHOD_END;
} // _testInterpolate


// ----------------------------------------------------------------------
void
pylith::meshio::TestOutputSolnPoints::_calcField(pylith::topology::Field* field,
                                                 const OutputSolnPointsData& data) { // _calcField
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(field);

    const pylith::topology::Mesh& mesh = field->mesh();
    pylith::topology::Stratum verticesStratum(mesh.dmMesh(), topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();

    pylith::topology::VecVisitorMesh fieldVisitor(*field);
    PylithScalar* fieldArray = fieldVisitor.localArray();CPPUNIT_ASSERT(fieldArray);

    pylith::topology::CoordsVisitor coordsVisitor(mesh.dmMesh());
    PylithScalar* coordsArray = coordsVisitor.localArray();CPPUNIT_ASSERT(coordsArray);

    const PetscInt fiberDim = data.fiberDim;
    const PetscInt spaceDim = data.spaceDim;

    for (PetscInt v = vStart, index = 0; v < vEnd; ++v) {
        const PetscInt foff = fieldVisitor.sectionOffset(v);
        CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));

        const PetscInt coff = coordsVisitor.sectionOffset(v);
        CPPUNIT_ASSERT_EQUAL(spaceDim, coordsVisitor.sectionDof(v));

        for (PetscInt d = 0; d < fiberDim; ++d, ++index) {
            PylithScalar value = 0.0;
            for (PetscInt iv = 0; iv < spaceDim; ++iv) {
                const PetscInt ic = d*spaceDim + iv;
                value += data.coefs[ic]*coordsArray[coff+iv];
            } // for
            fieldArray[foff+d] = value;
        } // for
    } // for

    PYLITH_METHOD_END;
} // _calcField


// End of file
