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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestOutputSolnPoints.hh" // Implementation of class methods

#include "pylith/meshio/OutputSolnPoints.hh"

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "data/OutputSolnPointsDataTri3.hh"
#include "data/OutputSolnPointsDataQuad4.hh"
#include "data/OutputSolnPointsDataTet4.hh"
#include "data/OutputSolnPointsDataHex8.hh"

#include <string.h> // USES strcmp()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestOutputSolnPoints );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestOutputSolnPoints::testConstructor(void)
{ // testConstructor
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
// Test setupInterpolator for tri3 mesh.
void
pylith::meshio::TestOutputSolnPoints::testSetupInterpolatorTri3(void)
{ // testSetupInterpolatorTri3
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataTri3 data;

    _testSetupInterpolator(data);

    PYLITH_METHOD_END;
} // testSetupInterpolatorTri3


// ----------------------------------------------------------------------
// Test interpolation for tri3 mesh.
void
pylith::meshio::TestOutputSolnPoints::testInterpolateTri3(void)
{ // testInterpolateTri3
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataTri3 data;

    _testInterpolate(data);

    PYLITH_METHOD_END;
} // testInterpolateTri3


// ----------------------------------------------------------------------
// Test setupInterpolator for quad4 mesh.
void
pylith::meshio::TestOutputSolnPoints::testSetupInterpolatorQuad4(void)
{ // testSetupInterpolatorQuad4
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataQuad4 data;

    _testSetupInterpolator(data);

    PYLITH_METHOD_END;
} // testSetupInterpolatorQuad4


// ----------------------------------------------------------------------
// Test interpolation for quad4 mesh.
void
pylith::meshio::TestOutputSolnPoints::testInterpolateQuad4(void)
{ // testInterpolateQuad4
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataQuad4 data;

    _testInterpolate(data);

    PYLITH_METHOD_END;
} // testInterpolateQuad4


// ----------------------------------------------------------------------
// Test setupInterpolator for tet4 mesh.
void
pylith::meshio::TestOutputSolnPoints::testSetupInterpolatorTet4(void)
{ // testSetupInterpolatorTet4
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataTet4 data;

    _testSetupInterpolator(data);

    PYLITH_METHOD_END;
} // testSetupInterpolatorTet4


// ----------------------------------------------------------------------
// Test interpolation for tet4 mesh.
void
pylith::meshio::TestOutputSolnPoints::testInterpolateTet4(void)
{ // testInterpolateTet4
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataTet4 data;

    _testInterpolate(data);

    PYLITH_METHOD_END;
} // testInterpolateTet4


// ----------------------------------------------------------------------
// Test setupInterpolator for hex8 mesh.
void
pylith::meshio::TestOutputSolnPoints::testSetupInterpolatorHex8(void)
{ // testSetupInterpolatorHex8
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataHex8 data;

    _testSetupInterpolator(data);

    PYLITH_METHOD_END;
} // testSetupInterpolatorHex8


// ----------------------------------------------------------------------
// Test interpolation for hex8 mesh.
void
pylith::meshio::TestOutputSolnPoints::testInterpolateHex8(void)
{ // testInterpolateHex8
    PYLITH_METHOD_BEGIN;

    OutputSolnPoints output;
    OutputSolnPointsDataHex8 data;

    _testInterpolate(data);

    PYLITH_METHOD_END;
} // testInterpolateHex8


// ----------------------------------------------------------------------
// Test setupInterpolator().
void
pylith::meshio::TestOutputSolnPoints::_testSetupInterpolator(const OutputSolnPointsData& data)
{ // _testSetupInterpolator
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
    cs.initialize();
    mesh.coordsys(&cs);
    MeshIOAscii iohandler;
    iohandler.filename(data.meshFilename);
    iohandler.read(&mesh);

    OutputSolnPoints output;
    CPPUNIT_ASSERT(data.points);
    output.setupInterpolator(&mesh, data.points, numPoints, spaceDim, data.names, numPoints, normalizer);

    PetscDM dmMesh = output.pointsMesh().dmMesh(); CPPUNIT_ASSERT(dmMesh);

    // Check vertices
    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();
    CPPUNIT_ASSERT_EQUAL(numVerticesE, verticesStratum.size());
    for (PetscInt v=vStart, index = 0; v < vEnd; ++v, ++index) {
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

        err = DMPlexGetConeSize(dmMesh, c, &coneSize); PYLITH_CHECK_ERROR(err);
        err = DMPlexGetCone(dmMesh, c, &cone); PYLITH_CHECK_ERROR(err);

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
pylith::meshio::TestOutputSolnPoints::_testInterpolate(const OutputSolnPointsData& data)
{ // _testInterpolate
    PYLITH_METHOD_BEGIN;
#if 0

    const int numPoints = data.numPoints;
    const int spaceDim = data.spaceDim;

    const int numVerticesE = numPoints;
    const int numCellsE = numPoints;
    const int numCornersE = 1;

    topology::Mesh mesh;
    spatialdata::geocoords::CSCart cs;
    spatialdata::units::Nondimensional normalizer;

    cs.setSpaceDim(spaceDim);
    cs.initialize();
    mesh.coordsys(&cs);
    MeshIOAscii iohandler;
    iohandler.filename(data.meshFilename);
    iohandler.read(&mesh);

    OutputSolnPoints output;
    CPPUNIT_ASSERT(data.points);
    output.setupInterpolator(&mesh, data.points, numPoints, spaceDim, data.names, numPoints, normalizer);

    // Create field with data.
    const char* fieldName = "data_field";

    // Create field for interpolated data.

    PetscDM pointsMeshDM = output.pointsMesh().dmMesh(); CPPUNIT_ASSERT(pointsMeshDM);
    pylith::topology::Field fieldInterp(pointsMeshDM);
    fieldInterp.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
    fieldInterp.allocate();
    fieldInterp.label(field.label());
    fieldInterp.vectorFieldType(field.vectorFieldType());
    fieldInterp.scale(field.scale());
    fieldInterp.zeroAll();

    const char* context = field.c_str();
    fieldInterp.createScatter(*_pointsMesh, context);

    PetscVec fieldInterpVec = fieldInterp.vector(context); assert(fieldInterpVec);
    err = DMInterpolationSetDof(output._interpolator, fiberDim); PYLITH_CHECK_ERROR(err);
    err = DMInterpolationEvaluate(output._interpolator, pointsMeshDM, field.localVector(), fieldInterpVec); PYLITH_CHECK_ERROR(err);

    // Check interpolated field
    topology::Stratum verticesStratum(pointsMeshDM, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();
    CPPUNIT_ASSERT_EQUAL(numVerticesE, verticesStratum.size());
    for (PetscInt v=vStart, index = 0; v < vEnd; ++v, ++index) {
        const int vertexE = numCellsE + index;
	// :TODO: STUFF GOES HERE
    } // for

#endif
    PYLITH_METHOD_END;
} // _testInterpolate


// End of file
