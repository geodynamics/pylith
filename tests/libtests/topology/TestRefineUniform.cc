// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestRefineUniform.hh" // Implementation of class methods

#include "pylith/topology/RefineUniform.hh" // USES RefineUniform
#include "tests/src/FaultCohesiveStub.hh" // USES FaultCohesiveStub

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "pylith/utils/array.hh" // USES int_array

#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::logic_error

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
// Setup testing data.
pylith::topology::TestRefineUniform::TestRefineUniform(TestRefineUniform_Data* data) :
    _data(data) {
    REQUIRE(_data);
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::topology::TestRefineUniform::~TestRefineUniform(void) {
    delete _data;_data = NULL;
} // destructor


#include <iostream>
// ------------------------------------------------------------------------------------------------
// Test refine().
void
pylith::topology::TestRefineUniform::testRefine(void) {
    PYLITH_METHOD_BEGIN;
    REQUIRE(_data);

    Mesh mesh(_data->cellDim);
    _initializeMesh(&mesh);

    RefineUniform refiner;
    Mesh newMesh(_data->cellDim);
    refiner.refine(&newMesh, mesh, _data->refineLevel);

    // Check mesh dimension
    REQUIRE(_data->cellDim == newMesh.getDimension());

    const PetscDM dmNewMesh = newMesh.getDM();REQUIRE(dmNewMesh);

    // Check vertices
    pylith::topology::Stratum verticesStratum(dmNewMesh, topology::Stratum::DEPTH, 0);
    REQUIRE(_data->numVertices == verticesStratum.size());
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();

    pylith::topology::CoordsVisitor coordsVisitor(dmNewMesh);
    const int spaceDim = _data->spaceDim;
    for (PetscInt v = vStart; v < vEnd; ++v) {
        CHECK(spaceDim == coordsVisitor.sectionDof(v));
    } // for

    // Check cells
    pylith::topology::Stratum cellsStratum(dmNewMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();
    const PetscInt numCells = cellsStratum.size();

    REQUIRE(_data->numCells+_data->numCellsCohesive == numCells);
    // Normal cells
    for (PetscInt c = cStart; c < _data->numCells; ++c) {
        DMPolytopeType ct;
        PetscInt *closure = NULL;
        PetscInt closureSize, numCorners = 0;

        PylithCallPetscRequire(DMPlexGetCellType(dmNewMesh, c, &ct));
        PylithCallPetscRequire(DMPlexGetTransitiveClosure(dmNewMesh, c, PETSC_TRUE, &closureSize, &closure));
        for (PetscInt p = 0; p < closureSize*2; p += 2) {
            const PetscInt point = closure[p];
            if ((point >= vStart) && (point < vEnd)) {
                closure[numCorners++] = point;
            } // if
        } // for
        PylithCallPetscRequire(DMPlexInvertCell(ct, closure));
        CHECK(_data->numCorners == numCorners);
        PylithCallPetscRequire(DMPlexRestoreTransitiveClosure(dmNewMesh, c, PETSC_TRUE, &closureSize, &closure));
    } // for

    // Cohesive cells
    for (PetscInt c = _data->numCells; c < cEnd; ++c) {
        DMPolytopeType ct;
        PetscInt *closure = NULL;
        PetscInt closureSize, numCorners = 0;

        PylithCallPetscRequire(DMPlexGetCellType(dmNewMesh, c, &ct));
        PylithCallPetscRequire(DMPlexGetTransitiveClosure(dmNewMesh, c, PETSC_TRUE, &closureSize, &closure));
        for (PetscInt p = 0; p < closureSize*2; p += 2) {
            const PetscInt point = closure[p];
            if ((point >= vStart) && (point < vEnd)) {
                closure[numCorners++] = point;
            } // if
        } // for
        PylithCallPetscRequire(DMPlexInvertCell(ct, closure));
        CHECK(_data->numCornersCohesive == numCorners);
        PylithCallPetscRequire(DMPlexRestoreTransitiveClosure(dmNewMesh, c, PETSC_TRUE, &closureSize, &closure));
    } // for

    // check materials
    PetscInt matId = 0;
    PetscInt matIdSum = 0; // Use sum of material ids as simple checksum.
    for (PetscInt c = cStart; c < cEnd; ++c) {
        PylithCallPetscRequire(DMGetLabelValue(dmNewMesh, pylith::topology::Mesh::cells_label_name, c, &matId));
        matIdSum += matId;
    } // for
    CHECK(_data->matIdSum == matIdSum);

    // Check vertex groups
    pylith::string_vector vertexGroupNames;
    pylith::meshio::MeshBuilder::getVertexGroupNames(&vertexGroupNames, newMesh);
    for (size_t iGroup = 0; iGroup < _data->numVertexGroups; ++iGroup) {
        INFO("Checking vertex group '"<<_data->vertexGroupNames[iGroup]<<"'.");
        int_array points;
        pylith::meshio::MeshBuilder::getVertexGroup(&points, newMesh, _data->vertexGroupNames[iGroup]);
        REQUIRE(_data->vertexGroupSizes[iGroup] == points.size());
        for (size_t iPoint = 0; iPoint < points.size(); ++iPoint) {
            CHECK(points[iPoint] >= 0);
            CHECK(points[iPoint] < _data->numVertices);
        } // for
    } // for

    // Check face groups
    pylith::topology::Stratum facesStratum(dmNewMesh, topology::Stratum::HEIGHT, 1);
    const PetscInt fStart = facesStratum.begin();
    const PetscInt fEnd = facesStratum.end();
    pylith::string_vector faceGroupNames;
    pylith::meshio::MeshBuilder::getFaceGroupNames(&faceGroupNames, newMesh);
    REQUIRE(_data->numFaceGroups == faceGroupNames.size());
    for (size_t iGroup = 0; iGroup < _data->numFaceGroups; ++iGroup) {
        INFO("Checking face group '"<<_data->faceGroupNames[iGroup]<<"'.");
        PetscInt numFaces = 0;
        const PetscInt labelValue = 1;
        PylithCallPetscRequire(DMGetStratumSize(dmNewMesh, _data->faceGroupNames[iGroup], labelValue, &numFaces));
        REQUIRE(_data->faceGroupSizes[iGroup] == size_t(numFaces));
        PetscIS facesIS = NULL;
        const PetscInt *faces = NULL;
        PylithCallPetscRequire(DMGetStratumIS(dmNewMesh, _data->faceGroupNames[iGroup], labelValue, &facesIS));
        PylithCallPetscRequire(ISGetIndices(facesIS, &faces));
        for (PetscInt iFace = 0; iFace < numFaces; ++iFace) {
            CHECK((faces[iFace] >= fStart && faces[iFace] < fEnd));
        } // for
        PylithCallPetscRequire(ISRestoreIndices(facesIS, &faces));
        PylithCallPetscRequire(ISDestroy(&facesIS));

    } // for

    PYLITH_METHOD_END;
} // testRefine


// ------------------------------------------------------------------------------------------------
void
pylith::topology::TestRefineUniform::_initializeMesh(Mesh* const mesh) {
    PYLITH_METHOD_BEGIN;
    REQUIRE(_data);
    REQUIRE(mesh);

    pylith::meshio::MeshIOAscii iohandler;
    iohandler.setFilename(_data->filename);
    iohandler.read(mesh);

    // Adjust topology if necessary.
    if (_data->faultA) {
        faults::FaultCohesiveStub faultA;
        faultA.setCohesiveLabelValue(100);
        faultA.setSurfaceLabelName(_data->faultA);
        faultA.adjustTopology(mesh);
    } // if

    if (_data->faultB) {
        faults::FaultCohesiveStub faultB;
        faultB.setCohesiveLabelValue(101);
        faultB.setSurfaceLabelName(_data->faultB);
        faultB.adjustTopology(mesh);
    } // if

    PYLITH_METHOD_END;
} // _initializeMesh


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::topology::TestRefineUniform_Data::TestRefineUniform_Data(void) :
    filename(NULL),
    refineLevel(0),
    faultA(NULL),
    faultB(NULL),
    isSimplexMesh(true),
    numVertices(0),
    spaceDim(0),
    cellDim(0),
    numCells(0),
    numCorners(0),
    numCellsCohesive(0),
    numCornersCohesive(0),
    matIdSum(0),
    vertexGroupSizes(NULL),
    vertexGroupNames(NULL),
    numVertexGroups(0),
    faceGroupSizes(NULL),
    faceGroupNames(NULL),
    numFaceGroups(0) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::topology::TestRefineUniform_Data::~TestRefineUniform_Data(void) {}


// End of file
