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

#include "TestFaultKin.hh" // Implementation of class methods

#include "pylith/bc/BoundaryCondition.cc" // HASA BoundaryCondition
#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/KinSrc.hh" // USES KinSrc
#include "pylith/materials/Elasticity.hh" // USES Elasticity
#include "pylith/materials/RheologyElasticity.hh" // USES RheologyElasticity
#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // pythia::journal

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::mmstests::TestFaultKin::setUp(void) {
    MMSTest::setUp();

    _fault = new pylith::faults::FaultCohesiveKin;CPPUNIT_ASSERT(_fault);
    _data = NULL;
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate testing data.
void
pylith::mmstests::TestFaultKin::tearDown(void) {
    for (size_t i = 0; i < _materials.size(); ++i) {
        delete _materials[i];_materials[i] = NULL;
    } // for
    delete _fault;_fault = NULL;
    for (size_t i = 0; i < _bcs.size(); ++i) {
        delete _bcs[i];_bcs[i] = NULL;
    } // for
    delete _data;_data = NULL;

    MMSTest::tearDown();
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
// Initialize objects for test.
void
pylith::mmstests::TestFaultKin::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    pylith::meshio::MeshIOAscii iohandler;
    CPPUNIT_ASSERT(_data->meshFilename);
    iohandler.filename(_data->meshFilename);
    iohandler.read(_mesh);CPPUNIT_ASSERT(_mesh);

    CPPUNIT_ASSERT_MESSAGE("Test mesh does not contain any cells.",
                           pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    CPPUNIT_ASSERT_MESSAGE("Test mesh does not contain any vertices.",
                           pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

    // Set up coordinates.
    _mesh->setCoordSys(_data->cs);
    CPPUNIT_ASSERT(_data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->normalizer);

    // Set up materials
    for (size_t iMat = 0; iMat < _materials.size(); ++iMat) {
        CPPUNIT_ASSERT(_materials[iMat]);
        _materials[iMat]->setAuxiliaryFieldDB(_data->matAuxDB);

        for (int iSubfield = 0; iSubfield < _data->matNumAuxSubfields; ++iSubfield) {
            const pylith::topology::FieldBase::Discretization& info = _data->matAuxDiscretizations[iSubfield];
            _materials[iMat]->setAuxiliarySubfieldDiscretization(_data->matAuxSubfields[iSubfield], info.basisOrder, info.quadOrder,
                                                                 _data->spaceDim, info.cellBasis, info.feSpace, info.isBasisContinuous);
        } // for
    } // for

    // Set up fault
    CPPUNIT_ASSERT(_fault);
    CPPUNIT_ASSERT(_data->kinsrc);
    _fault->adjustTopology(_mesh);
    const int numRuptures = 1;
    const char* ruptureNames[1] = { "rupture" };
    pylith::faults::KinSrc* ruptures[1] = { _data->kinsrc };
    _fault->setEqRuptures(ruptureNames, numRuptures, ruptures, numRuptures);
    _data->kinsrc->auxFieldDB(_data->faultAuxDB);
    for (int i = 0; i < _data->faultNumAuxSubfields; ++i) {
        const pylith::topology::FieldBase::Discretization& info = _data->faultAuxDiscretizations[i];
        _fault->setAuxiliarySubfieldDiscretization(_data->faultAuxSubfields[i], info.basisOrder, info.quadOrder,
                                                   _data->spaceDim-1, info.cellBasis, info.feSpace, info.isBasisContinuous);
    } // for

    // Set up problem.
    CPPUNIT_ASSERT(_problem);
    CPPUNIT_ASSERT(_data->normalizer);
    _problem->setNormalizer(*_data->normalizer);
    _problem->setGravityField(_data->gravityField);
    _problem->setMaterials(&_materials[0], _materials.size());
    pylith::faults::FaultCohesive* ics[1] = { _fault };
    _problem->setInterfaces(ics, 1);
    _problem->setBoundaryConditions(&_bcs[0], _bcs.size());
    _problem->setStartTime(_data->startTime);
    _problem->setEndTime(_data->endTime);
    _problem->setInitialTimeStep(_data->timeStep);

    // Set up solution field.
    CPPUNIT_ASSERT( (!_data->isExplicit && 2 == _data->numSolnSubfields) ||
                    (_data->isExplicit && 3 == _data->numSolnSubfields) );
    CPPUNIT_ASSERT(_data->solnDiscretizations);

    CPPUNIT_ASSERT(!_solution);
    _solution = new pylith::topology::Field(*_mesh);CPPUNIT_ASSERT(_solution);
    _solution->setLabel("solution");
    pylith::problems::SolutionFactory factory(*_solution, *_data->normalizer);
    int iField = 0;
    factory.addDisplacement(_data->solnDiscretizations[iField++]);
    if (_data->isExplicit) {
        factory.addVelocity(_data->solnDiscretizations[iField++]);
    } // if
    factory.addLagrangeMultiplierFault(_data->solnDiscretizations[iField++]);
    _problem->setSolution(_solution);

    pylith::testing::MMSTest::_initialize();

    PYLITH_METHOD_END;
} // _initialize


// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::mmstests::TestFaultKin_Data::TestFaultKin_Data(void) :
    spaceDim(0),
    meshFilename(NULL),
    cs(NULL),
    gravityField(NULL),
    normalizer(new spatialdata::units::Nondimensional),

    startTime(0.0),
    endTime(0.0),
    timeStep(0.0),
    isExplicit(false),

    numSolnSubfields(0),
    solnDiscretizations(NULL),

    rheology(NULL),
    matNumAuxSubfields(0),
    matAuxSubfields(NULL),
    matAuxDiscretizations(NULL),
    matAuxDB(new spatialdata::spatialdb::UserFunctionDB),

    kinsrc(NULL),
    faultNumAuxSubfields(0),
    faultAuxSubfields(NULL),
    faultAuxDiscretizations(NULL),
    faultAuxDB(new spatialdata::spatialdb::UserFunctionDB) {
    CPPUNIT_ASSERT(normalizer);

    CPPUNIT_ASSERT(matAuxDB);
    matAuxDB->setLabel("material auxiliary field spatial database");

    CPPUNIT_ASSERT(faultAuxDB);
    faultAuxDB->setLabel("fault auxiliary field spatial database");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::mmstests::TestFaultKin_Data::~TestFaultKin_Data(void) {
    delete cs;cs = NULL;
    delete gravityField;gravityField = NULL;
    delete normalizer;normalizer = NULL;
    delete rheology;rheology = NULL;
    delete matAuxDB;matAuxDB = NULL;
    delete kinsrc;kinsrc = NULL;
    delete faultAuxDB;faultAuxDB = NULL;
} // destructor


// End of file
