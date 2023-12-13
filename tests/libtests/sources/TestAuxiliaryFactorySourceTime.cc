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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestAuxiliaryFactorySourceTime.hh" // Implementation of class methods

#include "pylith/sources/AuxiliaryFactorySourceTime.hh" // Test subject
#include "tests/src/FieldTester.hh" // USES FieldTester

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::sources::TestAuxiliaryFactorySourceTime::TestAuxiliaryFactorySourceTime(TestAuxiliaryFactorySourceTime_Data* data) :
    _data(data) {
    PYLITH_METHOD_BEGIN;

    assert(_data->normalizer);
    _data->normalizer->setLengthScale(1.0e+03);
    _data->normalizer->setTimeScale(2.0);
    _data->normalizer->setDensityScale(3.0e+3);
    _data->normalizer->setPressureScale(2.25e+10);

    pylith::topology::Field::SubfieldInfo info;
    pylith::string_vector componentNames;

    // momentTensor
    componentNames.resize(6);
    componentNames[0] = "moment_tensor_xx";
    componentNames[1] = "moment_tensor_yy";
    componentNames[2] = "moment_tensor_zz";
    componentNames[3] = "moment_tensor_xy";
    componentNames[4] = "moment_tensor_yz";
    componentNames[5] = "moment_tensor_xz";

    info.description = pylith::topology::Field::Description(
        "moment_tensor",
        "moment_tensor",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::TENSOR,
        _data->normalizer->getPressureScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, _data->auxDim, _data->auxDim, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, false
        );
    info.index = 0;
    _data->subfields["moment_tensor"] = info;
    if (2 == _data->auxDim) {
        _data->subfields["moment_tensor"].description.numComponents = 4;
        _data->subfields["moment_tensor"].description.vectorFieldType = pylith::topology::Field::OTHER;
    } // if

    // timeDelay
    componentNames.resize(1);
    componentNames[0] = "time_delay";
    info.description = pylith::topology::Field::Description(
        "time_delay",
        "time_delay",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getTimeScale(),
        pylith::topology::FieldQuery::validatorNonnegative
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _data->auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 1;
    _data->subfields["time_delay"] = info;

    // centerFrequency
    componentNames.resize(1);
    componentNames[0] = "center_frequency";
    info.description = pylith::topology::Field::Description(
        "center_frequency",
        "center_frequency",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        1.0 / _data->normalizer->getTimeScale(),
        pylith::topology::FieldQuery::validatorNonnegative
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, _data->auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 2;
    _data->subfields["center_frequency"] = info;

    _initialize();

    PYLITH_METHOD_END;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::TestAuxiliaryFactorySourceTime::~TestAuxiliaryFactorySourceTime(void) {
    PYLITH_METHOD_BEGIN;

    delete _factory;_factory = NULL;
    delete _data;_data = NULL;
    delete _mesh;_mesh = NULL;
    delete _auxiliaryField;_auxiliaryField = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test adding moment tensor and time delay subfields.
void
pylith::sources::TestAuxiliaryFactorySourceTime::testAdd(void) {
    PYLITH_METHOD_BEGIN;
    assert(_factory);
    assert(_data);

    CHECK(!_auxiliaryField->hasSubfield("moment_tensor"));
    CHECK(!_auxiliaryField->hasSubfield("time_delay"));
    CHECK(!_auxiliaryField->hasSubfield("center_frequency"));

    _factory->addMomentTensor();
    _factory->addTimeDelay();
    _factory->addCenterFrequency();

    assert(_data->normalizer);

    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["moment_tensor"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["time_delay"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["center_frequency"]);

    PYLITH_METHOD_END;
} // testAdd


// ------------------------------------------------------------------------------------------------
// Test setValues().
void
pylith::sources::TestAuxiliaryFactorySourceTime::testSetValuesFromDB(void) {
    PYLITH_METHOD_BEGIN;

    assert(_factory);

    _factory->addMomentTensor();
    _factory->addTimeDelay();
    _factory->addCenterFrequency();
    _auxiliaryField->subfieldsSetup();
    _auxiliaryField->createDiscretization();
    _auxiliaryField->allocate();

    assert(_data);
    assert(_data->normalizer);
    _factory->setValuesFromDB();
    pylith::testing::FieldTester::checkFieldWithDB(*_auxiliaryField, _data->auxiliaryDB, _data->normalizer->getLengthScale());

    PYLITH_METHOD_END;
} // testSetValues


// ------------------------------------------------------------------------------------------------
// Initialze mesh, coordinate system, auxiliary field, and factory.
void
pylith::sources::TestAuxiliaryFactorySourceTime::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);

    pylith::meshio::MeshIOAscii iohandler;
    assert(_data->meshFilename);
    iohandler.setFilename(_data->meshFilename);
    _mesh = new pylith::topology::Mesh();assert(_mesh);
    iohandler.read(_mesh);

    assert(pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    assert(pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

    // Setup coordinates.
    _mesh->setCoordSys(_data->cs);
    assert(_data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->normalizer);

    _auxiliaryField = new pylith::topology::Field(*_mesh);assert(_auxiliaryField);
    _auxiliaryField->setLabel("auxiliary");

    _factory = new AuxiliaryFactorySourceTime();
    assert(_data->auxiliaryDB);
    _factory->setQueryDB(_data->auxiliaryDB);
    typedef std::map<std::string, pylith::topology::Field::SubfieldInfo>::const_iterator subfield_iter;
    for (subfield_iter iter = _data->subfields.begin(); iter != _data->subfields.end(); ++iter) {
        const char* subfieldName = iter->first.c_str();
        const pylith::topology::Field::Discretization& fe = iter->second.fe;
        _factory->setSubfieldDiscretization(subfieldName, fe.basisOrder, fe.quadOrder, fe.dimension, fe.isFaultOnly, fe.cellBasis,
                                            fe.feSpace, fe.isBasisContinuous);
    } // for
    assert(_data->normalizer);
    _factory->initialize(_auxiliaryField, *_data->normalizer, _data->dimension);

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
pylith::sources::TestAuxiliaryFactorySourceTime_Data::TestAuxiliaryFactorySourceTime_Data(void) :
    meshFilename(NULL),
    cs(NULL),
    normalizer(new spatialdata::units::Nondimensional),
    auxiliaryDB(new spatialdata::spatialdb::UserFunctionDB) {}


// ------------------------------------------------------------------------------------------------
pylith::sources::TestAuxiliaryFactorySourceTime_Data::~TestAuxiliaryFactorySourceTime_Data(void) {
    delete cs;cs = NULL;
    delete normalizer;normalizer = NULL;
    delete auxiliaryDB;auxiliaryDB = NULL;
} // destructor
