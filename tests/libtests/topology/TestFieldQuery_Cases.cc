// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
//

#include <portinfo>

#include "TestFieldQuery.hh" // Implementation of class methods

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "pylith/scales/Scales.hh" // USES Scales

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace topology {
        class TestFieldQuery_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestFieldQuery_Cases {
public:

    // Data factory methods
    static TestFieldQuery_Data* Tri(void);

    static TestFieldQuery_Data* Quad(void);

    static TestFieldQuery_Data* Tet(void);

    static TestFieldQuery_Data* Hex(void);

private:

    static const char* disp_units(void) {
        return "m";
    }

    static const char* temp_units(void) {
        return "K";
    }

    static double disp_2d_x(const double x,
                            const double y) {
        return 1.0 + 2.4*x + 3.2*y;
    } // disp_x

    static double disp_2d_y(const double x,
                            const double y) {
        return 0.4 - 0.1*x + 0.6*y;
    } // disp_y

    static double temp_2d(const double x,
                          const double y) {
        return 20.0 + 0.1*x*x + 0.3*x*y -0.2*y*y;
    } // temp

    static double disp_3d_x(const double x,
                            const double y,
                            const double z) {
        return 1.0 + 2.4*x + 3.2*y - 0.7*z;
    } // disp_x

    static double disp_3d_y(const double x,
                            const double y,
                            const double z) {
        return 0.4 - 0.1*x + 0.6*y + 0.3*z;
    } // disp_y

    static double temp_3d(const double x,
                          const double y,
                          const double z) {
        return 20.0 + 0.1*x*x + 0.3*x*y -0.2*y*y - 3.0*y*z;
    } // temp

    static
    void _setData(TestFieldQuery_Data* data);

    static
    void _setData2D(TestFieldQuery_Data* data);

    static
    void _setData3D(TestFieldQuery_Data* data);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestFieldQuery::Quad::testConstructor", "[TestFieldQuery][Tri][testConstructor]") {
    pylith::topology::TestFieldQuery(pylith::topology::TestFieldQuery_Cases::Tri()).testConstructor();
}
TEST_CASE("TestFieldQuery::Quad::testSetQuery", "[TestFieldQuery][Tri][testSetQuery]") {
    pylith::topology::TestFieldQuery(pylith::topology::TestFieldQuery_Cases::Tri()).testSetQuery();
}
TEST_CASE("TestFieldQuery::Quad::testOpenClose", "[TestFieldQuery][Tri][testOpenClose]") {
    pylith::topology::TestFieldQuery(pylith::topology::TestFieldQuery_Cases::Tri()).testOpenClose();
}
TEST_CASE("TestFieldQuery::Quad::testQuery", "[TestFieldQuery][Tri][testQuery]") {
    pylith::topology::TestFieldQuery(pylith::topology::TestFieldQuery_Cases::Tri()).testQuery();
}
TEST_CASE("TestFieldQuery::Quad::testQueryNull", "[TestFieldQuery][Tri][testQueryNull]") {
    pylith::topology::TestFieldQuery(pylith::topology::TestFieldQuery_Cases::Tri()).testQueryNull();
}
TEST_CASE("TestFieldQuery::Quad::testValidatorPositive", "[TestFieldQuery][Tri][testValidatorPositive]") {
    pylith::topology::TestFieldQuery(pylith::topology::TestFieldQuery_Cases::Tri()).testValidatorPositive();
}
TEST_CASE("TestFieldQuery::Quad::testValidatorNonnegative", "[TestFieldQuery][Tri][testValidatorNonnegative]") {
    pylith::topology::TestFieldQuery(pylith::topology::TestFieldQuery_Cases::Tri()).testValidatorNonnegative();
}

// ------------------------------------------------------------------------------------------------
void
pylith::topology::TestFieldQuery_Cases::_setData(TestFieldQuery_Data* data) {
    assert(data->scales);
    data->scales->setLengthScale(0.02);
    data->scales->setTimeScale(10.0);
    data->scales->setRigidityScale(1.0e+6);

    data->numAuxSubfields = 2;
    static const char* auxSubfields[2] = { "displacement", "temperature" };
    data->auxSubfields = const_cast<const char**>(auxSubfields);
    static const char* displacementComponents[2] = { "displacement_x", "displacement_y" };
    static const char* temperatureComponents[1] = {"temperature"};
    static const pylith::topology::Field::Description auxDescriptions[2] = {
        pylith::topology::Field::Description(
            "displacement", // label
            "displacement", // alias
            pylith::string_vector(displacementComponents, displacementComponents+2),
            2,
            pylith::topology::Field::VECTOR, // vectorFieldType
            1000.0),
        pylith::topology::Field::Description(
            "temperature", // label
            "temperature", // alias
            pylith::string_vector(temperatureComponents, temperatureComponents+1),
            1,
            pylith::topology::Field::SCALAR, // vectorFieldType
            1.0),
    };
    data->auxDescriptions = const_cast<pylith::topology::Field::Description*>(auxDescriptions);

    static const pylith::topology::Field::Discretization auxDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 2), // displacement
        pylith::topology::Field::Discretization(2, 2), // temperature
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(auxDiscretizations);

} // _setData


// ------------------------------------------------------------------------------------------------
void
pylith::topology::TestFieldQuery_Cases::_setData2D(TestFieldQuery_Data* data) {
    _setData(data);
    assert(data->cs);
    data->cs->setSpaceDim(2);

    assert(data->auxDB);
    data->auxDB->setDescription("Auxiliary db 2D");
    data->auxDB->setCoordSys(*data->cs);
    data->auxDB->addValue("displacement_x", disp_2d_x, disp_units());
    data->auxDB->addValue("displacement_y", disp_2d_y, disp_units());
    data->auxDB->addValue("temperature", temp_2d, temp_units());

} // _setData2D


// ------------------------------------------------------------------------------------------------
pylith::topology::TestFieldQuery_Data*
pylith::topology::TestFieldQuery_Cases::Tri(void) {
    TestFieldQuery_Data* data = new TestFieldQuery_Data();assert(data);

    _setData2D(data);

    // Mesh information.
    const size_t numVertices = 4;
    const size_t spaceDim = 2;
    const size_t cellDim = 2;
    const size_t numCells = 2;
    const size_t numCorners = 3;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::TRIANGLE;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        -1.0,  0.0,
        +1.0, -2.0,
        +0.8, +1.2,
        +3.0,  0.0,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0, 1, 2,
        3, 2, 1,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);

    return data;
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::topology::TestFieldQuery_Data*
pylith::topology::TestFieldQuery_Cases::Quad(void) {
    TestFieldQuery_Data* data = new TestFieldQuery_Data();assert(data);

    _setData2D(data);

    // Mesh information.
    const size_t numVertices = 6;
    const size_t spaceDim = 2;
    const size_t cellDim = 2;
    const size_t numCells = 2;
    const size_t numCorners = 4;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::QUADRILATERAL;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        0.0, 0.0,
        1.0, 0.0,
        0.0, 1.0,
        1.0, 1.0,
        1.8, 0.4,
        1.9, 1.5,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0, 1, 3, 2,
        5, 3, 1, 4,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);

    return data;
} // Quad


// ------------------------------------------------------------------------------------------------
void
pylith::topology::TestFieldQuery_Cases::_setData3D(TestFieldQuery_Data* data) {
    _setData(data);

    assert(data->cs);
    data->cs->setSpaceDim(3);

    assert(data->auxDB);
    data->auxDB->setDescription("Auxiliary db 3D");
    data->auxDB->setCoordSys(*data->cs);
    data->auxDB->addValue("displacement_x", disp_3d_x, disp_units());
    data->auxDB->addValue("displacement_y", disp_3d_y, disp_units());
    data->auxDB->addValue("temperature", temp_3d, temp_units());

} // _setData3D


// ------------------------------------------------------------------------------------------------
pylith::topology::TestFieldQuery_Data*
pylith::topology::TestFieldQuery_Cases::Tet(void) {
    TestFieldQuery_Data* data = new TestFieldQuery_Data();assert(data);

    _setData3D(data);

    // Mesh information.
    const size_t numVertices = 5;
    const size_t spaceDim = 3;
    const size_t cellDim = 3;
    const size_t numCells = 2;
    const size_t numCorners = 4;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::TETRAHEDRON;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        -2.0,  0.0,  0.0,
        +0.0, -1.0,  0.0,
        +2.0,  0.0,  0.0,
        +0.0, +2.0,  0.0,
        +0.0, +1.0, +2.0,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0, 1, 3, 4,
        4, 1, 3, 2,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);

    return data;
} // Tet


// ------------------------------------------------------------------------------------------------
pylith::topology::TestFieldQuery_Data*
pylith::topology::TestFieldQuery_Cases::Hex(void) {
    TestFieldQuery_Data* data = new TestFieldQuery_Data();assert(data);

    _setData3D(data);

    // Mesh information.
    const size_t numVertices = 12;
    const size_t spaceDim = 3;
    const size_t cellDim = 3;
    const size_t numCells = 2;
    const size_t numCorners = 8;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::HEXAHEDRON;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        +0.0, +0.0, +0.0,
        +0.0, +0.0, +2.0,
        +0.0, +0.0, +4.5,
        +1.8, +0.1, +0.1,
        +1.9, +0.0, +2.0,
        +1.8, -0.1, +3.9,
        +0.0, +3.0, +0.2,
        +0.1, +2.8, +1.9,
        +0.0, +3.2, +5.0,
        +2.0, +3.0, +0.0,
        +2.1, +2.9, +2.0,
        +2.0, +3.0, +4.0,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0,  3,  9,  6,  1,  4, 10,  7,
        10, 11,  8,  7,  4,  5,  2,  1,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);

    return data;
} // Hex


// End of file
