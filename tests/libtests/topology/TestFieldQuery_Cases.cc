// -*- C++ -*-
//
// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
//

#include <portinfo>

#include "TestFieldQuery.hh" // Implementation of class methods

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// -----------------------------------------------------------------------------
namespace pylith {
    namespace topology {
        class TestFieldQuery_DispTemp;

        class TestFieldQuery_DispTemp2D;
        class TestFieldQuery_Tri;
        class TestFieldQuery_Quad;

        class TestFieldQuery_DispTemp3D;
        class TestFieldQuery_Tet;
        class TestFieldQuery_Hex;
    } // topologyfwd
} // pylith

// ---------------------------------------------------------------------
class pylith::topology::TestFieldQuery_DispTemp : public pylith::topology::TestFieldQuery {
protected:

    void setUp(void) {
        TestFieldQuery::setUp();

        CPPUNIT_ASSERT(_data->normalizer);
        _data->normalizer->setLengthScale(1.0);
        _data->normalizer->setTimeScale(10.0);
        _data->normalizer->setPressureScale(0.1);
        _data->normalizer->setDensityScale(2.0);

        _data->numAuxSubfields = 2;
        static const char* auxSubfields[2] = { "displacement", "temperature" };
        _data->auxSubfields = const_cast<const char**>(auxSubfields);
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
        _data->auxDescriptions = const_cast<pylith::topology::Field::Description*>(auxDescriptions);

        static const pylith::topology::Field::Discretization auxDiscretizations[2] = {
            pylith::topology::Field::Discretization(1, 2), // displacement
            pylith::topology::Field::Discretization(2, 2), // temperature
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(auxDiscretizations);

    } // setUp

}; // TestFieldQuery_DispTemp

// ---------------------------------------------------------------------
class pylith::topology::TestFieldQuery_DispTemp2D : public pylith::topology::TestFieldQuery_DispTemp {
    // Spatial database user functions for auxiliary subfields.
    static const char* disp_units(void) {
        return "m";
    }

    static const char* temp_units(void) {
        return "K";
    }

    static double disp_x(const double x,
                         const double y) {
        return 1.0 + 2.4*x + 3.2*y;
    } // disp_x

    static double disp_y(const double x,
                         const double y) {
        return 0.4 - 0.1*x + 0.6*y;
    } // disp_y

    static double temp(const double x,
                       const double y) {
        return 20.0 + 0.1*x*x + 0.3*x*y -0.2*y*y;
    } // temp

protected:

    void setUp(void) {
        TestFieldQuery_DispTemp::setUp();

        CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(2);

        CPPUNIT_ASSERT(_data->auxDB);
        _data->auxDB->setLabel("Auxiliary db 2D");
        _data->auxDB->setCoordSys(*_data->cs);
        _data->auxDB->addValue("displacement_x", disp_x, disp_units());
        _data->auxDB->addValue("displacement_y", disp_y, disp_units());
        _data->auxDB->addValue("temperature", temp, temp_units());

    } // setUp

}; // TestFieldQuery_DispTemp2D

// ---------------------------------------------------------------------
class pylith::topology::TestFieldQuery_Tri : public pylith::topology::TestFieldQuery_DispTemp2D {
    // Spatial database user functions for auxiliary subfields.

protected:

    CPPUNIT_TEST_SUB_SUITE(TestFieldQuery_Tri, TestFieldQuery_DispTemp2D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFieldQuery_DispTemp2D::setUp();

        // Mesh information.
        _data->cellDim = 2;
        _data->numVertices = 4;
        _data->numCells = 2;
        _data->numCorners = 3;
        static const int _cells[2*3] = {
            0, 1, 2,
            3, 2, 1,
        };
        _data->cells = const_cast<int*>(_cells);
        static const PylithScalar _coordinates[4*2] = {
            -1.0,  0.0,
            +1.0, -2.0,
            +0.8, +1.2,
            +3.0,  0.0,
        };
        _data->coordinates = const_cast<PylithScalar*>(_coordinates);

    } // setUp

}; // TestFieldQuery_Tri
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::topology::TestFieldQuery_Tri);

// ---------------------------------------------------------------------
class pylith::topology::TestFieldQuery_Quad : public pylith::topology::TestFieldQuery_DispTemp2D {
    // Spatial database user functions for auxiliary subfields.

protected:

    CPPUNIT_TEST_SUB_SUITE(TestFieldQuery_Quad, TestFieldQuery_DispTemp2D);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFieldQuery_DispTemp2D::setUp();

        // Mesh information.
        _data->cellDim = 2;
        _data->numVertices = 6;
        _data->numCells = 2;
        _data->numCorners = 4;
        static const int _cells[2*4] = {
            0, 1, 3, 2,
            5, 3, 1, 4,
        };
        _data->cells = const_cast<int*>(_cells);
        static const PylithScalar _coordinates[6*2] = {
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0,
            1.0, 1.0,
            1.8, 0.4,
            1.9, 1.5,
        };
        _data->coordinates = const_cast<PylithScalar*>(_coordinates);

    } // setUp

}; // TestFieldQuery_Quad
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::topology::TestFieldQuery_Quad);

// ---------------------------------------------------------------------
class pylith::topology::TestFieldQuery_DispTemp3D : public pylith::topology::TestFieldQuery_DispTemp {
    // Spatial database user functions for auxiliary subfields.
    static const char* disp_units(void) {
        return "m";
    }

    static const char* temp_units(void) {
        return "K";
    }

    static double disp_x(const double x,
                         const double y,
                         const double z) {
        return 1.0 + 2.4*x + 3.2*y - 0.7*z;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double z) {
        return 0.4 - 0.1*x + 0.6*y + 0.3*z;
    } // disp_y

    static double temp(const double x,
                       const double y,
                       const double z) {
        return 20.0 + 0.1*x*x + 0.3*x*y -0.2*y*y - 3.0*y*z;
    } // temp

protected:

    void setUp(void) {
        TestFieldQuery_DispTemp::setUp();

        CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(3);

        CPPUNIT_ASSERT(_data->auxDB);
        _data->auxDB->setLabel("Auxiliary db 3D");
        _data->auxDB->setCoordSys(*_data->cs);
        _data->auxDB->addValue("displacement_x", disp_x, disp_units());
        _data->auxDB->addValue("displacement_y", disp_y, disp_units());
        _data->auxDB->addValue("temperature", temp, temp_units());

    } // setUp

}; // TestFieldQuery_DispTemp3D

// ---------------------------------------------------------------------
class pylith::topology::TestFieldQuery_Tet : public pylith::topology::TestFieldQuery_DispTemp3D {
    // Spatial database user functions for auxiliary subfields.

protected:

    CPPUNIT_TEST_SUB_SUITE(TestFieldQuery_Tet, TestFieldQuery);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFieldQuery_DispTemp3D::setUp();

        // Mesh information.
        _data->cellDim = 3;
        _data->numVertices = 5;
        _data->numCells = 2;
        _data->numCorners = 4;
        static const int _cells[2*4] = {
            0, 1, 3, 4,
            4, 1, 3, 2,
        };
        _data->cells = const_cast<int*>(_cells);
        static const PylithScalar _coordinates[5*3] = {
            -2.0,  0.0,  0.0,
            +0.0, -1.0,  0.0,
            +2.0,  0.0,  0.0,
            +0.0, +2.0,  0.0,
            +0.0, +1.0, +2.0,
        };
        _data->coordinates = const_cast<PylithScalar*>(_coordinates);

    } // setUp

}; // TestFieldQuery_Tet
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::topology::TestFieldQuery_Tet);

// ---------------------------------------------------------------------
class pylith::topology::TestFieldQuery_Hex : public pylith::topology::TestFieldQuery_DispTemp3D {
    // Spatial database user functions for auxiliary subfields.

protected:

    CPPUNIT_TEST_SUB_SUITE(TestFieldQuery_Hex, TestFieldQuery);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFieldQuery_DispTemp3D::setUp();

        // Mesh information.
        _data->cellDim = 3;
        _data->numVertices = 12;
        _data->numCells = 2;
        _data->numCorners = 8;
        static const int _cells[2*8] = {
            0,  3,  9,  6,  1,  4, 10,  7,
            10, 11,  8,  7,  4,  5,  2,  1,
        };
        _data->cells = const_cast<int*>(_cells);
        static const PylithScalar _coordinates[12*3] = {
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
        _data->coordinates = const_cast<PylithScalar*>(_coordinates);

    } // setUp

}; // TestFieldQuery_Hex
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::topology::TestFieldQuery_Hex);

// End of file
