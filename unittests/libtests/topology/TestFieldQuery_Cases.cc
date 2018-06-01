// -*- C++ -*-
//
// -----------------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
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

        class TestFieldQuery_Quad;
    } // topologyfwd
} // pylith


// ---------------------------------------------------------------------
class pylith::topology::TestFieldQuery_DispTemp : public pylith::topology::TestFieldQuery {

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
        TestFieldQuery::setUp();

        CPPUNIT_ASSERT(_data->normalizer);
        _data->normalizer->lengthScale(1000.0);
        _data->normalizer->timeScale(10.0);
        _data->normalizer->pressureScale(0.1);
        _data->normalizer->densityScale(2.0);

        _data->numAuxSubfields = 2;
        static const char* auxSubfields[2] = { "displacement", "temperature" };
        _data->auxSubfields = const_cast<const char**>(auxSubfields);
        static const char* displacementComponents[2] = { "displacement_x", "displacement_y" };
        static const char* temperatureComponents[1] = {"temperature"};
        static const pylith::topology::Field::Description auxDescriptions[2] = {
            {
                "displacement", // label
                "displacement", // alias
                pylith::topology::Field::VECTOR, // vectorFieldType
                pylith::string_vector(displacementComponents, displacementComponents+2),
                2,
                1000.0,
                NULL,
            },{
                "temperature", // label
                "temperature", // alias
                pylith::topology::Field::SCALAR, // vectorFieldType
                pylith::string_vector(temperatureComponents, temperatureComponents+1),
                1,
                1.0,
                NULL,
            },
        };
        _data->auxDescriptions = const_cast<pylith::topology::Field::Description*>(auxDescriptions);

        static const pylith::topology::Field::Discretization auxDiscretizations[2] = {
            {1, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // displacement
            {2, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // temperature
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(auxDiscretizations);

        CPPUNIT_ASSERT(_data->auxDB);
        _data->auxDB->addValue("displacement_x", disp_x, disp_units());
        _data->auxDB->addValue("displacement_y", disp_y, disp_units());
        _data->auxDB->addValue("temperature", temp, temp_units());
    }   // setUp

}; // TestFieldQuery_DispTemp

// ---------------------------------------------------------------------
class pylith::topology::TestFieldQuery_Quad : public pylith::topology::TestFieldQuery_DispTemp {

    // Spatial database user functions for auxiliary subfields.

protected:
    CPPUNIT_TEST_SUB_SUITE(TestFieldQuery_Quad, TestFieldQuery);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFieldQuery_DispTemp::setUp();

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

        CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(2);
        _data->cs->initialize();

        _data->auxDB->coordsys(*_data->cs);
    }   // setUp


};  // TestFieldQuery_Quad
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::topology::TestFieldQuery_Quad);


// End of file
