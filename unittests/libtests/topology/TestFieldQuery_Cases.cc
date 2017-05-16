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

// -----------------------------------------------------------------------------
namespace pylith {
    namespace topology {

        // ---------------------------------------------------------------------
        class TestFieldQuery_Quad : public TestFieldQuery {

            CPPUNIT_TEST_SUB_SUITE(TestFieldQuery_Quad, TestFieldQuery);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestFieldQuery::setUp();

                // Mesh information.
                _data->cellDim = 2;
                _data->numVertices = 4;
                _data->numCells = 1;
                _data->numCorners = 4;
                static const int _cells[1*4] = {
                    0, 1, 2, 3,
                };
                _data->cells = const_cast<int*>(_cells);
                static const PylithScalar _coordinates[4*2] = {
                    0.0, 0.0,
                    1.0, 0.0,
                    0.0, 1.0,
                    1.0, 1.0,
                };
                _data->coordinates = const_cast<PylithScalar*>(_coordinates);

                // Subfield A
                _data->descriptionA.label = "displacement";
                _data->descriptionA.vectorFieldType = FieldBase::VECTOR;
                _data->descriptionA.scale = 2.0;
                _data->descriptionA.numComponents = 2;
                _data->descriptionA.componentNames.resize(2);
                _data->descriptionA.componentNames[0] = "displacement_x";
                _data->descriptionA.componentNames[1] = "displacement_y";
                _data->descriptionA.validator = NULL;

                _data->discretizationA.basisOrder = 1;
                _data->discretizationA.quadOrder = 1;
                _data->discretizationA.isBasisContinuous = true;
                _data->discretizationA.feSpace = pylith::topology::FieldBase::POLYNOMIAL_SPACE;
                static const PylithScalar _subfieldAValues[4*2] = {
                    1.1, 1.2,
                    2.1, 2.2,
                    3.1, 3.2,
                    4.1, 4.2,
                };
                _data->subfieldAValues = const_cast<PylithScalar*>(_subfieldAValues);

                // Subfield B
                _data->descriptionB.label = "fluid_pressure";
                _data->descriptionB.vectorFieldType = FieldBase::SCALAR;
                _data->descriptionB.scale = 0.1;
                _data->descriptionB.numComponents = 1;
                _data->descriptionB.componentNames.resize(1);
                _data->descriptionB.componentNames[0] = "fluid_pressure";
                _data->descriptionB.validator = NULL;

                _data->discretizationB.basisOrder = 1;
                _data->discretizationB.quadOrder = 1;
                _data->discretizationB.isBasisContinuous = true;
                _data->discretizationB.feSpace = pylith::topology::FieldBase::POLYNOMIAL_SPACE;
                static const PylithScalar _subfieldBValues[4*1] = {
                    1.3,
                    2.3,
                    3.3,
                    4.3,
                };
                _data->subfieldBValues = const_cast<PylithScalar*>(_subfieldBValues);
            }   // setUp


        };  // TestFieldQuery_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestFieldQuery_Quad);

    }   // topology
}   // pylith


// End of file
