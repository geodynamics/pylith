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
// Copyright (c) 2010-2016 University of California, Davis
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

            CPPUNIT_TEST_SUB_SUITE( TestFieldQuery_Quad, TestFieldQuery );
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
                _data->subfieldAName = "displacement";
                _data->subfieldAType = FieldBase::VECTOR;
                _data->subfieldAScale = 2.0;
                _data->subfieldANumComponents = 2;
                static const char* _subfieldAComponents[2] = {"displacement_x", "displacement_y"};
                _data->subfieldAComponents = const_cast<const char**>(_subfieldAComponents);
                static const PylithScalar _subfieldAValues[4*2] = {
                    1.1, 1.2,
                    2.1, 2.2,
                    3.1, 3.2,
                    4.1, 4.2,
                };
                _data->subfieldAValues = const_cast<PylithScalar*>(_subfieldAValues);
                _data->subfieldABasisOrder = 1;
                _data->subfieldAQuadOrder = 1;

                // Subfield B
                _data->subfieldBName = "fluid_pressure";
                _data->subfieldAType = FieldBase::SCALAR;
                _data->subfieldBScale = 0.1;
                _data->subfieldBNumComponents = 1;
                static const char* _subfieldBComponents[1] = {"pressure"};
                _data->subfieldBComponents = const_cast<const char**>(_subfieldBComponents);
                static const PylithScalar _subfieldBValues[4*1] = {
                    1.3,
                    2.3,
                    3.3,
                    4.3,
                };
                _data->subfieldBValues = const_cast<PylithScalar*>(_subfieldBValues);
                _data->subfieldBBasisOrder = 1;
                _data->subfieldBQuadOrder = 1;
            }   // setUp


        };  // TestFieldQuery_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION( TestFieldQuery_Quad );

    }   // topology
}   // pylith


// End of file
