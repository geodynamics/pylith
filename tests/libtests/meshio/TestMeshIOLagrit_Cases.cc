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

#include "TestMeshIOLagrit.hh" // Implementation of class methods

// ----------------------------------------------------------------------
namespace pylith {
    namespace meshio {

	// --------------------------------------------------------------
	class TestMeshIOLagrit_Tet : public TestMeshIOLagrit {
	    CPPUNIT_TEST_SUB_SUITE(TestMeshIOLagrit_Tet, TestMeshIOLagrit);
	    CPPUNIT_TEST_SUITE_END();

	public :
	    void setUp(void) {
                TestMeshIOLagrit::setUp();
		_data = new TestMeshIOLagrit_Data();CPPUNIT_ASSERT(_data);
		_data->numVertices = 12;
		_data->spaceDim = 3;
		_data->numCells = 12;
		_data->cellDim = 3;
		_data->numCorners = 4;

		static const PylithScalar vertices[12*3] = {
		    0.00000E+000,  -5.00000E-001,  -5.00000E-001,
		    0.00000E+000,  -5.00000E-001,   5.00000E-001,
		    1.00000E+000,  -5.00000E-001,  -5.00000E-001,
		    1.00000E+000,  -5.00000E-001,   5.00000E-001,
		    0.00000E+000,   5.00000E-001,  -5.00000E-001,
		    0.00000E+000,   5.00000E-001,   5.00000E-001,
		    1.00000E+000,   5.00000E-001,  -5.00000E-001,
		    1.00000E+000,   5.00000E-001,   5.00000E-001,
		    -1.00000E+000,  -5.00000E-001,  -5.00000E-001,
		    -1.00000E+000,  -5.00000E-001,   5.00000E-001,
		    -1.00000E+000,   5.00000E-001,  -5.00000E-001,
		    -1.00000E+000,   5.00000E-001,   5.00000E-001
		};
		_data->vertices = const_cast<PylithScalar*>(vertices);

		static const PylithInt cells[12*4] = {
		    9,      5,      8,     10,
		    8,      9,      1,      5,
		    1,      3,      2,      4,
		    8,      5,      4,     10,
		    8,      1,      4,      5,
		    5,      3,      6,      7,
		    4,      3,      2,      6,
		    5,      3,      4,      6,
		    1,      3,      4,      5,
		    9,      5,     10,     11,
		    8,      1,      0,      4,
		    0,      1,      2,      4
		};
		_data->cells = const_cast<PylithInt*>(cells);
		static const PylithInt materialIds[12] = {
		    2, 2, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1
		};
		_data->materialIds = const_cast<PylithInt*>(materialIds);
		
		_data->numGroups = 7;
		static const PylithInt groupSizes[7] = {
		    4, 4, 4, 6, 6, 6, 6,
		};
		_data->groupSizes = const_cast<PylithInt*>(groupSizes);
		static const PylithInt groups[4+4+4+6+6+6+6] = {
		    0,  1,  4,  5,
		    8,  9, 10, 11,
		    2,  3,  6,  7,
		    0,  1,  2,  3,  8,  9,
		    4,  5,  6,  7, 10, 11,
		    0,  2,  4,  6,  8, 10,
		    1,  3,  5,  7,  9, 11
		};
		_data->groups = const_cast<PylithInt*>(groups);
		static const char* groupNames[7] = {
		    "fault",
		    "xm",
		    "xp",
		    "ym",
		    "yp",
		    "zm",
		    "zp"
		};
		_data->groupNames = const_cast<char**>(groupNames);
		static const char* groupTypes[7] = {
		    "vertex",
		    "vertex",
		    "vertex",
		    "vertex",
		    "vertex",
		    "vertex",
		    "vertex",
		};
		_data->groupTypes = const_cast<char**>(groupTypes);
	    } // setUp
	}; // class TestMeshIOLagrit_Tet

        // --------------------------------------------------------------
        class TestMeshIOLagrit_Tet_Ascii : public TestMeshIOLagrit_Tet {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOLagrit_Tet_Ascii, TestMeshIOLagrit_Tet);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOLagrit_Tet::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filenameGmv = "data/cube2_ascii.gmv";
                _data->filenamePset = "data/cube2_ascii.pset";
		_data->ioInt32 = false;
		_data->isRecordHeader32Bit = false;
            } // setUp
        }; // class TestMeshIOLagrit_Tet_Ascii
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOLagrit_Tet_Ascii);

	
        // --------------------------------------------------------------
        class TestMeshIOLagrit_Tet_Binary : public TestMeshIOLagrit_Tet {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOLagrit_Tet_Binary, TestMeshIOLagrit_Tet);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOLagrit_Tet::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filenameGmv = "data/cube2_binary.gmv";
                _data->filenamePset = "data/cube2_binary.pset";
		_data->ioInt32 = true;
		_data->isRecordHeader32Bit = true;
            } // setUp
        }; // class TestMeshIOLagrit_Tet_Ascii
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOLagrit_Tet_Binary);

	
        // --------------------------------------------------------------
        class TestMeshIOLagrit_Tet_Binary32 : public TestMeshIOLagrit_Tet {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOLagrit_Tet_Binary32, TestMeshIOLagrit_Tet);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOLagrit_Tet::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filenameGmv = "data/cube2_binary_32on64.gmv";
                _data->filenamePset = "data/cube2_binary_32on64.pset";
		_data->ioInt32 = true;
		_data->isRecordHeader32Bit = false;
            } // setUp
        }; // class TestMeshIOLagrit_Tet4_Binary32
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOLagrit_Tet_Binary32);

	
        // --------------------------------------------------------------
        class TestMeshIOLagrit_Tet_Binary64 : public TestMeshIOLagrit_Tet {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOLagrit_Tet_Binary64, TestMeshIOLagrit_Tet);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOLagrit_Tet::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filenameGmv = "data/cube2_binary_64.gmv";
                _data->filenamePset = "data/cube2_binary_64.pset";
		_data->ioInt32 = false;
		_data->isRecordHeader32Bit = false;
            } // setUp
        }; // class TestMeshIOLagrit_Tet_Binary64
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOLagrit_Tet_Binary64);

	
    } // meshio
} // pylith


// End of file 
