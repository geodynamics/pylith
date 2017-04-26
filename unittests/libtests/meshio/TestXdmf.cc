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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestXdmf.hh" // Implementation of class methods

#include "pylith/meshio/Xdmf.hh" // USES Xdmf

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include <cstring> // USES strcmp()

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestXdmf::setUp(void)
{ // setUp
    _data = new TestXdmf_Data;
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::meshio::TestXdmf::tearDown(void)
{ // tearDown
    delete _data; _data = NULL;
} // tearDown


// ----------------------------------------------------------------------
// Test constructor.
void
pylith::meshio::TestXdmf::testConstructor(void)
{ // testConstructor
    PYLITH_METHOD_BEGIN;

    Xdmf one;

    PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test write().
void
pylith::meshio::TestXdmf::testWrite(void)
{ // testWrite
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->filenameHDF5);
    CPPUNIT_ASSERT(_data->filenameXdmf);

    Xdmf metafile;
    metafile.write(_data->filenameXdmf, _data->filenameHDF5);

    const std::string filenameE = std::string("data/") + std::string(_data->filenameXdmf);

    std::ifstream fileInE(filenameE.c_str());
    CPPUNIT_ASSERT(fileInE.is_open());

    std::ifstream fileIn(_data->filenameXdmf);
    CPPUNIT_ASSERT(fileIn.is_open());

    const int maxLen = 256;
    char line[maxLen];
    char lineE[maxLen];

    int i = 1;
    while(!fileInE.eof()) {
        fileInE.getline(lineE, maxLen);
        fileIn.getline(line, maxLen);
        if (0 != strcmp(line, lineE)) {
            std::cerr << "Line " << i << " of file '" << _data->filenameXdmf << "' is incorrect." << std::endl;
            CPPUNIT_ASSERT(false);
        } // if
        ++i;
    } // while

    fileInE.close();
    fileIn.close();

    PYLITH_METHOD_END;
} // testWrite


// ----------------------------------------------------------------------
// Constructor
pylith::meshio::TestXdmf_Data::TestXdmf_Data(void) :
    filenameHDF5(NULL),
    filenameXdmf(NULL)
{ // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::TestXdmf_Data::~TestXdmf_Data(void)
{ // destructor
} // destructor


// End of file 
