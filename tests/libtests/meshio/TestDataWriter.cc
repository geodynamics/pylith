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

#include "TestDataWriter.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::TestDataWriter_Data::TestDataWriter_Data(void) :
    meshFilename(NULL),
    faultLabel(NULL),
    faultId(-1),
    spaceDim(0),
    lengthScale(1.0),
    time(0.0),
    timeFormat(NULL),
    vertexNumPoints(0),
    vertexValues(NULL),
    vertexNumDOF(0),
    cellNumPoints(0),
    cellValues(NULL),
    cellNumDOF(0) {}


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::TestDataWriter_Data::~TestDataWriter_Data(void) {}


// End of file
