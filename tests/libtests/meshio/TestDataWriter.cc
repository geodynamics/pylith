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
    vertexScalarValues(NULL),
    vertexVectorValues(NULL),
    vertexTensorValues(NULL),
    vertexOtherValues(NULL),
    vertexScalarNumComponents(0),
    vertexVectorNumComponents(0),
    vertexTensorNumComponents(0),
    vertexOtherNumComponents(0),
    cellNumPoints(0),
    cellScalarValues(NULL),
    cellVectorValues(NULL),
    cellTensorValues(NULL),
    cellOtherValues(NULL),
    cellScalarNumComponents(0),
    cellVectorNumComponents(0),
    cellTensorNumComponents(0),
    cellOtherNumComponents(0)
{ // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::TestDataWriter_Data::~TestDataWriter_Data(void)
{ // destructor
} // destructor


// End of file
