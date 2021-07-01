// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include "OutputSolnPointsData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnPointsData::OutputSolnPointsData(void) :
    meshFilename(0),
    spaceDim(0),
    numPoints(0),
    points(0),
    names(0),
    fiberDim(0),
    coefs(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSolnPointsData::~OutputSolnPointsData(void)
{ // destructor
} // destructor


// End of file
