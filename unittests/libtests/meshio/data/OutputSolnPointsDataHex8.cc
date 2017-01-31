// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include "OutputSolnPointsDataHex8.hh"

const char* pylith::meshio::OutputSolnPointsDataHex8::_meshFilename =
  "data/hex8.mesh";

const int pylith::meshio::OutputSolnPointsDataHex8::_spaceDim = 3;
const int pylith::meshio::OutputSolnPointsDataHex8::_numPoints = 4;
const PylithScalar pylith::meshio::OutputSolnPointsDataHex8::_points[12] = {
  0.1, 0.4, 0.2,// interior points
  0.3, 0.1, 0.8,
  0.001, 0.75, -0.0,// edge
  0.999, 0.0, 0.0, // vertex
};
const char* pylith::meshio::OutputSolnPointsDataHex8::_names[4] = {
    "AA",
    "BB",
    "CC",
    "DD",
};

const int pylith::meshio::OutputSolnPointsDataHex8::_fiberDim = 1;

const PylithScalar pylith::meshio::OutputSolnPointsDataHex8::_coefs[1*3] = { 
    -3.0, 6.0, 12.0,
};

// ----------------------------------------------------------------------
pylith::meshio::OutputSolnPointsDataHex8::OutputSolnPointsDataHex8(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  spaceDim = _spaceDim;
  numPoints = _numPoints;
  points = const_cast<PylithScalar*>(_points);
  names = const_cast<const char**>(_names);
  
  fiberDim = _fiberDim;
  coefs = const_cast<PylithScalar*>(_coefs);

} // constructor

pylith::meshio::OutputSolnPointsDataHex8::~OutputSolnPointsDataHex8(void)
{}


// End of file
