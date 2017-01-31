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

#include "OutputSolnPointsDataQuad4.hh"

const char* pylith::meshio::OutputSolnPointsDataQuad4::_meshFilename =
  "data/quad4.mesh";

const int pylith::meshio::OutputSolnPointsDataQuad4::_spaceDim = 2;
const int pylith::meshio::OutputSolnPointsDataQuad4::_numPoints = 7;
const PylithScalar pylith::meshio::OutputSolnPointsDataQuad4::_points[] = {
    0.0,  0.1, // interior points
    0.3,  0.4,
   -0.6, -0.7,
   -1.0,  0.9,
   -0.3,  0.8,
    0.3,  0.99999, // edge point
   -0.999999,  0.99999, // vertex point
};
const char* pylith::meshio::OutputSolnPointsDataQuad4::_names[7] = {
    "AA",
    "BB",
    "CC",
    "DD",
    "EE",
    "FF",
    "GG",
};

const int pylith::meshio::OutputSolnPointsDataQuad4::_fiberDim = 3;

const PylithScalar pylith::meshio::OutputSolnPointsDataQuad4::_coefs[3*2] = {
    2.0, -3.0, 10.0,
    -1.0, 6.0, 13.0,
};

// ----------------------------------------------------------------------
pylith::meshio::OutputSolnPointsDataQuad4::OutputSolnPointsDataQuad4(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  spaceDim = _spaceDim;
  numPoints = _numPoints;
  points = const_cast<PylithScalar*>(_points);
  names = const_cast<const char**>(_names);

  fiberDim = _fiberDim;
  coefs = const_cast<PylithScalar*>(_coefs);

} // constructor

pylith::meshio::OutputSolnPointsDataQuad4::~OutputSolnPointsDataQuad4(void)
{}


// End of file
