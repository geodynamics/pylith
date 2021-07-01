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

#include "OutputSolnPointsDataQuad4.hh"

const char* pylith::meshio::OutputSolnPointsDataQuad4::_meshFilename =
  "data/mesh_quad4.exo";

const int pylith::meshio::OutputSolnPointsDataQuad4::_spaceDim = 2;
const int pylith::meshio::OutputSolnPointsDataQuad4::_numPoints = 9;
const PylithScalar pylith::meshio::OutputSolnPointsDataQuad4::_points[9*2] = {
    -2466.89062768, -3338.70224665, // random interior points
    -2573.79512745,  3527.28189631,
    626.05179461, -3749.09327302,
    -3378.87781999,  2950.4481556 ,
    -2791.09686116,   956.10709564,
    1356.65774294, -1018.21897423,
    -3740.56959013,   983.93227394,
    +3999.0, -200.0, // edge
    -4000.0, -4000.0, // vertex
};
const char* pylith::meshio::OutputSolnPointsDataQuad4::_names[9] = {
    "AA",
    "BB",
    "CC",
    "DD",
    "EE",
    "FF",
    "GG",
    "II",
    "JJ",
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
