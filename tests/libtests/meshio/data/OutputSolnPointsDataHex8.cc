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

#include "OutputSolnPointsDataHex8.hh"

const char* pylith::meshio::OutputSolnPointsDataHex8::_meshFilename =
  "data/mesh_hex8.exo";

const int pylith::meshio::OutputSolnPointsDataHex8::_spaceDim = 3;
const int pylith::meshio::OutputSolnPointsDataHex8::_numPoints = 12;
const PylithScalar pylith::meshio::OutputSolnPointsDataHex8::_points[12*3] = {
    // random interior points
    19689.91631843,  18343.86802356, -35239.73155509,
    -7018.54864543,  -3408.40367508, -27864.2279251,
    -9932.13792313,  17125.97147268, -1059.89238979,
    -8400.92028815, -33469.44233322, -6750.00853956,
    -22055.34076694, -11371.10294793, -17933.1315953,
    3843.82866346,   -709.66642423, -11472.36248362,
    -4587.37956464, -31694.65154831, -34577.58065838,
    27845.04819411,   7575.97018408, -27891.47911585,
    -8533.19891649, -36938.5637657, -10686.81635502,
    30721.80825612,  -3822.81946493, -3741.87309063,
    -40000.0, +40000.0, -2300.0, // edge
    +40000.0, -40000.0, +0.0, // vertex
};
const char* pylith::meshio::OutputSolnPointsDataHex8::_names[12] = {
    "AA",
    "BB",
    "CC",
    "DD",
    "EE",
    "FF",
    "GG",
    "HH",
    "II",
    "JJ",
    "KK",
    "LL",
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
