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

#include "OutputSolnPointsDataTri3.hh"

const char* pylith::meshio::OutputSolnPointsDataTri3::_meshFilename =
  "data/mesh_tri3.exo";

const int pylith::meshio::OutputSolnPointsDataTri3::_spaceDim = 2;
const int pylith::meshio::OutputSolnPointsDataTri3::_numPoints = 9;
const PylithScalar pylith::meshio::OutputSolnPointsDataTri3::_points[9*2] = {
    97.16316825,   421.37654702, // random interior points
    -1563.74874936, -3404.05692152,
     -695.88887916, -1040.46076016,
     -148.63433749, -2994.98857429,
    -3769.33953777,  2270.48242702,
     1440.63801728, -2212.39736902,
    -1162.09140595,  1542.16123633,
    -4000.0, 500.0, // edge
    +4000.0, -4000.0, // vertex
};
const char* pylith::meshio::OutputSolnPointsDataTri3::_names[9] = {
    "AA",
    "BB",
    "CC",
    "DD",
    "EE",
    "FF",
    "GG",
    "HH",
    "II",
};

const int pylith::meshio::OutputSolnPointsDataTri3::_fiberDim = 1;

const PylithScalar pylith::meshio::OutputSolnPointsDataTri3::_coefs[1*2] = {
    10.0, -12.0,
};

// ----------------------------------------------------------------------
pylith::meshio::OutputSolnPointsDataTri3::OutputSolnPointsDataTri3(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  spaceDim = _spaceDim;
  numPoints = _numPoints;
  points = const_cast<PylithScalar*>(_points);
  names = const_cast<const char**>(_names);

  fiberDim = _fiberDim;
  coefs = const_cast<PylithScalar*>(_coefs);

} // constructor

pylith::meshio::OutputSolnPointsDataTri3::~OutputSolnPointsDataTri3(void)
{}


// End of file
