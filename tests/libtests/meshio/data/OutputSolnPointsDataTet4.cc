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

#include "OutputSolnPointsDataTet4.hh"

const char* pylith::meshio::OutputSolnPointsDataTet4::_meshFilename =
  "data/mesh_tet4.exo";

const int pylith::meshio::OutputSolnPointsDataTet4::_spaceDim = 3;
const int pylith::meshio::OutputSolnPointsDataTet4::_numPoints = 10;
const PylithScalar pylith::meshio::OutputSolnPointsDataTet4::_points[10*3] = {
    // random interior points
    -6000.15881915,  -3580.23292515, -31990.60993758,
    -13360.21306575, -15520.37170445, -10060.94019049,
    -13700.63488871, -22720.9494787, -3190.61843645,
    15800.48184977,  17340.62982775, -2650.08094591,
    15750.68552586,  24000.35915383, -16020.69611792,
    -7060.50564414,   9380.7932739, -27110.94878431,
    -17110.06012685,  32450.33022291, -27920.33736877,
    -22060.4290158 ,  38590.77984257, -11200.16910736,
    -40000.0, +40000.0, -2300.0, // edge
    +40000.0, -40000.0, +0.0, // vertex
};
const char* pylith::meshio::OutputSolnPointsDataTet4::_names[10] = {
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
};

const int pylith::meshio::OutputSolnPointsDataTet4::_fiberDim = 2;

const PylithScalar pylith::meshio::OutputSolnPointsDataTet4::_coefs[2*3] = {
    2.0, -3.0, 10.0,
    5.0, 5.0, -8.0,
};

// ----------------------------------------------------------------------
pylith::meshio::OutputSolnPointsDataTet4::OutputSolnPointsDataTet4(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  spaceDim = _spaceDim;
  numPoints = _numPoints;
  points = const_cast<PylithScalar*>(_points);
  names = const_cast<const char**>(_names);

  fiberDim = _fiberDim;
  coefs = const_cast<PylithScalar*>(_coefs);

} // constructor

pylith::meshio::OutputSolnPointsDataTet4::~OutputSolnPointsDataTet4(void)
{}


// End of file
