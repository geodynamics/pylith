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
// Copyright (c) 2010-2015 University of California, Davis
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

const int pylith::meshio::OutputSolnPointsDataHex8::_numVertices = 12;
const int pylith::meshio::OutputSolnPointsDataHex8::_fiberDim = 1;

const PylithScalar pylith::meshio::OutputSolnPointsDataHex8::_field[] = { 
  1.0,
  1.1,
  1.2,
  1.3,
  1.4,
  1.5,
  1.6,
  1.7,
  1.8,
  1.9,
  2.0,
  2.1,
};

const PylithScalar pylith::meshio::OutputSolnPointsDataHex8::_fieldInterp[] = {
  0.0,
  0.0,
  0.0,
  0.0,
};


// ----------------------------------------------------------------------
pylith::meshio::OutputSolnPointsDataHex8::OutputSolnPointsDataHex8(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  spaceDim = _spaceDim;
  numPoints = _numPoints;
  points = const_cast<PylithScalar*>(_points);

  numVertices = _numVertices;
  fiberDim = _fiberDim;
  field = const_cast<PylithScalar*>(_field);

  fieldInterp = const_cast<PylithScalar*>(_fieldInterp);

} // constructor

pylith::meshio::OutputSolnPointsDataHex8::~OutputSolnPointsDataHex8(void)
{}


// End of file
