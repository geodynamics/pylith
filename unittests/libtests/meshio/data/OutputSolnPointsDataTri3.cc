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

#include "OutputSolnPointsDataTri3.hh"

const char* pylith::meshio::OutputSolnPointsDataTri3::_meshFilename =
  "data/tri3.mesh";

const int pylith::meshio::OutputSolnPointsDataTri3::_spaceDim = 2;
const int pylith::meshio::OutputSolnPointsDataTri3::_numPoints = 5;
const PylithScalar pylith::meshio::OutputSolnPointsDataTri3::_points[] = {
   0.1, 0.4, // interior points
  -0.2, 0.3,
   0.3, 0.1,
   0.001, 0.75, // edge
   0.999, 0.0, // vertex
};

const int pylith::meshio::OutputSolnPointsDataTri3::_numVertices = 4;
const int pylith::meshio::OutputSolnPointsDataTri3::_fiberDim = 1;

const PylithScalar pylith::meshio::OutputSolnPointsDataTri3::_field[] = {
  1.0,
  1.1,
  1.2,
  1.3,
};

const PylithScalar pylith::meshio::OutputSolnPointsDataTri3::_fieldInterp[] = {
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
};


// ----------------------------------------------------------------------
pylith::meshio::OutputSolnPointsDataTri3::OutputSolnPointsDataTri3(void)
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

pylith::meshio::OutputSolnPointsDataTri3::~OutputSolnPointsDataTri3(void)
{}


// End of file
