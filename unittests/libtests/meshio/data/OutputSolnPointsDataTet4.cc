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

#include "OutputSolnPointsDataTet4.hh"

const char* pylith::meshio::OutputSolnPointsDataTet4::_meshFilename =
  "data/tet4.mesh";

const int pylith::meshio::OutputSolnPointsDataTet4::_spaceDim = 3;
const int pylith::meshio::OutputSolnPointsDataTet4::_numPoints = 5;
const PylithScalar pylith::meshio::OutputSolnPointsDataTet4::_points[] = {
   0.1, 0.4, 0.1, // interior points
  -0.2, 0.3, 0.2,
   0.3, 0.1, -0.0,
   0.001, 0.75, 0.0001, // edge
   0.999, 0.0, 0.0, // vertex
};

const int pylith::meshio::OutputSolnPointsDataTet4::_numVertices = 5;
const int pylith::meshio::OutputSolnPointsDataTet4::_fiberDim = 2;

const PylithScalar pylith::meshio::OutputSolnPointsDataTet4::_field[] = {
  1.0, 2.5,
  1.1, 2.4,
  1.2, 2.3,
  1.3, 2.2,
  1.4, 2.1,
};

const PylithScalar pylith::meshio::OutputSolnPointsDataTet4::_fieldInterp[] = {
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
};


// ----------------------------------------------------------------------
pylith::meshio::OutputSolnPointsDataTet4::OutputSolnPointsDataTet4(void)
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

pylith::meshio::OutputSolnPointsDataTet4::~OutputSolnPointsDataTet4(void)
{}


// End of file
