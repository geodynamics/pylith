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

#include "GeomDataTet3D.hh"

const int pylith::feassemble::GeomDataTet3D::_cellDim = 3;

const int pylith::feassemble::GeomDataTet3D::_spaceDim = 3;

const PylithScalar pylith::feassemble::GeomDataTet3D::_gravityVec[] = {
  0.0, 0.0, -9.80665 };

const int pylith::feassemble::GeomDataTet3D::_numCorners = 4;

const int pylith::feassemble::GeomDataTet3D::_numLocs = 2;

const PylithScalar pylith::feassemble::GeomDataTet3D::_vertices[] = {
  -1.3, -0.8, 0.2,
  -1.0, 2.4, -0.3,
  2.1, -0.7, 0.1,
  -0.1, 0.2, 3.0
};

const PylithScalar pylith::feassemble::GeomDataTet3D::_locations[] = {
  0.345, 0.397, 0.319,
  0.459, 0.727, 0.693
};

// Reference cell has volume of 2.0, so divide by 2.0;
const PylithScalar pylith::feassemble::GeomDataTet3D::_jacobian[] = {
  3.4/2.0, 0.3/2.0, 1.2/2.0, 0.1/2.0, 3.2/2.0, 1.0/2.0, -0.1/2.0, -0.5/2.0, 2.8/2.0,
  3.4/2.0, 0.3/2.0, 1.2/2.0, 0.1/2.0, 3.2/2.0, 1.0/2.0, -0.1/2.0, -0.5/2.0, 2.8/2.0
};

// Reference cell has volume of 2.0, so divide by 8.0;
const PylithScalar pylith::feassemble::GeomDataTet3D::_jacobianDet[] = {
  32.374/8.0,
  32.374/8.0
};

pylith::feassemble::GeomDataTet3D::GeomDataTet3D(void)
{ // constructor
  cellDim = _cellDim;
  spaceDim = _spaceDim;
  numCorners = _numCorners;
  numLocs = _numLocs;
  gravityVec = const_cast<PylithScalar*>(_gravityVec);
  vertices = const_cast<PylithScalar*>(_vertices);
  locations = const_cast<PylithScalar*>(_locations);
  jacobian = const_cast<PylithScalar*>(_jacobian);
  jacobianDet = const_cast<PylithScalar*>(_jacobianDet);
} // constructor

pylith::feassemble::GeomDataTet3D::~GeomDataTet3D(void)
{}


// End of file
