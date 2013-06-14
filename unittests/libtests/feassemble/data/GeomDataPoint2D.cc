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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "GeomDataPoint2D.hh"

const int pylith::feassemble::GeomDataPoint2D::_cellDim = 0;

const int pylith::feassemble::GeomDataPoint2D::_spaceDim = 2;

const PylithScalar pylith::feassemble::GeomDataPoint2D::_gravityVec[] = {
  0.0, -9.80665 };

const int pylith::feassemble::GeomDataPoint2D::_numCorners = 1;

const int pylith::feassemble::GeomDataPoint2D::_numLocs = 2;

const PylithScalar pylith::feassemble::GeomDataPoint2D::_vertices[] = {
  1.3, 5.4,
  4.1, 7.5
};

const PylithScalar pylith::feassemble::GeomDataPoint2D::_locations[] = {
  0.0,
  0.0
};

const PylithScalar pylith::feassemble::GeomDataPoint2D::_jacobian[] = {
  1.0,
  1.0
};

const PylithScalar pylith::feassemble::GeomDataPoint2D::_jacobianDet[] = {
  1.0,
  1.0
};

pylith::feassemble::GeomDataPoint2D::GeomDataPoint2D(void)
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

pylith::feassemble::GeomDataPoint2D::~GeomDataPoint2D(void)
{}


// End of file
