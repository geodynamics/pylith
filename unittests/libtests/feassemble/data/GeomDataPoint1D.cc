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

#include "GeomDataPoint1D.hh"

const int pylith::feassemble::GeomDataPoint1D::_cellDim = 0;

const int pylith::feassemble::GeomDataPoint1D::_spaceDim = 1;

const PylithScalar pylith::feassemble::GeomDataPoint1D::_gravityVec[] = {
  -9.80665 };

const int pylith::feassemble::GeomDataPoint1D::_numCorners = 1;

const int pylith::feassemble::GeomDataPoint1D::_numLocs = 2;

const PylithScalar pylith::feassemble::GeomDataPoint1D::_vertices[] = {
  1.2,
  4.5
};

const PylithScalar pylith::feassemble::GeomDataPoint1D::_locations[] = {
  0.0,
  0.0
};

const PylithScalar pylith::feassemble::GeomDataPoint1D::_jacobian[] = {
  1.0,
  1.0
};

const PylithScalar pylith::feassemble::GeomDataPoint1D::_jacobianDet[] = {
  1.0,
  1.0
};

pylith::feassemble::GeomDataPoint1D::GeomDataPoint1D(void)
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

pylith::feassemble::GeomDataPoint1D::~GeomDataPoint1D(void)
{}


// End of file
