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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "GeomDataPoint3D.hh"

const int pylith::feassemble::GeomDataPoint3D::_cellDim = 0;

const int pylith::feassemble::GeomDataPoint3D::_spaceDim = 3;

const double pylith::feassemble::GeomDataPoint3D::_gravityVec[] = {
  0.0, 0.0, -9.80665 };

const int pylith::feassemble::GeomDataPoint3D::_numCorners = 1;

const int pylith::feassemble::GeomDataPoint3D::_numLocs = 2;

const double pylith::feassemble::GeomDataPoint3D::_vertices[] = {
  1.22, 4.35, 6.56,
  4.45, 5.62, 2.55
};

const double pylith::feassemble::GeomDataPoint3D::_locations[] = {
  0.0,
  0.0
};

const double pylith::feassemble::GeomDataPoint3D::_jacobian[] = {
  1.0,
  1.0
};

const double pylith::feassemble::GeomDataPoint3D::_jacobianDet[] = {
  1.0,
  1.0
};

pylith::feassemble::GeomDataPoint3D::GeomDataPoint3D(void)
{ // constructor
  cellDim = _cellDim;
  spaceDim = _spaceDim;
  numCorners = _numCorners;
  numLocs = _numLocs;
  gravityVec = const_cast<double*>(_gravityVec);
  vertices = const_cast<double*>(_vertices);
  locations = const_cast<double*>(_locations);
  jacobian = const_cast<double*>(_jacobian);
  jacobianDet = const_cast<double*>(_jacobianDet);
} // constructor

pylith::feassemble::GeomDataPoint3D::~GeomDataPoint3D(void)
{}


// End of file
