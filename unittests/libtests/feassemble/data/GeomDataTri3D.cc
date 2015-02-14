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

#include "GeomDataTri3D.hh"

const int pylith::feassemble::GeomDataTri3D::_cellDim = 2;

const int pylith::feassemble::GeomDataTri3D::_spaceDim = 3;

const PylithScalar pylith::feassemble::GeomDataTri3D::_gravityVec[] = {
  0.0, 0.0, -9.80665 };

const int pylith::feassemble::GeomDataTri3D::_numCorners = 3;

const int pylith::feassemble::GeomDataTri3D::_numLocs = 2;

const PylithScalar pylith::feassemble::GeomDataTri3D::_vertices[] = {
  1.2, 1.3, -0.4,
  -0.3, -1.7, 0.1,
  -0.1, 0.4, 1.4
};

const PylithScalar pylith::feassemble::GeomDataTri3D::_locations[] = {
  0.345, 0.397,
  0.459, 0.727
};

// Reference cell has area of 2.0, so divide by 2.0;
const PylithScalar pylith::feassemble::GeomDataTri3D::_jacobian[] = {
  -1.5/2.0, -1.3/2.0, -3.0/2.0, -0.9/2.0, 0.5/2.0, 1.8/2.0,
  -1.5/2.0, -1.3/2.0, -3.0/2.0, -0.9/2.0, 0.5/2.0, 1.8/2.0
};

// Reference cell has area of 2.0, so divide by 4.0;
const PylithScalar pylith::feassemble::GeomDataTri3D::_jacobianDet[] = {
  5.933590818383081/4.0,
  5.933590818383081/4.0
};

pylith::feassemble::GeomDataTri3D::GeomDataTri3D(void)
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

pylith::feassemble::GeomDataTri3D::~GeomDataTri3D(void)
{}


// End of file
