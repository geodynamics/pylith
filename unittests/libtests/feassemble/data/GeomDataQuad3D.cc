// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include "GeomDataQuad3D.hh"

const int pylith::feassemble::GeomDataQuad3D::_cellDim = 2;

const int pylith::feassemble::GeomDataQuad3D::_spaceDim = 3;

const int pylith::feassemble::GeomDataQuad3D::_numCorners = 4;

const int pylith::feassemble::GeomDataQuad3D::_numLocs = 5;

const double pylith::feassemble::GeomDataQuad3D::_vertices[] = {
  1.1, -0.9, 0.1,
  1.0, 0.7, 0.8,
  -0.4, 0.6, 1.2,
  -0.3, -1.3, -0.1
};

const double pylith::feassemble::GeomDataQuad3D::_locations[] = {
  0.0, 0.0,
  1.0, 0.0,
  0.0, 1.0,
  1.0, 1.0,
  0.3, 0.8
};

const double pylith::feassemble::GeomDataQuad3D::_jacobian[] = {
  -0.1, -1.4, 1.6, -0.4, 0.7, -0.2,
  -0.1, -1.4, 1.6, -0.1, 0.7, 0.4,
  -0.1, -1.4, 1.9, -0.4, 1.3, -0.2,
  -0.1, -1.4, 1.9, -0.1, 1.3, 0.4,
  -0.1, -1.4, 1.84, -0.31, 1.18, -0.02
};

pylith::feassemble::GeomDataQuad3D::GeomDataQuad3D(void)
{ // constructor
  cellDim = _cellDim;
  spaceDim = _spaceDim;
  numCorners = _numCorners;
  numLocs = _numLocs;
  vertices = const_cast<double*>(_vertices);
  locations = const_cast<double*>(_locations);
  jacobian = const_cast<double*>(_jacobian);
} // constructor

pylith::feassemble::GeomDataQuad3D::~GeomDataQuad3D(void)
{}


// End of file
