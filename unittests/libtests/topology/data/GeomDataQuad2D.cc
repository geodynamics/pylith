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

#include "GeomDataQuad2D.hh"

const int pylith::topology::GeomDataQuad2D::_cellDim = 2;

const int pylith::topology::GeomDataQuad2D::_spaceDim = 2;

const int pylith::topology::GeomDataQuad2D::_numCorners = 4;

const int pylith::topology::GeomDataQuad2D::_numLocs = 5;

const double pylith::topology::GeomDataQuad2D::_vertices[] = {
  0.3, 0.1,
  0.8, -0.2,
  0.7, 1.2,
  -0.1, 1.6
};

const double pylith::topology::GeomDataQuad2D::_locations[] = {
  0.0, 0.0,
  1.0, 0.0,
  0.0, 1.0,
  1.0, 1.0,
  0.4, 0.7
};

const double pylith::topology::GeomDataQuad2D::_jacobian[] = {
  0.5, -0.4, -0.3, 1.5,
  0.5, -0.1, -0.3, 1.4,
  0.8, -0.4, -0.4, 1.5,
  0.8, -0.1, -0.4, 1.4,
  0.71, -0.28, -0.37, 1.46
};

pylith::topology::GeomDataQuad2D::GeomDataQuad2D(void)
{ // constructor
  cellDim = _cellDim;
  spaceDim = _spaceDim;
  numCorners = _numCorners;
  numLocs = _numLocs;
  vertices = const_cast<double*>(_vertices);
  locations = const_cast<double*>(_locations);
  jacobian = const_cast<double*>(_jacobian);
} // constructor

pylith::topology::GeomDataQuad2D::~GeomDataQuad2D(void)
{}


// End of file
