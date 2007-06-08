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

const int pylith::feassemble::GeomDataQuad2D::_cellDim = 2;

const int pylith::feassemble::GeomDataQuad2D::_spaceDim = 2;

const int pylith::feassemble::GeomDataQuad2D::_numCorners = 4;

const int pylith::feassemble::GeomDataQuad2D::_numLocs = 5;

const double pylith::feassemble::GeomDataQuad2D::_vertices[] = {
  0.3, 0.1,
  0.8, -0.2,
  0.7, 1.2,
  -0.1, 1.6
};

const double pylith::feassemble::GeomDataQuad2D::_locations[] = {
  0.0, 0.0,
  1.0, 0.0,
  0.0, 1.0,
  1.0, 1.0,
  0.4, 0.7
};

const double pylith::feassemble::GeomDataQuad2D::_jacobian[] = {
  0.5, -0.4, -0.3, 1.5,
  0.5, -0.1, -0.3, 1.4,
  0.8, -0.4, -0.4, 1.5,
  0.8, -0.1, -0.4, 1.4,
  0.71, -0.28, -0.37, 1.46
};

const double pylith::feassemble::GeomDataQuad2D::_jacobianDet[] = {
  0.63,
  0.67,
  1.04,
  1.08,
  0.933
};

pylith::feassemble::GeomDataQuad2D::GeomDataQuad2D(void)
{ // constructor
  cellDim = _cellDim;
  spaceDim = _spaceDim;
  numCorners = _numCorners;
  numLocs = _numLocs;
  vertices = const_cast<double*>(_vertices);
  locations = const_cast<double*>(_locations);
  jacobian = const_cast<double*>(_jacobian);
  jacobianDet = const_cast<double*>(_jacobianDet);
} // constructor

pylith::feassemble::GeomDataQuad2D::~GeomDataQuad2D(void)
{}


// End of file
