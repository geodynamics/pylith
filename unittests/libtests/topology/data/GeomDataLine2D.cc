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

#include "GeomDataLine2D.hh"

const int pylith::topology::GeomDataLine2D::_cellDim = 1;

const int pylith::topology::GeomDataLine2D::_spaceDim = 2;

const int pylith::topology::GeomDataLine2D::_numCorners = 2;

const int pylith::topology::GeomDataLine2D::_numLocs = 2;

const double pylith::topology::GeomDataLine2D::_vertices[] = {
  1.2, 2.4,
  4.5, -1.4
};

const double pylith::topology::GeomDataLine2D::_locations[] = {
  0.345,
  0.459
};

const double pylith::topology::GeomDataLine2D::_jacobian[] = {
  3.3, -3.8,
  3.3, -3.8
};

pylith::topology::GeomDataLine2D::GeomDataLine2D(void)
{ // constructor
  cellDim = _cellDim;
  spaceDim = _spaceDim;
  numCorners = _numCorners;
  numLocs = _numLocs;
  vertices = const_cast<double*>(_vertices);
  locations = const_cast<double*>(_locations);
  jacobian = const_cast<double*>(_jacobian);
} // constructor

pylith::topology::GeomDataLine2D::~GeomDataLine2D(void)
{}


// End of file
