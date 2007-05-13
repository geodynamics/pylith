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

#include "GeomDataLine3D.hh"

const int pylith::topology::GeomDataLine3D::_cellDim = 1;

const int pylith::topology::GeomDataLine3D::_spaceDim = 3;

const int pylith::topology::GeomDataLine3D::_numCorners = 2;

const int pylith::topology::GeomDataLine3D::_numLocs = 2;

const double pylith::topology::GeomDataLine3D::_vertices[] = {
  1.2, 2.3, 2.4,
  4.5, -1.5, 5.3
};

const double pylith::topology::GeomDataLine3D::_locations[] = {
  0.345,
  0.459
};

const double pylith::topology::GeomDataLine3D::_jacobian[] = {
  3.3, -3.8, 2.9,
  3.3, -3.8, 2.9
};

pylith::topology::GeomDataLine3D::GeomDataLine3D(void)
{ // constructor
  cellDim = _cellDim;
  spaceDim = _spaceDim;
  numCorners = _numCorners;
  numLocs = _numLocs;
  vertices = const_cast<double*>(_vertices);
  locations = const_cast<double*>(_locations);
  jacobian = const_cast<double*>(_jacobian);
} // constructor

pylith::topology::GeomDataLine3D::~GeomDataLine3D(void)
{}


// End of file
