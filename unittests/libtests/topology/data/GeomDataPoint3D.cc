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

#include "GeomDataPoint3D.hh"

const int pylith::topology::GeomDataPoint3D::_cellDim = 0;

const int pylith::topology::GeomDataPoint3D::_spaceDim = 3;

const int pylith::topology::GeomDataPoint3D::_numCorners = 1;

const int pylith::topology::GeomDataPoint3D::_numLocs = 2;

const double pylith::topology::GeomDataPoint3D::_vertices[] = {
  1.22, 4.35, 6.56,
  4.45, 5.62, 2.55
};

const double pylith::topology::GeomDataPoint3D::_locations[] = {
  0.0,
  0.0
};

const double pylith::topology::GeomDataPoint3D::_jacobian[] = {
  1.0,
  1.0
};

pylith::topology::GeomDataPoint3D::GeomDataPoint3D(void)
{ // constructor
  cellDim = _cellDim;
  spaceDim = _spaceDim;
  numCorners = _numCorners;
  numLocs = _numLocs;
  vertices = const_cast<double*>(_vertices);
  locations = const_cast<double*>(_locations);
  jacobian = const_cast<double*>(_jacobian);
} // constructor

pylith::topology::GeomDataPoint3D::~GeomDataPoint3D(void)
{}


// End of file
