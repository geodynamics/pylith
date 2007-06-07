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

#include "GeomDataPoint1D.hh"

const int pylith::feassemble::GeomDataPoint1D::_cellDim = 0;

const int pylith::feassemble::GeomDataPoint1D::_spaceDim = 1;

const int pylith::feassemble::GeomDataPoint1D::_numCorners = 1;

const int pylith::feassemble::GeomDataPoint1D::_numLocs = 2;

const double pylith::feassemble::GeomDataPoint1D::_vertices[] = {
  1.2,
  4.5
};

const double pylith::feassemble::GeomDataPoint1D::_locations[] = {
  0.0,
  0.0
};

const double pylith::feassemble::GeomDataPoint1D::_jacobian[] = {
  1.0,
  1.0
};

pylith::feassemble::GeomDataPoint1D::GeomDataPoint1D(void)
{ // constructor
  cellDim = _cellDim;
  spaceDim = _spaceDim;
  numCorners = _numCorners;
  numLocs = _numLocs;
  vertices = const_cast<double*>(_vertices);
  locations = const_cast<double*>(_locations);
  jacobian = const_cast<double*>(_jacobian);
} // constructor

pylith::feassemble::GeomDataPoint1D::~GeomDataPoint1D(void)
{}


// End of file
