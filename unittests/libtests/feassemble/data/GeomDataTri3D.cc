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

#include "GeomDataTri3D.hh"

const int pylith::feassemble::GeomDataTri3D::_cellDim = 2;

const int pylith::feassemble::GeomDataTri3D::_spaceDim = 3;

const int pylith::feassemble::GeomDataTri3D::_numCorners = 3;

const int pylith::feassemble::GeomDataTri3D::_numLocs = 2;

const double pylith::feassemble::GeomDataTri3D::_vertices[] = {
  1.2, 1.3, -0.4,
  -0.3, -1.7, 0.1,
  -0.1, 0.4, 1.4
};

const double pylith::feassemble::GeomDataTri3D::_locations[] = {
  0.345, 0.397,
  0.459, 0.727
};

const double pylith::feassemble::GeomDataTri3D::_jacobian[] = {
  -1.5, -1.3, -3.0, -0.9, 0.5, 1.8,
  -1.5, -1.3, -3.0, -0.9, 0.5, 1.8
};

const double pylith::feassemble::GeomDataTri3D::_jacobianDet[] = {
  5.933590818383081,
  5.933590818383081
};

pylith::feassemble::GeomDataTri3D::GeomDataTri3D(void)
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

pylith::feassemble::GeomDataTri3D::~GeomDataTri3D(void)
{}


// End of file
