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

#include "AbsorbingDampersDataLine2.hh"

const char* pylith::bc::AbsorbingDampersDataLine2::_meshFilename = 
  "data/line2.mesh";

const int pylith::bc::AbsorbingDampersDataLine2::_numBasis = 1;
const int pylith::bc::AbsorbingDampersDataLine2::_numQuadPts = 1;
const double pylith::bc::AbsorbingDampersDataLine2::_quadPts[] = {
  0.0,
};
const double pylith::bc::AbsorbingDampersDataLine2::_quadWts[] = {
  1.0,
};
const double pylith::bc::AbsorbingDampersDataLine2::_basis[] = {
  1.0,
};
const double pylith::bc::AbsorbingDampersDataLine2::_basisDerivRef[] = {
  1.0,
};

const char* pylith::bc::AbsorbingDampersDataLine2::_spatialDBFilename = "data/elasticstrain1d.spatialdb";
const int pylith::bc::AbsorbingDampersDataLine2::_id = 2;
const char* pylith::bc::AbsorbingDampersDataLine2::_label = "bc0";

const double pylith::bc::AbsorbingDampersDataLine2::_dt =   1.00000000e-02;
const double pylith::bc::AbsorbingDampersDataLine2::_fieldTpdt[] = {
  1.20000000e+00,
  1.20000000e+00,
  1.70000000e+00,
};
const double pylith::bc::AbsorbingDampersDataLine2::_fieldT[] = {
  1.10000000e+00,
  1.10000000e+00,
  1.50000000e+00,
};
const double pylith::bc::AbsorbingDampersDataLine2::_fieldTmdt[] = {
  1.00000000e+00,
  1.00000000e+00,
  1.30000000e+00,
};

const int pylith::bc::AbsorbingDampersDataLine2::_spaceDim = 1;
const int pylith::bc::AbsorbingDampersDataLine2::_cellDim = 0;
const int pylith::bc::AbsorbingDampersDataLine2::_numVertices = 2;
const int pylith::bc::AbsorbingDampersDataLine2::_numCells = 2;
const int pylith::bc::AbsorbingDampersDataLine2::_numCorners = 1;
const int pylith::bc::AbsorbingDampersDataLine2::_cells[] = {
  2,
  4,
};


const double pylith::bc::AbsorbingDampersDataLine2::_dampingConsts[] = {
  2500.0*6000.0,
  2500.0*6000.0,
};
const double pylith::bc::AbsorbingDampersDataLine2::_valsResidual[] = {
  0.0,
  0.0,
  0.0,
};
const double pylith::bc::AbsorbingDampersDataLine2::_valsJacobian[] = {
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
};

pylith::bc::AbsorbingDampersDataLine2::AbsorbingDampersDataLine2(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);

  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  quadPts = const_cast<double*>(_quadPts);
  quadWts = const_cast<double*>(_quadWts);
  basis = const_cast<double*>(_basis);
  basisDerivRef = const_cast<double*>(_basisDerivRef);

  spatialDBFilename = const_cast<char*>(_spatialDBFilename);
  id = _id;
  label = const_cast<char*>(_label);

  dt = _dt;
  fieldTpdt = const_cast<double*>(_fieldTpdt);
  fieldT = const_cast<double*>(_fieldT);
  fieldTmdt = const_cast<double*>(_fieldTmdt);

  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numVertices = _numVertices;
  numCells = _numCells;
  numCorners = _numCorners;
  cells = const_cast<int*>(_cells);

  dampingConsts = const_cast<double*>(_dampingConsts);
  valsResidual = const_cast<double*>(_valsResidual);
  valsJacobian = const_cast<double*>(_valsJacobian);
} // constructor

pylith::bc::AbsorbingDampersDataLine2::~AbsorbingDampersDataLine2(void)
{}


// End of file
