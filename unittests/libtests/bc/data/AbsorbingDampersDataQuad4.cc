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

#include "AbsorbingDampersDataQuad4.hh"

const char* pylith::bc::AbsorbingDampersDataQuad4::_meshFilename = 
  "data/quad4.mesh";

const int pylith::bc::AbsorbingDampersDataQuad4::_numBasis = 2;
const int pylith::bc::AbsorbingDampersDataQuad4::_numQuadPts = 1;
const double pylith::bc::AbsorbingDampersDataQuad4::_quadPts[] = {
  0.0,
};
const double pylith::bc::AbsorbingDampersDataQuad4::_quadWts[] = {
  2.0,
};
const double pylith::bc::AbsorbingDampersDataQuad4::_basis[] = {
  0.5,
  0.5,
};
const double pylith::bc::AbsorbingDampersDataQuad4::_basisDerivRef[] = {
  -0.5,
   0.5,
};

const char* pylith::bc::AbsorbingDampersDataQuad4::_spatialDBFilename = "data/elasticplanestrain.spatialdb";
const int pylith::bc::AbsorbingDampersDataQuad4::_id = 2;
const char* pylith::bc::AbsorbingDampersDataQuad4::_label = "bc2";

const double pylith::bc::AbsorbingDampersDataQuad4::_dt =   0.25;
const double pylith::bc::AbsorbingDampersDataQuad4::_fieldTmdt[] = {
  1.0,  2.4,
  1.1,  1.8,
  1.2,  2.4,
  1.3,  2.2,
  1.4,  2.4,
  1.5,  1.6,
};
const double pylith::bc::AbsorbingDampersDataQuad4::_fieldT[] = {
  1.1,  2.0,
  1.3,  2.1,
  1.5,  2.2,
  1.7,  2.3,
  1.9,  2.4,
  2.1,  2.5,
};
const double pylith::bc::AbsorbingDampersDataQuad4::_fieldTIncr[] = {
  1.2,  1.6,
  1.5,  2.4,
  1.8,  2.0,
  2.1,  2.4,
  2.4,  2.4,
  2.7,  3.4,
};

const int pylith::bc::AbsorbingDampersDataQuad4::_spaceDim = 2;
const int pylith::bc::AbsorbingDampersDataQuad4::_cellDim = 1;
const int pylith::bc::AbsorbingDampersDataQuad4::_numVertices = 4;
const int pylith::bc::AbsorbingDampersDataQuad4::_numCells = 2;
const int pylith::bc::AbsorbingDampersDataQuad4::_numCorners = 2;
const int pylith::bc::AbsorbingDampersDataQuad4::_cells[] = {
  3, 2,
  6, 7,
};


const double pylith::bc::AbsorbingDampersDataQuad4::_dampingConsts[] = {
  1.25e+07, 7.5e+06,
  1.25e+07, 7.5e+06,
};
const double pylith::bc::AbsorbingDampersDataQuad4::_valsResidual[] = {
  -3.75000000e+06,    7.50000000e+05,
  -3.75000000e+06,    7.50000000e+05,
  0.0, 0.0, 
  0.0, 0.0,
  -1.37500000e+07,   -6.75000000e+06,
  -1.37500000e+07,   -6.75000000e+06,
};
const double pylith::bc::AbsorbingDampersDataQuad4::_valsJacobian[] = {
  1.25e+07, 0.0, // 0x
  1.25e+07, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 7.5e+06, // 0y
  0.0, 7.5e+06,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  1.25e+07, 0.0, // 1x
  1.25e+07, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 7.5e+06, // 1y
  0.0, 7.5e+06,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 2x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 2y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 3x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 3y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 4x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  1.25e+07, 0.0,
  1.25e+07, 0.0,
  0.0, 0.0, // 4y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 7.5e+06,
  0.0, 7.5e+06,
  0.0, 0.0, // 5x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  1.25e+07, 0.0,
  1.25e+07, 0.0,
  0.0, 0.0, // 5y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 7.5e+06,
  0.0, 7.5e+06,
};

pylith::bc::AbsorbingDampersDataQuad4::AbsorbingDampersDataQuad4(void)
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
  fieldTIncr = const_cast<double*>(_fieldTIncr);
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

pylith::bc::AbsorbingDampersDataQuad4::~AbsorbingDampersDataQuad4(void)
{}


// End of file
