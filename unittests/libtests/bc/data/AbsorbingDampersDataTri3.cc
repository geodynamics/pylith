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

#include "AbsorbingDampersDataTri3.hh"

const char* pylith::bc::AbsorbingDampersDataTri3::_meshFilename = 
  "data/tri3.mesh";

const int pylith::bc::AbsorbingDampersDataTri3::_numBasis = 2;
const int pylith::bc::AbsorbingDampersDataTri3::_numQuadPts = 1;
const double pylith::bc::AbsorbingDampersDataTri3::_quadPts[] = {
  0.0,
};
const double pylith::bc::AbsorbingDampersDataTri3::_quadWts[] = {
  2.0,
};
const double pylith::bc::AbsorbingDampersDataTri3::_basis[] = {
  0.5,
  0.5,
};
const double pylith::bc::AbsorbingDampersDataTri3::_basisDerivRef[] = {
  -0.5,
   0.5,
};

const char* pylith::bc::AbsorbingDampersDataTri3::_spatialDBFilename = "data/elasticplanestrain.spatialdb";
const int pylith::bc::AbsorbingDampersDataTri3::_id = 2;
const char* pylith::bc::AbsorbingDampersDataTri3::_label = "bc";

const double pylith::bc::AbsorbingDampersDataTri3::_dt =   0.25;
const double pylith::bc::AbsorbingDampersDataTri3::_fieldTpdt[] = {
  1.0,
  1.1,
  1.2,
  1.3,
};
const double pylith::bc::AbsorbingDampersDataTri3::_fieldT[] = {
  1.1,
  1.3,
  1.5,
  1.7,
};
const double pylith::bc::AbsorbingDampersDataTri3::_fieldTmdt[] = {
  1.2,
  1.5,
  1.8,
  2.1,
};

const int pylith::bc::AbsorbingDampersDataTri3::_spaceDim = 2;
const int pylith::bc::AbsorbingDampersDataTri3::_cellDim = 1;
const int pylith::bc::AbsorbingDampersDataTri3::_numVertices = 2;
const int pylith::bc::AbsorbingDampersDataTri3::_numCells = 1;
const int pylith::bc::AbsorbingDampersDataTri3::_numCorners = 2;
const int pylith::bc::AbsorbingDampersDataTri3::_cells[] = {
  3, 5,
};


const double pylith::bc::AbsorbingDampersDataTri3::_dampingConsts[] = {
  12.5e+6, 7.5e+6,
};
const double pylith::bc::AbsorbingDampersDataTri3::_valsResidual[] = {
  0.0, 0.0,
  2.4e+07, 6.0e+06,
  0.0, 0.0,
  2.4e+07, 6.0e+06
};
const double pylith::bc::AbsorbingDampersDataTri3::_valsJacobian[] = {
  0.0, 0.0, // 0   // FIX THESE
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 1
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 2
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 3
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
};

pylith::bc::AbsorbingDampersDataTri3::AbsorbingDampersDataTri3(void)
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

pylith::bc::AbsorbingDampersDataTri3::~AbsorbingDampersDataTri3(void)
{}


// End of file
