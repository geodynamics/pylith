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

/* Mesh: hex8b.mesh
 *
 * Neumann BC on faces defined by vertices 0, 2, 4, 6, 8, 10.
 *
 * Constant value of 4 Pa horizontal shear applied on both faces.
 * Since each face has an area of 2 m^2, we should obtain the following
 * values:
 * 2 N in x-direction:  Vertices 0, 2, 8, 10
 * 4 N in x-direction:  Vertices 4, 6.
 */

#include "NeumannDataHex8.hh"

const char* pylith::bc::NeumannDataHex8::_meshFilename = 
  "data/hex8b.mesh";
const int pylith::bc::NeumannDataHex8::_spaceDim = 3;
const int pylith::bc::NeumannDataHex8::_cellDim = 2;
const int pylith::bc::NeumannDataHex8::_numBasis = 4;
const int pylith::bc::NeumannDataHex8::_numQuadPts = 4;
const double pylith::bc::NeumannDataHex8::_quadPts[] = {
  -0.57735027, -0.57735027,
  +0.57735027, -0.57735027,
  +0.57735027, +0.57735027,
  -0.57735027, +0.57735027,
};
const double pylith::bc::NeumannDataHex8::_quadWts[] = {
  1.0, 1.0, 1.0, 1.0
};
const double pylith::bc::NeumannDataHex8::_basis[] = {
  0.62200847,  0.16666667,  0.16666667,  0.0446582,
  0.16666667,  0.62200847,  0.0446582,   0.16666667,
  0.16666667,  0.0446582,   0.62200847,  0.16666667,
  0.0446582,   0.16666667,  0.16666667,  0.62200847,
};
const double pylith::bc::NeumannDataHex8::_basisDeriv[] = {
  -0.39433757, -0.39433757,
  +0.39433757, -0.10566243,
  +0.10566243, +0.10566243,
  -0.10566243, +0.39433757,

  -0.39433757, -0.10566243,
  +0.39433757, -0.39433757,
  +0.10566243, +0.39433757,
  -0.10566243, +0.10566243,

  -0.10566243, -0.10566243,
  +0.10566243, -0.39433757,
  +0.39433757, +0.39433757,
  -0.39433757, +0.10566243,

  -0.10566243, -0.39433757,
  +0.10566243, -0.10566243,
  +0.39433757, +0.10566243,
  -0.39433757, +0.39433757,
};
const double pylith::bc::NeumannDataHex8::_verticesRef[] = {
  -1.0, -1.0,
  +1.0, -1.0,
  +1.0, +1.0,
  -1.0, +1.0
};

const int pylith::bc::NeumannDataHex8::_id = 0;
const char* pylith::bc::NeumannDataHex8::_label = "tractionVerts";
const char* pylith::bc::NeumannDataHex8::_dbFilename =
  "data/hex8b.spatialdb";

const int pylith::bc::NeumannDataHex8::_numBoundaryCells = 2;
const int pylith::bc::NeumannDataHex8::_numVertices = 12;
const int pylith::bc::NeumannDataHex8::_numCorners[] = { 4, 4 };
const int pylith::bc::NeumannDataHex8::_cells[] = { 0,  4,  6,  2,
						    4,  8, 10,  6 };
const double pylith::bc::NeumannDataHex8::_tractionsCell[] = { 1.0, 0.0, 0.0,
							       1.0, 0.0, 0.0,
							       1.0, 0.0, 0.0,
							       1.0, 0.0, 0.0,
							       1.0, 0.0, 0.0,
							       1.0, 0.0, 0.0,
							       1.0, 0.0, 0.0,
							       1.0, 0.0, 0.0};
const double pylith::bc::NeumannDataHex8::_valsResidual[] = { 2.0, 0.0, 0.0,
							      0.0, 0.0, 0.0,
							      2.0, 0.0, 0.0,
							      0.0, 0.0, 0.0,
							      4.0, 0.0, 0.0,
							      0.0, 0.0, 0.0,
							      4.0, 0.0, 0.0,
							      0.0, 0.0, 0.0,
							      2.0, 0.0, 0.0,
							      0.0, 0.0, 0.0,
							      2.0, 0.0, 0.0,
							      0.0, 0.0, 0.0};


pylith::bc::NeumannDataHex8::NeumannDataHex8(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  quadPts = const_cast<double*>(_quadPts);
  quadWts = const_cast<double*>(_quadWts);
  basis = const_cast<double*>(_basis);
  basisDeriv = const_cast<double*>(_basisDeriv);
  verticesRef = const_cast<double*>(_verticesRef);
  id = _id;
  label = const_cast<char*>(_label);
  dbFilename = const_cast<char*>(_dbFilename);
  numBoundaryCells = _numBoundaryCells;
  numVertices = _numVertices;
  numCorners = const_cast<int*>(_numCorners);
  cells = const_cast<int*>(_cells);
  tractionsCell = const_cast<double*>(_tractionsCell);
  valsResidual = const_cast<double*>(_valsResidual);

} // constructor

pylith::bc::NeumannDataHex8::~NeumannDataHex8(void)
{}


// End of file
