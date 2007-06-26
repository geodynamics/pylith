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

/* Original mesh
 *
 * Cells are 0-3, vertices are 4-9.
 *
 *
 *         9
 *        / \
 *       /   \
 *      /     \
 *     /       \
 *    8---------5
 *     \       /|\
 *      \     / | \
 *       \   /  |  \
 *        \ /   |   \
 *         4    |    7
 *          \   |   /
 *           \  |  /
 *            \ | /
 *             \|/
 *              6
 *
 *
 * After adding cohesive elements
 *
 * Cells are 0-3, 13-14, vertices are 4-12.
 *
 *         9
 *        / \
 *       /   \
 *      /     \
 *     /       \
 *    8---------  5
 * 15 |        11/|
 *   14--------10 |
 *     \       /| |\
 *      \     / | | \
 *       \   /  | |  \
 *        \ /   | |   \
 *         4    | |    7
 *          \   | |   /
 *           \  | |  /
 *            \ | | /
 *             \| |/
 *             12-6
 *               13
 */


#include "CohesiveKinDataTri3d.hh"

const char* pylith::faults::CohesiveKinDataTri3d::_meshFilename =
  "data/tri3d.mesh";

const int pylith::faults::CohesiveKinDataTri3d::_spaceDim = 2;

const int pylith::faults::CohesiveKinDataTri3d::_cellDim = 1;

const int pylith::faults::CohesiveKinDataTri3d::_numBasis = 2;

const int pylith::faults::CohesiveKinDataTri3d::_numQuadPts = 1;

const double pylith::faults::CohesiveKinDataTri3d::_quadPts[] = {
  0.0,
};

const double pylith::faults::CohesiveKinDataTri3d::_quadWts[] = {
  2.0,
};

const double pylith::faults::CohesiveKinDataTri3d::_basis[] = {
  0.5,
  0.5
};

const double pylith::faults::CohesiveKinDataTri3d::_basisDeriv[] = {
  -0.5,
   0.5
};

const double pylith::faults::CohesiveKinDataTri3d::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveKinDataTri3d::_id = 10;

const char* pylith::faults::CohesiveKinDataTri3d::_label = "fault";

const char* pylith::faults::CohesiveKinDataTri3d::_finalSlipFilename = 
  "data/tri3d_finalslip.spatialdb";

const char* pylith::faults::CohesiveKinDataTri3d::_slipTimeFilename = 
  "data/tri3d_sliptime.spatialdb";

const char* pylith::faults::CohesiveKinDataTri3d::_peakRateFilename = 
  "data/tri3d_peakrate.spatialdb";

const char* pylith::faults::CohesiveKinDataTri3d::_matPropsFilename = 
  "data/bulkprops_2d.spatialdb";

const double pylith::faults::CohesiveKinDataTri3d::_fieldT[] = {
  6.1, 8.1,
  6.2, 8.2,
  6.3, 8.3,
  6.4, 8.4,
  6.5, 8.5,
  6.6, 8.6,
  6.7, 8.7,
  6.8, 8.8, // 11
  6.9, 8.9,
  6.0, 8.0, // 13
  7.1, 9.1,
  7.2, 9.2, // 15
};

const int pylith::faults::CohesiveKinDataTri3d::_numConstraintVert = 3;

const double pylith::faults::CohesiveKinDataTri3d::_orientation[] = {
  +0.70710678118654757, -0.70710678118654757,  
  -0.70710678118654757, -0.70710678118654757,
  0.0, -1.0,  -1.0,  0.0,
 +1.0,  0.0,   0.0, -1.0
};

const int pylith::faults::CohesiveKinDataTri3d::_constraintVertices[] = {
  11, 13, 15
};

const int pylith::faults::CohesiveKinDataTri3d::_constraintCells[] = {
  16, 16, 17
};

const double pylith::faults::CohesiveKinDataTri3d::_valsResidual[] = {
  0.0,  0.0,
 -1.4142135623730949, -11.030865786510143, // 5
 -8.0,  -6.0, // 6
  0.0,  0.0,
 +7.2,  -9.2, // 8
  0.0,  0.0,
 +1.4142135623730949, +11.030865786510143, // 10
  1.05057813143, 0.0456773100622, // 11
 +8.0, +6.0, // 12
  0.989535448086, 0.0824612873405, // 13
 -7.2, +9.2, // 14
  0.90435792846,  0.10852295130, // 15
};

const double pylith::faults::CohesiveKinDataTri3d::_valsJacobian[] = {
  0.0, 0.0, // 4x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 4y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 5x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -0.70710678118654757, +0.70710678118654757, // 11
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 5y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +0.70710678118654757, +0.70710678118654757, // 11
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 6x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 13
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 6y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 13
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 7x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 7y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 8x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, // ??
  0.0, 0.0, // 8y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // ??
  0.0, 0.0, // 9x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 9y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 10x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +0.70710678118654757, -0.70710678118654757, // 11
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 10y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -0.70710678118654757, -0.70710678118654757, // 11
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 11x
 -0.70710678118654757, +0.70710678118654757, // 5
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +0.70710678118654757, -0.70710678118654757, // 10
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 11y
 +0.70710678118654757, +0.70710678118654757, // 5
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -0.70710678118654757, -0.70710678118654757, // 10
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 12x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, // 13
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 12y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, // 13
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 13x
  0.0, 0.0,
  0.0,+1.0, // 6
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, // 12
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 13y
  0.0, 0.0,
 +1.0, 0.0, // 6
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, // 12
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 14x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 15
  0.0, 0.0, // 14y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, // 15
  0.0, 0.0, // 15x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, // 8
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 14
  0.0, 0.0,
  0.0, 0.0, // 15y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 8
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, // 14
  0.0, 0.0,
};

const double pylith::faults::CohesiveKinDataTri3d::_pseudoStiffness = 2.4;

pylith::faults::CohesiveKinDataTri3d::CohesiveKinDataTri3d(void)
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
  finalSlipFilename = const_cast<char*>(_finalSlipFilename);
  slipTimeFilename = const_cast<char*>(_slipTimeFilename);
  peakRateFilename = const_cast<char*>(_peakRateFilename);
  matPropsFilename = const_cast<char*>(_matPropsFilename);
  fieldT = const_cast<double*>(_fieldT);
  orientation = const_cast<double*>(_orientation);
  constraintVertices = const_cast<int*>(_constraintVertices);
  constraintCells = const_cast<int*>(_constraintCells);
  valsResidual = const_cast<double*>(_valsResidual);
  valsJacobian = const_cast<double*>(_valsJacobian);
  pseudoStiffness = _pseudoStiffness;
  numConstraintVert = _numConstraintVert;  
} // constructor

pylith::faults::CohesiveKinDataTri3d::~CohesiveKinDataTri3d(void)
{}


// End of file
