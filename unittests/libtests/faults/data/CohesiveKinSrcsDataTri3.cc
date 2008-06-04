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
 * Cells are 0-1, vertices are 2-5.
 *
 *              3
 *             /|\
 *            / | \
 *           /  |  \
 *          /   |   \
 *         2    |    5
 *          \   |   /
 *           \  |  /
 *            \ | /
 *             \|/
 *              4
 *
 *
 * After adding cohesive elements
 *
 * Cells are 0-1, 8, vertices are 2-7.
 *
 *              6 -7- 3
 *             /|     |\
 *            / |     | \
 *           /  |     |  \
 *          /   |     |   \
 *         2    |     |    5
 *          \   |     |   /
 *           \  |     |  /
 *            \ |     | /
 *             \|     |/
 *              8 -9- 4
 *              6 -8- 3
 *             /|     |\
 *            / |     | \
 *           /  |     |  \
 *          /   |     |   \
 *         2    |     |    5
 *          \   |     |   /
 *           \  |     |  /
 *            \ |     | /
 *             \|     |/
 *              7 -9- 4
 */

#include "CohesiveKinSrcsDataTri3.hh"

const char* pylith::faults::CohesiveKinSrcsDataTri3::_meshFilename =
  "data/tri3.mesh";

const int pylith::faults::CohesiveKinSrcsDataTri3::_spaceDim = 2;

const int pylith::faults::CohesiveKinSrcsDataTri3::_cellDim = 1;

const int pylith::faults::CohesiveKinSrcsDataTri3::_numBasis = 2;

const int pylith::faults::CohesiveKinSrcsDataTri3::_numQuadPts = 1;

const double pylith::faults::CohesiveKinSrcsDataTri3::_quadPts[] = {
  0.0,
};

const double pylith::faults::CohesiveKinSrcsDataTri3::_quadWts[] = {
  2.0,
};

const double pylith::faults::CohesiveKinSrcsDataTri3::_basis[] = {
  0.5,
  0.5
};

const double pylith::faults::CohesiveKinSrcsDataTri3::_basisDeriv[] = {
  -0.5,
   0.5
};

const double pylith::faults::CohesiveKinSrcsDataTri3::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveKinSrcsDataTri3::_id = 10;

const char* pylith::faults::CohesiveKinSrcsDataTri3::_label = "fault";

const char* pylith::faults::CohesiveKinSrcsDataTri3::_finalSlipFilename = 
  "data/tri3_finalslip.spatialdb";

const char* pylith::faults::CohesiveKinSrcsDataTri3::_slipTimeFilename = 
  "data/tri3_sliptime.spatialdb";

const char* pylith::faults::CohesiveKinSrcsDataTri3::_peakRateFilename = 
  "data/tri3_peakrate.spatialdb";

const char* pylith::faults::CohesiveKinSrcsDataTri3::_matPropsFilename = 
  "data/bulkprops_2d.spatialdb";

const double pylith::faults::CohesiveKinSrcsDataTri3::_fieldT[] = {
  8.1, 9.1,
  8.2, 9.2,
  8.3, 9.3,
  8.4, 9.4,
  8.5, 9.5,
  8.7, 9.7,
  8.6, 9.6, // 8
  8.8, 9.8, // 9
};

const int pylith::faults::CohesiveKinSrcsDataTri3::_numConstraintVert = 2;

const double pylith::faults::CohesiveKinSrcsDataTri3::_orientation[] = {
  0.0, -1.0,  -1.0, 0.0,
  0.0, -1.0,  -1.0, 0.0
};

const double pylith::faults::CohesiveKinSrcsDataTri3::_area[] = {
  1.0,
  1.0,
};

const int pylith::faults::CohesiveKinSrcsDataTri3::_constraintVertices[] = {
  8, 9
};

const int pylith::faults::CohesiveKinSrcsDataTri3::_constraintCells[] = {
  11, 11
};

const double pylith::faults::CohesiveKinSrcsDataTri3::_valsResidual[] = {
  0.0,  0.0,
  0.0,  0.0,
  0.0,  0.0,
  0.0,  0.0,
  0.0,  0.0,
  0.0,  0.0,
  1.42290872258, 0.06186559663, // 8
  1.25514570482, 0.10459547540, // 9
};

const double pylith::faults::CohesiveKinSrcsDataTri3::_valsResidualIncr[] = {
  0.0,  0.0,
 -9.6, -8.6, // 3
 -9.8, -8.8, // 4
  0.0,  0.0,
 +9.6, +8.6, // 6
 +9.8, +8.8, // 7
  0.02583782954, 0.00112338389, // 8
  0.02698044341, 0.00224837028, // 9
};

const double pylith::faults::CohesiveKinSrcsDataTri3::_valsJacobian[] = {
  0.0, 0.0, // 2x
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 3x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 8
  0.0, 0.0,
  0.0, 0.0, // 3y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, //  8
  0.0, 0.0,
  0.0, 0.0, // 4x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, //  9
  0.0, 0.0, // 4y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, //  9
  0.0, 0.0, // 5x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0,
  0.0, 0.0, // 6x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, //  8
  0.0, 0.0,
  0.0, 0.0, // 6y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  8
  0.0, 0.0,
  0.0, 0.0, // 7x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, //  9
  0.0, 0.0, // 7y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  9

  0.0, 0.0, // 8x
  0.0,+1.0, //  3
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, //  6
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 8y
 +1.0, 0.0, //  3
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  6
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,

  0.0, 0.0, // 9x
  0.0, 0.0,
  0.0,+1.0, //  4
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, //  7
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 9y
  0.0, 0.0,
 +1.0, 0.0, //  4
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  7
  0.0, 0.0,
  0.0, 0.0,
};

const double pylith::faults::CohesiveKinSrcsDataTri3::_pseudoStiffness = 2.4;

pylith::faults::CohesiveKinSrcsDataTri3::CohesiveKinSrcsDataTri3(void)
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
  area = const_cast<double*>(_area);
  constraintVertices = const_cast<int*>(_constraintVertices);
  constraintCells = const_cast<int*>(_constraintCells);
  valsResidual = const_cast<double*>(_valsResidual);
  valsResidualIncr = const_cast<double*>(_valsResidualIncr);
  valsJacobian = const_cast<double*>(_valsJacobian);
  pseudoStiffness = _pseudoStiffness;
  numConstraintVert = _numConstraintVert;  
} // constructor

pylith::faults::CohesiveKinSrcsDataTri3::~CohesiveKinSrcsDataTri3(void)
{}


// End of file
