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

/* Original mesh using Sieve labels.
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
 * After adding cohesive elements using Sieve labels.
 *
 * Cells are 0-1, 8, vertices are 2-7.
 *
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

#include "CohesiveDynLDataTri3.hh"

const char* pylith::faults::CohesiveDynLDataTri3::_meshFilename =
  "data/tri3.mesh";

const int pylith::faults::CohesiveDynLDataTri3::_spaceDim = 2;

const int pylith::faults::CohesiveDynLDataTri3::_cellDim = 1;

const int pylith::faults::CohesiveDynLDataTri3::_numBasis = 2;

const int pylith::faults::CohesiveDynLDataTri3::_numQuadPts = 1;

const double pylith::faults::CohesiveDynLDataTri3::_quadPts[] = {
  0.0,
};

const double pylith::faults::CohesiveDynLDataTri3::_quadWts[] = {
  2.0,
};

const double pylith::faults::CohesiveDynLDataTri3::_basis[] = {
  0.5,
  0.5
};

const double pylith::faults::CohesiveDynLDataTri3::_basisDeriv[] = {
  -0.5,
   0.5
};

const double pylith::faults::CohesiveDynLDataTri3::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveDynLDataTri3::_id = 10;

const char* pylith::faults::CohesiveDynLDataTri3::_label = "fault";

const char* pylith::faults::CohesiveDynLDataTri3::_initialTractFilename = 
  "data/tri3_initialtract.spatialdb";

const double pylith::faults::CohesiveDynLDataTri3::_fieldT[] = {
  8.1, 9.1,
  8.2, 9.2, // 3
  8.3, 9.3, // 4
  8.4, 9.4,
  8.5, 9.5, // 6
  8.7, 9.7, // 7
  8.6, 9.6, // 8
  8.8, 9.8, // 9
};

// :TODO: Make sensible values for Jacobian for DOF on positive and
// negative sides of the fault. Add semi-random values for other DOF.
const double pylith::faults::CohesiveDynLDataTri3::_jacobian[] = {
  1.0, 0.0, // 2x (values for row associated with vertex with label 2, x DOF)
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 1.0, // 2y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 3x (values for row associated with vertex label 3, x DOF)
  1.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 8 (column for DOF associated with vertex label 8)
  0.0, 0.0,
  0.0, 0.0, // 3y
  0.0, 1.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, //  8
  0.0, 0.0,
  0.0, 0.0, // 4x
  0.0, 0.0,
  1.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, //  9
  0.0, 0.0, // 4y
  0.0, 0.0,
  0.0, 1.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, //  9
  0.0, 0.0, // 5x
  0.0, 0.0,
  0.0, 0.0,
  1.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 5y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 1.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 6x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  1.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, //  8
  0.0, 0.0,
  0.0, 0.0, // 6y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 1.0,
  0.0, 0.0,
 -1.0, 0.0, //  8
  0.0, 0.0,
  0.0, 0.0, // 7x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  1.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, //  9
  0.0, 0.0, // 7y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 1.0,
  0.0, 0.0,
 -1.0, 0.0, //  9

  0.0, 0.0, // 8x (rows associated with Lagrange multiplier vertex label 8)
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

// ----------------------------------------------------------------------
// Computed values
// ----------------------------------------------------------------------

const double pylith::faults::CohesiveDynLDataTri3::_orientation[] = {
  0.0, -1.0,  -1.0, 0.0,
  0.0, -1.0,  -1.0, 0.0
};

const double pylith::faults::CohesiveDynLDataTri3::_area[] = {
  1.0,
  1.0,
};

const double pylith::faults::CohesiveDynLDataTri3::_initialTractions[] = {
  1.0, -2.0,
  1.1, -2.1,
};


const int pylith::faults::CohesiveDynLDataTri3::_numConstraintVert = 2;
const int pylith::faults::CohesiveDynLDataTri3::_constraintVertices[] = {
  8, 9
};
const int pylith::faults::CohesiveDynLDataTri3::_constraintCells[] = {
  11, 11
};

// ----------------------------------------------------------------------
// Stick case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynLDataTri3::_fieldIncrStick[] = {
  8.1, 9.1,
  8.2, 9.2, // 3
  8.3, 9.3, // 4
  8.4, 9.4,
  8.5, 9.5, // 6
  8.7, 9.7, // 7
  8.6, 9.6, // 8
  8.8, 9.8, // 9
};

// No change in fieldIncr
// Zero slip

// ----------------------------------------------------------------------
// Slip case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynLDataTri3::_fieldIncrSlip[] = {
  8.1, 9.1,
  8.2, 9.2, // 3
  8.3, 9.3, // 4
  8.4, 9.4,
  8.5, 9.5, // 6
  8.7, 9.7, // 7
  8.6, 9.6, // 8
  8.8, 9.8, // 9
};

// Output
// :TODO: Update Lagrange multiplier values
const double pylith::faults::CohesiveDynLDataTri3::_fieldIncrSlipE[] = {
  8.1, 9.1,
  8.2, 9.2, // 3
  8.3, 9.3, // 4
  8.4, 9.4,
  8.5, 9.5, // 6
  8.7, 9.7, // 7
  8.6, 9.6, // 8
  8.8, 9.8, // 9
};

// :TODO: Update slip values based on changes in Lagrange multiplier values
const double pylith::faults::CohesiveDynLDataTri3::_slipSlipE[] = {
  0.0, 0.0,
  0.0, 0.0,
};

// ----------------------------------------------------------------------
// Open case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynLDataTri3::_fieldIncrOpen[] = {
  8.1, 9.1,
  8.2, 9.2, // 3
  8.3, 9.3, // 4
  8.4, 9.4,
  8.5, 9.5, // 6
  8.7, 9.7, // 7
  8.6, 9.6, // 8
  8.8, 9.8, // 9
};

// Output
const double pylith::faults::CohesiveDynLDataTri3::_fieldIncrOpenE[] = {
  8.1, 9.1,
  8.2, 9.2, // 3
  8.3, 9.3, // 4
  8.4, 9.4,
  8.5, 9.5, // 6
  8.7, 9.7, // 7
  8.6, 9.6, // 8
  8.8, 9.8, // 9
};

const double pylith::faults::CohesiveDynLDataTri3::_slipOpenE[] = {
  0.0, 0.0,
  0.0, 0.0,
};

// ----------------------------------------------------------------------
pylith::faults::CohesiveDynLDataTri3::CohesiveDynLDataTri3(void)
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
  initialTractFilename = const_cast<char*>(_initialTractFilename);

  fieldT = const_cast<double*>(_fieldT);
  jacobian = const_cast<double*>(_jacobian);
  orientation = const_cast<double*>(_orientation);
  area = const_cast<double*>(_area);
  initialTractions = const_cast<double*>(_initialTractions);

  constraintVertices = const_cast<int*>(_constraintVertices);
  constraintCells = const_cast<int*>(_constraintCells);
  numConstraintVert = _numConstraintVert;  

  // Stick
  fieldIncrStick = const_cast<double*>(_fieldIncrStick);

  // Slip
  fieldIncrSlip = const_cast<double*>(_fieldIncrSlip);
  fieldIncrSlipE = const_cast<double*>(_fieldIncrSlipE);
  slipSlipE = const_cast<double*>(_slipSlipE);

  // Open
  fieldIncrOpen = const_cast<double*>(_fieldIncrOpen);
  fieldIncrOpenE = const_cast<double*>(_fieldIncrOpenE);
  slipOpenE = const_cast<double*>(_slipOpenE);
} // constructor

pylith::faults::CohesiveDynLDataTri3::~CohesiveDynLDataTri3(void)
{}


// End of file
