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
 * Cells are 0-1, vertices are 2-7.
 *
 *       3 -------- 5 -------- 7
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       2 -------- 4 -------- 6
 *
 * After adding cohesive elements
 *
 * Cells are 0-1,10 vertices are 2-9.
 *
 *       3 -------- 5 -11-- 9 -------- 7
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       2 -------- 4 -10-- 8 -------- 6
 */

#include "CohesiveDynLDataQuad4.hh"

const char* pylith::faults::CohesiveDynLDataQuad4::_meshFilename =
  "data/quad4.mesh";

const int pylith::faults::CohesiveDynLDataQuad4::_spaceDim = 2;

const int pylith::faults::CohesiveDynLDataQuad4::_cellDim = 1;

const int pylith::faults::CohesiveDynLDataQuad4::_numBasis = 2;

const int pylith::faults::CohesiveDynLDataQuad4::_numQuadPts = 1;

const double pylith::faults::CohesiveDynLDataQuad4::_quadPts[] = {
  0.0,
};

const double pylith::faults::CohesiveDynLDataQuad4::_quadWts[] = {
  2.0,
};

const double pylith::faults::CohesiveDynLDataQuad4::_basis[] = {
  0.5,
  0.5
};

const double pylith::faults::CohesiveDynLDataQuad4::_basisDeriv[] = {
  -0.5,
   0.5
};

const double pylith::faults::CohesiveDynLDataQuad4::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveDynLDataQuad4::_id = 10;

const char* pylith::faults::CohesiveDynLDataQuad4::_label = "fault";

const char* pylith::faults::CohesiveDynLDataQuad4::_initialTractFilename = 
  "data/quad4_initialtract.spatialdb";

const double pylith::faults::CohesiveDynLDataQuad4::_fieldT[] = {
  8.1, 9.1,
  8.2, 9.2,
  8.3, 9.3, // 4
  8.4, 9.4, // 5
  8.5, 9.5,
  8.6, 9.6,
  8.7, 9.7, // 8
  8.9, 9.9, // 9
  8.8, 9.8, // 10
  8.0, 9.0, // 11
};

// :TODO: Make sensible values for Jacobian for DOF on positive and
// negative sides of the fault. Add semi-random values for other DOF.
const double pylith::faults::CohesiveDynLDataQuad4::_jacobian[] = {
  1.0, 0.0, // 2x
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 3x
  1.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 3y
  0.0, 1.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 4x
  0.0, 0.0,
  1.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, // 10
  0.0, 0.0,
  0.0, 0.0, // 4y
  0.0, 0.0,
  0.0, 1.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  10
  0.0, 0.0,
  0.0, 0.0, // 5x
  0.0, 0.0,
  0.0, 0.0,
  1.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, //  11
  0.0, 0.0, // 5y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 1.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  11
  0.0, 0.0, // 6x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  1.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 6y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 1.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 7x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  1.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 7y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 1.0,
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
  1.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, //  10
  0.0, 0.0,
  0.0, 0.0, // 8y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 1.0,
  0.0, 0.0,
 +1.0, 0.0, // 10
  0.0, 0.0,
  0.0, 0.0, // 9x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  1.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 11
  0.0, 0.0, // 9y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 1.0,
  0.0, 0.0,
 +1.0, 0.0, // 11
  0.0, 0.0, // 10x
  0.0, 0.0,
  0.0,-1.0, //  4
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, //  8
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 10y
  0.0, 0.0,
 -1.0, 0.0, //  4
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, //  8
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 11x
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, //  5
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, //  9
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 11y
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  5
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 9
  0.0, 0.0,
  0.0, 0.0,
};

// ----------------------------------------------------------------------
// Computed values
// ----------------------------------------------------------------------

const double pylith::faults::CohesiveDynLDataQuad4::_orientation[] = {
  0.0,  -1.0,  -1.0, 0.0,
  0.0,  -1.0,  -1.0, 0.0
};

const double pylith::faults::CohesiveDynLDataQuad4::_area[] = {
  1.0,
  1.0,
};

const double pylith::faults::CohesiveDynLDataQuad4::_initialTractions[] = {
  1.0, -2.0,
  1.1, -2.1,
};


const int pylith::faults::CohesiveDynLDataQuad4::_numConstraintVert = 2;
const int pylith::faults::CohesiveDynLDataQuad4::_constraintVertices[] = {
  10, 11
};
const int pylith::faults::CohesiveDynLDataQuad4::_constraintCells[] = {
  13, 13
};

// ----------------------------------------------------------------------
// Stick case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynLDataQuad4::_fieldIncrStick[] = {
  1.1, 29.1,
  1.2, 29.2,
  1.3, 29.3, // 4
  1.4, 29.4, // 5
  1.5, 29.5,
  1.6, 29.6,
  1.7, 29.7, // 8
  1.9, 29.9, // 9
  1.8, -29.8, // 10
  1.0, -29.0, // 11
};

// No change in fieldIncr
// Zero slip

// ----------------------------------------------------------------------
// Slip case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynLDataQuad4::_fieldIncrSlip[] = {
  9.1, 10.1,
  9.2, 10.2,
  9.3, 10.3, // 4
  9.4, 10.4, // 5
  9.5, 10.5,
  9.6, 10.6,
  9.7, 10.7, // 8
  9.9, 10.9, // 9
  9.8, -10.8, // 10
  9.0, -10.0, // 11
};

// Output
// :TODO: Update Lagrange multiplier values
const double pylith::faults::CohesiveDynLDataQuad4::_fieldIncrSlipE[] = {
  9.1, 10.1,
  9.2, 10.2,
  9.3, 10.3, // 4
  9.4, 10.4, // 5
  9.5, 10.5,
  9.6, 10.6,
  9.7, 10.7, // 8
  9.9, 10.9, // 9
  -7.0, -10.8, // 10
  -6.14, -10.0, // 11
};

// :TODO: Update slip values based on changes in Lagrange multiplier values
const double pylith::faults::CohesiveDynLDataQuad4::_slipSlipE[] = {
  33.6, 0.0,
  30.28, 0.0,
};

// ----------------------------------------------------------------------
// Open case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynLDataQuad4::_fieldIncrOpen[] = {
  9.1, 10.1,
  9.2, 10.2,
  9.3, 10.3, // 4
  9.4, 10.4, // 5
  9.5, 10.5,
  9.6, 10.6,
  9.7, 10.7, // 8
  9.9, 10.9, // 9
  9.8, 10.8, // 10
  9.0, 10.0, // 11
};

// Output
const double pylith::faults::CohesiveDynLDataQuad4::_fieldIncrOpenE[] = {
  9.1, 10.1,
  9.2, 10.2,
  9.3, 10.3, // 4
  9.4, 10.4, // 5
  9.5, 10.5,
  9.6, 10.6,
  9.7, 10.7, // 8
  9.9, 10.9, // 9
  -8.8, -9.8, // 10
  -8.0, -9.0, // 11
};

const double pylith::faults::CohesiveDynLDataQuad4::_slipOpenE[] = {
  37.2, 41.2,
  34.0, 38.0,
};

// ----------------------------------------------------------------------
pylith::faults::CohesiveDynLDataQuad4::CohesiveDynLDataQuad4(void)
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

pylith::faults::CohesiveDynLDataQuad4::~CohesiveDynLDataQuad4(void)
{}


// End of file
