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
 * Cells are 0-1, vertices are 2-6.
 *
 * 2   3,4,5  6
 *
 *     ^^^^^ Face in x-y plane
 *
 * After adding cohesive elements
 *
 * Cells are 0-1,10, vertices are 2-9.
 *
 * 2   3,4,5  7,8,9   6
 *             10,11,12
 *     ^^^^^^^^^^^^ Cohesive element in x-y plane.
 */

#include "CohesiveDynDataTet4.hh"

const char* pylith::faults::CohesiveDynDataTet4::_meshFilename =
  "data/tet4.mesh";

const int pylith::faults::CohesiveDynDataTet4::_spaceDim = 3;

const int pylith::faults::CohesiveDynDataTet4::_cellDim = 2;

const int pylith::faults::CohesiveDynDataTet4::_numBasis = 3;

const int pylith::faults::CohesiveDynDataTet4::_numQuadPts = 1;

const double pylith::faults::CohesiveDynDataTet4::_quadPts[] = {
  -3.33333333e-01,  -3.33333333e-01,
};

const double pylith::faults::CohesiveDynDataTet4::_quadWts[] = {
  2.0,
};

const double pylith::faults::CohesiveDynDataTet4::_basis[] = {
  3.33333333e-01,  3.33333333e-01,
  3.33333333e-01,};

const double pylith::faults::CohesiveDynDataTet4::_basisDeriv[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.00000000e+00,
};

const double pylith::faults::CohesiveDynDataTet4::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
};

const int pylith::faults::CohesiveDynDataTet4::_id = 10;

const char* pylith::faults::CohesiveDynDataTet4::_label = "fault";

const char* pylith::faults::CohesiveDynDataTet4::_initialTractFilename = 
  "data/tet4_initialtract.spatialdb";

const double pylith::faults::CohesiveDynDataTet4::_fieldT[] = {
  7.1, 8.1, 9.1,
  7.2, 8.2, 9.2, // 3
  7.3, 8.3, 9.3, // 4
  7.4, 8.4, 9.4, // 5
  7.5, 8.5, 9.5,
  7.6, 8.6, 9.6, // 7
  7.8, 8.8, 9.8, // 8
  7.0, 8.0, 9.0, // 9
  7.7, 8.7, 9.7, // 10
  7.9, 8.9, 9.9, // 11
  7.1, 8.1, 9.1, // 12
};

// :TODO: Make sensible values for Jacobian for DOF on positive and
// negative sides of the fault. Add semi-random values for other DOF.
const double pylith::faults::CohesiveDynDataTet4::_jacobian[] = {
    1,  0.1,  0.2,
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,    1,  1.1,
  1.2,  1.3,  1.4,
  1.5,  1.6,  1.7,
  1.8,  1.9,    2,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
    3,  3.1,  3.2,
  3.3,    1,  3.4,
  3.5,  3.6,  3.7,
  3.8,  3.9,    4,
  4.1,  4.2,  4.3,
  4.4,  4.5,  4.6,
  4.7,  4.8,  4.9,
    5,  5.1,  5.2,
  5.3,  5.4,  5.5,
  5.6,  5.7,  5.8,
  5.9,    6,  6.1,
  6.2,  6.3,  6.4,
  6.5,  6.6,    1,
  6.7,  6.8,  6.9,
    7,  7.1,  7.2,
  7.3,  7.4,  7.5,
  7.6,  7.7,  7.8,
  7.9,    8,  8.1,
  8.2,  8.3,  8.4,
  8.5,  8.6,  8.7,
  8.8,  8.9,    9,
  9.1,  9.2,  9.3,
  9.4,  9.5,  9.6,
  9.7,  9.8,  9.9,
    1,   10, 10.1,
 10.2, 10.3, 10.4,
 10.5, 10.6, 10.7,
 10.8, 10.9,   11,
 11.1, 11.2, 11.3,
 11.4, 11.5, 11.6,
 11.7, 11.8, 11.9,
   12, 12.1,   -1,
 12.2, 12.3, 12.4,
 12.5, 12.6, 12.7,
 12.8, 12.9,   13,
 13.1,    1, 13.2,
 13.3, 13.4, 13.5,
 13.6, 13.7, 13.8,
 13.9,   14, 14.1,
 14.2, 14.3, 14.4,
 14.5, 14.6, 14.7,
 14.8, 14.9,   15,
   -1, 15.1, 15.2,
 15.3, 15.4, 15.5,
 15.6, 15.7, 15.8,
 15.9,   16, 16.1,
 16.2, 16.3,    1,
 16.4, 16.5, 16.6,
 16.7, 16.8, 16.9,
   17, 17.1, 17.2,
 17.3, 17.4, 17.5,
 17.6, 17.7, 17.8,
 17.9,   18, 18.1,
 18.2,   -1, 18.3,
 18.4, 18.5, 18.6,
 18.7, 18.8, 18.9,
   19, 19.1, 19.2,
 19.3, 19.4, 19.5,
    1, 19.6, 19.7,
 19.8, 19.9,   20,
 20.1, 20.2, 20.3,
 20.4, 20.5, 20.6,
 20.7, 20.8, 20.9,
   21, 21.1, 21.2,
 21.3, 21.4, 21.5,
 21.6, 21.7,   -1,
 21.8, 21.9,   22,
 22.1, 22.2, 22.3,
 22.4, 22.5, 22.6,
 22.7,    1, 22.8,
 22.9,   23, 23.1,
 23.2, 23.3, 23.4,
 23.5, 23.6, 23.7,
 23.8, 23.9,   24,
 24.1, 24.2, 24.3,
 24.4, 24.5, 24.6,
   -1, 24.7, 24.8,
 24.9,   25, 25.1,
 25.2, 25.3, 25.4,
 25.5, 25.6, 25.7,
 25.8, 25.9,    1,
   26, 26.1, 26.2,
 26.3, 26.4, 26.5,
 26.6, 26.7, 26.8,
 26.9,   27, 27.1,
 27.2, 27.3, 27.4,
 27.5, 27.6, 27.7,
 27.8,   -1, 27.9,
   28, 28.1, 28.2,
 28.3, 28.4, 28.5,
 28.6, 28.7, 28.8,
 28.9,   29, 29.1,
    1, 29.2, 29.3,
 29.4, 29.5, 29.6,
 29.7, 29.8, 29.9,
   30, 30.1, 30.2,
 30.3, 30.4, 30.5,
 30.6, 30.7, 30.8,
 30.9,   31, 31.1,
 31.2, 31.3,   -1,
 31.4, 31.5, 31.6,
 31.7, 31.8, 31.9,
   32, 32.1, 32.2,
 32.3,    1, 32.4,
 32.5, 32.6, 32.7,
 32.8, 32.9,   33,
 33.1, 33.2, 33.3,
 33.4, 33.5, 33.6,
 33.7, 33.8, 33.9,
   34, 34.1, 34.2,
   -1, 34.3, 34.4,
 34.5, 34.6, 34.7,
 34.8, 34.9,   35,
 35.1, 35.2, 35.3,
 35.4, 35.5,    1,
 35.6, 35.7, 35.8,
 35.9,   36, 36.1,
 36.2, 36.3, 36.4,
 36.5, 36.6, 36.7,
 36.8, 36.9,   37,
 37.1, 37.2, 37.3,
 37.4,   -1, 37.5,
 37.6, 37.7, 37.8,
 37.9,   38, 38.1,
 38.2, 38.3, 38.4,
 38.5, 38.6, 38.7,
    1, 38.8, 38.9,
   39, 39.1, 39.2,
 39.3, 39.4, 39.5,
 39.6, 39.7, 39.8,
 39.9,   40, 40.1,
 40.2, 40.3, 40.4,
 40.5, 40.6, 40.7,
 40.8, 40.9,   41,
 41.1, 41.2, 41.3,
 41.4, 41.5, 41.6,
 41.7, 41.8, 41.9,
   42,    1, 42.1,
 42.2, 42.3, 42.4,
 42.5, 42.6, 42.7,
 42.8, 42.9,   43,
 43.1, 43.2, 43.3,
 43.4, 43.5, 43.6,
 43.7, 43.8, 43.9,
   44, 44.1, 44.2,
 44.3, 44.4, 44.5,
 44.6, 44.7, 44.8,
 44.9,   45, 45.1,
 45.2, 45.3,    1,
 45.4, 45.5, 45.6,
 45.7, 45.8, 45.9,
   46, 46.1, 46.2,
 46.3, 46.4, 46.5,
 46.6, 46.7, 46.8,
 46.9,   47, 47.1,
 47.2, 47.3, 47.4,
 47.5, 47.6, 47.7,
 47.8, 47.9,   48,
 48.1, 48.2, 48.3,
 48.4, 48.5, 48.6,
    1, 48.7, 48.8,
 48.9,   49, 49.1,
 49.2, 49.3, 49.4,
 49.5, 49.6,    1,
 49.7, 49.8, 49.9,
   50, 50.1, 50.2,
 50.3, 50.4, 50.5,
 50.6, 50.7, 50.8,
 50.9,   51, 51.1,
 51.2, 51.3, 51.4,
 51.5, 51.6, 51.7,
 51.8,    1, 51.9,
   52, 52.1, 52.2,
 52.3, 52.4, 52.5,
    1, 52.6, 52.7,
 52.8, 52.9,   53,
 53.1, 53.2, 53.3,
 53.4, 53.5, 53.6,
 53.7, 53.8, 53.9,
   54, 54.1, 54.2,
 54.3, 54.4, 54.5,
 54.6, 54.7, 54.8,
 54.9,   55,    1,
 55.1, 55.2, 55.3,
 55.4, 55.5, 55.6,
 55.7,    1, 55.8,
 55.9,   56, 56.1,
 56.2, 56.3, 56.4,
 56.5, 56.6, 56.7,
 56.8, 56.9,   57,
 57.1, 57.2, 57.3,
 57.4, 57.5, 57.6,
 57.7, 57.8, 57.9,
   58, 58.1, 58.2,
    1, 58.3, 58.4,
 58.5, 58.6, 58.7,
 58.8, 58.9,   59,
 59.1, 59.2,    1,
 59.3, 59.4, 59.5,
 59.6, 59.7, 59.8,
 59.9,   60, 60.1,
 60.2, 60.3, 60.4,
 60.5, 60.6, 60.7,
 60.8, 60.9,   61,
 61.1, 61.2, 61.3,
 61.4,    1, 61.5,
 61.6, 61.7, 61.8,
 61.9,   62, 62.1,
    1, 62.2, 62.3,
 62.4, 62.5, 62.6,
 62.7, 62.8, 62.9,
   63, 63.1, 63.2,
 63.3, 63.4, 63.5,
 63.6, 63.7, 63.8,
 63.9,   64, 64.1,
 64.2, 64.3, 64.4,
 64.5, 64.6,    1,
 64.7, 64.8, 64.9,
   65, 65.1, 65.2,
 65.3,    1, 65.4,
 65.5, 65.6, 65.7,
 65.8, 65.9,   66,
 66.1, 66.2, 66.3,
 66.4, 66.5, 66.6,
 66.7, 66.8, 66.9,
   67, 67.1, 67.2,
 67.3, 67.4, 67.5,
 67.6, 67.7, 67.8,
    1, 67.9,   68,
 68.1, 68.2, 68.3,
 68.4, 68.5, 68.6,
 68.7, 68.8,    1,
 68.9,   69, 69.1,
 69.2, 69.3, 69.4,
 69.5, 69.6, 69.7,
 69.8, 69.9,   70,
 70.1, 70.2, 70.3,
 70.4, 70.5, 70.6,
 70.7, 70.8, 70.9,
   71,    1, 71.1,
 71.2, 71.3, 71.4,
 71.5, 71.6, 71.7,
    1, 71.8, 71.9,
   72, 72.1, 72.2,
 72.3, 72.4, 72.5,
 72.6, 72.7, 72.8,
 72.9,   73, 73.1,
 73.2, 73.3, 73.4,
 73.5, 73.6, 73.7,
 73.8, 73.9,   74,
 74.1, 74.2,    1,
 74.3, 74.4, 74.5,
 74.6, 74.7, 74.8,
 74.9,    1,   75,
 75.1, 75.2, 75.3,
 75.4,   -1, 75.5,
 75.6, 75.7, 75.8,
 75.9,   76, 76.1,
 76.2, 76.3, 76.4,
 76.5,    1, 76.6,
 76.7, 76.8, 76.9,
   77, 77.1, 77.2,
 77.3, 77.4, 77.5,
 77.6, 77.7, 77.8,
 77.9,   78, 78.1,
 78.2, 78.3, 78.4,
 78.5, 78.6,   -1,
 78.7, 78.8, 78.9,
   79, 79.1, 79.2,
 79.3, 79.4, 79.5,
 79.6, 79.7,    1,
 79.8, 79.9,   80,
 80.1, 80.2, 80.3,
 80.4, 80.5, 80.6,
 80.7, 80.8, 80.9,
   81, 81.1, 81.2,
 81.3, 81.4, 81.5,
   -1, 81.6, 81.7,
 81.8, 81.9,   82,
 82.1, 82.2, 82.3,
 82.4, 82.5, 82.6,
    1, 82.7, 82.8,
 82.9,   83, 83.1,
 83.2, 83.3, 83.4,
 83.5, 83.6, 83.7,
 83.8, 83.9,   84,
 84.1, 84.2, 84.3,
 84.4, 84.5, 84.6,
 84.7, 84.8, 84.9,
   85,   -1, 85.1,
 85.2, 85.3, 85.4,
 85.5, 85.6, 85.7,
 85.8, 85.9,   86,
 86.1,    1, 86.2,
 86.3, 86.4, 86.5,
 86.6, 86.7, 86.8,
 86.9,   87, 87.1,
 87.2, 87.3, 87.4,
 87.5, 87.6, 87.7,
 87.8, 87.9,   88,
 88.1, 88.2,   -1,
 88.3, 88.4, 88.5,
 88.6, 88.7, 88.8,
 88.9,   89, 89.1,
 89.2, 89.3,    1,
 89.4, 89.5, 89.6,
 89.7, 89.8, 89.9,
   90, 90.1, 90.2,
 90.3, 90.4, 90.5,
 90.6, 90.7, 90.8,
 90.9,   91, 91.1,
   -1, 91.2, 91.3,
 91.4, 91.5, 91.6,
 91.7, 91.8, 91.9,
   92, 92.1, 92.2,
    1, 92.3, 92.4,
 92.5, 92.6, 92.7,
 92.8, 92.9,   93,
 93.1, 93.2, 93.3,
 93.4, 93.5, 93.6,
 93.7, 93.8, 93.9,
   94, 94.1, 94.2,
 94.3, 94.4, 94.5,
 94.6,   -1, 94.7,
 94.8, 94.9,   95,
 95.1, 95.2, 95.3,
 95.4, 95.5, 95.6,
 95.7,    1, 95.8,
 95.9,   96, 96.1,
 96.2, 96.3, 96.4,
 96.5, 96.6, 96.7,
 96.8, 96.9,   97,
 97.1, 97.2, 97.3,
 97.4, 97.5, 97.6,
 97.7, 97.8,   -1,
 97.9,   98, 98.1,
 98.2, 98.3, 98.4,
 98.5, 98.6, 98.7,
 98.8, 98.9,    1,
   99, 99.1, 99.2,
 99.3, 99.4, 99.5,
 99.6, 99.7, 99.8,
 99.9,  100,100.1,
100.2,100.3,100.4,
100.5,100.6,100.7,
   -1,100.8,100.9,
  101,101.1,101.2,
101.3,101.4,101.5,
101.6,101.7,101.8,
    1,101.9,  102,
102.1,102.2,102.3,
102.4,102.5,102.6,
102.7,102.8,102.9,
};

// ----------------------------------------------------------------------
// Computed values
// ----------------------------------------------------------------------

const double pylith::faults::CohesiveDynDataTet4::_orientation[] = {
  0.0, -1.0, 0.0,    0.0, 0.0, +1.0,    -1.0, 0.0, 0.0,
  0.0, -1.0, 0.0,    0.0, 0.0, +1.0,    -1.0, 0.0, 0.0,
  0.0, -1.0, 0.0,    0.0, 0.0, +1.0,    -1.0, 0.0, 0.0,
};

const double pylith::faults::CohesiveDynDataTet4::_area[] = {
  1.0/3.0, 
  1.0/3.0, 
  1.0/3.0,
};

const double pylith::faults::CohesiveDynDataTet4::_forcesInitial[] = {
  3.1/3.0, -1.1/3.0, +2.1/3.0,
  3.1/3.0, -1.1/3.0, +2.1/3.0,
  3.1/3.0, -1.1/3.0, +2.1/3.0,
};


const int pylith::faults::CohesiveDynDataTet4::_numConstraintVert = 3;
const int pylith::faults::CohesiveDynDataTet4::_constraintVertices[] = {
  10, 11, 12
};

// ----------------------------------------------------------------------
// Stick case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynDataTet4::_fieldIncrStick[] = {
  1.1, 2.1, 35.1,
  1.2, 2.2, 35.2, // 3
  1.3, 2.3, 35.3, // 4
  1.4, 2.4, 35.4, // 5
  1.5, 2.5, 35.5,
  1.6, 2.6, 35.6, // 7
  1.8, 2.8, 35.8, // 8
  1.0, 2.0, 35.0, // 9
  1.7, 2.7, -35.7, // 10
  1.9, 2.9, -35.9, // 11
  1.1, 2.1, -35.1, // 12
};

// No change in fieldIncr
// Zero slip

// ----------------------------------------------------------------------
// Slip case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynDataTet4::_fieldIncrSlip[] = {
  8.1, 9.1, 10.1,
  8.2, 9.2, 10.2, // 3
  8.3, 9.3, 10.3, // 4
  8.4, 9.4, 10.4, // 5
  8.5, 9.5, 10.5,
  8.6, 9.6, 10.6, // 7
  8.8, 9.8, 10.8, // 8
  8.0, 9.0, 10.0, // 9
  8.7, 9.7, -10.7, // 10
  8.9, 9.9, -10.9, // 11
  8.1, 9.1, -10.1, // 12
};

// Output
const double pylith::faults::CohesiveDynDataTet4::_fieldIncrSlipE[] = {
   8.100000000000,   9.100000000000,  10.100000000000,
   8.200000000000,   8.391727731714,  10.956284985259,
   8.300000000000,   8.791277340217,  10.815192651071,
   8.400000000000,   9.060755249362,  10.758883268626,
   8.500000000000,   9.500000000000,  10.500000000000,
   8.600000000000,  10.408272268286,   9.843715014741,
   8.800000000000,  10.308722659783,  10.284807348929,
   8.000000000000,   9.339244750638,   9.641116731374,
  -7.300777685147,  -8.252092036995, -10.700000000000,
  -7.500201409882,  -8.452606339630, -10.900000000000,
  -6.702681322117,  -7.650402548711, -10.100000000000,
};

const double pylith::faults::CohesiveDynDataTet4::_slipSlipE[] = {
  -1.616544536572,  -1.512569970519,   0.000000000000,
  -1.017445319566,  -1.030385302142,   0.000000000000,
  -0.678489501275,  -0.717766537252,   0.000000000000,
};

// ----------------------------------------------------------------------
// Open case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynDataTet4::_fieldIncrOpen[] = {
  8.1, 9.1, 10.1,
  8.2, 9.2, 10.2, // 3
  8.3, 9.3, 10.3, // 4
  8.4, 9.4, 10.4, // 5
  8.5, 9.5, 10.5,
  8.6, 9.6, 10.6, // 7
  8.8, 9.8, 10.8, // 8
  8.0, 9.0, 10.0, // 9
  8.7, 9.7, 10.7, // 10
  8.9, 9.9, 10.9, // 11
  8.1, 9.1, 10.1, // 12
};

// Output
const double pylith::faults::CohesiveDynDataTet4::_fieldIncrOpenE[] = {
   8.100000000000,   9.100000000000,  10.100000000000,
   8.200000000000,   8.707485448105,  11.300792216319,
   8.300000000000,   9.089936432069,  11.135493258511,
   8.400000000000,   9.353421152116,  11.068689140292,
   8.500000000000,   9.500000000000,  10.500000000000,
   8.600000000000,  10.092514551895,   9.499207783681,
   8.800000000000,  10.010063567931,   9.964506741489,
   8.000000000000,   9.046578847884,   9.331310859708,
  -7.700000000000,  -8.700000000000,  -9.700000000000,
  -7.900000000000,  -8.900000000000,  -9.900000000000,
  -7.100000000000,  -8.100000000000,  -9.100000000000,
};

const double pylith::faults::CohesiveDynDataTet4::_slipOpenE[] = {
  -0.985029103789,  -2.201584432637,   0.000000000000,
  -0.420127135861,  -1.670986517021,   0.000000000000,
  -0.093157695768,  -1.337378280584,   0.000000000000,
};

// ----------------------------------------------------------------------
pylith::faults::CohesiveDynDataTet4::CohesiveDynDataTet4(void)
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
  forcesInitial = const_cast<double*>(_forcesInitial);

  constraintVertices = const_cast<int*>(_constraintVertices);
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

pylith::faults::CohesiveDynDataTet4::~CohesiveDynDataTet4(void)
{}


// End of file
