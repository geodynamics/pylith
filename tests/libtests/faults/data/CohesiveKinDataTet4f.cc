// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
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
 * 2   7,9,11   3,4,5  6
 *      34,35,36
 *     ^^^^^^^^^^^^ Cohesive element in x-y plane.
 */

#include "CohesiveKinDataTet4f.hh"

const char* pylith::faults::CohesiveKinDataTet4f::_meshFilename =
  "data/tet4f.mesh";

const int pylith::faults::CohesiveKinDataTet4f::_spaceDim = 3;

const int pylith::faults::CohesiveKinDataTet4f::_cellDim = 2;

const int pylith::faults::CohesiveKinDataTet4f::_numBasis = 3;

const int pylith::faults::CohesiveKinDataTet4f::_numQuadPts = 3;

const PylithScalar pylith::faults::CohesiveKinDataTet4f::_quadPts[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
};

const PylithScalar pylith::faults::CohesiveKinDataTet4f::_quadWts[] = {
  2.0/3.0, 2.0/3.0, 2.0/3.0,
};

const PylithScalar pylith::faults::CohesiveKinDataTet4f::_basis[] = {
  1.0, 0.0, 0.0,
  0.0, 1.0, 0.0,
  0.0, 0.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveKinDataTet4f::_basisDeriv[] = {
 -0.50000000e+00, -0.50000000e+00,
  0.50000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.50000000e+00,
 -0.50000000e+00, -0.50000000e+00,
  0.50000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.50000000e+00,
 -0.50000000e+00, -0.50000000e+00,
  0.50000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.50000000e+00,
};

const PylithScalar pylith::faults::CohesiveKinDataTet4f::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
};

const int pylith::faults::CohesiveKinDataTet4f::_id = 10;

const char* pylith::faults::CohesiveKinDataTet4f::_label = "fault";

const char* pylith::faults::CohesiveKinDataTet4f::_finalSlipFilename = 
  "data/tet4_finalslip.spatialdb";

const char* pylith::faults::CohesiveKinDataTet4f::_slipTimeFilename = 
  "data/tet4_sliptime.spatialdb";

const char* pylith::faults::CohesiveKinDataTet4f::_riseTimeFilename = 
  "data/tet4_risetime.spatialdb";

const PylithScalar pylith::faults::CohesiveKinDataTet4f::_fieldT[] = {
  7.1, 8.1, 9.1,
  7.2, 8.2, 9.2, // 3
  7.3, 8.3, 9.3, // 4
  7.4, 8.4, 9.4, // 5
  7.5, 8.5, 9.5,
  7.6, 8.6, 9.6, // 7
  7.8, 8.8, 9.8, // 8
  7.0, 8.0, 9.0, // 9
  7.7, 8.7, 9.7, // 34
  7.9, 8.9, 9.9, // 35
  7.1, 8.1, 9.1, // 36
};


const PylithScalar pylith::faults::CohesiveKinDataTet4f::_fieldIncr[] = {
  3.1, 4.1, 5.1,
  3.2, 4.2, 5.2, // 3
  3.3, 4.3, 5.3, // 4
  3.4, 4.4, 5.4, // 5
  3.5, 4.5, 5.5,
  3.6, 4.6, 5.6, // 7
  3.8, 4.8, 5.8, // 8
  3.0, 4.0, 5.0, // 9
  3.7, 4.7, 5.7, // 34
  3.9, 4.9, 5.9, // 35
  3.1, 4.1, 5.1, // 36
};

const PylithScalar pylith::faults::CohesiveKinDataTet4f::_jacobianLumped[] = {
  1.1, 1.1, 1.1,
  1.2, 1.2, 1.2, // 3
  1.3, 1.3, 1.3, // 4
  1.4, 1.4, 1.4, // 5
  1.5, 1.5, 1.5,
  1.6, 1.6, 1.6, // 7
  1.8, 1.8, 1.8, // 8
  1.0, 1.0, 1.0, // 9
  1.0, 1.0, 1.0, // 34
  1.0, 1.0, 1.0, // 35
  1.0, 1.0, 1.0, // 36
};

const PylithScalar pylith::faults::CohesiveKinDataTet4f::_orientation[] = {
  0.0, -1.0, 0.0,    0.0, 0.0, +1.0,    -1.0, 0.0, 0.0,
  0.0, -1.0, 0.0,    0.0, 0.0, +1.0,    -1.0, 0.0, 0.0,
  0.0, -1.0, 0.0,    0.0, 0.0, +1.0,    -1.0, 0.0, 0.0,
};

const PylithScalar pylith::faults::CohesiveKinDataTet4f::_area[] = {
  1.0/3.0,
  1.0/3.0,
  1.0/3.0,
};

const int pylith::faults::CohesiveKinDataTet4f::_numFaultVertices = 3;
const int pylith::faults::CohesiveKinDataTet4f::_verticesFault[] = {
   2,  3,  4
};
const int pylith::faults::CohesiveKinDataTet4f::_edgesLagrange[] = {
  34, 35, 36
};
const int pylith::faults::CohesiveKinDataTet4f::_verticesNegative[] = {
   4,  5,  6
};
const int pylith::faults::CohesiveKinDataTet4f::_verticesPositive[] = {
   8,  9,  10
};

const int pylith::faults::CohesiveKinDataTet4f::_numCohesiveCells = 1;
const int pylith::faults::CohesiveKinDataTet4f::_cellMappingFault[] = {
  5
};
const int pylith::faults::CohesiveKinDataTet4f::_cellMappingCohesive[] = {
  2
};


const PylithScalar pylith::faults::CohesiveKinDataTet4f::_residual[] = {
  0.0,  0.0,  0.0,
   +7.7/3.0,  +8.7/3.0,  +9.7/3.0, // 3
  +7.9/3.0,  +8.9/3.0,  +9.9/3.0, // 4
  +7.1/3.0,  +8.1/3.0,  +9.1/3.0, // 5
  0.0,  0.0,  0.0,
  -7.7/3.0,  -8.7/3.0,  -9.7/3.0, // 7
  -7.9/3.0,  -8.9/3.0,  -9.9/3.0, // 8
  -7.1/3.0,  -8.1/3.0,  -9.1/3.0, // 9
  -1.0/3.0*(7.6-7.2 +  0.07938069066),
  -1.0/3.0*(8.6-8.2 +  1.82575588523),
  -1.0/3.0*(9.6-9.2 +  0.55566483464), // 34
  -1.0/3.0*(7.8-7.3 +  0.14140241667),
  -1.0/3.0*(8.8-8.3 +  1.69682900001),
  -1.0/3.0*(9.8-9.3 +  0.56560966667), // 35
  -1.0/3.0*(7.0-7.4 +  0.18205179147),
  -1.0/3.0*(8.0-8.4 +  1.51709826228),
  -1.0/3.0*(9.0-9.4 +  0.54615537442), // 36
};

const PylithScalar pylith::faults::CohesiveKinDataTet4f::_jacobian[] = {
  0.0, 0.0, 0.0, // 2x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 2y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 2z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 3x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 -1.0/3.0, 0.0, 0.0, // 10
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 3y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,-1.0/3.0, 0.0, // 10
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 3z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,-1.0/3.0, // 10
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 4x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 -1.0/3.0, 0.0, 0.0, // 11
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 4y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,-1.0/3.0, 0.0, // 11
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 4z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,-1.0/3.0, // 11
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 5x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 -1.0/3.0, 0.0, 0.0, // 12
  0.0, 0.0, 0.0, // 5y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,-1.0/3.0, 0.0, // 12
  0.0, 0.0, 0.0, // 5z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,-1.0/3.0, // 12
  0.0, 0.0, 0.0, // 6x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 6y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 6z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 7x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 +1.0/3.0, 0.0, 0.0, // 10
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 7y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,+1.0/3.0, 0.0, // 10
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 7z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,+1.0/3.0, // 10
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 8x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 +1.0/3.0, 0.0, 0.0, // 11
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 8y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,+1.0/3.0, 0.0, // 11
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 8z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,+1.0/3.0, // 11
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 9x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 +1.0/3.0, 0.0, 0.0, // 12
  0.0, 0.0, 0.0, // 9y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,+1.0/3.0, 0.0, // 12
  0.0, 0.0, 0.0, // 9z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,+1.0/3.0, // 12
  0.0, 0.0, 0.0, // 10x
 -1.0/3.0, 0.0, 0.0, // 3
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 +1.0/3.0, 0.0, 0.0, // 7
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 10y
  0.0,-1.0/3.0, 0.0, // 3
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,+1.0/3.0, 0.0, // 7
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 10z
  0.0, 0.0,-1.0/3.0, // 3
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,+1.0/3.0, // 7
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 11x
  0.0, 0.0, 0.0,
 -1.0/3.0, 0.0, 0.0, // 4
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 +1.0/3.0, 0.0, 0.0, // 8
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 11y
  0.0, 0.0, 0.0,
  0.0,-1.0/3.0, 0.0, // 4
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,+1.0/3.0, 0.0, // 8
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 11z
  0.0, 0.0, 0.0,
  0.0, 0.0,-1.0/3.0, // 4
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,+1.0/3.0, // 8
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 12x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 -1.0/3.0, 0.0, 0.0, // 5
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 +1.0/3.0, 0.0, 0.0, // 9
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 12y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,-1.0/3.0, 0.0, // 5
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,+1.0/3.0, 0.0, // 9
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 12z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,-1.0/3.0, // 5
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,+1.0/3.0, // 9
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
};

pylith::faults::CohesiveKinDataTet4f::CohesiveKinDataTet4f(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  quadPts = const_cast<PylithScalar*>(_quadPts);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDeriv = const_cast<PylithScalar*>(_basisDeriv);
  verticesRef = const_cast<PylithScalar*>(_verticesRef);
  id = _id;
  label = const_cast<char*>(_label);
  finalSlipFilename = const_cast<char*>(_finalSlipFilename);
  slipTimeFilename = const_cast<char*>(_slipTimeFilename);
  riseTimeFilename = const_cast<char*>(_riseTimeFilename);
  fieldT = const_cast<PylithScalar*>(_fieldT);
  fieldIncr = const_cast<PylithScalar*>(_fieldIncr);
  jacobianLumped = const_cast<PylithScalar*>(_jacobianLumped);
  orientation = const_cast<PylithScalar*>(_orientation);
  area = const_cast<PylithScalar*>(_area);
  residual = const_cast<PylithScalar*>(_residual);
  jacobian = const_cast<PylithScalar*>(_jacobian);
  verticesFault = const_cast<int*>(_verticesFault);
  edgesLagrange = const_cast<int*>(_edgesLagrange);
  verticesNegative = const_cast<int*>(_verticesNegative);
  verticesPositive = const_cast<int*>(_verticesPositive);
  numFaultVertices = _numFaultVertices;  
  cellMappingFault = const_cast<int*>(_cellMappingFault);
  cellMappingCohesive = const_cast<int*>(_cellMappingCohesive);
  numCohesiveCells = _numCohesiveCells;  
} // constructor

pylith::faults::CohesiveKinDataTet4f::~CohesiveKinDataTet4f(void)
{}


// End of file
