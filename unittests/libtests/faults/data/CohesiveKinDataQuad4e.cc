// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/* Original mesh
 *
 * Cells are 0-3, vertices are 4-12.
 *
 *      10 --------11 --------12
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       5 -------- 7 -------- 9
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       4 -------- 6 -------- 8
 *
 * After adding cohesive elements
 *
 * Cells are 0-3,16-17 vertices are 4-15.
 *
 *      10 --------11-18-15 --------12
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       5 -------- 7-17-14 -------- 9
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       4 -------- 6-16-13 -------- 8
 */


#include "CohesiveKinDataQuad4e.hh"

const char* pylith::faults::CohesiveKinDataQuad4e::_meshFilename =
  "data/quad4e.mesh";

const int pylith::faults::CohesiveKinDataQuad4e::_spaceDim = 2;

const int pylith::faults::CohesiveKinDataQuad4e::_cellDim = 1;

const int pylith::faults::CohesiveKinDataQuad4e::_numBasis = 2;

const int pylith::faults::CohesiveKinDataQuad4e::_numQuadPts = 1;

const double pylith::faults::CohesiveKinDataQuad4e::_quadPts[] = {
  0.0,
};

const double pylith::faults::CohesiveKinDataQuad4e::_quadWts[] = {
  2.0,
};

const double pylith::faults::CohesiveKinDataQuad4e::_basis[] = {
  0.5,
  0.5
};

const double pylith::faults::CohesiveKinDataQuad4e::_basisDeriv[] = {
  -0.5,
   0.5
};

const double pylith::faults::CohesiveKinDataQuad4e::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveKinDataQuad4e::_id = 10;

const char* pylith::faults::CohesiveKinDataQuad4e::_label = "fault";

const char* pylith::faults::CohesiveKinDataQuad4e::_finalSlipFilename = 
  "data/quad4e_finalslip.spatialdb";

const char* pylith::faults::CohesiveKinDataQuad4e::_slipTimeFilename = 
  "data/quad4e_sliptime.spatialdb";

const char* pylith::faults::CohesiveKinDataQuad4e::_riseTimeFilename = 
  "data/quad4e_risetime.spatialdb";

const double pylith::faults::CohesiveKinDataQuad4e::_fieldT[] = {
  3.1, 5.1,
  3.2, 5.2,
  3.3, 5.3, // 6
  3.4, 5.4, // 7
  3.5, 5.5,
  3.6, 5.6,
  3.7, 5.7,
  3.8, 5.8, // 11
  3.9, 5.9,
  3.0, 5.0, // 13
  4.2, 6.2, // 14
  4.4, 6.4, // 15
  4.1, 6.1, // 16
  4.3, 6.3, // 17
  4.5, 6.5, // 18
};



const double pylith::faults::CohesiveKinDataQuad4e::_fieldIncr[] = {
  6.1, 4.1,
  6.2, 4.2,
  6.3, 4.3, // 6
  6.4, 4.4, // 7
  6.5, 4.5,
  6.6, 4.6,
  6.7, 4.7,
  6.8, 4.8, // 11
  6.9, 4.9,
  6.0, 4.0, // 13
  5.2, 3.2, // 14
  5.4, 3.4, // 15
  5.1, 3.1, // 16
  5.3, 3.3, // 17
  5.5, 3.5, // 18
};


const double pylith::faults::CohesiveKinDataQuad4e::_jacobianLumped[] = {
  1.1, 7.1,
  1.2, 7.2,
  1.3, 7.3, // 6
  1.4, 7.4, // 7
  1.5, 7.5,
  1.6, 7.6,
  1.7, 7.7,
  1.8, 7.8, // 11
  1.9, 7.9,
  1.0, 7.0, // 13
  2.2, 3.2, // 14
  2.4, 3.4, // 15
  1.0, 1.0, // 16
  2.0, 2.0, // 17
  1.0, 1.0, // 18
};



const double pylith::faults::CohesiveKinDataQuad4e::_orientation[] = {
  0.0, -1.0,  -1.0, 0.0,
  0.0, -1.0,  -1.0, 0.0,
  0.0, -1.0,  -1.0, 0.0,
};

const double pylith::faults::CohesiveKinDataQuad4e::_area[] = {
  1.0,
  2.0,
  1.0,
};

const int pylith::faults::CohesiveKinDataQuad4e::_numFaultVertices = 3;
const int pylith::faults::CohesiveKinDataQuad4e::_verticesFault[] = {
  3, 2, 4
};
const int pylith::faults::CohesiveKinDataQuad4e::_verticesLagrange[] = {
  17, 16, 18
};
const int pylith::faults::CohesiveKinDataQuad4e::_verticesNegative[] = {
  7, 6, 11
};
const int pylith::faults::CohesiveKinDataQuad4e::_verticesPositive[] = {
  14, 13, 15
};

const int pylith::faults::CohesiveKinDataQuad4e::_numCohesiveCells = 2;
const int pylith::faults::CohesiveKinDataQuad4e::_cellMappingFault[] = {
  0, 1
};
const int pylith::faults::CohesiveKinDataQuad4e::_cellMappingCohesive[] = {
  19, 20
};


const double pylith::faults::CohesiveKinDataQuad4e::_residual[] = {
  0.0,  0.0,
  0.0,  0.0,
 -4.2, -6.2, // 6
 -8.6,-12.6, // 7
  0.0,  0.0,
  0.0,  0.0,
  0.0,  0.0,
 -4.4, -6.4, // 11
  0.0,  0.0,
 +4.2, +6.2, // 13
 +8.6,+12.6, // 15
 +4.4, +6.4, // 17
  0.5*(3.0-3.3 + 4.2-3.4) + 0.5*(0.14794836271 + 0.08241148423),
  0.5*(5.0-5.3 + 6.2-5.4) + 0.5*(1.77538035254 + 1.89546413727), // 16
  0.5*(3.0-3.3 + 4.2-3.4) + 0.5*(0.14794836271 + 0.08241148423) +
  0.5*(4.2-3.4 + 4.4-3.8) + 0.5*(0.08241148423 + 0.19186497837),
  0.5*(5.0-5.3 + 6.2-5.4) + 0.5*(1.77538035254 + 1.89546413727) +
  0.5*(6.2-5.4 + 6.4-5.8) + 0.5*(1.89546413727 + 1.59887481971), // 17
  0.5*(4.2-3.4 + 4.4-3.8) + 0.5*(0.08241148423 + 0.19186497837),
  0.5*(6.2-5.4 + 6.4-5.8) + 0.5*(1.89546413727 + 1.59887481971), // 18
};

const double pylith::faults::CohesiveKinDataQuad4e::_residualIncr[] = {
  0.0,  0.0,
  0.0,  0.0,
 -4.2, -6.2, // 6
 -8.6,-12.6, // 7
  0.0,  0.0,
  0.0,  0.0,
  0.0,  0.0,
 -4.4, -6.4, // 11
  0.0,  0.0,
 +4.2, +6.2, // 13
 +8.6,+12.6, // 15
 +4.4, +6.4, // 17
  0.5*(3.0-3.3 + 4.2-3.4) + 0.5*(0.14794836271 + 0.08241148423),
  0.5*(5.0-5.3 + 6.2-5.4) + 0.5*(1.77538035254 + 1.89546413727), // 16
  0.5*(3.0-3.3 + 4.2-3.4) + 0.5*(0.14794836271 + 0.08241148423) +
  0.5*(4.2-3.4 + 4.4-3.8) + 0.5*(0.08241148423 + 0.19186497837),
  0.5*(5.0-5.3 + 6.2-5.4) + 0.5*(1.77538035254 + 1.89546413727) +
  0.5*(6.2-5.4 + 6.4-5.8) + 0.5*(1.89546413727 + 1.59887481971), // 17
  0.5*(4.2-3.4 + 4.4-3.8) + 0.5*(0.08241148423 + 0.19186497837),
  0.5*(6.2-5.4 + 6.4-5.8) + 0.5*(1.89546413727 + 1.59887481971), // 18
};

const double pylith::faults::CohesiveKinDataQuad4e::_jacobian[] = {
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
  0.0, 0.0,
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +0.5, 0.0, // 16
 +0.5, 0.0, // 17
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+0.5, // 16
  0.0,+0.5, // 17
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
 +0.5, 0.0, // 16
 +1.0, 0.0, // 17
 +0.5, 0.0, // 18
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
  0.0,+0.5, // 16
  0.0,+1.0, // 17
  0.0,+0.5, // 18
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 11x
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
  0.0, 0.0,
 +0.5, 0.0, // 17
 +0.5, 0.0, // 18
  0.0, 0.0, // 11y
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
  0.0, 0.0,
  0.0,+0.5, // 17
  0.0,+0.5, // 18
  0.0, 0.0, // 12x
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
  0.0, 0.0,
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 13x
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
 -0.5, 0.0, // 16
 -0.5, 0.0, // 17
  0.0, 0.0,
  0.0, 0.0, // 13y
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
  0.0,-0.5, // 16
  0.0,-0.5, // 17
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
  0.0, 0.0,
 -0.5, 0.0, // 16
 -1.0, 0.0, // 17
 -0.5, 0.0, // 18
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
  0.0, 0.0,
  0.0,-0.5, // 16
  0.0,-1.0, // 17
  0.0,-0.5, // 18
  0.0, 0.0, // 15x
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
  0.0, 0.0,
 -0.5, 0.0, // 17
 -0.5, 0.0, // 18
  0.0, 0.0, // 15y
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
  0.0, 0.0,
  0.0,-0.5, // 17
  0.0,-0.5, // 18
  0.0, 0.0, // 16x
  0.0, 0.0,
 +0.5, 0.0, // 6
 +0.5, 0.0, // 7
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -0.5, 0.0, // 13
 -0.5, 0.0, // 14
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 16y
  0.0, 0.0,
  0.0,+0.5, // 6
  0.0,+0.5, // 7
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-0.5, // 13
  0.0,-0.5, // 14
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 17x
  0.0, 0.0,
 +0.5, 0.0, // 6
 +1.0, 0.0, // 7
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +0.5, 0.0, // 11
  0.0, 0.0,
 -0.5, 0.0, // 13
 -1.0, 0.0, // 14
 -0.5, 0.0, // 15
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 17y
  0.0, 0.0,
  0.0,+0.5, // 6
  0.0,+1.0, // 7
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+0.5, // 11
  0.0, 0.0,
  0.0,-0.5, // 13
  0.0,-1.0, // 14
  0.0,-0.5, // 15
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 18x
  0.0, 0.0,
  0.0, 0.0,
 +0.5, 0.0, // 7
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +0.5, 0.0, // 11
  0.0, 0.0,
  0.0, 0.0,
 -0.5, 0.0, // 14
 -0.5, 0.0, // 15
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 18y
  0.0, 0.0,
  0.0, 0.0,
  0.0,+0.5, // 7
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+0.5, // 11
  0.0, 0.0,
  0.0, 0.0,
  0.0,-0.5, // 14
  0.0,-0.5, // 15
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
};

pylith::faults::CohesiveKinDataQuad4e::CohesiveKinDataQuad4e(void)
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
  riseTimeFilename = const_cast<char*>(_riseTimeFilename);
  fieldT = const_cast<double*>(_fieldT);
  fieldIncr = const_cast<double*>(_fieldIncr);
  jacobianLumped = const_cast<double*>(_jacobianLumped);
  orientation = const_cast<double*>(_orientation);
  area = const_cast<double*>(_area);
  residual = const_cast<double*>(_residual);
  residualIncr = const_cast<double*>(_residualIncr);
  jacobian = const_cast<double*>(_jacobian);
  verticesFault = const_cast<int*>(_verticesFault);
  verticesLagrange = const_cast<int*>(_verticesLagrange);
  verticesNegative = const_cast<int*>(_verticesNegative);
  verticesPositive = const_cast<int*>(_verticesPositive);
  numFaultVertices = _numFaultVertices;  
  cellMappingFault = const_cast<int*>(_cellMappingFault);
  cellMappingCohesive = const_cast<int*>(_cellMappingCohesive);
  numCohesiveCells = _numCohesiveCells;  
} // constructor

pylith::faults::CohesiveKinDataQuad4e::~CohesiveKinDataQuad4e(void)
{}


// End of file
