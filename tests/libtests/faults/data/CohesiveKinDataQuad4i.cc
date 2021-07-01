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
 */


#include "CohesiveKinDataQuad4i.hh"

const char* pylith::faults::CohesiveKinDataQuad4i::_meshFilename = "data/quad4i.mesh";

const int pylith::faults::CohesiveKinDataQuad4i::_spaceDim = 2;

const int pylith::faults::CohesiveKinDataQuad4i::_cellDim = 1;

const int pylith::faults::CohesiveKinDataQuad4i::_numBasis = 2;

const int pylith::faults::CohesiveKinDataQuad4i::_numQuadPts = 2;

const PylithScalar pylith::faults::CohesiveKinDataQuad4i::_quadPts[2] = {
  -1.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveKinDataQuad4i::_quadWts[2] = {
  1.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveKinDataQuad4i::_basis[2*2] = {
  1.0, 0.0,
  0.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveKinDataQuad4i::_basisDeriv[2*2] = {
  -0.5, 0.5,
  -0.5, 0.5,
};

const PylithScalar pylith::faults::CohesiveKinDataQuad4i::_verticesRef[2] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveKinDataQuad4i::_id = 100;

const char* pylith::faults::CohesiveKinDataQuad4i::_label = "fault";
const char* pylith::faults::CohesiveKinDataQuad4i::_edge = "edge";

const char* pylith::faults::CohesiveKinDataQuad4i::_finalSlipFilename = 
  "data/quad4i_finalslip.spatialdb";

const char* pylith::faults::CohesiveKinDataQuad4i::_slipTimeFilename = 
  "data/quad4i_sliptime.spatialdb";

const char* pylith::faults::CohesiveKinDataQuad4i::_riseTimeFilename = 
  "data/quad4i_risetime.spatialdb";

const PylithScalar pylith::faults::CohesiveKinDataQuad4i::_fieldT[(14+2)*2] = {
  3.1, 5.1, // 8
  3.2, 5.2, // 9
  3.3, 5.3, // 10
  3.4, 5.4, // 11
  3.5, 5.5, // 12
  3.6, 5.6, // 13
  3.7, 5.7, // 14
  3.8, 5.8, // 15
  3.9, 5.9, // 16
  3.0, 5.0, // 17
  4.2, 6.2, // 18
  4.4, 6.4, // 19
  4.1, 6.1, // 20
  4.3, 6.3, // 21
  4.5, 6.5, // 41
  4.7, 6.6, // 42
};



const PylithScalar pylith::faults::CohesiveKinDataQuad4i::_fieldIncr[(14+2)*2] = {
  6.1, 4.1, // 8
  6.2, 4.2, // 9
  6.3, 4.3, // 10
  6.4, 4.4, // 11
  6.5, 4.5, // 12
  6.6, 4.6, // 13
  6.7, 4.7, // 14
  6.8, 4.8, // 15
  6.9, 4.9, // 16
  6.0, 4.0, // 17
  5.2, 3.2, // 18
  5.4, 3.4, // 19
  5.1, 3.1, // 20
  5.3, 3.3, // 21
  5.5, 3.5, // 41
  5.7, 3.7, // 42
};


const PylithScalar pylith::faults::CohesiveKinDataQuad4i::_jacobianLumped[(14+2)*2] = {
  1.1, 7.1, // 8
  1.2, 7.2, // 9
  1.3, 7.3, // 10
  1.4, 7.4, // 11
  1.5, 7.5, // 12
  1.6, 7.6, // 13
  1.7, 7.7, // 14
  1.8, 7.8, // 15
  1.9, 7.9, // 16
  1.0, 7.0, // 17
  2.2, 3.2, // 18
  2.4, 3.4, // 19
  2.6, 3.5, // 20
  2.8, 3.6, // 21
  3.0, 3.7, // 41
  1.0, 1.0, // 42
};



const PylithScalar pylith::faults::CohesiveKinDataQuad4i::_orientation[3*2*2] = {
  0.0, +1.0,   1.0, 0.0,
  0.0, +1.0,   1.0, 0.0,
  0.0, +1.0,   1.0, 0.0,
};

const PylithScalar pylith::faults::CohesiveKinDataQuad4i::_area[3] = {
   5.0,
  10.0,
   5.0,
};

const int pylith::faults::CohesiveKinDataQuad4i::_numFaultVertices = 3;
const int pylith::faults::CohesiveKinDataQuad4i::_verticesFault[3] = {
  4,  5,  6,
};
const int pylith::faults::CohesiveKinDataQuad4i::_edgesLagrange[3] = {
  43, 41, 42,
};
const int pylith::faults::CohesiveKinDataQuad4i::_verticesNegative[3] = {
  12, 15, 18,
};
const int pylith::faults::CohesiveKinDataQuad4i::_verticesPositive[3] = {
  12, 20, 21,
};

const int pylith::faults::CohesiveKinDataQuad4i::_numCohesiveCells = 2;
const int pylith::faults::CohesiveKinDataQuad4i::_cellMappingFault[2] = {
  7, 8,
};
const int pylith::faults::CohesiveKinDataQuad4i::_cellMappingCohesive[2] = {
  6, 7,
};


const PylithScalar pylith::faults::CohesiveKinDataQuad4i::_residual[(14+2)*2] = {
  0.0,  0.0, // 8
  0.0,  0.0, // 9
  0.0,  0.0, // 10
  0.0,  0.0, // 11
  0.0,  0.0, // 12
  0.0,  0.0, // 13
  0.0,  0.0, // 14
 +10.0*4.5, +10.0*6.5, // 15
  0.0,  0.0, // 16
  0.0,  0.0, // 17
 +5.0*4.7, +5.0*6.6, // 18
  0.0,  0.0, // 19
 -10.0*4.5, -10.0*6.5, // 20
 -5.0*4.7, -5.0*6.6, // 21
 -10.0*(4.1-3.8 - 0.14794836271),
 -10.0*(6.1-5.8 - 1.77538035254), // 41
  -5.0*(4.3-4.2 - 0.19186497837),
  -5.0*(6.3-6.2 - 1.59887481971), // 42
};

const PylithScalar pylith::faults::CohesiveKinDataQuad4i::_jacobian[(14+2)*2*(14+2)*2] = {
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0,
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0,
-10.0, 0.0, // 41
  0.0, 0.0,
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
  0.0, 0.0,
  0.0,-10.0, // 41
  0.0, 0.0,
  0.0, 0.0, // 16x
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
  0.0, 0.0,
  0.0, 0.0, // 16y
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
  0.0, 0.0,
  0.0, 0.0, // 17x
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
  0.0, 0.0,
  0.0, 0.0, // 17y
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
  0.0, 0.0,
  0.0, 0.0, // 18x
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
 -5.0, 0.0, // 42
  0.0, 0.0, // 18y
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
  0.0,-5.0, // 42
  0.0, 0.0, // 19x
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
  0.0, 0.0,
  0.0, 0.0, // 19y
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
  0.0, 0.0,
  0.0, 0.0, // 20x
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
+10.0, 0.0, // 41
  0.0, 0.0,
  0.0, 0.0, // 20y
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
  0.0,+10.0, // 41
  0.0, 0.0,
  0.0, 0.0, // 21x
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
 +5.0, 0.0, // 42
  0.0, 0.0, // 21y
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
  0.0,+5.0, // 42
  0.0, 0.0, // 41x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
-10.0, 0.0, // 15
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +10.0, 0.0, // 20
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 41y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-10.0, // 15
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+10.0, // 20
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 42x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -5.0, 0.0, // 18
  0.0, 0.0,
  0.0, 0.0,
 +5.0, 0.0, // 21
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 42y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-5.0, // 18
  0.0, 0.0,
  0.0, 0.0,
  0.0,+5.0, // 21
  0.0, 0.0,
  0.0, 0.0,
};

pylith::faults::CohesiveKinDataQuad4i::CohesiveKinDataQuad4i(void)
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
  edge = const_cast<char*>(_edge);
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

pylith::faults::CohesiveKinDataQuad4i::~CohesiveKinDataQuad4i(void)
{}


// End of file
