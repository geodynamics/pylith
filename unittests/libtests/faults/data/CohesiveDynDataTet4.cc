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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
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
 * Cells are 0-1,2, vertices are 3-10,11-13.
 *
 * 3   4,5,6  8,9,10   7
 *             34,35,36
 *     ^^^^^^^^^^^^ Cohesive element in x-y plane.
 */

#include "CohesiveDynDataTet4.hh"

const char* pylith::faults::CohesiveDynDataTet4::_meshFilename =
  "data/tet4.mesh";

const int pylith::faults::CohesiveDynDataTet4::_spaceDim = 3;

const int pylith::faults::CohesiveDynDataTet4::_cellDim = 2;

const int pylith::faults::CohesiveDynDataTet4::_numBasis = 3;

const int pylith::faults::CohesiveDynDataTet4::_numQuadPts = 3;

const PylithScalar pylith::faults::CohesiveDynDataTet4::_quadPts[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
};

const PylithScalar pylith::faults::CohesiveDynDataTet4::_quadWts[] = {
  2.0/3.0, 2.0/3.0, 2.0/3.0,
};

const PylithScalar pylith::faults::CohesiveDynDataTet4::_basis[] = {
  1.0, 0.0, 0.0,
  0.0, 1.0, 0.0,
  0.0, 0.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveDynDataTet4::_basisDeriv[] = {
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

const PylithScalar pylith::faults::CohesiveDynDataTet4::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
};

const int pylith::faults::CohesiveDynDataTet4::_id = 10;

const char* pylith::faults::CohesiveDynDataTet4::_label = "fault";

const char* pylith::faults::CohesiveDynDataTet4::_initialTractFilename = 
  "data/tet4_initialtract.spatialdb";

const PylithScalar pylith::faults::CohesiveDynDataTet4::_fieldT[] = {
  7.1, 8.1, 9.1,
  7.2, 8.2, 9.2, // 4
  7.3, 8.3, 9.3, // 5
  7.4, 8.4, 9.4, // 6
  7.5, 8.5, 9.5,
  7.2, 8.2, 9.2, // 8
  7.3, 8.3, 9.3, // 9
  7.4, 8.4, 9.4, // 10
 -7.7, 18.7, 19.7, // 34
 -7.9, 18.9, 19.9, // 35
 -7.1, 18.1, 19.1, // 36
};

const PylithScalar pylith::faults::CohesiveDynDataTet4::_jacobian[] = {
  1.0,  0.1,  0.2, // 2x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  1.2,  1.3,  1.4,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  2.0,  4.1,  3.2, // 2y
  2.3,  4.4,  3.5,
  2.6,  4.7,  3.8,
  2.9,  4.0,  3.1,
  2.2,  4.3,  3.4,
  2.5,  4.6,  3.7,
  2.8,  4.9,  3.0,
  2.1,  4.2,  3.3,
  2.4,  4.5,  3.6,
  2.7,  4.8,  3.9,
  2.0,  4.1,  3.2,
  1.0,  7.1,  2.2, // 2z
  1.3,  7.4,  2.5,
  1.6,  7.7,  2.8,
  1.9,  7.0,  2.1,
  1.2,  7.3,  2.4,
  1.5,  7.6,  2.7,
  1.8,  7.9,  2.0,
  1.1,  7.2,  2.3,
  1.4,  7.5,  2.6,
  1.7,  7.8,  2.9,
  1.0,  7.1,  2.2,
  1.0,  0.1,  0.2, // 3x
 +4.0, -1.1, -1.2, // 3
 -1.3, -1.4, -1.5, // 4
 -1.6, -1.7, -1.8, // 5
  0.3,  0.4,  0.5,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  2.0,  4.1,  3.2, // 3y
 -1.1, +4.1, -2.3, // 3
 -2.4, -2.5, -2.6, // 4
 -2.7, -2.8, -2.9, // 5
  2.3,  4.4,  3.5,
  2.5,  4.6,  3.7,
  2.8,  4.9,  3.0,
  2.1,  4.2,  3.3,
  2.4,  4.5,  3.6,
  2.7,  4.8,  3.9,
  2.0,  4.1,  3.2,
  1.0,  7.1,  2.2, // 3z
 -1.2, -2.3, +4.2, // 3
 -1.0, -1.1, -1.2, // 4
 -1.3, -1.4, -1.5, // 5
  1.5,  7.6,  2.7,
  1.3,  7.4,  2.5,
  1.8,  7.9,  2.0,
  1.1,  7.2,  2.3,
  1.4,  7.5,  2.6,
  1.7,  7.8,  2.9,
  1.0,  7.1,  2.2,
  1.0,  0.1,  0.2, // 4x
 -1.3, -2.4, -1.0, // 3
 +4.3, -0.2, -0.3, // 4
 -0.4, -0.5, -0.6, // 5
  0.3,  0.4,  0.5,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  2.0,  4.1,  3.2, // 4y
 -1.4, -2.5, -1.1, // 3
 -0.2, +4.4, -0.9, // 4
 -0.8, -0.7, -0.5, // 5
  2.3,  4.4,  3.5,
  2.5,  4.6,  3.7,
  2.8,  4.9,  3.0,
  2.1,  4.2,  3.3,
  2.4,  4.5,  3.6,
  2.7,  4.8,  3.9,
  2.0,  4.1,  3.2,
  1.0,  7.1,  2.2, // 4z
 -1.5, -2.6, -1.2, // 3
 -0.3, -0.9, +4.5, // 4
 -1.1, -1.2, -1.3, // 5
  1.3,  7.4,  2.5,
  1.5,  7.6,  2.7,
  1.8,  7.9,  2.0,
  1.1,  7.2,  2.3,
  1.4,  7.5,  2.6,
  1.7,  7.8,  2.9,
  1.0,  7.1,  2.2,
  1.0,  0.1,  0.2, // 5x
 -1.6, -2.7, -1.3, // 3
 -0.4, -0.8, -1.1, // 4
 +4.6, -1.8, -1.5, // 5
  0.3,  0.4,  0.5,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  1.0,  0.1,  0.2, // 5y
 -1.7, -2.8, -1.4, // 3
 -0.5, -0.7, -1.2, // 4
 -1.8, +4.7, -1.1, // 5
  0.3,  0.4,  0.5,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  1.0,  7.1,  2.2, // 5z
 -1.8, -2.9, -1.5, // 3
 -0.6, -0.5, -1.3, // 4
 -1.5, -1.1, +4.8, // 5
  1.3,  7.4,  2.5,
  1.5,  7.6,  2.7,
  1.8,  7.9,  2.0,
  1.1,  7.2,  2.3,
  1.4,  7.5,  2.6,
  1.7,  7.8,  2.9,
  1.0,  7.1,  2.2,
  1.0,  0.1,  0.2, // 6x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  1.2,  1.3,  1.4,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  2.0,  4.1,  3.2, // 6y
  2.3,  4.4,  3.5,
  2.6,  4.7,  3.8,
  2.9,  4.0,  3.1,
  2.2,  4.3,  3.4,
  2.5,  4.6,  3.7,
  2.8,  4.9,  3.0,
  2.1,  4.2,  3.3,
  2.4,  4.5,  3.6,
  2.7,  4.8,  3.9,
  2.0,  4.1,  3.2,
  1.0,  7.1,  2.2, // 6z
  1.3,  7.4,  2.5,
  1.6,  7.7,  2.8,
  1.9,  7.0,  2.1,
  1.2,  7.3,  2.4,
  1.5,  7.6,  2.7,
  1.8,  7.9,  2.0,
  1.1,  7.2,  2.3,
  1.4,  7.5,  2.6,
  1.7,  7.8,  2.9,
  1.0,  7.1,  2.2,
  1.0,  0.1,  0.2, // 7x
  0.3,  0.4,  0.5,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
 +5.0, -1.1, -1.2, // 3
 -1.3, -1.4, -1.5, // 4
 -1.6, -1.7, -1.8, // 5
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  2.0,  4.1,  3.2, // 7y
  2.3,  4.4,  3.5,
  2.5,  4.6,  3.7,
  2.8,  4.9,  3.0,
  2.1,  4.2,  3.3,
 -1.1, +5.1, -2.3, // 7
 -2.4, -2.5, -2.6, // 8
 -2.7, -2.8, -2.9, // 9
  2.4,  4.5,  3.6,
  2.7,  4.8,  3.9,
  2.0,  4.1,  3.2,
  1.0,  7.1,  2.2, // 7z
  1.5,  7.6,  2.7,
  1.8,  7.9,  2.0,
  1.1,  7.2,  2.3,
  1.3,  7.4,  2.5,
 -1.2, -2.3, +5.2, // 7
 -1.0, -1.1, -1.2, // 8
 -1.3, -1.4, -1.5, // 9
  1.4,  7.5,  2.6,
  1.7,  7.8,  2.9,
  1.0,  7.1,  2.2,
  1.0,  0.1,  0.2, // 8x
  0.3,  0.4,  0.5,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
 -1.3, -2.4, -1.0, // 7
 +5.3, -0.2, -0.3, // 8
 -0.4, -0.5, -0.6, // 9
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  2.0,  4.1,  3.2, // 8y
  2.3,  4.4,  3.5,
  2.5,  4.6,  3.7,
  2.8,  4.9,  3.0,
  2.1,  4.2,  3.3,
 -1.4, -2.5, -1.1, // 7
 -0.2, +5.4, -0.9, // 8
 -0.8, -0.7, -0.5, // 9
  2.4,  4.5,  3.6,
  2.7,  4.8,  3.9,
  2.0,  4.1,  3.2,
  1.0,  7.1,  2.2, // 8z
  1.3,  7.4,  2.5,
  1.5,  7.6,  2.7,
  1.8,  7.9,  2.0,
  1.1,  7.2,  2.3,
 -1.5, -2.6, -1.2, // 7
 -0.3, -0.9, +5.5, // 8
 -1.1, -1.2, -1.3, // 9
  1.4,  7.5,  2.6,
  1.7,  7.8,  2.9,
  1.0,  7.1,  2.2,
  1.0,  0.1,  0.2, // 9x
  0.3,  0.4,  0.5,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
 -1.6, -2.7, -1.3, // 7
 -0.4, -0.8, -1.1, // 8
 +5.6, -1.8, -1.5, // 9
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  2.0,  4.1,  3.2, // 9y
  2.3,  4.4,  3.5,
  2.5,  4.6,  3.7,
  2.8,  4.9,  3.0,
  2.1,  4.2,  3.3,
 -1.7, -2.8, -1.4, // 7
 -0.5, -0.7, -1.2, // 8
 -1.8, +5.7, -1.1, // 9
  2.4,  4.5,  3.6,
  2.7,  4.8,  3.9,
  2.0,  4.1,  3.2,
  1.0,  7.1,  2.2, // 9z
  1.3,  7.4,  2.5,
  1.5,  7.6,  2.7,
  1.8,  7.9,  2.0,
  1.1,  7.2,  2.3,
 -1.8, -2.9, -1.5, // 7
 -0.6, -0.5, -1.3, // 8
 -1.5, -1.1, +5.8, // 9
  1.4,  7.5,  2.6,
  1.7,  7.8,  2.9,
  1.0,  7.1,  2.2,
  1.0,  0.1,  0.2, // 10x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  1.2,  1.3,  1.4,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  2.0,  4.1,  3.2, // 10y
  2.3,  4.4,  3.5,
  2.6,  4.7,  3.8,
  2.9,  4.0,  3.1,
  2.2,  4.3,  3.4,
  2.5,  4.6,  3.7,
  2.8,  4.9,  3.0,
  2.1,  4.2,  3.3,
  2.4,  4.5,  3.6,
  2.7,  4.8,  3.9,
  2.0,  4.1,  3.2,
  1.0,  7.1,  2.2, // 10z
  1.3,  7.4,  2.5,
  1.6,  7.7,  2.8,
  1.9,  7.0,  2.1,
  1.2,  7.3,  2.4,
  1.5,  7.6,  2.7,
  1.8,  7.9,  2.0,
  1.1,  7.2,  2.3,
  1.4,  7.5,  2.6,
  1.7,  7.8,  2.9,
  1.0,  7.1,  2.2,
  1.0,  0.1,  0.2, // 11x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  1.2,  1.3,  1.4,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  2.0,  4.1,  3.2, // 11y
  2.3,  4.4,  3.5,
  2.6,  4.7,  3.8,
  2.9,  4.0,  3.1,
  2.2,  4.3,  3.4,
  2.5,  4.6,  3.7,
  2.8,  4.9,  3.0,
  2.1,  4.2,  3.3,
  2.4,  4.5,  3.6,
  2.7,  4.8,  3.9,
  2.0,  4.1,  3.2,
  1.0,  7.1,  2.2, // 11z
  1.3,  7.4,  2.5,
  1.6,  7.7,  2.8,
  1.9,  7.0,  2.1,
  1.2,  7.3,  2.4,
  1.5,  7.6,  2.7,
  1.8,  7.9,  2.0,
  1.1,  7.2,  2.3,
  1.4,  7.5,  2.6,
  1.7,  7.8,  2.9,
  1.0,  7.1,  2.2,
  1.0,  0.1,  0.2, // 12x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  1.2,  1.3,  1.4,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  2.0,  4.1,  3.2, // 12y
  2.3,  4.4,  3.5,
  2.6,  4.7,  3.8,
  2.9,  4.0,  3.1,
  2.2,  4.3,  3.4,
  2.5,  4.6,  3.7,
  2.8,  4.9,  3.0,
  2.1,  4.2,  3.3,
  2.4,  4.5,  3.6,
  2.7,  4.8,  3.9,
  2.0,  4.1,  3.2,
  1.0,  7.1,  2.2, // 12z
  1.3,  7.4,  2.5,
  1.6,  7.7,  2.8,
  1.9,  7.0,  2.1,
  1.2,  7.3,  2.4,
  1.5,  7.6,  2.7,
  1.8,  7.9,  2.0,
  1.1,  7.2,  2.3,
  1.4,  7.5,  2.6,
  1.7,  7.8,  2.9,
  1.0,  7.1,  2.2,
};

// ----------------------------------------------------------------------
// Computed values
// ----------------------------------------------------------------------

const PylithScalar pylith::faults::CohesiveDynDataTet4::_orientation[] = {
  0.0, +1.0, 0.0,    0.0, 0.0, +1.0,    +1.0, 0.0, 0.0,
  0.0, +1.0, 0.0,    0.0, 0.0, +1.0,    +1.0, 0.0, 0.0,
  0.0, +1.0, 0.0,    0.0, 0.0, +1.0,    +1.0, 0.0, 0.0,
};

const PylithScalar pylith::faults::CohesiveDynDataTet4::_area[] = {
  1.0/3.0, 
  1.0/3.0, 
  1.0/3.0,
};

const PylithScalar pylith::faults::CohesiveDynDataTet4::_initialTractions[] = {
  // Fault coordinate frame
  +1.0, +2.0, -3.0,
  +1.1, +2.1, -3.1,
  +1.2, +2.2, -3.2,
};


const int pylith::faults::CohesiveDynDataTet4::_numConstraintEdges = 3;
const int pylith::faults::CohesiveDynDataTet4::_constraintEdges[] = {
  34, 35, 36
};
const int pylith::faults::CohesiveDynDataTet4::_negativeVertices[] = {
   4,  5,  6
};

// ----------------------------------------------------------------------
// Stick case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataTet4::_fieldIncrStick[] = {
  1.1, 2.1, 3.1,
  1.2, 2.2, 3.2, // 4
  1.3, 2.3, 3.3, // 5
  1.4, 2.4, 3.4, // 6
  1.5, 2.5, 3.5,
  1.2, 2.2, 3.2, // 8
  1.3, 2.3, 3.3, // 9
  1.4, 2.4, 3.4, // 10
 -81.7, 2.7, 3.7, // 34
 -81.9, 2.9, 3.9, // 35
 -81.1, 2.1, 3.1, // 36
};

// No slip
const PylithScalar pylith::faults::CohesiveDynDataTet4::_slipStickE[] = {
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
};

// ----------------------------------------------------------------------
// Slip case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataTet4::_fieldIncrSlip[] = {
  1.1, 2.1, 3.1,
  1.2, 2.2, 3.2, // 4
  1.3, 2.3, 3.3, // 5
  1.4, 2.4, 3.4, // 6
  1.5, 2.5, 3.5,
  1.2, 2.2, 3.2, // 7
  1.3, 2.3, 3.3, // 8
  1.4, 2.4, 3.4, // 10
 -4.7, 5.7, 6.7, // 34
 -4.9, 5.9, 6.9, // 35
 -4.1, 5.1, 6.1, // 36
};

// Output
const PylithScalar pylith::faults::CohesiveDynDataTet4::_fieldIncrSlipE[] = {
   1.100000000000,   2.100000000000,   3.100000000000,
   1.200000000000,   2.200000916680,   3.200000185945,
   1.300000000000,   2.299999954294,   3.300000084646,
   1.400000000000,   2.400000465197,   3.400000369651,
   1.500000000000,   2.500000000000,   3.500000000000,
   1.200000000000,   2.199999083320,   3.199999814055,
   1.300000000000,   2.300000045706,   3.299999915354,
   1.400000000000,   2.399999534803,   3.399999630349,
  -4.700000000000, -13.650158708856, -14.236237291549,
  -4.900000000000, -13.683824215808, -14.263164878373,
  -4.100000000000, -13.548480328050, -14.156107942537,
};

const PylithScalar pylith::faults::CohesiveDynDataTet4::_slipSlipE[] = {
  -0.000001833360,  -0.000000371890,   0.000000000000,
   0.000000091413,  -0.000000169291,   0.000000000000,
  -0.000000930394,  -0.000000739303,   0.000000000000,
};

// ----------------------------------------------------------------------
// Open case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataTet4::_fieldIncrOpen[] = {
  1.1, 2.1, 3.1,
  1.2, 2.2, 3.2, // 4
  1.3, 2.3, 3.3, // 5
  1.4, 2.4, 3.4, // 6
  1.5, 2.5, 3.5,
  1.2, 2.2, 3.2, // 8
  1.3, 2.3, 3.3, // 9
  1.4, 2.4, 3.4, // 10
 +80.7, 2.7, 3.7, // 34
 +80.9, 2.9, 3.9, // 35
 +80.1, 2.1, 3.1, // 36
};

// Output
const PylithScalar pylith::faults::CohesiveDynDataTet4::_fieldIncrOpenE[] = {
   1.100000000000,   2.100000000000,   3.100000000000,
   1.199999161233,   2.200004202086,   3.200002477003,
   1.299998110594,   2.300001999917,   3.300002482177,
   1.400000000000,   2.400002588150,   3.400002474071,
   1.500000000000,   2.500000000000,   3.500000000000,
   1.200000838767,   2.199995797914,   3.199997522997,
   1.300001889406,   2.299998000083,   3.299997517823,
   1.400000000000,   2.399997411850,   3.399997525929,
   7.700000000000, -18.700000000000, -19.700000000000,
   7.900000000000, -18.900000000000, -19.900000000000,
   7.100000000000, -18.100000000000, -19.100000000000,
};

const PylithScalar pylith::faults::CohesiveDynDataTet4::_slipOpenE[] = {
  -0.000008404171,  -0.000004954006,   0.000001677534,
  -0.000003999835,  -0.000004964354,   0.000003778811,
  -0.000005176299,  -0.000004948142,   0.000000000000,
};

// ----------------------------------------------------------------------
pylith::faults::CohesiveDynDataTet4::CohesiveDynDataTet4(void)
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
  initialTractFilename = const_cast<char*>(_initialTractFilename);

  fieldT = const_cast<PylithScalar*>(_fieldT);
  jacobian = const_cast<PylithScalar*>(_jacobian);
  orientation = const_cast<PylithScalar*>(_orientation);
  area = const_cast<PylithScalar*>(_area);
  initialTractions = const_cast<PylithScalar*>(_initialTractions);

  constraintEdges = const_cast<int*>(_constraintEdges);
  negativeVertices = const_cast<int*>(_negativeVertices);
  numConstraintEdges = _numConstraintEdges;  

  // Stick
  fieldIncrStick = const_cast<PylithScalar*>(_fieldIncrStick);
  slipStickE = const_cast<PylithScalar*>(_slipStickE);

  // Slip
  fieldIncrSlip = const_cast<PylithScalar*>(_fieldIncrSlip);
  fieldIncrSlipE = const_cast<PylithScalar*>(_fieldIncrSlipE);
  slipSlipE = const_cast<PylithScalar*>(_slipSlipE);

  // Open
  fieldIncrOpen = const_cast<PylithScalar*>(_fieldIncrOpen);
  fieldIncrOpenE = const_cast<PylithScalar*>(_fieldIncrOpenE);
  slipOpenE = const_cast<PylithScalar*>(_slipOpenE);
} // constructor

pylith::faults::CohesiveDynDataTet4::~CohesiveDynDataTet4(void)
{}


// End of file
