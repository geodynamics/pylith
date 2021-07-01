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

/* Mesh: tet4.mesh
 *
 * Neumann BC on face defined by vertices 1, 2, 3.
 * This face has an area of sqrt(3)/2.
 * The unit vectors should be:
 * Normal:       ( 1/sqrt(3),  1/sqrt(3),  1/sqrt(3))
 * Horiz-shear:  (-1/sqrt(2),  1/sqrt(2),  0)
 * Vert-shear:   (-1/sqrt(6), -1/sqrt(6),  2/sqrt(6))
 * Applying unit values of horiz-shear, vert-shear, and normal
 * tractions should give x, y, and z tractions that are the sum of each
 * column, yielding values of approximately (-0.538005, 0.876209, 1.39385).
 * Integrating these tractions over the area of the triangle and distributing
 * values to the face nodes, we should get x, y, and z forces of about
 * ( -1.5530860877e-01,  2.529396817e-01,  4.0236892706e-01)
 * at vertices 1, 2, and 3.
 *
 */

#include "NeumannDataTet4.hh"

const char* pylith::bc::NeumannDataTet4::_meshFilename = 
  "data/tet4.mesh";

const int pylith::bc::NeumannDataTet4::_numBasis = 3;
const int pylith::bc::NeumannDataTet4::_numQuadPts = 1;
const PylithScalar pylith::bc::NeumannDataTet4::_quadPts[] = {
  -0.3333333333333333, -0.3333333333333333
};
const PylithScalar pylith::bc::NeumannDataTet4::_quadWts[] = {
  2.0,
};
const PylithScalar pylith::bc::NeumannDataTet4::_basis[] = {
  0.3333333333333333,
  0.3333333333333333,
  0.3333333333333333,
};
const PylithScalar pylith::bc::NeumannDataTet4::_basisDerivRef[] = {
 -0.5, -0.5,
  0.5,  0.0,
  0.0,  0.5,
};

const char* pylith::bc::NeumannDataTet4::_spatialDBFilename =
  "data/tet4_tractions.spatialdb";
const int pylith::bc::NeumannDataTet4::_id = 0;
const char* pylith::bc::NeumannDataTet4::_label = "bc3";

const int pylith::bc::NeumannDataTet4::_spaceDim = 3;
const int pylith::bc::NeumannDataTet4::_cellDim = 2;
const int pylith::bc::NeumannDataTet4::_numVertices = 3;
const int pylith::bc::NeumannDataTet4::_numCells = 1;
const int pylith::bc::NeumannDataTet4::_numCorners = 3;
/* Now vertices are renumbered in the submesh */
const int pylith::bc::NeumannDataTet4::_cells[] = {
  1 /*3*/, 2 /*4*/, 3 /*5*/,
};

const PylithScalar pylith::bc::NeumannDataTet4::_tractionsCell[] = {
  -0.5380048025,  0.87620875991,  1.3938468501
};
const PylithScalar pylith::bc::NeumannDataTet4::_valsResidual[] = {
  0.0,               0.0,               0.0,
  -1.5530860877e-01,  2.529396817e-01,  4.0236892706e-01,
  -1.5530860877e-01,  2.529396817e-01,  4.0236892706e-01,
  -1.5530860877e-01,  2.529396817e-01,  4.0236892706e-01,
  0.0,               0.0,               0.0,
};


pylith::bc::NeumannDataTet4::NeumannDataTet4(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);

  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  quadPts = const_cast<PylithScalar*>(_quadPts);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDerivRef = const_cast<PylithScalar*>(_basisDerivRef);

  spatialDBFilename = const_cast<char*>(_spatialDBFilename);
  id = _id;
  label = const_cast<char*>(_label);

  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numVertices = _numVertices;
  numCells = _numCells;
  numCorners = _numCorners;
  cells = const_cast<int*>(_cells);

  tractionsCell = const_cast<PylithScalar*>(_tractionsCell);
  valsResidual = const_cast<PylithScalar*>(_valsResidual);

} // constructor

pylith::bc::NeumannDataTet4::~NeumannDataTet4(void)
{}


// End of file
