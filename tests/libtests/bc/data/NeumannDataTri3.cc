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

/* Mesh: tri3.mesh
 *
 * Neumann BC on edge defined by vertices 1, 3.
 * This edge has a length of sqrt(2).
 * The unit vectors should be:
 * Normal: ( 1/sqrt(2), -1/sqrt(2))
 * Shear:  ( 1/sqrt(2),  1/sqrt(2))
 * Applying shear and normal tractions of (1,  1) should give x and y
 * tractions of (sqrt(2), 0).
 * Integrating these tractions over the length of the edge and distributing
 * values to the edge nodes, we should get x and y forces of (1.0, 0.0)
 * at vertices 1 and 3.
 *
 */

#include "NeumannDataTri3.hh"

const char* pylith::bc::NeumannDataTri3::_meshFilename = 
  "data/tri3.mesh";

const int pylith::bc::NeumannDataTri3::_numBasis = 2;
const int pylith::bc::NeumannDataTri3::_numQuadPts = 1;
const PylithScalar pylith::bc::NeumannDataTri3::_quadPts[] = {
  0.0,
};
const PylithScalar pylith::bc::NeumannDataTri3::_quadWts[] = {
  2.0,
};
const PylithScalar pylith::bc::NeumannDataTri3::_basis[] = {
  0.5,
  0.5,
};
const PylithScalar pylith::bc::NeumannDataTri3::_basisDerivRef[] = {
 -0.5,
  0.5,
};

const char* pylith::bc::NeumannDataTri3::_spatialDBFilename =
  "data/tri3_tractions.spatialdb";
const int pylith::bc::NeumannDataTri3::_id = 0;
const char* pylith::bc::NeumannDataTri3::_label = "bc";

const int pylith::bc::NeumannDataTri3::_spaceDim = 2;
const int pylith::bc::NeumannDataTri3::_cellDim = 1;
const int pylith::bc::NeumannDataTri3::_numVertices = 2;
const int pylith::bc::NeumannDataTri3::_numCells = 1;
const int pylith::bc::NeumannDataTri3::_numCorners = 2;
/* Now vertices are renumbered in the submesh */
const int pylith::bc::NeumannDataTri3::_cells[] = {
  1 /*3*/, 2 /*5*/,
};

const PylithScalar pylith::bc::NeumannDataTri3::_tractionsCell[] = {
  1.4142135624,  0.0,
};
const PylithScalar pylith::bc::NeumannDataTri3::_valsResidual[] = {
  0.0,  0.0,
  1.0,  0.0,
  0.0,  0.0,
  1.0,  0.0,
};


pylith::bc::NeumannDataTri3::NeumannDataTri3(void)
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

pylith::bc::NeumannDataTri3::~NeumannDataTri3(void)
{}


// End of file
