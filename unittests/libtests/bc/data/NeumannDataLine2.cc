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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/* Mesh: line2.mesh
 *
 * Neumann BC on vertex 0.
 * For this problem, the traction should be the same as a point force.
 * For a point force of 1 applied at vertex 0, the applied traction
 * and the residual should both have a value of 1 for this vertex.
 *
 */

#include "NeumannDataLine2.hh"

const char* pylith::bc::NeumannDataLine2::_meshFilename = 
  "data/line2.mesh";

const int pylith::bc::NeumannDataLine2::_numBasis = 1;
const int pylith::bc::NeumannDataLine2::_numQuadPts = 1;
const PylithScalar pylith::bc::NeumannDataLine2::_quadPts[] = {
  0.0,
};
const PylithScalar pylith::bc::NeumannDataLine2::_quadWts[] = {
  1.0,
};
const PylithScalar pylith::bc::NeumannDataLine2::_basis[] = {
  1.0,
};
const PylithScalar pylith::bc::NeumannDataLine2::_basisDerivRef[] = {
  1.0,
};

const char* pylith::bc::NeumannDataLine2::_spatialDBFilename =
  "data/line2_tractions.spatialdb";
const int pylith::bc::NeumannDataLine2::_id = 0;
const char* pylith::bc::NeumannDataLine2::_label = "bc1";

const int pylith::bc::NeumannDataLine2::_spaceDim = 1;
const int pylith::bc::NeumannDataLine2::_cellDim = 0;
const int pylith::bc::NeumannDataLine2::_numVertices = 1;
const int pylith::bc::NeumannDataLine2::_numCells = 1;
const int pylith::bc::NeumannDataLine2::_numCorners = 1;
/* Now vertices are renumbered in the submesh */
const int pylith::bc::NeumannDataLine2::_cells[] = {
  1 /*2*/,
};

const PylithScalar pylith::bc::NeumannDataLine2::_tractionsCell[] = {
  1.0,
};
const PylithScalar pylith::bc::NeumannDataLine2::_valsResidual[] = {
  1.0,
  0.0,
  0.0,
};


pylith::bc::NeumannDataLine2::NeumannDataLine2(void)
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

pylith::bc::NeumannDataLine2::~NeumannDataLine2(void)
{}


// End of file
