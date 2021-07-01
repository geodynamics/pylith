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

/* Mesh: quad4.mesh
 *
 * Neumann BC on edge defined by vertices 0, 2, 4.
 * A linear variation in the normal traction is specified, ranging from
 * a value of zero at x = -1 (vertex 0) to one at x = 1 (vertex 4).
 * The traction values at the integration points should be:
 * cell 0:  -0.1056624327, -0.3943375673
 * cell 1:  -0.6056624327, -0.8943375673
 * The Jacobian determinant along the element edges is 1/4.
 * The shape function values at the two element integration points are:
 * 0.21132486541, 0.78867513459
 * Integrating over the cell edges, the nodal forces (y-direction)
 * should be:
 * vertex 0:  -0.08333333333
 * vertex 2:  -0.5
 * vertex 4:  -0.41666666667
 *
 */

#include "NeumannDataQuad4.hh"

const char* pylith::bc::NeumannDataQuad4::_meshFilename = 
  "data/quad4.mesh";

const int pylith::bc::NeumannDataQuad4::_numBasis = 2;
const int pylith::bc::NeumannDataQuad4::_numQuadPts = 2;
const PylithScalar pylith::bc::NeumannDataQuad4::_quadPts[] = {
  -0.57735027,
   0.57735027,
};
const PylithScalar pylith::bc::NeumannDataQuad4::_quadWts[] = {
  1.0,
  1.0,
};
const PylithScalar pylith::bc::NeumannDataQuad4::_basis[] = {
  0.78867513459,
  0.21132486541,
  0.21132486541,
  0.78867513459,
};
const PylithScalar pylith::bc::NeumannDataQuad4::_basisDerivRef[] = {
 -0.5,
  0.5,
 -0.5,
  0.5,
};

const char* pylith::bc::NeumannDataQuad4::_spatialDBFilename =
  "data/quad4_tractions.spatialdb";
const int pylith::bc::NeumannDataQuad4::_id = 0;
const char* pylith::bc::NeumannDataQuad4::_label = "bc3";

const int pylith::bc::NeumannDataQuad4::_spaceDim = 2;
const int pylith::bc::NeumannDataQuad4::_cellDim = 1;
const int pylith::bc::NeumannDataQuad4::_numVertices = 3;
const int pylith::bc::NeumannDataQuad4::_numCells = 2;
const int pylith::bc::NeumannDataQuad4::_numCorners = 2;
/* Now vertices are renumbered in the submesh */
const int pylith::bc::NeumannDataQuad4::_cells[] = {
  2 /*2*/, 3 /*4*/,
  3 /*4*/, 4 /*6*/,
};

const PylithScalar pylith::bc::NeumannDataQuad4::_tractionsCell[] = {
  0.0, -0.1056624327,
  0.0, -0.3943375673,
  0.0, -0.6056624327,
  0.0, -0.8943375673,
};
const PylithScalar pylith::bc::NeumannDataQuad4::_valsResidual[] = {
  0.0, -0.08333333333,
  0.0,  0.0,
  0.0, -0.5,
  0.0,  0.0,
  0.0, -0.41666666667,
  0.0,  0.0,
};


pylith::bc::NeumannDataQuad4::NeumannDataQuad4(void)
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

pylith::bc::NeumannDataQuad4::~NeumannDataQuad4(void)
{}


// End of file
