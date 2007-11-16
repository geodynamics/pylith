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
const double pylith::bc::NeumannDataTri3::_quadPts[] = {
  0.0,
};
const double pylith::bc::NeumannDataTri3::_quadWts[] = {
  2.0,
};
const double pylith::bc::NeumannDataTri3::_basis[] = {
  0.5,
  0.5,
};
const double pylith::bc::NeumannDataTri3::_basisDerivRef[] = {
 -0.5,
  0.5,
};

const char* pylith::bc::NeumannDataTri3::_spatialDBFilename =
  "data/tri3-tractions.spatialdb";
const int pylith::bc::NeumannDataTri3::_id = 0;
const char* pylith::bc::NeumannDataTri3::_label = "bc";

const int pylith::bc::NeumannDataTri3::_spaceDim = 2;
const int pylith::bc::NeumannDataTri3::_cellDim = 1;

const int pylith::bc::NeumannDataTri3::_numBoundaryVertices = 2;
const int pylith::bc::NeumannDataTri3::_numBoundaryCells = 1;
const int pylith::bc::NeumannDataTri3::_numCorners = 2;
const double pylith::bc::NeumannDataTri3::_cellVertices[] = { 0.0,-1.0,
							      1.0, 0.0};
const double pylith::bc::NeumannDataTri3::_tractionsCell[] = {
  1.4142135624,  0.0,
};
const double pylith::bc::NeumannDataTri3::_valsResidual[] = {
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
  quadPts = const_cast<double*>(_quadPts);
  quadWts = const_cast<double*>(_quadWts);
  basis = const_cast<double*>(_basis);
  basisDerivRef = const_cast<double*>(_basisDerivRef);

  spatialDBFilename = const_cast<char*>(_spatialDBFilename);
  id = _id;
  label = const_cast<char*>(_label);

  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numBoundaryVertices = _numBoundaryVertices;
  numBoundaryCells = _numBoundaryCells;
  numCorners = _numCorners;

  cellVertices = const_cast<double*>(_cellVertices);
  tractionsCell = const_cast<double*>(_tractionsCell);
  valsResidual = const_cast<double*>(_valsResidual);

} // constructor

pylith::bc::NeumannDataTri3::~NeumannDataTri3(void)
{}


// End of file
