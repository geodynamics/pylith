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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#if !defined(pylith_feassemble_quadraturedata_hh)
#define pylith_feassemble_quadraturedata_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
  namespace feassemble {
     class QuadratureData;
  } // pylith
} // feassemble

class pylith::feassemble::QuadratureData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  QuadratureData(void);

  /// Destructor
  ~QuadratureData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  int numVertices; ///< Number of vertices
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int numCells; ///< Number of cells (=1)
  int cellDim; ///< Number of dimensions associated with cell
  int numBasis; ///< Number of vertices in cell
  int numQuadPts; ///< Number of quadrature points

  PylithScalar* vertices; ///< Pointer to coordinates of vertices
  int* cells; ///< Pointer to indices of vertices in cells

  PylithScalar* verticesRef; ///< Coordinates of vertices in ref cell
  PylithScalar* quadPtsRef; ///< Coordinates of quad pts in ref cell
  PylithScalar* quadWts; ///< Weights of quadrature points
  PylithScalar* quadPts; ///< Coordinates of quad pts in cell

  PylithScalar* basis; ///< Basis fns at quadrature points
  PylithScalar* basisDerivRef; ///< Derivatices of basis fns at quad pts (cell)
  PylithScalar* basisDeriv; ///< Derivatices of basis fns at quad pts (global)
  PylithScalar* jacobian; ///< Jacobian at quadrature points
  PylithScalar* jacobianDet; ///< Determinant of quadrature points
  PylithScalar* jacobianInv; ///< Inverse of Jacobian at quadruature points

};

#endif // pylith_feassemble_quadraturedata_hh

// End of file
