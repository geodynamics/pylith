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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#if !defined(pylith_bc_neumanndata_hh)
#define pylith_bc_neumanndata_hh

namespace pylith {
  namespace bc {
     class NeumannData;
  } // pylith
} // bc

class pylith::bc::NeumannData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  NeumannData(void);

  /// Destructor
  ~NeumannData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  char* meshFilename; ///< Filename for input mesh.

  /// @name Quadrature information
  //@{
  int numBasis; ///< Number of basis functions for surface cells.
  int numQuadPts; ///< Number of quadrature points per boundary cell.
  double* quadPts; ///< Coordinates of quadrature points in ref cell.
  double* quadWts; ///< Weights of quadrature points.
  double* basis; ///< Cell basis functions at quad points.
  double* basisDerivRef; ///< Derivatives of basis functions at quad points.
  //@}

  /// @name Parameter information
  //@{
  char* spatialDBFilename; ///< Filename of simple spatial database.
  int id; ///< Boundary condition identifier
  char* label; ///< Label for boundary condition group
  //@}

  /// @name Boundary mesh information
  //@{
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int cellDim; ///< Number of dimensions associated with cell
  int numVertices; ///< Number of vertices
  int numCells; ///< Number of cells
  int numCorners; ///< Number of vertices in cell
  int* cells; ///< Indices of vertices in cells
  //@}

  /// @name Calculated values.
  //@{
  double* tractionsCell; ///< Expected traction values at quadrature points.
  double* valsResidual; ///< Expected residual at each vertex.
  //@}


};

#endif // pylith_bc_neumanndata_hh

// End of file
