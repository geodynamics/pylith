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

#if !defined(pylith_bc_absorbingdampersdata_hh)
#define pylith_bc_absorbingdampersdata_hh

namespace pylith {
  namespace bc {
     class AbsorbingDampersData;
  } // pylith
} // bc

class pylith::bc::AbsorbingDampersData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  AbsorbingDampersData(void);

  /// Destructor
  ~AbsorbingDampersData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  char* meshFilename; ///< Name of file with input mesh

  /// @name Quadrature information
  //@{
  int numBasis; ///< Number of basis functions for cell
  int numQuadPts; ///< Number of quadrature points
  double* quadPts; ///< Coordinates of quad pts in ref cell
  double* quadWts; ///< Weights of quadrature points
  double* basis; ///< Basis fns at quadrature points
  double* basisDerivRef; ///< Derivatives of basis fns at quad pts
  //@}

  /// @name Parameter information
  //@{
  const char* spatialDBFilename; ///< Filename for database of parameters.
  int id; ///< Identifier for boundary condition
  char* label; ///< Label for boundary condition
  //@}

  /// @name Input fields
  //@{
  double dt; ///< Time step
  double* fieldTIncr; ///< Input increment field for time to to t+dt.
  double* fieldT; ///< Input field at time t.
  double* fieldTmdt; ///< Input field at time t-dt.
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
  double* dampingConsts; ///< Expected values from initialization.
  double* valsResidual; ///< Expected values from residual calculation.
  double* valsJacobian; ///< Expected values from Jacobian calculation.
  //@}
};

#endif // pylith_bc_absorbingdampersdata_hh


// End of file
