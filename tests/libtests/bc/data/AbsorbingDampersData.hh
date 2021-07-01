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

#if !defined(pylith_bc_absorbingdampersdata_hh)
#define pylith_bc_absorbingdampersdata_hh

#include "pylith/utils/types.hh" // USES PylithScalar

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

  /// @name Scales information for nondimensionalization.
  //@{
  PylithScalar lengthScale; ///< Length scale.
  PylithScalar pressureScale; ///< Pressure scale.
  PylithScalar timeScale; ///< Time scale.
  PylithScalar densityScale; ///< Density scale.
  //@}

  /// @name Quadrature information
  //@{
  int numBasis; ///< Number of basis functions for cell
  int numQuadPts; ///< Number of quadrature points
  PylithScalar* quadPts; ///< Coordinates of quad pts in ref cell
  PylithScalar* quadWts; ///< Weights of quadrature points
  PylithScalar* basis; ///< Basis fns at quadrature points
  PylithScalar* basisDerivRef; ///< Derivatives of basis fns at quad pts
  //@}

  /// @name Parameter information
  //@{
  const char* spatialDBFilename; ///< Filename for database of parameters.
  int id; ///< Identifier for boundary condition
  char* label; ///< Label for boundary condition
  //@}

  /// @name Input fields
  //@{
  PylithScalar dt; ///< Time step
  PylithScalar* fieldTIncr; ///< Input increment field for time to to t+dt.
  PylithScalar* fieldT; ///< Input field at time t.
  PylithScalar* fieldTmdt; ///< Input field at time t-dt.
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
  PylithScalar* dampingConsts; ///< Expected values from initialization.
  PylithScalar* valsResidual; ///< Expected values from residual calculation.
  PylithScalar* valsJacobian; ///< Expected values from Jacobian calculation.
  //@}
};

#endif // pylith_bc_absorbingdampersdata_hh


// End of file
