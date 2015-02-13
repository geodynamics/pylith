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

#if !defined(pylith_feassemble_integratordata_hh)
#define pylith_feassemble_integratordata_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
  namespace feassemble {
     class IntegratorData;
  } // pylith
} // feassemble

class pylith::feassemble::IntegratorData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  IntegratorData(void);

  /// Destructor
  ~IntegratorData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  /// @name Mesh information
  //@{
  int spaceDim; ///< Number of dimensions in vertex coordinates
  PylithScalar* gravityVec; ///< Constant gravity vector (for now)
  int cellDim; ///< Number of dimensions associated with cell
  int numVertices; ///< Number of vertices
  int numCells; ///< Number of cells
  PylithScalar* vertices; ///< Coordinates of vertices
  int* cells; ///< Indices of vertices in cells
  PylithScalar* verticesRef; ///< Coordinates of vertices in ref cell (dual basis)
  //@}

  /// @name Quadrature information
  //@{
  int numBasis; ///< Number of vertices in cell
  int numQuadPts; ///< Number of quadrature points
  PylithScalar* quadPts; ///< Coordinates of quad pts in ref cell
  PylithScalar* quadWts; ///< Weights of quadrature points
  PylithScalar* basis; ///< Basis fns at quadrature points
  PylithScalar* basisDerivRef; ///< Derivatives of basis fns at quad pts
  //@}

  /// @name Material information
  //@{
  char* matType; ///< String corresponding to material type.
  char* matDBFilename; ///< Filename for database of material properties.
  int matId; ///< Material identifier.
  char* matLabel; ///< Label of material.
  //@}

  /// @name Scales information for nondimensionalization.
  //@{
  PylithScalar lengthScale; ///< Length scale.
  PylithScalar pressureScale; ///< Pressure scale.
  PylithScalar timeScale; ///< Time scale.
  PylithScalar densityScale; ///< Density scale.
  //@}  

  /// @name Input fields
  //@{
  PylithScalar dt; ///< Time step
  PylithScalar* fieldTIncr; ///< Input field increment for time t to time t+dt.
  PylithScalar* fieldT; ///< Input field at time t.
  PylithScalar* fieldTmdt; ///< Input field at time t-dt.
  //@}

  /// @name Calculated values.
  //@{
  PylithScalar* valsResidual; ///< Expected values from residual calculation.
  PylithScalar* valsJacobian; ///< Expected values from Jacobian calculation.
  //@}
};

#endif // pylith_feassemble_integratordata_hh

// End of file
