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

#if !defined(pylith_faults_cohesivekindata_hh)
#define pylith_faults_cohesivekindata_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
  namespace faults {
     class CohesiveKinData;
  } // pylith
} // faults

class pylith::faults::CohesiveKinData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  CohesiveKinData(void);

  /// Destructor
  ~CohesiveKinData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  char* meshFilename; ///< Filename for input mesh

  /// @name Scales information for nondimensionalization.
  //@{
  PylithScalar lengthScale; ///< Length scale.
  PylithScalar pressureScale; ///< Pressure scale.
  PylithScalar timeScale; ///< Time scale.
  PylithScalar densityScale; ///< Density scale.
  //@}

  /// @name Quadrature information
  //@{
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int cellDim; ///< Number of dimensions associated with cell
  int numBasis; ///< Number of vertices in cell
  int numQuadPts; ///< Number of quadrature points
  PylithScalar* quadPts; ///< Coordinates of quad pts in ref cell
  PylithScalar* quadWts; ///< Weights of quadrature points
  PylithScalar* basis; ///< Basis fns at quadrature points
  PylithScalar* basisDeriv; ///< Derivatives of basis fns at quad pts
  PylithScalar* verticesRef; ///< Coordinates of vertices in ref cell (dual basis)
  //@}

  /// @name Fault information
  //@{
  int id; ///< Fault material identifier
  char* label; ///< Label for fault
  char* edge; ///< Label for fault edge
  char* finalSlipFilename; ///< Name of db for final slip
  char* slipTimeFilename; ///< Name of db for slip time
  char* riseTimeFilename; ///< Name of db for rise time
  char* matPropsFilename; ///< Name of db for bulk material properties
  //@}

  /// @name Input fields
  //@{
  PylithScalar* fieldT; ///< Field over domain at time t.
  PylithScalar* fieldIncr; ///< Solution increment field over domain at time t.
  PylithScalar* jacobianLumped; ///< Lumped Jacobian.
  //@}

  /// @name Calculated values.
  //@{
  PylithScalar* orientation; ///< Expected values for fault orientation.
  PylithScalar* area; ///< Expected values for fault area.
  PylithScalar* residual; ///< Expected values from residual calculation.

  PylithScalar* jacobian; ///< Expected values from Jacobian calculation.

  /// Expected values for solution increment after adjustment.
  PylithScalar* fieldIncrAdjusted;

  int* verticesFault; ///< Expected points for Fault vertices
  int* edgesLagrange; ///< Expected points for Lagrange multipliers
  int* verticesPositive; ///< Expected points for vertices on + side of fault.
  int* verticesNegative; ///< Expected points for vertices on - side of fault.
  int numFaultVertices; ///< Number of fault vertices

  int numCohesiveCells;
  int* cellMappingFault;
  int* cellMappingCohesive;
  //@}

};

#endif // pylith_faults_cohesivekindata_hh

// End of file
