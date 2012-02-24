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

#if !defined(pylith_faults_cohesivekindata_hh)
#define pylith_faults_cohesivekindata_hh

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

  /// @name Quadrature information
  //@{
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int cellDim; ///< Number of dimensions associated with cell
  int numBasis; ///< Number of vertices in cell
  int numQuadPts; ///< Number of quadrature points
  double* quadPts; ///< Coordinates of quad pts in ref cell
  double* quadWts; ///< Weights of quadrature points
  double* basis; ///< Basis fns at quadrature points
  double* basisDeriv; ///< Derivatives of basis fns at quad pts
  double* verticesRef; ///< Coordinates of vertices in ref cell (dual basis)
  //@}

  /// @name Fault information
  //@{
  int id; ///< Fault material identifier
  char* label; ///< Label for fault
  char* finalSlipFilename; ///< Name of db for final slip
  char* slipTimeFilename; ///< Name of db for slip time
  char* riseTimeFilename; ///< Name of db for rise time
  char* matPropsFilename; ///< Name of db for bulk material properties
  //@}

  /// @name Input fields
  //@{
  double* fieldT; ///< Field over domain at time t.
  double* fieldIncr; ///< Solution increment field over domain at time t.
  double* jacobianLumped; ///< Lumped Jacobian.
  //@}

  /// @name Calculated values.
  //@{
  double* orientation; ///< Expected values for fault orientation.
  double* area; ///< Expected values for fault area.
  double* residual; ///< Expected values from residual calculation.

  /// Expected values from residual calculation using solution increment.
  double* residualIncr;

  double* jacobian; ///< Expected values from Jacobian calculation.

  /// Expected values for solution increment after adjustment.
  double* fieldIncrAdjusted;

  int* verticesFault; ///< Expected points for Fault vertices
  int* verticesLagrange; ///< Expected points for Lagrange vertices
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
