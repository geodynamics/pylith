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

#if !defined(pylith_faults_cohesivedynldata_hh)
#define pylith_faults_cohesivedynldata_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
  namespace faults {
     class CohesiveDynData;
  } // pylith
} // faults

class pylith::faults::CohesiveDynData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  CohesiveDynData(void);

  /// Destructor
  ~CohesiveDynData(void);

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
  char* initialTractFilename; ///< Name of db for initial tractions.
  //@}

  /// @name Input fields
  //@{
  PylithScalar* fieldT; ///< Solution field at time t.
  PylithScalar* fieldIncrStick; ///< Soln increment field at time t for stick case.
  PylithScalar* fieldIncrSlip; ///< Soln increment field at time t for slipping case.
  PylithScalar* fieldIncrOpen; ///< Soln increment field at time t for opening case.
  PylithScalar* jacobian; ///< Jacobian sparse matrix.
  //@}

  /// @name Calculated values.
  //@{
  PylithScalar* orientation; ///< Expected values for fault orientation.
  PylithScalar* area; ///< Expected values for fault area.
  PylithScalar* initialTractions; ///< Expected values for initial tractions.
  PylithScalar* slipStickE; ///< Expected values for slip for sticking case.
  PylithScalar* fieldIncrSlipE; ///< Expected values for solution increment for slipping case.
  PylithScalar* slipSlipE; ///< Expected values for slip for slipping case.
  PylithScalar* fieldIncrOpenE; ///< Expected values for solution increment for opening case.
  PylithScalar* slipOpenE; ///< Expected values for slip for opening case.

  int* constraintEdges; ///< Expected points for constraint edges
  int* negativeVertices; ///< Expected points for negative side fault vertices
  int numConstraintEdges; ///< Number of constraint edges
  //@}

};

#endif // pylith_faults_cohesivedynldata_hh

// End of file
