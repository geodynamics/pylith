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

#if !defined(pylith_feassemble_integratordata_hh)
#define pylith_feassemble_integratordata_hh

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
  char* meshFilename; ///< Name of mesh file.
  //@}

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
  //@}

  /// @name Material information
  //@{
  double* matType; ///< String corresponding to material type.
  double* matDBFilename; ///< Filename for database of material properties.
  double* matId; ///< Material identifier.
  double* matLabel; ///< Label of material.
  //@}

  /// @name Input fields
  //@{
  double* fieldTpdt; ///< Input field at time t+dt.
  double* fieldT; ///< Input field at time t.
  double* fieldTmdt; ///< Input field at time t-dt.
  //@}

  /// @name Calculated values.
  //@{
  double* valsResidual; ///< Expected values from residual calculation.
  double* valsJacobian; ///< Expected values from Jacobian calculation.
  //@}
};

#endif // pylith_feassemble_integratordata_hh

// End of file
