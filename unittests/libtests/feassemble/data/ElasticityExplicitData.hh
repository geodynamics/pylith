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

#if !defined(pylith_feassemble_elasticityexplicitdata_hh)
#define pylith_feassemble_elasticityexplicitdata_hh

#include "IntegratorData.hh" // ISA IntegratorData

namespace pylith {
  namespace feassemble {
     class ElasticityExplicitData;
  } // pylith
} // feassemble

class pylith::feassemble::ElasticityExplicitData : public IntegratorData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  ElasticityExplicitData(void);

  /// Destructor
  ~ElasticityExplicitData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  /// @name Calculated values.
  //@{
  double* valsResidualLumped; ///< Expected values from residual calculation (lumped Jacobian).
  double* valsJacobianLumped; ///< Expected values from lumped Jacobian calculation.
  //@}
};

#endif // pylith_feassemble_elasticityexplicitdata_hh

// End of file
