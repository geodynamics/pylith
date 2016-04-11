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

#if !defined(pylith_materials_testisotropiclinearelasticityplanestrain_data_hh)
#define pylith_materials_testisotropiclinearelasticityplanestrain_data_hh

#include "TestMaterialNew_Data.hh" // ISA TestMaterialNew_Data

namespace pylith {
  namespace materials {
     class TestIsotropicLinearElasticityPlaneStrain_Data;
  } // pylith
} // materials

class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_Data : public TestMaterialNew_Data
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  TestIsotropicLinearElasticityPlaneStrain_Data(const bool useInertia,
						const bool useBodyForce,
						const bool useInitialState);

  /// Destructor
  ~TestIsotropicLinearElasticityPlaneStrain_Data(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  // SPECIFIC TO MATERIAL, VALUES DEPEND ON TEST CASE

  PetscPointFunc _kernelsRHSResidual[2*2]; ///< FE kernels for RHS residual, G(t,s).
  PetscPointJac _kernelsRHSJacobian[2*2*4]; ///< FE kernels for RHS Jacobian, G(t,s).
  PetscPointFunc _kernelsLHSResidual[2*2]; ///< FE kernels for LHS residual, F(t,s,\dot{s}).
  PetscPointJac _kernelsLHSJacobianImplicit[2*2*4]; ///< FE kernels for LHS Jacobian, F(t,s,\dot{s}) with implicit time-stepping.
  PetscPointJac _kernelsLHSJacobianExplicit[2*2*4];///< FE kernels for LHS Jacobian, F(t,s,\dot{s}) with

  bool useInertia; ///< Flag indicating test case uses inertia.
  bool useBodyForce; ///< Flag indicating test case uses body force.
  bool useInitialState; ///< Flag indicating test case uses initial state.

};

#endif // pylith_materials_testisotropiclinearelasticityplanestrain_data_hh

// End of file
