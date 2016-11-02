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

#include "TestIsotropicLinearElasticityPlaneStrain_Data.hh"

extern "C" {
#include "pylith/fekernels/dispvel.h" // USES DispVel kernels
#include "pylith/fekernels/elasticity.h" // USES Elasticity kernels
#include "pylith/fekernels/linearelasticityplanestrain.h" // USES IsotropicLinearElasticityPlaneStrain kernels
}

// ----------------------------------------------------------------------
// Constructor
pylith::materials::TestIsotropicLinearElasticityPlaneStrain_Data::TestIsotropicLinearElasticityPlaneStrain_Data(const bool withInertia,
												      const bool withBodyForce,
												      const bool withReferenceState) :
  useInertia(withInertia),
  useBodyForce(withBodyForce),
  useReferenceState(withReferenceState)
{ // constructor
  dimension = 2;
  numSolnFields = 2;

  // RHS Residual (nfields*2)
  // disp
  _kernelsRHSResidual[0*2+0] = pylith_fekernels_DispVel_g0u;
  _kernelsRHSResidual[0*2+1] = NULL;
  // vel
  _kernelsRHSResidual[1*2+0] = NULL;
  _kernelsRHSResidual[1*2+1] = pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1v;

  // RHS Jacobian (nfields*nfields*4)
  // disp/disp
  _kernelsRHSJacobian[(0*2+0)*4+0] = NULL;
  _kernelsRHSJacobian[(0*2+0)*4+1] = NULL;
  _kernelsRHSJacobian[(0*2+0)*4+2] = NULL;
  _kernelsRHSJacobian[(0*2+0)*4+3] = NULL;

  // disp/vel
  _kernelsRHSJacobian[(0*2+1)*4+0] = pylith_fekernels_DispVel_Jg0uv;
  _kernelsRHSJacobian[(0*2+1)*4+1] = NULL;
  _kernelsRHSJacobian[(0*2+1)*4+2] = NULL;
  _kernelsRHSJacobian[(0*2+1)*4+3] = NULL;

  // vel/disp
  _kernelsRHSJacobian[(1*2+0)*4+0] = NULL;
  _kernelsRHSJacobian[(1*2+0)*4+1] = NULL;
  _kernelsRHSJacobian[(1*2+0)*4+2] = NULL;
  _kernelsRHSJacobian[(1*2+0)*4+3] = pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jg3vu;

  // vel/vel
  _kernelsRHSJacobian[(1*2+1)*4+0] = NULL;
  _kernelsRHSJacobian[(1*2+1)*4+1] = NULL;
  _kernelsRHSJacobian[(1*2+1)*4+2] = NULL;
  _kernelsRHSJacobian[(1*2+1)*4+3] = NULL;

  // LHS Residual (nfields*2)
  // disp
  _kernelsLHSResidual[0*2+0] = pylith_fekernels_DispVel_f0u;
  _kernelsLHSResidual[0*2+1] = NULL;
  // vel
  _kernelsLHSResidual[1*2+0] = NULL;
  _kernelsLHSResidual[1*2+1] = NULL;

  // LHS Jacobian Implicit (nfields*nfields*4)
  // disp/disp
  _kernelsLHSJacobianImplicit[(0*2+0)*4+0] = pylith_fekernels_DispVel_Jf0uu_implicit;
  _kernelsLHSJacobianImplicit[(0*2+0)*4+1] = NULL;
  _kernelsLHSJacobianImplicit[(0*2+0)*4+2] = NULL;
  _kernelsLHSJacobianImplicit[(0*2+0)*4+3] = NULL;

  // disp/vel
  _kernelsLHSJacobianImplicit[(0*2+1)*4+0] = NULL;
  _kernelsLHSJacobianImplicit[(0*2+1)*4+1] = NULL;
  _kernelsLHSJacobianImplicit[(0*2+1)*4+2] = NULL;
  _kernelsLHSJacobianImplicit[(0*2+1)*4+3] = NULL;

  // vel/disp
  _kernelsLHSJacobianImplicit[(1*2+0)*4+0] = NULL;
  _kernelsLHSJacobianImplicit[(1*2+0)*4+1] = NULL;
  _kernelsLHSJacobianImplicit[(1*2+0)*4+2] = NULL;
  _kernelsLHSJacobianImplicit[(1*2+0)*4+3] = NULL;

  // vel/vel
  _kernelsLHSJacobianImplicit[(1*2+1)*4+0] = NULL;
  _kernelsLHSJacobianImplicit[(1*2+1)*4+1] = NULL;
  _kernelsLHSJacobianImplicit[(1*2+1)*4+2] = NULL;
  _kernelsLHSJacobianImplicit[(1*2+1)*4+3] = NULL;

  // LHS Jacobian Explicit (nfields*nfields*4)
  // disp/disp
  _kernelsLHSJacobianExplicit[(0*2+0)*4+0] = pylith_fekernels_DispVel_Jf0uu_explicit;
  _kernelsLHSJacobianExplicit[(0*2+0)*4+1] = NULL;
  _kernelsLHSJacobianExplicit[(0*2+0)*4+2] = NULL;
  _kernelsLHSJacobianExplicit[(0*2+0)*4+3] = NULL;

  // disp/vel
  _kernelsLHSJacobianExplicit[(0*2+1)*4+0] = NULL;
  _kernelsLHSJacobianExplicit[(0*2+1)*4+1] = NULL;
  _kernelsLHSJacobianExplicit[(0*2+1)*4+2] = NULL;
  _kernelsLHSJacobianExplicit[(0*2+1)*4+3] = NULL;

  // vel/disp
  _kernelsLHSJacobianExplicit[(1*2+0)*4+0] = NULL;
  _kernelsLHSJacobianExplicit[(1*2+0)*4+1] = NULL;
  _kernelsLHSJacobianExplicit[(1*2+0)*4+2] = NULL;
  _kernelsLHSJacobianExplicit[(1*2+0)*4+3] = NULL;

  // vel/vel
  _kernelsLHSJacobianExplicit[(1*2+1)*4+0] = NULL;
  _kernelsLHSJacobianExplicit[(1*2+1)*4+1] = NULL;
  _kernelsLHSJacobianExplicit[(1*2+1)*4+2] = NULL;
  _kernelsLHSJacobianExplicit[(1*2+1)*4+3] = NULL;

  kernelsRHSResidual = const_cast<PetscPointFunc*>(_kernelsRHSResidual);
  kernelsRHSJacobian = const_cast<PetscPointJac*>(_kernelsRHSJacobian);
  kernelsLHSResidual = const_cast<PetscPointFunc*>( _kernelsLHSResidual);
  kernelsLHSJacobianImplicit = const_cast<PetscPointJac*>(_kernelsLHSJacobianImplicit);
  kernelsLHSJacobianExplicit = const_cast<PetscPointJac*>( _kernelsLHSJacobianExplicit);
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::materials::TestIsotropicLinearElasticityPlaneStrain_Data::~TestIsotropicLinearElasticityPlaneStrain_Data(void)
{ // destructor
} // destructor

// End of file
