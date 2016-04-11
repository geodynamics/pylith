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

#if !defined(pylith_materials_testmaterialnew_data_hh)
#define pylith_materials_testmaterialnew_data_hh

#include "pylith/utils/types.hh" // HASA PylithScalar
#include "pylith/topology/Field.hh" // HASA FieldBase::Discretization

#include "petscds.h" // USES PetscPointFunc, PetsPointJac

namespace pylith {
  namespace materials {
     class TestMaterialNew_Data;
  } // pylith
} // materials

class pylith::materials::TestMaterialNew_Data
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  TestMaterialNew_Data(void);

  /// Destructor
  ~TestMaterialNew_Data(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  // GENERAL, VALUES DEPEND ON TEST CASE
  const char* meshFilename; ///< Name of file with ASCII mesh.
  const char* materialLabel; ///< Label defining cells associated with material.
  int materialId; ///< Material id.
  const char* boundaryLabel; ///< Group defining domain boundary.

  PylithReal lengthScale; ///< Length scale for nondimensionalization.
  PylithReal timeScale; ///< Time scale for nondimensionalization.
  PylithReal pressureScale; ///< Pressure scale for nondimensionalization.
  PylithReal densityScale; ///< Density scale for nondimensionalization.

  PylithReal t; ///< Time for solution in simulation.
  PylithReal dt; ///< Time step in simulation.
  PylithReal tshift; ///< Time shift for LHS Jacobian.
  
  topology::Field::DiscretizeInfo* solnDiscretizations; ///< Discretizations for solution fields.
  const char* solnDBFilename; ///< Name of file with data for solution.
  const char* pertDBFilename; ///< Name of file with data for perturbation.

  int numAuxFields; ///< Number of auxiliary fields.
  const char** auxFields; ///< Names of auxiliary fields.
  topology::Field::DiscretizeInfo* auxDiscretizations; ///< Discretizations for auxiliary fields.  
  const char* auxDBFilename; ///< Name of file with data for auxFieldsDB.

  // GENERAL, VALUES DEPEND ONLY ON MATERIAL
  int dimension; ///< Dimension of material.
  int numSolnFields; ///< Number of solution fields.

  static const int numKernelsResidual;
  static const int numKernelsJacobian;
  PetscPointFunc* kernelsRHSResidual; ///< FE kernels for RHS residual, G(t,s).
  PetscPointJac* kernelsRHSJacobian; ///< FE kernels for RHS Jacobian, G(t,s).
  PetscPointFunc* kernelsLHSResidual; ///< FE kernels for LHS residual, F(t,s,\dot{s}).
  PetscPointJac* kernelsLHSJacobianImplicit; ///< FE kernels for LHS Jacobian, F(t,s,\dot{s}) with implicit time-stepping.
  PetscPointJac* kernelsLHSJacobianExplicit;///< FE kernels for LHS Jacobian, F(t,s,\dot{s}) with expicit time-stepping.
  
};

#endif // pylith_materials_testmaterialnew_data_hh

// End of file
