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

#if !defined(pylith_materials_isotropiclinearelasticityplanestraindata_hh)
#define pylith_materials_isotropiclinearelasticityplanestraindata_hh

#include "pylith/utils/types.hh" // HASA PylithScalar
#include "pylith/topology/Field.hh" // HASA FieldBase::Discretization

#include "petscds.h" // USES PetscPointFunc, PetsPointJac

namespace pylith {
  namespace materials {
     class IsotropicLinearElasticityPlaneStrainData;
  } // pylith
} // materials

class pylith::materials::IsotropicLinearElasticityPlaneStrainData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  IsotropicLinearElasticityPlaneStrainData(void);

  /// Destructor
  ~IsotropicLinearElasticityPlaneStrainData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  // GENERAL, VALUES DEPEND ON TEST CASE
  char* meshFilename; ///< Name of file with ASCII mesh.
  char* materialLabel; ///< Label defining cells associated with material.
  int materialId; ///< Material id.
  char* boundaryLabel; ///< Group defining domain boundary.

  PylithReal lengthScale; ///< Length scale for nondimensionalization.
  PylithReal timeScale; ///< Time scale for nondimensionalization.
  PylithReal pressureScale; ///< Pressure scale for nondimensionalization.
  PylithReal densityScale; ///< Density scale for nondimensionalization.

  PylithReal t; ///< Time in simulation.
  PylithReal dt; ///< Time step in simulation.
  
  topology::Field::DiscretizeInfo* solnDiscretizations; ///< Discretizations for solution fields.
  char* soln1DBFilename; ///< Name of file with data for solution1.
  char* soln2DBFilename; ///< Name of file with data for solution2.

  int numAuxFields; ///< Number of auxiliary fields.
  char** auxFields; ///< Names of auxiliary fields.
  topology::Field::DiscretizeInfo* auxDiscretizations; ///< Discretizations for auxiliary fields.  
  char* auxDBFilename; ///< Name of file with data for auxFieldsDB.

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

  
  // SPECIFIC TO MATERIAL, VALUES DEPEND ON TEST CASE
  bool useInertia; ///< Test uses inertia.
  bool useBodyForce; ///< Test uses body force.

};

#endif // pylith_materials_materialdata_hh

// End of file
