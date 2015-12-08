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

#if !defined(pylith_materials_isotropiclinearelasticityplanestraindata_tri3_hh)
#define pylith_materials_isotropiclinearelasticityplanestraindata_tri3_hh

#include "IsotropicLinearElasticityPlaneStrainData.hh"

namespace pylith {
  namespace materials {
    class IsotropicLinearElasticityPlaneStrainData_Tri3;
  } // pylith
} // materials

class pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3 : public IsotropicLinearElasticityPlaneStrainData
{

public: 

  /// Constructor
  IsotropicLinearElasticityPlaneStrainData_Tri3(void);

  /// Destructor
  ~IsotropicLinearElasticityPlaneStrainData_Tri3(void);

  /** Function for computing displacement solution to match LHS and RHS residual functions.
   *
   * @param[in] dim Spatial dimension.
   * @param[in] t Current time.
   * @param[in] x Coordinates (nondimensioned) of point location for query.
   * @param[in] nvalues Size of values array.
   * @param[out] values Array of values to be returned.
   * @param[in] context Query context.
   * @returns PETSc error code (0 for success).
   */
  static
  PetscErrorCode _querySolutionDisplacement(PylithInt dim, 
					    PylithReal t, 
					    const PylithReal x[],
					    PylithInt nvalues,
					    PylithScalar* values, 
					    void* context);

  /** Function for computing velocity solution to match LHS and RHS residual functions.
   *
   * @param[in] dim Spatial dimension.
   * @param[in] t Current time.
   * @param[in] x Coordinates (nondimensioned) of point location for query.
   * @param[in] nvalues Size of values array.
   * @param[out] values Array of values to be returned.
   * @param[in] context Query context.
   * @returns PETSc error code (0 for success).
   */
  static
  PetscErrorCode _querySolutionVelocity(PylithInt dim, 
					PylithReal t, 
					const PylithReal x[],
					PylithInt nvalues,
					PylithScalar* values, 
					void* context);

private:

  static const char* _filenameMesh;
  static const char* _label;
  static const int _id;
  static const int _dimension;

  static const bool _useInertia;
  static const bool _useBodyForce;

  static const int _numSolnFields;
  static const pylith::topology::Field::DiscretizeInfo _discretizations[1];

  static const char* _filenameAuxFieldsDB;

  static const PetscPointFunc _kernelsRHSResidual[2*2];
  static const PetscPointJac _kernelsRHSJacobian[2*2*4];
  static const PetscPointFunc _kernelsLHSResidual[2*2];
  static const PetscPointJac _kernelsLHSJacobianImplicit[2*2*4];
  static const PetscPointJac _kernelsLHSJacobianExplicit[2*2*4];

  static const PylithReal _lengthScale;
  static const PylithReal _timeScale;
  static const PylithReal _densityScale;
  static const PylithReal _pressureScale;

};

#endif // pylith_materials_isotropiclinearelasticityplanestraindata_tri3_hh

// End of file
