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

#include "IsotropicLinearElasticityPlaneStrainData_Tri3.hh"

extern "C" {
#include "pylith/fekernels/dispvel.h" // USES DispVel kernels
#include "pylith/fekernels/elasticity.h" // USES Elasticity kernels
#include "pylith/fekernels/linearelasticityplanestrain.h" // USES IsotropicLinearElasticityPlaneStrain kernels
}

const char* pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_filenameMesh = "data/tri3.mesh";
const char* pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_label = "IsotropicLinearElascitity";
const int pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_id = 0;
const int pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_dimension = 2;

const bool pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_useInertia = false;
const bool pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_useBodyForce = false;

const int pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_numSolnFields = 1;
const pylith::topology::Field::DiscretizeInfo pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_discretizations[1] = {
  {1, 1, true},
};

const PetscPointFunc pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_kernelsRHSResidual[2*2] =
  {
    // displacement
    NULL,
    pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1,
    // velocity
    pylith_fekernels_DispVel_g0,
    NULL,
  };

const PetscPointJac pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_kernelsRHSJacobian[2*2*4] = 
  {
    // disp/disp
    NULL,
    NULL,
    NULL,
    pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jg3_uu,
    // disp/vel
    NULL,
    NULL,
    NULL,
    NULL,
    // vel/disp
    NULL,
    NULL,
    NULL,
    NULL,
    // vel/vel
    pylith_fekernels_DispVel_Jg0_vv,
    NULL,
    NULL,
    NULL,
};

const PetscPointFunc pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_kernelsLHSResidual[2*2] = 
  {
    // displacement
    NULL,
    NULL,
    // velocity
    pylith_fekernels_DispVel_f0,
    NULL,
  };

const PetscPointJac pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_kernelsLHSJacobianImplicit[2*2*4] = 
  {
    // disp/disp
    NULL,
    NULL,
    NULL,
    NULL,
    // disp/vel
    NULL,
    NULL,
    NULL,
    NULL,
    // vel/disp
    pylith_fekernels_DispVel_Jf0_vu_implicit,
    NULL,
    NULL,
    NULL,
    // vel/vel
    NULL,
    NULL,
    NULL,
    NULL,
  };
const PetscPointJac pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_kernelsLHSJacobianExplicit[2*2*4] = 
  {
    // disp/disp
    NULL,
    NULL,
    NULL,
    NULL,
    // disp/vel
    NULL,
    NULL,
    NULL,
    NULL,
    // vel/disp
    pylith_fekernels_DispVel_Jf0_vu_explicit,
    NULL,
    NULL,
    NULL,
    // vel/vel
    NULL,
    NULL,
    NULL,
    NULL,
  };

const char* pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_filenameAuxFieldsDB = "data/matinitialize.spatialdb";

const PylithReal pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_lengthScale =   1.00000000e+03;
const PylithReal pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_timeScale =   2.00000000e+00;
const PylithReal pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_densityScale =   9.00000000e+04;
const PylithReal pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_pressureScale =   2.25000000e+10;


pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::IsotropicLinearElasticityPlaneStrainData_Tri3(void)
{ // constructor
  filenameMesh = const_cast<char*>(_filenameMesh);
  label = const_cast<char*>(_label);
  id = _id;
  dimension = _dimension;

  useInertia = _useInertia;
  useBodyForce = _useBodyForce;

  numSolnFields = _numSolnFields;
  discretizations = const_cast<pylith::topology::Field::DiscretizeInfo*>(_discretizations);

  filenameAuxFieldsDB = const_cast<char*>(_filenameAuxFieldsDB);

  kernelsRHSResidual = const_cast<PetscPointFunc*>(_kernelsRHSResidual);
  kernelsRHSJacobian = const_cast<PetscPointJac*>(_kernelsRHSJacobian);
  kernelsLHSResidual = const_cast<PetscPointFunc*>( _kernelsLHSResidual);
  kernelsLHSJacobianImplicit = const_cast<PetscPointJac*>(_kernelsLHSJacobianImplicit);
  kernelsLHSJacobianExplicit = const_cast<PetscPointJac*>( _kernelsLHSJacobianExplicit);

  querySolutionDisplacement = _querySolutionDisplacement;
  querySolutionDisplacementDot = _querySolutionDisplacementDot;
  querySolutionVelocity = _querySolutionVelocity;

  lengthScale = _lengthScale;
  timeScale = _timeScale;
  densityScale = _densityScale;
  pressureScale = _pressureScale;
} // constructor

pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::~IsotropicLinearElasticityPlaneStrainData_Tri3(void)
{}

// Function for computing displacement solution to match LHS and RHS residual functions.
PetscErrorCode
pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_querySolutionDisplacement(PylithInt dim, 
											     PylithReal t, 
											     const PylithReal x[],
											     PylithInt nvalues,
											     PylithScalar* values, 
											     void* context)
{ // _querySolutionDisplacement
  PYLITH_METHOD_BEGIN;

  const int _dim = 2;
  const int _nvalues = _dim;

  assert(x);
  assert(values);
  assert(context);
  assert(_nvalues == nvalues);
  assert(_dim == dim);

  const pylith::topology::FieldQuery::DBQueryContext* queryctx = (pylith::topology::FieldQuery::DBQueryContext*)context;assert(queryctx);

  // Tell database which values we want.
  const int numDBValues = _dim;
  PylithReal dbValues[numDBValues];
  const char* dbValueNames[numDBValues] = {"displacement-x", "displacement-y"};
  queryctx->db->queryVals(dbValueNames, numDBValues);

  // Dimensionalize query location coordinates.
  assert(queryctx->lengthScale > 0);
  double xDim[_dim];
  for (int i=0; i < _dim; ++i) {
    xDim[i] = x[i] * queryctx->lengthScale;
  } // for

  assert(queryctx->cs);
  const int err = queryctx->db->query(dbValues, numDBValues, xDim, _dim, queryctx->cs);
  if (err) {
    std::ostringstream msg;
    msg << "Could not find displacement field at (";
    for (int i=0; i < _dim; ++i)
      msg << "  " << xDim[i];
    msg << ") in using spatial database '" << queryctx->db->label() << "'.";
    (*PetscErrorPrintf)(msg.str().c_str());
    PYLITH_METHOD_RETURN(1);
  } // if
  
  for (int i=0; i < _dim; ++i) {
    values[i] = dbValues[i] / queryctx->valueScale;
  } // for

  PYLITH_METHOD_RETURN(0);
} // _querySolutionDisplacement


// Function for computing velocity solution to match LHS and RHS residual functions.
PetscErrorCode
pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_querySolutionVelocity(PylithInt dim, 
											 PylithReal t, 
											 const PylithReal x[],
											 PylithInt nvalues,
											 PylithScalar* values, 
											 void* context)
{ // _querySolutionVelocity
  PYLITH_METHOD_BEGIN;

  const int _dim = 2;
  const int _nvalues = _dim;

  assert(x);
  assert(values);
  assert(context);
  assert(_nvalues == nvalues);
  assert(_dim == dim);

  const pylith::topology::FieldQuery::DBQueryContext* queryctx = (pylith::topology::FieldQuery::DBQueryContext*)context;assert(queryctx);

  // Tell database which values we want.
  const int numDBValues = _dim;
  PylithReal dbValues[numDBValues];
  const char* dbValueNames[numDBValues] = {"velocity-x", "velocity-y"};
  queryctx->db->queryVals(dbValueNames, numDBValues);

  // Dimensionalize query location coordinates.
  assert(queryctx->lengthScale > 0);
  double xDim[_dim];
  for (int i=0; i < _dim; ++i) {
    xDim[i] = x[i] * queryctx->lengthScale;
  } // for

  assert(queryctx->cs);
  const int err = queryctx->db->query(dbValues, numDBValues, xDim, _dim, queryctx->cs);
  if (err) {
    std::ostringstream msg;
    msg << "Could not find velocity field at (";
    for (int i=0; i < _dim; ++i)
      msg << "  " << xDim[i];
    msg << ") in using spatial database '" << queryctx->db->label() << "'.";
    (*PetscErrorPrintf)(msg.str().c_str());
    PYLITH_METHOD_RETURN(1);
  } // if
  
  for (int i=0; i < _dim; ++i) {
    values[i] = dbValues[i] / queryctx->valueScale;
  } // for

  PYLITH_METHOD_RETURN(0);
} // _querySolutionVelocity


// Function for computing displacement solution to match LHS and RHS residual functions.
PetscErrorCode
pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_querySolutionDisplacementDot(PylithInt dim, 
												PylithReal t, 
												const PylithReal x[],
												PylithInt nvalues,
												PylithScalar* values, 
												void* context)
{ // _querySolutionDisplacementDot
  PYLITH_METHOD_BEGIN;

  const int _dim = 2;
  const int _nvalues = _dim;

  assert(x);
  assert(values);
  assert(context);
  assert(_nvalues == nvalues);
  assert(_dim == dim);

  const pylith::topology::FieldQuery::DBQueryContext* queryctx = (pylith::topology::FieldQuery::DBQueryContext*)context;assert(queryctx);

  // Tell database which values we want.
  const int numDBValues = _dim;
  PylithReal dbValues[numDBValues];
  const char* dbValueNames[numDBValues] = {"displacement-dot-x", "displacement-dot-y"};
  queryctx->db->queryVals(dbValueNames, numDBValues);

  // Dimensionalize query location coordinates.
  assert(queryctx->lengthScale > 0);
  double xDim[_dim];
  for (int i=0; i < _dim; ++i) {
    xDim[i] = x[i] * queryctx->lengthScale;
  } // for

  assert(queryctx->cs);
  const int err = queryctx->db->query(dbValues, numDBValues, xDim, _dim, queryctx->cs);
  if (err) {
    std::ostringstream msg;
    msg << "Could not find time derivative of displacement field at (";
    for (int i=0; i < _dim; ++i)
      msg << "  " << xDim[i];
    msg << ") in using spatial database '" << queryctx->db->label() << "'.";
    (*PetscErrorPrintf)(msg.str().c_str());
    PYLITH_METHOD_RETURN(1);
  } // if
  
  for (int i=0; i < _dim; ++i) {
    values[i] = dbValues[i] / queryctx->valueScale;
  } // for

  PYLITH_METHOD_RETURN(0);
} // _querySolutionDisplacementDot


// Function for computing time derivative of velocity solution to match LHS and RHS residual functions.
PetscErrorCode
pylith::materials::IsotropicLinearElasticityPlaneStrainData_Tri3::_querySolutionVelocityDot(PylithInt dim, 
											    PylithReal t, 
											    const PylithReal x[],
											    PylithInt nvalues,
											    PylithScalar* values, 
											    void* context)
{ // _querySolutionVelocityDot
  PYLITH_METHOD_BEGIN;

  const int _dim = 2;
  const int _nvalues = _dim;

  assert(x);
  assert(values);
  assert(context);
  assert(_nvalues == nvalues);
  assert(_dim == dim);

  const pylith::topology::FieldQuery::DBQueryContext* queryctx = (pylith::topology::FieldQuery::DBQueryContext*)context;assert(queryctx);

  // Tell database which values we want.
  const int numDBValues = _dim;
  PylithReal dbValues[numDBValues];
  const char* dbValueNames[numDBValues] = {"velocity-dot-x", "velocity-dot-y"};
  queryctx->db->queryVals(dbValueNames, numDBValues);

  // Dimensionalize query location coordinates.
  assert(queryctx->lengthScale > 0);
  double xDim[_dim];
  for (int i=0; i < _dim; ++i) {
    xDim[i] = x[i] * queryctx->lengthScale;
  } // for

  assert(queryctx->cs);
  const int err = queryctx->db->query(dbValues, numDBValues, xDim, _dim, queryctx->cs);
  if (err) {
    std::ostringstream msg;
    msg << "Could not find time derivative of velocity field at (";
    for (int i=0; i < _dim; ++i)
      msg << "  " << xDim[i];
    msg << ") in using spatial database '" << queryctx->db->label() << "'.";
    (*PetscErrorPrintf)(msg.str().c_str());
    PYLITH_METHOD_RETURN(1);
  } // if
  
  for (int i=0; i < _dim; ++i) {
    values[i] = dbValues[i] / queryctx->valueScale;
  } // for

  PYLITH_METHOD_RETURN(0);
} // _querySolutionVelocityDot


// End of file
