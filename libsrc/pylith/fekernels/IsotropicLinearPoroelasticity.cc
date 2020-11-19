/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University of Chicago
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2017 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/IsotropicLinearPoroelasticity.hh"
#include "pylith/fekernels/Poroelasticity.hh" // USES Poroelasticity kernels
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include <cassert> // USES assert()

// =====================================================================================================================
// Kernels for isotropic, linear poroelasticity plane strain.
// =====================================================================================================================
// ----------------------------------------------------------------------

// ================================= MMS =======================================

// ----------------------------------------------------------------------
// f0u function for quadratic space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0_mms_ql_u(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar f0u[]) {
 const PylithInt _dim = 2;

 // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 // Incoming re-packed auxiliary field.

 // Poroelasticity

 // IsotropicLinearPoroelasticity
 const PylithInt i_shearModulus = numA - 5;
 const PylithInt i_drainedBulkModulus = numA - 4;
 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar shearModulus = a[aOff[i_shearModulus]];
 const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];

 const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

  for (PylithInt d = 0; d < _dim-1; ++d) {
    f0u[d] -= 2.0*shearModulus - biotCoefficient*t;
  }
  f0u[_dim-1] -= 2.0*lambda + 4.0*shearModulus - biotCoefficient*t;
} // f0_quadratic_linear_u

// ----------------------------------------------------------------------
// f0p function for quadratic space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0_mms_ql_p(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar f0p[]) {

 const PylithInt _dim = 2;

  // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];

 PylithScalar       sum    = 0.0;

 for (PylithInt d = 0; d < _dim; ++d) {
   sum += x[d];
 }
 f0p[0] += s_t ? biotCoefficient*s_t[sOff[i_trace_strain]] : 0.0;
 f0p[0] += s_t ? s_t[sOff[i_pressure]]/biotModulus     : 0.0;
 f0p[0] -= sum/biotModulus;
} // f0_quadratic_linear_p

// ----------------------------------------------------------------------
// f0u function for quadratic space and trigonometric time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0_mms_qt_u(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar f0u[]) {

 const PylithInt _dim = 2;

 // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 // Incoming re-packed auxiliary field.

 // Poroelasticity

 // IsotropicLinearPoroelasticity
 const PylithInt i_shearModulus = numA - 5;
 const PylithInt i_drainedBulkModulus = numA - 4;
 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar shearModulus = a[aOff[i_shearModulus]];
 const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];

 const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

  for (PylithInt d = 0; d < _dim-1; ++d) {
    f0u[d] -= 2.0*shearModulus - biotCoefficient*PetscCosReal(t);
  }
  f0u[_dim-1] -= 2.0*lambda + 4.0*shearModulus - biotCoefficient*PetscCosReal(t);
} // f0_quadratic_trig_u

// ----------------------------------------------------------------------
// f0p function for quadratic space and trigonometric time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0_mms_qt_p(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar f0p[]) {

 const PylithInt _dim = 2;

 // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];

 PylithScalar       sum    = 0.0;

 for (PylithInt d = 0; d < _dim; ++d) {
   sum += x[d];
 }
 f0p[0] += s_t ? biotCoefficient*s_t[sOff[i_trace_strain]] : 0.0;
 f0p[0] += s_t ? s_t[sOff[i_pressure]]/biotModulus     : 0.0;
 f0p[0] += PetscSinReal(t)*sum/biotModulus;
} // f0_quadratic_trig_p

// ----------------------------------------------------------------------
// f0u function for trigonometric space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0_mms_tl_u(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar f0u[]) {

 const PylithInt _dim = 2;

 // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 // Incoming re-packed auxiliary field.

 // Poroelasticity

 // IsotropicLinearPoroelasticity
 const PylithInt i_shearModulus = numA - 5;
 const PylithInt i_drainedBulkModulus = numA - 4;
 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar shearModulus = a[aOff[i_shearModulus]];
 const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];

 const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

 for (PylithInt d = 0; d < _dim-1; ++d) {
   f0u[d] += PetscSqr(2.*PETSC_PI)*PetscSinReal(2.*PETSC_PI*x[d])*(2.*shearModulus + lambda) + 2.0*(shearModulus + lambda) - 2.*PETSC_PI*biotCoefficient*PetscSinReal(2.*PETSC_PI*x[d])*t;
 }
 f0u[_dim-1] += PetscSqr(2.*PETSC_PI)*PetscSinReal(2.*PETSC_PI*x[_dim-1])*(2.*shearModulus + lambda) - 2.*PETSC_PI*biotCoefficient*PetscSinReal(2.*PETSC_PI*x[_dim-1])*t;
} // f0_trig_linear_u

// ----------------------------------------------------------------------
// f0p function for trigonometric space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0_mms_tl_p(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar f0p[]) {

 const PylithInt _dim = 2;

 // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 // Poroelasticity
 const PylithInt i_fluidViscosity = 2;

 // IsotropicLinearPoroelasticity
 const PylithInt i_isotropicPerm = numA - 1;
 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];
 const PylithScalar kappa = a[aOff[i_isotropicPerm]] / a[aOff[i_fluidViscosity]];
 PylithScalar       sum    = 0.0;

 for (PylithInt d = 0; d < _dim; ++d) {
   sum += PetscCosReal(2.*PETSC_PI*x[d]);
 }
 f0p[0] += s_t ? biotCoefficient*s_t[i_trace_strain] : 0.0;
 f0p[0] += s_t ? s_t[sOff[i_pressure]]/biotModulus     : 0.0;
 f0p[0] -= sum/biotModulus - 4*PetscSqr(PETSC_PI)*kappa*sum*t;
} // f0_quadratic_trig_p

// ================================= STD =======================================


// ================================= LHS =======================================


// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_quasistatic(const PylithInt dim,
                                                                  const PylithInt numS,
                                                                  const PylithInt numA,
                                                                  const PylithInt sOff[],
                                                                  const PylithInt sOff_x[],
                                                                  const PylithScalar s[],
                                                                  const PylithScalar s_t[],
                                                                  const PylithScalar s_x[],
                                                                  const PylithInt aOff[],
                                                                  const PylithInt aOff_x[],
                                                                  const PylithScalar a[],
                                                                  const PylithScalar a_t[],
                                                                  const PylithScalar a_x[],
                                                                  const PylithReal t,
                                                                  const PylithScalar x[],
                                                                  const PylithInt numConstants,
                                                                  const PylithScalar constants[],
                                                                  PylithScalar f0p[]) {
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Incoming re-packed auxiliary field.

    // Poroelasticity

    // IsotropicLinearPoroelasticity
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0p[0] += biotCoefficient*trace_strain_t;
    f0p[0] += pressure_t/biotModulus;
} // f0p_QS

// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_inertia(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar f0p[]) {

    const PylithInt _dim = 2;

    // Incoming re-packed solution field.
    const PylithInt i_poro_pres = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity

    // IsotropicLinearPoroelasticity
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar poro_pres_t = s_t[sOff[i_poro_pres]];
    const PylithScalar* vel_x = &s_x[sOff[i_velocity]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
      trace_strain_t += vel_x[d*_dim+d];
    }

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0p[0] += biotCoefficient*trace_strain_t;
    f0p[0] += poro_pres_t/biotModulus;
} // f0p_DYN

// -----------------------------------------------------------------------------
// Jf0pe function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf0pe(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal utshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jf0[]) {

    const PylithInt _dim = 2;

    // Incoming re-packed auxiliary field.
    // IsotropicLinearPoroelasticity
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    Jf0[0] = utshift * biotCoefficient;
} // Jf0pe

// -----------------------------------------------------------------------------
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf0pp(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal utshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jf0[]) {

    const PylithInt _dim = 2;

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticity
    const PylithInt i_biotModulus = numA - 2;
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    Jf0[0] = utshift / biotModulus;
} // Jf0pp


// ============================== RHS ==========================================

// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1p_gravity(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar g1p[]) {

    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticity
    const PylithInt i_isotropicPerm = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];
    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar isotropicPerm = a[aOff[i_isotropicPerm]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    for (PylithInt d = 0; d < _dim; ++d) {
        g1p[d] -= (isotropicPerm / fluidViscosity) * (pressure_x[d] - fluidDensity*gravityField[d]);
    } // for

} // g1p_gravity

// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1p_gravity_tensor_permeability(const PylithInt dim,
                                                                                  const PylithInt numS,
                                                                                  const PylithInt numA,
                                                                                  const PylithInt sOff[],
                                                                                  const PylithInt sOff_x[],
                                                                                  const PylithScalar s[],
                                                                                  const PylithScalar s_t[],
                                                                                  const PylithScalar s_x[],
                                                                                  const PylithInt aOff[],
                                                                                  const PylithInt aOff_x[],
                                                                                  const PylithScalar a[],
                                                                                  const PylithScalar a_t[],
                                                                                  const PylithScalar a_x[],
                                                                                  const PylithReal t,
                                                                                  const PylithScalar x[],
                                                                                  const PylithInt numConstants,
                                                                                  const PylithScalar constants[],
                                                                                  PylithScalar g1p[]) {

    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticity
    const PylithInt i_tensorPerm = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];
    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar* vectorPerm = &a[aOff[i_tensorPerm]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    PylithScalar tensorPerm[4];
    tensorPerm[0] = vectorPerm[0];
    tensorPerm[1] = vectorPerm[3];
    tensorPerm[2] = vectorPerm[1];
    tensorPerm[3] = vectorPerm[3];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            g1p[i] -= (tensorPerm[i*_dim+j] / fluidViscosity) * (pressure_x[j] - fluidDensity*gravityField[j]);
        } // for
    } // for

} // g1p_gravity_tensor_permeability


// ----------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITHOUT gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1p(const PylithInt dim,
                                                      const PylithInt numS,
                                                      const PylithInt numA,
                                                      const PylithInt sOff[],
                                                      const PylithInt sOff_x[],
                                                      const PylithScalar s[],
                                                      const PylithScalar s_t[],
                                                      const PylithScalar s_x[],
                                                      const PylithInt aOff[],
                                                      const PylithInt aOff_x[],
                                                      const PylithScalar a[],
                                                      const PylithScalar a_t[],
                                                      const PylithScalar a_x[],
                                                      const PylithReal t,
                                                      const PylithScalar x[],
                                                      const PylithInt numConstants,
                                                      const PylithScalar constants[],
                                                      PylithScalar g1p[]) {

    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticity
    const PylithInt i_isotropicPerm = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar isotropicPerm = a[aOff[i_isotropicPerm]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    for (PylithInt d = 0; d < _dim; ++d) {
        g1p[d] -= (isotropicPerm / fluidViscosity) * pressure_x[d];
    } // for
} // g1p

// ----------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITHOUT gravity, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1p_tensor_permeability(const PylithInt dim,
                                                     const PylithInt numS,
                                                     const PylithInt numA,
                                                     const PylithInt sOff[],
                                                     const PylithInt sOff_x[],
                                                     const PylithScalar s[],
                                                     const PylithScalar s_t[],
                                                     const PylithScalar s_x[],
                                                     const PylithInt aOff[],
                                                     const PylithInt aOff_x[],
                                                     const PylithScalar a[],
                                                     const PylithScalar a_t[],
                                                     const PylithScalar a_x[],
                                                     const PylithReal t,
                                                     const PylithScalar x[],
                                                     const PylithInt numConstants,
                                                     const PylithScalar constants[],
                                                     PylithScalar g1p[]) {

    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticity
    const PylithInt i_tensorPerm = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar* vectorPerm = &a[aOff[i_tensorPerm]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    PylithScalar tensorPerm[4];
    tensorPerm[0] = vectorPerm[0];
    tensorPerm[1] = vectorPerm[3];
    tensorPerm[2] = vectorPerm[1];
    tensorPerm[3] = vectorPerm[3];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            g1p[i] -= (tensorPerm[i*_dim+j] / fluidViscosity) * (pressure_x[j]);
        } // for
    } // for
} // g1p_tensor_permeability

// ----------------------------------------------------------------------
// g1u function for isotropic linear Poroelasticity WITHOUT reference stress and strain.
// Quasi - Static Case
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1u(const PylithInt dim,
                                                      const PylithInt numS,
                                                      const PylithInt numA,
                                                      const PylithInt sOff[],
                                                      const PylithInt sOff_x[],
                                                      const PylithScalar s[],
                                                      const PylithScalar s_t[],
                                                      const PylithScalar s_x[],
                                                      const PylithInt aOff[],
                                                      const PylithInt aOff_x[],
                                                      const PylithScalar a[],
                                                      const PylithScalar a_t[],
                                                      const PylithScalar a_x[],
                                                      const PylithReal t,
                                                      const PylithScalar x[],
                                                      const PylithInt numConstants,
                                                      const PylithScalar constants[],
                                                      PylithScalar g1[]) {

    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar pressure = s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticity
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    for (PylithInt c = 0; c < _dim; ++c) {
      for (PylithInt d = 0; d < _dim; ++d) {
        g1[c*_dim+d] += shearModulus * (disp_x[c*_dim+d] + disp_x[d*_dim+c]);
      } // for
      g1[c*_dim+c] += (drainedBulkModulus - (2.0*shearModulus)/3.0) * trace_strain;
      g1[c*_dim+c] -= biotCoefficient*pressure;
    } // for
} // g1u

// ----------------------------------------------------------------------
// g1u function for isotropic linear Poroelasticity plane strain with reference stress and strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1u_refstate(const PylithInt dim,
                                                               const PylithInt numS,
                                                               const PylithInt numA,
                                                               const PylithInt sOff[],
                                                               const PylithInt sOff_x[],
                                                               const PylithScalar s[],
                                                               const PylithScalar s_t[],
                                                               const PylithScalar s_x[],
                                                               const PylithInt aOff[],
                                                               const PylithInt aOff_x[],
                                                               const PylithScalar a[],
                                                               const PylithScalar a_t[],
                                                               const PylithScalar a_x[],
                                                               const PylithReal t,
                                                               const PylithScalar x[],
                                                               const PylithInt numConstants,
                                                               const PylithScalar constants[],
                                                               PylithScalar g1[]) {

   const PylithInt _dim = 2;

   // Incoming solution fields.
   const PylithInt i_disp = 0;
   const PylithInt i_pressure = 1;
   const PylithInt i_trace_strain = 2;

   const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
   const PylithScalar pressure = s[sOff[i_pressure]];
   const PylithScalar trace_strain = s[sOff[i_trace_strain]];

   // Incoming auxiliary fields.

   // IsotropicLinearPoroelasticity
   const PylithInt i_rstress = numA - 7;
   const PylithInt i_rstrain = numA - 6;
   const PylithInt i_shearModulus = numA - 5;
   const PylithInt i_drainedBulkModulus = numA - 4;
   const PylithInt i_biotCoefficient = numA - 3;
   const PylithInt i_biotModulus = numA - 2;

   const PylithScalar shearModulus = a[aOff[i_shearModulus]];
   const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
   const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
   const PylithScalar biotModulus = a[aOff[i_biotModulus]];
   const PylithScalar* refStress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
   const PylithScalar* refStrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

   const PylithScalar ref_trace_strain = refStrain[0] + refStrain[1] + refStrain[2];
   const PylithScalar meanrstress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
   const PylithScalar meanStress = meanrstress + drainedBulkModulus * (trace_strain - ref_trace_strain);
   const PylithScalar alphaPres = biotCoefficient * pressure;
   const PylithReal traceTerm = (-2.0/3.0)*shearModulus * trace_strain;

   // Convert reference vectors to refrence tensors
   PylithScalar refStressTensor[_dim*_dim];
   PylithScalar refStrainTensor[_dim*_dim];
   PylithInt refTensorPos[9] = {0, 3, 5, 3, 1, 4, 5, 4, 2};

   for (PylithInt i = 0; i < _dim; ++i) {
     for (PylithInt j = 0; j < _dim; ++j) {
       refStressTensor[i*_dim+j] = refStress[refTensorPos[i*_dim+j]];
       refStrainTensor[i*_dim+j] = refStrain[refTensorPos[i*_dim+j]];
     } // for
   } // for

   for (PylithInt i = 0; i < _dim; ++i) {
       g1[i*_dim+i] -= (meanStress - alphaPres);
       g1[i*_dim+i] -= refStress[i*_dim+i] - meanrstress + traceTerm;
       for (PylithInt j = 0; j < _dim; ++j) {
           g1[i*_dim+j] -= shearModulus * (disp_x[i*_dim+j] + disp_x[j*_dim+i]) - refStrain[i*_dim+j];
       } // for
   } // for
} // g1u_refstate


// ----------------------------------------------------------------------
// g1v function for isotropic linear Poroelasticity WITHOUT reference stress and strain.
// Dynamic Case
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1v(const PylithInt dim,
                                                      const PylithInt numS,
                                                      const PylithInt numA,
                                                      const PylithInt sOff[],
                                                      const PylithInt sOff_x[],
                                                      const PylithScalar s[],
                                                      const PylithScalar s_t[],
                                                      const PylithScalar s_x[],
                                                      const PylithInt aOff[],
                                                      const PylithInt aOff_x[],
                                                      const PylithScalar a[],
                                                      const PylithScalar a_t[],
                                                      const PylithScalar a_x[],
                                                      const PylithReal t,
                                                      const PylithScalar x[],
                                                      const PylithInt numConstants,
                                                      const PylithScalar constants[],
                                                      PylithScalar g1[]) {

   const PylithInt _dim = 2;

   // Incoming solution fields.
   const PylithInt i_disp = 0;
   const PylithInt i_pressure = 1;

   const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
   const PylithScalar pressure = s[sOff[i_pressure]];

   PylithScalar trace_strain = 0.0;
   for (PylithInt d = 0; d < _dim; ++d) {
     trace_strain += disp_x[d*_dim+d];
   }

   // Incoming auxiliary fields.

   // IsotropicLinearPoroelasticity
   const PylithInt i_shearModulus = numA - 5;
   const PylithInt i_drainedBulkModulus = numA - 4;
   const PylithInt i_biotCoefficient = numA - 3;
   const PylithInt i_biotModulus = numA - 2;

   const PylithScalar shearModulus = a[aOff[i_shearModulus]];
   const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
   const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
   const PylithScalar biotModulus = a[aOff[i_biotModulus]];

   for (PylithInt c = 0; c < _dim; ++c) {
     for (PylithInt d = 0; d < _dim; ++d) {
       g1[c*dim+d] += shearModulus * (disp_x[c*_dim+d] + disp_x[d*_dim+c]);
     } // for
     g1[c*dim+c] += (drainedBulkModulus - (2.0*shearModulus)/3.0) * trace_strain;
     g1[c*dim+c] -= biotCoefficient*pressure;
   } // for
} // g1v

// ----------------------------------------------------------------------
// g1v function for isotropic linear Poroelasticity plane strain with reference stress and strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1v_refstate(const PylithInt dim,
                                                               const PylithInt numS,
                                                               const PylithInt numA,
                                                               const PylithInt sOff[],
                                                               const PylithInt sOff_x[],
                                                               const PylithScalar s[],
                                                               const PylithScalar s_t[],
                                                               const PylithScalar s_x[],
                                                               const PylithInt aOff[],
                                                               const PylithInt aOff_x[],
                                                               const PylithScalar a[],
                                                               const PylithScalar a_t[],
                                                               const PylithScalar a_x[],
                                                               const PylithReal t,
                                                               const PylithScalar x[],
                                                               const PylithInt numConstants,
                                                               const PylithScalar constants[],
                                                               PylithScalar g1[]) {

    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar pressure = s[sOff[i_pressure]];

    PylithScalar trace_strain = 0.0;
    for (PylithInt d = 0; d < _dim; ++d) {
     trace_strain += disp_x[d*_dim+d];
    }

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticity
    const PylithInt i_rstress = numA - 7;
    const PylithInt i_rstrain = numA - 6;
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar* refStress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refStrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithScalar ref_trace_strain = refStrain[0] + refStrain[1] + refStrain[2];
    const PylithScalar meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
    const PylithScalar meanStress = meanRefStress + drainedBulkModulus * (trace_strain - ref_trace_strain);
    const PylithScalar alphaPres = biotCoefficient * pressure;
    const PylithReal traceTerm = (-2.0/3.0)*shearModulus * trace_strain;

    // Convert reference vectors to refrence tensors
    PylithScalar refStressTensor[_dim*_dim];
    PylithScalar refStrainTensor[_dim*_dim];
    PylithInt refTensorPos[9] = {0, 3, 5, 3, 1, 4, 5, 4, 2};

    for (PylithInt i = 0; i < _dim; ++i) {
      for (PylithInt j = 0; j < _dim; ++j) {
        refStressTensor[i*_dim+j] = refStress[refTensorPos[i*_dim+j]];
        refStrainTensor[i*_dim+j] = refStrain[refTensorPos[i*_dim+j]];
      } // for
    } // for

    for (PylithInt i = 0; i < _dim; ++i) {
        g1[i*_dim+i] -= (meanStress - alphaPres);
        g1[i*_dim+i] -= refStressTensor[i*_dim+i] - meanRefStress + traceTerm;
        for (PylithInt j = 0; j < _dim; ++j) {
            g1[i*_dim+j] -= shearModulus * (disp_x[i*_dim+j] + disp_x[j*_dim+i]) - refStrainTensor[i*_dim+j];
        } // for
    } // for
} // g1v_refstate

// RHS Jacobians

// -----------------------------------------------------------------------------
// Jg2ue function for isotropic linear poroelasticity.

void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg2ue(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal utshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jg2[]) {

    const PylithInt _dim = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    for (PylithInt d = 0; d < _dim; ++d) {
        Jg2[d*_dim+d] += drainedBulkModulus - (2.0*shearModulus) / 3.0;
    } // for
} // Jg2ue


// -----------------------------------------------------------------------------
// Jg2up function for isotropic linear poroelasticity.

void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg2up(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal utshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jg2[]) {

    const PylithInt _dim = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    for (PylithInt d = 0; d < _dim; ++d) {
        Jg2[d*_dim+d] -= biotCoefficient;
    } // for
} // Jg2up

// -----------------------------------------------------------------------------
// Jg2vp function for isotropic linear poroelasticity.
// vp refers to dynamic formulation (velocity / pressure)
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg2vp(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal utshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jg2[]) {

    const PylithInt _dim = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    for (PylithInt d = 0; d < _dim; ++d) {
        Jg2[d*_dim+d] -= biotCoefficient ;
    } // for
} // Jg2vp

// ---------------------------------------------------------------------------------------------------------------------
/* Jg3_uu entry function for isotropic linear poroelasticity WITHOUT reference stress and reference strain.
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg3uu(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal s_tshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jg3[]) {

    const PylithInt _dim = 2;

    // index of Incoming auxiliary fields.
    // Poroelasticity

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    for (PylithInt c = 0; c < _dim; ++c) {
      for (PylithInt d = 0; d < _dim; ++d) {
        Jg3[((c*_dim + c)*_dim + d)*_dim + d] += shearModulus;
        Jg3[((c*_dim + d)*_dim + d)*_dim + c] += shearModulus;
      } // for
    } // for
} // Jg3uu

// ----------------------------------------------------------------------
/* Jg3pp entry function for isotropic linear poroelasticity. */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg3pp(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal utshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jg3[]) {

    const PylithInt _dim = 2;

    // index of Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_isotropicPermeability = numA - 1;

    const PylithScalar isotropicPermeablity = a[aOff[i_isotropicPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    for (PylithInt d = 0; d < _dim; ++d) {
      Jg3[d*_dim+d] -= isotropicPermeablity/fluidViscosity;
    }
} // Jg3pp

// ----------------------------------------------------------------------
/* Jg3pp entry function for isotropic linear poroelasticity, permeability in tensor form */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg3pp_tensor_permeability(const PylithInt dim,
                                                                            const PylithInt numS,
                                                                            const PylithInt numA,
                                                                            const PylithInt sOff[],
                                                                            const PylithInt sOff_x[],
                                                                            const PylithScalar s[],
                                                                            const PylithScalar s_t[],
                                                                            const PylithScalar s_x[],
                                                                            const PylithInt aOff[],
                                                                            const PylithInt aOff_x[],
                                                                            const PylithScalar a[],
                                                                            const PylithScalar a_t[],
                                                                            const PylithScalar a_x[],
                                                                            const PylithReal t,
                                                                            const PylithReal utshift,
                                                                            const PylithScalar x[],
                                                                            const PylithInt numConstants,
                                                                            const PylithScalar constants[],
                                                                            PylithScalar Jg3[]) {

    const PylithInt _dim = 2;

    // index of Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_tensorPerm = numA - 1;

    const PylithScalar* vectorPerm = &a[aOff[i_tensorPerm]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    PylithScalar tensorPerm[4];
    tensorPerm[0] = vectorPerm[0];
    tensorPerm[1] = vectorPerm[3];
    tensorPerm[2] = vectorPerm[1];
    tensorPerm[3] = vectorPerm[3];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            Jg3[i*_dim+j] -= tensorPerm[i*_dim+j]/fluidViscosity;
        } // for
    } // for

} // Jg3pp_tensorPerm


// ---------------------------------------------------------------------------------------------------------------------
/* Jg3_vu entry function for isotropic linear poroelasticity WITHOUT reference stress and reference strain.
 *
 * stress_ij = C_ijkl strain_kl
 *
 * stress_11 = C1111 strain_11 + C1122 strain_22, C1111=lambda+2mu, C1122=lambda.
 *
 * stress_12 = C1212 strain_12 + C1221 strain_21. C1212 = C1221 from symmetry, so C1212 = C1221 = shearModulus.
 *
 * For reference:
 *
 * Isotropic:
 *  C_ijkl = bulkModulus * delta_ij * delta_kl
 *   + shearModulus * (delta_ik*delta_jl + delta_il*delta*jk - 2/3*delta_ij*delta_kl)
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg3vu(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal s_tshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jg3[]) {
    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_disp = 0;
//    const PylithInt i_poro_pres = 1;

    // Incoming auxiliary fields.

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
//    const PylithScalar poro_pres = s[sOff[i_poro_pres]];

    PylithReal strainTrace = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
      strainTrace += disp_x[d*_dim+d];
    } // for

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    const PylithScalar meanStress = drainedBulkModulus * strainTrace;
    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar lambda2mu = lambda + 2.0*shearModulus;

    const PylithReal C1111 = lambda2mu;// Get auxiliary factory associated with physics.
    const PylithReal C2222 = lambda2mu;
    const PylithReal C1122 = lambda;
    const PylithReal C1212 = shearModulus;

    /* j(f,g,df,dg) = C(f,df,g,dg)
     *
     * 0: j0000 = C1111
     * 1: j0001 = C1112 = 0
     * 4: j0100 = C1121, symmetry C1112 = 0
     * 5: j0101 = C1122
     *
     * 2: j0010 = C1211 = 0
     * 3: j0011 = C1212
     * 6: j0110 = C1221, symmetry C1212
     * 7: j0111 = C1222 = 0
     *
     * 8: j1000 = C2111 = 0
     * 9: j1001 = C2112, symmetry C1212
     * 12: j1100 = C2121, symmetry C1212
     * 13: j1101 = C2122, symmetry C1222 = 0
     *
     * 10: j1010 = C2211, symmetry C1122
     * 11: j1011 = C2212, symmetry C1222 = 0
     * 14: j1110 = C2221, symmetry C1222 = 0
     * 15: j1111 = C2222
     */

    Jg3[ 0] -= C1111; // j0000
    Jg3[ 3] -= C1212; // j0011
    Jg3[ 5] -= C1122; // j0101
    Jg3[ 6] -= C1212; // j0110, C1221
    Jg3[ 9] -= C1212; // j1001, C2112
    Jg3[10] -= C1122; // j1010, C2211
    Jg3[12] -= C1212; // j1100, C2121
    Jg3[15] -= C2222; // j1111

} // Jg3vu

// ========================== Helper Kernels ===================================


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate stress vector for isotropic linear
 * poroelasticity WITHOUT a reference stress and strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim), ...]
 * Auxiliary fields: [density(1), ..., shear_modulus(1), bulk_modulus(1)]
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::cauchyStress(const PylithInt dim,
                                                               const PylithInt numS,
                                                               const PylithInt numA,
                                                               const PylithInt sOff[],
                                                               const PylithInt sOff_x[],
                                                               const PylithScalar s[],
                                                               const PylithScalar s_t[],
                                                               const PylithScalar s_x[],
                                                               const PylithInt aOff[],
                                                               const PylithInt aOff_x[],
                                                               const PylithScalar a[],
                                                               const PylithScalar a_t[],
                                                               const PylithScalar a_x[],
                                                               const PylithReal t,
                                                               const PylithScalar x[],
                                                               const PylithInt numConstants,
                                                               const PylithScalar constants[],
                                                               PylithScalar stressVector[]) {

    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar pressure = s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    // Incoming auxiliary fields.
    // Poroelasticity

    // IsotropicLinearPoroelasticity
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    // Create and populate stress tensor

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            stressTensor[i*_dim+j] += shearModulus * (disp_x[i*_dim+j] + disp_x[j*_dim+i]);
      } // for
      stressTensor[i*_dim+i] += (drainedBulkModulus - (2.0*shearModulus)/3.0) * trace_strain;
      stressTensor[i*_dim+i] -= biotCoefficient*pressure;
    } // for

    // Construct stress vector
    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar stress_zz = 0.5*lambda/(lambda+shearModulus) * (stressTensor[0*_dim+0] + stressTensor[1*_dim+1]);

    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stress_zz;
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
} // cauchyStress

// ---------------------------------------------------------------------------------------------------------------------
/* Calculate stress vector for isotropic linear
 * poroelasticity WITH a reference stress/strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim), ...]
 * Auxiliary fields: [density(1), ..., refstress(4), refstrain(4), shear_modulus(1), bulk_modulus(1)]
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::cauchyStress_refstate(const PylithInt dim,
                                                                        const PylithInt numS,
                                                                        const PylithInt numA,
                                                                        const PylithInt sOff[],
                                                                        const PylithInt sOff_x[],
                                                                        const PylithScalar s[],
                                                                        const PylithScalar s_t[],
                                                                        const PylithScalar s_x[],
                                                                        const PylithInt aOff[],
                                                                        const PylithInt aOff_x[],
                                                                        const PylithScalar a[],
                                                                        const PylithScalar a_t[],
                                                                        const PylithScalar a_x[],
                                                                        const PylithReal t,
                                                                        const PylithScalar x[],
                                                                        const PylithInt numConstants,
                                                                        const PylithScalar constants[],
                                                                        PylithScalar stressVector[]) {

   const PylithInt _dim = 2;

   // Incoming solution fields.
   const PylithInt i_disp = 0;
   const PylithInt i_pressure = 1;
   const PylithInt i_trace_strain = 2;

   const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
   const PylithScalar pressure = s[sOff[i_pressure]];
   const PylithScalar trace_strain = s[sOff[i_trace_strain]];

   // Incoming auxiliary fields.

   // IsotropicLinearPoroelasticity
   const PylithInt i_rstress = numA - 7;
   const PylithInt i_rstrain = numA - 6;
   const PylithInt i_shearModulus = numA - 5;
   const PylithInt i_drainedBulkModulus = numA - 4;
   const PylithInt i_biotCoefficient = numA - 3;
   const PylithInt i_biotModulus = numA - 2;

   const PylithScalar shearModulus = a[aOff[i_shearModulus]];
   const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
   const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
   const PylithScalar biotModulus = a[aOff[i_biotModulus]];
   const PylithScalar* refStress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
   const PylithScalar* refStrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

   const PylithScalar ref_trace_strain = refStrain[0] + refStrain[1] + refStrain[2];
   const PylithScalar meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
   const PylithScalar meanStress = meanRefStress + drainedBulkModulus * (trace_strain - ref_trace_strain);
   const PylithScalar alphaPres = biotCoefficient * pressure;
   const PylithReal traceTerm = (-2.0/3.0)*shearModulus * trace_strain;

   // Convert reference vectors to refrence tensors
   PylithScalar refStressTensor[_dim*_dim];
   PylithScalar refStrainTensor[_dim*_dim];
   PylithInt refTensorPos[9] = {0, 3, 5, 3, 1, 4, 5, 4, 2};

   for (PylithInt i = 0; i < _dim; ++i) {
     for (PylithInt j = 0; j < _dim; ++j) {
       refStressTensor[i*_dim+j] = refStress[refTensorPos[i*_dim+j]];
       refStrainTensor[i*_dim+j] = refStrain[refTensorPos[i*_dim+j]];
     } // for
   } // for

   // Create and populate stress tensor

   PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };

   for (PylithInt i = 0; i < _dim; ++i) {
       stressTensor[i*_dim+i] -= (meanStress - alphaPres);
       stressTensor[i*_dim+i] -= refStressTensor[i*_dim+i] - meanRefStress + traceTerm;
       for (PylithInt j = 0; j < _dim; ++j) {
           stressTensor[i*_dim+j] -= shearModulus * (disp_x[i*_dim+j] + disp_x[j*_dim+i]) - refStrainTensor[i*_dim+j];
       } // for
   } // for

   // Generate stress vector

    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar stress_zz = refStress[2] +  0.5*lambda/(lambda+shearModulus) *
                                       (stressTensor[0*_dim+0]-refStress[0] + stressTensor[1*_dim+1]-refStress[1]);

    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stress_zz;
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
} // cauchyStress_refstate

// =====================================================================================================================
// Kernels for isotropic, linear poroelasticity in 3D.
// =====================================================================================================================

// ----------------------------------------------------------------------

// ================================= MMS =======================================

// ----------------------------------------------------------------------
// f0u function for quadratic space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::f0_mms_ql_u(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar f0u[]) {
 const PylithInt _dim = 3;

 // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 // Incoming re-packed auxiliary field.

 // Poroelasticity

 // IsotropicLinearPoroelasticity
 const PylithInt i_shearModulus = numA - 5;
 const PylithInt i_drainedBulkModulus = numA - 4;
 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar shearModulus = a[aOff[i_shearModulus]];
 const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];
 const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

  for (PylithInt d = 0; d < _dim-1; ++d) {
    f0u[d] -= 2.0*shearModulus - biotCoefficient*t;
  }
  f0u[_dim-1] -= 2.0*lambda + 4.0*shearModulus - biotCoefficient*t;
} // f0_quadratic_linear_u

// ----------------------------------------------------------------------
// f0p function for quadratic space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::f0_mms_ql_p(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar f0p[]) {

 const PylithInt _dim = 3;

  // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];

 PylithScalar       sum    = 0.0;

 for (PylithInt d = 0; d < _dim; ++d) {
   sum += x[d];
 }
 f0p[0] += s_t ? biotCoefficient*s_t[sOff[i_trace_strain]] : 0.0;
 f0p[0] += s_t ? s_t[sOff[i_pressure]]/biotModulus     : 0.0;
 f0p[0] -= sum/biotModulus;
} // f0_quadratic_linear_p

// ----------------------------------------------------------------------
// f0u function for quadratic space and trigonometric time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::f0_mms_qt_u(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar f0u[]) {

 const PylithInt _dim = 3;

 // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 // Incoming re-packed auxiliary field.

 // Poroelasticity

 // IsotropicLinearPoroelasticity
 const PylithInt i_shearModulus = numA - 5;
 const PylithInt i_drainedBulkModulus = numA - 4;
 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar shearModulus = a[aOff[i_shearModulus]];
 const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];
 const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

  for (PylithInt d = 0; d < _dim-1; ++d) {
    f0u[d] -= 2.0*shearModulus - biotCoefficient*PetscCosReal(t);
  }
  f0u[_dim-1] -= 2.0*lambda + 4.0*shearModulus - biotCoefficient*PetscCosReal(t);
} // f0_quadratic_trig_u

// ----------------------------------------------------------------------
// f0p function for quadratic space and trigonometric time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::f0_mms_qt_p(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar f0p[]) {

 const PylithInt _dim = 3;

 // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];

 PylithScalar       sum    = 0.0;

 for (PylithInt d = 0; d < _dim; ++d) {
   sum += x[d];
 }
 f0p[0] += s_t ? biotCoefficient*s_t[sOff[i_trace_strain]] : 0.0;
 f0p[0] += s_t ? s_t[sOff[i_pressure]]/biotModulus     : 0.0;
 f0p[0] += PetscSinReal(t)*sum/biotModulus;
} // f0_quadratic_trig_p

// ----------------------------------------------------------------------
// f0u function for trigonometric space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::f0_mms_tl_u(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar f0u[]) {

 const PylithInt _dim = 3;

 // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 // Incoming re-packed auxiliary field.

 // Poroelasticity

 // IsotropicLinearPoroelasticity
 const PylithInt i_shearModulus = numA - 5;
 const PylithInt i_drainedBulkModulus = numA - 4;
 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar shearModulus = a[aOff[i_shearModulus]];
 const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];
 const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

 for (PylithInt d = 0; d < _dim-1; ++d) {
   f0u[d] += PetscSqr(2.*PETSC_PI)*PetscSinReal(2.*PETSC_PI*x[d])*(2.*shearModulus + lambda) + 2.0*(shearModulus + lambda) - 2.*PETSC_PI*biotCoefficient*PetscSinReal(2.*PETSC_PI*x[d])*t;
 }
 f0u[_dim-1] += PetscSqr(2.*PETSC_PI)*PetscSinReal(2.*PETSC_PI*x[_dim-1])*(2.*shearModulus + lambda) - 2.*PETSC_PI*biotCoefficient*PetscSinReal(2.*PETSC_PI*x[_dim-1])*t;
} // f0_trig_linear_u

// ----------------------------------------------------------------------
// f0p function for trigonometric space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::f0_mms_tl_p(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar f0p[]) {

 const PylithInt _dim = 3;

 // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 // Poroelasticity
 const PylithInt i_fluidViscosity = 2;

 // IsotropicLinearPoroelasticity
 const PylithInt i_isotropicPerm = numA - 1;
 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];
 const PylithScalar kappa = a[aOff[i_isotropicPerm]] / a[aOff[i_fluidViscosity]];
 PylithScalar       sum    = 0.0;

 for (PylithInt d = 0; d < _dim; ++d) {
   sum += PetscCosReal(2.*PETSC_PI*x[d]);
 }
 f0p[0] += s_t ? biotCoefficient*s_t[i_trace_strain] : 0.0;
 f0p[0] += s_t ? s_t[sOff[i_pressure]]/biotModulus     : 0.0;
 f0p[0] -= sum/biotModulus - 4*PetscSqr(PETSC_PI)*kappa*sum*t;
} // f0_quadratic_trig_p

// ================================= STD =======================================


// ================================= LHS =======================================


// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_quasistatic(const PylithInt dim,
                                                                  const PylithInt numS,
                                                                  const PylithInt numA,
                                                                  const PylithInt sOff[],
                                                                  const PylithInt sOff_x[],
                                                                  const PylithScalar s[],
                                                                  const PylithScalar s_t[],
                                                                  const PylithScalar s_x[],
                                                                  const PylithInt aOff[],
                                                                  const PylithInt aOff_x[],
                                                                  const PylithScalar a[],
                                                                  const PylithScalar a_t[],
                                                                  const PylithScalar a_x[],
                                                                  const PylithReal t,
                                                                  const PylithScalar x[],
                                                                  const PylithInt numConstants,
                                                                  const PylithScalar constants[],
                                                                  PylithScalar f0p[]) {
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Incoming re-packed auxiliary field.

    // Poroelasticity

    // IsotropicLinearPoroelasticity
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0p[0] += biotCoefficient*trace_strain_t;
    f0p[0] += pressure_t/biotModulus;
} // f0p_QS

// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_inertia(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar f0p[]) {

    const PylithInt _dim = 3;

    // Incoming re-packed solution field.
    const PylithInt i_poro_pres = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity

    // IsotropicLinearPoroelasticity
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar poro_pres_t = s_t[sOff[i_poro_pres]];
    const PylithScalar* vel_x = &s_x[sOff[i_velocity]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
      trace_strain_t += vel_x[d*_dim+d];
    }

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0p[0] += biotCoefficient*trace_strain_t;
    f0p[0] += poro_pres_t/biotModulus;
} // f0p_DYN

// -----------------------------------------------------------------------------
// Jf0pe function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf0pe(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal utshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jf0[]) {

    const PylithInt _dim = 3;

    // Incoming re-packed auxiliary field.
    // IsotropicLinearPoroelasticity
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    Jf0[0] = utshift * biotCoefficient;
} // Jf0pe

// -----------------------------------------------------------------------------
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf0pp(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal utshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jf0[]) {

    const PylithInt _dim = 3;

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticity
    const PylithInt i_biotModulus = numA - 2;
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    Jf0[0] = utshift / biotModulus;
} // Jf0pp


// ============================== RHS ==========================================

// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::g1p_gravity(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
                                                              const PylithScalar s[],
                                                              const PylithScalar s_t[],
                                                              const PylithScalar s_x[],
                                                              const PylithInt aOff[],
                                                              const PylithInt aOff_x[],
                                                              const PylithScalar a[],
                                                              const PylithScalar a_t[],
                                                              const PylithScalar a_x[],
                                                              const PylithReal t,
                                                              const PylithScalar x[],
                                                              const PylithInt numConstants,
                                                              const PylithScalar constants[],
                                                              PylithScalar g1p[]) {

    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticity
    const PylithInt i_isotropicPerm = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];
    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar isotropicPerm = a[aOff[i_isotropicPerm]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    for (PylithInt d = 0; d < _dim; ++d) {
        g1p[d] -= (isotropicPerm / fluidViscosity) * (pressure_x[d] - fluidDensity*gravityField[d]);
    } // for

} // g1p_gravity

// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::g1p_gravity_tensor_permeability(const PylithInt dim,
                                                                                  const PylithInt numS,
                                                                                  const PylithInt numA,
                                                                                  const PylithInt sOff[],
                                                                                  const PylithInt sOff_x[],
                                                                                  const PylithScalar s[],
                                                                                  const PylithScalar s_t[],
                                                                                  const PylithScalar s_x[],
                                                                                  const PylithInt aOff[],
                                                                                  const PylithInt aOff_x[],
                                                                                  const PylithScalar a[],
                                                                                  const PylithScalar a_t[],
                                                                                  const PylithScalar a_x[],
                                                                                  const PylithReal t,
                                                                                  const PylithScalar x[],
                                                                                  const PylithInt numConstants,
                                                                                  const PylithScalar constants[],
                                                                                  PylithScalar g1p[]) {

    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticity
    const PylithInt i_tensorPerm = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];
    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar* vectorPerm = &a[aOff[i_tensorPerm]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    PylithScalar tensorPerm[9];
    tensorPerm[0] = vectorPerm[0];
    tensorPerm[1] = vectorPerm[3];
    tensorPerm[2] = vectorPerm[5];
    tensorPerm[3] = vectorPerm[3];
    tensorPerm[4] = vectorPerm[1];
    tensorPerm[5] = vectorPerm[4];
    tensorPerm[6] = vectorPerm[5];
    tensorPerm[7] = vectorPerm[4];
    tensorPerm[8] = vectorPerm[2];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            g1p[i] -= (tensorPerm[i*_dim+j] / fluidViscosity) * (pressure_x[j] - fluidDensity*gravityField[j]);
        } // for
    } // for

} // g1p_gravity_tensor_permeability


// ----------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITHOUT gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::g1p(const PylithInt dim,
                                                      const PylithInt numS,
                                                      const PylithInt numA,
                                                      const PylithInt sOff[],
                                                      const PylithInt sOff_x[],
                                                      const PylithScalar s[],
                                                      const PylithScalar s_t[],
                                                      const PylithScalar s_x[],
                                                      const PylithInt aOff[],
                                                      const PylithInt aOff_x[],
                                                      const PylithScalar a[],
                                                      const PylithScalar a_t[],
                                                      const PylithScalar a_x[],
                                                      const PylithReal t,
                                                      const PylithScalar x[],
                                                      const PylithInt numConstants,
                                                      const PylithScalar constants[],
                                                      PylithScalar g1p[]) {

    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticity
    const PylithInt i_isotropicPerm = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar isotropicPerm = a[aOff[i_isotropicPerm]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    for (PylithInt d = 0; d < _dim; ++d) {
        g1p[d] -= (isotropicPerm / fluidViscosity) * pressure_x[d];
    } // for
} // g1p

// ----------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITHOUT gravity, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::g1p_tensor_permeability(const PylithInt dim,
                                                     const PylithInt numS,
                                                     const PylithInt numA,
                                                     const PylithInt sOff[],
                                                     const PylithInt sOff_x[],
                                                     const PylithScalar s[],
                                                     const PylithScalar s_t[],
                                                     const PylithScalar s_x[],
                                                     const PylithInt aOff[],
                                                     const PylithInt aOff_x[],
                                                     const PylithScalar a[],
                                                     const PylithScalar a_t[],
                                                     const PylithScalar a_x[],
                                                     const PylithReal t,
                                                     const PylithScalar x[],
                                                     const PylithInt numConstants,
                                                     const PylithScalar constants[],
                                                     PylithScalar g1p[]) {

    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticity
    const PylithInt i_tensorPerm = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar* vectorPerm = &a[aOff[i_tensorPerm]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    PylithScalar tensorPerm[9];
    tensorPerm[0] = vectorPerm[0];
    tensorPerm[1] = vectorPerm[3];
    tensorPerm[2] = vectorPerm[5];
    tensorPerm[3] = vectorPerm[3];
    tensorPerm[4] = vectorPerm[1];
    tensorPerm[5] = vectorPerm[4];
    tensorPerm[6] = vectorPerm[5];
    tensorPerm[7] = vectorPerm[4];
    tensorPerm[8] = vectorPerm[2];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            g1p[i] -= (tensorPerm[i*_dim+j] / fluidViscosity) * (pressure_x[j]);
        } // for
    } // for
} // g1p_tensor_permeability

// ----------------------------------------------------------------------
// g1u function for isotropic linear Poroelasticity WITHOUT reference stress and strain.
// Quasi - Static Case
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::g1u(const PylithInt dim,
                                                      const PylithInt numS,
                                                      const PylithInt numA,
                                                      const PylithInt sOff[],
                                                      const PylithInt sOff_x[],
                                                      const PylithScalar s[],
                                                      const PylithScalar s_t[],
                                                      const PylithScalar s_x[],
                                                      const PylithInt aOff[],
                                                      const PylithInt aOff_x[],
                                                      const PylithScalar a[],
                                                      const PylithScalar a_t[],
                                                      const PylithScalar a_x[],
                                                      const PylithReal t,
                                                      const PylithScalar x[],
                                                      const PylithInt numConstants,
                                                      const PylithScalar constants[],
                                                      PylithScalar g1[]) {

    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar pressure = s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticity
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    for (PylithInt c = 0; c < _dim; ++c) {
      for (PylithInt d = 0; d < _dim; ++d) {
        g1[c*_dim+d] += shearModulus * (disp_x[c*_dim+d] + disp_x[d*_dim+c]);
      } // for
      g1[c*_dim+c] += (drainedBulkModulus - (2.0*shearModulus)/3.0) * trace_strain;
      g1[c*_dim+c] -= biotCoefficient*pressure;
    } // for
} // g1u

// ----------------------------------------------------------------------
// g1u function for isotropic linear Poroelasticity plane strain with reference stress and strain.
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::g1u_refstate(const PylithInt dim,
                                                               const PylithInt numS,
                                                               const PylithInt numA,
                                                               const PylithInt sOff[],
                                                               const PylithInt sOff_x[],
                                                               const PylithScalar s[],
                                                               const PylithScalar s_t[],
                                                               const PylithScalar s_x[],
                                                               const PylithInt aOff[],
                                                               const PylithInt aOff_x[],
                                                               const PylithScalar a[],
                                                               const PylithScalar a_t[],
                                                               const PylithScalar a_x[],
                                                               const PylithReal t,
                                                               const PylithScalar x[],
                                                               const PylithInt numConstants,
                                                               const PylithScalar constants[],
                                                               PylithScalar g1[]) {

   const PylithInt _dim = 3;

   // Incoming solution fields.
   const PylithInt i_disp = 0;
   const PylithInt i_pressure = 1;
   const PylithInt i_trace_strain = 2;

   const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
   const PylithScalar pressure = s[sOff[i_pressure]];
   const PylithScalar trace_strain = s[sOff[i_trace_strain]];

   // Incoming auxiliary fields.

   // IsotropicLinearPoroelasticity
   const PylithInt i_rstress = numA - 7;
   const PylithInt i_rstrain = numA - 6;
   const PylithInt i_shearModulus = numA - 5;
   const PylithInt i_drainedBulkModulus = numA - 4;
   const PylithInt i_biotCoefficient = numA - 3;
   const PylithInt i_biotModulus = numA - 2;

   const PylithScalar shearModulus = a[aOff[i_shearModulus]];
   const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
   const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
   const PylithScalar biotModulus = a[aOff[i_biotModulus]];
   const PylithScalar* refStress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
   const PylithScalar* refStrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

   const PylithScalar ref_trace_strain = refStrain[0] + refStrain[1] + refStrain[2];
   const PylithScalar meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
   const PylithScalar meanStress = meanRefStress + drainedBulkModulus * (trace_strain - ref_trace_strain);
   const PylithScalar alphaPres = biotCoefficient * pressure;
   const PylithReal traceTerm = (-2.0/3.0)*shearModulus * trace_strain;

   // Convert reference vectors to refrence tensors
   PylithScalar refStressTensor[_dim*_dim];
   PylithScalar refStrainTensor[_dim*_dim];
   PylithInt refTensorPos[9] = {0, 3, 5, 3, 1, 4, 5, 4, 2};

   for (PylithInt i = 0; i < _dim; ++i) {
     for (PylithInt j = 0; j < _dim; ++j) {
       refStressTensor[i*_dim+j] = refStress[refTensorPos[i*_dim+j]];
       refStrainTensor[i*_dim+j] = refStrain[refTensorPos[i*_dim+j]];
     } // for
   } // for

   for (PylithInt i = 0; i < _dim; ++i) {
       g1[i*_dim+i] -= (meanStress - alphaPres);
       g1[i*_dim+i] -= refStressTensor[i*_dim+i] - meanRefStress + traceTerm;
       for (PylithInt j = 0; j < _dim; ++j) {
           g1[i*_dim+j] -= shearModulus * (disp_x[i*_dim+j] + disp_x[j*_dim+i]) - refStrainTensor[i*_dim+j];
       } // for
   } // for
} // g1u_refstate


// ----------------------------------------------------------------------
// g1v function for isotropic linear Poroelasticity WITHOUT reference stress and strain.
// Dynamic Case
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::g1v(const PylithInt dim,
                                                      const PylithInt numS,
                                                      const PylithInt numA,
                                                      const PylithInt sOff[],
                                                      const PylithInt sOff_x[],
                                                      const PylithScalar s[],
                                                      const PylithScalar s_t[],
                                                      const PylithScalar s_x[],
                                                      const PylithInt aOff[],
                                                      const PylithInt aOff_x[],
                                                      const PylithScalar a[],
                                                      const PylithScalar a_t[],
                                                      const PylithScalar a_x[],
                                                      const PylithReal t,
                                                      const PylithScalar x[],
                                                      const PylithInt numConstants,
                                                      const PylithScalar constants[],
                                                      PylithScalar g1[]) {

   const PylithInt _dim = 3;

   // Incoming solution fields.
   const PylithInt i_disp = 0;
   const PylithInt i_pressure = 1;

   const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
   const PylithScalar pressure = s[sOff[i_pressure]];

   PylithScalar trace_strain = 0.0;
   for (PylithInt d = 0; d < _dim; ++d) {
     trace_strain += disp_x[d*_dim+d];
   }

   // Incoming auxiliary fields.

   // IsotropicLinearPoroelasticity
   const PylithInt i_shearModulus = numA - 5;
   const PylithInt i_drainedBulkModulus = numA - 4;
   const PylithInt i_biotCoefficient = numA - 3;
   const PylithInt i_biotModulus = numA - 2;

   const PylithScalar shearModulus = a[aOff[i_shearModulus]];
   const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
   const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
   const PylithScalar biotModulus = a[aOff[i_biotModulus]];

   for (PylithInt c = 0; c < _dim; ++c) {
     for (PylithInt d = 0; d < _dim; ++d) {
       g1[c*dim+d] += shearModulus * (disp_x[c*_dim+d] + disp_x[d*_dim+c]);
     } // for
     g1[c*dim+c] += (drainedBulkModulus - (2.0*shearModulus)/3.0) * trace_strain;
     g1[c*dim+c] -= biotCoefficient*pressure;
   } // for
} // g1v

// ----------------------------------------------------------------------
// g1v function for isotropic linear Poroelasticity plane strain with reference stress and strain.
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::g1v_refstate(const PylithInt dim,
                                                               const PylithInt numS,
                                                               const PylithInt numA,
                                                               const PylithInt sOff[],
                                                               const PylithInt sOff_x[],
                                                               const PylithScalar s[],
                                                               const PylithScalar s_t[],
                                                               const PylithScalar s_x[],
                                                               const PylithInt aOff[],
                                                               const PylithInt aOff_x[],
                                                               const PylithScalar a[],
                                                               const PylithScalar a_t[],
                                                               const PylithScalar a_x[],
                                                               const PylithReal t,
                                                               const PylithScalar x[],
                                                               const PylithInt numConstants,
                                                               const PylithScalar constants[],
                                                               PylithScalar g1[]) {

    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar pressure = s[sOff[i_pressure]];

    PylithScalar trace_strain = 0.0;
    for (PylithInt d = 0; d < _dim; ++d) {
     trace_strain += disp_x[d*_dim+d];
    }

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticity
    const PylithInt i_rstress = numA - 7;
    const PylithInt i_rstrain = numA - 6;
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];
    const PylithScalar* refStress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refStrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithScalar ref_trace_strain = refStrain[0] + refStrain[1] + refStrain[2];
    const PylithScalar meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
    const PylithScalar meanStress = meanRefStress + drainedBulkModulus * (trace_strain - ref_trace_strain);
    const PylithScalar alphaPres = biotCoefficient * pressure;
    const PylithReal traceTerm = (-2.0/3.0)*shearModulus * trace_strain;

    // Convert reference vectors to refrence tensors
    PylithScalar refStressTensor[_dim*_dim];
    PylithScalar refStrainTensor[_dim*_dim];
    PylithInt refTensorPos[9] = {0, 3, 5, 3, 1, 4, 5, 4, 2};

    for (PylithInt i = 0; i < _dim; ++i) {
      for (PylithInt j = 0; j < _dim; ++j) {
        refStressTensor[i*_dim+j] = refStress[refTensorPos[i*_dim+j]];
        refStrainTensor[i*_dim+j] = refStrain[refTensorPos[i*_dim+j]];
      } // for
    } // for

    for (PylithInt i = 0; i < _dim; ++i) {
        g1[i*_dim+i] -= (meanStress - alphaPres);
        g1[i*_dim+i] -= refStrainTensor[i*_dim+i] - meanRefStress + traceTerm;
        for (PylithInt j = 0; j < _dim; ++j) {
            g1[i*_dim+j] -= shearModulus * (disp_x[i*_dim+j] + disp_x[j*_dim+i]) - refStrainTensor[i*_dim+j];
        } // for
    } // for
} // g1v_refstate

// RHS Jacobians

// -----------------------------------------------------------------------------
// Jg2ue function for isotropic linear poroelasticity.

void
pylith::fekernels::IsotropicLinearPoroelasticity3D::Jg2ue(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal utshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jg2[]) {

    const PylithInt _dim = 3;

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    for (PylithInt d = 0; d < _dim; ++d) {
        Jg2[d*_dim+d] += drainedBulkModulus - (2.0*shearModulus) / 3.0;
    } // for
} // Jg2ue


// -----------------------------------------------------------------------------
// Jg2up function for isotropic linear poroelasticity.

void
pylith::fekernels::IsotropicLinearPoroelasticity3D::Jg2up(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal utshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jg2[]) {

    const PylithInt _dim = 3;

    // Isotropic Linear Poroelasticity
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    for (PylithInt d = 0; d < _dim; ++d) {
        Jg2[d*_dim+d] -= biotCoefficient;
    } // for
} // Jg2up

// -----------------------------------------------------------------------------
// Jg2vp function for isotropic linear poroelasticity.
// vp refers to dynamic formulation (velocity / pressure)
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::Jg2vp(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal utshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jg2[]) {

    const PylithInt _dim = 3;

    // Isotropic Linear Poroelasticity
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    for (PylithInt d = 0; d < _dim; ++d) {
        Jg2[d*_dim+d] -= biotCoefficient ;
    } // for
} // Jg2vp

// ---------------------------------------------------------------------------------------------------------------------
/* Jg3_uu entry function for isotropic linear poroelasticity WITHOUT reference stress and reference strain.
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::Jg3uu(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal s_tshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jg3[]) {

    const PylithInt _dim = 3;

    // index of Incoming auxiliary fields.
    // Poroelasticity

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithInt Nc = _dim;

    for (PylithInt c = 0; c < Nc; ++c) {
      for (PylithInt d = 0; d < _dim; ++d) {
        Jg3[((c*Nc + c)*_dim + d)*_dim + d] += shearModulus;
        Jg3[((c*Nc + d)*_dim + d)*_dim + c] += shearModulus;
      } // for
    } // for
} // Jg3uu

// ----------------------------------------------------------------------
/* Jg3pp entry function for isotropic linear poroelasticity. */
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::Jg3pp(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal utshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jg3[]) {

    const PylithInt _dim = 3;

    // index of Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_isotropicPermeability = numA - 1;

    const PylithScalar isotropicPermeablity = a[aOff[i_isotropicPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    for (PylithInt d = 0; d < _dim; ++d) {
      Jg3[d*_dim+d] -= isotropicPermeablity/fluidViscosity;
    }
} // Jg3pp

// ----------------------------------------------------------------------
/* Jg3pp entry function for isotropic linear poroelasticity, permeability in tensor form */
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::Jg3pp_tensor_permeability(const PylithInt dim,
                                                                            const PylithInt numS,
                                                                            const PylithInt numA,
                                                                            const PylithInt sOff[],
                                                                            const PylithInt sOff_x[],
                                                                            const PylithScalar s[],
                                                                            const PylithScalar s_t[],
                                                                            const PylithScalar s_x[],
                                                                            const PylithInt aOff[],
                                                                            const PylithInt aOff_x[],
                                                                            const PylithScalar a[],
                                                                            const PylithScalar a_t[],
                                                                            const PylithScalar a_x[],
                                                                            const PylithReal t,
                                                                            const PylithReal utshift,
                                                                            const PylithScalar x[],
                                                                            const PylithInt numConstants,
                                                                            const PylithScalar constants[],
                                                                            PylithScalar Jg3[]) {

    const PylithInt _dim = 3;

    // index of Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_tensorPerm = numA - 1;

    const PylithScalar* vectorPerm = &a[aOff[i_tensorPerm]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    PylithScalar tensorPerm[9];
    tensorPerm[0] = vectorPerm[0];
    tensorPerm[1] = vectorPerm[3];
    tensorPerm[2] = vectorPerm[5];
    tensorPerm[3] = vectorPerm[3];
    tensorPerm[4] = vectorPerm[1];
    tensorPerm[5] = vectorPerm[4];
    tensorPerm[6] = vectorPerm[5];
    tensorPerm[7] = vectorPerm[4];
    tensorPerm[8] = vectorPerm[2];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            Jg3[i*_dim+j] -= tensorPerm[i*_dim+j]/fluidViscosity;
        } // for
    } // for

} // Jg3pp_tensorPerm


// ---------------------------------------------------------------------------------------------------------------------
/* Jg3_vu entry function for isotropic linear poroelasticity WITHOUT reference stress and reference strain.
 *
 * stress_ij = C_ijkl strain_kl
 *
 * stress_11 = C1111 strain_11 + C1122 strain_22, C1111=lambda+2mu, C1122=lambda.
 *
 * stress_12 = C1212 strain_12 + C1221 strain_21. C1212 = C1221 from symmetry, so C1212 = C1221 = shearModulus.
 *
 * For reference:
 *
 * Isotropic:
 *  C_ijkl = bulkModulus * delta_ij * delta_kl
 *   + shearModulus * (delta_ik*delta_jl + delta_il*delta*jk - 2/3*delta_ij*delta_kl)
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::Jg3vu(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithReal t,
                                                        const PylithReal s_tshift,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar Jg3[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_disp = 0;
//    const PylithInt i_poro_pres = 1;

    // Incoming auxiliary fields.

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
//    const PylithScalar poro_pres = s[sOff[i_poro_pres]];

    PylithReal strainTrace = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
      strainTrace += disp_x[d*_dim+d];
    } // for

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    const PylithScalar meanStress = drainedBulkModulus * strainTrace;
    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar lambda2mu = lambda + 2.0*shearModulus;

    const PylithReal C1111 = lambda2mu;// Get auxiliary factory associated with physics.
    const PylithReal C2222 = lambda2mu;
    const PylithReal C1122 = lambda;
    const PylithReal C1212 = shearModulus;

    /* j(f,g,df,dg) = C(f,df,g,dg)
     *
     * 0: j0000 = C1111
     * 1: j0001 = C1112 = 0
     * 4: j0100 = C1121, symmetry C1112 = 0
     * 5: j0101 = C1122
     *
     * 2: j0010 = C1211 = 0
     * 3: j0011 = C1212
     * 6: j0110 = C1221, symmetry C1212
     * 7: j0111 = C1222 = 0
     *
     * 8: j1000 = C2111 = 0
     * 9: j1001 = C2112, symmetry C1212
     * 12: j1100 = C2121, symmetry C1212
     * 13: j1101 = C2122, symmetry C1222 = 0
     *
     * 10: j1010 = C2211, symmetry C1122
     * 11: j1011 = C2212, symmetry C1222 = 0
     * 14: j1110 = C2221, symmetry C1222 = 0
     * 15: j1111 = C2222
     */

    Jg3[ 0] -= C1111; // j0000
    Jg3[ 3] -= C1212; // j0011
    Jg3[ 5] -= C1122; // j0101
    Jg3[ 6] -= C1212; // j0110, C1221
    Jg3[ 9] -= C1212; // j1001, C2112
    Jg3[10] -= C1122; // j1010, C2211
    Jg3[12] -= C1212; // j1100, C2121
    Jg3[15] -= C2222; // j1111

} // Jg3vu

// ========================== Helper Kernels ===================================


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate stress vector for isotropic linear
 * poroelasticity WITHOUT a reference stress and strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim), ...]
 * Auxiliary fields: [density(1), ..., shear_modulus(1), bulk_modulus(1)]
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::cauchyStress(const PylithInt dim,
                                                               const PylithInt numS,
                                                               const PylithInt numA,
                                                               const PylithInt sOff[],
                                                               const PylithInt sOff_x[],
                                                               const PylithScalar s[],
                                                               const PylithScalar s_t[],
                                                               const PylithScalar s_x[],
                                                               const PylithInt aOff[],
                                                               const PylithInt aOff_x[],
                                                               const PylithScalar a[],
                                                               const PylithScalar a_t[],
                                                               const PylithScalar a_x[],
                                                               const PylithReal t,
                                                               const PylithScalar x[],
                                                               const PylithInt numConstants,
                                                               const PylithScalar constants[],
                                                               PylithScalar stressVector[]) {

    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar pressure = s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    // Incoming auxiliary fields.
    // Poroelasticity

    // IsotropicLinearPoroelasticity
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    // Create and populate stress tensor

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            stressTensor[i*_dim+j] += shearModulus * (disp_x[i*_dim+j] + disp_x[j*_dim+i]);
      } // for
      stressTensor[i*_dim+i] += (drainedBulkModulus - (2.0*shearModulus)/3.0) * trace_strain;
      stressTensor[i*_dim+i] -= biotCoefficient*pressure;
    } // for

    // Construct stress vector
    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar stress_zz = 0.5*lambda/(lambda+shearModulus) * (stressTensor[0*_dim+0] + stressTensor[1*_dim+1]);

    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stress_zz;
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
} // cauchyStress

// ---------------------------------------------------------------------------------------------------------------------
/* Calculate stress vector for isotropic linear
 * poroelasticity WITH a reference stress/strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim), ...]
 * Auxiliary fields: [density(1), ..., refstress(4), refstrain(4), shear_modulus(1), bulk_modulus(1)]
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity3D::cauchyStress_refstate(const PylithInt dim,
                                                                        const PylithInt numS,
                                                                        const PylithInt numA,
                                                                        const PylithInt sOff[],
                                                                        const PylithInt sOff_x[],
                                                                        const PylithScalar s[],
                                                                        const PylithScalar s_t[],
                                                                        const PylithScalar s_x[],
                                                                        const PylithInt aOff[],
                                                                        const PylithInt aOff_x[],
                                                                        const PylithScalar a[],
                                                                        const PylithScalar a_t[],
                                                                        const PylithScalar a_x[],
                                                                        const PylithReal t,
                                                                        const PylithScalar x[],
                                                                        const PylithInt numConstants,
                                                                        const PylithScalar constants[],
                                                                        PylithScalar stressVector[]) {

   const PylithInt _dim = 3;

   // Incoming solution fields.
   const PylithInt i_disp = 0;
   const PylithInt i_pressure = 1;
   const PylithInt i_trace_strain = 2;

   const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
   const PylithScalar pressure = s[sOff[i_pressure]];
   const PylithScalar trace_strain = s[sOff[i_trace_strain]];

   // Incoming auxiliary fields.

   // IsotropicLinearPoroelasticity
   const PylithInt i_rstress = numA - 7;
   const PylithInt i_rstrain = numA - 6;
   const PylithInt i_shearModulus = numA - 5;
   const PylithInt i_drainedBulkModulus = numA - 4;
   const PylithInt i_biotCoefficient = numA - 3;
   const PylithInt i_biotModulus = numA - 2;

   const PylithScalar shearModulus = a[aOff[i_shearModulus]];
   const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
   const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
   const PylithScalar biotModulus = a[aOff[i_biotModulus]];
   const PylithScalar* refStress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
   const PylithScalar* refStrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

   const PylithScalar ref_trace_strain = refStrain[0] + refStrain[1] + refStrain[2];
   const PylithScalar meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
   const PylithScalar meanStress = meanRefStress + drainedBulkModulus * (trace_strain - ref_trace_strain);
   const PylithScalar alphaPres = biotCoefficient * pressure;
   const PylithReal traceTerm = (-2.0/3.0)*shearModulus * trace_strain;

   // Convert reference vectors to refrence tensors
   PylithScalar refStressTensor[_dim*_dim];
   PylithScalar refStrainTensor[_dim*_dim];
   PylithInt refTensorPos[9] = {0, 3, 5, 3, 1, 4, 5, 4, 2};

   for (PylithInt i = 0; i < _dim; ++i) {
     for (PylithInt j = 0; j < _dim; ++j) {
       refStressTensor[i*_dim+j] = refStress[refTensorPos[i*_dim+j]];
       refStrainTensor[i*_dim+j] = refStrain[refTensorPos[i*_dim+j]];
     } // for
   } // for

   // Create and populate stress tensor

   PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

   for (PylithInt i = 0; i < _dim; ++i) {
       stressTensor[i*_dim+i] -= (meanStress - alphaPres);
       stressTensor[i*_dim+i] -= refStressTensor[i*_dim+i] - meanRefStress + traceTerm;
       for (PylithInt j = 0; j < _dim; ++j) {
           stressTensor[i*_dim+j] -= shearModulus * (disp_x[i*_dim+j] + disp_x[j*_dim+i]) - refStrainTensor[i*_dim+j];
       } // for
   } // for

   // Generate stress vector

    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar stress_zz = refStress[2] +  0.5*lambda/(lambda+shearModulus) *
                                       (stressTensor[0*_dim+0]-refStress[0] + stressTensor[1*_dim+1]-refStress[1]);

    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stress_zz;
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
} // cauchyStress_refstate

// End of file
