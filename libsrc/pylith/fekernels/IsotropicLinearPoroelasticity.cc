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
// Kernels for isotropic, linear poroelasticity
// =====================================================================================================================
// ----------------------------------------------------------------------

// ================================= MMS =======================================

// ----------------------------------------------------------------------
// f0u function for quadratic space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticity::f0_mms_ql_u(const PylithInt dim,
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

 // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 // Incoming re-packed auxiliary field.

 // Poroelasticity

 // IsotropicLinearPoroelasticity
 const PylithInt i_shearModulus = numA - 5;
 const PylithInt i_undrainedBulkModulus = numA - 4;
 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar shearModulus = a[aOff[i_shearModulus]];
 const PylithScalar undrainedBulkModulus = a[aOff[i_undrainedBulkModulus]];
 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];
 const PylithScalar drainedBulkModulus = undrainedBulkModulus - biotCoefficient*biotCoefficient*biotModulus;
 const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

  for (PylithInt d = 0; d < dim-1; ++d) {
    f0u[d] -= 2.0*shearModulus - biotCoefficient*t;
  }
  f0u[dim-1] -= 2.0*lambda + 4.0*shearModulus - biotCoefficient*t;
} // f0_quadratic_linear_u

// ----------------------------------------------------------------------
// f0p function for quadratic space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticity::f0_mms_ql_p(const PylithInt dim,
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
  // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];

 PylithScalar       sum    = 0.0;

 for (PylithInt d = 0; d < dim; ++d) {
   sum += x[d];
 }
 f0p[0] += s_t ? biotCoefficient*s_t[sOff[i_trace_strain]] : 0.0;
 f0p[0] += s_t ? s_t[sOff[i_pressure]]/biotModulus     : 0.0;
 f0p[0] -= sum/biotModulus;
} // f0_quadratic_linear_p

// ----------------------------------------------------------------------
// f0u function for quadratic space and trigonometric time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticity::f0_mms_qt_u(const PylithInt dim,
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

 // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 // Incoming re-packed auxiliary field.

 // Poroelasticity

 // IsotropicLinearPoroelasticity
 const PylithInt i_shearModulus = numA - 5;
 const PylithInt i_undrainedBulkModulus = numA - 4;
 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar shearModulus = a[aOff[i_shearModulus]];
 const PylithScalar undrainedBulkModulus = a[aOff[i_undrainedBulkModulus]];
 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];
 const PylithScalar drainedBulkModulus = undrainedBulkModulus - biotCoefficient*biotCoefficient*biotModulus;
 const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

  for (PylithInt d = 0; d < dim-1; ++d) {
    f0u[d] -= 2.0*shearModulus - biotCoefficient*PetscCosReal(t);
  }
  f0u[dim-1] -= 2.0*lambda + 4.0*shearModulus - biotCoefficient*PetscCosReal(t);
} // f0_quadratic_trig_u

// ----------------------------------------------------------------------
// f0p function for quadratic space and trigonometric time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticity::f0_mms_qt_p(const PylithInt dim,
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
 // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];

 PylithScalar       sum    = 0.0;

 for (PylithInt d = 0; d < dim; ++d) {
   sum += x[d];
 }
 f0p[0] += s_t ? biotCoefficient*s_t[sOff[i_trace_strain]] : 0.0;
 f0p[0] += s_t ? s_t[sOff[i_pressure]]/biotModulus     : 0.0;
 f0p[0] += PetscSinReal(t)*sum/biotModulus;
} // f0_quadratic_trig_p

// ----------------------------------------------------------------------
// f0u function for trigonometric space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticity::f0_mms_tl_u(const PylithInt dim,
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
 // Incoming re-packed solution field.
 const PylithInt i_pressure = 1;
 const PylithInt i_trace_strain = 2;

 // Incoming re-packed auxiliary field.

 // Poroelasticity

 // IsotropicLinearPoroelasticity
 const PylithInt i_shearModulus = numA - 5;
 const PylithInt i_undrainedBulkModulus = numA - 4;
 const PylithInt i_biotCoefficient = numA - 3;
 const PylithInt i_biotModulus = numA - 2;

 const PylithScalar shearModulus = a[aOff[i_shearModulus]];
 const PylithScalar undrainedBulkModulus = a[aOff[i_undrainedBulkModulus]];
 const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
 const PylithScalar biotModulus = a[aOff[i_biotModulus]];
 const PylithScalar drainedBulkModulus = undrainedBulkModulus - biotCoefficient*biotCoefficient*biotModulus;
 const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

 for (PylithInt d = 0; d < dim-1; ++d) {
   f0u[d] += PetscSqr(2.*PETSC_PI)*PetscSinReal(2.*PETSC_PI*x[d])*(2.*shearModulus + lambda) + 2.0*(shearModulus + lambda) - 2.*PETSC_PI*biotCoefficient*PetscSinReal(2.*PETSC_PI*x[d])*t;
 }
 f0u[dim-1] += PetscSqr(2.*PETSC_PI)*PetscSinReal(2.*PETSC_PI*x[dim-1])*(2.*shearModulus + lambda) - 2.*PETSC_PI*biotCoefficient*PetscSinReal(2.*PETSC_PI*x[dim-1])*t;
} // f0_trig_linear_u

// ----------------------------------------------------------------------
// f0p function for trigonometric space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticity::f0_mms_tl_p(const PylithInt dim,
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

 for (PylithInt d = 0; d < dim; ++d) {
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
pylith::fekernels::IsotropicLinearPoroelasticity::f0p_quasistatic(const PylithInt dim,
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
pylith::fekernels::IsotropicLinearPoroelasticity::f0p_inertia(const PylithInt dim,
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

    for (PylithInt d = 0; d < dim; ++d) {
      trace_strain_t += vel_x[d*dim+d];
    }

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0p[0] += biotCoefficient*trace_strain_t;
    f0p[0] += poro_pres_t/biotModulus;
} // f0p_DYN

// -----------------------------------------------------------------------------
// Jf0pe function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticity::Jf0pe(const PylithInt dim,
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

    // Incoming re-packed auxiliary field.
    // IsotropicLinearPoroelasticity
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    Jf0[0] = utshift * biotCoefficient;
} // Jf0pe

// -----------------------------------------------------------------------------
void
pylith::fekernels::IsotropicLinearPoroelasticity::Jf0pp(const PylithInt dim,
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
pylith::fekernels::IsotropicLinearPoroelasticity::g1p_gravity(const PylithInt dim,
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

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 3;

    // IsotropicLinearPoroelasticity
    const PylithInt i_isotropicPerm = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];
    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar isotropicPerm = a[aOff[i_isotropicPerm]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    for (PylithInt d = 0; d < dim; ++d) {
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
pylith::fekernels::IsotropicLinearPoroelasticity::g1p_gravity_tensor_permeability(const PylithInt dim,
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
    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 3;

    // IsotropicLinearPoroelasticity
    const PylithInt i_tensorPerm = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];
    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar* tensorPerm = &a[aOff[i_tensorPerm]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    for (PylithInt i = 0; i < dim; ++i) {
        for (PylithInt j = 0; j < dim; j++) {
            g1p[i] -= (tensorPerm[i*dim+j] / fluidViscosity) * (pressure_x[j] - fluidDensity*gravityField[j]);
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
pylith::fekernels::IsotropicLinearPoroelasticity::g1p(const PylithInt dim,
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

    for (PylithInt d = 0; d < dim; ++d) {
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
pylith::fekernels::IsotropicLinearPoroelasticity::g1p_tensor_permeability(const PylithInt dim,
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
    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticity
    const PylithInt i_tensorPerm = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar* tensorPerm = &a[aOff[i_tensorPerm]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    for (PylithInt i = 0; i < dim; ++i) {
        for (PylithInt j = 0; j < dim; j++) {
            g1p[i] -= (tensorPerm[i*dim+j] / fluidViscosity) * (pressure_x[j]);
        } // for
    } // for
} // g1p_tensor_permeability

// ----------------------------------------------------------------------
// g1u function for isotropic linear Poroelasticity WITHOUT reference stress and strain.
// Quasi - Static Case
void
pylith::fekernels::IsotropicLinearPoroelasticity::g1u(const PylithInt dim,
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
    const PylithInt i_undrainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar undrainedBulkModulus = a[aOff[i_undrainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];
    const PylithScalar drainedBulkModulus = undrainedBulkModulus - biotCoefficient*biotCoefficient*biotModulus;

    for (PylithInt c = 0; c < dim; ++c) {
      for (PylithInt d = 0; d < dim; ++d) {
        g1[c*dim+d] += shearModulus * (disp_x[c*dim+d] + disp_x[d*dim+c]);
      } // for
      g1[c*dim+c] += (drainedBulkModulus - (2.0*shearModulus)/3.0) * trace_strain;
      g1[c*dim+c] -= biotCoefficient*pressure;
    } // for
} // g1u

// ----------------------------------------------------------------------
// g1u function for isotropic linear Poroelasticity plane strain with reference stress and strain.
void
pylith::fekernels::IsotropicLinearPoroelasticity::g1u_refstate(const PylithInt dim,
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_poro_pres = 1;
    const PylithInt i_trace_strain = 2;

    // Incoming auxiliary fields.

    // Poroelasticity

    // IsotropicLinearPoroelasticity
    const PylithInt i_rstress = numA - 7;
    const PylithInt i_rstrain = numA - 6;
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_undrainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithInt _numS = 3; // Number passed on to stress kernels.
    const PylithInt sOffCouple[3] = { sOff[i_disp], sOff[i_poro_pres], sOff_x[i_trace_strain] };
    const PylithInt sOffCouple_x[3] = { sOff_x[i_disp], sOff_x[i_poro_pres], sOff_x[i_trace_strain] };

    const PylithInt numAMean = 5; // Number passed to mean stress kernel.
    const PylithInt aOffMean[5] = { aOff[i_undrainedBulkModulus], aOff[i_biotCoefficient], aOff[i_biotModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[5] = { aOff_x[i_undrainedBulkModulus], aOff_x[i_biotCoefficient], aOff_x[i_biotModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    PylithScalar stressTensor[dim*dim];

    for (PylithInt d = 0; d < dim*dim; ++d) {
      stressTensor[d] = 0.0;
    }

    meanStress_refstate(dim, _numS, numAMean,sOffCouple, sOffCouple_x, s, s_t,
                        s_x,aOffMean, aOffMean_x, a, a_t, a_x,t, x, numConstants,
                        constants, stressTensor);

    deviatoricStress_refstate(dim, _numS, numADev,sOffCouple, sOffCouple_x, s,
                              s_t, s_x,aOffDev, aOffDev_x, a, a_t, a_x,t, x,
                              numConstants, constants, stressTensor);

    for (PylithInt d = 0; d < dim*dim; ++d) {
      g1[d] -= stressTensor[d];
    } // for
} // g1u_refstate


// ----------------------------------------------------------------------
// g1v function for isotropic linear Poroelasticity WITHOUT reference stress and strain.
// Dynamic Case
void
pylith::fekernels::IsotropicLinearPoroelasticity::g1v(const PylithInt dim,
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
   // Incoming solution fields.
   const PylithInt i_disp = 0;
   const PylithInt i_pressure = 1;

   const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
   const PylithScalar pressure = s[sOff[i_pressure]];

   PylithScalar trace_strain = 0.0;
   for (PylithInt d = 0; d < dim; ++d) {
     trace_strain += disp_x[d*dim+d];
   }

   // Incoming auxiliary fields.

   // IsotropicLinearPoroelasticity
   const PylithInt i_shearModulus = numA - 5;
   const PylithInt i_undrainedBulkModulus = numA - 4;
   const PylithInt i_biotCoefficient = numA - 3;
   const PylithInt i_biotModulus = numA - 2;

   const PylithScalar shearModulus = a[aOff[i_shearModulus]];
   const PylithScalar undrainedBulkModulus = a[aOff[i_undrainedBulkModulus]];
   const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
   const PylithScalar biotModulus = a[aOff[i_biotModulus]];
   const PylithScalar drainedBulkModulus = undrainedBulkModulus - biotCoefficient*biotCoefficient*biotModulus;

   for (PylithInt c = 0; c < dim; ++c) {
     for (PylithInt d = 0; d < dim; ++d) {
       g1[c*dim+d] += shearModulus * (disp_x[c*dim+d] + disp_x[d*dim+c]);
     } // for
     g1[c*dim+c] += (drainedBulkModulus - (2.0*shearModulus)/3.0) * trace_strain;
     g1[c*dim+c] -= biotCoefficient*pressure;
   } // for
} // g1v

// ----------------------------------------------------------------------
// g1v function for isotropic linear Poroelasticity plane strain with reference stress and strain.
void
pylith::fekernels::IsotropicLinearPoroelasticity::g1v_refstate(const PylithInt dim,
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_poro_pres = 1; ///SHOULDN'T THIS BE EQUAL TO 1 ??? (JOSIMAR)
    const PylithInt i_velocity = 2;

    // Incoming auxiliary fields.

    // Poroelasticity

    // IsotropicLinearPoroelasticity
    const PylithInt i_rstress = numA - 7;
    const PylithInt i_rstrain = numA - 6;
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_undrainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithInt _numS = 3; // Number passed on to stress kernels.
    const PylithInt sOffCouple[3] = { sOff[i_disp], sOff[i_poro_pres], sOff_x[i_velocity] };
    const PylithInt sOffCouple_x[3] = { sOff_x[i_disp], sOff_x[i_poro_pres], sOff_x[i_velocity] };

    const PylithInt numAMean = 5; // Number passed to mean stress kernel.
    const PylithInt aOffMean[5] = { aOff[i_undrainedBulkModulus], aOff[i_biotCoefficient], aOff[i_biotModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[5] = { aOff_x[i_undrainedBulkModulus], aOff_x[i_biotCoefficient], aOff_x[i_biotModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    PylithScalar stressTensor[dim*dim];

    for (PylithInt d = 0; d < dim*dim; ++d) {
      stressTensor[d] = 0.0;
    }

    meanStress_refstate(dim, _numS, numAMean,sOffCouple, sOffCouple_x, s, s_t,
                        s_x,aOffMean, aOffMean_x, a, a_t, a_x,t, x, numConstants,
                        constants, stressTensor);

    deviatoricStress_refstate(dim, _numS, numADev,sOffCouple, sOffCouple_x, s,
                              s_t, s_x,aOffDev, aOffDev_x, a, a_t, a_x,t, x,
                              numConstants, constants, stressTensor);

    for (PylithInt d = 0; d < dim*dim; ++d) {
      g1[d] -= stressTensor[d];
    } // for
} // g1v_refstate

// RHS Jacobians

// -----------------------------------------------------------------------------
// Jg2ue function for isotropic linear poroelasticity.

void
pylith::fekernels::IsotropicLinearPoroelasticity::Jg2ue(const PylithInt dim,
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

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_undrainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar undrainedBulkModulus = a[aOff[i_undrainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    const PylithScalar drainedBulkModulus = undrainedBulkModulus - biotCoefficient*biotCoefficient*biotModulus;

    for (PylithInt d = 0; d < dim; ++d) {
        Jg2[d*dim+d] += drainedBulkModulus - (2.0*shearModulus) / 3.0;
    } // for
} // Jg2ue


// -----------------------------------------------------------------------------
// Jg2up function for isotropic linear poroelasticity.

void
pylith::fekernels::IsotropicLinearPoroelasticity::Jg2up(const PylithInt dim,
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

    // Isotropic Linear Poroelasticity
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    for (PylithInt d = 0; d < dim; ++d) {
        Jg2[d*dim+d] -= biotCoefficient;
    } // for
} // Jg2up

// -----------------------------------------------------------------------------
// Jg2vp function for isotropic linear poroelasticity.
// vp refers to dynamic formulation (velocity / pressure)
void
pylith::fekernels::IsotropicLinearPoroelasticity::Jg2vp(const PylithInt dim,
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

    // Isotropic Linear Poroelasticity
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    for (PylithInt d = 0; d < dim; ++d) {
        Jg2[d*dim+d] -= biotCoefficient ;
    } // for
} // Jg2vp

// ---------------------------------------------------------------------------------------------------------------------
/* Jg3_uu entry function for isotropic linear poroelasticity WITHOUT reference stress and reference strain.
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity::Jg3uu(const PylithInt dim,
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
    // index of Incoming auxiliary fields.
    // Poroelasticity

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithInt Nc = dim;

    for (PylithInt c = 0; c < Nc; ++c) {
      for (PylithInt d = 0; d < dim; ++d) {
        Jg3[((c*Nc + c)*dim + d)*dim + d] += shearModulus;
        Jg3[((c*Nc + d)*dim + d)*dim + c] += shearModulus;
      } // for
    } // for
} // Jg3uu

// ----------------------------------------------------------------------
/* Jg3pp entry function for isotropic linear poroelasticity. */
void
pylith::fekernels::IsotropicLinearPoroelasticity::Jg3pp(const PylithInt dim,
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

    // index of Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_isotropicPermeability = numA - 1;

    const PylithScalar isotropicPermeablity = a[aOff[i_isotropicPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    for (PylithInt d = 0; d < dim; ++d) {
      Jg3[d*dim+d] -= isotropicPermeablity/fluidViscosity;
    }
} // Jg3pp

// ----------------------------------------------------------------------
/* Jg3pp entry function for isotropic linear poroelasticity, permeability in tensor form */
void
pylith::fekernels::IsotropicLinearPoroelasticity::Jg3pp_tensor_permeability(const PylithInt dim,
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
    // index of Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_tensorPermeability = numA - 1;

    const PylithScalar* tensorPermeablity = &a[aOff[i_tensorPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    for (PylithInt i = 0; i < dim; ++i) {
        for (PylithInt j = 0; j < dim; j++) {
            Jg3[i*dim+j] -= tensorPermeablity[i*dim+j]/fluidViscosity;
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
pylith::fekernels::IsotropicLinearPoroelasticity::Jg3vu(const PylithInt dim,
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
    //const PylithInt _dim = dim;

    // Incoming solution field.
    const PylithInt i_disp = 0;
//    const PylithInt i_poro_pres = 1;

    // Incoming auxiliary fields.

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_undrainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 6;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
//    const PylithScalar poro_pres = s[sOff[i_poro_pres]];

    PylithReal strainTrace = 0.0;

    for (PylithInt d = 0; d < dim; ++d) {
      strainTrace += disp_x[d*dim+d];
    } // for

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar undrainedBulkModulus = a[aOff[i_undrainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    const PylithScalar drainedBulkModulus = undrainedBulkModulus - biotCoefficient*biotCoefficient*biotModulus;
    const PylithScalar meanStress = drainedBulkModulus * strainTrace;
    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar lambda2mu = lambda + 2.0*shearModulus;

    if (dim == 2) {

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

    } else if (dim == 3) {

        // All other values are either zero or equal to one of these.
        const PylithReal C1111 = lambda2mu;
        const PylithReal C1122 = lambda;
        const PylithReal C1212 = shearModulus;
        /* j(f,g,df,dg) = C(f,df,g,dg)
         *
         *  0:  j0000 = C1111 = lambda + 2.0*shearModulus
         *  1:  j0001 = C1112 = 0
         *  2:  j0002 = C1113 = 0
         *  3:  j0010 = C1211 = 0
         *  4:  j0011 = C1212 = shearModulus
         *  5:  j0012 = C1213 = 0
         *  6:  j0020 = C1311 = 0
         *  7:  j0021 = C1312 = 0
         *  8:  j0022 = C1313 = shearModulus
         *  9:  j0100 = C1121 = 0
         * 10:  j0101 = C1122 = lambda
         * 11:  j0102 = C1123 = 0
         * 12:  j0110 = C1221 = shearModulus
         * 13:  j0111 = C1222 = 0
         * 14:  j0112 = C1223 = 0
         * 15:  j0120 = C1321 = 0
         * 16:  j0121 = C1322 = 0
         * 17:  j0122 = C1323 = 0
         * 18:  j0200 = C1131 = 0
         * 19:  j0201 = C1132 = 0
         * 20:  j0202 = C1133 = lambda
         * 21:  j0210 = C1231 = 0
         * 22:  j0211 = C1232 = 0
         * 23:  j0212 = C1233 = 0
         * 24:  j0220 = C1331 = shearModulus
         * 25:  j0221 = C1332 = 0
         * 26:  j0222 = C1333 = 0
         * 27:  j1000 = C2111 = 0
         * 28:  j1001 = C2112 = shearModulus
         * 29:  j1002 = C2113 = 0
         * 30:  j1010 = C2211 = lambda
         * 31:  j1011 = C2212 = 0
         * 32:  j1012 = C2213 = 0
         * 33:  j1020 = C2311 = 0
         * 34:  j1021 = C2312 = 0
         * 35:  j1022 = C2313 = 0
         * 36:  j1100 = C2121 = shearModulus
         * 37:  j1101 = C2122 = 0
         * 38:  j1102 = C2123 = 0
         * 39:  j1110 = C2221 = 0
         * 40:  j1111 = C2222 = lambda + 2.0*shearModulus
         * 41:  j1112 = C2223 = 0
         * 42:  j1120 = C2321 = 0
         * 43:  j1121 = C2322 = 0
         * 44:  j1122 = C2323 = shearModulus
         * 45:  j1200 = C2131 = 0
         * 46:  j1201 = C2132 = 0
         * 47:  j1202 = C2133 = 0
         * 48:  j1210 = C2231 = 0
         * 49:  j1211 = C2232 = 0
         * 50:  j1212 = C2233 = lambda
         * 51:  j1220 = C2331 = 0
         * 52:  j1221 = C2332 = shearModulus
         * 53:  j1222 = C2333 = 0
         * 54:  j2000 = C3111 = 0
         * 55:  j2001 = C3112 = 0
         * 56:  j2002 = C3113 = shearModulus
         * 57:  j2010 = C3211 = 0
         * 58:  j2011 = C3212 = 0
         * 59:  j2012 = C3213 = 0
         * 60:  j2020 = C3311 = lambda
         * 61:  j2021 = C3312 = 0
         * 62:  j2022 = C3313 = 0
         * 63:  j2100 = C3121 = 0
         * 64:  j2101 = C3122 = 0
         * 65:  j2102 = C3123 = 0
         * 66:  j2110 = C3221 = 0
         * 67:  j2111 = C3222 = 0
         * 68:  j2112 = C3223 = shearModulus
         * 69:  j2120 = C3321 = 0
         * 70:  j2121 = C3322 = lambda
         * 71:  j2122 = C3323 = 0
         * 72:  j2200 = C3131 = shearModulus
         * 73:  j2201 = C3132 = 0
         * 74:  j2202 = C3133 = 0
         * 75:  j2210 = C3231 = 0
         * 76:  j2211 = C3232 = shearModulus
         * 77:  j2212 = C3233 = 0
         * 78:  j2220 = C3331 = 0
         * 79:  j2221 = C3332 = 0
         * 80:  j2222 = C3333 = lambda + 2.0*shearModulus
         */

        // Nonzero Jacobian entries.
        Jg3[ 0] -= C1111; // j0000
        Jg3[ 4] -= C1212; // j0011
        Jg3[ 8] -= C1212; // j0022
        Jg3[10] -= C1122; // j0101
        Jg3[12] -= C1212; // j0110
        Jg3[20] -= C1122; // j0202
        Jg3[24] -= C1212; // j0220
        Jg3[28] -= C1212; // j1001
        Jg3[30] -= C1122; // j1010
        Jg3[36] -= C1212; // j1100
        Jg3[40] -= C1111; // j1111
        Jg3[44] -= C1212; // j1122
        Jg3[50] -= C1122; // j1212
        Jg3[52] -= C1212; // j1221
        Jg3[56] -= C1212; // j2002
        Jg3[60] -= C1122; // j2020
        Jg3[68] -= C1212; // j2112
        Jg3[70] -= C1122; // j2121
        Jg3[72] -= C1212; // j2200
        Jg3[76] -= C1212; // j2211
        Jg3[80] -= C1111; // j2222

    } //elseif
} // Jg3vu





// ========================== Helper Kernels ===================================

// ----------------------------------------------------------------------
/* Calculate mean stress for isotropic linear
 * poroelasticity WITHOUT reference stress and strain.
 *
 * meanStress = bulkModulus * strain_kk
 * alphaPres = biotCoefficient * p
 * stress += (meanStress - alphaPres) * delta_ij
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity::meanStress(const PylithInt dim,
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
                                                             PylithScalar stress[]) {


    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Incoming auxiliary field.

    // IsotropicLinearPoroelasticity
    const PylithInt i_undrainedBulkModulus = 0;
    const PylithInt i_biotCoefficient = 1;
    const PylithInt i_biotModulus = 2;

    const PylithScalar pressure = s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];
;
    const PylithScalar undrainedBulkModulus = a[aOff[i_undrainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];
    const PylithScalar drainedBulkModulus = undrainedBulkModulus - biotCoefficient*biotCoefficient*biotModulus;

    for (PylithInt i = 0; i < dim; ++i) {
        stress[i*dim+i] += (drainedBulkModulus*trace_strain - biotCoefficient*pressure);
    } // for
} // meanStress

// ----------------------------------------------------------------------
/* Calculate deviatoric stress for isotropic linear
 * poroelasticity WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * i==j
 * devStress_ii = 2*shearModulus*strain_ii - 2/3*shearModulus*strain_kk
 *
 * i!=j
 * devStress_ij = 2*shearModulus*strain_ij
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity::deviatoricStress(const PylithInt dim,
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
                                                                   PylithScalar stress[]) {
    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithInt i_trace_strain = 2;

    // Incoming auxiliary field.

    // IsotropicLinearPoroelasticity
    const PylithInt i_shearModulus = 0;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    for (PylithInt i = 0; i < dim; ++i) {
      for (PylithInt j = 0; j < dim; ++j) {
        stress[i*dim+j] += shearModulus * (disp_x[i*dim+j] + disp_x[j*dim+i]);
      } // for j
      stress[i*dim+i] -= (2.0/3.0)*shearModulus * trace_strain;
    } // for i

} // deviatoricStress

// ----------------------------------------------------------------------
/* Calculate mean stress for isotropic linear
* poroelasticity WITH reference stress and reference strain.
 *
 * meanStress = meanRefStress + bulkModulus * (strain_kk - refstrain_kk)
 * alphaPres = biotCoefficient * p
 * stress += (meanStress - alphaPres) * delta_ij
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity::meanStress_refstate(const PylithInt dim,
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
                                                                      PylithScalar stress[]) {


    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithInt i_poro_pres = 1;
    const PylithInt i_trace_strain = 2;

    // Incoming auxiliary field.

    // IsotropicLinearPoroelasticity
    const PylithInt i_undrainedBulkModulus = 0;
    const PylithInt i_biotCoefficient = 1;
    const PylithInt i_biotModulus = 2;
    const PylithInt i_rstress = 3;
    const PylithInt i_rstrain = 4;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar poro_pres = s[sOff[i_poro_pres]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    const PylithScalar undrainedBulkModulus = a[aOff[i_undrainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithScalar ref_trace_strain = refstrain[0] + refstrain[1] + refstrain[2];

    const PylithScalar drainedBulkModulus = undrainedBulkModulus - biotCoefficient*biotCoefficient*biotModulus;
    const PylithScalar meanrstress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithScalar meanStress = meanrstress + drainedBulkModulus * (trace_strain - ref_trace_strain);
    const PylithScalar alphaPres = biotCoefficient * poro_pres;

    for (PylithInt i = 0; i < dim; ++i) {
        stress[i*dim+i] += (meanStress - alphaPres);
    } // for
} // meanStress_refstate

// ----------------------------------------------------------------------
/* Calculate deviatoric stress for isotropic linear
* poroelasticity WITH reference stress and reference strain.
*
* devStress_ij = stress_ij - meanStress*delta_ij
*
* i==j
* devStress_ii = refstress_ii - meanRefstress + 2*shearModulus*(strain_ii - refstrain_ii) - 2/3*shearModulus*(strain_kk - refstrain_kk)
*
* i!=j
* devStress_ij = refstress_ij + 2*shearModulus*(strain_ij - refstrain_ij)
*/
void
pylith::fekernels::IsotropicLinearPoroelasticity::deviatoricStress_refstate(const PylithInt dim,
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
                                                                            PylithScalar stress[]) {

    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithInt i_trace_strain = 2;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_rstress = 1;
    const PylithInt i_rstrain = 2;

    // IsotropicLinearPoroelasticity
    const PylithInt i_shearModulus = 0;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithReal refstrainTrace = refstrain[0] + refstrain[1] + refstrain[2];
    const PylithReal meanrstress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithReal traceTerm = (-2.0/3.0)*shearModulus * trace_strain;

    for (PylithInt i = 0; i < dim; ++i) {
      for (PylithInt j = 0; j< dim; ++j) {
        stress[i*dim+j] += shearModulus * (disp_x[i*dim+j] + disp_x[j*dim+i]) - refstrain[i*dim+j];
        if (i == j) {
          stress[i*dim+j] += refstress[i*dim+j] - meanrstress + traceTerm;
        } // if
      } // for j
    } // for i
} // deviatoricStress_refstate

// ----------------------------------------------------------------------
/* Calculate mean stress for isotropic linear
 * poroelasticity WITHOUT reference stress and strain.
 *
 * meanStress = bulkModulus * strain_kk
 * alphaPres = biotCoefficient * p
 * stress += (meanStress - alphaPres) * delta_ij
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity::meanStress_inertia(const PylithInt dim,
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
                                                                     PylithScalar stress[]) {
    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithInt i_pressure = 1;
    //const PylithInt i_velocity = 2;

    // Incoming auxiliary field.

    // IsotropicLinearPoroelasticity
    const PylithInt i_undrainedBulkModulus = 0;
    const PylithInt i_biotCoefficient = 1;
    const PylithInt i_biotModulus = 2;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar pressure = s[sOff[i_pressure]];

    PylithScalar trace_strain = 0.0;

    for (PylithInt d = 0; d < dim; ++d) {
        trace_strain += disp_x[d*dim+d];
    }

    const PylithScalar undrainedBulkModulus = a[aOff[i_undrainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];
    const PylithScalar drainedBulkModulus = undrainedBulkModulus - biotCoefficient*biotCoefficient*biotModulus;

    for (PylithInt d = 0; d < dim; ++d) {
        stress[d*dim+d] += (drainedBulkModulus*trace_strain - biotCoefficient*pressure);
    } // for

} // meanStress_dyn

// ----------------------------------------------------------------------
/* Calculate deviatoric stress for isotropic linear
 * poroelasticity WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * i==j
 * devStress_ii = 2*shearModulus*strain_ii - 2/3*shearModulus*strain_kk
 *
 * i!=j
 * devStress_ij = 2*shearModulus*strain_ij
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity::deviatoricStress_inertia(const PylithInt dim,
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
                                                                           PylithScalar stress[]) {
    // Incoming solution field.
    const PylithInt i_disp = 0;

    // Incoming auxiliary field.

    // IsotropicLinearPoroelasticity
    const PylithInt i_shearModulus = 0;
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    PylithScalar trace_strain = 0.0;

    for (PylithInt d = 0; d < dim; ++d) {
      trace_strain += disp_x[d*dim+d];
    }

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    for (PylithInt i = 0; i < dim; ++i) {
      for (PylithInt j = 0; j < dim; ++j) {
        stress[i*dim+j] += shearModulus * (disp_x[i*dim+j] + disp_x[j*dim+i]);
      } // for j
      stress[i*dim+i] -= (2.0/3.0)*shearModulus * trace_strain;
    } // for i
} // deviatoricStress_dyn

// ----------------------------------------------------------------------
/* Calculate mean stress for isotropic linear
* poroelasticity WITH reference stress and reference strain.
 *
 * meanStress = meanRefStress + bulkModulus * (strain_kk - refstrain_kk)
 * alphaPres = biotCoefficient * p
 * stress += (meanStress - alphaPres) * delta_ij
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity::meanStress_inertia_refstate(const PylithInt dim,
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
                                                                              PylithScalar stress[]) {


    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithInt i_poro_pres = 1;
    const PylithInt i_velocity = 2;

    // Incoming auxiliary field.

    // IsotropicLinearPoroelasticity
    const PylithInt i_undrainedBulkModulus = 0;
    const PylithInt i_biotCoefficient = 1;
    const PylithInt i_biotModulus = 2;
    const PylithInt i_rstress = 3;
    const PylithInt i_rstrain = 4;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar poro_pres = s[sOff[i_poro_pres]];
    const PylithScalar trace_strain = disp_x[0*dim+0] + disp_x[1*dim+1] + disp_x[2*dim+2];

    const PylithScalar undrainedBulkModulus = a[aOff[i_undrainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithScalar ref_trace_strain = refstrain[0] + refstrain[1] + refstrain[2];

    const PylithScalar drainedBulkModulus = undrainedBulkModulus - biotCoefficient*biotCoefficient*biotModulus;
    const PylithScalar meanrstress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithScalar meanStress = meanrstress + drainedBulkModulus * (trace_strain - ref_trace_strain);
    const PylithScalar alphaPres = biotCoefficient * poro_pres;

    for (PylithInt i = 0; i < dim; ++i) {
        stress[i*dim+i] += (meanStress - alphaPres);
    } // for
} // meanStress_dyn_refstate

// ----------------------------------------------------------------------
/* Calculate deviatoric stress for isotropic linear
* poroelasticity WITH reference stress and reference strain.
*
* devStress_ij = stress_ij - meanStress*delta_ij
*
* i==j
* devStress_ii = refstress_ii - meanRefstress + 2*shearModulus*(strain_ii - refstrain_ii) - 2/3*shearModulus*(strain_kk - refstrain_kk)
*
* i!=j
* devStress_ij = refstress_ij + 2*shearModulus*(strain_ij - refstrain_ij)
*/
void
pylith::fekernels::IsotropicLinearPoroelasticity::deviatoricStress_inertia_refstate(const PylithInt dim,
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
                                                                                    PylithScalar stress[]) {

    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithInt i_velocity = 2;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_rstress = 1;
    const PylithInt i_rstrain = 2;

    // IsotropicLinearPoroelasticity
    const PylithInt i_shearModulus = 0;

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar trace_strain = disp_x[0*dim+0] + disp_x[1*dim+1] + disp_x[2*dim+2];
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithReal refstrainTrace = refstrain[0] + refstrain[1] + refstrain[2];
    const PylithReal meanrstress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithReal traceTerm = (-2.0/3.0)*shearModulus * trace_strain;

    for (PylithInt i = 0; i < dim; ++i) {
      for (PylithInt j = 0; j< dim; ++j) {
        stress[i*dim+j] += shearModulus * (disp_x[i*dim+j] + disp_x[j*dim+i]) - refstrain[i*dim+j];
        if (i == j){
          stress[i*dim+j] += refstress[i*dim+j] - meanrstress + traceTerm;
        } // if
      } // for j
    } // for i
} // deviatoricStress_dyn_refstate


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
pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress(const PylithInt dim,
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_poro_pres = 1; ///SHOULDN'T THIS BE EQUAL TO 1 ??? (JOSIMAR)
    //const PylithInt i_trace_strain = 2;

    // Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_rstress = numA - 7;
    const PylithInt i_rstrain = numA - 6;

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_poissonsRatio = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3 ;

    //const PylithInt _numS = 3; // Number passed on to stress kernels.
    //const PylithInt sOffCouple[3] = { sOff[i_disp], sOff[i_poro_pres], sOff[i_trace_strain]};
    //const PylithInt sOffCouple_x[3] = { sOff_x[i_disp], sOff_x[i_poro_pres], sOff_x[i_trace_strain]};

    const PylithInt _numS = 2; // Number passed on to stress kernels.
    const PylithInt sOffCouple[2] = { sOff[i_disp], sOff[i_poro_pres]};
    const PylithInt sOffCouple_x[2] = { sOff_x[i_disp], sOff_x[i_poro_pres]};

    const PylithInt numAMean = 5; // Number passed to mean stress kernel.
    const PylithInt aOffMean[5] = { aOff[i_shearModulus], aOff[i_poissonsRatio], aOff[i_biotCoefficient], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[5] = { aOff_x[i_shearModulus], aOff_x[i_poissonsRatio], aOff_x[i_biotCoefficient], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    if (dim == 1) {

        PylithScalar stressTensor[1] = { 0.0};
        meanStress(dim, _numS, numAMean, sOffCouple, sOffCouple_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                   t, x, numConstants, constants, stressTensor);
        deviatoricStress(dim, _numS, numADev, sOffCouple, sOffCouple_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                         t, x, numConstants, constants, stressTensor);
        stressVector[0] = stressTensor[0]; // stress_xx

    } else if (dim == 2) {

        PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
        meanStress(dim, _numS, numAMean, sOffCouple, sOffCouple_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                   t, x, numConstants, constants, stressTensor);
        deviatoricStress(dim, _numS, numADev, sOffCouple, sOffCouple_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                         t, x, numConstants, constants, stressTensor);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar poissonsRatio = a[aOff[i_poissonsRatio]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        const PylithScalar drainedBulkModulus = shearModulus * (2.0 * (1.0 + poissonsRatio) ) / (3.0 * (1.0 - 2.0*poissonsRatio) );
        const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;
        const PylithScalar stress_zz = 0.5*lambda/(lambda+shearModulus) * (stressTensor[0*dim+0] + stressTensor[1*dim+1]);

        stressVector[0] = stressTensor[0*dim+0]; // stress_xx
        stressVector[1] = stressTensor[1*dim+1]; // stress_yy
        stressVector[2] = stress_zz;
        stressVector[3] = stressTensor[0*dim+1]; // stress_xy

    } else if (dim == 3) {

        PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        meanStress(dim, _numS, numAMean, sOffCouple, sOffCouple_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                   t, x, numConstants, constants, stressTensor);
        deviatoricStress(dim, _numS, numADev, sOffCouple, sOffCouple_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                         t, x, numConstants, constants, stressTensor);

        stressVector[0] = stressTensor[0*dim+0]; // stress_xx
        stressVector[1] = stressTensor[1*dim+1]; // stress_yy
        stressVector[2] = stressTensor[2*dim+2]; // stress_zz
        stressVector[3] = stressTensor[0*dim+1]; // stress_xy
        stressVector[4] = stressTensor[1*dim+2]; // stress_yz
        stressVector[5] = stressTensor[0*dim+2]; // stress_xz
    } //elseif

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
pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress_refstate(const PylithInt dim,
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
    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_poro_pres = 1; ///SHOULDN'T THIS BE EQUAL TO 1 ??? (JOSIMAR)
    //const PylithInt i_trace_strain = 2;

    // Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_rstress = numA - 7;
    const PylithInt i_rstrain = numA - 6;

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_poissonsRatio = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3 ;

    //const PylithInt _numS = 3; // Number passed on to stress kernels.
    //const PylithInt sOffCouple[3] = { sOff[i_disp], sOff[i_poro_pres], sOff[i_trace_strain]};
    //const PylithInt sOffCouple_x[3] = { sOff_x[i_disp], sOff_x[i_poro_pres], sOff_x[i_trace_strain]};

    const PylithInt _numS = 2; // Number passed on to stress kernels.
    const PylithInt sOffCouple[2] = { sOff[i_disp], sOff[i_poro_pres]};
    const PylithInt sOffCouple_x[2] = { sOff_x[i_disp], sOff_x[i_poro_pres]};

    const PylithInt numAMean = 5; // Number passed to mean stress kernel.
    const PylithInt aOffMean[5] = { aOff[i_shearModulus], aOff[i_poissonsRatio], aOff[i_biotCoefficient], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[5] = { aOff[i_shearModulus], aOff_x[i_poissonsRatio], aOff_x[i_biotCoefficient], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    if (dim == 1) {

        PylithScalar stressTensor[1] = { 0.0};
        meanStress_refstate(dim, _numS, numAMean, sOffCouple, sOffCouple_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                   t, x, numConstants, constants, stressTensor);
        deviatoricStress_refstate(dim, _numS, numADev, sOffCouple, sOffCouple_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                         t, x, numConstants, constants, stressTensor);
        stressVector[0] = stressTensor[0]; // stress_xx

    } else if (dim == 2) {

        PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
        meanStress_refstate(dim, _numS, numAMean, sOffCouple, sOffCouple_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                            t, x, numConstants, constants, stressTensor);
        deviatoricStress_refstate(dim, _numS, numADev, sOffCouple, sOffCouple_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                                  t, x, numConstants, constants, stressTensor);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar poissonsRatio = a[aOff[i_poissonsRatio]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        const PylithScalar drainedBulkModulus = shearModulus * (2.0 * (1.0 + poissonsRatio) ) / (3.0 * (1.0 - 2.0*poissonsRatio) );
        const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;
        const PylithScalar* rstress = &a[aOff[i_rstress]];
        const PylithScalar stress_zz = rstress[2] +
                                       0.5*lambda/(lambda+shearModulus) *
                                       (stressTensor[0*dim+0]-rstress[0] + stressTensor[1*dim+1]-rstress[1]);

        stressVector[0] = stressTensor[0*dim+0]; // stress_xx
        stressVector[1] = stressTensor[1*dim+1]; // stress_yy
        stressVector[2] = stress_zz;
        stressVector[3] = stressTensor[0*dim+1]; // stress_xy

    } else if (dim == 3) {

        PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        meanStress_refstate(dim, _numS, numAMean, sOffCouple, sOffCouple_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                            t, x, numConstants, constants, stressTensor);
        deviatoricStress_refstate(dim, _numS, numADev, sOffCouple, sOffCouple_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                                  t, x, numConstants, constants, stressTensor);
        stressVector[0] = stressTensor[0*dim+0]; // stress_xx
        stressVector[1] = stressTensor[1*dim+1]; // stress_yy
        stressVector[2] = stressTensor[2*dim+2]; // stress_zz
        stressVector[3] = stressTensor[0*dim+1]; // stress_xy
        stressVector[4] = stressTensor[1*dim+2]; // stress_yz
        stressVector[5] = stressTensor[0*dim+2]; // stress_xz
    } //elseif

} // cauchyStress_refstate


// End of file
