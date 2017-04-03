/* -*- C -*-
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
 * Copyright (c) 2010-2015 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/elasticity.h"

/* ======================================================================
 * Generic elasticity kernels for inertia and body forces.
 *
 * Solution fields: [disp(dim), vel(dim)]
 *
 * ======================================================================
 */

/* ---------------------------------------------------------------------- */
/* f0 function for generic elasticity terms (inertia).
 */
void
pylith_fekernels_Elasticity_f0v_inertia(const PylithInt dim,
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
                                        PylithScalar f0[])
{ /* Elasticity_f0v_inertia */
    const PylithInt _numS = 2;
    const PylithInt i_vel = 1;
    const PylithScalar* vel_t = &s_t[sOff[i_vel]]; /* acceleration */

    const PylithInt _numA = 1;
    const PylithInt i_density = 0;
    const PylithScalar density = a[aOff[i_density]];

    PylithInt i;

    assert(_numS == numS);
    assert(_numA == numA);
    assert(sOff);
    assert(s_t);
    assert(aOff);
    assert(a);

    for (i=0; i < dim; ++i) {
        f0[i] += vel_t[i] * density;
    } /* for */
} /* Elasticity_f0v_inertia */


/* ---------------------------------------------------------------------- */
/* g0 function for generic elasticity terms (body forces).
 */
void
pylith_fekernels_Elasticity_g0v_grav(const PylithInt dim,
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
                                     PylithScalar g0[])
{ /* Elasticity_g0v_grav */
    const PylithInt _numS = 0;

    const PylithInt _numA = 2;
    const PylithInt i_density = 0;
    const PylithInt i_gravityField = 1;
    const PylithScalar density = a[aOff[i_density]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    PylithInt i;

    assert(_numS == numS);
    assert(_numA == numA);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_gravityField] >= 0);
    assert(a);

    for (i=0; i < dim; ++i) {
        g0[i] += density*gravityField[i];
    } /* for */
} /* Elasticity_g0v_grav */


/* ---------------------------------------------------------------------- */
/* g0 function for generic elasticity terms (body forces).
 */
void
pylith_fekernels_Elasticity_g0v_bodyforce(const PylithInt dim,
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
                                          PylithScalar g0[])
{ /* Elasticity_g0v_bodyforce */
    const PylithInt _numS = 2;

    const PylithInt _numA = 1;
    const PylithInt i_bodyForce = 0;
    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

    PylithInt i;

    assert(_numS == numS);
    assert(_numA == numA);
    assert(aOff);
    assert(aOff[i_bodyForce] >= 0);
    assert(a);

    for (i=0; i < dim; ++i) {
        g0[i] += bodyForce[i];
    } /* for */
} /* Elasticity_g0v_bodyforce */


/* ---------------------------------------------------------------------- */
/** Jf0 function for generic elasticity terms (inertia) with implicit time stepping.
 */
void
pylith_fekernels_Elasticity_Jf0vv_inertiaimplicit(const PylithInt dim,
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
                                                  PylithScalar Jf0[])
{ /* Elasticity_Jf0vv_inertiaimplicit */
    const PylithInt _numS = 2;

    const PylithInt _numA = 1;
    const PylithInt i_density = 0;
    const PylithScalar density = a[aOff[i_density]];

    PylithInt i, j;

    assert(_numS == numS);
    assert(_numA == numA);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(a);

    for (i=0; i < dim; ++i) {
        for (j=0; j < dim; ++j) {
            Jf0[i*dim+j] += utshift * density;
        } /* for */
    } /* for */
} /* Elasticity_Jf0vv_inertiaimplicit */




/** Jf0 function for generic elasticity terms (inertia) with explicit time stepping.
 */
void
pylith_fekernels_Elasticity_Jf0vv_inertiaexplicit(const PylithInt dim,
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
                                                  PylithScalar Jf0[])
{ /* Elasticity_Jf0vv_inertiaexplicit */
    const PylithInt _numS = 2;
    const PylithInt i_disp = 0;
    const PylithInt i_vel = 1;

    const PylithInt _numA = 1;
    const PylithInt i_density = 0;
    const PylithScalar density = a[aOff[i_density]];


    assert(_numS == numS);
    assert(_numA == numA);
    assert(aOff);
    assert(a);

    Jf0[i_disp*_numS+i_vel] += density;
} /* Elasticity_Jf0vv_inertiaexplicit */


/** Function to compute Maxwell viscous strain coefficient.
 */
PylithScalar
pylith_fekernels_Elasticity_Maxwell_VisStrain_Coeff(const PylithScalar dt,
						    const PylithScalar maxwellTime)
{ // VisStrain_Coeff
#if 0
  if (maxwellTime <= 0.0)
    throw std::runtime_error("Maxwell time must be greater than 0.");

  // Define cutoff values
  const PylithScalar timeFrac = 1.0e-10;

  // Compute viscous strain parameter.  The ratio of dt and
  // maxwellTime should never approach timeFrac for any reasonable
  // computation, but I have put in alternative solutions just in
  // case.

  // For now, assume time step size is reasonable to avoid if statements.
  PylithScalar dq = 0.0;

  // Use series expansion if dt is very small, since default solution
  // blows up otherwise.
  if (dt < timeFrac*maxwellTime) {
    PylithScalar fSign = 1.0;
    PylithScalar factorial = 1.0;
    PylithScalar fraction = 1.0;
    dq = 1.0;

    const int numTerms = 5;
    for (int iTerm=2; iTerm <= numTerms; ++iTerm) {
      factorial *= iTerm;
      fSign *= -1.0;
      fraction *= dt / maxwellTime;
      dq += fSign * fraction / factorial;
    } // for
    PetscLogFlops(8*(numTerms-1));
  } else if (maxwellTime < timeFrac*dt) {
    // Throw away exponential term if maxwellTime is very small.
    dq = maxwellTime / dt;
    PetscLogFlops(1);
  } else{
    // Default solution.
    dq = maxwellTime*(1.0-exp(-dt/maxwellTime))/dt;
    PetscLogFlops(6);
  } // else
#endif

  PylithScalar dq = maxwellTime*(1.0-exp(-dt/maxwellTime))/dt;

  return dq;
} // VisStrain_Coeff


/* End of file */
