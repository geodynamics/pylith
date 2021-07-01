/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University at Buffalo
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2021 University of California, Davis
 *
 * See LICENSE.md for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/Viscoelasticity.hh"

#include <cassert> // USES assert()

/* ======================================================================
 * Generic viscoelastic functions and kernels.
 * ======================================================================
 */

// ----------------------------------------------------------------------
// Function to compute Maxwell viscous strain coefficient.
PylithScalar
pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(const PylithScalar dt,
                                                              const PylithScalar maxwellTime) {
#if 0
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
        for (int iTerm = 2; iTerm <= numTerms; ++iTerm) {
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
    } else {
        // Default solution.
        dq = maxwellTime*(1.0-exp(-dt/maxwellTime))/dt;
        PetscLogFlops(6);
    } // else
#endif

    PylithScalar dq = maxwellTime*(1.0-exp(-dt/maxwellTime))/dt;

    return dq;
} // maxwellViscousStrainCoef


// Compute 2D scalar product of two tensors represented as vectors.
// 6 FLOPs per call.
PylithScalar
pylith::fekernels::Viscoelasticity::scalarProduct2D(const PylithScalar* tensor1, const PylithScalar* tensor2)
{ // scalarProduct2D
	const PylithScalar scalarProduct = tensor1[0] * tensor2[0] + tensor1[1] * tensor2[1] + 2.0 * tensor1[2] * tensor2[2];

	return scalarProduct;
} // scalarProduct2D


// Compute 2DPS scalar product of two tensors represented as vectors.
// In this case all 3 normal components are present, but only one shear comp.
// 8 FLOPs per call.
PylithScalar
pylith::fekernels::Viscoelasticity::scalarProduct2DPS(const PylithScalar* tensor1, const PylithScalar* tensor2)
{ // scalarProduct2DPS
	const PylithScalar scalarProduct = tensor1[0] * tensor2[0] + tensor1[1] * tensor2[1] + tensor1[2] * tensor2[2] + 
		2.0 * tensor1[3] * tensor2[3];

	return scalarProduct;
} // scalarProduct2DPS


// Compute 3D scalar product of two tensors represented as vectors.
// 12 FLOPs per call.
PylithScalar
pylith::fekernels::Viscoelasticity::scalarProduct3D(const PylithScalar* tensor1, const PylithScalar* tensor2)
{ // scalarProduct3D
	const PylithScalar scalarProduct = tensor1[0] * tensor2[0] + tensor1[1] * tensor2[1] + tensor1[2] * tensor2[2] +
		2.0 * (tensor1[3] * tensor2[3] + tensor1[4] * tensor2[4] + tensor1[5] * tensor2[5]);

	return scalarProduct;
} // scalarProduct3D


// End of file
