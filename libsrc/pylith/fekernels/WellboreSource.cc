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
 * Copyright (c) 2010-2015 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/WellboreSource.hh"

#include <cassert> // USES assert()
#include <iostream> // debugging.

// ----------------------------------------------------------------------
// f0 entry function for displacement equation: f0u = \dot{u}.
void
pylith::fekernels::WellboreSource::f0p(const PylithInt dim,
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
                                       PylithScalar f0[]) {
    assert(sOff);
    assert(s);
    assert(s_t);
    assert(f0);

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;

    // Incoming re-packed auxiliary field.
    const PylithInt i_fluidDensity = 0;
    const PylithInt i_fluidViscosity = 1;
    const PylithInt i_isotropicPermeability = 2;
    const PylithInt i_wellboreRadius = 3;
    const PylithInt i_wellboreLength = 4;
    const PylithInt i_wellborePressure = 5;
    const PylithInt i_wellboreCharacter = 6;
    const PylithInt i_elementDimensions = 7;

    const PylithScalar pressure = s[sOff[i_pressure]];

    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
    const PylithScalar wellboreRadius = a[aOff[i_wellboreRadius]];
    const PylithScalar wellboreLength = a[aOff[i_wellboreLength]];
    const PylithScalar wellborePressure = a[aOff[i_wellborePressure]];
    const PylithScalar wellboreCharacter = a[aOff[i_wellboreCharacter]];
    const PylithScalar* elementDimensions = &a[aOff[i_elementDimensions]];

    PylithScalar r_e = 0.28 * ( PetscSqrtReal(PetscSqrtReal(isotropicPermeability / isotropicPermeability) * PetscSqr(elementDimensions[0]) + PetscSqrtReal(isotropicPermeability/isotropicPermeability) * PetscSqr(elementDimensions[1]) ) /
                                ( PetscPowReal(isotropicPermeability/isotropicPermeability, 0.25) + PetscPowReal(isotropicPermeability/isotropicPermeability, 0.25) ) );
    // Well Index, L**3
    PylithScalar WI = 2.0*PETSC_PI * ( PetscSqrtReal(isotropicPermeability*isotropicPermeability)*wellboreLength / PetscLogReal(r_e / wellboreRadius) );

    f0[0] += WI * wellboreCharacter * (1.0 / fluidViscosity) * (pressure - wellborePressure);
} // f0p


// ----------------------------------------------------------------------
// Jf0 function for pressure equation
void
pylith::fekernels::WellboreSource::Jf0pp(const PylithInt dim,
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
                                         PylithScalar Jf0[]) {
    // Incoming re-packed auxiliary field.
    const PylithInt i_fluidDensity = 0;
    const PylithInt i_fluidViscosity = 1;
    const PylithInt i_isotropicPermeability = 2;
    const PylithInt i_wellboreRadius = 3;
    const PylithInt i_wellboreLength = 4;
    const PylithInt i_wellboreCharacter = 6;
    const PylithInt i_elementDimensions = 7;

    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
    const PylithScalar wellboreRadius = a[aOff[i_wellboreRadius]];
    const PylithScalar wellboreLength = a[aOff[i_wellboreLength]];
    const PylithScalar wellboreCharacter = a[aOff[i_wellboreCharacter]];
    const PylithScalar* elementDimensions = &a[aOff[i_elementDimensions]];

    PylithScalar r_e = 0.28 * ( PetscSqrtReal(PetscSqrtReal(isotropicPermeability / isotropicPermeability) * PetscSqr(elementDimensions[0]) + PetscSqrtReal(isotropicPermeability/isotropicPermeability) * PetscSqr(elementDimensions[1]) ) /
                                ( PetscPowReal(isotropicPermeability/isotropicPermeability, 0.25) + PetscPowReal(isotropicPermeability/isotropicPermeability, 0.25) ) );
    PylithScalar W = 2.0*PETSC_PI * ( PetscSqrtReal(isotropicPermeability*isotropicPermeability)*wellboreLength / PetscLogReal(r_e / wellboreRadius) );

    Jf0[0] += W * wellboreCharacter * (1.0 / fluidViscosity);

} // Jf0uu


// End of file
