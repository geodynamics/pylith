// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

namespace pylith {
    namespace feassemble {
        class AuxiliaryFactory; ///< Creates auxiliary subfields.

        class PhysicsImplementation; ///< Abstract base class for constraints and integrators.

        class DSLabelAccess; ///< Utility class for accessing PetscDMLabel, PetscDS, and PetscWeakForm.
        class FEKernelKey; ///< Utility class for managing keys for pointwise functions in finite-element integrations.
        class Integrator; ///< Abstract base class for finite-element integration.
        class IntegratorDomain; ///< Abstract base class for finite-element integration over portions on the domain.
        class IntegratorBoundary; ///< Abstract base class for finite-element integration over a boundary.
        class IntegratorInterface; ///< Abstract base class for finite-element integration over an interior interface.
        class IntegrationData; ///< Data used in finite-element integration (residual, solution, t, dt, ...)
        class InterfacePatches; ///< Interface integration patches.
        class UpdateStateVars; ///< Manager for updating state variables.
        class JacobianValues; ///< Manager for setting Jacobian values without finite-element integration.

        class Constraint; ///< Abstract base class for finite-element constraints.
        class ConstraintSpatialDB; ///< Finite-element constraints via auxiliary field from spatial database.
        class ConstraintUserFn; ///< Finite-element constraints via user-specified function (testing).
        class ConstraintSimple; ///< Finite-element simple constraints via user-specified function (testing).

    } // feassemble
} // pylith

// End of file
