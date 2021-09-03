// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/** @file libsrc/feassemble/feassemblefwd.hh
 *
 * @brief Forward declarations for PyLith feassemble objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_feassemble_feassemblefwd_hh)
#define pylith_feassemble_feassemblefwd_hh

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
        class InterfacePatches; ///< Interface integration patches.
        class UpdateStateVars; ///< Manager for updating state variables.
        class JacobianValues; ///< Manager for setting Jacobian values without finite-element integration.

        class Constraint; ///< Abstract base class for finite-element constraints.
        class ConstraintSpatialDB; ///< Finite-element constraints via auxiliary field from spatial database.
        class ConstraintUserFn; ///< Finite-element constraints via user-specified function (testing).
        class ConstraintSimple; ///< Finite-element simple constraints via user-specified function (testing).

    } // feassemble
} // pylith

#endif // pylith_feassemble_feassemblefwd_hh

// End of file
