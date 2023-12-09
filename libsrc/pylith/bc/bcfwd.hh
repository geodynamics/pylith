// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file libsrc/bc/bcfwd.hh
 *
 * @brief Forward declarations for PyLith boundary condition objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_bc_bcfwd_hh)
#define pylith_bc_bcfwd_hh

namespace pylith {
    namespace bc {
        class BoundaryCondition;
        class DiagnosticFieldFactory;
        class TimeDependentAuxiliaryFactory;

        class ConstraintBoundary;
        class Dirichlet;
        class DirichletTimeDependent;
        class DirichletUserFn;

        class NeumannTimeDependent;
        class NeumannUserFn;
        class AbsorbingDampers;

        class AbsorbingDampersAuxiliaryFactory;

    } // bc
} // pylith

#endif // pylith_bc_bcfwd_hh

// End of file
