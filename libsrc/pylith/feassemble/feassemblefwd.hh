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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
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
        class Integrator; ///< Abstract base class for finite-element integration.
        class IntegratorDomain; ///< Abstract base class for finite-element integration over portions on the domain.
        class IntegratorBoundary; ///< Abstract base class for finite-element integration over a boundary.
        class IntegratorInterface; ///< Abstract base class for finite-element integration over an interior interface.

        class Constraint; ///< Abstract base class for finite-element constraints.
        class ConstraintBoundary; ///< Abstract base class for finite-element constraints on a boundary.

        class IntegratorPointwise; ///< Integration of terms in governing equation.
        class ConstraintPointwise; ///< Constrained degrees of freedom.

        class AuxiliaryFactory; ///< Creates auxiliary subfields.

        class Observer; ///< Observer of subject.
        class Observers; ///< Subject being Observers.
        class IntegratorObserver; //< Observes integrator.

    } // feassemble
} // pylith

#endif // pylith_feassemble_feassemblefwd_hh

// End of file
