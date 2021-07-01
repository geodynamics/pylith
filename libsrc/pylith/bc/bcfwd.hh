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
        class TimeDependentAuxiliaryFactory;

        class ConstraintBoundary;
        class Dirichlet;
        class DirichletTimeDependent;
        class DirichletUserFn;

        class NeumannTimeDependent;
        class AbsorbingDampers;

        class AbsorbingDampersAuxiliaryFactory;

    } // bc
} // pylith

#endif // pylith_bc_bcfwd_hh

// End of file
