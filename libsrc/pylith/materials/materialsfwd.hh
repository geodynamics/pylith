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

/** @file libsrc/materials/materialsfwd.hh
 *
 * @brief Forward declarations for PyLith materials objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_materials_materialsfwd_hh)
#define pylith_materials_materialsfwd_hh

namespace pylith {
    namespace materials {
        // New stuff
        class Material;

        class Elasticity;
        class AuxiliaryFactoryElasticity;
        class DerivedFactoryElasticity;

        class RheologyElasticity;
        class IsotropicLinearElasticity;
        class IsotropicLinearMaxwell;
        class IsotropicLinearGenMaxwell;
        class IsotropicPowerLaw;
        class IsotropicLinearIncompElasticity;
        class AuxiliaryFactoryElastic;
        class AuxiliaryFactoryViscoelastic;

        class IncompressibleElasticity;

        class RheologyIncompressibleElasticity;
        class IsotropicLinearIncompElasticity;

        class Poroelasticity;
        class AuxiliaryFactoryPoroelasticity;
        class DerivedFactoryPoroelasticity;

        class RheologyPoroelasticity;
        class IsotropicLinearPoroelasticity;
        class AuxiliaryFactoryPoroelastic;

        class Query;

    } // materials
} // pylith

#endif // pylith_materials_materialsfwd_hh

// End of file
