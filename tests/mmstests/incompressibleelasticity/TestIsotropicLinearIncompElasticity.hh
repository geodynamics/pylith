// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

/**
 * @file tests/mmstests/elasticity/TestIsotropicLinearIncompElasticity.hh
 *
 * @brief C++ TestIsotropicLinearIncompElasticity object
 *
 * C++ unit testing for IsotropicLinearIncompElasticity.
 */

#if !defined(pylith_mmstests_testisotropiclinearincompelasticity_hh)
#define pylith_mmstests_testisotropiclinearincompelasticity_hh

#include "TestIncompressibleElasticity.hh" // ISA TestIncompressibleElasticity

#include "pylith/materials/IsotropicLinearIncompElasticity.hh" // HOLDSA IsotropicLinearIncompElasticity

namespace pylith {
    namespace mmstests {
        class TestIsotropicLinearIncompElasticity;
    } // tests/mmstests
} // pylith

/// C++ unit testing for IsotropicLinearIncompElasticity
class pylith::mmstests::TestIsotropicLinearIncompElasticity : public pylith::mmstests::TestIncompressibleElasticity {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::materials::IsotropicLinearIncompElasticity* _rheology; ///< Rheology for testing.

}; // class TestIsotropicLinearIncompElasticity

#endif // pylith_mmstests_testisotropiclinearincompelasticity_hh

// End of file
