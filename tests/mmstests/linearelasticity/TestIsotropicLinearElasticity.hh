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
 * @file tests/mmstests/elasticity/TestIsotropicLinearElasticity.hh
 *
 * @brief C++ TestIsotropicLinearElasticity object
 *
 * C++ unit testing for IsotropicLinearElasticity.
 */

#if !defined(pylith_mmstests_testisotropiclinearelasticity_hh)
#define pylith_mmstests_testisotropiclinearelasticity_hh

#include "TestElasticity.hh" // ISA TestElasticity

#include "pylith/materials/IsotropicLinearElasticity.hh" // HOLDSA IsotropicLinearElasticity

namespace pylith {
    namespace mmstests {
        class TestIsotropicLinearElasticity;
    } // tests/mmstests
} // pylith

/// C++ unit testing for IsotropicLinearElasticity
class pylith::mmstests::TestIsotropicLinearElasticity : public pylith::mmstests::TestElasticity {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::materials::IsotropicLinearElasticity* _rheology; ///< Rheology for testing.

}; // class TestIsotropicLinearElasticity

#endif // pylith_mmstests_testisotropiclinearelasticity_hh

// End of file
